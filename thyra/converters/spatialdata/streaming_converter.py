# thyra/converters/spatialdata/streaming_converter.py

"""Streaming SpatialData converter with direct Zarr write.

This converter processes MSI data in a memory-efficient streaming manner:
- Two-pass approach: count non-zeros, then write directly to Zarr
- Writes directly to final output without scipy matrix in memory
- Memory stays bounded regardless of dataset size (~200MB for any size)
- Supports both CSR (row-wise) and CSC (column-wise) sparse formats
"""

import gc
import logging
import shutil
import tempfile
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
import zarr
from numpy.typing import NDArray
from scipy import sparse
from tqdm import tqdm

from .base_spatialdata_converter import SPATIALDATA_AVAILABLE, BaseSpatialDataConverter

if SPATIALDATA_AVAILABLE:
    import geopandas as gpd
    import xarray as xr
    from anndata import AnnData
    from spatialdata import SpatialData
    from spatialdata.models import Image2DModel, ShapesModel, TableModel
    from spatialdata.transformations import Affine, Identity

logger = logging.getLogger(__name__)


@contextmanager
def _suppress_upstream_warnings():
    """Suppress known upstream warnings from ome_zarr, zarr v3, and spatialdata."""
    with warnings.catch_warnings():
        # ome_zarr passes deprecated kwargs to da.to_zarr()
        warnings.filterwarnings(
            "ignore",
            message="Passing storage-related arguments",
            category=FutureWarning,
        )
        # zarr v3 does not recognize parquet files inside the store
        warnings.filterwarnings(
            "ignore", message="Object at.*is not recognized", category=UserWarning
        )
        # zarr v3 consolidated metadata not in spec yet
        warnings.filterwarnings(
            "ignore",
            message="Consolidated metadata is currently not",
            category=UserWarning,
        )
        yield


class StreamingSpatialDataConverter(BaseSpatialDataConverter):
    """Memory-efficient streaming converter for MSI data to SpatialData format.

    Uses a two-pass approach to keep memory bounded during processing:
    - Pass 1: Count non-zeros to build sparse matrix structure (indptr)
    - Pass 2: Write indices and data directly to Zarr

    For large datasets (>50GB), uses CSC format with memory-mapped files
    to handle datasets of any size with ~200MB RAM. For smaller datasets,
    uses the simpler CSR approach.

    Key advantages:
    - Memory stays bounded (~200MB) regardless of dataset size
    - Writes directly to Zarr without scipy matrix in memory
    - CSC format enables efficient ion image extraction
    - Handles 1M+ pixels without memory issues
    """

    # Threshold in GB above which PCS (Pre-calculated Scatter) method is used
    # for memory-efficient CSC conversion. Below this, standard COO approach is used.
    # Set to 30 GB to catch continuous mode datasets that would cause memory spikes
    # (the standard method can use 2-3x the dataset size in peak memory).
    PCS_SIZE_THRESHOLD_GB: float = 30.0

    def __init__(
        self,
        *args,
        chunk_size: int = 5000,
        temp_dir: Optional[Path] = None,
        use_csc: Union[bool, Literal["auto"]] = "auto",
        **kwargs,
    ):
        """Initialize streaming converter.

        Args:
            *args: Arguments passed to BaseSpatialDataConverter
            chunk_size: Number of spectra to process before writing to disk.
                Larger values use more memory but reduce I/O overhead.
                Default: 5000 spectra per chunk.
            temp_dir: Directory for temporary Zarr files. If None, uses
                system temp directory. Cleaned up after conversion.
            use_csc: Controls CSC format conversion method:
                - "auto" (default): Use PCS method for large datasets (>50GB),
                  standard method for smaller datasets
                - True: Always use PCS method for memory-efficient CSC
                - False: Use CSR format instead
            **kwargs: Keyword arguments passed to BaseSpatialDataConverter

        Note:
            Intensity thresholding (filtering noise below a minimum value) is
            handled at the reader level via the `intensity_threshold` parameter
            passed to the reader constructor.
        """
        kwargs["handle_3d"] = False  # Force 2D mode for now
        super().__init__(*args, **kwargs)

        self._chunk_size = chunk_size
        self._temp_dir = temp_dir
        self._cleanup_temp = temp_dir is None
        self._use_csc = use_csc
        self._zarr_store: Optional[zarr.Group] = None
        self._temp_path: Optional[Path] = None

    def _suppress_reader_progress(self) -> None:
        """Suppress progress output from reader during streaming passes."""
        setattr(self.reader, "_quiet_mode", True)

    def _estimate_output_size_gb(self) -> float:
        """Estimate the output dataset size in GB.

        Uses metadata to estimate the dense matrix size:
        n_pixels * n_mz_bins * 4 bytes (float32)

        Returns:
            Estimated size in gigabytes
        """
        metadata = self.reader.get_essential_metadata()

        # Get dimensions
        n_x, n_y, n_z = metadata.dimensions
        n_pixels = n_x * n_y * n_z

        # Estimate number of m/z bins after resampling
        if self._resampling_config:
            if isinstance(self._resampling_config, dict):
                n_mz_bins = self._resampling_config.get("target_bins") or 10000
            else:
                # ResamplingConfig dataclass
                n_mz_bins = (
                    getattr(self._resampling_config, "target_bins", None) or 10000
                )
        else:
            # Estimate from mass range with ~0.01 Da resolution
            min_mass, max_mass = metadata.mass_range
            n_mz_bins = int((max_mass - min_mass) / 0.01)

        # Dense matrix size in bytes (float32 = 4 bytes)
        dense_bytes = n_pixels * n_mz_bins * 4

        # Convert to GB
        size_gb = dense_bytes / (1024**3)

        logger.info(
            f"Estimated output size: {size_gb:.1f} GB "
            f"({n_pixels:,} pixels x {n_mz_bins:,} m/z bins)"
        )

        return size_gb

    def _should_use_pcs(self) -> bool:
        """Determine whether to use PCS method based on settings and dataset size.

        Returns:
            True if PCS method should be used, False otherwise
        """
        if self._use_csc is True:
            return True
        elif self._use_csc is False:
            return False
        else:  # "auto"
            estimated_size = self._estimate_output_size_gb()
            use_pcs = estimated_size > self.PCS_SIZE_THRESHOLD_GB
            if use_pcs:
                logger.info(
                    f"Dataset size ({estimated_size:.1f} GB) exceeds threshold "
                    f"({self.PCS_SIZE_THRESHOLD_GB} GB) - using PCS method"
                )
            else:
                logger.info(
                    f"Dataset size ({estimated_size:.1f} GB) below threshold "
                    f"({self.PCS_SIZE_THRESHOLD_GB} GB) - using standard method"
                )
            return use_pcs

    def convert(self) -> bool:
        """Stream-convert MSI data to SpatialData format.

        Overrides the base convert() method. Writes directly to final output
        Zarr without creating scipy sparse matrix, keeping memory bounded.

        For large datasets (>50GB by default), uses the Pre-calculated Scatter
        (PCS) approach for CSC format:
        - O(N) time complexity (no sorting required)
        - ~200 MB memory regardless of dataset size
        - CSC format for efficient ion image access

        For smaller datasets, uses the standard COO approach which is simpler.

        Returns:
            True if conversion was successful, False otherwise
        """
        try:
            # Initialize (loads metadata, mass axis, etc.)
            self._initialize_conversion()

            # Decide which method to use based on settings and dataset size
            if self._should_use_pcs():
                # Use no-cache approach for memory-efficient CSC conversion
                # This processes spectra twice but eliminates the ~200GB cache file
                result = self._convert_to_csc_no_cache()
                logger.info(
                    f"Zero-copy CSC conversion complete: {result['total_nnz']:,} non-zeros"
                )
            else:
                # Use standard COO approach for smaller datasets
                self._setup_temp_storage()
                coo_result = self._stream_build_coo()
                data_structures = self._create_data_structures_from_coo(coo_result)
                self._finalize_data(data_structures)
                self._save_output(data_structures)
                self._cleanup_temp_storage()
                logger.info("Zero-copy COO conversion complete")

            return True

        except Exception as e:
            logger.error(f"Error during zero-copy conversion: {e}")
            import traceback

            logger.error(f"Detailed traceback:\n{traceback.format_exc()}")
            return False

        finally:
            self.reader.close()

    def _setup_temp_storage(self) -> None:
        """Set up temporary Zarr storage for COO components."""
        if self._temp_dir is None:
            self._temp_path = Path(tempfile.mkdtemp(prefix="streaming_coo_"))
            self._cleanup_temp = True
        else:
            self._temp_path = self._temp_dir
            self._cleanup_temp = False

        logger.info(f"Temp storage: {self._temp_path}")

        # Create Zarr store for COO chunks
        zarr_path = self._temp_path / "coo_chunks.zarr"
        self._zarr_store = zarr.open_group(str(zarr_path), mode="w")

    def _cleanup_temp_storage(self) -> None:
        """Clean up temporary storage."""
        if self._cleanup_temp and self._temp_path is not None:
            try:
                shutil.rmtree(self._temp_path, ignore_errors=True)
                logger.debug(f"Cleaned up temp storage: {self._temp_path}")
            except Exception as e:
                logger.warning(f"Failed to cleanup temp storage: {e}")

    def _stream_build_coo(self) -> Dict[str, Any]:
        """Stream-build CSR matrix using two-pass direct Zarr write.

        This approach uses bounded memory by:
        1. Pass 1: Count nnz per row to size arrays and build indptr
        2. Pass 2: Write indices and data directly to Zarr

        Returns:
            Dictionary with matrix info and metadata
        """
        if self._dimensions is None:
            raise ValueError("Dimensions not initialized")
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis not initialized")
        if self._zarr_store is None:
            raise ValueError("Zarr store not initialized")

        n_x, n_y, n_z = self._dimensions
        n_rows = n_x * n_y * n_z
        n_cols = len(self._common_mass_axis)

        logger.info(
            f"Streaming CSR build (two-pass): {n_rows:,} pixels x {n_cols:,} cols"
        )

        self._suppress_reader_progress()
        total_spectra = self._get_total_spectra_count()

        # Pass 1: Count non-zeros and compute TIC/average spectrum
        pass1_result = self._coo_pass1_count_nonzeros(
            n_x, n_y, n_z, n_rows, n_cols, total_spectra
        )

        # Setup Zarr arrays
        indices_arr, data_arr = self._coo_setup_zarr_arrays(
            n_rows, n_cols, pass1_result["total_nnz"], pass1_result["indptr"]
        )

        # Pass 2: Write data to Zarr (position-aware using indptr)
        self._coo_pass2_write_data(
            indices_arr, data_arr, pass1_result["indptr"], total_spectra
        )

        logger.info(
            f"Pass 2 complete: {pass1_result['total_nnz']:,} non-zeros written to Zarr"
        )

        avg_spectrum = pass1_result["total_intensity"] / max(
            pass1_result["pixel_count"], 1
        )

        return {
            "total_nnz": pass1_result["total_nnz"],
            "n_rows": n_rows,
            "n_cols": n_cols,
            "tic_values": pass1_result["tic_values"],
            "avg_spectrum": avg_spectrum,
            "pixel_count": pass1_result["pixel_count"],
            "avg_spectrum_per_region": pass1_result.get("avg_spectrum_per_region"),
        }

    def _init_region_accumulators(self, n_cols: int) -> Tuple[
        Optional[Dict[tuple, int]],
        Optional[Dict[int, NDArray[np.float64]]],
        Optional[Dict[int, int]],
    ]:
        """Initialise per-region accumulation structures.

        Returns:
            (region_map, region_total, region_count) -- all None when no
            region map is available.
        """
        region_map = self._region_map if hasattr(self, "_region_map") else None
        if region_map is None:
            return None, None, None
        unique_regions = sorted(set(region_map.values()))
        region_total = {r: np.zeros(n_cols, dtype=np.float64) for r in unique_regions}
        region_count: Dict[int, int] = {r: 0 for r in unique_regions}
        return region_map, region_total, region_count

    @staticmethod
    def _accumulate_region(
        region_map: Dict[tuple, int],
        region_total: Dict[int, NDArray[np.float64]],
        region_count: Dict[int, int],
        x: int,
        y: int,
        mz_indices: NDArray[np.int_],
        resampled_ints: NDArray[np.float64],
    ) -> None:
        """Accumulate spectrum into the appropriate region bucket."""
        rn = region_map.get((x, y), -1)
        if rn in region_total:
            np.add.at(region_total[rn], mz_indices, resampled_ints)
            region_count[rn] += 1

    @staticmethod
    def _compute_region_averages(
        region_total: Optional[Dict[int, NDArray[np.float64]]],
        region_count: Optional[Dict[int, int]],
    ) -> Optional[Dict[str, NDArray[np.float64]]]:
        """Compute per-region mean spectra from accumulators."""
        if region_total is None or region_count is None:
            return None
        return {
            str(r): total / max(region_count.get(r, 0), 1)
            for r, total in region_total.items()
        }

    def _coo_pass1_count_nonzeros(
        self,
        n_x: int,
        n_y: int,
        n_z: int,
        n_rows: int,
        n_cols: int,
        total_spectra: int,
    ) -> Dict[str, Any]:
        """Pass 1: Count non-zeros per row and compute TIC/total intensity.

        Args:
            n_x: Number of pixels in x dimension.
            n_y: Number of pixels in y dimension.
            n_z: Number of pixels in z dimension.
            n_rows: Total number of rows (pixels).
            n_cols: Number of columns (m/z bins).
            total_spectra: Total number of spectra to process.

        Returns:
            Dictionary with nnz_per_row, total_nnz, tic_values, total_intensity,
            pixel_count, and indptr.
        """
        logger.info(
            f"Pass 1: Counting non-zeros per row ({total_spectra:,} spectra)..."
        )

        tic_values = np.zeros((n_y, n_x), dtype=np.float64)
        total_intensity = np.zeros(n_cols, dtype=np.float64)
        nnz_per_row = np.zeros(n_rows, dtype=np.int64)
        total_nnz = 0
        pixel_count = 0

        region_map, region_total, region_count = self._init_region_accumulators(n_cols)

        with tqdm(
            total=total_spectra, desc="Pass 1: Counting", unit="spectrum"
        ) as pbar:
            for coords, mzs, intensities in self.reader.iter_spectra(
                batch_size=self._buffer_size
            ):
                x, y, z = coords
                pixel_idx = z * (n_x * n_y) + y * n_x + x

                mz_indices, resampled_ints = self._process_spectrum(mzs, intensities)
                nnz = len(mz_indices)

                nnz_per_row[pixel_idx] = nnz
                total_nnz += nnz

                if 0 <= y < n_y and 0 <= x < n_x:
                    tic_values[y, x] = float(np.sum(resampled_ints))

                if nnz > 0:
                    np.add.at(total_intensity, mz_indices, resampled_ints)

                    if region_map is not None:
                        self._accumulate_region(
                            region_map,
                            region_total,
                            region_count,
                            x,
                            y,
                            mz_indices,
                            resampled_ints,
                        )

                pixel_count += 1
                pbar.update(1)

        logger.info(f"Pass 1 complete: {total_nnz:,} total non-zeros")

        # Build indptr from nnz counts
        indptr = np.zeros(n_rows + 1, dtype=np.int64)
        indptr[1:] = np.cumsum(nnz_per_row)
        del nnz_per_row
        gc.collect()

        return {
            "total_nnz": total_nnz,
            "tic_values": tic_values,
            "total_intensity": total_intensity,
            "pixel_count": pixel_count,
            "indptr": indptr,
            "avg_spectrum_per_region": self._compute_region_averages(
                region_total, region_count
            ),
        }

    def _coo_setup_zarr_arrays(
        self,
        n_rows: int,
        n_cols: int,
        total_nnz: int,
        indptr: NDArray[np.int64],
    ) -> Tuple[Any, Any]:
        """Setup CSR component arrays in Zarr store.

        Args:
            n_rows: Number of rows in the matrix.
            n_cols: Number of columns in the matrix.
            total_nnz: Total number of non-zero entries.
            indptr: Row pointers array.

        Returns:
            Tuple of (indices_arr, data_arr) Zarr arrays.
        """
        if self._zarr_store is None:
            raise RuntimeError("Zarr store not initialized")
        X_group = self._zarr_store.create_group("X")
        X_group.attrs["encoding-type"] = "csr_matrix"
        X_group.attrs["encoding-version"] = "0.1.0"
        X_group.attrs["shape"] = [n_rows, n_cols]

        # Use int64 for indptr if total_nnz exceeds int32 max
        indptr_dtype = np.int64 if total_nnz > np.iinfo(np.int32).max else np.int32
        indptr_arr = X_group.create_array("indptr", data=indptr.astype(indptr_dtype))
        indptr_arr.attrs["encoding-type"] = "array"
        indptr_arr.attrs["encoding-version"] = "0.2.0"

        chunk_size_zarr = min(total_nnz, 1000000)
        indices_arr = X_group.create_array(
            "indices",
            shape=(total_nnz,),
            dtype=np.int32,
            chunks=(chunk_size_zarr,),
        )
        indices_arr.attrs["encoding-type"] = "array"
        indices_arr.attrs["encoding-version"] = "0.2.0"
        data_arr = X_group.create_array(
            "data",
            shape=(total_nnz,),
            dtype=np.float64,
            chunks=(chunk_size_zarr,),
        )
        data_arr.attrs["encoding-type"] = "array"
        data_arr.attrs["encoding-version"] = "0.2.0"

        return indices_arr, data_arr

    def _coo_pass2_write_data(
        self,
        indices_arr: Any,
        data_arr: Any,
        indptr: NDArray[np.int64],
        total_spectra: int,
    ) -> None:
        """Pass 2: Write spectrum data to correct CSR positions.

        Each spectrum's data must be written at the position indicated by
        the indptr for its pixel_idx (row-major grid position), NOT in
        iteration order. The reader may yield spectra in frame order which
        differs from row-major pixel order.

        Args:
            indices_arr: Zarr array for column indices.
            data_arr: Zarr array for data values.
            indptr: Row pointers from Pass 1 (indexed by pixel position).
            total_spectra: Total number of spectra to process.
        """
        if self._dimensions is None:
            raise ValueError("Dimensions not initialized")

        n_x, n_y, n_z = self._dimensions

        logger.info("Pass 2: Writing data to Zarr (position-aware)...")

        if hasattr(self.reader, "reset"):
            self.reader.reset()
        self._suppress_reader_progress()

        # Track write position per row, starting at each row's indptr offset
        write_pos = indptr[:-1].copy()

        # Buffer writes for efficiency: collect (zarr_position, col_idx, value)
        buf_positions: list = []
        buf_indices: list = []
        buf_data: list = []
        buf_size = 0
        flush_threshold = 5_000_000  # flush every ~5M entries

        with tqdm(total=total_spectra, desc="Pass 2: Writing", unit="spectrum") as pbar:
            for coords, mzs, intensities in self.reader.iter_spectra(
                batch_size=self._buffer_size
            ):
                x, y, z = coords
                pixel_idx = z * (n_x * n_y) + y * n_x + x

                mz_indices, resampled_ints = self._process_spectrum(mzs, intensities)
                nnz = len(mz_indices)

                if nnz > 0:
                    pos = write_pos[pixel_idx]
                    positions = np.arange(pos, pos + nnz)
                    buf_positions.append(positions)
                    buf_indices.append(mz_indices.astype(np.int32))
                    buf_data.append(resampled_ints.astype(np.float64))
                    write_pos[pixel_idx] += nnz
                    buf_size += nnz

                if buf_size >= flush_threshold:
                    self._flush_positioned_to_zarr(
                        buf_positions,
                        buf_indices,
                        buf_data,
                        indices_arr,
                        data_arr,
                    )
                    buf_positions = []
                    buf_indices = []
                    buf_data = []
                    buf_size = 0

                pbar.update(1)

        # Flush remaining
        if buf_positions:
            self._flush_positioned_to_zarr(
                buf_positions,
                buf_indices,
                buf_data,
                indices_arr,
                data_arr,
            )

    def _flush_positioned_to_zarr(
        self,
        buf_positions: list,
        buf_indices: list,
        buf_data: list,
        indices_arr: Any,
        data_arr: Any,
    ) -> None:
        """Flush buffered data to correct positions in Zarr arrays.

        Sorts entries by position so contiguous ranges can be written
        efficiently in bulk.

        Args:
            buf_positions: List of position arrays (where to write).
            buf_indices: List of column index arrays.
            buf_data: List of data value arrays.
            indices_arr: Zarr array for column indices.
            data_arr: Zarr array for data values.
        """
        if not buf_positions:
            return

        all_pos = np.concatenate(buf_positions)
        all_idx = np.concatenate(buf_indices)
        all_dat = np.concatenate(buf_data)

        # Sort by position for efficient sequential writes
        order = np.argsort(all_pos)
        all_pos = all_pos[order]
        all_idx = all_idx[order]
        all_dat = all_dat[order]

        # Write contiguous runs in bulk
        n = len(all_pos)
        i = 0
        while i < n:
            start_pos = all_pos[i]
            # Find end of contiguous run
            j = i + 1
            while j < n and all_pos[j] == all_pos[j - 1] + 1:
                j += 1
            # Write the contiguous block
            indices_arr[int(start_pos) : int(start_pos) + (j - i)] = all_idx[i:j]
            data_arr[int(start_pos) : int(start_pos) + (j - i)] = all_dat[i:j]
            i = j

        del all_pos, all_idx, all_dat
        gc.collect()

    def _process_spectrum(
        self, mzs: np.ndarray, intensities: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Process a single spectrum - resample and return indices/values.

        Note: Intensity thresholding is handled at the reader level before data
        reaches this method. This method only handles resampling and zero filtering.

        Args:
            mzs: Mass values (already filtered by reader if threshold is set)
            intensities: Intensity values (already filtered by reader if threshold is set)

        Returns:
            Tuple of (mz_indices, resampled_intensities) with zeros filtered out
        """
        if not self._resampling_config:
            # No resampling - map m/z values to indices directly
            mz_indices = self._map_mass_to_indices(mzs)
            return mz_indices, intensities

        # Check for optimized nearest-neighbor path
        if hasattr(self, "_resampling_method"):
            from ...resampling import ResamplingMethod

            if self._resampling_method == ResamplingMethod.NEAREST_NEIGHBOR:
                mz_indices, resampled = self._nearest_neighbor_resample(
                    mzs, intensities
                )
                return mz_indices, resampled

        # Fallback: general resampling with zero filtering
        resampled_ints = self._resample_spectrum(mzs, intensities)
        if self._common_mass_axis is None:
            raise RuntimeError("Common mass axis not initialized")
        mz_indices = np.arange(len(self._common_mass_axis))
        mask = resampled_ints != 0
        return mz_indices[mask], resampled_ints[mask]

    def _create_data_structures_from_coo(
        self, coo_result: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Create SpatialData structures from CSR components in Zarr.

        Reads CSR components (indptr, indices, data) from Zarr and creates
        a scipy sparse matrix. This is more memory-efficient than the old
        COO approach because we read pre-built CSR components.

        Args:
            coo_result: Dictionary with matrix info and metadata

        Returns:
            Data structures dictionary for SpatialData creation
        """
        logger.info("Reading CSR components from Zarr...")

        n_rows = coo_result["n_rows"]
        n_cols = coo_result["n_cols"]

        # Read CSR components directly from Zarr
        if self._zarr_store is None:
            raise RuntimeError("Zarr store not initialized")
        X_group = self._zarr_store["X"]
        if not isinstance(X_group, zarr.Group):
            raise TypeError("Expected zarr.Group for X")
        indptr = np.asarray(X_group["indptr"])
        indices = np.asarray(X_group["indices"])
        data = np.asarray(X_group["data"])

        logger.info(f"Loaded CSR components: {len(data):,} entries")

        # For large datasets (>2.1B entries), ensure 64-bit indices for scipy
        if len(data) > np.iinfo(np.int32).max:
            logger.info("Large dataset detected, using 64-bit sparse matrix indices")
            indptr = indptr.astype(np.int64)
            indices = indices.astype(np.int64)

        # Create CSR matrix directly (no COO intermediate)
        sparse_matrix: Union[sparse.csr_matrix, sparse.csc_matrix] = sparse.csr_matrix(
            (data, indices, indptr),
            shape=(n_rows, n_cols),
            dtype=np.float64,
        )

        del indptr, indices, data
        gc.collect()

        # Convert to CSC if needed
        if self._sparse_format == "csc":
            logger.info("Converting CSR to CSC format...")
            sparse_matrix = sparse_matrix.tocsc()
            gc.collect()

        logger.info(
            f"Created sparse matrix: {sparse_matrix.shape}, {sparse_matrix.nnz:,} nnz"
        )

        # Build data structures
        if self._dimensions is None:
            raise ValueError("Dimensions not initialized")

        n_x, n_y, n_z = self._dimensions

        # Create slice data structure (similar to 2D converter)
        slice_id = f"{self.dataset_id}_z0"

        tables: Dict[str, Any] = {}
        shapes: Dict[str, Any] = {}
        images: Dict[str, Any] = {}

        return {
            "mode": "2d_slices",
            "slices_data": {
                slice_id: {
                    "sparse_matrix": sparse_matrix,
                    "coords_df": self._create_coordinates_dataframe_for_slice(0),
                    "tic_values": coo_result["tic_values"],
                }
            },
            "tables": tables,
            "shapes": shapes,
            "images": images,
            "var_df": self._create_mass_dataframe(),
            "avg_spectrum": coo_result["avg_spectrum"],
            "pixel_count": coo_result["pixel_count"],
        }

    def _create_data_structures(self) -> Dict[str, Any]:
        """Not used in streaming mode - required by ABC."""
        raise NotImplementedError("Streaming converter uses _stream_build_coo instead")

    def _create_coordinates_dataframe_for_slice(self, z_value: int) -> pd.DataFrame:
        """Create a coordinates dataframe for a single Z-slice.

        Args:
            z_value: Z-index of the slice

        Returns:
            DataFrame with pixel coordinates
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")

        n_x, n_y, _ = self._dimensions

        # Pre-allocate arrays for better performance
        pixel_count = n_x * n_y
        y_values: NDArray[np.int32] = np.repeat(np.arange(n_y, dtype=np.int32), n_x)
        x_values: NDArray[np.int32] = np.tile(np.arange(n_x, dtype=np.int32), n_y)
        instance_ids: NDArray[np.int32] = np.arange(pixel_count, dtype=np.int32)

        # Create DataFrame in one operation
        coords_df = pd.DataFrame(
            {
                "y": y_values,
                "x": x_values,
                "instance_id": instance_ids,
                "region": f"{self.dataset_id}_z{z_value}_pixels",
            }
        )

        # Set index efficiently
        coords_df["instance_id"] = coords_df["instance_id"].astype(str)
        coords_df.set_index("instance_id", inplace=True)

        # Add spatial coordinates in a vectorized operation
        coords_df["spatial_x"] = coords_df["x"] * self.pixel_size_um
        coords_df["spatial_y"] = coords_df["y"] * self.pixel_size_um

        # Always add per-pixel region numbers for a consistent schema.
        region_map = getattr(self, "_region_map", None)
        if region_map is not None:
            region_numbers = np.full(pixel_count, -1, dtype=np.int32)
            for i in range(pixel_count):
                key = (int(x_values[i]), int(y_values[i]))
                if key in region_map:
                    region_numbers[i] = region_map[key]
            coords_df["region_number"] = region_numbers
        else:
            coords_df["region_number"] = np.ones(pixel_count, dtype=np.int32)

        return coords_df

    def _finalize_data(self, data_structures: Dict[str, Any]) -> None:
        """Finalize data by creating tables, shapes, and images.

        Args:
            data_structures: Data structures containing processed data
        """
        if not SPATIALDATA_AVAILABLE:
            raise ImportError("SpatialData dependencies not available")

        avg_spectrum = data_structures["avg_spectrum"]
        self._non_empty_pixel_count = data_structures["pixel_count"]

        # Process each slice
        for slice_id, slice_data in data_structures["slices_data"].items():
            try:
                sparse_matrix = slice_data["sparse_matrix"]
                coords_df = slice_data["coords_df"]

                # Create AnnData for this slice
                adata = AnnData(
                    X=sparse_matrix,
                    obs=coords_df,
                    var=data_structures["var_df"],
                )

                # Add average spectrum to .uns
                adata.uns["average_spectrum"] = avg_spectrum

                # Add per-region mean spectra for multi-region datasets
                avg_per_region = data_structures.get("avg_spectrum_per_region")
                if avg_per_region is not None:
                    adata.uns["average_spectrum_per_region"] = avg_per_region

                # Add MSI metadata to .uns
                self._add_metadata_to_uns(adata)

                # Make sure region column exists and is correct
                region_key = f"{slice_id}_pixels"
                if "region" not in adata.obs.columns:
                    adata.obs["region"] = pd.Categorical([region_key] * len(adata))
                elif not isinstance(adata.obs["region"].dtype, pd.CategoricalDtype):
                    adata.obs["region"] = pd.Categorical(adata.obs["region"])

                # Make sure instance_key is a string column
                adata.obs["instance_key"] = adata.obs.index.astype(str)

                # Create table model
                table = TableModel.parse(
                    adata,
                    region=region_key,
                    region_key="region",
                    instance_key="instance_key",
                )

                # Add to tables and create shapes
                data_structures["tables"][slice_id] = table
                data_structures["shapes"][region_key] = self._create_pixel_shapes(
                    adata, is_3d=False
                )

                # Create TIC image for this slice
                tic_values = slice_data["tic_values"]
                y_size, x_size = tic_values.shape

                # Add channel dimension to make it (c, y, x) as required by SpatialData
                tic_values_with_channel = tic_values.reshape(1, y_size, x_size)

                # When alignment exists, use raster indices + Affine transform
                if self._tic_to_image_matrix is not None:
                    tic_image = xr.DataArray(
                        tic_values_with_channel,
                        dims=("c", "y", "x"),
                        coords={
                            "c": [0],
                            "y": np.arange(y_size),
                            "x": np.arange(x_size),
                        },
                    )
                    transform = Affine(
                        self._tic_to_image_matrix,
                        input_axes=("x", "y"),
                        output_axes=("x", "y"),
                    )
                else:
                    tic_image = xr.DataArray(
                        tic_values_with_channel,
                        dims=("c", "y", "x"),
                        coords={
                            "c": [0],
                            "y": np.arange(y_size) * self.pixel_size_um,
                            "x": np.arange(x_size) * self.pixel_size_um,
                        },
                    )
                    transform = Identity()

                # Create Image2DModel for the TIC image
                data_structures["images"][f"{slice_id}_tic"] = Image2DModel.parse(
                    tic_image,
                    transformations={
                        self.dataset_id: transform,
                        "global": transform,
                    },
                )

            except Exception as e:
                logger.error(f"Error processing slice {slice_id}: {e}")
                import traceback

                logger.debug(f"Detailed traceback:\n{traceback.format_exc()}")
                raise

        # Add optical images if available
        self._add_optical_images(data_structures)

    # ========================================================================
    # No-Cache CSC Conversion (Optimized)
    # Two-pass approach without disk caching - processes spectra twice
    # but eliminates ~200GB cache file I/O
    # ========================================================================

    def _convert_to_csc_no_cache(self) -> Dict[str, Any]:
        """Convert MSI data to CSC sparse format without disk caching.

        This optimized method processes spectra twice but eliminates the
        large cache file (~200GB for big datasets):

        1. Pre-scan: Count entries per column + compute TIC + average spectrum
           - Light pass: just resampling and counting, no disk I/O
        2. Allocate memory-mapped files for CSC arrays
        3. Main pass: Process spectra again and scatter directly to CSC
        4. Write CSC arrays to Zarr

        This works because nearest-neighbor resampling is deterministic -
        same input always produces same output. The 2x CPU cost of resampling
        is far less than the disk I/O cost of caching.

        Returns:
            Dictionary with conversion statistics
        """
        from uuid import uuid4

        if self._dimensions is None:
            raise ValueError("Dimensions not initialized")
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis not initialized")

        n_x, n_y, n_z = self._dimensions
        n_rows = n_x * n_y * n_z
        n_cols = len(self._common_mass_axis)

        logger.info(
            f"Streaming CSC (no-cache): {n_rows:,} pixels x {n_cols:,} m/z bins"
        )

        # Create temp directory for memmap files only (no cache file)
        temp_dir = self.output_path.parent / f".streaming_csc_{uuid4().hex[:8]}"
        temp_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Step 1: Pre-scan - count entries per column (no caching)
            logger.info("Step 1/3: Pre-scan (counting entries per column)...")
            prescan_result = self._prescan_count_columns(n_rows, n_cols, n_x, n_y)

            col_counts = prescan_result["col_counts"]
            total_nnz = prescan_result["total_nnz"]
            tic_values = prescan_result["tic_values"]
            avg_spectrum = prescan_result["avg_spectrum"]
            pixel_count = prescan_result["pixel_count"]
            avg_per_region = prescan_result.get("avg_spectrum_per_region")

            if total_nnz == 0:
                logger.warning("No non-zero entries found!")

            # Build indptr from col_counts
            indptr = np.zeros(n_cols + 1, dtype=np.int64)
            indptr[1:] = np.cumsum(col_counts)

            # Step 2: Allocate memory-mapped files for CSC arrays
            logger.info(f"Step 2/3: Allocating memmap ({total_nnz:,} entries)...")
            mm_indices, mm_data = self._allocate_csc_memmap_arrays(total_nnz, temp_dir)

            # Step 3: Main pass - process and scatter directly to CSC
            logger.info("Step 3/3: Processing spectra and scattering to CSC...")
            self._scatter_spectra_direct(mm_indices, mm_data, indptr, n_x, pixel_count)

            # Write CSC arrays to Zarr
            logger.info("Writing CSC arrays to Zarr...")
            self._write_csc_arrays_to_zarr(
                mm_indices,
                mm_data,
                indptr,
                n_rows,
                n_cols,
                total_nnz,
                tic_values,
                avg_spectrum,
                pixel_count,
                avg_per_region,
            )

            # Cleanup memmap references before deleting files
            del mm_indices, mm_data
            gc.collect()

            logger.info(f"Streaming CSC (no-cache) complete: {total_nnz:,} non-zeros")

            return {
                "total_nnz": total_nnz,
                "n_rows": n_rows,
                "n_cols": n_cols,
                "pixel_count": pixel_count,
            }

        finally:
            # Clean up temp directory
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)

    def _prescan_count_columns(
        self, n_rows: int, n_cols: int, n_x: int, n_y: int
    ) -> Dict[str, Any]:
        """Pre-scan spectra to count entries per column without caching.

        This is a lightweight pass that:
        1. Counts how many entries each m/z column will have (for CSC indptr)
        2. Computes TIC values per pixel
        3. Accumulates total intensity for average spectrum

        No data is cached to disk - we'll reprocess spectra in the main pass.

        Args:
            n_rows: Total number of pixels
            n_cols: Number of m/z bins
            n_x, n_y: Spatial dimensions

        Returns:
            Dictionary with col_counts, total_nnz, tic_values, avg_spectrum, pixel_count
        """
        # Allocate counting arrays (very small memory footprint)
        col_counts = np.zeros(n_cols, dtype=np.int64)
        total_intensity = np.zeros(n_cols, dtype=np.float64)
        tic_values = np.zeros((n_y, n_x), dtype=np.float64)

        total_nnz = 0
        pixel_count = 0

        region_map, region_total, region_count = self._init_region_accumulators(n_cols)

        self._suppress_reader_progress()
        total_spectra = self._get_total_spectra_count()

        with tqdm(
            total=total_spectra,
            desc="Pre-scan",
            unit="spectrum",
        ) as pbar:
            for coords, mzs, intensities in self.reader.iter_spectra(
                batch_size=self._buffer_size
            ):
                x, y, z = coords

                # Process spectrum (same resampling as main pass)
                mz_indices, resampled_ints = self._process_spectrum(mzs, intensities)

                nnz = len(mz_indices)
                if nnz > 0:
                    # Count entries per column (vectorized)
                    np.add.at(col_counts, mz_indices, 1)
                    total_nnz += nnz

                    # Accumulate for average spectrum (vectorized)
                    np.add.at(total_intensity, mz_indices, resampled_ints)

                    # TIC
                    if 0 <= y < n_y and 0 <= x < n_x:
                        tic_values[y, x] = resampled_ints.sum()

                    if region_map is not None:
                        self._accumulate_region(
                            region_map,
                            region_total,
                            region_count,
                            x,
                            y,
                            mz_indices,
                            resampled_ints,
                        )

                pixel_count += 1
                pbar.update(1)

        # Compute average spectrum
        avg_spectrum = total_intensity / max(pixel_count, 1)

        logger.info(
            f"  Pre-scan complete: {total_nnz:,} entries across {n_cols:,} columns"
        )

        return {
            "col_counts": col_counts,
            "total_nnz": total_nnz,
            "tic_values": tic_values,
            "avg_spectrum": avg_spectrum,
            "pixel_count": pixel_count,
            "avg_spectrum_per_region": self._compute_region_averages(
                region_total, region_count
            ),
        }

    def _scatter_spectra_direct(
        self,
        mm_indices: np.memmap,
        mm_data: np.memmap,
        indptr: NDArray[np.int64],
        n_x: int,
        pixel_count: int,
    ) -> None:
        """Process spectra and scatter directly to CSC arrays.

        This is the main pass that processes spectra again (same resampling
        as pre-scan) and scatters values directly to their CSC positions.

        Args:
            mm_indices: Memory-mapped array for row indices
            mm_data: Memory-mapped array for values
            indptr: Column pointers array (from pre-scan)
            n_x: Number of columns in spatial grid
            pixel_count: Number of spectra to process
        """
        # Current write position for each column
        write_pos = indptr[:-1].copy()

        # Periodic flush to reduce memory pressure
        flush_interval = 100_000
        spectra_since_flush = 0

        # Reset reader for second pass
        # For real readers (ImzML, Bruker), iter_spectra() is a generator factory
        # that creates a fresh iterator each time - no need to recreate the reader.
        # For mock readers with random data, we need to reset the random seed.
        if hasattr(self.reader, "reset"):
            self.reader.reset()

        self._suppress_reader_progress()

        with tqdm(
            total=pixel_count,
            desc="Scatter to CSC",
            unit="spectrum",
        ) as pbar:
            for coords, mzs, intensities in self.reader.iter_spectra(
                batch_size=self._buffer_size
            ):
                x, y, z = coords

                # Process spectrum (deterministic - same result as pre-scan)
                mz_indices, resampled_ints = self._process_spectrum(mzs, intensities)

                if len(mz_indices) > 0:
                    # Compute row index (row-major order)
                    row_idx = y * n_x + x

                    # Vectorized scatter
                    destinations = write_pos[mz_indices]
                    mm_indices[destinations] = row_idx
                    mm_data[destinations] = resampled_ints

                    # Increment write positions
                    write_pos[mz_indices] += 1

                pbar.update(1)
                spectra_since_flush += 1

                # Periodic flush
                if spectra_since_flush >= flush_interval:
                    mm_indices.flush()
                    mm_data.flush()
                    spectra_since_flush = 0

        # Final flush
        mm_indices.flush()
        mm_data.flush()

        logger.info("  Scatter complete, memmap flushed")

    def _allocate_csc_memmap_arrays(
        self, total_nnz: int, temp_dir: Path
    ) -> Tuple[np.memmap, np.memmap]:
        """Allocate memory-mapped files for CSC indices and data arrays.

        Uses numpy memmap so the OS handles memory management via virtual memory.
        This keeps RAM usage minimal regardless of array size.

        Args:
            total_nnz: Total number of non-zero entries
            temp_dir: Directory for temporary files

        Returns:
            Tuple of (mm_indices, mm_data) memory-mapped arrays
        """
        # Ensure at least 1 element for empty datasets
        size = max(total_nnz, 1)

        # CSC indices (row indices for each non-zero) - int32 sufficient for rows
        mm_indices = np.memmap(
            temp_dir / "csc_indices.bin",
            dtype=np.int32,
            mode="w+",
            shape=(size,),
        )

        # CSC data (values for each non-zero)
        mm_data = np.memmap(
            temp_dir / "csc_data.bin",
            dtype=np.float64,
            mode="w+",
            shape=(size,),
        )

        logger.info(
            f"  Allocated memmap: {size * 4 / (1024**3):.2f} GB indices + "
            f"{size * 8 / (1024**3):.2f} GB data"
        )

        return mm_indices, mm_data

    def _write_csc_arrays_to_zarr(
        self,
        mm_indices: np.memmap,
        mm_data: np.memmap,
        indptr: NDArray[np.int64],
        n_rows: int,
        n_cols: int,
        total_nnz: int,
        tic_values: NDArray[np.float64],
        avg_spectrum: NDArray[np.float64],
        pixel_count: int,
        avg_spectrum_per_region: dict[str, NDArray[np.float64]] | None = None,
    ) -> None:
        """Write CSC arrays to Zarr store with SpatialData-compatible structure.

        Creates the complete Zarr directory structure including:
        - CSC sparse matrix (indptr, indices, data)
        - Observation metadata (coordinates, region keys)
        - Variable metadata (m/z values)
        - Essential metadata and average spectrum

        Reads from memmap sequentially and writes in aligned chunks for efficiency.

        Args:
            mm_indices: Memory-mapped CSC indices array.
            mm_data: Memory-mapped CSC data array.
            indptr: Column pointers for CSC format.
            n_rows: Number of rows in the matrix.
            n_cols: Number of columns in the matrix.
            total_nnz: Total non-zero count.
            tic_values: TIC image array.
            avg_spectrum: Average spectrum.
            pixel_count: Number of non-empty pixels.
            avg_spectrum_per_region: Per-region mean spectra, or None.
        """
        from datetime import datetime

        from ... import __version__

        slice_id = f"{self.dataset_id}_z0"
        region_key = f"{slice_id}_pixels"
        if self._dimensions is None:
            raise ValueError("Dimensions not initialized")
        n_x, n_y, n_z = self._dimensions

        # Clean output directory
        if self.output_path.exists():
            shutil.rmtree(self.output_path)

        # Create Zarr store
        store = zarr.open_group(str(self.output_path), mode="w")

        # Root attributes
        store.attrs["spatialdata_attrs"] = {
            "version": "0.2",
            "spatialdata_software_version": "0.6.1",
        }
        store.attrs["pixel_size_x_um"] = float(self.pixel_size_um)
        store.attrs["pixel_size_y_um"] = float(self.pixel_size_um)
        store.attrs["pixel_size_units"] = "micrometers"
        store.attrs["coordinate_system"] = "physical_micrometers"
        store.attrs["msi_converter_version"] = __version__
        store.attrs["conversion_timestamp"] = datetime.now().isoformat()

        # Create table structure
        tables_group = store.create_group("tables")
        table_group = tables_group.create_group(slice_id)

        table_group.attrs["encoding-type"] = "anndata"
        table_group.attrs["encoding-version"] = "0.1.0"
        table_group.attrs["spatialdata-encoding-type"] = "ngff:regions_table"
        table_group.attrs["region"] = region_key
        table_group.attrs["region_key"] = "region"
        table_group.attrs["instance_key"] = "instance_key"
        table_group.attrs["version"] = "0.2"

        # Empty groups required by AnnData
        for group_name in ["layers", "obsm", "obsp", "varm", "varp"]:
            g = table_group.create_group(group_name)
            g.attrs["encoding-type"] = "dict"
            g.attrs["encoding-version"] = "0.1.0"

        raw_arr = table_group.create_array("raw", data=np.array(False))
        raw_arr.attrs["encoding-type"] = "array"
        raw_arr.attrs["encoding-version"] = "0.2.0"

        # X group (CSC matrix)
        X_group = table_group.create_group("X")
        X_group.attrs["encoding-type"] = "csc_matrix"
        X_group.attrs["encoding-version"] = "0.1.0"
        X_group.attrs["shape"] = [n_rows, n_cols]

        # Write indptr (small, do it directly) - use int64
        indptr_arr = X_group.create_array(
            "indptr",
            data=indptr.astype(np.int64),
            chunks=(min(len(indptr), 100_000),),
        )
        indptr_arr.attrs["encoding-type"] = "array"
        indptr_arr.attrs["encoding-version"] = "0.2.0"

        # Zarr chunk settings
        z_chunk = 1_000_000
        actual_nnz = max(total_nnz, 1)

        indices_arr = X_group.create_array(
            "indices",
            shape=(actual_nnz,),
            dtype=np.int32,
            chunks=(min(actual_nnz, z_chunk),),
        )
        indices_arr.attrs["encoding-type"] = "array"
        indices_arr.attrs["encoding-version"] = "0.2.0"
        data_arr = X_group.create_array(
            "data",
            shape=(actual_nnz,),
            dtype=np.float64,
            chunks=(min(actual_nnz, z_chunk),),
        )
        data_arr.attrs["encoding-type"] = "array"
        data_arr.attrs["encoding-version"] = "0.2.0"

        # Sequential transfer from memmap to Zarr (aligned to chunks)
        read_buffer_size = z_chunk * 50  # ~200 MB buffer

        with tqdm(
            total=total_nnz,
            desc="Step 4/4: Writing to Zarr",
            unit="entries",
            unit_scale=True,
        ) as pbar:
            for start in range(0, total_nnz, read_buffer_size):
                end = min(start + read_buffer_size, total_nnz)
                indices_arr[start:end] = mm_indices[start:end]
                data_arr[start:end] = mm_data[start:end]
                pbar.update(end - start)

        # obs (coordinates)
        str_dtype = np.dtypes.StringDType()
        obs_group = table_group.create_group("obs")
        obs_group.attrs["encoding-type"] = "dataframe"
        obs_group.attrs["encoding-version"] = "0.2.0"
        obs_group.attrs["_index"] = "instance_id"
        obs_group.attrs["column-order"] = [
            "y",
            "x",
            "region",
            "spatial_x",
            "spatial_y",
            "instance_key",
        ]

        y_values = np.repeat(np.arange(n_y, dtype=np.int32), n_x)
        x_values = np.tile(np.arange(n_x, dtype=np.int32), n_y)
        instance_ids = np.array([str(i) for i in range(n_rows)], dtype=str_dtype)
        spatial_x = x_values.astype(np.float64) * self.pixel_size_um
        spatial_y = y_values.astype(np.float64) * self.pixel_size_um

        a = obs_group.create_array("y", data=y_values)
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"
        a = obs_group.create_array("x", data=x_values)
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"
        a = obs_group.create_array("spatial_x", data=spatial_x)
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"
        a = obs_group.create_array("spatial_y", data=spatial_y)
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"
        a = obs_group.create_array("instance_id", data=instance_ids)
        a.attrs["encoding-type"] = "string-array"
        a.attrs["encoding-version"] = "0.2.0"
        a = obs_group.create_array("instance_key", data=instance_ids)
        a.attrs["encoding-type"] = "string-array"
        a.attrs["encoding-version"] = "0.2.0"

        # Region as categorical
        region_group = obs_group.create_group("region")
        region_group.attrs["encoding-type"] = "categorical"
        region_group.attrs["encoding-version"] = "0.2.0"
        region_group.attrs["ordered"] = False
        a = region_group.create_array(
            "categories", data=np.array([region_key], dtype=str_dtype)
        )
        a.attrs["encoding-type"] = "string-array"
        a.attrs["encoding-version"] = "0.2.0"
        a = region_group.create_array("codes", data=np.zeros(n_rows, dtype=np.int8))
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"

        # var (mass axis)
        var_group = table_group.create_group("var")
        var_group.attrs["encoding-type"] = "dataframe"
        var_group.attrs["encoding-version"] = "0.2.0"
        var_group.attrs["_index"] = "_index"
        var_group.attrs["column-order"] = ["mz"]

        mz_values = self._common_mass_axis
        if mz_values is None:
            raise RuntimeError("Common mass axis not initialized")
        mz_index = np.array([f"mz_{i}" for i in range(n_cols)], dtype=str_dtype)
        a = var_group.create_array("_index", data=mz_index)
        a.attrs["encoding-type"] = "string-array"
        a.attrs["encoding-version"] = "0.2.0"
        a = var_group.create_array("mz", data=mz_values)
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"

        # uns (metadata)
        uns_group = table_group.create_group("uns")
        uns_group.attrs["encoding-type"] = "dict"
        uns_group.attrs["encoding-version"] = "0.1.0"

        sd_attrs = uns_group.create_group("spatialdata_attrs")
        sd_attrs.attrs["encoding-type"] = "dict"
        sd_attrs.attrs["encoding-version"] = "0.1.0"
        a = sd_attrs.create_array("region", data=np.array(region_key, dtype=str_dtype))
        a.attrs["encoding-type"] = "string"
        a.attrs["encoding-version"] = "0.2.0"
        a = sd_attrs.create_array(
            "region_key", data=np.array("region", dtype=str_dtype)
        )
        a.attrs["encoding-type"] = "string"
        a.attrs["encoding-version"] = "0.2.0"
        a = sd_attrs.create_array(
            "instance_key",
            data=np.array("instance_key", dtype=str_dtype),
        )
        a.attrs["encoding-type"] = "string"
        a.attrs["encoding-version"] = "0.2.0"

        em_group = uns_group.create_group("essential_metadata")
        em_group.attrs["encoding-type"] = "dict"
        em_group.attrs["encoding-version"] = "0.1.0"
        a = em_group.create_array(
            "dimensions", data=np.asarray(np.array(self._dimensions))
        )
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"
        a = em_group.create_array(
            "mass_range",
            data=np.array([mz_values.min(), mz_values.max()]),
        )
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"
        a = em_group.create_array(
            "source_path",
            data=np.array(str(self.reader.data_path), dtype=str_dtype),
        )
        a.attrs["encoding-type"] = "string"
        a.attrs["encoding-version"] = "0.2.0"
        a = em_group.create_array(
            "spectrum_type",
            data=np.array("processed", dtype=str_dtype),
        )
        a.attrs["encoding-type"] = "string"
        a.attrs["encoding-version"] = "0.2.0"

        a = uns_group.create_array("average_spectrum", data=avg_spectrum)
        a.attrs["encoding-type"] = "array"
        a.attrs["encoding-version"] = "0.2.0"

        # Per-region mean spectra for multi-region datasets
        if avg_spectrum_per_region is not None:
            pr_group = uns_group.create_group("average_spectrum_per_region")
            pr_group.attrs["encoding-type"] = "dict"
            pr_group.attrs["encoding-version"] = "0.1.0"
            for region_key_str, region_avg in avg_spectrum_per_region.items():
                ra = pr_group.create_array(region_key_str, data=region_avg)
                ra.attrs["encoding-type"] = "array"
                ra.attrs["encoding-version"] = "0.2.0"

        # Create empty images and shapes groups (will be populated below)
        store.create_group("images")
        store.create_group("shapes")

        # Add TIC image and pixel shapes using SpatialData
        logger.info("  Adding TIC image and pixel shapes...")
        self._add_tic_image_and_shapes_to_store(
            tic_values, n_rows, n_x, n_y, slice_id, region_key
        )

        # Consolidate metadata after all elements are written
        logger.info("  Consolidating metadata...")
        with _suppress_upstream_warnings():
            zarr.consolidate_metadata(str(self.output_path))

    def _add_tic_image_and_shapes_to_store(
        self,
        tic_values: NDArray[np.float64],
        n_rows: int,
        n_x: int,
        n_y: int,
        slice_id: str,
        region_key: str,
    ) -> None:
        """Add TIC image and pixel shapes to the Zarr store using SpatialData.

        This method is called after the main table has been written to Zarr.
        It uses SpatialData's models to create properly formatted images and
        shapes, then writes them to the existing store.

        Args:
            tic_values: 2D array of TIC values (n_y, n_x).
            n_rows: Total number of pixels.
            n_x: Number of pixels in x dimension.
            n_y: Number of pixels in y dimension.
            slice_id: Identifier for this slice (e.g., "msi_dataset_z0").
            region_key: Region key for shapes (e.g., "msi_dataset_z0_pixels").
        """
        if not SPATIALDATA_AVAILABLE:
            logger.warning("SpatialData not available, skipping TIC image and shapes")
            return

        # === Create TIC Image ===
        y_size, x_size = tic_values.shape

        # Add channel dimension (c, y, x) as required by SpatialData
        tic_values_3d = tic_values.reshape(1, y_size, x_size)

        # When alignment exists, use raster indices + Affine transform
        if self._tic_to_image_matrix is not None:
            tic_xarray = xr.DataArray(
                tic_values_3d,
                dims=("c", "y", "x"),
                coords={
                    "c": [0],
                    "y": np.arange(y_size),
                    "x": np.arange(x_size),
                },
            )
            tic_transform = Affine(
                self._tic_to_image_matrix,
                input_axes=("x", "y"),
                output_axes=("x", "y"),
            )
        else:
            tic_xarray = xr.DataArray(
                tic_values_3d,
                dims=("c", "y", "x"),
                coords={
                    "c": [0],
                    "y": np.arange(y_size) * self.pixel_size_um,
                    "x": np.arange(x_size) * self.pixel_size_um,
                },
            )
            tic_transform = Identity()

        tic_image = Image2DModel.parse(
            tic_xarray,
            transformations={
                self.dataset_id: tic_transform,
                "global": tic_transform,
            },
        )

        # === Create Pixel Shapes ===
        gdf = self._create_streaming_pixel_shapes(n_rows, n_x, n_y)

        shape_transform = Identity()
        shapes = ShapesModel.parse(
            gdf,
            transformations={
                self.dataset_id: shape_transform,
                "global": shape_transform,
            },
        )

        # === Write to existing store using SpatialData ===
        # Open the store and add elements
        try:
            with _suppress_upstream_warnings():
                # Also suppress table-annotating-missing-shapes warning (shapes
                # are created and written immediately after this read).
                warnings.filterwarnings(
                    "ignore",
                    message="The table is annotating.*which is not present",
                    category=UserWarning,
                )
                sdata = SpatialData.read(str(self.output_path))

                # Add TIC image
                tic_name = f"{slice_id}_tic"
                sdata.images[tic_name] = tic_image

                # Add shapes
                sdata.shapes[region_key] = shapes

                # Write elements to disk
                sdata.write_element(tic_name, overwrite=True)
                sdata.write_element(region_key, overwrite=True)

                logger.info(
                    f"  Added TIC image '{tic_name}' ({x_size}x{y_size}) and "
                    f"{n_rows:,} pixel shapes"
                )

                # Add optical images if available
                self._add_optical_images_to_sdata(sdata, slice_id)

        except Exception as e:
            logger.warning(f"Failed to add TIC image and shapes: {e}")
            import traceback

            logger.debug(f"Traceback:\n{traceback.format_exc()}")

    def _create_streaming_pixel_shapes(
        self, n_rows: int, n_x: int, n_y: int
    ) -> "gpd.GeoDataFrame":
        """Create pixel shape geometries for the streaming converter.

        Args:
            n_rows: Total number of pixels.
            n_x: Number of pixels in x dimension.
            n_y: Number of pixels in y dimension.

        Returns:
            GeoDataFrame with pixel box geometries.
        """
        y_indices = np.repeat(np.arange(n_y), n_x)
        x_indices = np.tile(np.arange(n_x), n_y)

        from shapely import box as shapely_box_vectorized
        from shapely.geometry import box as shapely_box_single

        valid_indices: Optional[List[int]] = None

        if (
            self._alignment_result is not None
            and self._alignment_result.region_mappings
        ):
            # Use optical alignment - transform raster coords to image pixels.
            default_half_pixel = self._alignment_result.region_mappings[
                0
            ].get_half_pixel_size()

            valid_geometries: List[Any] = []
            valid_indices = []
            for i in range(n_rows):
                rx, ry = int(x_indices[i]), int(y_indices[i])
                img_coords = self._alignment_result.transform_point(rx, ry)

                if img_coords is not None:
                    ix, iy = img_coords
                    half_pixel = self._alignment_result.get_half_pixel_size(rx, ry)
                    if half_pixel is None:
                        half_pixel = default_half_pixel
                    half_x, half_y = half_pixel
                    valid_geometries.append(
                        shapely_box_single(
                            ix - half_x, iy - half_y, ix + half_x, iy + half_y
                        )
                    )
                    valid_indices.append(i)

            n_skipped = n_rows - len(valid_indices)
            if n_skipped > 0:
                logger.info(
                    f"Created {len(valid_geometries)} shapes using optical "
                    f"alignment (skipped {n_skipped} empty grid positions)"
                )

            geometries = valid_geometries
        else:
            # Physical micrometer coordinates
            half_pixel_um = self.pixel_size_um / 2
            spatial_x = x_indices * self.pixel_size_um
            spatial_y = y_indices * self.pixel_size_um

            geometries = shapely_box_vectorized(
                spatial_x - half_pixel_um,
                spatial_y - half_pixel_um,
                spatial_x + half_pixel_um,
                spatial_y + half_pixel_um,
            )

        # Create GeoDataFrame with string indices matching obs
        if valid_indices is not None:
            instance_ids = [str(i) for i in valid_indices]
        else:
            instance_ids = [str(i) for i in range(n_rows)]
        return gpd.GeoDataFrame(geometry=geometries, index=instance_ids)

    def _add_optical_images_to_sdata(self, sdata: "SpatialData", slice_id: str) -> None:
        """Add optical images directly to SpatialData store.

        The primary alignment image (from .mis <ImageFile>) gets an Identity
        transform and is loaded first. Other images get a Scale transform
        mapping their pixel coordinates to the primary image's space.

        Args:
            sdata: SpatialData object to add images to
            slice_id: Slice identifier for transformations
        """
        if not self._include_optical:
            return

        optical_paths = self.reader.get_optical_image_paths()
        if not optical_paths:
            logger.debug("No optical images found")
            return

        import tifffile

        logger.info(f"  Adding {len(optical_paths)} optical image(s)...")

        # Load primary image first so we know its dimensions for scaling others
        primary_paths = [p for p in optical_paths if self._is_primary_optical(p)]
        other_paths = [p for p in optical_paths if not self._is_primary_optical(p)]

        for tiff_path in primary_paths + other_paths:
            try:
                # Generate image name
                image_name = self._generate_optical_image_name(tiff_path)

                logger.info(f"    Loading: {tiff_path.name} as '{image_name}'")

                # Read TIFF
                with tifffile.TiffFile(tiff_path) as tif:
                    img_data = tif.pages[0].asarray()

                    # Convert to (c, y, x) format
                    if img_data.ndim == 2:
                        img_data = img_data[np.newaxis, :, :]
                    elif img_data.ndim == 3:
                        img_data = np.moveaxis(img_data, -1, 0)
                    else:
                        logger.warning(
                            f"Unexpected dimensions {img_data.ndim} for {tiff_path.name}"
                        )
                        continue

                    n_channels, y_size, x_size = img_data.shape

                    # Create xarray DataArray
                    optical_xarray = xr.DataArray(
                        img_data,
                        dims=("c", "y", "x"),
                        coords={
                            "c": list(range(n_channels)),
                            "y": np.arange(y_size),
                            "x": np.arange(x_size),
                        },
                    )

                    # Determine transform: primary gets Identity,
                    # others get Scale to align with primary
                    is_primary = self._is_primary_optical(tiff_path)
                    if is_primary:
                        transform = Identity()
                        self._primary_optical_dims = (x_size, y_size)
                        logger.info(f"    Primary alignment image: {x_size}x{y_size}")
                    elif self._primary_optical_dims is not None:
                        transform = self._compute_optical_scale_transform(
                            x_size, y_size
                        )
                    else:
                        transform = Identity()

                    optical_image = Image2DModel.parse(
                        optical_xarray,
                        transformations={
                            self.dataset_id: transform,
                            "global": transform,
                        },
                    )

                    # Add to SpatialData and write
                    sdata.images[image_name] = optical_image
                    sdata.write_element(image_name, overwrite=True)

                    logger.info(
                        f"    Added optical image ({x_size}x{y_size}, {n_channels}ch)"
                    )

            except Exception as e:
                logger.warning(f"Failed to load optical image {tiff_path.name}: {e}")
