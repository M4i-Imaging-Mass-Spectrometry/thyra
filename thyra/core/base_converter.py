import logging
from abc import ABC, abstractmethod
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pandas import DataFrame
from scipy import sparse
from tqdm import tqdm

from .base_reader import BaseMSIReader


class PixelSizeSource(Enum):
    """Enum to track how pixel size was determined."""

    DEFAULT = "default"  # Using default 1.0 (fallback)
    USER_PROVIDED = "manual"  # User explicitly provided via parameter
    AUTO_DETECTED = "automatic"  # Detected from metadata


class BaseMSIConverter(ABC):
    """Base class for MSI data converters with shared functionality.

    Implements common processing steps while allowing format-specific
    customization.
    """

    def __init__(
        self,
        reader: BaseMSIReader,
        output_path: Union[str, Path, PathLike[str]],
        dataset_id: str = "msi_dataset",
        pixel_size_um: float = 1.0,
        pixel_size_source: PixelSizeSource = PixelSizeSource.DEFAULT,
        compression_level: int = 5,
        handle_3d: bool = False,
        **kwargs: Any,
    ):
        """Initialize the MSI converter.

        Args:
            reader: MSI data reader instance
            output_path: Path for output file
            dataset_id: Identifier for the dataset
            pixel_size_um: Size of each pixel in micrometers
            pixel_size_source: How pixel size was determined
            compression_level: Compression level for output
            handle_3d: Whether to process as 3D data
            **kwargs: Additional keyword arguments
        """
        self.reader = reader
        self.output_path = Path(output_path)
        self.dataset_id = dataset_id
        self.pixel_size_um = pixel_size_um
        self.pixel_size_source = pixel_size_source
        self.compression_level = compression_level
        self.handle_3d = handle_3d
        self.options: Dict[str, Any] = kwargs
        self._common_mass_axis: Optional[NDArray[np.float64]] = None
        self._dimensions: Optional[Tuple[int, int, int]] = None
        self._metadata: Optional[dict[str, Any]] = None
        from ..config import DEFAULT_BUFFER_SIZE

        self._buffer_size = DEFAULT_BUFFER_SIZE

        # Essential metadata properties (loaded during initialization)
        self._coordinate_bounds: Optional[Tuple[float, float, float, float]] = None
        self._n_spectra: Optional[int] = None
        self._estimated_memory_gb: Optional[float] = None

    def convert(self) -> bool:
        """Template method defining the conversion workflow.

        Returns:
        --------
        bool: True if conversion was successful, False otherwise.
        """
        try:
            self._initialize_conversion()
            data_structures = self._create_data_structures()
            self._process_spectra(data_structures)
            self._finalize_data(data_structures)
            success = self._save_output(data_structures)

            return success
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            import traceback

            logging.error(f"Detailed traceback:\n{traceback.format_exc()}")
            return False
        finally:
            self.reader.close()

    def _initialize_conversion(self) -> None:
        """Initialize conversion by loading essential metadata first, then other data."""
        logging.info("Loading essential dataset information...")
        try:
            # Load essential metadata first (fast, single query for Bruker)
            essential = self.reader.get_essential_metadata()

            self._dimensions = essential.dimensions
            if any(d <= 0 for d in self._dimensions):
                raise ValueError(
                    f"Invalid dimensions: {self._dimensions}. All dimensions "
                    f"must be positive."
                )

            # Store essential metadata for use throughout conversion
            self._coordinate_bounds = essential.coordinate_bounds
            self._n_spectra = essential.n_spectra
            self._estimated_memory_gb = essential.estimated_memory_gb

            # Override pixel size only if using default value and metadata
            # is available
            if (
                self.pixel_size_source == PixelSizeSource.DEFAULT
                and essential.pixel_size
            ):
                old_size = self.pixel_size_um
                self.pixel_size_um = essential.pixel_size[0]
                self.pixel_size_source = PixelSizeSource.AUTO_DETECTED
                logging.info(
                    f"Auto-detected pixel size: {self.pixel_size_um} um "
                    f"(was default: {old_size} um)"
                )
            elif self.pixel_size_source == PixelSizeSource.USER_PROVIDED:
                logging.info(
                    f"Using user-specified pixel size: {self.pixel_size_um} um"
                )

            # Load mass axis separately (still expensive operation)
            self._common_mass_axis = self.reader.get_common_mass_axis()
            if len(self._common_mass_axis) == 0:
                raise ValueError(
                    "Common mass axis is empty. Cannot proceed with " "conversion."
                )

            # Only load comprehensive metadata if needed (lazy loading)
            self._metadata = None  # Will be loaded on demand

            logging.info(f"Dataset dimensions: {self._dimensions}")
            logging.info(f"Coordinate bounds: {self._coordinate_bounds}")
            logging.info(f"Total spectra: {self._n_spectra}")
            logging.info(f"Estimated memory: {self._estimated_memory_gb:.2f} GB")
            logging.info(f"Common mass axis length: {len(self._common_mass_axis)}")
        except Exception as e:
            logging.error(f"Error during initialization: {e}")
            raise

    @abstractmethod
    def _create_data_structures(self) -> Any:
        """Create format-specific data structures.

        Returns:
        --------
        Any: Format-specific data structures to be used in subsequent steps.
        """
        pass

    def _process_spectra(self, data_structures: Any) -> None:
        """Process all spectra from the reader and integrate into data structures.

        Parameters:
        -----------
        data_structures: Format-specific data containers created by
            _create_data_structures.
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized.")

        total_spectra = self._get_total_spectra_count()
        logging.info(
            f"Converting {total_spectra} spectra to "
            f"{self.__class__.__name__.replace('Converter', '')} format..."
        )

        setattr(self.reader, "_quiet_mode", True)

        # Process spectra with unified progress tracking
        with tqdm(
            total=total_spectra, desc="Converting spectra", unit="spectrum"
        ) as pbar:
            for coords, mzs, intensities in self.reader.iter_spectra(
                batch_size=self._buffer_size
            ):
                self._process_single_spectrum(data_structures, coords, mzs, intensities)
                pbar.update(1)

    def _get_total_spectra_count(self) -> int:
        """Get the total number of spectra for progress tracking.

        Uses cached essential metadata for efficient access.
        """
        # Use cached spectra count from essential metadata
        if self._n_spectra is not None:
            return self._n_spectra

        # Fallback: try reader-specific methods
        if hasattr(self.reader, "n_spectra"):
            return int(self.reader.n_spectra)

        # For ImzML readers, count coordinates
        if hasattr(self.reader, "parser") and self.reader.parser is not None:
            if hasattr(self.reader.parser, "coordinates"):
                return len(self.reader.parser.coordinates)

        # For Bruker readers, try frame count methods
        if hasattr(self.reader, "_get_frame_count"):
            return int(self.reader._get_frame_count())

        # Final fallback: calculate from dimensions
        if self._dimensions is not None:
            total_pixels = (
                self._dimensions[0] * self._dimensions[1] * self._dimensions[2]
            )
            logging.warning(
                f"Could not determine exact spectra count, estimating "
                f"{total_pixels} from dimensions"
            )
            return total_pixels

        # Should not reach here if initialization was successful
        raise ValueError(
            "Cannot determine spectra count - conversion not properly " "initialized"
        )

    def _process_single_spectrum(
        self,
        data_structures: Any,
        coords: Tuple[int, int, int],
        mzs: NDArray[np.float64],
        intensities: NDArray[np.float64],
    ) -> None:
        """Process a single spectrum.

        Args:
            data_structures: Format-specific data containers
            coords: (x, y, z) coordinates
            mzs: m/z values
            intensities: Intensity values
        """
        # Default implementation - to be overridden by subclasses if needed
        pass

    def _finalize_data(self, data_structures: Any) -> None:
        """Perform any final processing on the data structures before saving.

        Args:
            data_structures: Format-specific data containers
        """
        # Default implementation - to be overridden by subclasses if needed
        pass

    def _get_comprehensive_metadata(self) -> Dict[str, Any]:
        """Lazy load comprehensive metadata when needed."""
        if self._metadata is None:
            logging.info("Loading comprehensive metadata...")
            comprehensive = self.reader.get_comprehensive_metadata()
            self._metadata = comprehensive.raw_metadata
        return self._metadata

    @abstractmethod
    def _save_output(self, data_structures: Any) -> bool:
        """Save the processed data to the output format.

        Args:
            data_structures: Format-specific data containers

        Returns:
            True if saving was successful, False otherwise
        """
        pass

    def add_metadata(self, metadata: Any) -> None:
        """Add comprehensive metadata to the output.

        Base implementation provides common metadata structure.
        Subclasses should override to add format-specific metadata storage.

        Args:
            metadata: Any object that can store metadata
        """
        # Get comprehensive metadata for complete information
        comprehensive_metadata = self.reader.get_comprehensive_metadata()

        # Create structured metadata dict that subclasses can use
        self._structured_metadata = {
            # Conversion metadata
            "conversion_info": {
                "dataset_id": self.dataset_id,
                "pixel_size_um": self.pixel_size_um,
                "handle_3d": self.handle_3d,
                "compression_level": self.compression_level,
                "converter_class": self.__class__.__name__,
                "conversion_timestamp": pd.Timestamp.now().isoformat(),
            },
            # Essential metadata for quick access
            "essential_metadata": {
                "dimensions": comprehensive_metadata.essential.dimensions,
                "coordinate_bounds": (
                    comprehensive_metadata.essential.coordinate_bounds
                ),
                "mass_range": comprehensive_metadata.essential.mass_range,
                "pixel_size": comprehensive_metadata.essential.pixel_size,
                "n_spectra": comprehensive_metadata.essential.n_spectra,
                "estimated_memory_gb": (
                    comprehensive_metadata.essential.estimated_memory_gb
                ),
                "source_path": comprehensive_metadata.essential.source_path,
                "is_3d": comprehensive_metadata.essential.is_3d,
                "has_pixel_size": (comprehensive_metadata.essential.has_pixel_size),
            },
            # Format-specific metadata from source
            "format_specific_metadata": comprehensive_metadata.format_specific,
            "acquisition_parameters": (comprehensive_metadata.acquisition_params),
            "instrument_information": comprehensive_metadata.instrument_info,
            "raw_metadata": comprehensive_metadata.raw_metadata,
            # Processing statistics
            "processing_stats": {
                "total_grid_pixels": (
                    self._dimensions[0] * self._dimensions[1] * self._dimensions[2]
                    if self._dimensions
                    else 0
                ),
                "coordinate_bounds": self._coordinate_bounds,
                "estimated_memory_gb": self._estimated_memory_gb,
            },
        }

        # Subclasses should override to add this structured metadata to
        # their outputs
        logging.info(f"Base metadata structure prepared for {self.__class__.__name__}")

        # Default implementation does nothing - subclasses should override
        pass

    # --- Common Utility Methods ---

    def _create_sparse_matrix(self) -> sparse.lil_matrix:
        """Create sparse matrix for storing intensity values.

        Returns:
        --------
        sparse.lil_matrix: Sparse matrix for storing intensity values
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized.")
        n_x, n_y, n_z = self._dimensions
        n_pixels = n_x * n_y * n_z
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized.")
        n_masses = len(self._common_mass_axis)

        logging.info(
            f"Creating sparse matrix for {n_pixels} pixels and "
            f"{n_masses} mass values"
        )

        return sparse.lil_matrix((n_pixels, n_masses), dtype=np.float64)

    def _create_coordinates_dataframe(self) -> pd.DataFrame:
        """Create a DataFrame containing pixel coordinates.

        Returns:
        --------
        pd.DataFrame: DataFrame with pixel coordinates
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized.")
        n_x, n_y, n_z = self._dimensions

        coords = []
        for z in range(n_z):
            for y in range(n_y):
                for x in range(n_x):
                    pixel_idx = z * (n_y * n_x) + y * n_x + x
                    coords.append(
                        {
                            "z": z,
                            "y": y,
                            "x": x,
                            "pixel_id": str(
                                pixel_idx
                            ),  # Convert to string for compatibility
                        }
                    )

        coords_df: pd.DataFrame = pd.DataFrame(coords)
        coords_df.set_index("pixel_id", inplace=True)

        # Add spatial coordinates
        coords_df["spatial_x"] = coords_df["x"] * self.pixel_size_um
        coords_df["spatial_y"] = coords_df["y"] * self.pixel_size_um
        coords_df["spatial_z"] = coords_df["z"] * self.pixel_size_um

        return coords_df

    def _create_mass_dataframe(self) -> pd.DataFrame:
        """Create a DataFrame containing mass values.

        Returns:
        --------
        pd.DataFrame: DataFrame with mass values
        """
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized.")
        var_df: DataFrame = pd.DataFrame({"mz": self._common_mass_axis})
        # Convert to string index for compatibility
        var_df["mz_str"] = var_df["mz"].astype(str)
        var_df.set_index("mz_str", inplace=True)

        return var_df

    def _get_pixel_index(self, x: int, y: int, z: int) -> int:
        """Convert 3D coordinates to a flat array index.

        Args:
            x: X coordinate
            y: Y coordinate
            z: Z coordinate

        Returns:
            Flat index
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized.")
        n_x, n_y, _ = self._dimensions
        return z * (n_y * n_x) + y * n_x + x

    def _map_mass_to_indices(self, mzs: NDArray[np.float64]) -> NDArray[np.int_]:
        """Map m/z values to indices in the common mass axis with high accuracy.

        Args:
            mzs: Array of m/z values

        Returns:
            Array of indices in common mass axis
        """
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized.")

        if mzs.size == 0:
            return np.array([], dtype=int)

        # Use searchsorted for exact mapping
        indices = np.searchsorted(self._common_mass_axis, mzs)

        # Ensure indices are within bounds
        indices = np.clip(indices, 0, len(self._common_mass_axis) - 1)

        # For complete accuracy, validate the indices
        # Very small tolerance threshold for floating point differences
        max_diff = 1e-6
        mask = np.abs(self._common_mass_axis[indices] - mzs) <= max_diff

        return indices[mask]

    def _add_to_sparse_matrix(
        self,
        sparse_matrix: sparse.lil_matrix,
        pixel_idx: int,
        mz_indices: NDArray[np.int_],
        intensities: NDArray[np.float64],
    ) -> None:
        """Add intensity values to a sparse matrix efficiently.

        Args:
            sparse_matrix: Target sparse matrix
            pixel_idx: Flat pixel index
            mz_indices: Indices in common mass axis
            intensities: Intensity values
        """
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized.")

        if mz_indices.size == 0 or intensities.size == 0:
            logging.debug(
                f"Empty data for pixel {pixel_idx}: {mz_indices.size} "
                f"indices, "
                f"{intensities.size} intensities"
            )
            return

        n_masses = len(self._common_mass_axis)

        # Filter out invalid indices and zero intensities in a single pass
        valid_mask = (mz_indices < n_masses) & (intensities > 0)

        logging.info(
            f"Pixel {pixel_idx}: {len(mz_indices)} input indices, "
            f"{np.sum(valid_mask)} valid after filtering"
        )
        logging.info(
            f"  Index bounds check: "
            f"{np.sum(mz_indices < n_masses)}/{len(mz_indices)}"
        )
        logging.info(
            f"  Intensity > 0 check: " f"{np.sum(intensities > 0)}/{len(intensities)}"
        )

        if not np.any(valid_mask):
            logging.info(f"No valid data to store for pixel {pixel_idx}")
            return

        # Extract valid values
        valid_indices = mz_indices[valid_mask]
        valid_intensities = intensities[valid_mask]

        logging.info(
            f"Storing {len(valid_indices)} values for pixel {pixel_idx}, "
            f"intensity sum: {np.sum(valid_intensities):.2e}"
        )

        # Use bulk assignment for better performance
        sparse_matrix[pixel_idx, valid_indices] = valid_intensities

        logging.info(f"After storage - matrix nnz: {sparse_matrix.nnz}")
