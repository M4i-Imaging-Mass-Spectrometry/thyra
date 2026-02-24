# thyra/converters/spatialdata/base_spatialdata_converter.py

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from ...alignment import AreaAlignmentResult, TeachingPointAlignment
from ...core.base_converter import BaseMSIConverter, PixelSizeSource
from ...core.base_reader import BaseMSIReader
from ...metadata.types import ComprehensiveMetadata, EssentialMetadata
from ...resampling import ResamplingDecisionTree, ResamplingMethod
from ...resampling.types import ResamplingConfig

# Check SpatialData availability (defer imports to avoid issues)
SPATIALDATA_AVAILABLE = False
_import_error_msg = None
try:
    import geopandas as gpd
    import tifffile
    import xarray as xr
    from anndata import AnnData
    from shapely.geometry import box
    from spatialdata import SpatialData
    from spatialdata.models import Image2DModel, ShapesModel, TableModel
    from spatialdata.transformations import Affine, Identity

    SPATIALDATA_AVAILABLE = True
except (ImportError, NotImplementedError) as e:
    _import_error_msg = str(e)
    logging.warning(f"SpatialData dependencies not available: {e}")
    SPATIALDATA_AVAILABLE = False

    # Create dummy classes for registration
    AnnData = None
    SpatialData = None
    TableModel = None
    ShapesModel = None
    Image2DModel = None
    Affine = None
    Identity = None
    box = None
    gpd = None


class BaseSpatialDataConverter(BaseMSIConverter, ABC):
    """Base converter for MSI data to SpatialData format with shared functionality."""

    def __init__(
        self,
        reader: BaseMSIReader,
        output_path: Path,
        dataset_id: str = "msi_dataset",
        pixel_size_um: float = 1.0,
        pixel_size_source: PixelSizeSource = PixelSizeSource.DEFAULT,
        handle_3d: bool = False,
        pixel_size_detection_info: Optional[Dict[str, Any]] = None,
        resampling_config: Optional[Union[Dict[str, Any], ResamplingConfig]] = None,
        sparse_format: str = "csc",
        include_optical: bool = True,
        **kwargs: Any,
    ) -> None:
        """Initialize the base SpatialData converter.

        Args:
            reader: MSI data reader
            output_path: Path for output file
            dataset_id: Identifier for the dataset
            pixel_size_um: Size of each pixel in micrometers
            pixel_size_source: How pixel size was determined
            handle_3d: Whether to process as 3D data (True) or 2D slices
                (False)
            pixel_size_detection_info: Optional metadata about pixel size
                detection
            resampling_config: Optional resampling configuration dict
            sparse_format: Sparse matrix format ('csc' or 'csr', default: 'csc')
            include_optical: Whether to include optical images in output
                (default: True)
            **kwargs: Additional keyword arguments

        Raises:
            ImportError: If SpatialData dependencies are not available
            ValueError: If pixel_size_um is not positive or dataset_id is
                empty
        """
        # Check if SpatialData is available
        if not SPATIALDATA_AVAILABLE:
            error_msg = (
                f"SpatialData dependencies not available: "
                f"{_import_error_msg}. "
                f"Please install required packages or fix dependency "
                f"conflicts."
            )
            raise ImportError(error_msg)

        # Validate inputs
        if pixel_size_um <= 0:
            raise ValueError(f"pixel_size_um must be positive, got {pixel_size_um}")
        if not dataset_id.strip():
            raise ValueError("dataset_id cannot be empty")

        # Extract pixel_size_detection_info from kwargs if provided
        kwargs_filtered = dict(kwargs)
        if (
            pixel_size_detection_info is None
            and "pixel_size_detection_info" in kwargs_filtered
        ):
            pixel_size_detection_info = kwargs_filtered.pop("pixel_size_detection_info")

        super().__init__(
            reader,
            output_path,
            dataset_id=dataset_id,
            pixel_size_um=pixel_size_um,
            pixel_size_source=pixel_size_source,
            handle_3d=handle_3d,
            **kwargs_filtered,
        )

        self._non_empty_pixel_count: int = 0
        self._pixel_size_detection_info = pixel_size_detection_info
        self._resampling_config = resampling_config
        self._sparse_format = sparse_format.lower()
        self._include_optical = include_optical
        if self._sparse_format not in ("csc", "csr"):
            raise ValueError(
                f"sparse_format must be 'csc' or 'csr', got '{sparse_format}'"
            )

        # Set up resampling if enabled
        if self._resampling_config:
            self._setup_resampling()
            # Note: _build_resampled_mass_axis() will be called in _initialize_conversion()
            # after reader metadata is fully loaded

        # Cache for dense mass axis indices (to avoid repeated np.arange calls)
        self._cached_mass_axis_indices: Optional[NDArray[np.int_]] = None

        # Optical-MSI alignment (computed from FlexImaging Area definitions)
        self._alignment_result: Optional[AreaAlignmentResult] = None
        # Affine matrix mapping TIC raster indices to optical image pixels
        self._tic_to_image_matrix: Optional[NDArray[np.float64]] = None

        # Metadata caches (populated lazily during conversion)
        self._essential_metadata_cached: Optional[EssentialMetadata] = None
        self._comprehensive_metadata_cached: Optional[ComprehensiveMetadata] = None
        self._spectrum_metadata_cached: Optional[Dict[str, Any]] = None
        self._resampling_metadata_cached: Optional[Dict[str, Any]] = None

    def _setup_resampling(self) -> None:  # noqa: C901
        """Set up resampling configuration and strategy."""
        if not self._resampling_config:
            return

        # Access method - handle both dict and dataclass
        # Check for ResamplingConfig dataclass first (more specific)
        if isinstance(self._resampling_config, ResamplingConfig):
            method = self._resampling_config.method
            axis_type = self._resampling_config.axis_type
        elif isinstance(self._resampling_config, dict):
            method_raw = self._resampling_config.get("method", "auto")
            axis_type_raw = self._resampling_config.get("axis_type", "auto")
            # Convert string to enum if needed
            if isinstance(method_raw, str):
                method_map = {
                    "nearest_neighbor": ResamplingMethod.NEAREST_NEIGHBOR,
                    "tic_preserving": ResamplingMethod.TIC_PRESERVING,
                }
                method = method_map.get(method_raw, None)
            else:
                method = method_raw
            # Convert axis_type string to enum if needed
            if isinstance(axis_type_raw, str):
                from ...resampling.types import AxisType

                axis_type_map = {
                    "constant": AxisType.CONSTANT,
                    "linear_tof": AxisType.LINEAR_TOF,
                    "reflector_tof": AxisType.REFLECTOR_TOF,
                    "orbitrap": AxisType.ORBITRAP,
                    "fticr": AxisType.FTICR,
                }
                axis_type = axis_type_map.get(axis_type_raw, None)
            else:
                axis_type = axis_type_raw

        # If method is None or "auto", use DecisionTree to determine strategy
        if method is None:
            try:
                # Get metadata from reader for instrument detection
                metadata = self._get_reader_metadata_for_resampling()
                tree = ResamplingDecisionTree()
                detected_method = tree.select_strategy(metadata)
                logging.info(f"Auto-detected resampling method: {detected_method}")
                self._resampling_method = detected_method
            except NotImplementedError as e:
                logging.error(f"Auto-detection failed: {e}")
                logging.info("Falling back to nearest_neighbor for resampling")
                self._resampling_method = ResamplingMethod.NEAREST_NEIGHBOR
        else:
            # Use provided method directly (already an enum)
            self._resampling_method = method

        logging.info(f"Using resampling method: {self._resampling_method}")

        # Store axis_type override if provided (will be used in _build_resampled_mass_axis)
        self._manual_axis_type = axis_type

        # Store resampling parameters - handle both dict and dataclass
        if isinstance(self._resampling_config, ResamplingConfig):
            # ResamplingConfig dataclass
            self._target_bins = self._resampling_config.target_bins
            self._min_mz = self._resampling_config.min_mz
            self._max_mz = self._resampling_config.max_mz
            self._width_at_mz = self._resampling_config.mass_width_da
            self._reference_mz = self._resampling_config.reference_mz
        elif isinstance(self._resampling_config, dict):
            self._target_bins = self._resampling_config.get("target_bins", None)
            self._min_mz = self._resampling_config.get("min_mz")
            self._max_mz = self._resampling_config.get("max_mz")
            self._width_at_mz = self._resampling_config.get("width_at_mz")
            self._reference_mz = self._resampling_config.get("reference_mz", 1000.0)

    def _get_cached_metadata_for_resampling(self) -> Dict[str, Any]:
        """Get cached metadata for resampling decision tree to avoid multiple reader calls."""
        if self._resampling_metadata_cached is not None:
            return self._resampling_metadata_cached

        # If not cached yet, extract and cache it
        return self._get_reader_metadata_for_resampling()

    def _get_reader_metadata_for_resampling(self) -> Dict[str, Any]:
        """Extract metadata from reader for resampling decision tree."""
        metadata: Dict[str, Any] = {}

        # Extract different types of metadata
        self._extract_essential_metadata(metadata)
        self._extract_comprehensive_metadata(metadata)
        self._extract_spectrum_metadata(metadata)

        # Cache for later reuse
        self._resampling_metadata_cached = metadata
        return metadata

    def _extract_essential_metadata(self, metadata: Dict[str, Any]) -> None:
        """Extract essential metadata for resampling decisions."""
        try:
            # Use cached essential metadata if available
            if self._essential_metadata_cached is not None:
                essential = self._essential_metadata_cached
            else:
                essential = self.reader.get_essential_metadata()
                self._essential_metadata_cached = essential

            if hasattr(essential, "source_path"):
                metadata["source_path"] = str(essential.source_path)

            # Add essential metadata for resampling decisions
            metadata["essential_metadata"] = {
                "spectrum_type": getattr(essential, "spectrum_type", None),
                "dimensions": essential.dimensions,
                "mass_range": essential.mass_range,
                "source_path": str(essential.source_path),
                "total_peaks": getattr(essential, "total_peaks", None),
                "n_spectra": getattr(essential, "n_spectra", None),
            }
        except Exception as e:
            logging.debug(f"Could not extract essential metadata: {e}")

    def _extract_comprehensive_metadata(self, metadata: Dict[str, Any]) -> None:
        """Extract comprehensive metadata including Bruker GlobalMetadata."""
        try:
            # Use cached comprehensive metadata if available
            if self._comprehensive_metadata_cached is not None:
                comp_meta = self._comprehensive_metadata_cached
            else:
                comp_meta = self.reader.get_comprehensive_metadata()
                self._comprehensive_metadata_cached = comp_meta

            self._extract_bruker_metadata(metadata, comp_meta)
            self._extract_instrument_info(metadata, comp_meta)
        except Exception as e:
            logging.debug(f"Could not extract comprehensive metadata: {e}")

    def _extract_bruker_metadata(self, metadata: Dict[str, Any], comp_meta) -> None:
        """Extract Bruker GlobalMetadata from comprehensive metadata."""
        if (
            hasattr(comp_meta, "raw_metadata")
            and "global_metadata" in comp_meta.raw_metadata
        ):
            metadata["GlobalMetadata"] = comp_meta.raw_metadata["global_metadata"]
            logging.debug(
                f"Extracted Bruker GlobalMetadata with keys: "
                f"{list(metadata['GlobalMetadata'].keys())}"
            )

    def _extract_instrument_info(self, metadata: Dict[str, Any], comp_meta) -> None:
        """Extract instrument_info for fallback detection."""
        if hasattr(comp_meta, "instrument_info"):
            metadata["instrument_info"] = comp_meta.instrument_info
            logging.debug(f"Extracted instrument_info: {comp_meta.instrument_info}")

        # Extract format_specific for FlexImaging detection
        if hasattr(comp_meta, "format_specific"):
            metadata["format_specific"] = comp_meta.format_specific
            logging.debug(f"Extracted format_specific: {comp_meta.format_specific}")

        # Extract acquisition_params for additional detection
        if hasattr(comp_meta, "acquisition_params"):
            metadata["acquisition_params"] = comp_meta.acquisition_params
            logging.debug(
                f"Extracted acquisition_params: {comp_meta.acquisition_params}"
            )

    def _extract_spectrum_metadata(self, metadata: Dict[str, Any]) -> None:
        """Extract ImzML-specific spectrum metadata."""
        try:
            if hasattr(self.reader, "get_spectrum_metadata"):
                # Use cached spectrum metadata if available
                if self._spectrum_metadata_cached is not None:
                    spec_meta = self._spectrum_metadata_cached
                else:
                    spec_meta = self.reader.get_spectrum_metadata()
                    self._spectrum_metadata_cached = spec_meta

                if spec_meta:
                    metadata.update(spec_meta)
        except Exception as e:
            logging.debug(f"Could not extract spectrum metadata: {e}")

    def _add_metadata_to_uns(self, adata) -> None:
        """Add MSI metadata to AnnData .uns for preservation in SpatialData.

        This stores comprehensive metadata including:
        - Essential metadata (dimensions, mass range, source)
        - Format-specific metadata (FlexImaging areas, teaching points, etc.)
        - Acquisition parameters
        - Instrument information
        - Raw metadata (complete original metadata for future use)
        """
        try:
            comp_meta = self.reader.get_comprehensive_metadata()

            self._store_format_specific(adata, comp_meta)
            self._store_acquisition_params(adata, comp_meta)
            self._store_instrument_info(adata, comp_meta)
            self._store_essential_metadata(adata, comp_meta)
            self._store_raw_metadata(adata, comp_meta)

            logging.debug("Added MSI metadata to AnnData .uns")
        except Exception as e:
            logging.debug(f"Could not add metadata to .uns: {e}")

    def _store_format_specific(self, adata, comp_meta) -> None:
        """Store format-specific metadata (FlexImaging areas, teaching points)."""
        if hasattr(comp_meta, "format_specific") and comp_meta.format_specific:
            adata.uns["format_specific"] = comp_meta.format_specific

    def _store_acquisition_params(self, adata, comp_meta) -> None:
        """Store acquisition parameters."""
        if hasattr(comp_meta, "acquisition_params") and comp_meta.acquisition_params:
            adata.uns["acquisition_params"] = comp_meta.acquisition_params

    def _store_instrument_info(self, adata, comp_meta) -> None:
        """Store instrument information."""
        if hasattr(comp_meta, "instrument_info") and comp_meta.instrument_info:
            adata.uns["instrument_info"] = comp_meta.instrument_info

    def _store_essential_metadata(self, adata, comp_meta) -> None:
        """Store essential metadata (convert tuples to lists for Zarr)."""
        if not hasattr(comp_meta, "essential"):
            return

        essential = comp_meta.essential
        dims = essential.dimensions
        mrange = essential.mass_range
        adata.uns["essential_metadata"] = {
            "source_path": str(essential.source_path),
            "dimensions": list(dims) if dims else None,
            "mass_range": list(mrange) if mrange else None,
            "spectrum_type": getattr(essential, "spectrum_type", None),
        }

    def _store_raw_metadata(self, adata, comp_meta) -> None:
        """Store raw metadata (complete original data for future use)."""
        if hasattr(comp_meta, "raw_metadata") and comp_meta.raw_metadata:
            adata.uns["raw_metadata"] = self._serialize_for_zarr(comp_meta.raw_metadata)

    def _serialize_for_zarr(self, obj):
        """Recursively convert tuples to lists for Zarr serialization."""
        if isinstance(obj, dict):
            return {k: self._serialize_for_zarr(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._serialize_for_zarr(item) for item in obj]
        elif hasattr(obj, "__dict__"):
            # Convert dataclass/object to dict
            return self._serialize_for_zarr(vars(obj))
        else:
            return obj

    def _calculate_bins_from_width(
        self, min_mz: float, max_mz: float, axis_type
    ) -> int:
        """Calculate optimal number of bins from desired width at reference m/z.

        Args:
            min_mz: Minimum m/z of the mass range
            max_mz: Maximum m/z of the mass range
            axis_type: The axis type (determines physics-based spacing)

        Returns:
            Calculated number of bins
        """
        # Calculate bins based on axis type physics
        if hasattr(axis_type, "value"):
            axis_name = axis_type.value
        else:
            axis_name = str(axis_type).split(".")[-1].lower()

        if self._width_at_mz is None:
            # Use axis-type-specific defaults for optimal resolution
            if axis_name == "linear_tof":
                # LINEAR_TOF (FlexImaging): Use SCiLS-like default ~17 mDa at m/z 300
                # This matches typical MALDI-TOF resolution and gives ~30k bins
                width_at_mz = 0.017  # 17.0 mDa in Da
                reference_mz = 300.0
            else:
                # Default for other axis types: 5.0 mDa at m/z 1000
                width_at_mz = 0.005  # 5.0 mDa in Da
                reference_mz = 1000.0
        else:
            width_at_mz = self._width_at_mz
            reference_mz = self._reference_mz

        logging.info(
            f"Calculating bins for {width_at_mz*1000:.1f} mDa width at m/z {reference_mz:.1f}"
        )

        if axis_name == "reflector_tof":
            # REFLECTOR_TOF: constant relative resolution (width ∝ m/z)
            # relative_resolution = reference_mz / width_at_mz
            # For logarithmic spacing: bins ≈ ln(max_mz/min_mz) * (reference_mz / width_at_mz)
            relative_resolution = reference_mz / width_at_mz
            bins = int(np.log(max_mz / min_mz) * relative_resolution)

        elif axis_name == "linear_tof":
            # LINEAR_TOF: bin width = k * sqrt(m/z), where k = width_at_mz / sqrt(reference_mz)
            # Number of bins: n = (2/k) * (sqrt(max_mz) - sqrt(min_mz))
            # This matches SCiLS Lab's "Linear TOF" mass axis calculation
            k = width_at_mz / np.sqrt(reference_mz)
            bins = int((2.0 / k) * (np.sqrt(max_mz) - np.sqrt(min_mz)))

        elif axis_name == "orbitrap":
            # ORBITRAP: width ∝ m/z^1.5
            # For 1/sqrt spacing: bins ≈ (1/sqrt(min_mz) - 1/sqrt(max_mz)) *
            # (reference_mz^1.5 / width_at_mz)
            scaling_factor = (reference_mz**1.5) / width_at_mz
            bins = int((1 / np.sqrt(min_mz) - 1 / np.sqrt(max_mz)) * scaling_factor)

        else:
            # LINEAR/CONSTANT: uniform spacing
            # bins = (max_mz - min_mz) / width_at_mz
            bins = int((max_mz - min_mz) / width_at_mz)

        # Ensure minimum bin count
        bins = max(100, bins)

        logging.info(f"Calculated {bins} bins for {axis_name} axis type")
        return bins

    def _get_reference_params(self, axis_type) -> Tuple[float, float]:
        """Get reference width and m/z for the given axis type.

        Returns axis-type-specific defaults if user didn't specify.
        """
        if self._width_at_mz is not None:
            return self._width_at_mz, self._reference_mz

        # Get axis name for default selection
        if hasattr(axis_type, "value"):
            axis_name = axis_type.value
        else:
            axis_name = str(axis_type).split(".")[-1].lower()

        # Use axis-type-specific defaults
        if axis_name == "linear_tof":
            # LINEAR_TOF (FlexImaging): Use SCiLS-like default ~17 mDa at m/z 300
            return 0.017, 300.0
        # Default for other axis types: 5.0 mDa at m/z 1000
        return 0.005, 1000.0

    def _build_resampled_mass_axis(self) -> None:
        """Build resampled mass axis using physics-based generators."""
        from ...resampling.common_axis import CommonAxisBuilder

        # Use cached essential metadata to avoid reader calls
        # This should be called after _initialize_conversion() has loaded metadata
        if self._essential_metadata_cached is None:
            # Cache essential metadata for reuse
            self._essential_metadata_cached = self.reader.get_essential_metadata()

        mass_range = self._essential_metadata_cached.mass_range
        min_mz = mass_range[0] if self._min_mz is None else self._min_mz
        max_mz = mass_range[1] if self._max_mz is None else self._max_mz

        # Check if manual axis type override was provided
        if hasattr(self, "_manual_axis_type") and self._manual_axis_type is not None:
            axis_type = self._manual_axis_type
            logging.info(f"Using manually specified axis type: {axis_type}")
        else:
            # Get metadata for axis type selection (minimize reader calls)
            metadata = self._get_cached_metadata_for_resampling()
            tree = ResamplingDecisionTree()
            axis_type = tree.select_axis_type(metadata)
            logging.info(f"Auto-detected axis type: {axis_type}")

        # Calculate bins based on width if specified OR if no bins were specified (default)
        if self._width_at_mz is not None or self._target_bins is None:
            # Either user specified width OR using default (calculate from 5mDa@1000)
            target_bins = self._calculate_bins_from_width(min_mz, max_mz, axis_type)
        else:
            # User explicitly specified bin count
            target_bins = self._target_bins

        logging.info(
            f"Building resampled mass axis: {min_mz:.2f} - {max_mz:.2f} m/z, "
            f"{target_bins} bins"
        )

        # Build the physics-based axis
        builder = CommonAxisBuilder()

        # Determine reference parameters for physics generators
        reference_width, reference_mz = self._get_reference_params(axis_type)

        if hasattr(axis_type, "value") and axis_type.value != "constant":
            # Use physics-based generator with reference parameters
            mass_axis = builder.build_physics_axis(
                min_mz=min_mz,
                max_mz=max_mz,
                num_bins=target_bins,
                axis_type=axis_type,
                reference_mz=reference_mz,
                reference_width=reference_width,
            )
            logging.info(
                f"Built physics-based {axis_type} mass axis with "
                f"{len(mass_axis.mz_values)} points"
            )
        else:
            # Fall back to uniform axis
            mass_axis = builder.build_uniform_axis(min_mz, max_mz, target_bins)
            logging.info(
                f"Built uniform mass axis with " f"{len(mass_axis.mz_values)} points"
            )

        # Override the parent's common mass axis
        self._common_mass_axis = mass_axis.mz_values.astype(np.float64)
        assert self._common_mass_axis is not None

        # Cache the mass axis indices array to avoid repeated np.arange() calls
        self._cached_mass_axis_indices = np.arange(
            len(self._common_mass_axis), dtype=np.int_
        )

        # Calculate bin sizes for informative logging
        bin_widths = np.diff(self._common_mass_axis)
        min_bin_size = np.min(bin_widths) * 1000  # Convert to mDa
        max_bin_size = np.max(bin_widths) * 1000  # Convert to mDa

        logging.info(
            f"Resampled mass axis created: {len(self._common_mass_axis)} bins, "
            f"range {self._common_mass_axis[0]:.2f}-{self._common_mass_axis[-1]:.2f} m/z, "
            f"bin sizes {min_bin_size:.2f}-{max_bin_size:.2f} mDa ({axis_type})"
        )

    def _initialize_conversion(self) -> None:
        """Override parent initialization to preserve resampled mass axis."""
        logging.info("Loading essential dataset information...")
        try:
            # Load essential metadata first (fast, single query for Bruker)
            essential = self.reader.get_essential_metadata()
            # Cache for reuse during resampling setup
            self._essential_metadata_cached = essential

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

            # Override pixel size only if using default and metadata is available
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

            # Handle mass axis setup
            config_status = "SET" if self._resampling_config else "NOT SET"
            logging.info(f"Mass axis mode: resampling_config={config_status}")
            if self._resampling_config:
                # Build resampled mass axis now that reader metadata is loaded
                logging.info(
                    "Building RESAMPLED mass axis (resampling enabled) - "
                    "will NOT iterate through all spectra"
                )
                self._build_resampled_mass_axis()
                assert self._common_mass_axis is not None
                logging.info(
                    f"Built resampled mass axis with "
                    f"{len(self._common_mass_axis)} bins"
                )
            else:
                # No resampling - load raw mass axis as usual
                if self.reader.has_shared_mass_axis:
                    logging.info(
                        "Loading RAW mass axis (no resampling) - "
                        "continuous mode, reading m/z from first spectrum only"
                    )
                else:
                    logging.warning(
                        "Building RAW mass axis (no resampling) - "
                        "processed mode, iterating ALL spectra to collect unique m/z values. "
                        "This is slow for large datasets!"
                    )
                self._common_mass_axis = self.reader.get_common_mass_axis()
                if len(self._common_mass_axis) == 0:
                    raise ValueError(
                        "Common mass axis is empty. Cannot proceed with " "conversion."
                    )
                logging.info(
                    f"Using raw mass axis with "
                    f"{len(self._common_mass_axis)} unique m/z values"
                )

            # Only load comprehensive metadata if needed (lazy loading)
            self._metadata = None  # Will be loaded on demand

            # Compute optical alignment for FlexImaging data
            self._compute_optical_alignment()
            self._build_tic_to_image_affine()

            logging.info(f"Dataset dimensions: {self._dimensions}")
            logging.info(f"Coordinate bounds: {self._coordinate_bounds}")
            logging.info(f"Total spectra: {self._n_spectra}")
            logging.info(f"Estimated memory: {self._estimated_memory_gb:.2f} GB")
            assert self._common_mass_axis is not None
            logging.info(f"Common mass axis length: {len(self._common_mass_axis)}")
        except Exception as e:
            logging.error(f"Error during initialization: {e}")
            raise

    def _map_mass_to_indices(self, mzs: NDArray[np.float64]) -> NDArray[np.int_]:
        """Override mass mapping to handle resampling with interpolation."""
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized.")

        if mzs.size == 0:
            return np.array([], dtype=int)

        # If resampling is enabled, we need to interpolate instead of exact
        # matching
        if self._resampling_config:
            return self._resample_spectrum_to_indices(mzs)
        else:
            # Use parent's exact matching for non-resampled data
            return super()._map_mass_to_indices(mzs)

    def _resample_spectrum_to_indices(
        self, mzs: NDArray[np.float64]
    ) -> NDArray[np.int_]:
        """Map spectrum m/z values to resampled mass axis indices using interpolation."""
        # For resampled data, we want to return ALL indices in the resampled
        # axis. The actual resampling/interpolation will be handled in the
        # processing
        assert self._common_mass_axis is not None
        return np.arange(len(self._common_mass_axis), dtype=np.int_)

    def _process_single_spectrum(
        self,
        data_structures: Any,
        coords: Tuple[int, int, int],
        mzs: NDArray[np.float64],
        intensities: NDArray[np.float64],
    ) -> None:
        """Override spectrum processing to handle resampling with sparse optimization."""
        # Only log detailed per-spectrum info at DEBUG level
        logging.debug(
            f"Processing spectrum at {coords}: {len(mzs)} peaks, "
            f"intensity sum: {np.sum(intensities):.2e}"
        )

        if self._resampling_config:
            # OPTIMIZATION: Use sparse resampling directly for nearest_neighbor
            if (
                hasattr(self, "_resampling_method")
                and self._resampling_method == ResamplingMethod.NEAREST_NEIGHBOR
            ):
                # Sparse path - much faster, no dense array allocation
                mz_indices, resampled_intensities = self._nearest_neighbor_resample(
                    mzs, intensities
                )
                logging.debug(
                    f"Resampled (sparse): {len(resampled_intensities)} non-zero bins, "
                    f"sum: {np.sum(resampled_intensities):.2e}"
                )
            else:
                # Dense path for other methods (TIC-preserving, etc.)
                resampled_intensities = self._resample_spectrum(mzs, intensities)
                # Use cached indices instead of creating new array every time
                assert self._cached_mass_axis_indices is not None
                mz_indices = self._cached_mass_axis_indices
                logging.debug(
                    f"Resampled: {len(resampled_intensities)} values, "
                    f"sum: {np.sum(resampled_intensities):.2e}"
                )

            # Call the specific converter's processing with resampled data
            self._process_resampled_spectrum(
                data_structures, coords, mz_indices, resampled_intensities
            )
        else:
            # Use standard processing for non-resampled data
            mz_indices = self._map_mass_to_indices(mzs)
            logging.debug(
                f"Mapped to {len(mz_indices)} indices, "
                f"intensity sum: {np.sum(intensities):.2e}"
            )
            self._process_resampled_spectrum(
                data_structures, coords, mz_indices, intensities
            )

    def _resample_spectrum(
        self, mzs: NDArray[np.float64], intensities: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Resample a single spectrum onto the common mass axis using TIC-preserving method.

        Note: Nearest neighbor resampling is handled directly in _process_single_spectrum
        for performance optimization (returns sparse data).
        """
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized")

        # This method is only called for TIC-preserving resampling (dense path)
        return self._tic_preserving_resample(mzs, intensities)

    def _nearest_neighbor_resample(
        self, mzs: NDArray[np.float64], intensities: NDArray[np.float64]
    ) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
        """Resample spectrum using nearest neighbor interpolation.

        Maps each m/z value to its nearest bin in the common mass axis and
        accumulates intensities. Returns only non-zero bins for efficiency.

        Args:
            mzs: Original m/z values from spectrum
            intensities: Corresponding intensity values

        Returns:
            Tuple of (bin_indices, accumulated_intensities) containing only non-zero bins
        """
        if mzs.size == 0:
            return np.array([], dtype=np.int_), np.array([], dtype=np.float64)

        # Ensure common mass axis is initialized
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized")
        # Find insertion points using vectorized binary search
        indices = np.searchsorted(self._common_mass_axis, mzs)

        # OPTIMIZED: Handle boundary and nearest neighbor in one pass
        # Clip to valid range
        indices_clipped = np.clip(indices, 0, len(self._common_mass_axis) - 1)

        # For non-boundary points, check if left is closer
        # Only check where we're not at the left edge
        check_left = indices > 0
        if np.any(check_left):
            # Get distances only for points that need checking
            mz_values = self._common_mass_axis[indices_clipped[check_left]]
            mz_values_left = self._common_mass_axis[indices_clipped[check_left] - 1]
            mz_query = mzs[check_left]

            # Use left if it's closer
            use_left = np.abs(mz_values_left - mz_query) < np.abs(mz_values - mz_query)
            indices_clipped[check_left] = np.where(
                use_left, indices_clipped[check_left] - 1, indices_clipped[check_left]
            )

        # OPTIMIZATION: Use pandas-style groupby or bincount for accumulation
        # bincount without minlength is fastest
        accumulated = np.bincount(indices_clipped, weights=intensities)

        # Extract non-zero bins for sparse representation
        nonzero_mask = accumulated > 0
        nonzero_indices = np.where(nonzero_mask)[0].astype(np.int_)
        nonzero_values = accumulated[nonzero_mask]

        return nonzero_indices, nonzero_values.astype(np.float64)

    def _tic_preserving_resample(
        self, mzs: NDArray[np.float64], intensities: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Resample using TIC-preserving linear interpolation - optimized."""
        assert self._common_mass_axis is not None
        if mzs.size == 0:
            return np.zeros(len(self._common_mass_axis))

        # OPTIMIZED: Check if already sorted to avoid unnecessary sorting
        if np.all(mzs[:-1] <= mzs[1:]):
            # Already sorted - use directly
            mzs_sorted = mzs
            intensities_sorted = intensities
        else:
            # Need to sort for interpolation
            sort_indices = np.argsort(mzs)
            mzs_sorted = mzs[sort_indices]
            intensities_sorted = intensities[sort_indices]

        # Interpolate onto the common mass axis (np.interp is highly optimized)
        assert self._common_mass_axis is not None
        resampled = np.interp(
            self._common_mass_axis,
            mzs_sorted,
            intensities_sorted,
            left=0,
            right=0,
        )

        return resampled

    def _process_resampled_spectrum(
        self,
        data_structures: Any,
        coords: Tuple[int, int, int],
        mz_indices: NDArray[np.int_],
        intensities: NDArray[np.float64],
    ) -> None:
        """Process a spectrum with resampled intensities - to be overridden by subclasses."""
        # This method should be overridden by specific converters (2D/3D)
        pass

    def _create_sparse_matrix(self) -> Dict[str, Any]:
        """Create COO arrays for storing intensity values.

        Returns:
            Dictionary containing pre-allocated COO arrays and metadata

        Raises:
            ValueError: If dimensions or common mass axis are not initialized
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized")

        n_x, n_y, n_z = self._dimensions
        n_pixels = n_x * n_y * n_z
        n_masses = len(self._common_mass_axis)

        # Get exact number of peaks from cached metadata (no iteration needed!)
        assert self._essential_metadata_cached is not None
        exact_nnz = self._essential_metadata_cached.total_peaks

        logging.info(
            f"Pre-allocating COO arrays: {n_pixels:,} pixels x {n_masses:,} m/z bins"
        )
        logging.info(f"Exact non-zero values from metadata: {exact_nnz:,}")

        # Pre-allocate arrays with exact size (uint32 for indices, all values are positive)
        return {
            "rows": np.empty(exact_nnz, dtype=np.uint32),
            "cols": np.empty(exact_nnz, dtype=np.uint32),
            "data": np.empty(exact_nnz, dtype=np.float64),
            "current_idx": 0,
            "n_rows": n_pixels,
            "n_cols": n_masses,
        }

    def _create_coordinates_dataframe(self) -> pd.DataFrame:
        """Create coordinates dataframe with pixel positions.

        Returns:
            DataFrame with pixel coordinates

        Raises:
            ValueError: If dimensions are not initialized
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")

        n_x, n_y, n_z = self._dimensions

        # Pre-allocate arrays for better performance
        coords_data = []

        pixel_idx = 0
        for z in range(n_z):
            for y in range(n_y):
                for x in range(n_x):
                    coords_data.append(
                        {
                            "x": x,
                            "y": y,
                            "z": z if n_z > 1 else 0,
                            "instance_id": str(pixel_idx),
                            "region": f"{self.dataset_id}_pixels",
                            "spatial_x": x * self.pixel_size_um,
                            "spatial_y": y * self.pixel_size_um,
                            "spatial_z": (z * self.pixel_size_um if n_z > 1 else 0.0),
                        }
                    )
                    pixel_idx += 1

        coords_df = pd.DataFrame(coords_data)
        coords_df.set_index("instance_id", inplace=True)
        return coords_df

    def _create_mass_dataframe(self) -> pd.DataFrame:
        """Create m/z dataframe for variable metadata.

        Returns:
            DataFrame with m/z values

        Raises:
            ValueError: If common mass axis is not initialized
        """
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized")

        return pd.DataFrame(
            {"mz": self._common_mass_axis},
            index=[f"mz_{i}" for i in range(len(self._common_mass_axis))],
        )

    def _get_pixel_index(self, x: int, y: int, z: int) -> int:
        """Calculate linear pixel index from 3D coordinates.

        Args:
            x: X coordinate
            y: Y coordinate
            z: Z coordinate

        Returns:
            Linear pixel index

        Raises:
            ValueError: If dimensions are not initialized
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")

        n_x, n_y, _ = self._dimensions
        return z * (n_x * n_y) + y * n_x + x

    def _add_to_sparse_matrix(
        self,
        coo_arrays: Dict[str, Any],
        pixel_idx: int,
        mz_indices: NDArray[np.int_],
        intensities: NDArray[np.float64],
    ) -> None:
        """Add intensity data to COO arrays efficiently.

        Args:
            coo_arrays: Dictionary with pre-allocated COO arrays
            pixel_idx: Linear pixel index
            mz_indices: Indices for mass values
            intensities: Intensity values to add
        """
        # Filter out zero intensities to maintain sparsity
        nonzero_mask = intensities != 0.0
        if not np.any(nonzero_mask):
            return

        valid_mz_indices = mz_indices[nonzero_mask]
        valid_intensities = intensities[nonzero_mask]
        n = len(valid_intensities)

        # Check if we need to resize arrays
        current_idx = coo_arrays["current_idx"]
        end_idx = current_idx + n
        current_size = len(coo_arrays["rows"])

        if end_idx > current_size:
            # Resize arrays by 50% to accommodate more data
            new_size = int(current_size * 1.5)
            if new_size < end_idx:
                new_size = end_idx + int(current_size * 0.1)  # Ensure it fits

            logging.info(
                f"Resizing COO arrays from {current_size:,} to {new_size:,} elements"
            )

            coo_arrays["rows"] = np.resize(coo_arrays["rows"], new_size)
            coo_arrays["cols"] = np.resize(coo_arrays["cols"], new_size)
            coo_arrays["data"] = np.resize(coo_arrays["data"], new_size)

        # Direct array assignment (vectorized, very fast)
        coo_arrays["rows"][current_idx:end_idx] = pixel_idx
        coo_arrays["cols"][current_idx:end_idx] = valid_mz_indices
        coo_arrays["data"][current_idx:end_idx] = valid_intensities

        coo_arrays["current_idx"] = end_idx

    def _create_pixel_shapes(
        self, adata: AnnData, is_3d: bool = False
    ) -> "ShapesModel":
        """Create geometric shapes for pixels with proper transformations.

        When optical alignment is available (FlexImaging with Area definitions),
        shapes are created in optical image pixel coordinates for proper overlay.
        Otherwise, shapes use physical (micrometer) coordinates.

        Args:
            adata: AnnData object containing coordinates
            is_3d: Whether to handle as 3D data

        Returns:
            SpatialData shapes model

        Raises:
            ImportError: If required SpatialData dependencies are not available
        """
        if not SPATIALDATA_AVAILABLE:
            raise ImportError("SpatialData dependencies not available")

        geometries = []

        # Track valid indices (for alignment mode where we skip empty positions)
        valid_indices: Optional[List[int]] = None

        if self._alignment_result is not None:
            # Use optical alignment - transform raster coords to image pixels
            # Only create shapes for positions that have actual spectra (in pos_to_region)
            raster_x: NDArray[np.int_] = adata.obs["x"].values
            raster_y: NDArray[np.int_] = adata.obs["y"].values

            # Default half-pixel size (fallback if region lookup fails)
            default_half_pixel = 10.0
            if self._alignment_result.region_mappings:
                # Use first region as default
                rm = self._alignment_result.region_mappings[0]
                default_half_pixel = rm.get_half_pixel_size()

            valid_indices = []
            for i in range(len(adata)):
                rx, ry = int(raster_x[i]), int(raster_y[i])
                img_coords = self._alignment_result.transform_point(rx, ry)

                if img_coords is not None:
                    ix, iy = img_coords
                    # Get region-specific half-pixel size (square pixels)
                    half_pixel = self._alignment_result.get_half_pixel_size(rx, ry)
                    if half_pixel is None:
                        half_pixel = default_half_pixel

                    pixel_box = box(
                        ix - half_pixel,
                        iy - half_pixel,
                        ix + half_pixel,
                        iy + half_pixel,
                    )
                    geometries.append(pixel_box)
                    valid_indices.append(i)
                # Skip positions without spectra - empty grid cells

            n_skipped = len(adata) - len(valid_indices)
            if n_skipped > 0:
                logging.info(
                    f"Created {len(geometries)} shapes using optical alignment "
                    f"(skipped {n_skipped} empty grid positions)"
                )
            else:
                logging.info(
                    f"Created {len(geometries)} shapes using optical alignment "
                    f"(image pixel coordinates)"
                )
        else:
            # Standard physical coordinates (micrometers)
            x_coords: NDArray[np.float64] = adata.obs["spatial_x"].values
            y_coords: NDArray[np.float64] = adata.obs["spatial_y"].values
            half_pixel = self.pixel_size_um / 2

            for i in range(len(adata)):
                x, y = x_coords[i], y_coords[i]
                pixel_box = box(
                    x - half_pixel, y - half_pixel, x + half_pixel, y + half_pixel
                )
                geometries.append(pixel_box)

        # Create GeoDataFrame with appropriate index
        if valid_indices is not None:
            # Use filtered indices for alignment mode
            filtered_index = adata.obs.index[valid_indices]
            gdf = gpd.GeoDataFrame(geometry=geometries, index=filtered_index)
        else:
            gdf = gpd.GeoDataFrame(geometry=geometries, index=adata.obs.index)

        # Set up transform
        transform = Identity()
        transformations = {self.dataset_id: transform, "global": transform}

        # Parse shapes
        shapes = ShapesModel.parse(gdf, transformations=transformations)
        return shapes

    def _compute_optical_alignment(self) -> None:
        """Compute optical-MSI alignment from reader metadata.

        For FlexImaging data with Area definitions, this computes the
        transformation that maps MSI raster coordinates to optical image
        pixel coordinates.
        """
        # Check if reader has FlexImaging-specific metadata
        if not hasattr(self.reader, "mis_metadata"):
            logging.debug("Reader does not have mis_metadata, skipping alignment")
            return

        mis_metadata = getattr(self.reader, "mis_metadata", {})
        areas = mis_metadata.get("areas", [])

        if not areas:
            logging.debug("No Area definitions found, skipping alignment")
            return

        # Get required data for alignment
        positions = getattr(self.reader, "_positions", [])
        header = getattr(self.reader, "_header", {})

        if not positions:
            logging.warning("No position data available for alignment")
            return

        first_raster_x = header.get("first_raster_x", 0)
        first_raster_y = header.get("first_raster_y", 0)

        # Compute area-based alignment
        try:
            aligner = TeachingPointAlignment()
            self._alignment_result = aligner.compute_area_alignment(
                areas=areas,
                poslog_positions=positions,
                first_raster_x=first_raster_x,
                first_raster_y=first_raster_y,
            )
            logging.info(
                f"Computed optical alignment with "
                f"{len(self._alignment_result.region_mappings)} region mappings"
            )
        except Exception as e:
            logging.warning(f"Failed to compute optical alignment: {e}")
            self._alignment_result = None

    def _build_tic_to_image_affine(self) -> None:
        """Build affine matrix mapping TIC raster-index coords to image pixels.

        When optical alignment is available, this creates a 3x3 affine matrix
        that transforms TIC image coordinates (integer raster indices) into
        optical image pixel coordinates, so the TIC overlays correctly on the
        optical image in SpatialData.

        The matrix encodes: image_pixel = scale * raster_index + offset
        where offset places the first pixel center at image_min + half_pixel.
        """
        if self._alignment_result is None:
            return
        if not self._alignment_result.region_mappings:
            return

        rm = self._alignment_result.region_mappings[0]
        scale = rm._get_pixel_scale()
        half = scale / 2.0

        # Translation: raster index (0,0) maps to (image_min + half, image_min + half)
        tx = rm.image_min_x + half
        ty = rm.image_min_y + half

        # 3x3 affine: [[scale, 0, tx], [0, scale, ty], [0, 0, 1]]
        self._tic_to_image_matrix = np.array(
            [[scale, 0, tx], [0, scale, ty], [0, 0, 1]], dtype=np.float64
        )
        logging.info(
            f"Built TIC-to-image affine: scale={scale:.2f}, "
            f"offset=({tx:.1f}, {ty:.1f})"
        )

    def _add_optical_images(self, data_structures: Dict[str, Any]) -> None:
        """Load and add optical images from the reader to data structures.

        Finds optical TIFF files associated with the MSI data and adds them
        as image layers in the SpatialData output.

        Args:
            data_structures: Data structures dict to add images to
        """
        if not self._include_optical:
            return

        optical_paths = self.reader.get_optical_image_paths()
        if not optical_paths:
            logging.debug("No optical images found")
            return

        logging.info(f"Found {len(optical_paths)} optical image(s)")

        for tiff_path in optical_paths:
            try:
                self._load_single_optical_image(tiff_path, data_structures)
            except Exception as e:
                logging.warning(f"Failed to load optical image {tiff_path.name}: {e}")

    def _load_single_optical_image(
        self, tiff_path: Path, data_structures: Dict[str, Any]
    ) -> None:
        """Load a single optical TIFF and add it to data structures.

        Args:
            tiff_path: Path to the TIFF file
            data_structures: Data structures dict to add the image to
        """
        # Generate a clean name for the image layer
        image_name = self._generate_optical_image_name(tiff_path)

        logging.info(f"Loading optical image: {tiff_path.name} as '{image_name}'")

        # Read TIFF using tifffile (handles large files efficiently)
        with tifffile.TiffFile(tiff_path) as tif:
            # Read the first page/frame
            img_data = tif.pages[0].asarray()

            # Get image dimensions
            if img_data.ndim == 2:
                # Grayscale: (y, x) -> (c, y, x)
                img_data = img_data[np.newaxis, :, :]
                n_channels = 1
            elif img_data.ndim == 3:
                # RGB/RGBA: (y, x, c) -> (c, y, x)
                img_data = np.moveaxis(img_data, -1, 0)
                n_channels = img_data.shape[0]
            else:
                logging.warning(
                    f"Unexpected image dimensions {img_data.ndim} for {tiff_path.name}"
                )
                return

            y_size, x_size = img_data.shape[1], img_data.shape[2]

            # Create xarray DataArray
            optical_image = xr.DataArray(
                img_data,
                dims=("c", "y", "x"),
                coords={
                    "c": np.arange(n_channels),
                    "y": np.arange(y_size),
                    "x": np.arange(x_size),
                },
                attrs={
                    "source_file": tiff_path.name,
                    "original_path": str(tiff_path),
                },
            )

            # Create Image2DModel
            transform = Identity()
            data_structures["images"][image_name] = Image2DModel.parse(
                optical_image,
                transformations={
                    self.dataset_id: transform,
                    "global": transform,
                },
            )

            logging.info(
                f"Added optical image '{image_name}': {x_size}x{y_size} "
                f"({n_channels} channel{'s' if n_channels > 1 else ''})"
            )

    def _generate_optical_image_name(self, tiff_path: Path) -> str:
        """Generate a clean name for an optical image layer.

        Args:
            tiff_path: Path to the TIFF file

        Returns:
            Clean name for the image layer (e.g., 'optical_0000', 'optical_deriv')
        """
        stem = tiff_path.stem.lower()

        # Extract meaningful suffix from filename
        if "_0000" in stem:
            suffix = "highres"
        elif "_0001" in stem:
            suffix = "derived"
        elif "deriv" in stem:
            suffix = "overview"
        else:
            # Use stem with special chars replaced
            suffix = stem.replace(" ", "_").replace("-", "_")
            # Truncate if too long
            if len(suffix) > 30:
                suffix = suffix[:30]

        return f"{self.dataset_id}_optical_{suffix}"

    def _save_output(self, data_structures: Dict[str, Any]) -> bool:
        """Save the data to SpatialData format.

        Args:
            data_structures: Data structures to save

        Returns:
            True if saving was successful, False otherwise
        """
        if not SPATIALDATA_AVAILABLE:
            raise ImportError("SpatialData dependencies not available")

        try:
            # Create SpatialData object with images included
            sdata = SpatialData(
                tables=data_structures["tables"],
                shapes=data_structures["shapes"],
                images=data_structures["images"],
            )

            # Add metadata
            self.add_metadata(sdata)

            # Write to disk
            sdata.write(str(self.output_path))
            logging.info(f"Successfully saved SpatialData to {self.output_path}")
            return True
        except Exception as e:
            logging.error(f"Error saving SpatialData: {e}")
            import traceback

            logging.debug(f"Detailed traceback:\n{traceback.format_exc()}")
            return False

    def add_metadata(self, metadata: "SpatialData") -> None:
        """Add comprehensive metadata to the SpatialData object.

        Args:
            metadata: SpatialData object to add metadata to
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")

        # Call parent to prepare structured metadata
        super().add_metadata(metadata)

        # Get comprehensive metadata object for detailed access
        comprehensive_metadata_obj = self.reader.get_comprehensive_metadata()

        # Setup attributes and add pixel size metadata
        self._setup_spatialdata_attrs(metadata, comprehensive_metadata_obj)

        # Add comprehensive dataset metadata if supported
        self._add_comprehensive_metadata(metadata)

    def _setup_spatialdata_attrs(
        self, metadata: "SpatialData", comprehensive_metadata_obj
    ) -> None:
        """Setup SpatialData attributes with pixel size and metadata."""
        if not hasattr(metadata, "attrs") or metadata.attrs is None:
            metadata.attrs = {}

        logging.info("Adding comprehensive metadata to SpatialData.attrs")

        # Create pixel size attributes
        pixel_size_attrs = self._create_pixel_size_attrs()

        # Add comprehensive metadata sections
        self._add_comprehensive_sections(pixel_size_attrs, comprehensive_metadata_obj)

        # Update SpatialData attributes
        metadata.attrs.update(pixel_size_attrs)

    def _create_pixel_size_attrs(self) -> Dict[str, Any]:
        """Create pixel size and conversion metadata attributes."""
        # Import version dynamically
        try:
            from ... import __version__

            version = __version__
        except ImportError:
            version = "unknown"

        # Base pixel size metadata
        pixel_size_attrs = {
            "pixel_size_x_um": float(self.pixel_size_um),
            "pixel_size_y_um": float(self.pixel_size_um),
            "pixel_size_units": "micrometers",
            "coordinate_system": "physical_micrometers",
            "msi_converter_version": version,
            "conversion_timestamp": pd.Timestamp.now().isoformat(),
        }

        # Add pixel size detection provenance if available
        if self._pixel_size_detection_info is not None:
            pixel_size_attrs["pixel_size_detection_info"] = dict(
                self._pixel_size_detection_info
            )
            logging.info(
                f"Added pixel size detection info: "
                f"{self._pixel_size_detection_info}"
            )

        # Add conversion metadata
        assert self._dimensions is not None
        pixel_size_attrs["msi_dataset_info"] = {
            "dataset_id": self.dataset_id,
            "total_grid_pixels": self._dimensions[0]
            * self._dimensions[1]
            * self._dimensions[2],
            "non_empty_pixels": self._non_empty_pixel_count,
            "dimensions_xyz": list(self._dimensions),
        }

        return pixel_size_attrs

    def _add_comprehensive_sections(
        self, pixel_size_attrs: Dict[str, Any], comprehensive_metadata_obj
    ) -> None:
        """Add comprehensive metadata sections to attributes."""
        if comprehensive_metadata_obj.format_specific:
            pixel_size_attrs["format_specific_metadata"] = (
                comprehensive_metadata_obj.format_specific
            )

        if comprehensive_metadata_obj.acquisition_params:
            pixel_size_attrs["acquisition_parameters"] = (
                comprehensive_metadata_obj.acquisition_params
            )

        if comprehensive_metadata_obj.instrument_info:
            pixel_size_attrs["instrument_information"] = (
                comprehensive_metadata_obj.instrument_info
            )

    def _add_comprehensive_metadata(self, metadata: "SpatialData") -> None:
        """Add comprehensive dataset metadata if SpatialData supports it."""
        if not hasattr(metadata, "metadata"):
            return

        # Start with structured metadata from base class
        metadata_dict = self._structured_metadata.copy()

        # Add SpatialData-specific enhancements
        metadata_dict["non_empty_pixels"] = self._non_empty_pixel_count  # type: ignore[assignment]
        metadata_dict.update(
            {
                "spatialdata_specific": {
                    "zarr_compression_level": self.compression_level,
                    "tables_count": len(getattr(metadata, "tables", {})),
                    "shapes_count": len(getattr(metadata, "shapes", {})),
                    "images_count": len(getattr(metadata, "images", {})),
                },
            }
        )

        # Add pixel size detection provenance if available
        if self._pixel_size_detection_info is not None:
            metadata_dict["pixel_size_provenance"] = self._pixel_size_detection_info

        # Add conversion options used
        metadata_dict["conversion_options"] = {
            "handle_3d": self.handle_3d,
            "pixel_size_um": self.pixel_size_um,
            "dataset_id": self.dataset_id,
            **self.options,
        }

        metadata.metadata = metadata_dict

        logging.info(
            f"Comprehensive metadata persisted to SpatialData with "
            f"{len(metadata_dict)} top-level sections"
        )

    @abstractmethod
    def _create_data_structures(self) -> Dict[str, Any]:
        """Create data structures for the specific converter type."""
        pass

    @abstractmethod
    def _finalize_data(self, data_structures: Dict[str, Any]) -> None:
        """Finalize data structures for the specific converter type."""
        pass
