# thyra/convert.py
import logging
import traceback
import warnings
from pathlib import Path
from typing import Any, Dict, Literal, Optional, Tuple, Union

from .core.base_converter import PixelSizeSource
from .core.registry import detect_format, get_converter_class, get_reader_class

logger = logging.getLogger(__name__)

warnings.filterwarnings(
    "ignore",
    message=r"Accession IMS:1000046.*",  # ignore UserWarning
    category=UserWarning,
    module=r"pyimzml.ontology.ontology",
)

# warnings.filterwarnings(
#     "ignore",
#     category=CryptographyDeprecationWarning
# )


def _validate_paths_parameters(
    input_path: Union[str, Path], output_path: Union[str, Path]
) -> bool:
    """Validate path parameters."""
    if not input_path or not isinstance(input_path, (str, Path)):
        logger.error("Input path must be a valid string or Path object")
        return False

    if not output_path or not isinstance(output_path, (str, Path)):
        logger.error("Output path must be a valid string or Path object")
        return False

    return True


def _validate_string_parameters(format_type: str, dataset_id: str) -> bool:
    """Validate string parameters."""
    if not isinstance(format_type, str) or not format_type.strip():
        logger.error("Format type must be a non-empty string")
        return False

    if not isinstance(dataset_id, str) or not dataset_id.strip():
        logger.error("Dataset ID must be a non-empty string")
        return False

    return True


def _validate_numeric_parameters(
    pixel_size_um: Optional[float], handle_3d: bool
) -> bool:
    """Validate numeric and boolean parameters."""
    if pixel_size_um is not None and (
        not isinstance(pixel_size_um, (int, float)) or pixel_size_um <= 0
    ):
        logger.error("Pixel size must be a positive number")
        return False

    return True


def _validate_input_parameters(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    format_type: str,
    dataset_id: str,
    pixel_size_um: Optional[float],
    handle_3d: bool,
) -> bool:
    """Validate all input parameters for convert_msi function."""
    return (
        _validate_paths_parameters(input_path, output_path)
        and _validate_string_parameters(format_type, dataset_id)
        and _validate_numeric_parameters(pixel_size_um, handle_3d)
    )


def _validate_paths(input_path: Path, output_path: Path) -> bool:
    """Validate that input exists and output doesn't exist."""
    if not input_path.exists():
        logger.error(f"Input path does not exist: {input_path}")
        return False

    if output_path.exists():
        logger.error(f"Destination {output_path} already exists.")
        return False

    return True


def _create_reader(
    input_path: Path, reader_options: Optional[Dict[str, Any]] = None
) -> Tuple[Any, str]:
    """Create and return a reader for the input format.

    Args:
        input_path: Path to the input MSI data
        reader_options: Optional format-specific reader options (e.g., calibration settings)

    Returns:
        Tuple of (reader instance, detected format string)
    """
    input_format = detect_format(input_path)
    logger.info(f"Detected format: {input_format}")
    reader_class = get_reader_class(input_format)
    logger.info(f"Using reader: {reader_class.__name__}")

    # Pass reader options to the reader if provided
    options = reader_options or {}
    return reader_class(input_path, **options), input_format


def _determine_pixel_size(
    reader: Any, pixel_size_um: Optional[float], input_format: str
) -> Tuple[float, PixelSizeSource, Dict[str, Any]]:
    """Determine pixel size either from metadata or user input."""
    if pixel_size_um is not None:
        # Manual pixel size was provided
        pixel_size_detection_info = {
            "method": "manual",
            "source_format": input_format,
            "detection_successful": False,
            "note": "Pixel size manually specified via --pixel-size parameter",
        }
        return pixel_size_um, PixelSizeSource.USER_PROVIDED, pixel_size_detection_info

    # Attempt automatic detection
    logger.info("Attempting automatic pixel size detection...")
    essential_metadata = reader.get_essential_metadata()

    if essential_metadata.pixel_size is None:
        logger.error("Pixel size not found in metadata")
        logger.error("Use --pixel-size parameter (e.g., --pixel-size 25)")
        raise ValueError("Pixel size not found in metadata")

    final_pixel_size = essential_metadata.pixel_size[0]  # Use X size
    logger.info(f"Detected pixel size: {final_pixel_size:.1f} um")

    pixel_size_detection_info = {
        "method": "automatic",
        "detected_x_um": float(essential_metadata.pixel_size[0]),
        "detected_y_um": float(essential_metadata.pixel_size[1]),
        "source_format": input_format,
        "detection_successful": True,
        "note": "Pixel size automatically detected from source metadata",
    }

    return final_pixel_size, PixelSizeSource.AUTO_DETECTED, pixel_size_detection_info


def _should_use_streaming(streaming: Union[bool, Literal["auto"]], reader: Any) -> bool:
    """Determine if streaming converter should be used."""
    if streaming is True:
        return True
    if streaming != "auto":
        return False

    # Auto-detect based on estimated dataset size (>10GB)
    try:
        essential_meta = reader.get_essential_metadata()
        dims = essential_meta.dimensions
        n_pixels = dims[0] * dims[1] * dims[2]
        # Rough estimate: assume average 10k peaks per spectrum, 8 bytes each
        estimated_gb = (n_pixels * 10000 * 8) / (1024**3)
        if estimated_gb > 10:
            logger.info(
                f"Auto-detected large dataset (~{estimated_gb:.1f} GB), "
                "using streaming converter"
            )
            return True
    except Exception as e:
        logger.debug(f"Could not estimate dataset size for auto-streaming: {e}")
    return False


def _create_converter(
    format_type: str,
    reader: Any,
    output_path: Path,
    dataset_id: str,
    pixel_size_um: float,
    pixel_size_source: PixelSizeSource,
    handle_3d: bool,
    pixel_size_detection_info: Dict[str, Any],
    resampling_config: Optional[Dict[str, Any]] = None,
    sparse_format: str = "csc",
    include_optical: bool = True,
    streaming: Union[bool, Literal["auto"]] = "auto",
    **kwargs: Any,
) -> Any:
    """Create and return a converter for the specified format."""
    converter_kwargs = {
        "dataset_id": dataset_id,
        "pixel_size_um": pixel_size_um,
        "pixel_size_source": pixel_size_source,
        "handle_3d": handle_3d,
        "pixel_size_detection_info": pixel_size_detection_info,
        "resampling_config": resampling_config,
        "include_optical": include_optical,
        **kwargs,
    }

    # Try streaming converter if requested
    if (
        _should_use_streaming(streaming, reader)
        and "spatialdata" in format_type.lower()
    ):
        try:
            from .converters.spatialdata import StreamingSpatialDataConverter

            logger.info("Using streaming converter for memory-efficient processing")
            return StreamingSpatialDataConverter(
                reader, output_path, **converter_kwargs
            )
        except ImportError as e:
            logger.warning(f"Streaming converter not available: {e}")
            logger.warning("Falling back to standard converter")

    try:
        converter_class = get_converter_class(format_type.lower())
        logger.info(f"Using converter: {converter_class.__name__}")
    except ValueError as e:
        if "spatialdata" in format_type.lower():
            logger.error(
                "SpatialData converter is not available due to " "dependency issues."
            )
            logger.error("This is commonly caused by zarr version incompatibility.")
            logger.error("Try upgrading your dependencies:")
            logger.error("  pip install --upgrade anndata spatialdata zarr")
            logger.error("Or create a fresh environment with compatible versions.")
            raise ValueError("SpatialData converter unavailable") from e
        else:
            raise e
    return converter_class(
        reader,
        output_path,
        dataset_id=dataset_id,
        pixel_size_um=pixel_size_um,
        pixel_size_source=pixel_size_source,
        handle_3d=handle_3d,
        pixel_size_detection_info=pixel_size_detection_info,
        resampling_config=resampling_config,
        sparse_format=sparse_format,
        include_optical=include_optical,
        **kwargs,
    )


def _perform_conversion_with_cleanup(converter: Any, reader: Any) -> bool:
    """Perform the conversion and handle reader cleanup."""
    try:
        logger.info("Starting conversion...")
        result = converter.convert()
        logger.info(f"Conversion {'completed successfully' if result else 'failed'}")
        return bool(result)
    finally:
        if hasattr(reader, "close"):
            reader.close()


def convert_msi(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    format_type: str = "spatialdata",
    dataset_id: str = "msi_dataset",
    pixel_size_um: Optional[float] = None,
    handle_3d: bool = False,
    resampling_config: Optional[Dict[str, Any]] = None,
    reader_options: Optional[Dict[str, Any]] = None,
    sparse_format: str = "csc",
    include_optical: bool = True,
    streaming: Union[bool, Literal["auto"]] = "auto",
    region: Optional[int] = None,
    **kwargs: Any,
) -> bool:
    """Convert MSI data to the specified format.

    Provides automatic pixel size detection from metadata or
    accepts user-specified values.

    Args:
        input_path: Path to input MSI data file or directory
        output_path: Path for output file
        format_type: Output format type (default: "spatialdata")
        dataset_id: Identifier for the dataset
        pixel_size_um: Pixel size in micrometers (None for auto)
        handle_3d: Whether to process as 3D data (default: False)
        resampling_config: Optional resampling configuration
        reader_options: Optional format-specific reader options:
            - intensity_threshold: float - Minimum intensity
              to include. Default: None (no filtering).
            - use_recalibrated_state: bool - For Bruker data,
              use active/recalibrated calibration (default True).
        sparse_format: Sparse matrix format ('csc' or 'csr')
        include_optical: Include optical images (default: True)
        streaming: Use streaming converter for large datasets.
            - "auto": Auto-detect based on dataset size >10GB (default)
            - True: Force streaming converter
            - False: Force standard converter
        region: For multi-region datasets (e.g. Bruker timsTOF),
            select a specific region number. None (default)
            converts all regions. Passed to the reader as
            reader_options["region"].
        **kwargs: Additional keyword arguments

    Returns:
        True if conversion was successful, False otherwise
    """
    # Validate input parameters
    if not _validate_input_parameters(
        input_path,
        output_path,
        format_type,
        dataset_id,
        pixel_size_um,
        handle_3d,
    ):
        return False

    # Convert to Path objects and validate
    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    logger.info(f"Processing input file: {input_path}")

    if not _validate_paths(input_path, output_path):
        return False

    # Merge region into reader_options if provided
    if region is not None:
        reader_options = dict(reader_options or {})
        reader_options["region"] = region

    try:
        # Create reader with format-specific options
        reader, input_format = _create_reader(input_path, reader_options)

        # Determine pixel size
        final_pixel_size, pixel_size_source, pixel_size_detection_info = (
            _determine_pixel_size(reader, pixel_size_um, input_format)
        )

        # Create converter
        converter = _create_converter(
            format_type,
            reader,
            output_path,
            dataset_id,
            final_pixel_size,
            pixel_size_source,
            handle_3d,
            pixel_size_detection_info,
            resampling_config,
            sparse_format,
            include_optical=include_optical,
            streaming=streaming,
            **kwargs,
        )

        # Perform conversion with cleanup
        return _perform_conversion_with_cleanup(converter, reader)

    except Exception as e:
        logger.error(f"Error during conversion: {e}")
        logger.error(f"Detailed traceback:\n{traceback.format_exc()}")
        return False
