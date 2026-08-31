"""Metadata handling for MSI data formats.

This module provides a comprehensive metadata system for mass spectrometry
imaging (MSI) data,
featuring two-phase extraction (essential and comprehensive) optimized for
interpolation and
analysis workflows.

Architecture:
- Two-phase extraction: Fast essential metadata for processing, comprehensive
  for analysis
- Format-specific extractors: Optimized for ImzML and Bruker formats
- Structured data types: Type-safe dataclasses for metadata representation
- Caching and lazy loading: Efficient metadata access patterns

Core Components:
- EssentialMetadata: Critical metadata for processing decisions and
  interpolation setup
- ComprehensiveMetadata: Complete metadata including format-specific details
- MetadataExtractor: Abstract base class for format-specific extractors
- Format-specific extractors: ImzMLMetadataExtractor, BrukerMetadataExtractor

Example usage:
    >>> from thyra.metadata import EssentialMetadata, \
    ...     ComprehensiveMetadata
    >>> from thyra.metadata.extractors import ImzMLMetadataExtractor
    >>> from pyimzml.ImzMLParser import ImzMLParser
    >>>
    >>> # Create extractor and get essential metadata (fast)
    >>> parser = ImzMLParser("data.imzML")
    >>> extractor = ImzMLMetadataExtractor(parser, Path("data.imzML"))
    >>> essential = extractor.get_essential()
    >>>
    >>> # Access key properties for interpolation
    >>> print(f"Dimensions: {essential.dimensions}")
    >>> print(f"Pixel size: {essential.pixel_size}")
    >>> print(f"Mass range: {essential.mass_range}")
    >>> print(f"Is 3D: {essential.is_3d}")
    >>>
    >>> # Get comprehensive metadata when needed (slower)
    >>> comprehensive = extractor.get_comprehensive()
    >>> print(f"Instrument info: {comprehensive.instrument_info}")
"""

# Format-specific extractors (re-export from extractors submodule)
from .extractors import (
    BrukerMetadataExtractor,
    ImzMLMetadataExtractor,
    get_extractor_for_format,
    list_supported_formats,
)

# Core data types
from .types import ComprehensiveMetadata, EssentialMetadata

# Read-side compatibility for stores written by earlier versions
from .uns_compat import sanitize_uns_string_arrays

# Base classes - import delayed to avoid circular imports


# Public API
__all__ = [
    # Core data types
    "EssentialMetadata",
    "ComprehensiveMetadata",
    # Format-specific extractors
    "ImzMLMetadataExtractor",
    "BrukerMetadataExtractor",
    # Utility functions
    "get_extractor_for_format",
    "list_supported_formats",
    "sanitize_uns_string_arrays",
]

# Version information for the metadata system
__version__ = "2.0.0"
__metadata_api_version__ = "2.0"


def create_extractor(format_name: str, *args, **kwargs):
    """Factory function to create a metadata extractor for a given format.

    This is a convenience function that combines format detection with
    extractor creation.

    Args:
        format_name: Format identifier (e.g., 'imzml', 'bruker', 'tsf', 'tdf')
        *args: Positional arguments to pass to the extractor constructor
        **kwargs: Keyword arguments to pass to the extractor constructor

    Returns:
        Initialized metadata extractor instance

    Raises:
        ValueError: If format is not supported

    Example:
        >>> # For ImzML
        >>> from pyimzml.ImzMLParser import ImzMLParser
        >>> parser = ImzMLParser("data.imzML")
        >>> extractor = create_extractor('imzml', parser, Path("data.imzML"))
        >>>
        >>> # For Bruker
        >>> import sqlite3
        >>> conn = sqlite3.connect("data.d/analysis.tsf")
        >>> extractor = create_extractor('bruker', conn, Path("data.d"))
    """
    extractor_class = get_extractor_for_format(format_name)
    return extractor_class(*args, **kwargs)


def get_metadata_summary(extractor) -> dict:
    """Get a summary of metadata information for quick inspection.

    This function provides a convenient way to get key metadata information
    without needing to access individual properties.

    Args:
        extractor: Initialized metadata extractor instance

    Returns:
        Dictionary with summary information

    Example:
        >>> summary = get_metadata_summary(extractor)
        >>> print(f"Dataset: {summary['dimensions']} pixels, "
        ...       f"{summary['mass_range']} m/z")
    """
    essential = extractor.get_essential()

    return {
        "dimensions": essential.dimensions,
        "coordinate_bounds": essential.coordinate_bounds,
        "mass_range": essential.mass_range,
        "pixel_size": essential.pixel_size,
        "n_spectra": essential.n_spectra,
        "estimated_memory_gb": essential.estimated_memory_gb,
        "source_path": essential.source_path,
        "is_3d": essential.is_3d,
        "has_pixel_size": essential.has_pixel_size,
    }
