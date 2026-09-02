# thyra/metadata/extractors/__init__.py
"""Format-specific metadata extractors for MSI data.

This module provides specialized metadata extractors for different MSI data
formats,
each optimized for efficient two-phase metadata extraction (essential and
comprehensive).

Available Extractors:
- ImzMLMetadataExtractor: For ImzML format files
- BrukerMetadataExtractor: For Bruker TSF/TDF format files
- WatersMetadataExtractor: For Waters .raw imaging files

Example usage:
    >>> from thyra.metadata.extractors import \
    ...     ImzMLMetadataExtractor, BrukerMetadataExtractor
    >>> from pyimzml.ImzMLParser import ImzMLParser
    >>>
    >>> # For ImzML files
    >>> parser = ImzMLParser("data.imzML")
    >>> extractor = ImzMLMetadataExtractor(parser, Path("data.imzML"))
    >>> essential = extractor.get_essential()
    >>>
    >>> # For Bruker files
    >>> import sqlite3
    >>> conn = sqlite3.connect("data.d/analysis.tsf")
    >>> extractor = BrukerMetadataExtractor(conn, Path("data.d"))
    >>> essential = extractor.get_essential()
"""

from .bruker_extractor import BrukerMetadataExtractor
from .imzml_extractor import ImzMLMetadataExtractor
from .mzpeak_extractor import MzPeakMetadataExtractor
from .phi_extractor import PhiMetadataExtractor
from .solarix_extractor import SolarixMetadataExtractor
from .waters_extractor import WatersMetadataExtractor

# Public API
__all__ = [
    "BrukerMetadataExtractor",
    "ImzMLMetadataExtractor",
    "MzPeakMetadataExtractor",
    "PhiMetadataExtractor",
    "SolarixMetadataExtractor",
    "WatersMetadataExtractor",
]

# Format mapping for dynamic extraction
FORMAT_EXTRACTORS = {
    "imzml": ImzMLMetadataExtractor,
    "bruker": BrukerMetadataExtractor,
    "tsf": BrukerMetadataExtractor,
    "tdf": BrukerMetadataExtractor,
    "solarix": SolarixMetadataExtractor,
    "waters": WatersMetadataExtractor,
    "phi": PhiMetadataExtractor,
    "mzpeak": MzPeakMetadataExtractor,
}


def get_extractor_for_format(format_name: str):
    """Get the appropriate metadata extractor class for a given format.

    Args:
        format_name: Format identifier (e.g., 'imzml', 'bruker', 'tsf', 'tdf')

    Returns:
        Metadata extractor class for the specified format

    Raises:
        ValueError: If format is not supported

    Example:
        >>> extractor_class = get_extractor_for_format('imzml')
        >>> extractor = extractor_class(parser, path)
    """
    format_lower = format_name.lower()
    if format_lower not in FORMAT_EXTRACTORS:
        supported_formats = ", ".join(FORMAT_EXTRACTORS.keys())
        raise ValueError(
            f"Unsupported format '{format_name}'. "
            f"Supported formats: {supported_formats}"
        )
    return FORMAT_EXTRACTORS[format_lower]


def list_supported_formats() -> list[str]:
    """List all supported metadata extraction formats.

    Returns:
        List of supported format names
    """
    return list(FORMAT_EXTRACTORS.keys())
