# thyra/metadata/core/base_extractor.py
import logging
from abc import ABC, abstractmethod
from typing import Any, Optional

from ..metadata.types import ComprehensiveMetadata, EssentialMetadata

logger = logging.getLogger(__name__)


class MetadataExtractor(ABC):
    """Abstract base class for format-specific metadata extractors."""

    def __init__(self, data_source: Any):
        """Initialize metadata extractor with data source.

        Args:
            data_source: Format-specific data source (parser, connection, etc.)
        """
        self.data_source = data_source
        self._essential_cache: Optional[EssentialMetadata] = None
        self._comprehensive_cache: Optional[ComprehensiveMetadata] = None

    @abstractmethod
    def _extract_essential_impl(self) -> EssentialMetadata:
        """Format-specific implementation of essential metadata extraction.

        This method should be optimized for speed and extract only the minimum
        metadata needed for processing decisions and interpolation setup.

        Returns:
            EssentialMetadata: Critical metadata for processing
        """
        pass

    @abstractmethod
    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        """Format-specific implementation of comprehensive metadata extraction.

        This method can be slower and should extract all available metadata
        including format-specific details, acquisition parameters, etc.

        Returns:
            ComprehensiveMetadata: Complete metadata including
            format-specific details
        """
        pass

    def get_essential(self) -> EssentialMetadata:
        """Get essential metadata (cached after first call).

        Returns:
            EssentialMetadata: Critical metadata for processing
        """
        if self._essential_cache is None:
            logger.info("Extracting essential metadata...")
            self._essential_cache = self._extract_essential_impl()
            logger.debug(
                f"Essential metadata extracted: "
                f"{self._essential_cache.dimensions} dimensions, "
                f"{self._essential_cache.n_spectra} spectra"
            )
        return self._essential_cache

    def get_comprehensive(self) -> ComprehensiveMetadata:
        """Get comprehensive metadata (cached after first call).

        Returns:
            ComprehensiveMetadata: Complete metadata including
            format-specific details
        """
        if self._comprehensive_cache is None:
            logger.info("Extracting comprehensive metadata...")
            # Ensure essential metadata is loaded first
            self.get_essential()
            self._comprehensive_cache = self._extract_comprehensive_impl()
            logger.debug(
                f"Comprehensive metadata extracted with "
                f"{len(self._comprehensive_cache.raw_metadata)} raw entries"
            )
        return self._comprehensive_cache

    def clear_cache(self) -> None:
        """Clear cached metadata to force re-extraction."""
        self._essential_cache = None
        self._comprehensive_cache = None
        logger.debug("Metadata cache cleared")
