"""Two phases, because one of them costs a pass over the data.

Metadata splits by price, not by subject. ``get_essential()`` returns what
the conversion needs to set itself up -- grid, bounds, mass range, pixel
size, counts -- and the caller pays for it on every run.
``get_comprehensive()`` returns everything the source records, for
provenance, and on several formats that means parsing far more of the
file. Keeping them apart is what lets ``preview_msi`` describe an
acquisition without decoding it.

The split is a price, not a promise about cost: "essential" is cheap for a
format whose header states the counts and expensive for one that does not.
imzML scanned every spectrum in processed mode to find the true mass range,
and PHI has to aggregate the whole ion-event stream to learn which pixels
fired at all -- which is why both offer a way to skip the counting and say
so in the result rather than guessing (``n_spectra_counted``, issue #240;
``mass_range_known``, issue #360).  On imzML that way is a whole second
extractor, :class:`ImzMLHeaderExtractor`, because the expensive half is
the spectrum list itself and not one query within it.

Both results are cached on the instance, so the expensive phase runs once
per reader. :meth:`MetadataExtractor.clear_cache` exists for the case
where the underlying source changed beneath a live extractor.
"""

import logging
from abc import ABC, abstractmethod
from typing import Any, Optional

from ..metadata.personal_data import without_personal_data
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

        The vendor dictionaries come back without the fields that name a
        person, and with a path the vendor recorded reduced to its file
        name (see :mod:`thyra.metadata.personal_data`). That happens here,
        for every format at once, so that no extractor has to remember it.

        Returns:
            ComprehensiveMetadata: Complete metadata including
            format-specific details
        """
        if self._comprehensive_cache is None:
            logger.info("Extracting comprehensive metadata...")
            # Ensure essential metadata is loaded first
            self.get_essential()
            self._comprehensive_cache = without_personal_data(
                self._extract_comprehensive_impl()
            )
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
