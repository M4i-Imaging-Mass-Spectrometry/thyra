# thyra/core/base_reader.py
import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Generator, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..metadata.types import ComprehensiveMetadata, EssentialMetadata

if TYPE_CHECKING:
    from .base_extractor import MetadataExtractor

logger = logging.getLogger(__name__)


class BaseMSIReader(ABC):
    """Abstract base class for reading MSI data formats."""

    def __init__(
        self,
        data_path: Path,
        intensity_threshold: Optional[float] = None,
        **kwargs: object,
    ):
        """Initialize the reader with the path to the data.

        Args:
            data_path: Path to the data file or directory
            intensity_threshold: Minimum intensity value to include.
                Values below this threshold are filtered out during iteration.
                Useful for removing detector noise in continuous mode data.
                Default: None (no filtering, include all values).
            **kwargs: Additional reader-specific parameters
        """
        self.data_path = Path(data_path)
        self._intensity_threshold = intensity_threshold
        self._metadata_extractor: Optional["MetadataExtractor"] = None

        if intensity_threshold is not None:
            logger.info(
                f"Intensity threshold active: values < {intensity_threshold} "
                "will be filtered out"
            )

    @abstractmethod
    def _create_metadata_extractor(self) -> "MetadataExtractor":
        """Create format-specific metadata extractor."""
        pass

    @property
    def metadata_extractor(self) -> "MetadataExtractor":
        """Lazy-loaded metadata extractor."""
        if self._metadata_extractor is None:
            self._metadata_extractor = self._create_metadata_extractor()
        return self._metadata_extractor

    def get_essential_metadata(self) -> EssentialMetadata:
        """Get essential metadata for processing."""
        return self.metadata_extractor.get_essential()

    def get_comprehensive_metadata(self) -> ComprehensiveMetadata:
        """Get complete metadata."""
        return self.metadata_extractor.get_comprehensive()

    def get_optical_image_paths(self) -> List[Path]:
        """Get paths to optical/microscopy images associated with this data.

        Returns list of TIFF file paths that contain optical images of the
        sample. These images can be stored alongside MSI data in SpatialData
        output for multimodal analysis.

        Default implementation returns empty list. Subclasses should override
        to return paths to optical images specific to their format.

        Returns:
            List of paths to TIFF files, empty if no optical images available.
        """
        return []

    @abstractmethod
    def get_common_mass_axis(self) -> NDArray[np.float64]:
        """Return the common mass axis for all spectra.

        This method must always return a valid array. If no common mass
        axis can be created, implementations should raise an exception.
        """
        pass

    @property
    def has_shared_mass_axis(self) -> bool:
        """Check if all spectra share the same m/z axis.

        For continuous ImzML data, all spectra have identical m/z values,
        so get_common_mass_axis() only needs to read the first spectrum.
        For processed/centroid data, each spectrum may have different m/z
        values, requiring iteration through all spectra.

        Returns:
            True if all spectra share the same m/z axis (continuous mode),
            False if each spectrum has different m/z values (processed mode).
        """
        # Default implementation returns False (conservative assumption)
        # Subclasses should override if they can detect shared mass axis
        return False

    @abstractmethod
    def iter_spectra(self, batch_size: Optional[int] = None) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Iterate through spectra with optional batch processing.

        Args:
            batch_size: Optional batch size for spectrum iteration

        Yields:
            Tuple containing:

                - Coordinates (x, y, z) using 0-based indexing
                - m/z values array

                - Intensity values array

        Note:
            Subclasses should apply intensity threshold filtering by calling
            _apply_intensity_filter() on the intensities before yielding.
        """
        pass

    def _apply_intensity_filter(
        self,
        mzs: NDArray[np.float64],
        intensities: NDArray[np.float64],
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Apply intensity threshold filtering to spectrum data.

        Args:
            mzs: m/z values array
            intensities: Intensity values array

        Returns:
            Tuple of (filtered_mzs, filtered_intensities) with values below
            threshold removed. Returns original arrays if no threshold is set.
        """
        if self._intensity_threshold is None:
            return mzs, intensities

        mask = intensities >= self._intensity_threshold
        return mzs[mask], intensities[mask]

    def get_peak_counts_per_pixel(self) -> Optional[NDArray[np.int32]]:
        """Get per-pixel peak counts for CSR indptr construction.

        This method enables optimized streaming conversion by providing
        pre-computed peak counts, avoiding the need for a separate counting pass.

        Returns:
            Array of size n_pixels where arr[pixel_idx] = peak_count.
            pixel_idx = z * (n_x * n_y) + y * n_x + x
            Returns None if not supported/available for this reader.

        Note:
            Override in subclass to enable optimized streaming conversion.
            The default implementation returns None, which causes the
            streaming converter to fall back to a two-pass approach.

        Warning:
            When intensity_threshold is set, the actual peak counts after
            filtering may be lower than the values returned here, since this
            method typically returns pre-computed counts from metadata that
            don't account for intensity filtering. The streaming converter
            handles this gracefully by using a two-pass approach.
        """
        return None

    @staticmethod
    def map_mz_to_common_axis(
        mzs: NDArray[np.float64],
        intensities: NDArray[np.float64],
        common_axis: NDArray[np.float64],
    ) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
        """Map m/z values to indices in the common mass axis with high accuracy.

        This method ensures exact mapping of m/z values to the common mass axis
        without interpolation, preserving the original intensity values.

        Args:
            mzs: NDArray[np.float64] - Array of m/z values
            intensities: NDArray[np.float64] - Array of intensity values
            common_axis: NDArray[np.float64] - Common mass axis (sorted array
            of unique m/z values)

        Returns:
            Tuple of (indices in common mass axis, corresponding intensities)
        """
        if mzs.size == 0 or intensities.size == 0:
            return np.array([], dtype=int), np.array([])

        # Use searchsorted to find indices in common mass axis
        indices = np.searchsorted(common_axis, mzs)

        # Ensure indices are within bounds
        indices = np.clip(indices, 0, len(common_axis) - 1)

        # Verify that we're actually finding the right m/z values
        max_diff = 1e-6  # A very small tolerance threshold for floating
        # point differences
        indices_valid = np.abs(common_axis[indices] - mzs) <= max_diff

        # Return only the valid indices and their corresponding intensities
        return indices[indices_valid], intensities[indices_valid]

    def get_region_map(self) -> Optional[dict]:
        """Get per-pixel region mapping for multi-region datasets.

        Returns a dictionary mapping normalized (0-based) (x, y) coordinate
        tuples to integer region numbers. This enables the converter to annotate
        each pixel with its acquisition region in obs["region_number"].

        Default implementation returns None (single-region or no region info).
        Subclasses should override when region information is available.

        Returns:
            Dict mapping (x, y) tuples to region numbers, or None if
            region information is not available.
        """
        return None

    def get_region_info(self) -> Optional[list]:
        """Get summary information about acquisition regions.

        Returns a list of dictionaries, each describing one region with at
        minimum: {"region_number": int, "n_spectra": int}. Additional keys
        (e.g. "name") are format-specific and optional.

        Default implementation returns None (single-region or no region info).
        Subclasses should override when region information is available.

        Returns:
            List of region summary dicts, or None if region information
            is not available.
        """
        return None

    def reset(self) -> None:
        """Reset the reader to allow iterating from the beginning again.

        This method should reset any internal state so that iter_spectra()
        can be called again to iterate from the first spectrum.

        Default implementation does nothing (assumes reader state is stateless).
        Subclasses should override if they maintain iteration state.
        """
        pass

    @abstractmethod
    def close(self) -> None:
        """Close all open file handles."""
        pass

    def __enter__(self) -> "BaseMSIReader":
        """Context manager entry."""
        return self

    def __exit__(
        self,
        exc_type: Optional[type],
        exc_val: Optional[Exception],
        exc_tb: Optional[object],
    ) -> None:
        """Context manager exit with cleanup."""
        self.close()
