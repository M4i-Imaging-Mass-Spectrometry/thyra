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
    from .frames import FrameScans
    from .mobility import MobilityAxis
    from .msms import FragmentationSchedule

from .mass_axis import validate_max_mass_axis_length

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
            **kwargs: Additional reader-specific parameters.
                ``max_mass_axis_length`` caps how many unique m/z values a
                raw axis may reach before the build gives up. It lives here
                rather than on one reader because five of them build such an
                axis through the same accumulator (issue #294); it defaults
                to ``None``, no cap, everywhere except imzML, which sets its
                own documented default. A reader that does not build a raw
                axis simply never consults it.
        """
        self.data_path = Path(data_path)
        self._intensity_threshold = intensity_threshold
        self._metadata_extractor: Optional["MetadataExtractor"] = None
        # Validated here so a bad value fails while the caller is still
        # looking at their own arguments, which is what the imzML reader
        # already did with it.
        self.max_mass_axis_length: Optional[int] = validate_max_mass_axis_length(
            kwargs.get("max_mass_axis_length")
        )

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

        Returns the image files that hold optical images of the sample --
        TIFF, JPEG, PNG or BMP. These images can be stored alongside MSI data
        in SpatialData output for multimodal analysis.

        Default implementation returns empty list. Subclasses should override
        to return paths to optical images specific to their format.

        Returns:
            Paths to optical images, empty if none are available. When one of
            them is the alignment image (see
            :meth:`get_primary_optical_image_path`) it comes first.
        """
        return []

    def get_primary_optical_image_path(self) -> Optional[Path]:
        """Get the optical image the source designates as the alignment image.

        Some formats say which of an acquisition's images the registration is
        stated against -- FlexImaging names it in the .mis ``<ImageFile>``
        element -- and that one gets the identity transform while the others
        are scaled into its pixel space. Without it the converter has no way
        to tell an alignment scan from a slide overview that happens to sit
        in the same folder.

        Default implementation returns None, which leaves every optical image
        equal. Override alongside :meth:`get_optical_image_paths` when the
        format designates one; it must be a member of that list.

        Returns:
            Path to the designated image, or None if the format does not
            designate one or the designated file is not present.
        """
        return None

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
    def iter_spectra(self) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Iterate through spectra with optional batch processing.

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

    def get_mass_axis_annotations(
        self,
    ) -> Optional[dict]:
        """Get extra per-channel columns to store alongside the m/z axis.

        Returns a mapping of column name to an array with one entry per
        entry of :meth:`get_common_mass_axis`. These are written into the
        table's ``var`` next to ``mz``.

        This exists so a format whose native axis is not m/z can keep that
        axis in the output. Time-of-flight instruments measure flight time
        and derive m/z from a calibration, so storing the flight time makes
        the conversion reversible without the reader: a later recalibration
        can be applied to the stored times directly.

        Annotations are dropped if their length does not match the axis the
        converter actually writes, which is what happens when resampling is
        enabled and the axis is rebuilt.

        Returns:
            Mapping of column name to per-channel values, or None.
        """
        return None

    # ------------------------------------------------------------------
    # One read per frame for every table
    #
    # A source whose summed spectrum, mobility point cloud and precursor
    # spectra are all functions of one raw read per pixel can hand that
    # read over once, as a record, so a converter that writes several
    # tables reads the source once per pass rather than once per table.
    # See ``thyra.core.frames`` and design decision D5.
    # ------------------------------------------------------------------

    @property
    def has_frame_scans(self) -> bool:
        """Whether :meth:`iter_frame_scans` is available on this source."""
        return False

    def iter_frame_scans(self) -> Generator["FrameScans", None, None]:
        """Iterate the frames as records that derive every view of a pixel.

        Yields one :class:`~thyra.core.frames.FrameScans` per frame, in
        the same order and with the same coordinates as
        :meth:`iter_spectra`, including frames whose summed spectrum is
        empty (their :meth:`~thyra.core.frames.FrameScans.spectrum` is
        ``None``), so a consumer can still count them for the sinks that
        do not depend on the summed spectrum.

        Raises:
            NotImplementedError: When :attr:`has_frame_scans` is False.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not hand its frames over as records"
        )

    # ------------------------------------------------------------------
    # Ion mobility
    #
    # Mobility is a second coordinate on a feature, next to m/z. It never
    # enters obs, a coordinate system or a transform. ``iter_spectra`` keeps
    # yielding the mobility-summed spectrum, so nothing downstream learns
    # about mobility unless it asks through the methods below.
    # ------------------------------------------------------------------

    @property
    def has_ion_mobility(self) -> bool:
        """Whether the source carries a mobility dimension at all."""
        return False

    @property
    def has_shared_mobility_axis(self) -> bool:
        """Whether every pixel carries the same (m/z, mobility) feature pairs.

        True for a continuous imzML export whose mobility array is the same
        for every spectrum (TIMSImaging, TIMSCONVERT in continuous mode):
        the pairs are then a feature list the converter can write as a
        mobility-resolved table directly. False when each pixel has its
        own mobility values, which need a common mobility grid first.
        """
        return False

    def get_mobility_axis(self) -> Optional["MobilityAxis"]:
        """Describe the mobility dimension, or ``None`` when there is none."""
        return None

    def get_shared_mobility_features(
        self,
    ) -> Optional[Tuple[NDArray[np.float64], NDArray[np.float64]]]:
        """The shared ``(mz, mobility)`` feature pairs, when the axis is shared.

        Both arrays have one entry per feature, in the source's order;
        duplicate m/z values are expected -- that is what mobility splits.
        ``None`` when :attr:`has_shared_mobility_axis` is False.
        """
        return None

    def iter_mobility_spectra(self) -> Generator[
        Tuple[
            Tuple[int, int, int],
            NDArray[np.float64],
            NDArray[np.float64],
            NDArray[np.float64],
        ],
        None,
        None,
    ]:
        """Iterate through spectra with their per-point mobility values.

        Yields ``((x, y, z), mzs, mobilities, intensities)`` per pixel: the
        raw (m/z, mobility) point cloud, unbinned, with the three arrays
        parallel. Readers whose source has no mobility dimension do not
        implement this.

        Raises:
            NotImplementedError: If the reader exposes no mobility dimension.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not expose an ion mobility dimension"
        )

    # ------------------------------------------------------------------
    # Fragmentation
    #
    # Whether the spectra are fragment spectra, and of what. Like mobility,
    # this changes nothing about what ``iter_spectra`` yields -- an MS/MS
    # acquisition still comes through as one spectrum per pixel. It is
    # recorded so a store can say that its m/z axis is fragment m/z, which
    # the axis itself cannot.
    # ------------------------------------------------------------------

    def get_fragmentation(self) -> Optional["FragmentationSchedule"]:
        """Describe the fragmentation, or ``None`` when the reader cannot.

        ``None`` means "not reported", which is not the same as "MS1": a
        reader that has no way to tell says nothing rather than claiming
        the acquisition was unfragmented.
        """
        return None

    def iter_precursor_spectra(self) -> Generator[
        Tuple[
            Tuple[int, int, int],
            int,
            NDArray[np.float64],
            NDArray[np.float64],
        ],
        None,
        None,
    ]:
        """Iterate the fragment spectra of each precursor separately.

        Yields ``((x, y, z), window_index, mzs, intensities)``: one
        spectrum per pixel *per precursor*, where ``window_index`` is the
        position of that precursor in
        :meth:`get_fragmentation`'s ``windows``. Only a source that
        separates its precursors within a pixel implements this -- Bruker
        PASEF gives each isolation window its own slice of the mobility
        ramp, so the split is a filter on the scan number rather than a
        deconvolution.

        A pixel-precursor pair with no ion current is not yielded.
        :meth:`iter_spectra` is unaffected and keeps yielding the summed
        spectrum.

        Raises:
            NotImplementedError: If the reader cannot separate precursors.
        """
        raise NotImplementedError(
            f"{type(self).__name__} cannot separate the precursors of a pixel"
        )

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
