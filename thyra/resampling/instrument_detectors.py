"""Instrument detection strategies for resampling decisions.

This module implements the Strategy pattern for instrument detection.
Each detector class encapsulates the logic for identifying a specific
instrument type and returning appropriate resampling parameters.
"""

import logging
from abc import ABC, abstractmethod
from typing import List, Optional

from .data_characteristics import DataCharacteristics
from .types import AxisType, ResamplingMethod

logger = logging.getLogger(__name__)


class InstrumentDetector(ABC):
    """Abstract base class for instrument detection strategies.

    Each detector is responsible for:
    1. Detecting if data matches a specific instrument type
    2. Returning the appropriate resampling method
    3. Returning the appropriate mass axis type
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable name for logging."""
        pass

    @abstractmethod
    def matches(self, characteristics: DataCharacteristics) -> bool:
        """Check if data characteristics match this instrument.

        Args:
            characteristics: Extracted data characteristics

        Returns:
            True if this detector matches the data
        """
        pass

    @abstractmethod
    def get_resampling_method(self) -> ResamplingMethod:
        """Get the recommended resampling method for this instrument."""
        pass

    @abstractmethod
    def get_axis_type(self) -> AxisType:
        """Get the recommended mass axis type for this instrument."""
        pass

    @property
    def source_grid_law(self) -> Optional[AxisType]:
        """Spacing law of the *source* m/z grid, when it is known.

        ``tic_preserving`` interpolates onto the target axis and then applies
        one global scaling factor. That composite operator is exact only when
        the source grid follows the same spacing law as the target axis; off
        the diagonal it distorts the ratio between two peaks by exactly
        ``(m_hi / m_lo) ** (p_target - p_source)``, with ``p`` = 0, 0.5, 1,
        1.5, 2 for constant / linear_tof / reflector_tof / orbitrap / fticr.

        SCiLS Lab gates its own TIC-preserving resampling on the same
        condition -- "If all axis types are identical, a TIC preserving
        resampling is applied, otherwise a linear interpolation is performed"
        (SCiLS Lab 2026b User Guide, p.80).

        ``None`` -- the default -- means Thyra does not know the source law.
        That is the honest answer for every detector that matches on
        something other than a vendor format whose grid Thyra itself lays
        out, and :class:`InstrumentDetectorChain` refuses ``TIC_PRESERVING``
        for such a detector.
        """
        return None


class CentroidImzMLDetector(InstrumentDetector):
    """Detector for centroid ImzML data.

    Centroid data has discrete peaks and benefits from nearest-neighbor
    resampling with reflector TOF axis spacing (constant relative resolution).
    """

    @property
    def name(self) -> str:
        """Return detector name."""
        return "ImzML Centroid"

    def matches(self, characteristics: DataCharacteristics) -> bool:
        """Check if data is centroid spectrum type."""
        return characteristics.is_centroid_data

    def get_resampling_method(self) -> ResamplingMethod:
        """Return nearest-neighbor for centroid data."""
        return ResamplingMethod.NEAREST_NEIGHBOR

    def get_axis_type(self) -> AxisType:
        """Return reflector TOF axis for constant relative resolution."""
        return AxisType.REFLECTOR_TOF


class RapiflexDetector(InstrumentDetector):
    """Detector for Bruker Rapiflex MALDI-TOF data.

    Rapiflex profile data uses TIC-preserving resampling with equidistant
    (constant) axis spacing, matching SCiLS Lab convention.
    """

    @property
    def name(self) -> str:
        """Return detector name."""
        return "Rapiflex MALDI-TOF"

    def matches(self, characteristics: DataCharacteristics) -> bool:
        """Check if data came from Bruker flexImaging/Rapiflex.

        Both branches below are written by one producer, the Rapiflex
        metadata extractor, so a match here means the spectra arrive on the
        uniform-in-m/z grid ``RapiflexReader`` builds.

        Peak density is deliberately *not* consulted. It used to be: any
        profile data averaging more than
        ``Thresholds.PROFILE_PEAK_DENSITY`` points per spectrum was handed
        MALDI-TOF treatment regardless of what instrument produced it. That
        made the detector modality-blind -- a dense TOF-SIMS or Orbitrap
        profile acquisition would have been resampled by MALDI-TOF logic
        onto a MALDI-shaped axis, silently. SCiLS Lab does not guess
        modality either; its importer takes ``--project TIMSTOF|TOF|FT`` as
        an argument (2026b User Guide, p.81). Unknown-provenance profile
        data now falls through to :class:`DefaultDetector`, which bins
        counts rather than interpolating and so is safe for any modality.
        """
        # Direct Rapiflex format detection
        if characteristics.is_rapiflex_format:
            return True

        # Bruker MALDI-TOF detection
        return (
            characteristics.instrument_type == "MALDI-TOF"
            and characteristics.manufacturer == "Bruker"
        )

    def get_resampling_method(self) -> ResamplingMethod:
        """Return TIC-preserving for profile MALDI-TOF data."""
        return ResamplingMethod.TIC_PRESERVING

    def get_axis_type(self) -> AxisType:
        """Return constant axis matching SCiLS Lab convention."""
        return AxisType.CONSTANT

    @property
    def source_grid_law(self) -> Optional[AxisType]:
        """Report the constant law of the flexImaging source grid.

        ``RapiflexReader.get_common_mass_axis`` lays every spectrum out with
        ``np.linspace(mass_start, mass_end, n_points)``, so the source
        spacing is constant in m/z -- the same law as the ``constant`` target
        axis :meth:`get_axis_type` asks for. Source law equal to target law
        is what makes ``tic_preserving`` exact here, and it is the reason
        this is the only route on which the chain permits it.
        """
        return AxisType.CONSTANT


class TimsTOFDetector(InstrumentDetector):
    """Detector for Bruker timsTOF data.

    timsTOF produces centroid data with constant relative resolution,
    using nearest-neighbor resampling with reflector TOF axis spacing.
    """

    @property
    def name(self) -> str:
        """Return detector name."""
        return "timsTOF"

    def matches(self, characteristics: DataCharacteristics) -> bool:
        """Check if data is from timsTOF instrument."""
        return characteristics.is_timstof

    def get_resampling_method(self) -> ResamplingMethod:
        """Return nearest-neighbor for timsTOF centroid data."""
        return ResamplingMethod.NEAREST_NEIGHBOR

    def get_axis_type(self) -> AxisType:
        """Return reflector TOF axis for constant relative resolution."""
        return AxisType.REFLECTOR_TOF


class FTICRDetector(InstrumentDetector):
    """Detector for FT-ICR data.

    FT-ICR has quadratic mass axis spacing due to cyclotron frequency physics.
    Uses nearest-neighbor for centroid data.
    """

    @property
    def name(self) -> str:
        """Return detector name."""
        return "FT-ICR"

    def matches(self, characteristics: DataCharacteristics) -> bool:
        """Check if instrument type is FT-ICR."""
        return characteristics.instrument_type == "FT-ICR"

    def get_resampling_method(self) -> ResamplingMethod:
        """Return nearest-neighbor for FT-ICR data."""
        return ResamplingMethod.NEAREST_NEIGHBOR

    def get_axis_type(self) -> AxisType:
        """Return FTICR axis for quadratic mass scaling."""
        return AxisType.FTICR


class OrbitrapDetector(InstrumentDetector):
    """Detector for Orbitrap data.

    Orbitrap has 1/sqrt(m/z) resolution scaling.
    Uses nearest-neighbor for centroid data.
    """

    @property
    def name(self) -> str:
        """Return detector name."""
        return "Orbitrap"

    def matches(self, characteristics: DataCharacteristics) -> bool:
        """Check if instrument type is Orbitrap."""
        return characteristics.instrument_type == "Orbitrap"

    def get_resampling_method(self) -> ResamplingMethod:
        """Return nearest-neighbor for Orbitrap data."""
        return ResamplingMethod.NEAREST_NEIGHBOR

    def get_axis_type(self) -> AxisType:
        """Return Orbitrap axis for sqrt(m/z) resolution scaling."""
        return AxisType.ORBITRAP


class DefaultDetector(InstrumentDetector):
    """Fallback detector for unknown instruments.

    Uses conservative defaults: nearest-neighbor resampling with
    constant (equidistant) axis spacing.
    """

    @property
    def name(self) -> str:
        """Return detector name."""
        return "Unknown (default)"

    def matches(self, characteristics: DataCharacteristics) -> bool:
        """Always return True as fallback detector."""
        return True

    def get_resampling_method(self) -> ResamplingMethod:
        """Return nearest-neighbor as safe default."""
        return ResamplingMethod.NEAREST_NEIGHBOR

    def get_axis_type(self) -> AxisType:
        """Return constant axis as safe default."""
        return AxisType.CONSTANT


class InstrumentDetectorChain:
    """Chain of instrument detectors using Chain of Responsibility pattern.

    Iterates through detectors in priority order until one matches.
    The order is important - more specific detectors should come first.
    """

    def __init__(self, detectors: Optional[List[InstrumentDetector]] = None):
        """Initialize detector chain.

        Args:
            detectors: List of detectors in priority order.
                      If None, uses default detector chain.
        """
        if detectors is None:
            detectors = self._default_detectors()
        self.detectors = detectors

    @staticmethod
    def _default_detectors() -> List[InstrumentDetector]:
        """Create default detector chain in priority order."""
        return [
            # Specific instrument detectors first
            TimsTOFDetector(),
            RapiflexDetector(),
            FTICRDetector(),
            OrbitrapDetector(),
            # Generic spectrum type detector
            CentroidImzMLDetector(),
            # Fallback last
            DefaultDetector(),
        ]

    def detect(self, characteristics: DataCharacteristics) -> InstrumentDetector:
        """Find the first matching detector.

        Args:
            characteristics: Data characteristics to match

        Returns:
            The first detector that matches
        """
        for detector in self.detectors:
            if detector.matches(characteristics):
                logger.info(f"Detected instrument type: {detector.name}")
                return detector

        # Should never reach here due to DefaultDetector
        return DefaultDetector()

    def get_resampling_method(
        self, characteristics: DataCharacteristics
    ) -> ResamplingMethod:
        """Get resampling method for the detected instrument.

        ``TIC_PRESERVING`` additionally has to clear the matching-axis-law
        gate; see :meth:`_gate_tic_preserving`.

        Args:
            characteristics: Data characteristics to match

        Returns:
            Recommended resampling method for the instrument
        """
        detector = self.detect(characteristics)
        method = self._gate_tic_preserving(detector)
        logger.info(f"Selected resampling method: {method.name}")
        return method

    @staticmethod
    def _gate_tic_preserving(detector: InstrumentDetector) -> ResamplingMethod:
        """Allow ``TIC_PRESERVING`` only when source and target laws agree.

        SCiLS Lab applies TIC-preserving resampling to profile data only
        when all the mass axes being combined are of the same type, and
        linear interpolation otherwise (2026b User Guide, p.80). That is
        also precisely the condition under which Thyra's operator --
        interpolate, then rescale by one global factor -- is exact. See
        :attr:`InstrumentDetector.source_grid_law` for the error off the
        diagonal.

        A detector that has not declared its source grid law does not clear
        the gate. Auto-selection therefore cannot reach the interpolating
        path on data whose acquisition Thyra has not actually identified.

        Args:
            detector: The detector that matched.

        Returns:
            The detector's method, or ``NEAREST_NEIGHBOR`` in its place.
        """
        method = detector.get_resampling_method()
        if method is not ResamplingMethod.TIC_PRESERVING:
            return method

        axis_type = detector.get_axis_type()
        source_law = detector.source_grid_law
        if source_law is axis_type:
            return method

        logger.info(
            "%s asks for TIC_PRESERVING onto a %s axis, but the source grid "
            "law is %s. TIC-preserving resampling is exact only when the two "
            "match, so NEAREST_NEIGHBOR is used instead.",
            detector.name,
            axis_type.name,
            "unknown" if source_law is None else source_law.name,
        )
        return ResamplingMethod.NEAREST_NEIGHBOR

    def get_axis_type(self, characteristics: DataCharacteristics) -> AxisType:
        """Get axis type for the detected instrument.

        Args:
            characteristics: Data characteristics to match

        Returns:
            Recommended axis type for the instrument
        """
        detector = self.detect(characteristics)
        axis_type = detector.get_axis_type()
        logger.info(f"Selected axis type: {axis_type.name}")
        return axis_type
