"""Data characteristics for resampling decisions."""

from dataclasses import dataclass
from typing import Any, Dict, Optional

from .constants import SpectrumType, Thresholds


@dataclass
class DataCharacteristics:
    """Detected characteristics of MSI data for resampling decisions.

    This dataclass consolidates all the information needed to make
    intelligent resampling decisions, extracted from metadata.
    """

    # Core data properties
    has_shared_mass_axis: bool = False  # True for continuous, False for processed
    spectrum_type: Optional[str] = None  # "centroid spectrum" or "profile spectrum"

    # Instrument identification
    instrument_type: Optional[str] = None  # "MALDI-TOF", "timsTOF", "FT-ICR", etc.
    instrument_name: Optional[str] = None  # Specific model name
    manufacturer: Optional[str] = None  # "Bruker", etc.

    # Software provenance
    software_source: Optional[str] = None  # "SCiLS Lab", etc.

    # Data statistics
    total_peaks: Optional[int] = None
    n_spectra: Optional[int] = None

    # Format-specific flags
    is_rapiflex_format: bool = False
    is_timstof: bool = False

    @property
    def needs_resampling(self) -> bool:
        """Determine if data needs mass axis alignment.

        Continuous data (shared mass axis) doesn't need resampling for alignment.
        Processed data (different m/z per spectrum) needs resampling.
        """
        return not self.has_shared_mass_axis

    @property
    def is_profile_data(self) -> bool:
        """Check if this is profile (continuous signal) data."""
        return self.spectrum_type == SpectrumType.PROFILE

    @property
    def is_centroid_data(self) -> bool:
        """Check if this is centroid (discrete peaks) data."""
        return self.spectrum_type == SpectrumType.CENTROID

    @property
    def avg_peaks_per_spectrum(self) -> Optional[float]:
        """Calculate average peaks per spectrum."""
        if self.n_spectra and self.n_spectra > 0 and self.total_peaks:
            return self.total_peaks / self.n_spectra
        return None

    @property
    def is_high_density_profile(self) -> bool:
        """Check if this appears to be high-density profile data.

        Profile data typically has >5000 points per spectrum,
        indicating continuous signal rather than centroid peaks.
        """
        avg = self.avg_peaks_per_spectrum
        return (
            self.is_profile_data
            and avg is not None
            and avg > Thresholds.PROFILE_PEAK_DENSITY
        )

    @property
    def is_maldi_tof(self) -> bool:
        """Check if this is MALDI-TOF data."""
        return (
            self.is_rapiflex_format
            or self.instrument_type == "MALDI-TOF"
            or (self.manufacturer == "Bruker" and self.is_high_density_profile)
        )

    @classmethod
    def from_metadata(cls, metadata: Dict[str, Any]) -> "DataCharacteristics":
        """Create DataCharacteristics from metadata dictionary.

        Args:
            metadata: Metadata dictionary with nested essential_metadata,
                     instrument_info, format_specific, etc.

        Returns:
            DataCharacteristics instance populated from metadata.
        """
        essential = metadata.get("essential_metadata", {})
        instrument_info = metadata.get("instrument_info", {})
        format_specific = metadata.get("format_specific", {})
        global_meta = metadata.get("GlobalMetadata", {})

        # Extract values with safe defaults
        spectrum_type = essential.get("spectrum_type") if essential else None
        total_peaks = essential.get("total_peaks") if essential else None
        n_spectra = essential.get("n_spectra") if essential else None

        # Instrument info
        instrument_type = (
            instrument_info.get("instrument_type") if instrument_info else None
        )
        manufacturer = instrument_info.get("manufacturer") if instrument_info else None
        instrument_name = global_meta.get("InstrumentName") if global_meta else None

        # Format detection
        is_rapiflex = (
            format_specific.get("format") == "Rapiflex" if format_specific else False
        )
        # Substring match: Bruker names timsTOF variants several ways
        # ("timsTOF Maldi 2", "timsTOF Pro 2", "timsTOF fleX MALDI-2",
        # "timsTOF SCP", etc.).  An exact match was missing every
        # variant except the original Maldi 2, so the InstrumentDetectorChain
        # silently fell through TimsTOFDetector and ended up at the
        # DefaultDetector returning CONSTANT.  This match is permissive
        # but case-sensitive on the brand "timsTOF" / "timstof".
        is_timstof = bool(instrument_name and "timstof" in instrument_name.lower())

        return cls(
            spectrum_type=spectrum_type,
            total_peaks=total_peaks,
            n_spectra=n_spectra,
            instrument_type=instrument_type,
            instrument_name=instrument_name,
            manufacturer=manufacturer,
            is_rapiflex_format=is_rapiflex,
            is_timstof=is_timstof,
        )
