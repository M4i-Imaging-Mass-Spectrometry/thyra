"""
Tests for ResamplingDecisionTree with Strategy pattern for instrument detection.
"""

import pytest

from thyra.resampling.constants import SpectrumType
from thyra.resampling.data_characteristics import DataCharacteristics
from thyra.resampling.decision_tree import ResamplingDecisionTree
from thyra.resampling.instrument_detectors import (
    CentroidImzMLDetector,
    DefaultDetector,
    FTICRDetector,
    InstrumentDetectorChain,
    OrbitrapDetector,
    RapiflexDetector,
    TimsTOFDetector,
)
from thyra.resampling.types import AxisType, ResamplingMethod


class TestResamplingDecisionTree:
    """Test ResamplingDecisionTree with Strategy pattern."""

    def setup_method(self):
        """Setup test fixtures."""
        self.tree = ResamplingDecisionTree()

    def test_no_metadata_raises_error(self):
        """Test that None metadata raises NotImplementedError."""
        with pytest.raises(NotImplementedError) as exc_info:
            self.tree.select_strategy(None)

        assert "metadata" in str(exc_info.value).lower()

    def test_empty_metadata_uses_default(self):
        """Test that empty metadata uses default detector."""
        # With the Strategy pattern, empty metadata falls through to DefaultDetector
        method = self.tree.select_strategy({})
        assert method == ResamplingMethod.NEAREST_NEIGHBOR

    def test_timstof_detection_from_bruker_metadata(self):
        """Test timsTOF detection from Bruker GlobalMetadata."""
        bruker_timstof_metadata = {
            "GlobalMetadata": {"InstrumentName": "timsTOF Maldi 2"}
        }
        method = self.tree.select_strategy(bruker_timstof_metadata)
        assert method == ResamplingMethod.NEAREST_NEIGHBOR

        axis_type = self.tree.select_axis_type(bruker_timstof_metadata)
        assert axis_type == AxisType.REFLECTOR_TOF

    def test_centroid_spectrum_detection(self):
        """Test centroid spectrum detection from essential_metadata."""
        metadata = {"essential_metadata": {"spectrum_type": SpectrumType.CENTROID}}
        method = self.tree.select_strategy(metadata)
        assert method == ResamplingMethod.NEAREST_NEIGHBOR

    def test_centroid_spectrum_detection_in_instrument_info(self):
        """Test centroid spectrum detection via instrument info."""
        metadata = {
            "essential_metadata": {"spectrum_type": SpectrumType.CENTROID},
            "instrument_info": {"instrument_type": "Q-TOF"},
        }
        method = self.tree.select_strategy(metadata)
        assert method == ResamplingMethod.NEAREST_NEIGHBOR

    def test_profile_spectrum_with_high_density_uses_tic_preserving(self):
        """Test that profile spectrum with high peak density uses TIC preserving."""
        metadata = {
            "essential_metadata": {
                "spectrum_type": SpectrumType.PROFILE,
                "total_peaks": 10000000,  # 10 million peaks
                "n_spectra": 1000,  # 10000 peaks per spectrum
            }
        }
        method = self.tree.select_strategy(metadata)
        assert method == ResamplingMethod.TIC_PRESERVING

    def test_rapiflex_format_detection(self):
        """Test Rapiflex format detection."""
        metadata = {"format_specific": {"format": "Rapiflex"}}
        method = self.tree.select_strategy(metadata)
        assert method == ResamplingMethod.TIC_PRESERVING

        axis_type = self.tree.select_axis_type(metadata)
        assert axis_type == AxisType.CONSTANT

    def test_bruker_maldi_tof_detection(self):
        """Test Bruker MALDI-TOF detection."""
        metadata = {
            "instrument_info": {
                "instrument_type": "MALDI-TOF",
                "manufacturer": "Bruker",
            }
        }
        method = self.tree.select_strategy(metadata)
        assert method == ResamplingMethod.TIC_PRESERVING

    def test_axis_type_selection_no_metadata(self):
        """Test axis type selection with no metadata uses default."""
        axis_type = self.tree.select_axis_type(None)
        assert axis_type == AxisType.CONSTANT


class TestInstrumentDetectorChain:
    """Test InstrumentDetectorChain behavior."""

    def setup_method(self):
        """Setup test fixtures."""
        self.chain = InstrumentDetectorChain()

    def test_default_chain_order(self):
        """Test that default chain has correct order."""
        detector_types = [type(d).__name__ for d in self.chain.detectors]
        assert detector_types == [
            "TimsTOFDetector",
            "RapiflexDetector",
            "FTICRDetector",
            "OrbitrapDetector",
            "CentroidImzMLDetector",
            "DefaultDetector",
        ]

    def test_timstof_takes_priority(self):
        """Test that TimsTOF detector takes priority over others."""
        characteristics = DataCharacteristics(is_timstof=True)
        detector = self.chain.detect(characteristics)
        assert isinstance(detector, TimsTOFDetector)

    def test_rapiflex_detected_before_centroid(self):
        """Test Rapiflex detection takes priority over generic centroid."""
        characteristics = DataCharacteristics(
            is_rapiflex_format=True,
            spectrum_type=SpectrumType.CENTROID,
        )
        detector = self.chain.detect(characteristics)
        assert isinstance(detector, RapiflexDetector)

    def test_fallback_to_default(self):
        """Test fallback to DefaultDetector when nothing matches."""
        characteristics = DataCharacteristics()
        detector = self.chain.detect(characteristics)
        assert isinstance(detector, DefaultDetector)


class TestInstrumentDetectors:
    """Test individual instrument detectors."""

    def test_timstof_detector(self):
        """Test TimsTOF detector matching."""
        detector = TimsTOFDetector()
        assert detector.name == "timsTOF"

        # Should match when is_timstof is True
        chars = DataCharacteristics(is_timstof=True)
        assert detector.matches(chars)

        # Should not match when is_timstof is False
        chars = DataCharacteristics(is_timstof=False)
        assert not detector.matches(chars)

        assert detector.get_resampling_method() == ResamplingMethod.NEAREST_NEIGHBOR
        assert detector.get_axis_type() == AxisType.REFLECTOR_TOF

    def test_rapiflex_detector(self):
        """Test Rapiflex detector matching."""
        detector = RapiflexDetector()
        assert detector.name == "Rapiflex MALDI-TOF"

        # Should match Rapiflex format
        chars = DataCharacteristics(is_rapiflex_format=True)
        assert detector.matches(chars)

        # Should match Bruker MALDI-TOF
        chars = DataCharacteristics(
            instrument_type="MALDI-TOF",
            manufacturer="Bruker",
        )
        assert detector.matches(chars)

        # Should match high-density profile data
        chars = DataCharacteristics(
            spectrum_type=SpectrumType.PROFILE,
            total_peaks=10000000,
            n_spectra=1000,
        )
        assert detector.matches(chars)

        assert detector.get_resampling_method() == ResamplingMethod.TIC_PRESERVING
        assert detector.get_axis_type() == AxisType.CONSTANT

    def test_centroid_imzml_detector(self):
        """Test CentroidImzML detector matching."""
        detector = CentroidImzMLDetector()

        # Should match centroid data
        chars = DataCharacteristics(spectrum_type=SpectrumType.CENTROID)
        assert detector.matches(chars)

        # Should not match profile data
        chars = DataCharacteristics(spectrum_type=SpectrumType.PROFILE)
        assert not detector.matches(chars)

        assert detector.get_resampling_method() == ResamplingMethod.NEAREST_NEIGHBOR
        assert detector.get_axis_type() == AxisType.REFLECTOR_TOF

    def test_fticr_detector(self):
        """Test FTICR detector matching."""
        detector = FTICRDetector()

        chars = DataCharacteristics(instrument_type="FT-ICR")
        assert detector.matches(chars)

        chars = DataCharacteristics(instrument_type="Orbitrap")
        assert not detector.matches(chars)

        assert detector.get_axis_type() == AxisType.FTICR

    def test_orbitrap_detector(self):
        """Test Orbitrap detector matching."""
        detector = OrbitrapDetector()

        chars = DataCharacteristics(instrument_type="Orbitrap")
        assert detector.matches(chars)

        chars = DataCharacteristics(instrument_type="FT-ICR")
        assert not detector.matches(chars)

        assert detector.get_axis_type() == AxisType.ORBITRAP

    def test_default_detector(self):
        """Test DefaultDetector always matches."""
        detector = DefaultDetector()
        assert detector.name == "Unknown (default)"

        # Should always match
        chars = DataCharacteristics()
        assert detector.matches(chars)

        chars = DataCharacteristics(
            instrument_type="Unknown Instrument",
            spectrum_type="unknown",
        )
        assert detector.matches(chars)

        assert detector.get_resampling_method() == ResamplingMethod.NEAREST_NEIGHBOR
        assert detector.get_axis_type() == AxisType.CONSTANT


class TestDataCharacteristics:
    """Test DataCharacteristics dataclass."""

    def test_from_metadata_essential(self):
        """Test creating from essential metadata."""
        metadata = {
            "essential_metadata": {
                "spectrum_type": SpectrumType.CENTROID,
                "total_peaks": 1000000,
                "n_spectra": 500,
            }
        }
        chars = DataCharacteristics.from_metadata(metadata)

        assert chars.spectrum_type == SpectrumType.CENTROID
        assert chars.total_peaks == 1000000
        assert chars.n_spectra == 500
        assert chars.is_centroid_data

    def test_from_metadata_instrument_info(self):
        """Test extracting instrument info from metadata."""
        metadata = {
            "instrument_info": {
                "instrument_type": "MALDI-TOF",
                "manufacturer": "Bruker",
            }
        }
        chars = DataCharacteristics.from_metadata(metadata)

        assert chars.instrument_type == "MALDI-TOF"
        assert chars.manufacturer == "Bruker"

    def test_from_metadata_global_metadata(self):
        """Test extracting from GlobalMetadata."""
        metadata = {"GlobalMetadata": {"InstrumentName": "timsTOF Maldi 2"}}
        chars = DataCharacteristics.from_metadata(metadata)

        assert chars.instrument_name == "timsTOF Maldi 2"
        assert chars.is_timstof

    def test_is_high_density_profile(self):
        """Test high density profile detection."""
        # High density profile (>5000 peaks per spectrum)
        chars = DataCharacteristics(
            spectrum_type=SpectrumType.PROFILE,
            total_peaks=10000000,
            n_spectra=1000,
        )
        assert chars.is_high_density_profile
        assert chars.avg_peaks_per_spectrum == 10000.0

        # Low density profile
        chars = DataCharacteristics(
            spectrum_type=SpectrumType.PROFILE,
            total_peaks=1000,
            n_spectra=1000,
        )
        assert not chars.is_high_density_profile

        # Centroid data should not be high density profile
        chars = DataCharacteristics(
            spectrum_type=SpectrumType.CENTROID,
            total_peaks=10000000,
            n_spectra=1000,
        )
        assert not chars.is_high_density_profile

    def test_needs_resampling(self):
        """Test needs_resampling property."""
        # Continuous data doesn't need resampling
        chars = DataCharacteristics(has_shared_mass_axis=True)
        assert not chars.needs_resampling

        # Processed data needs resampling
        chars = DataCharacteristics(has_shared_mass_axis=False)
        assert chars.needs_resampling

    def test_is_maldi_tof(self):
        """Test MALDI-TOF detection."""
        # Rapiflex format
        chars = DataCharacteristics(is_rapiflex_format=True)
        assert chars.is_maldi_tof

        # Explicit MALDI-TOF type
        chars = DataCharacteristics(instrument_type="MALDI-TOF")
        assert chars.is_maldi_tof

        # Bruker with high-density profile
        chars = DataCharacteristics(
            manufacturer="Bruker",
            spectrum_type=SpectrumType.PROFILE,
            total_peaks=10000000,
            n_spectra=1000,
        )
        assert chars.is_maldi_tof
