# tests/unit/metadata/core/test_metadata_types.py
import pytest

from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata


class TestEssentialMetadata:
    """Test EssentialMetadata dataclass functionality."""

    def test_creation_with_required_fields(self):
        """Test creation with all required fields."""
        metadata = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )

        assert metadata.dimensions == (10, 20, 1)
        assert metadata.coordinate_bounds == (0.0, 100.0, 0.0, 200.0)
        assert metadata.mass_range == (100.0, 1000.0)
        assert metadata.pixel_size == (25.0, 25.0)
        assert metadata.n_spectra == 200
        assert metadata.total_peaks == 400000
        assert metadata.source_path == "/path/to/data.imzML"

    def test_creation_with_optional_none_pixel_size(self):
        """Test creation with None pixel size."""
        metadata = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=None,
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )

        assert metadata.pixel_size is None

    def test_has_pixel_size_property(self):
        """Test has_pixel_size property."""
        # With pixel size
        metadata_with_size = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )
        assert metadata_with_size.has_pixel_size is True

        # Without pixel size
        metadata_without_size = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=None,
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )
        assert metadata_without_size.has_pixel_size is False

    def test_is_3d_property(self):
        """Test is_3d property."""
        # 2D dataset
        metadata_2d = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )
        assert metadata_2d.is_3d is False

        # 3D dataset
        metadata_3d = EssentialMetadata(
            dimensions=(10, 20, 5),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=1000,
            total_peaks=2000000,
            source_path="/path/to/data.imzML",
        )
        assert metadata_3d.is_3d is True

    def test_immutability(self):
        """Test that EssentialMetadata is frozen/immutable."""
        metadata = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )

        # Should not be able to modify fields
        with pytest.raises(AttributeError):
            metadata.dimensions = (5, 10, 1)

        with pytest.raises(AttributeError):
            metadata.n_spectra = 100


class TestComprehensiveMetadata:
    """Test ComprehensiveMetadata dataclass functionality."""

    def test_creation_with_essential_metadata(self):
        """Test creation with essential metadata."""
        essential = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )

        comprehensive = ComprehensiveMetadata(
            essential=essential,
            format_specific={"imzml_version": "1.1.0", "scan_mode": "profile"},
            acquisition_params={"laser_power": 50, "scan_rate": 1000},
            instrument_info={"vendor": "Bruker", "model": "timsTOF"},
            raw_metadata={"original_key": "original_value"},
        )

        assert comprehensive.essential is essential
        assert comprehensive.format_specific["imzml_version"] == "1.1.0"
        assert comprehensive.acquisition_params["laser_power"] == 50
        assert comprehensive.instrument_info["vendor"] == "Bruker"
        assert comprehensive.raw_metadata["original_key"] == "original_value"

    def test_convenience_properties(self):
        """Test convenience properties that delegate to essential metadata."""
        essential = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )

        comprehensive = ComprehensiveMetadata(
            essential=essential,
            format_specific={},
            acquisition_params={},
            instrument_info={},
            raw_metadata={},
        )

        # Test convenience properties
        assert comprehensive.dimensions == (10, 20, 1)
        assert comprehensive.pixel_size == (25.0, 25.0)
        assert comprehensive.coordinate_bounds == (0.0, 100.0, 0.0, 200.0)

    def test_comprehensive_with_none_pixel_size(self):
        """Test comprehensive metadata with None pixel size."""
        essential = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=None,
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )

        comprehensive = ComprehensiveMetadata(
            essential=essential,
            format_specific={},
            acquisition_params={},
            instrument_info={},
            raw_metadata={},
        )

        assert comprehensive.pixel_size is None

    def test_empty_dictionaries(self):
        """Test creation with empty dictionaries."""
        essential = EssentialMetadata(
            dimensions=(10, 20, 1),
            coordinate_bounds=(0.0, 100.0, 0.0, 200.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(25.0, 25.0),
            n_spectra=200,
            total_peaks=400000,
            source_path="/path/to/data.imzML",
        )

        comprehensive = ComprehensiveMetadata(
            essential=essential,
            format_specific={},
            acquisition_params={},
            instrument_info={},
            raw_metadata={},
        )

        assert comprehensive.format_specific == {}
        assert comprehensive.acquisition_params == {}
        assert comprehensive.instrument_info == {}
        assert comprehensive.raw_metadata == {}
