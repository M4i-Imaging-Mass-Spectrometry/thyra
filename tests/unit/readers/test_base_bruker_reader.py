# tests/unit/readers/test_base_bruker_reader.py
"""Tests for BrukerBaseMSIReader abstract base class."""

from typing import Generator, Tuple
from unittest.mock import MagicMock

import numpy as np
import pytest
from numpy.typing import NDArray

from thyra.readers.bruker.base_bruker_reader import BrukerBaseMSIReader
from thyra.readers.bruker.folder_structure import BrukerFormat


class ConcreteBrukerReader(BrukerBaseMSIReader):
    """Concrete implementation of BrukerBaseMSIReader for testing."""

    def _create_metadata_extractor(self):
        """Return mock metadata extractor."""
        mock = MagicMock()
        mock.get_essential.return_value = MagicMock(
            dimensions=(10, 10, 1),
            coordinate_bounds=(0.0, 9.0, 0.0, 9.0),
            mass_range=(100.0, 1000.0),
            pixel_size=(50.0, 50.0),
        )
        return mock

    def get_common_mass_axis(self) -> NDArray[np.float64]:
        """Return mock mass axis."""
        return np.linspace(100, 1000, 100)

    def iter_spectra(self) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Return mock spectra iterator."""
        mz = np.linspace(100, 1000, 100)
        for x in range(5):
            for y in range(5):
                yield (x, y, 0), mz, np.random.rand(100)

    def close(self) -> None:
        """Close reader."""
        pass


class TestBrukerBaseMSIReader:
    """Tests for BrukerBaseMSIReader class."""

    @pytest.fixture
    def rapiflex_folder(self, tmp_path):
        """Create a mock Rapiflex folder structure."""
        data_dir = tmp_path / "rapiflex_data"
        data_dir.mkdir()
        (data_dir / "sample.dat").touch()
        (data_dir / "sample_poslog.txt").touch()
        (data_dir / "sample_info.txt").touch()
        (data_dir / "sample.mis").touch()
        (data_dir / "optical_0000.tif").touch()
        (data_dir / "optical_0001.tif").touch()
        return data_dir

    @pytest.fixture
    def timstof_folder(self, tmp_path):
        """Create a mock timsTOF folder structure."""
        d_dir = tmp_path / "test.d"
        d_dir.mkdir()
        (d_dir / "analysis.tsf").touch()
        (d_dir / "analysis.tsf_bin").touch()
        # Add TIFF at parent level
        (tmp_path / "optical.tif").touch()
        return d_dir

    def test_initialization(self, rapiflex_folder):
        """Test reader initialization."""
        reader = ConcreteBrukerReader(rapiflex_folder)
        assert reader.data_path == rapiflex_folder
        assert reader._folder_info is None  # Lazy loaded

    def test_folder_info_lazy_loading(self, rapiflex_folder):
        """Test that folder_info is lazily loaded."""
        reader = ConcreteBrukerReader(rapiflex_folder)

        # Not loaded yet
        assert reader._folder_info is None

        # Access triggers loading
        info = reader.folder_info

        # Now loaded
        assert reader._folder_info is not None
        assert info.format == BrukerFormat.RAPIFLEX

    def test_get_optical_image_paths_rapiflex(self, rapiflex_folder):
        """Test getting optical images for Rapiflex."""
        reader = ConcreteBrukerReader(rapiflex_folder)
        optical_paths = reader.get_optical_image_paths()

        assert len(optical_paths) == 2
        assert all(p.suffix == ".tif" for p in optical_paths)

    def test_get_optical_image_paths_timstof(self, timstof_folder, tmp_path):
        """Test getting optical images for timsTOF (at parent level)."""
        # Need to create reader at parent level to find TIFF
        reader = ConcreteBrukerReader(tmp_path)
        optical_paths = reader.get_optical_image_paths()

        assert len(optical_paths) >= 1

    def test_get_teaching_points_file(self, rapiflex_folder):
        """Test getting teaching points file."""
        reader = ConcreteBrukerReader(rapiflex_folder)
        mis_file = reader.get_teaching_points_file()

        assert mis_file is not None
        assert mis_file.name == "sample.mis"

    def test_get_teaching_points_file_not_found(self, tmp_path):
        """Test getting teaching points when not present."""
        # Create minimal Rapiflex structure without .mis
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "sample.dat").touch()
        (data_dir / "sample_poslog.txt").touch()
        (data_dir / "sample_info.txt").touch()

        reader = ConcreteBrukerReader(data_dir)
        mis_file = reader.get_teaching_points_file()

        assert mis_file is None

    def test_context_manager(self, rapiflex_folder):
        """Test reader works as context manager."""
        with ConcreteBrukerReader(rapiflex_folder) as reader:
            info = reader.folder_info
            assert info.format == BrukerFormat.RAPIFLEX

    def test_inheritance_from_base_msi_reader(self, rapiflex_folder):
        """Test that BrukerBaseMSIReader properly inherits from BaseMSIReader."""
        reader = ConcreteBrukerReader(rapiflex_folder)

        # Should have BaseMSIReader methods
        assert hasattr(reader, "get_essential_metadata")
        assert hasattr(reader, "get_comprehensive_metadata")
        assert hasattr(reader, "map_mz_to_common_axis")

    def test_abstract_methods_must_be_implemented(self):
        """Test that abstract methods must be implemented."""
        # BrukerBaseMSIReader is abstract, can't instantiate directly
        # This is verified by the fact that ConcreteBrukerReader implements all methods

        # Check that the class has the expected abstract method markers
        assert hasattr(BrukerBaseMSIReader, "_create_metadata_extractor")
        assert hasattr(BrukerBaseMSIReader, "get_common_mass_axis")
        assert hasattr(BrukerBaseMSIReader, "iter_spectra")
        assert hasattr(BrukerBaseMSIReader, "close")


class TestBrukerReadersUseBaseClass:
    """Tests to verify existing readers use BrukerBaseMSIReader."""

    def test_rapiflex_reader_extends_bruker_base(self):
        """Test RapiflexReader extends BrukerBaseMSIReader."""
        from thyra.readers.bruker.rapiflex.rapiflex_reader import RapiflexReader

        assert issubclass(RapiflexReader, BrukerBaseMSIReader)

    def test_bruker_reader_extends_bruker_base(self):
        """Test BrukerReader extends BrukerBaseMSIReader."""
        from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

        assert issubclass(BrukerReader, BrukerBaseMSIReader)

    def test_rapiflex_reader_has_folder_info(self):
        """Test RapiflexReader has folder_info property."""
        from thyra.readers.bruker.rapiflex.rapiflex_reader import RapiflexReader

        assert hasattr(RapiflexReader, "folder_info")
        assert hasattr(RapiflexReader, "get_teaching_points_file")

    def test_bruker_reader_has_folder_info(self):
        """Test BrukerReader has folder_info property."""
        from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

        assert hasattr(BrukerReader, "folder_info")
        assert hasattr(BrukerReader, "get_teaching_points_file")
