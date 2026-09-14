# tests/unit/readers/test_bruker_folder_structure.py
"""Tests for BrukerFolderStructure abstraction."""

from pathlib import Path

import pytest

from thyra.errors import ConversionRefused
from thyra.readers.bruker.folder_structure import (
    BrukerFolderInfo,
    BrukerFolderStructure,
    BrukerFormat,
)


def _rapiflex_folder(tmp_path: Path, name: str = "data") -> Path:
    """The three files that make a folder look like a Rapiflex acquisition."""
    data_dir = tmp_path / name
    data_dir.mkdir()
    (data_dir / "sample.dat").touch()
    (data_dir / "sample_poslog.txt").touch()
    (data_dir / "sample_info.txt").touch()
    return data_dir


def _write_mis(path: Path, image_file=None, original_image=None) -> Path:
    """A minimal FlexImaging .mis naming (or not naming) an optical image."""
    body = ['<?xml version="1.0"?>', '<ImagingSequence flexImagingVersion="5.0.0">']
    if image_file is not None:
        body.append(f"  <ImageFile>{image_file}</ImageFile>")
    if original_image is not None:
        body.append(f"  <OriginalImage>{original_image}</OriginalImage>")
    body.append("</ImagingSequence>")
    path.write_text("\n".join(body) + "\n", encoding="utf-8")
    return path


class TestBrukerFormat:
    """Tests for BrukerFormat enum."""

    def test_format_values(self):
        """Test that format enum has expected values."""
        assert BrukerFormat.TIMSTOF.value == "timstof"
        assert BrukerFormat.RAPIFLEX.value == "rapiflex"
        assert BrukerFormat.UNKNOWN.value == "unknown"


class TestBrukerFolderStructure:
    """Tests for BrukerFolderStructure class."""

    def test_detect_timstof_tsf(self, tmp_path):
        """Test detection of timsTOF TSF format."""
        # Create .d directory with analysis.tsf
        d_dir = tmp_path / "test.d"
        d_dir.mkdir()
        (d_dir / "analysis.tsf").touch()
        (d_dir / "analysis.tsf_bin").touch()

        folder = BrukerFolderStructure(d_dir)
        info = folder.analyze()

        assert info.format == BrukerFormat.TIMSTOF
        assert info.data_path == d_dir

    def test_detect_timstof_tdf(self, tmp_path):
        """Test detection of timsTOF TDF format."""
        # Create .d directory with analysis.tdf
        d_dir = tmp_path / "test.d"
        d_dir.mkdir()
        (d_dir / "analysis.tdf").touch()
        (d_dir / "analysis.tdf_bin").touch()

        folder = BrukerFolderStructure(d_dir)
        info = folder.analyze()

        assert info.format == BrukerFormat.TIMSTOF
        assert info.data_path == d_dir

    def test_detect_rapiflex(self, tmp_path):
        """Test detection of Rapiflex format."""
        # Create Rapiflex folder structure
        data_dir = tmp_path / "rapiflex_data"
        data_dir.mkdir()
        (data_dir / "sample.dat").touch()
        (data_dir / "sample_poslog.txt").touch()
        (data_dir / "sample_info.txt").touch()

        folder = BrukerFolderStructure(data_dir)
        info = folder.analyze()

        assert info.format == BrukerFormat.RAPIFLEX
        assert info.data_path == data_dir

    def test_detect_unknown_format(self, tmp_path):
        """Test detection returns UNKNOWN for unrecognized format."""
        # Create empty directory
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()

        folder = BrukerFolderStructure(empty_dir)
        info = folder.analyze()

        assert info.format == BrukerFormat.UNKNOWN

    def test_detect_timstof_from_parent(self, tmp_path):
        """Test detection of timsTOF from parent folder."""
        # Create parent folder containing .d subfolder
        parent = tmp_path / "experiment"
        parent.mkdir()
        d_dir = parent / "data.d"
        d_dir.mkdir()
        (d_dir / "analysis.tsf").touch()

        folder = BrukerFolderStructure(parent)
        info = folder.analyze()

        assert info.format == BrukerFormat.TIMSTOF
        assert info.data_path == d_dir

    def test_detect_rapiflex_from_parent(self, tmp_path):
        """Test detection of Rapiflex from parent folder."""
        # Create parent folder containing Rapiflex subfolder
        parent = tmp_path / "experiment"
        parent.mkdir()
        data_dir = parent / "data"
        data_dir.mkdir()
        (data_dir / "sample.dat").touch()
        (data_dir / "sample_poslog.txt").touch()
        (data_dir / "sample_info.txt").touch()

        folder = BrukerFolderStructure(parent)
        info = folder.analyze()

        assert info.format == BrukerFormat.RAPIFLEX
        assert info.data_path == data_dir

    def test_find_optical_images(self, tmp_path):
        """Test finding optical TIFF images."""
        # Create Rapiflex structure with TIFF files
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "sample.dat").touch()
        (data_dir / "sample_poslog.txt").touch()
        (data_dir / "sample_info.txt").touch()
        (data_dir / "optical_0000.tif").touch()
        (data_dir / "optical_0001.tiff").touch()

        folder = BrukerFolderStructure(data_dir)
        info = folder.analyze()

        assert len(info.optical_images) == 2
        assert any("optical_0000.tif" in str(p) for p in info.optical_images)
        assert any("optical_0001.tiff" in str(p) for p in info.optical_images)

    def test_find_optical_images_at_parent(self, tmp_path):
        """Test finding optical images in parent folder."""
        # Create structure where TIFFs are at parent level
        parent = tmp_path / "experiment"
        parent.mkdir()
        (parent / "optical.tif").touch()

        d_dir = parent / "data.d"
        d_dir.mkdir()
        (d_dir / "analysis.tsf").touch()

        folder = BrukerFolderStructure(parent)
        info = folder.analyze()

        assert len(info.optical_images) == 1
        assert "optical.tif" in str(info.optical_images[0])

    def test_find_optical_images_non_tiff(self, tmp_path):
        """JPEG, PNG and BMP are optical images too, not just TIFF."""
        data_dir = _rapiflex_folder(tmp_path)
        for name in ("scan_0000.jpg", "scan_0001.jpeg", "map.png", "slide.bmp"):
            (data_dir / name).touch()
        (data_dir / "notes.txt").touch()

        info = BrukerFolderStructure(data_dir).analyze()

        assert sorted(p.name for p in info.optical_images) == [
            "map.png",
            "scan_0000.jpg",
            "scan_0001.jpeg",
            "slide.bmp",
        ]

    def test_find_optical_images_ignores_suffix_case(self, tmp_path):
        """An upper-case suffix is found, and found exactly once."""
        data_dir = _rapiflex_folder(tmp_path)
        (data_dir / "scan.JPG").touch()
        (data_dir / "overview.TIFF").touch()

        info = BrukerFolderStructure(data_dir).analyze()

        assert sorted(p.name for p in info.optical_images) == [
            "overview.TIFF",
            "scan.JPG",
        ]

    def test_mis_named_image_wins_over_the_glob(self, tmp_path):
        """A folder with several JPEGs: the .mis says which one is the image."""
        data_dir = _rapiflex_folder(tmp_path)
        for name in ("slide_overview.jpg", "sample_0000.jpg", "sample_0001.jpg"):
            (data_dir / name).touch()
        _write_mis(data_dir / "sample.mis", image_file="sample_0000.jpg")

        info = BrukerFolderStructure(data_dir).analyze()

        assert info.primary_optical_image == data_dir / "sample_0000.jpg"
        # First, so the converter sizes the others against it.
        assert info.optical_images[0] == data_dir / "sample_0000.jpg"
        # And still listed exactly once, not twice.
        assert sorted(p.name for p in info.optical_images) == [
            "sample_0000.jpg",
            "sample_0001.jpg",
            "slide_overview.jpg",
        ]

    def test_original_image_is_provenance_only(self, tmp_path):
        """<OriginalImage> is an absolute path from another machine: never opened."""
        data_dir = _rapiflex_folder(tmp_path)
        (data_dir / "sample_0000.jpg").touch()
        _write_mis(
            data_dir / "sample.mis",
            image_file="sample_0000.jpg",
            original_image="Z:\\FlexImaging\\Runs\\2026\\sample_0000.jpg",
        )

        info = BrukerFolderStructure(data_dir).analyze()

        # Resolved from <ImageFile>, against the acquisition folder.
        assert info.primary_optical_image == data_dir / "sample_0000.jpg"
        assert all(p.is_file() for p in info.optical_images)

    def test_no_image_file_element_means_no_primary(self, tmp_path):
        """<OriginalImage> alone names nothing Thyra will resolve."""
        data_dir = _rapiflex_folder(tmp_path)
        (data_dir / "sample_0000.jpg").touch()
        _write_mis(
            data_dir / "sample.mis",
            image_file=None,
            original_image="Z:\\FlexImaging\\Runs\\2026\\sample_0000.jpg",
        )

        info = BrukerFolderStructure(data_dir).analyze()

        assert info.primary_optical_image is None
        # The glob still finds what is actually there.
        assert [p.name for p in info.optical_images] == ["sample_0000.jpg"]

    def test_mis_names_an_image_that_is_not_there(self, tmp_path, thyra_logs):
        """Named but missing: say so, and fall back to whatever is present."""
        import logging

        data_dir = _rapiflex_folder(tmp_path)
        (data_dir / "slide_overview.jpg").touch()
        _write_mis(data_dir / "sample.mis", image_file="sample_0000.jpg")

        # Not caplog: setup_logging sets propagate = False on the `thyra`
        # logger process-wide, so once any test in the session has invoked the
        # CLI, caplog's root handler stops seeing Thyra records entirely.
        with thyra_logs("thyra", logging.WARNING) as records:
            info = BrukerFolderStructure(data_dir).analyze()

        assert info.primary_optical_image is None
        assert [p.name for p in info.optical_images] == ["slide_overview.jpg"]
        assert any("sample_0000.jpg" in record.getMessage() for record in records)

    def test_mis_image_file_with_a_windows_directory_part(self, tmp_path):
        """Only the filename is used, wherever Thyra runs."""
        data_dir = _rapiflex_folder(tmp_path)
        (data_dir / "sample_0000.jpg").touch()
        _write_mis(data_dir / "sample.mis", image_file=r"images\sample_0000.jpg")

        info = BrukerFolderStructure(data_dir).analyze()

        assert info.primary_optical_image == data_dir / "sample_0000.jpg"

    def test_find_teaching_points_file(self, tmp_path):
        """Test finding .mis teaching points file."""
        # Create Rapiflex structure with .mis file
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "sample.dat").touch()
        (data_dir / "sample_poslog.txt").touch()
        (data_dir / "sample_info.txt").touch()
        (data_dir / "sample.mis").touch()

        folder = BrukerFolderStructure(data_dir)
        info = folder.analyze()

        assert info.teaching_points_file is not None
        assert info.teaching_points_file.name == "sample.mis"

    def test_find_metadata_files_rapiflex(self, tmp_path):
        """Test finding Rapiflex metadata files."""
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "sample.dat").touch()
        (data_dir / "sample_poslog.txt").touch()
        (data_dir / "sample_info.txt").touch()
        (data_dir / "sample.mis").touch()

        folder = BrukerFolderStructure(data_dir)
        info = folder.analyze()

        assert "data" in info.metadata_files
        assert "poslog" in info.metadata_files
        assert "info" in info.metadata_files
        assert "mis" in info.metadata_files

    def test_find_metadata_files_timstof(self, tmp_path):
        """Test finding timsTOF metadata files."""
        d_dir = tmp_path / "test.d"
        d_dir.mkdir()
        (d_dir / "analysis.tsf").touch()
        (d_dir / "analysis.tsf_bin").touch()

        folder = BrukerFolderStructure(d_dir)
        info = folder.analyze()

        assert "tsf" in info.metadata_files
        assert "tsf_bin" in info.metadata_files

    def test_nonexistent_path_raises(self, tmp_path):
        """Test that nonexistent path raises ValueError."""
        nonexistent = tmp_path / "nonexistent"

        folder = BrukerFolderStructure(nonexistent)
        with pytest.raises(ConversionRefused, match="does not exist"):
            folder.analyze()

    def test_classmethod_detect_format(self, tmp_path):
        """Test classmethod detect_format."""
        d_dir = tmp_path / "test.d"
        d_dir.mkdir()
        (d_dir / "analysis.tdf").touch()

        fmt = BrukerFolderStructure.detect_format(d_dir)
        assert fmt == BrukerFormat.TIMSTOF

    def test_classmethod_is_bruker_data(self, tmp_path):
        """Test classmethod is_bruker_data."""
        # Valid Bruker data
        d_dir = tmp_path / "test.d"
        d_dir.mkdir()
        (d_dir / "analysis.tsf").touch()

        assert BrukerFolderStructure.is_bruker_data(d_dir) is True

        # Invalid data
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()

        assert BrukerFolderStructure.is_bruker_data(empty_dir) is False

    def test_analyze_caching(self, tmp_path):
        """Test that analyze() result is cached."""
        d_dir = tmp_path / "test.d"
        d_dir.mkdir()
        (d_dir / "analysis.tsf").touch()

        folder = BrukerFolderStructure(d_dir)
        info1 = folder.analyze()
        info2 = folder.analyze()

        # Should be the same object (cached)
        assert info1 is info2


class TestBrukerFolderInfo:
    """Tests for BrukerFolderInfo dataclass."""

    def test_dataclass_fields(self, tmp_path):
        """Test BrukerFolderInfo has expected fields."""
        info = BrukerFolderInfo(
            path=tmp_path,
            format=BrukerFormat.TIMSTOF,
            data_path=tmp_path / "data.d",
            optical_images=[tmp_path / "image.tif"],
            teaching_points_file=tmp_path / "points.mis",
            metadata_files={"tdf": tmp_path / "analysis.tdf"},
        )

        assert info.path == tmp_path
        assert info.format == BrukerFormat.TIMSTOF
        assert info.data_path == tmp_path / "data.d"
        assert len(info.optical_images) == 1
        assert info.teaching_points_file == tmp_path / "points.mis"
        assert "tdf" in info.metadata_files

    def test_default_values(self, tmp_path):
        """Test BrukerFolderInfo default values."""
        info = BrukerFolderInfo(
            path=tmp_path,
            format=BrukerFormat.UNKNOWN,
            data_path=tmp_path,
        )

        assert info.optical_images == []
        assert info.teaching_points_file is None
        assert info.metadata_files == {}
