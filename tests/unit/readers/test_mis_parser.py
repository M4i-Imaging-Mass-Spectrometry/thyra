"""Tests for the .mis file parser and discovery helpers."""

from pathlib import Path

from thyra.readers.bruker.mis_parser import find_mis_file_for_d_folder, parse_mis_file


def _write_mis(tmp_path: Path, name: str, raster: str = "5,5") -> Path:
    mis = tmp_path / name
    mis.write_text(
        f"""<?xml version="1.0"?>
<ImagingSequence>
<ImageFile>img.tif</ImageFile>
<Raster>{raster}</Raster>
<Area Name="01"><Point>10,20</Point><Point>30,40</Point></Area>
</ImagingSequence>
"""
    )
    return mis


def test_parse_mis_extracts_raster(tmp_path: Path) -> None:
    mis = _write_mis(tmp_path, "sample.mis", raster="5,5")
    data = parse_mis_file(mis)
    assert data["raster"] == [5, 5]


def test_parse_mis_handles_rectangular_raster(tmp_path: Path) -> None:
    mis = _write_mis(tmp_path, "sample.mis", raster="10,20")
    data = parse_mis_file(mis)
    assert data["raster"] == [10, 20]


def test_find_mis_prefers_matching_stem(tmp_path: Path) -> None:
    """When several .mis files sit in the parent, prefer the one whose stem
    matches the .d folder stem.
    """
    d_folder = tmp_path / "sample_A.d"
    d_folder.mkdir()
    _write_mis(tmp_path, "sample_A.mis", raster="5,5")
    _write_mis(tmp_path, "sample_B.mis", raster="50,50")

    found = find_mis_file_for_d_folder(d_folder)
    assert found is not None
    assert found.name == "sample_A.mis"


def test_find_mis_falls_back_to_any(tmp_path: Path) -> None:
    d_folder = tmp_path / "sample.d"
    d_folder.mkdir()
    other = _write_mis(tmp_path, "different_name.mis", raster="7,7")

    found = find_mis_file_for_d_folder(d_folder)
    assert found == other


def test_find_mis_returns_none_when_missing(tmp_path: Path) -> None:
    d_folder = tmp_path / "sample.d"
    d_folder.mkdir()
    assert find_mis_file_for_d_folder(d_folder) is None
