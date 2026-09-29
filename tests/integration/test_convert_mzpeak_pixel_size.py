"""The pixel size an mzPeak conversion writes into the store.

An archive made from an imzML of before 2017 carries the AREA of a pixel
under IMS:1000046. Thyra 4.2.0 wrote that number into the store as a
length: 10000 um for a 100 um pixel, and ``convert_msi`` returned ``True``.
These tests read the value back from the store's ``zarr.json``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

import pytest

from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra.convert import convert_msi

pytest.importorskip("spatialdata", reason="SpatialData not installed")

OLD_PIXEL_SIZE = {
    "name": "pixel size",
    "accession": "IMS:1000046",
    "value": 10000,
    "unit": None,
}
EXTENT_X = {
    "name": "max dimension x",
    "accession": "IMS:1000044",
    "value": 500,
    "unit": "UO:0000017",
}
COUNT_X = {
    "name": "max count of pixels x",
    "accession": "IMS:1000042",
    "value": 5,
    "unit": None,
}


def _archive(tmp_path: Path, parameters: List[dict]) -> Path:
    """A 5 x 3 archive with these scan settings."""
    return build_mzpeak(
        tmp_path / "source.mzpeak",
        grid_spectra(5, 3, n_points=6),
        pixel_size=None,
        footer_metadata=False,
        index_metadata={
            "scan_settings_list": [
                {
                    "id": "scanSettings1",
                    "source_file_refs": [],
                    "targets": [],
                    "parameters": parameters,
                }
            ]
        },
    )


def _convert(source: Path, output: Path, pixel_size_um: Optional[float] = None):
    """Convert, leaving the pixel size to the file unless one is given."""
    return convert_msi(
        str(source),
        str(output),
        format_type="spatialdata",
        dataset_id="area",
        pixel_size_um=pixel_size_um,
    )


def _stored_pixel_size(store: Path):
    """``(x, y)`` pixel size of the global coordinate system, in um."""
    attributes = json.loads((store / "zarr.json").read_text(encoding="utf-8"))
    system = attributes["attributes"]["coordinate_systems"]["global"]
    return (system["pixel_size_um_x"], system["pixel_size_um_y"])


def test_an_area_is_stored_as_the_side_of_the_pixel(tmp_path):
    """Pixel size 10000, 5 pixels on x, extent 500 um: 100 um in the store."""
    archive = _archive(tmp_path, [OLD_PIXEL_SIZE, EXTENT_X, COUNT_X])
    store = tmp_path / "area.zarr"

    assert _convert(archive, store) is True

    assert _stored_pixel_size(store) == pytest.approx((100.0, 100.0))


def test_an_untested_value_asks_for_the_pixel_size(tmp_path):
    """No count and extent: the conversion is refused, as for imzML."""
    archive = _archive(tmp_path, [OLD_PIXEL_SIZE])
    store = tmp_path / "untested.zarr"

    assert _convert(archive, store) is False

    assert not store.exists()


def test_a_stated_pixel_size_is_used(tmp_path):
    """``--pixel-size`` is the way through, and it is what gets stored."""
    archive = _archive(tmp_path, [OLD_PIXEL_SIZE])
    store = tmp_path / "stated.zarr"

    assert _convert(archive, store, pixel_size_um=100.0) is True

    assert _stored_pixel_size(store) == pytest.approx((100.0, 100.0))
