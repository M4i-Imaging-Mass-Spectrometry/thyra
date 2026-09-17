# tests/unit/converters/test_optical_crop_to_region.py
"""The alignment image is cropped to the one region a conversion covers.

Every ``.d`` cut from a FlexImaging slide names the same whole-slide scan in
its ``.mis``, so a project of eleven sections held eleven copies of one
photo, each placed so that its own Area landed on its own MSI and the other
ten sections' pixels were centimetres off. Measured on a real store: five
copies, composed placements up to 54,873 um apart.

The contract has two ends and both are pinned here, because getting the
crop's transform wrong moves every registered overlay by the crop's origin,
silently:

* **coordinates**: the crop's element carries a translation by its origin,
  composed in front of whatever the whole scan carried, so a raster position
  lands on the same scan pixel through the crop as it did through the whole
  image -- checked in optical-pixel mode and in micrometre mode;
* **bytes**: the pixels streamed into the store are the scan's own window.

And the two cases that must not change: a multi-region file converted whole,
and an Area that already covers the scan.

The geometry is the collaborator test's -- one Area at (100, 200)-(500, 600)
over a 4 x 3 raster -- on a scan large enough to have something outside it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
import tifffile
from spatialdata.transformations import Identity, Translation

from thyra.converters.spatialdata.optical_image import (
    CROP_MARGIN_FRACTION,
    CroppedOpticalSource,
    OpticalCrop,
    OpticalImages,
    OpticalTiffSource,
    probe_optical_source,
)

from .test_optical_images_collaborator import (
    AREA,
    IMAGE_FILE,
    N_X,
    N_Y,
    PIXEL_SIZE,
    _optical,
    _StubReader,
)

#: The scan: 700 rows by 800 columns, so the Area (100..500, 200..600) has
#: slide around it on every side and a 10% margin still fits.
SCAN_ROWS, SCAN_COLS = 700, 800
#: The Area padded by CROP_MARGIN_FRACTION of its own 400 x 400 box.
EXPECTED_CROP = OpticalCrop(
    x0=60,
    y0=160,
    width=480,
    height=480,
    full_width=SCAN_COLS,
    full_height=SCAN_ROWS,
    region_id=0,
    region_name="Area0",
)
AXES = ("x", "y")


@pytest.fixture
def scan(tmp_path: Path) -> Path:
    """A whole-slide scan, striped so a crop decodes only some strips."""
    rng = np.random.default_rng(23)
    path = tmp_path / IMAGE_FILE
    tifffile.imwrite(
        str(path),
        rng.integers(0, 256, size=(SCAN_ROWS, SCAN_COLS, 3), dtype=np.uint8),
        rowsperstrip=64,
        compression="lzw",
    )
    return path


def _scan_pixels(path: Path) -> np.ndarray:
    return np.moveaxis(tifffile.imread(str(path)), -1, 0)


def _aligned(tmp_path: Path, reader: _StubReader, **kwargs: Any) -> OpticalImages:
    """A collaborator that has done everything up to declaring its images."""
    optical = _optical(tmp_path, reader, **kwargs)
    optical.compute_alignment()
    optical.build_tic_to_image_affine()
    optical.add_images({})
    return optical


def _matrix(transform: Any) -> np.ndarray:
    return transform.to_affine_matrix(input_axes=AXES, output_axes=AXES)


def _apply(matrix: np.ndarray, x: float, y: float) -> np.ndarray:
    return (matrix @ np.array([x, y, 1.0]))[:2]


# --- The crop ---------------------------------------------------------------


def test_a_single_region_conversion_crops_the_alignment_image(
    tmp_path: Path, scan: Path
) -> None:
    optical = _aligned(tmp_path, _StubReader([scan]))

    assert optical.crop == EXPECTED_CROP
    pending = optical.pending["ds_optical_highres"]
    assert isinstance(pending.source, CroppedOpticalSource)
    assert pending.source.shape == (3, 480, 480)
    # The pyramid is the crop's, not the scan's: 480 px halves below the
    # 1000 px floor at once, so there is none.
    assert pending.scale_factors == []
    # The other images still scale into the WHOLE scan's grid.
    assert optical._primary_dims == (SCAN_COLS, SCAN_ROWS)


def test_the_margin_is_a_tenth_of_the_box_on_every_side() -> None:
    """Pinned as arithmetic, so the constant cannot drift from the geometry."""
    crop = OpticalCrop.around(
        (100.0, 500.0, 200.0, 600.0), SCAN_COLS, SCAN_ROWS, mapping=_mapping()
    )
    assert crop is not None
    pad = 400 * CROP_MARGIN_FRACTION
    assert (crop.x0, crop.y0) == (100 - pad, 200 - pad)
    assert (crop.width, crop.height) == (400 + 2 * pad, 400 + 2 * pad)


def test_the_crop_is_clamped_to_the_scan() -> None:
    crop = OpticalCrop.around(
        (-50.0, 120.0, 650.0, 900.0), SCAN_COLS, SCAN_ROWS, mapping=_mapping()
    )
    assert crop is not None
    assert (crop.x0, crop.y0) == (0, 625)
    assert crop.x0 + crop.width == 137
    assert crop.y0 + crop.height == SCAN_ROWS


# --- Both ends of the transform contract ------------------------------------


def test_the_crop_carries_the_translation_that_puts_it_back(
    tmp_path: Path, scan: Path
) -> None:
    """Optical-pixel mode: the scan's grid stays ``global``; the crop moves.

    A raster position lands on scan pixel ``T(r)`` through the affine. In
    the crop that pixel is ``T(r) - origin``, and the element's transform
    must take it back to ``T(r)`` -- otherwise every polygon and every
    overlay is off by the origin.
    """
    optical = _aligned(tmp_path, _StubReader([scan]))
    transform = optical.pending["ds_optical_highres"].transformations["global"]

    assert isinstance(transform, Translation)
    origin = np.array([EXPECTED_CROP.x0, EXPECTED_CROP.y0], dtype=float)
    element = _matrix(transform)
    np.testing.assert_allclose(
        element, [[1.0, 0.0, 60.0], [0.0, 1.0, 160.0], [0.0, 0.0, 1.0]]
    )

    tic = optical.tic_to_image
    assert tic is not None
    for rx, ry in [(0, 0), (N_X - 1, N_Y - 1), (2, 1)]:
        on_scan = _apply(tic, rx, ry)
        in_crop = on_scan - origin
        # End one: the crop pixel the raster position sits on.
        assert (in_crop >= 0).all() and in_crop[0] < 480 and in_crop[1] < 480
        # End two: the element puts that pixel back on the scan.
        np.testing.assert_allclose(_apply(element, *in_crop), on_scan)


def test_in_micrometre_mode_the_translation_comes_before_the_affine(
    tmp_path: Path, scan: Path
) -> None:
    """The opt-out route Ousia takes: photo into micrometres, crop first.

    Through the whole scan a raster position ``r`` came out at
    ``pixel_size * r`` micrometres. Through the crop it must still: the
    crop pixel ``T(r) - origin`` is translated back onto the scan and then
    carried by the same inverse affine as before.
    """
    optical = _aligned(tmp_path, _StubReader([scan]), apply_alignment=False)
    transform = optical.pending["ds_optical_highres"].transformations["global"]
    composed = _matrix(transform)

    tic = optical.tic_to_image
    assert tic is not None
    ps_x, ps_y = PIXEL_SIZE
    to_um = np.diag([ps_x, ps_y, 1.0]) @ np.linalg.inv(tic)
    shift = np.array([[1.0, 0.0, 60.0], [0.0, 1.0, 160.0], [0.0, 0.0, 1.0]])
    np.testing.assert_allclose(composed, to_um @ shift)

    origin = np.array([60.0, 160.0])
    for rx, ry in [(0, 0), (N_X - 1, N_Y - 1), (1, 2)]:
        in_crop = _apply(tic, rx, ry) - origin
        np.testing.assert_allclose(
            _apply(composed, *in_crop), [rx * ps_x, ry * ps_y], atol=1e-9
        )
    # Both coordinate systems get the same transform object, as ever.
    assert optical.pending["ds_optical_highres"].transformations["ds"] is transform


def test_the_streamed_pixels_are_the_scans_own_window(
    tmp_path: Path, scan: Path
) -> None:
    """The bytes end: what lands in the store is the window, nothing resampled."""
    pytest.importorskip("spatialdata")
    from spatialdata import SpatialData
    from spatialdata.transformations import get_transformation

    optical = _aligned(tmp_path, _StubReader([scan]))
    name = "ds_optical_highres"
    streamed = optical.pending[name]
    store = tmp_path / "out.zarr"
    SpatialData(images={name: streamed.placeholder()}).write(store)
    assert optical.stream_pending_pixels() == 1

    element = SpatialData.read(str(store)).images[name]
    crop = EXPECTED_CROP
    window = _scan_pixels(scan)[
        :, crop.y0 : crop.y0 + crop.height, crop.x0 : crop.x0 + crop.width
    ]
    np.testing.assert_array_equal(np.asarray(element.values), window)
    np.testing.assert_allclose(
        _matrix(get_transformation(element, "global")),
        [[1.0, 0.0, 60.0], [0.0, 1.0, 160.0], [0.0, 0.0, 1.0]],
    )


def test_a_cropped_source_decodes_only_the_window(scan: Path) -> None:
    """The source contract, band by band: rows re-based, columns cut."""
    source = CroppedOpticalSource(OpticalTiffSource.probe(scan), EXPECTED_CROP)
    whole = _scan_pixels(scan)

    assert source.shape == (3, 480, 480)
    assert source.row_bytes == 3 * 480
    rows_seen = []
    for first_row, band in source.bands(100):
        rows_seen.append(first_row)
        assert band.shape[0] == 3 and band.shape[2] == 480
        np.testing.assert_array_equal(
            band, whole[:, 160 + first_row : 160 + first_row + band.shape[1], 60:540]
        )
    assert rows_seen == [0, 100, 200, 300, 400]


def test_a_raster_format_is_cropped_after_the_whole_decode(tmp_path: Path) -> None:
    """A PNG decodes whole or not at all; the window is what comes out."""
    PILImage = pytest.importorskip("PIL.Image")
    rng = np.random.default_rng(5)
    pixels = rng.integers(0, 256, size=(SCAN_ROWS, SCAN_COLS, 3), dtype=np.uint8)
    path = tmp_path / "scan.png"
    PILImage.fromarray(pixels).save(path)

    source = CroppedOpticalSource(probe_optical_source(path), EXPECTED_CROP)
    bands = list(source.bands(50))

    assert len(bands) == 1
    first_row, band = bands[0]
    assert first_row == 0
    np.testing.assert_array_equal(band, np.moveaxis(pixels, -1, 0)[:, 160:640, 60:540])


# --- What the store says ------------------------------------------------------


def test_the_root_attr_records_the_crop(tmp_path: Path, scan: Path) -> None:
    optical = _aligned(tmp_path, _StubReader([scan]))

    assert optical.root_attr() == {
        "alignment_element": "ds_optical_highres",
        "elements": {
            "ds_optical_highres": {
                "source_file": IMAGE_FILE,
                "crop": {
                    "origin": [60, 160],
                    "size": [480, 480],
                    "full_size": [SCAN_COLS, SCAN_ROWS],
                    "region_id": 0,
                    "region_name": "Area0",
                    "margin_fraction": CROP_MARGIN_FRACTION,
                },
            }
        },
    }


# --- What must not change -----------------------------------------------------


def test_a_multi_region_file_converted_whole_keeps_the_whole_scan(
    tmp_path: Path, scan: Path
) -> None:
    """Fused conversions are not this change's business: same image, same transform."""
    reader = _StubReader(
        [scan],
        areas=[AREA, {"name": "Area1", "p1": (550, 200), "p2": (750, 600)}],
    )
    reader._positions = [
        {"region": r, "raster_x": x + 10 * r, "raster_y": y}
        for r in (0, 1)
        for y in range(N_Y)
        for x in range(N_X)
    ]
    optical = _aligned(tmp_path, reader)

    assert len(optical.alignment.region_mappings) == 2
    assert optical.crop is None
    pending = optical.pending["ds_optical_highres"]
    assert isinstance(pending.source, OpticalTiffSource)
    assert pending.source.shape == (3, SCAN_ROWS, SCAN_COLS)
    assert isinstance(pending.transformations["global"], Identity)
    assert "crop" not in optical.root_attr()["elements"]["ds_optical_highres"]


def test_a_selected_region_of_a_multi_region_file_is_cropped_to_it(
    tmp_path: Path, scan: Path
) -> None:
    """``region=`` names one region; the crop is where THAT region's raster lands.

    The reader still maps every Area, so the affine is the global stretch,
    and the box is taken from it rather than from the Area's corners: it is
    the affine that places the TIC, and the crop has to agree with the TIC.
    """
    reader = _StubReader(
        [scan],
        areas=[AREA, {"name": "Area1", "p1": (550, 200), "p2": (750, 600)}],
    )
    reader._positions = [
        {"region": r, "raster_x": x + 10 * r, "raster_y": y}
        for r in (0, 1)
        for y in range(N_Y)
        for x in range(N_X)
    ]
    # What the timsTOF reader does when region 1 is selected: the raster
    # origin becomes that region's own minimum.
    reader._selected_region = 1
    reader._header = {"first_raster_x": 10, "first_raster_y": 0}
    optical = _aligned(tmp_path, reader)

    crop = optical.crop
    assert crop is not None
    assert (crop.region_id, crop.region_name) == (1, "Area1")
    tic = optical.tic_to_image
    assert tic is not None
    # Region 1's raster is TIC indices 0..N_X-1 by 0..N_Y-1 (its own
    # origin), so its footprint under the affine is the box the crop pads.
    sx, sy = tic[0, 0], tic[1, 1]
    x_lo, x_hi = tic[0, 2] - sx / 2, tic[0, 2] + (N_X - 1) * sx + sx / 2
    y_lo, y_hi = tic[1, 2] - sy / 2, tic[1, 2] + (N_Y - 1) * sy + sy / 2
    expected = OpticalCrop.around(
        (x_lo, x_hi, y_lo, y_hi),
        SCAN_COLS,
        SCAN_ROWS,
        mapping=optical.alignment.region_mappings[1],
    )
    assert crop == expected
    transform = optical.pending["ds_optical_highres"].transformations["global"]
    np.testing.assert_allclose(
        _matrix(transform), [[1.0, 0.0, crop.x0], [0.0, 1.0, crop.y0], [0.0, 0.0, 1.0]]
    )


def test_an_area_covering_the_scan_leaves_it_uncropped(
    tmp_path: Path, scan: Path
) -> None:
    reader = _StubReader(
        [scan], areas=[{"name": "Whole", "p1": (0, 0), "p2": (SCAN_COLS, SCAN_ROWS)}]
    )
    optical = _aligned(tmp_path, reader)

    assert optical.crop is None
    assert isinstance(
        optical.pending["ds_optical_highres"].transformations["global"], Identity
    )


def test_no_alignment_means_no_crop(tmp_path: Path, scan: Path) -> None:
    """Without an affine there is nothing to crop against."""
    reader = _StubReader([scan], areas=[])
    optical = _aligned(tmp_path, reader)

    assert optical.tic_to_image is None
    assert optical.crop is None
    assert optical.pending["ds_optical_highres"].source.shape == (
        3,
        SCAN_ROWS,
        SCAN_COLS,
    )


# --- helpers ------------------------------------------------------------------


def _mapping(region_id: int = 0, name: str = "Area0") -> Any:
    from thyra.alignment.teaching_points import RegionMapping

    return RegionMapping(
        region_id=region_id,
        name=name,
        raster_min_x=0,
        raster_max_x=N_X - 1,
        raster_min_y=0,
        raster_max_y=N_Y - 1,
        image_min_x=100,
        image_max_x=500,
        image_min_y=200,
        image_max_y=600,
    )
