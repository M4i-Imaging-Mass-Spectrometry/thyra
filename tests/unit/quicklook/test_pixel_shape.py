"""Pixels that are not square are drawn as wide and as tall as they are.

The spacing along a scan line and the spacing between lines are set by
different things, so a store can state two pixel sizes. Both commands read
them from the root attributes ``pixel_size_x_um`` and ``pixel_size_y_um``.
These tests go through the commands, on stores built here.
"""

from __future__ import annotations

import base64
import io
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest
from click.testing import CliRunner
from PIL import Image

from tests.fixtures.hand_built_store import MZ, make_table, write_store
from thyra.quicklook.cli import ion_image_command, report_command
from thyra.quicklook.render import VIRIDIS

#: The raster of the stores below: columns, rows.
WIDTH, HEIGHT = 30, 12

#: Where the figure puts its raster: the margin on its left.
_FIGURE_LEFT = 16

#: A pixel is drawn to scale up to this ratio of its two sides.
_BOUND = 20


def _positions(width: int = WIDTH, height: int = HEIGHT) -> List[Tuple[int, int]]:
    """Every position of a full rectangle that does not start at zero."""
    return [(4 + x, 2 + y) for y in range(height) for x in range(width)]


def _store(
    folder: Path,
    attrs: Optional[Dict[str, Any]] = None,
    width: int = WIDTH,
    height: int = HEIGHT,
) -> Path:
    """A store whose every pixel holds the same spectrum, with a TIC image.

    One value everywhere makes the raster a single colour, so its edges
    can be found in the figure. The mean spectrum has one peak, so the
    report shows one ion image.
    """
    positions = _positions(width, height)
    matrix = np.tile(np.arange(float(MZ.size), 0.0, -1.0), (len(positions), 1))
    table = make_table(positions=positions, matrix=matrix)
    tic = np.zeros((1, 2 + height, 4 + width))
    tic[0, 2:, 4:] = matrix.sum(axis=1).reshape(height, width)
    folder.mkdir(parents=True, exist_ok=True)
    return write_store(
        folder / "scan.zarr", {"ds_z0": table}, images={"ds_z0_tic": [tic]}, attrs=attrs
    )


def _sizes(x: Any, y: Any) -> Dict[str, Any]:
    return {"pixel_size_x_um": x, "pixel_size_y_um": y}


def _ion_image_bytes(store: Path) -> bytes:
    target = store.parent / "ion.png"
    result = CliRunner().invoke(
        ion_image_command, [str(store), "505", "-o", str(target)]
    )
    assert result.exit_code == 0, result.output
    return target.read_bytes()


def _report(store: Path) -> str:
    target = store.parent / "report.html"
    result = CliRunner().invoke(report_command, [str(store), "-o", str(target)])
    assert result.exit_code == 0, result.output
    return target.read_text(encoding="utf-8")


def _raster_size(png: bytes) -> Tuple[int, int]:
    """Width and height of the one-colour raster inside an ion-image figure."""
    figure = np.asarray(Image.open(io.BytesIO(png)))
    filled = (figure == VIRIDIS[255]).all(axis=2)
    rows = np.flatnonzero(filled[:, _FIGURE_LEFT])
    assert rows.size and rows[-1] - rows[0] + 1 == rows.size
    line = filled[rows[0], _FIGURE_LEFT:]
    width = int(line.size if line.all() else np.argmin(line))
    assert filled[rows[0] : rows[-1] + 1, _FIGURE_LEFT : _FIGURE_LEFT + width].all()
    return width, int(rows.size)


def _embedded_sizes(page: str) -> List[Tuple[int, int]]:
    """The size of each image embedded in a report, in page order."""
    found = re.findall(r'src="data:image/png;base64,([^"]+)"', page)
    return [Image.open(io.BytesIO(base64.b64decode(data))).size for data in found]


def _in_proportion(size: Tuple[int, int], ratio: float) -> bool:
    """Whether one pixel more or less on the shorter side gives ``ratio``."""
    width, height = size
    if ratio >= 1:
        return abs(width / ratio - height) <= 1
    return abs(height * ratio - width) <= 1


class TestTheDrawnProportions:
    """A raster of 30 by 12 pixels of 70 by 150 um covers 2100 by 1800 um."""

    AREA = (WIDTH * 70.0) / (HEIGHT * 150.0)

    def test_the_ion_image_has_the_proportions_of_the_measured_area(self, tmp_path):
        size = _raster_size(_ion_image_bytes(_store(tmp_path, _sizes(70.0, 150.0))))

        assert _in_proportion(size, self.AREA), size
        assert not _in_proportion(size, WIDTH / HEIGHT)

    def test_the_images_of_the_report_have_them_too(self, tmp_path):
        """The TIC image and the ion image, in that order."""
        sizes = _embedded_sizes(_report(_store(tmp_path, _sizes(70.0, 150.0))))

        assert len(sizes) == 2
        for size in sizes:
            assert _in_proportion(size, self.AREA), size

    def test_pixels_wider_than_tall_are_drawn_wider(self, tmp_path):
        store = _store(tmp_path, _sizes(150.0, 70.0))
        area = (WIDTH * 150.0) / (HEIGHT * 70.0)

        assert _in_proportion(_raster_size(_ion_image_bytes(store)), area)
        for size in _embedded_sizes(_report(store)):
            assert _in_proportion(size, area), size

    def test_a_small_raster_is_still_enlarged_to_be_seen(self, tmp_path):
        width, height = _raster_size(
            _ion_image_bytes(_store(tmp_path, _sizes(70.0, 150.0)))
        )

        assert max(width, height) >= 480
        assert width % WIDTH == 0  # along its shorter side every block is alike


class TestThePixelSizeIsStated:
    def test_the_report_states_both_sides(self, tmp_path):
        page = _report(_store(tmp_path, _sizes(70.0, 150.0)))

        assert "Each pixel is 70 x 150 µm, drawn to scale." in page

    def test_the_ion_image_gains_a_line_of_text(self, tmp_path):
        """The text is drawn, so it is measured: the figure grows by a line."""
        square = Image.open(
            io.BytesIO(_ion_image_bytes(_store(tmp_path / "a", _sizes(70.0, 70.0))))
        )
        stated = Image.open(
            io.BytesIO(_ion_image_bytes(_store(tmp_path / "b", _sizes(70.0, 70.5))))
        )

        # 70 by 70.5 um makes this raster one pixel taller, so the rest of
        # what the figure gains in height is the line.
        assert 10 < stated.height - square.height < 30

    def test_square_pixels_are_not_remarked_on(self, tmp_path):
        assert "Each pixel is" not in _report(_store(tmp_path, _sizes(70.0, 70.0)))


#: Root attributes that state no usable pixel size. An older store has
#: none; the others are what a damaged or hand-edited store could hold.
_NO_PIXEL_SIZE = {
    "both missing": {},
    "one missing": {"pixel_size_x_um": 70.0},
    "zero": _sizes(70.0, 0.0),
    "negative": _sizes(-70.0, 150.0),
    "not a number": _sizes(float("nan"), 150.0),
    "infinite": _sizes(70.0, float("inf")),
    "text": _sizes("70", "150"),
    "a truth value": _sizes(True, 150.0),
    "a list": _sizes([70.0], 150.0),
    "a whole number past every float": _sizes(10**400, 150.0),
}


@pytest.fixture(scope="module")
def square(tmp_path_factory) -> Tuple[bytes, str]:
    """What a store with square pixels of 20 um gives."""
    store = _store(tmp_path_factory.mktemp("square"), _sizes(20.0, 20.0))
    return _ion_image_bytes(store), _report(store)


class TestWithoutAUsablePixelSize:
    @pytest.mark.parametrize("attrs", _NO_PIXEL_SIZE.values(), ids=_NO_PIXEL_SIZE)
    def test_the_pixels_are_drawn_square_and_nothing_fails(
        self, tmp_path, square, attrs
    ):
        store = _store(tmp_path, attrs)

        assert _ion_image_bytes(store) == square[0]
        assert _report(store) == square[1]

    def test_sizes_that_differ_by_float_noise_are_square(self, tmp_path, square):
        store = _store(tmp_path, _sizes(20.0, 20.0 * (1 + 1e-12)))

        assert _ion_image_bytes(store) == square[0]
        assert _report(store) == square[1]


class TestAnExtremeRatio:
    """Past 1 to 20 a pixel is drawn at 1 to 20, and the text says so."""

    def test_a_figure_stops_at_the_bound(self, tmp_path):
        store = _store(tmp_path, _sizes(1.0, 5000.0), width=600, height=2)

        assert _raster_size(_ion_image_bytes(store)) == (600, 2 * _BOUND)

    def test_it_stops_there_the_other_way_round(self, tmp_path):
        store = _store(tmp_path, _sizes(5000.0, 1.0), width=2, height=600)

        assert _raster_size(_ion_image_bytes(store)) == (2 * _BOUND, 600)

    def test_the_report_says_what_it_drew(self, tmp_path):
        page = _report(_store(tmp_path, _sizes(1.0, 5000.0), width=600, height=2))

        assert "Each pixel is 1 x 5000 µm, drawn 1 to 20: not to scale." in page

    @pytest.mark.parametrize("sizes", [(1.0, 1e9), (1e9, 1.0), (1e-300, 1e300)])
    @pytest.mark.parametrize("shape", [(600, 1), (1, 600), (3, 2)])
    def test_no_image_is_huge_and_none_has_an_empty_side(self, tmp_path, sizes, shape):
        width, height = shape
        store = _store(tmp_path, _sizes(*sizes), width=width, height=height)

        figure = Image.open(io.BytesIO(_ion_image_bytes(store)))
        drawn = _embedded_sizes(_report(store))

        # A raster large enough is drawn pixel for pixel on its shorter
        # side; a small one is enlarged to less than twice 480 pixels.
        assert max(figure.size) < max(_BOUND * max(shape), 2 * 480) + 200
        assert len(drawn) == 2
        for size in drawn:
            assert min(size) >= 1 and max(size) < 2 * 480
