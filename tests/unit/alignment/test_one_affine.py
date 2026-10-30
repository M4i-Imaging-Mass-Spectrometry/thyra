"""One affine places every region, with or without teaching points (#520, D75).

Without the D31 lattice the raster is stretched over the box of the mapped
Areas. Before #520 the TIC image used that box while each region's polygons
used its own Area, a lone region ignored where the raster started, and a
region with no Area got no polygon at all. These tests pin the one rule:
every measured pixel, of a region with an Area or without one, is the TIC
cell :attr:`AreaAlignmentResult.cell_to_image` draws.
"""

import logging
from typing import Any, Dict, List

import numpy as np
import pytest

from thyra.alignment.affine import AffineTransform, points_span_plane
from thyra.alignment.teaching_points import TeachingPointAlignment, fit_raster_lattice


def _positions(region: int, xs: range, ys: range) -> List[Dict[str, Any]]:
    return [{"region": region, "raster_x": x, "raster_y": y} for y in ys for x in xs]


def _cell(matrix: np.ndarray, nx: int, ny: int) -> List[tuple]:
    corners = np.array(
        [[nx, ny, 1.0], [nx + 1, ny, 1.0], [nx + 1, ny + 1, 1.0], [nx, ny + 1, 1.0]]
    )
    return [(float(x), float(y)) for x, y, _ in corners @ matrix.T]


def test_a_lone_region_is_placed_from_where_the_raster_starts() -> None:
    """F170: the raster starts three steps before the region's first spot.

    The Area spans x 40..100 for three spots, 20 px a step. TIC cell 3 holds
    the first spot, so it must land on the Area's first 20 px.
    """
    area = {"name": "A0", "p1": [40, 20], "p2": [100, 60]}
    positions = _positions(0, range(100, 103), range(10, 12))
    result = TeachingPointAlignment().compute_area_alignment(
        [area], positions, first_raster_x=97, first_raster_y=10
    )

    matrix = result.cell_to_image
    assert matrix is not None
    assert result.transform_point(3, 0) == pytest.approx((50.0, 30.0))
    np.testing.assert_allclose(
        result.cell_corners(3, 0), _cell(matrix, 3, 0), rtol=1e-6
    )
    assert result.cell_corners(0, 0) is None  # no spot at raster x 97


def test_every_region_is_placed_by_the_one_matrix() -> None:
    """F173: Area 0 drawn a step too wide no longer splits TIC from polygon."""
    areas = [
        {"name": "A0", "p1": [20, 20], "p2": [100, 40]},
        {"name": "A1", "p1": [20, 40], "p2": [80, 60]},
    ]
    positions = _positions(0, range(100, 103), range(10, 11)) + _positions(
        1, range(100, 103), range(11, 12)
    )
    result = TeachingPointAlignment().compute_area_alignment(areas, positions, 100, 10)

    matrix = result.cell_to_image
    assert matrix is not None
    for p in positions:
        nx, ny = p["raster_x"] - 100, p["raster_y"] - 10
        np.testing.assert_allclose(
            result.cell_corners(nx, ny), _cell(matrix, nx, ny), rtol=1e-6
        )
        centre = (matrix @ np.array([nx + 0.5, ny + 0.5, 1.0]))[:2]
        assert result.transform_point(nx, ny) == pytest.approx(tuple(centre))


@pytest.mark.parametrize("missing", ["no Area", "placeholder"])
def test_a_region_without_an_area_still_gets_its_cells(missing: str) -> None:
    """F171, F118: region 1 has no Area, or one the parse could not read."""
    areas: List[Dict[str, Any]] = [{"name": "A0", "p1": [20, 20], "p2": [80, 40]}]
    if missing == "placeholder":
        areas.append({"name": "A1", "unusable": "one point"})
    positions = _positions(0, range(100, 103), range(10, 11)) + _positions(
        1, range(100, 103), range(11, 12)
    )
    result = TeachingPointAlignment().compute_area_alignment(areas, positions, 100, 10)

    assert [m.region_id for m in result.region_mappings] == [0]
    matrix = result.cell_to_image
    assert matrix is not None
    np.testing.assert_allclose(
        result.cell_corners(0, 1), _cell(matrix, 0, 1), rtol=1e-6
    )


def test_an_area_with_no_width_places_nothing() -> None:
    """F169: a singular matrix put zero-area polygons in the store."""
    area = {"name": "A0", "unusable": "no width or no height"}
    result = TeachingPointAlignment().compute_area_alignment(
        [area], _positions(0, range(3), range(2)), 0, 0
    )
    assert result.cell_to_image is None
    assert result.cell_corners(0, 0) is None


# --- Teaching points that do not span the photo (F198) ----------------------


SPANNING = [[100, 100], [900, 100], [100, 700]]
DEGENERATE = {
    "collinear, diagonal": [[100, 100], [300, 300], [500, 500]],
    "two identical": [[100, 100], [100, 100], [100, 700]],
    "collinear, along x": [[0, 0], [10, 0], [20, 0]],
}


def test_spanning_points_span_the_plane() -> None:
    assert points_span_plane(SPANNING)


@pytest.mark.parametrize("points", DEGENERATE.values(), ids=DEGENERATE.keys())
def test_degenerate_points_do_not(points) -> None:
    assert not points_span_plane(points)
    with pytest.raises(ValueError, match="one line"):
        AffineTransform.from_points(points, SPANNING)
    with pytest.raises(ValueError, match="one line"):
        AffineTransform.from_points(SPANNING, points)


@pytest.mark.parametrize("points", DEGENERATE.values(), ids=DEGENERATE.keys())
def test_degenerate_teaching_points_fit_no_lattice(points, thyra_logs) -> None:
    """The fit used to spread three spots over 1e17 photo pixels."""
    teaching = [{"image": p, "stage": [p[0] * 8, p[1] * 8]} for p in points]
    area = {"name": "A0", "p1": [20, 20], "p2": [80, 60]}
    with thyra_logs("thyra.alignment", logging.WARNING) as records:
        fit = fit_raster_lattice(
            teaching,
            None,
            (20.0, 20.0),
            [area],
            _positions(0, range(3), range(2)),
            0,
            0,
        )
    assert fit is None
    assert any("cannot place the raster" in r.getMessage() for r in records)
