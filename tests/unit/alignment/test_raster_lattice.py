"""FlexImaging's raster lattice, placed on the alignment image (D31, #428).

FlexImaging puts every spot of a run on one lattice: one raster step apart in
the teaching-point frame, through the ``.mis`` reference point. An Area
outline only selects which nodes are measured. Thyra used to stretch each
region over its outline's bounding box instead, which sized pixels from 0.87
to 1.14 raster steps and moved spots by up to a step.

The synthetic tests build a stage that is rotated and mirrored against the
photo, measure every node inside each Area the way FlexImaging does, and
check the fit lands each spot exactly. The last tests use two public runs
whose flexImaging spot list says where each spot really is.
"""

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest
from shapely.geometry import Point, Polygon

from thyra.alignment.teaching_points import TeachingPointAlignment, fit_raster_lattice
from thyra.readers.bruker.mis_parser import parse_mis_file

FIXTURES = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "fixtures"
    / "fleximaging_msv000088438"
)

STEP = 100.0  # um
REFERENCE_NODE = (10, 6)


class _Stage:
    """A photo-to-stage affine: 8 um per pixel, rotated 0.7 degrees, y up."""

    def __init__(self, rotation_deg: float = 0.7, um_per_px: float = 8.0) -> None:
        th = np.radians(rotation_deg)
        rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        self.linear = um_per_px * rot @ np.diag([1.0, -1.0])
        self.offset = np.array([-20000.0, 15000.0])

    def to_stage(self, px: np.ndarray) -> np.ndarray:
        return px @ self.linear.T + self.offset

    def to_image(self, um: np.ndarray) -> np.ndarray:
        return (um - self.offset) @ np.linalg.inv(self.linear).T

    def teaching_points(self) -> List[Dict[str, Any]]:
        image = np.array([[100.0, 100.0], [5000.0, 120.0], [300.0, 4000.0]])
        stage = self.to_stage(image)
        return [{"image": list(i), "stage": list(s)} for i, s in zip(image, stage)]

    def node_image(
        self, x: np.ndarray, y: np.ndarray, ref_px: np.ndarray
    ) -> np.ndarray:
        """Image pixel of raster node (x, y) on FlexImaging's lattice."""
        ref = self.to_stage(ref_px[np.newaxis])[0]
        um = np.c_[
            ref[0] + STEP * (x - REFERENCE_NODE[0]),
            ref[1] - STEP * (y - REFERENCE_NODE[1]),
        ]
        return self.to_image(um)


def _acquisition(
    stage: _Stage,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], np.ndarray]:
    """Two Areas drawn on the photo, and every lattice node inside each.

    Area 0 is a rectangle whose edges stop a third of a step past the outer
    spots; Area 1 a pentagon. Both are drawn in photo pixels, axis-aligned
    there, while the lattice is rotated: the outline is not the raster.
    """
    ref_px = np.array([100.0, 100.0])
    nodes = np.array([(x, y) for x in range(0, 60) for y in range(0, 60)], dtype=float)
    node_px = stage.node_image(nodes[:, 0], nodes[:, 1], ref_px)

    def corner(x: float, y: float) -> List[float]:
        return list(stage.node_image(np.array([x]), np.array([y]), ref_px)[0])

    rect_lo, rect_hi = corner(12 - 0.33, 8 - 0.33), corner(18 + 0.33, 12 + 0.33)
    rect = [[rect_lo[0], rect_lo[1]], [rect_hi[0], rect_hi[1]]]
    pentagon = [
        corner(*p)
        for p in [(30.4, 20.2), (40.6, 21.5), (42.1, 30.3), (35.2, 36.4), (29.5, 29.0)]
    ]
    areas = [
        {
            "name": "rect",
            "type": 0,
            "points": rect,
            "p1": rect[0],
            "p2": rect[1],
            "raster": [STEP, STEP],
        },
        {
            "name": "pentagon",
            "type": 3,
            "points": pentagon,
            "p1": [min(v[0] for v in pentagon), min(v[1] for v in pentagon)],
            "p2": [max(v[0] for v in pentagon), max(v[1] for v in pentagon)],
            "raster": [STEP, STEP],
        },
    ]
    outlines = [
        Polygon(
            [
                (rect[0][0], rect[0][1]),
                (rect[1][0], rect[0][1]),
                (rect[1][0], rect[1][1]),
                (rect[0][0], rect[1][1]),
            ]
        ),
        Polygon(pentagon),
    ]
    positions = []
    for region, outline in enumerate(outlines):
        for (x, y), (px, py) in zip(nodes, node_px):
            if outline.contains(Point(px, py)):
                positions.append(
                    {"region": region, "raster_x": int(x), "raster_y": int(y)}
                )
    return areas, positions, ref_px


def _fit(stage: _Stage, **overrides: Any):
    areas, positions, ref_px = _acquisition(stage)
    first_x = min(p["raster_x"] for p in positions)
    first_y = min(p["raster_y"] for p in positions)
    kwargs: Dict[str, Any] = dict(
        teaching_points=stage.teaching_points(),
        reference_point=tuple(ref_px),
        raster_step=(STEP, STEP),
        areas=areas,
        positions=positions,
        first_raster_x=first_x,
        first_raster_y=first_y,
    )
    kwargs.update(overrides)
    return fit_raster_lattice(**kwargs), positions, ref_px, (first_x, first_y)


def test_every_spot_lands_on_its_lattice_node() -> None:
    stage = _Stage()
    fit, positions, ref_px, (fx, fy) = _fit(stage)

    assert fit is not None
    assert fit.reference_node == REFERENCE_NODE
    assert fit.axis_signs == (1, -1)
    assert (fit.spots_outside, fit.nodes_unmeasured) == (0, 0)
    xs = np.array([p["raster_x"] for p in positions], dtype=float)
    ys = np.array([p["raster_y"] for p in positions], dtype=float)
    truth = stage.node_image(xs, ys, ref_px)
    got = np.array([fit.cell_centre(x - fx, y - fy) for x, y in zip(xs, ys)])
    np.testing.assert_allclose(got, truth, atol=1e-6)


def test_a_cell_is_one_raster_step_whatever_the_outline() -> None:
    """The bug: cell size came from the outline. Here it is the step."""
    stage = _Stage()
    fit, _, _, _ = _fit(stage)

    assert fit is not None
    corners = np.array(fit.cell_corners(0, 0))
    step_px = STEP / 8.0
    assert np.linalg.norm(corners[1] - corners[0]) == pytest.approx(step_px)
    assert np.linalg.norm(corners[3] - corners[0]) == pytest.approx(step_px)
    # The centre is the middle of the cell, the SpatialData convention.
    np.testing.assert_allclose(corners.mean(axis=0), fit.cell_centre(0, 0))


def test_without_a_reference_point_the_first_teaching_point_is_used() -> None:
    """FlexImaging writes the reference point on the first teaching point."""
    stage = _Stage()
    fit, _, _, _ = _fit(stage, reference_point=None)

    assert fit is not None
    assert fit.reference_from == "first teaching point"
    assert fit.reference_node == REFERENCE_NODE


def test_a_mirrored_stage_is_found_from_the_areas() -> None:
    """Raster X running left on the stage still puts every spot inside."""
    stage = _Stage()
    areas, positions, ref_px = _acquisition(stage)
    mirrored = [dict(p, raster_x=100 - p["raster_x"]) for p in positions]
    fit = fit_raster_lattice(
        teaching_points=stage.teaching_points(),
        reference_point=tuple(ref_px),
        raster_step=(STEP, STEP),
        areas=areas,
        positions=mirrored,
        first_raster_x=min(p["raster_x"] for p in mirrored),
        first_raster_y=min(p["raster_y"] for p in mirrored),
    )

    assert fit is not None
    assert fit.axis_signs == (-1, -1)
    assert fit.spots_outside == 0


def test_the_alignment_result_places_measured_pixels_only() -> None:
    stage = _Stage()
    areas, positions, ref_px = _acquisition(stage)
    fx = min(p["raster_x"] for p in positions)
    fy = min(p["raster_y"] for p in positions)
    result = TeachingPointAlignment().compute_area_alignment(
        areas,
        positions,
        fx,
        fy,
        teaching_points=stage.teaching_points(),
        reference_point=tuple(ref_px),
        raster_step=(STEP, STEP),
    )

    assert result.lattice is not None
    np.testing.assert_allclose(result.cell_to_image, result.lattice.cell_to_image)
    p = positions[0]
    nx, ny = p["raster_x"] - fx, p["raster_y"] - fy
    assert result.transform_point(nx, ny) == pytest.approx(
        result.lattice.cell_centre(nx, ny)
    )
    assert result.cell_corners(nx, ny) == pytest.approx(
        result.lattice.cell_corners(nx, ny)
    )
    # A node no spot was measured on has no pixel.
    assert result.transform_point(59 - fx, 59 - fy) is None
    assert result.cell_corners(59 - fx, 59 - fy) is None


def test_fewer_than_three_teaching_points_fall_back(thyra_logs) -> None:
    stage = _Stage()
    areas, positions, ref_px = _acquisition(stage)
    with thyra_logs("thyra.alignment", logging.WARNING) as records:
        result = TeachingPointAlignment().compute_area_alignment(
            areas, positions, 0, 0, teaching_points=stage.teaching_points()[:2]
        )

    assert result.lattice is None
    assert result.cell_to_image is None
    assert any("fewer than three teaching points" in r.getMessage() for r in records)


def test_areas_with_different_steps_fall_back(thyra_logs) -> None:
    stage = _Stage()
    areas, positions, ref_px = _acquisition(stage)
    areas[1]["raster"] = [50.0, 50.0]
    with thyra_logs("thyra.alignment", logging.WARNING) as records:
        result = TeachingPointAlignment().compute_area_alignment(
            areas,
            positions,
            0,
            0,
            teaching_points=stage.teaching_points(),
            reference_point=tuple(ref_px),
        )

    assert result.lattice is None
    assert any("different raster steps" in r.getMessage() for r in records)


# --- Two public runs: flexImaging's own spot list is the answer -------------


def _spot_list(path: Path) -> Dict[Tuple[int, int, int], Tuple[float, float]]:
    """``(region, X, Y)`` -> stage position, from a flexImaging spot list."""
    spots = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        x, y, name = line.split()[:3]
        m = re.match(r"R(\d+)X(\d+)Y(\d+)", name)
        assert m, name
        spots[(int(m.group(1)), int(m.group(2)), int(m.group(3)))] = (
            float(x),
            float(y),
        )
    return spots


@pytest.mark.parametrize(
    ("stem", "node"),
    [
        ("20210920_vc_rugose_1mMTCA_gordon", (10, 6)),
        ("20210921_vc_rugose_tims_gordon", (9, 6)),
    ],
)
def test_public_runs_land_on_flexImagings_spot_list(
    stem: str, node: Tuple[int, int]
) -> None:
    mis = parse_mis_file(FIXTURES / f"{stem}.mis")
    spots = _spot_list(FIXTURES / f"{stem}_spot_list.txt")
    positions = [{"region": r, "raster_x": x, "raster_y": y} for r, x, y in spots]
    fx = min(p["raster_x"] for p in positions)
    fy = min(p["raster_y"] for p in positions)

    result = TeachingPointAlignment().compute_area_alignment(
        mis["areas"],
        positions,
        fx,
        fy,
        teaching_points=mis["teaching_points"],
        reference_point=tuple(mis["reference_point"]),
        raster_step=tuple(mis["raster"]),
    )

    assert result.lattice is not None
    assert result.lattice.reference_node == node
    assert result.lattice.spots_outside == 0
    # The teaching points map a stage position onto the photo.
    tp = mis["teaching_points"]
    image = np.array([[*t["image"], 1.0] for t in tp])
    stage = np.array([t["stage"] for t in tp], dtype=float)
    to_stage = np.linalg.solve(image, stage)
    to_image = np.linalg.inv(np.vstack([to_stage.T, [0.0, 0.0, 1.0]]))
    step_px = 1000.0 / float(np.hypot(to_stage[0, 0], to_stage[0, 1]))
    worst = 0.0
    for (r, x, y), (sx, sy) in spots.items():
        truth = (to_image @ np.array([sx, sy, 1.0]))[:2]
        got = np.array(result.transform_point(x - fx, y - fy))
        worst = max(worst, float(np.linalg.norm(got - truth)) / step_px)
    # Measured 0.0010 and 0.0014 raster steps (1 um); before D31, 0.57.
    assert worst < 0.005
