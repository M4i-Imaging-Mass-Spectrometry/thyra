# thyra/alignment/teaching_points.py
"""Teaching point alignment for FlexImaging optical-MSI registration.

This module handles the alignment between optical images and MSI data
using teaching point calibration data from FlexImaging .mis files.

Coordinate Systems:
- Image pixels: (x, y) in the optical reference image (origin at top-left)
- Stage coordinates (teaching): From teaching point calibration
- Stage coordinates (poslog): From acquisition position log
- MSI raster: (x, y) grid positions in the MSI dataset (0-based)

Key Challenge:
FlexImaging uses different coordinate frames for teaching (image calibration)
and acquisition (stage movement). These frames may have large offsets that
cannot be reliably computed without additional reference data.

The module provides:
1. Reliable image <-> teaching stage transformation from teaching points
2. Estimated MSI <-> image transformation (may need manual verification)
3. Methods to manually specify alignment offsets if automatic alignment fails
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .affine import AffineTransform

logger = logging.getLogger(__name__)


@dataclass
class TeachingPoint:
    """A single teaching point correspondence.

    Attributes:
        image_x: X coordinate in image pixels
        image_y: Y coordinate in image pixels
        stage_x: X coordinate in stage units (micrometers)
        stage_y: Y coordinate in stage units (micrometers)
    """

    image_x: int
    image_y: int
    stage_x: int
    stage_y: int

    @classmethod
    def from_dict(cls, data: Dict[str, Tuple[int, int]]) -> "TeachingPoint":
        """Create from parsed dictionary format.

        Args:
            data: Dictionary with 'image' and 'stage' tuples

        Returns:
            TeachingPoint instance
        """
        img = data["image"]
        stage = data["stage"]
        return cls(
            image_x=img[0],
            image_y=img[1],
            stage_x=stage[0],
            stage_y=stage[1],
        )


@dataclass
class RasterPosition:
    """A position from the poslog with raster and physical coordinates.

    Attributes:
        raster_x: Raster grid X coordinate
        raster_y: Raster grid Y coordinate
        phys_x: Physical stage X coordinate (micrometers)
        phys_y: Physical stage Y coordinate (micrometers)
    """

    raster_x: int
    raster_y: int
    phys_x: float
    phys_y: float


@dataclass
class AlignmentResult:
    """Result of teaching point alignment computation.

    Attributes:
        image_to_stage: Affine transform from image pixels to stage coords
        stage_to_image: Inverse transform (stage to image pixels)
        msi_to_image: Transform from MSI raster to image pixels
        image_to_msi: Transform from image pixels to MSI raster
        stage_offset: Estimated offset between teaching and poslog stages
        rmse: Root mean square error of teaching point fit
        warnings: List of alignment warnings
    """

    image_to_stage: AffineTransform
    stage_to_image: AffineTransform
    msi_to_image: Optional[AffineTransform] = None
    image_to_msi: Optional[AffineTransform] = None
    stage_offset: Optional[Tuple[float, float]] = None
    rmse: float = 0.0
    warnings: List[str] = field(default_factory=list)


@dataclass
class RegionMapping:
    """Mapping for a single acquisition region to image coordinates.

    The image coordinates use min/max bounds from the Area definition.
    The p1/p2 points in an Area are two corners of a bounding box - their
    order indicates scan direction but for coordinate mapping we use the
    actual min/max values to ensure monotonic mapping:
    - raster_min_x maps to image_min_x
    - raster_max_x maps to image_max_x
    - raster_min_y maps to image_min_y
    - raster_max_y maps to image_max_y

    Attributes:
        region_id: Region number (0, 1, 2, ...)
        name: Region name from Area definition
        raster_min_x: Minimum raster X in this region
        raster_max_x: Maximum raster X in this region
        raster_min_y: Minimum raster Y in this region
        raster_max_y: Maximum raster Y in this region
        image_min_x: Minimum image X for this region
        image_max_x: Maximum image X for this region
        image_min_y: Minimum image Y for this region
        image_max_y: Maximum image Y for this region
    """

    region_id: int
    name: str
    raster_min_x: int
    raster_max_x: int
    raster_min_y: int
    raster_max_y: int
    image_min_x: int
    image_max_x: int
    image_min_y: int
    image_max_y: int

    def _get_pixel_scale_x(self) -> float:
        """Get the X pixel scale (image pixels per raster step in X).

        Returns:
            Pixel scale along X axis
        """
        n_pixels = max(1, self.raster_max_x - self.raster_min_x + 1)
        image_width = self.image_max_x - self.image_min_x
        return image_width / n_pixels

    def _get_pixel_scale_y(self) -> float:
        """Get the Y pixel scale (image pixels per raster step in Y).

        Returns:
            Pixel scale along Y axis
        """
        n_pixels = max(1, self.raster_max_y - self.raster_min_y + 1)
        image_height = self.image_max_y - self.image_min_y
        return image_height / n_pixels

    def raster_to_image(self, raster_x: int, raster_y: int) -> Tuple[float, float]:
        """Convert raster coordinates to image pixel coordinates.

        Maps raster positions to pixel centers within the Area bounding box.
        The first pixel center is at image_min + half_pixel, the last at
        image_max - half_pixel, so pixels fill the box edge-to-edge.

        X and Y scales are computed independently so non-square aspect
        ratios are handled correctly.

        Args:
            raster_x: Original raster X coordinate (not normalized)
            raster_y: Original raster Y coordinate (not normalized)

        Returns:
            Tuple of (image_x, image_y) pixel center coordinates
        """
        scale_x = self._get_pixel_scale_x()
        scale_y = self._get_pixel_scale_y()

        # Pixel center = box_min + half_pixel + offset * scale
        image_x = (
            self.image_min_x + scale_x / 2.0 + (raster_x - self.raster_min_x) * scale_x
        )
        image_y = (
            self.image_min_y + scale_y / 2.0 + (raster_y - self.raster_min_y) * scale_y
        )

        return image_x, image_y

    def get_half_pixel_size(self) -> Tuple[float, float]:
        """Get the half-pixel size for this region in image coordinates.

        Returns independent X and Y half-sizes to handle non-square
        aspect ratios correctly.

        Returns:
            Tuple of (half_x, half_y) in image pixels
        """
        return self._get_pixel_scale_x() / 2, self._get_pixel_scale_y() / 2


#: Raster axis directions tried against the teaching-point frame, in order of
#: preference. FlexImaging's teaching frame has y pointing up the image while
#: raster Y counts down it, so (+1, -1) is what every measured run uses.
_LATTICE_AXIS_SIGNS: Tuple[Tuple[int, int], ...] = ((1, -1), (1, 1), (-1, -1), (-1, 1))

#: A spot this close to its Area's outline, in raster steps, still counts as
#: inside it. FlexImaging draws rectangles whose edges pass through the
#: outermost spots, so an exact test would reject them.
_AREA_EDGE_TOLERANCE_STEPS = 0.1

#: How far around the estimated reference node the integer search looks.
_NODE_SEARCH_RADIUS = 3


@dataclass(frozen=True)
class LatticeFit:
    """Where FlexImaging put each raster spot on its alignment image (D31).

    FlexImaging places every spot of a run on one lattice: one raster step
    apart in the teaching-point frame, through the ``.mis`` reference
    point. An Area outline only selects which nodes are measured, so the
    outline's size says nothing about the pitch.

    Attributes:
        cell_to_image: 3x3 affine from the TIC image's own coordinates to
            alignment-image pixels. Cell ``(i, j)`` spans ``[i, i + 1) x
            [j, j + 1)``, the SpatialData convention, so its centre
            ``(i + 0.5, j + 0.5)`` lands on the spot at raster
            ``(first_raster_x + i, first_raster_y + j)``.
        reference_node: Raster ``(X, Y)`` of the node on the reference point.
        axis_signs: Direction of raster X and Y in the teaching frame.
        raster_step_um: ``(x, y)`` lattice step in micrometres.
        spots_outside: Measured spots that fall outside their own Area.
        nodes_unmeasured: Lattice nodes inside an Area that were not measured.
        reference_from: ``"ReferencePoint"``, or ``"first teaching point"``
            when the ``.mis`` has none.
    """

    cell_to_image: np.ndarray
    reference_node: Tuple[int, int]
    axis_signs: Tuple[int, int]
    raster_step_um: Tuple[float, float]
    spots_outside: int
    nodes_unmeasured: int
    reference_from: str

    def cell_centre(self, norm_x: float, norm_y: float) -> Tuple[float, float]:
        """Image pixel of the spot in TIC cell ``(norm_x, norm_y)``."""
        x, y, _ = self.cell_to_image @ np.array([norm_x + 0.5, norm_y + 0.5, 1.0])
        return float(x), float(y)

    def cell_corners(self, norm_x: float, norm_y: float) -> List[Tuple[float, float]]:
        """The four image-pixel corners of TIC cell ``(norm_x, norm_y)``."""
        corners = np.array(
            [
                [norm_x, norm_y, 1.0],
                [norm_x + 1, norm_y, 1.0],
                [norm_x + 1, norm_y + 1, 1.0],
                [norm_x, norm_y + 1, 1.0],
            ]
        )
        out = corners @ self.cell_to_image.T
        return [(float(x), float(y)) for x, y, _ in out]


def _area_outline(area: Dict[str, Any]) -> Optional[List[Tuple[float, float]]]:
    """An Area's outline in image pixels: its polygon, or its rectangle."""
    points = area.get("points")
    if points and len(points) >= 3:
        return [(float(x), float(y)) for x, y in points]
    corners = points if points and len(points) == 2 else None
    if corners is None and "p1" in area and "p2" in area:
        corners = [area["p1"], area["p2"]]
    if corners is None:
        return None
    (x0, y0), (x1, y1) = corners
    return [
        (float(x0), float(y0)),
        (float(x1), float(y0)),
        (float(x1), float(y1)),
        (float(x0), float(y1)),
    ]


class _LatticeSearch:
    """The whole-node search behind :func:`fit_raster_lattice`.

    Holds the teaching-point affine, the reference point in stage
    micrometres, the measured spots and each measured region's Area, and
    scores a candidate ``(node, axis signs)`` by how many spots fall outside
    their Area and how many Area nodes went unmeasured.
    """

    def __init__(
        self,
        to_stage: np.ndarray,
        ref: np.ndarray,
        step: Tuple[float, float],
        positions: List[Dict[str, Any]],
        areas: List[Dict[str, Any]],
        tolerance_px: float,
    ) -> None:
        from shapely.geometry import Polygon

        self.to_stage = to_stage
        self.to_image = np.linalg.inv(to_stage)
        self.ref = ref
        self.step_x, self.step_y = step
        self.raster_x = np.array([p["raster_x"] for p in positions], dtype=np.float64)
        self.raster_y = np.array([p["raster_y"] for p in positions], dtype=np.float64)
        self.regions = np.array([p["region"] for p in positions], dtype=np.int64)
        # Each measured region's Area as drawn, and padded by the edge
        # tolerance for the inside test.
        self.exact: Dict[int, Any] = {}
        self.padded: Dict[int, Any] = {}
        for region in np.unique(self.regions).tolist():
            if 0 <= region < len(areas):
                vertices = _area_outline(areas[region])
                if vertices is not None:
                    self.exact[region] = Polygon(vertices)
                    self.padded[region] = self.exact[region].buffer(tolerance_px)

    def to_px(
        self, node: Tuple[int, int], signs: Tuple[int, int], xs: Any, ys: Any
    ) -> Tuple[Any, Any]:
        """Image pixels of raster nodes ``(xs, ys)`` on the lattice."""
        stage_x = self.ref[0] + signs[0] * self.step_x * (xs - node[0])
        stage_y = self.ref[1] + signs[1] * self.step_y * (ys - node[1])
        m = self.to_image
        return (
            m[0, 0] * stage_x + m[0, 1] * stage_y + m[0, 2],
            m[1, 0] * stage_x + m[1, 1] * stage_y + m[1, 2],
        )

    def outside(self, node: Tuple[int, int], signs: Tuple[int, int]) -> int:
        """Measured spots that fall outside their own (padded) Area."""
        import shapely

        count = 0
        for region, polygon in self.padded.items():
            sel = self.regions == region
            px, py = self.to_px(node, signs, self.raster_x[sel], self.raster_y[sel])
            count += int((~shapely.contains_xy(polygon, px, py)).sum())
        return count

    def unmeasured(self, node: Tuple[int, int], signs: Tuple[int, int]) -> int:
        """Lattice nodes inside an Area as drawn that no spot was measured on.

        Against the outline as drawn: FlexImaging measures the nodes inside
        it, and the edge tolerance would count nodes it never meant to.
        """
        import shapely

        count = 0
        for region, polygon in self.exact.items():
            sel = self.regions == region
            measured = set(
                zip(self.raster_x[sel].astype(int), self.raster_y[sel].astype(int))
            )
            # The Area's outline in raster coordinates bounds the nodes to try.
            verts = np.array(polygon.exterior.coords)
            stage = verts @ self.to_stage[:2, :2].T + self.to_stage[:2, 2]
            gx = node[0] + (stage[:, 0] - self.ref[0]) / (signs[0] * self.step_x)
            gy = node[1] + (stage[:, 1] - self.ref[1]) / (signs[1] * self.step_y)
            xs, ys = np.meshgrid(
                np.arange(np.floor(gx.min()), np.ceil(gx.max()) + 1),
                np.arange(np.floor(gy.min()), np.ceil(gy.max()) + 1),
            )
            xs, ys = xs.ravel(), ys.ravel()
            px, py = self.to_px(node, signs, xs, ys)
            inside = shapely.contains_xy(polygon, px, py)
            nodes = zip(xs[inside].astype(int), ys[inside].astype(int))
            count += sum(1 for xy in nodes if xy not in measured)
        return count

    def candidates(self, signs: Tuple[int, int]) -> List[Tuple[int, int]]:
        """Nodes near where each region's centre puts the reference node."""
        guesses = set()
        for region, polygon in self.padded.items():
            sel = self.regions == region
            centre = self.to_stage @ np.array(
                [polygon.centroid.x, polygon.centroid.y, 1.0]
            )
            gx = self.raster_x[sel].mean() - (centre[0] - self.ref[0]) / (
                signs[0] * self.step_x
            )
            gy = self.raster_y[sel].mean() - (centre[1] - self.ref[1]) / (
                signs[1] * self.step_y
            )
            guesses.add((int(round(gx)), int(round(gy))))
        span = range(-_NODE_SEARCH_RADIUS, _NODE_SEARCH_RADIUS + 1)
        return sorted(
            {(gx + dx, gy + dy) for gx, gy in guesses for dx in span for dy in span}
        )

    def best(self) -> Tuple[int, int, Tuple[int, int], Tuple[int, int]]:
        """``(spots outside, unmeasured nodes, node, signs)`` of the best fit.

        Fewest spots outside first, then fewest unmeasured nodes, then the
        usual axis directions.
        """
        scored = [
            (self.outside(node, signs), preference, node, signs)
            for preference, signs in enumerate(_LATTICE_AXIS_SIGNS)
            for node in self.candidates(signs)
        ]
        fewest = min(s[0] for s in scored)
        n_unmeasured, _, node, signs, n_outside = min(
            (self.unmeasured(node, signs), preference, node, signs, n_out)
            for n_out, preference, node, signs in scored
            if n_out == fewest
        )
        return n_outside, n_unmeasured, node, signs


def fit_raster_lattice(
    teaching_points: List[Dict[str, Any]],
    reference_point: Optional[Tuple[float, float]],
    raster_step: Tuple[float, float],
    areas: List[Dict[str, Any]],
    positions: List[Dict[str, Any]],
    first_raster_x: int,
    first_raster_y: int,
) -> Optional[LatticeFit]:
    """Place the run's raster on the alignment image the way FlexImaging does.

    The teaching points give the image-to-stage affine. The lattice runs
    through the reference point with the raster step. Which raster node
    sits on the reference point is a whole-number question: it is the node
    that puts every measured spot inside its own Area, and of those the one
    that leaves the fewest unmeasured nodes inside the Areas. On the public
    MassIVE MSV000088438 runs this reproduces flexImaging's own spot list
    within 1 um.

    Args:
        teaching_points: ``.mis`` teaching points, ``{"image", "stage"}``.
        reference_point: ``.mis`` ``<ReferencePoint>`` in image pixels, or
            ``None`` to use the first teaching point (FlexImaging puts it
            there in every file seen).
        raster_step: ``(x, y)`` raster step in micrometres.
        areas: ``.mis`` Areas, in file order. Region ``N`` is Area ``N``.
        positions: Measured spots, ``{"raster_x", "raster_y", "region"}``.
        first_raster_x: Raster X of TIC cell 0.
        first_raster_y: Raster Y of TIC cell 0.

    Returns:
        The fit, or ``None`` when the inputs cannot place a lattice: fewer
        than three teaching points, no positive step, no Area that any
        measured region maps to.
    """
    step_x, step_y = (float(raster_step[0]), float(raster_step[1]))
    if len(teaching_points) < 3 or step_x <= 0 or step_y <= 0 or not positions:
        return None

    points = [TeachingPoint.from_dict(tp) for tp in teaching_points]
    image_to_stage = AffineTransform.from_points(
        [(float(p.image_x), float(p.image_y)) for p in points],
        [(float(p.stage_x), float(p.stage_y)) for p in points],
    )
    if reference_point is None:
        reference_point = (float(points[0].image_x), float(points[0].image_y))
        reference_from = "first teaching point"
    else:
        reference_from = "ReferencePoint"
    to_stage = image_to_stage.matrix
    ref = to_stage @ np.array([reference_point[0], reference_point[1], 1.0])

    # The edge tolerance in image pixels: steps times pixels per micrometre.
    um_per_px = 0.5 * (image_to_stage.scale_x + image_to_stage.scale_y)
    tolerance_px = _AREA_EDGE_TOLERANCE_STEPS * min(step_x, step_y) / um_per_px
    search = _LatticeSearch(
        to_stage, ref, (step_x, step_y), positions, areas, tolerance_px
    )
    if not search.padded:
        return None
    n_outside, n_unmeasured, node, signs = search.best()

    sx, sy = signs
    lattice = np.array(
        [
            [sx * step_x, 0.0, ref[0] + sx * step_x * (first_raster_x - node[0] - 0.5)],
            [0.0, sy * step_y, ref[1] + sy * step_y * (first_raster_y - node[1] - 0.5)],
            [0.0, 0.0, 1.0],
        ]
    )
    fit = LatticeFit(
        cell_to_image=search.to_image @ lattice,
        reference_node=(int(node[0]), int(node[1])),
        axis_signs=(int(sx), int(sy)),
        raster_step_um=(step_x, step_y),
        spots_outside=int(n_outside),
        nodes_unmeasured=int(n_unmeasured),
        reference_from=reference_from,
    )
    logger.info(
        f"FlexImaging lattice: reference node {fit.reference_node} "
        f"(from {reference_from}), step {step_x:g} x {step_y:g} um, "
        f"axes {fit.axis_signs}, {n_outside} of {len(positions)} spots outside "
        f"their Area, {n_unmeasured} unmeasured nodes inside the Areas"
    )
    if n_outside:
        logger.warning(
            f"{n_outside} of {len(positions)} measured spots fall outside their "
            f".mis Area on the best lattice; check that the .mis belongs to this "
            f"run. The optical alignment may be off by whole raster steps."
        )
    if signs != _LATTICE_AXIS_SIGNS[0]:
        logger.warning(
            f"FlexImaging raster axes {signs} differ from the usual (1, -1); "
            "the optical alignment follows the Areas."
        )
    return fit


@dataclass
class AreaAlignmentResult:
    """Result of area-based alignment computation.

    This is the preferred alignment method when Area definitions are available
    in the .mis file, as it provides direct region-to-image mapping.

    Attributes:
        region_mappings: List of RegionMapping objects for each region
        first_raster_x: First raster X offset from header
        first_raster_y: First raster Y offset from header
        pos_to_region: Mapping from (raster_x, raster_y) to region_id
        lattice: The teaching-point lattice (D31). ``None`` only when the
            ``.mis`` lacks what it needs, and then each region is stretched
            over its Area's bounding box instead, which is approximate.
    """

    region_mappings: List[RegionMapping]
    first_raster_x: int
    first_raster_y: int
    pos_to_region: Dict[Tuple[int, int], int] = field(default_factory=dict)
    lattice: Optional[LatticeFit] = None

    @property
    def cell_to_image(self) -> Optional[np.ndarray]:
        """The lattice's TIC-cell-to-image affine, or ``None`` without one."""
        return self.lattice.cell_to_image if self.lattice is not None else None

    def transform_point(
        self, norm_x: int, norm_y: int
    ) -> Optional[Tuple[float, float]]:
        """Transform normalized raster coordinates to image coordinates.

        Args:
            norm_x: Normalized (0-based) raster X coordinate
            norm_y: Normalized (0-based) raster Y coordinate

        Returns:
            Tuple of (image_x, image_y) or None if no mapping exists
        """
        # Convert normalized to original raster coords
        orig_x = norm_x + self.first_raster_x
        orig_y = norm_y + self.first_raster_y

        # Find which region this belongs to
        region_id = self.pos_to_region.get((orig_x, orig_y))
        if region_id is None:
            return None

        if self.lattice is not None:
            return self.lattice.cell_centre(norm_x, norm_y)

        # Find the region mapping
        for mapping in self.region_mappings:
            if mapping.region_id == region_id:
                return mapping.raster_to_image(orig_x, orig_y)

        return None

    def cell_corners(
        self, norm_x: int, norm_y: int
    ) -> Optional[List[Tuple[float, float]]]:
        """The image-pixel outline of a measured pixel, or ``None``.

        Under the lattice this is the exact cell the TIC image draws, a
        parallelogram when the photo is rotated against the stage. Without
        one it is the axis-aligned box of the bounding-box stretch.
        """
        if (norm_x + self.first_raster_x, norm_y + self.first_raster_y) not in (
            self.pos_to_region
        ):
            return None
        if self.lattice is not None:
            return self.lattice.cell_corners(norm_x, norm_y)
        centre = self.transform_point(norm_x, norm_y)
        half = self.get_half_pixel_size(norm_x, norm_y)
        if centre is None or half is None:
            return None
        (cx, cy), (hx, hy) = centre, half
        return [
            (cx - hx, cy - hy),
            (cx + hx, cy - hy),
            (cx + hx, cy + hy),
            (cx - hx, cy + hy),
        ]

    def get_half_pixel_size(
        self, norm_x: int, norm_y: int
    ) -> Optional[Tuple[float, float]]:
        """Get the half-pixel size for a given position.

        Returns independent X and Y half-sizes to handle non-square
        aspect ratios. The scale may vary between regions due to
        Area definitions.

        Args:
            norm_x: Normalized (0-based) raster X coordinate
            norm_y: Normalized (0-based) raster Y coordinate

        Returns:
            Tuple of (half_x, half_y) in image pixels, or None if
            no mapping exists
        """
        # Convert normalized to original raster coords
        orig_x = norm_x + self.first_raster_x
        orig_y = norm_y + self.first_raster_y

        # Find which region this belongs to
        region_id = self.pos_to_region.get((orig_x, orig_y))
        if region_id is None:
            return None

        if self.lattice is not None:
            m = self.lattice.cell_to_image
            return (
                float(np.hypot(m[0, 0], m[1, 0])) / 2.0,
                float(np.hypot(m[0, 1], m[1, 1])) / 2.0,
            )

        # Find the region mapping and get its half-pixel size
        for mapping in self.region_mappings:
            if mapping.region_id == region_id:
                return mapping.get_half_pixel_size()

        return None


class TeachingPointAlignment:
    """Computes alignment between optical images and MSI data.

    This class handles the coordinate system transformations needed to
    align optical images with MSI raster data using FlexImaging teaching
    points.

    The workflow:
    1. Parse teaching points from .mis metadata
    2. Compute image -> stage affine transformation
    3. Determine offset between teaching stage and poslog stage coords
    4. Compute final image -> MSI raster transformation

    Example:
        >>> aligner = TeachingPointAlignment()
        >>> result = aligner.compute_alignment(
        ...     teaching_points=reader.mis_metadata['teaching_points'],
        ...     poslog_positions=reader._positions,
        ...     raster_step=(20.0, 20.0),
        ... )
        >>> # Transform optical image coordinates to MSI raster
        >>> msi_coords = result.image_to_msi.transform_point(img_x, img_y)
    """

    def __init__(self):
        """Initialize the alignment calculator."""
        pass

    def compute_alignment(
        self,
        teaching_points: List[Dict[str, Tuple[int, int]]],
        poslog_positions: Optional[List[Dict[str, Any]]] = None,
        raster_step: Tuple[float, float] = (20.0, 20.0),
        raster_offset: Tuple[int, int] = (0, 0),
        flip_poslog_x: bool = False,
        flip_poslog_y: bool = False,
    ) -> AlignmentResult:
        """Compute alignment transformations from teaching points.

        Args:
            teaching_points: List of teaching point dictionaries with
                'image' and 'stage' keys containing (x, y) tuples
            poslog_positions: Optional list of position dictionaries from
                poslog parsing, used to estimate stage coordinate offset
            raster_step: (step_x, step_y) raster step size in micrometers
            raster_offset: (offset_x, offset_y) offset of first raster position
            flip_poslog_x: If True, negate poslog X coordinates (for inverted
                stage X-axis relative to teaching points)
            flip_poslog_y: If True, negate poslog Y coordinates (for inverted
                stage Y-axis relative to teaching points)

        Returns:
            AlignmentResult with computed transformations
        """
        warnings: List[str] = []

        # Parse and validate teaching points
        if len(teaching_points) < 3:
            n_pts = len(teaching_points)
            raise ValueError(f"At least 3 teaching points required, got {n_pts}")

        points = [TeachingPoint.from_dict(tp) for tp in teaching_points]
        logger.info(f"Processing {len(points)} teaching points")

        # Compute transformations and RMSE
        image_to_stage = self._compute_image_to_stage(points)
        rmse = self._compute_alignment_rmse(points, image_to_stage)

        if rmse > 10.0:  # More than 10 um error
            warnings.append(
                f"Teaching point fit has high error (RMSE={rmse:.2f} um). "
                "Alignment may be inaccurate."
            )

        stage_to_image = image_to_stage.inverse()

        # Compute MSI transforms if poslog available
        msi_result = self._compute_msi_alignment(
            points,
            poslog_positions,
            image_to_stage,
            stage_to_image,
            raster_step,
            flip_poslog_x,
            flip_poslog_y,
        )

        if msi_result["warning"]:
            warnings.append(msi_result["warning"])

        return AlignmentResult(
            image_to_stage=image_to_stage,
            stage_to_image=stage_to_image,
            msi_to_image=msi_result["msi_to_image"],
            image_to_msi=msi_result["image_to_msi"],
            stage_offset=msi_result["stage_offset"],
            rmse=rmse,
            warnings=warnings,
        )

    def _compute_alignment_rmse(
        self, points: List[TeachingPoint], transform: AffineTransform
    ) -> float:
        """Compute RMSE for alignment quality assessment."""
        src_pts = [(float(p.image_x), float(p.image_y)) for p in points]
        dst_pts = [(float(p.stage_x), float(p.stage_y)) for p in points]
        rmse, _ = transform.compute_residuals(src_pts, dst_pts)
        logger.info(f"Image->Stage RMSE: {rmse:.4f} um")
        return rmse

    def _apply_coordinate_flips(
        self,
        positions: List[Dict[str, Any]],
        flip_x: bool,
        flip_y: bool,
    ) -> List[Dict[str, Any]]:
        """Apply coordinate flips to poslog positions if needed."""
        if not flip_x and not flip_y:
            return positions

        flipped = []
        for pos in positions:
            flipped_pos = pos.copy()
            if flip_x:
                flipped_pos["phys_x"] = -pos["phys_x"]
            if flip_y:
                flipped_pos["phys_y"] = -pos["phys_y"]
            flipped.append(flipped_pos)
        return flipped

    def _compute_msi_alignment(
        self,
        points: List[TeachingPoint],
        poslog_positions: Optional[List[Dict[str, Any]]],
        image_to_stage: AffineTransform,
        stage_to_image: AffineTransform,
        raster_step: Tuple[float, float],
        flip_poslog_x: bool,
        flip_poslog_y: bool,
    ) -> Dict[str, Any]:
        """Compute MSI-to-image alignment from poslog positions."""
        result: Dict[str, Any] = {
            "stage_offset": None,
            "msi_to_image": None,
            "image_to_msi": None,
            "warning": None,
        }

        if not poslog_positions:
            return result

        flipped_positions = self._apply_coordinate_flips(
            poslog_positions, flip_poslog_x, flip_poslog_y
        )

        stage_offset = self._estimate_stage_offset(
            points, flipped_positions, raster_step
        )
        result["stage_offset"] = stage_offset

        if stage_offset is None:
            result["warning"] = (
                "Could not determine stage coordinate offset. "
                "MSI-to-image transform may require manual calibration."
            )
            return result

        logger.info(
            f"Estimated stage offset: ({stage_offset[0]:.1f}, "
            f"{stage_offset[1]:.1f}) um"
        )

        first_phys = (
            float(flipped_positions[0]["phys_x"]),
            float(flipped_positions[0]["phys_y"]),
        )

        msi_to_image, image_to_msi = self._compute_msi_transforms(
            image_to_stage,
            stage_to_image,
            stage_offset,
            raster_step,
            first_phys,
            flip_poslog_x,
        )
        result["msi_to_image"] = msi_to_image
        result["image_to_msi"] = image_to_msi

        return result

    def _compute_image_to_stage(self, points: List[TeachingPoint]) -> AffineTransform:
        """Compute affine transformation from image pixels to stage coords.

        Args:
            points: List of teaching points

        Returns:
            AffineTransform from image to stage coordinates
        """
        src_points = [(float(p.image_x), float(p.image_y)) for p in points]
        dst_points = [(float(p.stage_x), float(p.stage_y)) for p in points]

        return AffineTransform.from_points(src_points, dst_points)

    def _estimate_stage_offset(
        self,
        teaching_points: List[TeachingPoint],
        poslog_positions: List[Dict[str, Any]],
        raster_step: Tuple[float, float],
    ) -> Optional[Tuple[float, float]]:
        """Estimate offset between teaching stage and poslog stage coords.

        This computes the translation offset between the two coordinate
        systems by comparing their centers. The assumption is that the
        MSI acquisition region is approximately centered on the optical
        image region defined by the teaching points.

        Coordinate relationship:
        - teaching_stage = poslog_physical - offset
        - poslog_physical = teaching_stage + offset

        Args:
            teaching_points: Teaching point data
            poslog_positions: Position log entries
            raster_step: Raster step size (step_x, step_y) in um

        Returns:
            Estimated (offset_x, offset_y) or None if cannot determine
        """
        if not poslog_positions:
            return None

        # Extract poslog physical coordinates
        phys_x = np.array([pos["phys_x"] for pos in poslog_positions])
        phys_y = np.array([pos["phys_y"] for pos in poslog_positions])

        # Get teaching point stage coordinate center
        teaching_x = [p.stage_x for p in teaching_points]
        teaching_y = [p.stage_y for p in teaching_points]
        teaching_center_x = float(np.mean(teaching_x))
        teaching_center_y = float(np.mean(teaching_y))

        # Get poslog physical coordinate center
        poslog_center_x = float(np.mean(phys_x))
        poslog_center_y = float(np.mean(phys_y))

        logger.debug(
            f"Teaching stage center: "
            f"({teaching_center_x:.1f}, {teaching_center_y:.1f})"
        )
        logger.debug(
            f"Poslog physical center: "
            f"({poslog_center_x:.1f}, {poslog_center_y:.1f})"
        )

        # The offset is the difference between poslog physical coords
        # and teaching stage coords for the same physical location
        # offset = poslog_physical - teaching_stage
        offset_x = poslog_center_x - teaching_center_x
        offset_y = poslog_center_y - teaching_center_y

        logger.info(
            f"Stage coordinate offset (poslog - teaching): "
            f"({offset_x:.1f}, {offset_y:.1f}) um"
        )

        return (offset_x, offset_y)

    def _compute_msi_transforms(
        self,
        image_to_stage: AffineTransform,
        stage_to_image: AffineTransform,
        stage_offset: Tuple[float, float],
        raster_step: Tuple[float, float],
        first_phys: Tuple[float, float],
        flip_poslog_x: bool = False,
    ) -> Tuple[AffineTransform, AffineTransform]:
        """Compute transformations between MSI raster and image coordinates.

        The chain of transformations:
        MSI raster (0-based) -> poslog physical -> teaching stage -> image

        The poslog Y-axis is inverted: larger raster Y maps to smaller
        physical Y (step_y is effectively negative).

        Args:
            image_to_stage: Transform from image pixels to teaching stage
            stage_to_image: Inverse transform
            stage_offset: (offset_x, offset_y) where offset = poslog - teaching
            raster_step: (step_x, step_y) raster step in um (magnitudes)
            first_phys: (phys_x, phys_y) physical position at raster (0, 0)
            flip_poslog_x: If True, negate X scale (poslog X is inverted)

        Returns:
            Tuple of (msi_to_image, image_to_msi) transforms
        """
        step_x, step_y = raster_step
        offset_x, offset_y = stage_offset
        first_phys_x, first_phys_y = first_phys

        # Determine X scale sign based on flip
        # If flip_poslog_x is True, the poslog coordinates were negated,
        # so the raster step direction is also reversed
        x_scale = -step_x if flip_poslog_x else step_x

        # Normalized raster (0-based) to teaching stage coordinates:
        # 1. normalized_raster -> physical:
        #    phys_x = x_scale * norm_x + first_phys_x
        #    phys_y = -step_y * norm_y + first_phys_y  (Y is inverted!)
        # 2. physical -> teaching:
        #    teaching = physical - offset
        #
        # Combined: normalized raster -> teaching stage
        #    teaching_x = x_scale * norm_x + first_phys_x - offset_x
        #    teaching_y = -step_y * norm_y + first_phys_y - offset_y

        tx = first_phys_x - offset_x
        ty = first_phys_y - offset_y

        logger.debug(
            f"MSI transform: scale=({x_scale}, {-step_y}), "
            f"translation=({tx:.1f}, {ty:.1f})"
        )

        # Build MSI -> teaching stage transform
        # Note: step_y is negated because Y-axis is inverted in poslog
        msi_to_teaching = AffineTransform.from_scale_translate(
            scale_x=x_scale,
            scale_y=-step_y,  # Negative because poslog Y is inverted
            tx=tx,
            ty=ty,
        )

        # Chain: MSI -> teaching stage -> image
        msi_to_image = msi_to_teaching.compose(stage_to_image)

        # Inverse: image -> MSI
        image_to_msi = msi_to_image.inverse()

        return msi_to_image, image_to_msi

    def validate_alignment(
        self,
        result: AlignmentResult,
        image_shape: Tuple[int, int],
        msi_shape: Tuple[int, int],
    ) -> List[str]:
        """Validate alignment by checking if coordinates map sensibly.

        Args:
            result: AlignmentResult to validate
            image_shape: (height, width) of optical image
            msi_shape: (height, width) of MSI raster

        Returns:
            List of validation warnings (empty if OK)
        """
        warnings = []
        img_h, img_w = image_shape
        msi_h, msi_w = msi_shape

        if result.msi_to_image is None:
            warnings.append("MSI-to-image transform not available")
            return warnings

        # Check corners of MSI raster map to within image bounds
        corners = [
            (0, 0),
            (msi_w - 1, 0),
            (0, msi_h - 1),
            (msi_w - 1, msi_h - 1),
        ]

        for cx, cy in corners:
            ix, iy = result.msi_to_image.transform_point(cx, cy)

            # Allow some margin outside image bounds
            margin = 0.2  # 20% margin
            if ix < -img_w * margin or ix > img_w * (1 + margin):
                warnings.append(
                    f"MSI corner ({cx}, {cy}) maps to image X={ix:.0f}, "
                    f"outside valid range [0, {img_w}]"
                )
            if iy < -img_h * margin or iy > img_h * (1 + margin):
                warnings.append(
                    f"MSI corner ({cx}, {cy}) maps to image Y={iy:.0f}, "
                    f"outside valid range [0, {img_h}]"
                )

        return warnings

    def compute_area_alignment(
        self,
        areas: List[Dict[str, Any]],
        poslog_positions: List[Dict[str, Any]],
        first_raster_x: int,
        first_raster_y: int,
        teaching_points: Optional[List[Dict[str, Any]]] = None,
        reference_point: Optional[Tuple[float, float]] = None,
        raster_step: Optional[Tuple[float, float]] = None,
    ) -> AreaAlignmentResult:
        """Compute alignment using Area definitions from .mis file.

        With three teaching points and a raster step, every spot is placed
        on FlexImaging's lattice (:func:`fit_raster_lattice`, D31). Without
        them, each region is stretched over its Area's bounding box, which
        is approximate: the outline is drawn by hand and is not the extent
        of the spots.

        Args:
            areas: List of Area dictionaries with 'name', 'p1', 'p2' keys
                where p1 and p2 are (x, y) tuples of image pixel coordinates,
                and optionally 'points', 'type' and 'raster'
            poslog_positions: List of position dictionaries from poslog parsing
            first_raster_x: First raster X offset from header
            first_raster_y: First raster Y offset from header
            teaching_points: ``.mis`` teaching points
            reference_point: ``.mis`` ``<ReferencePoint>`` in image pixels
            raster_step: Raster step in micrometres, used when the Areas
                state none of their own

        Returns:
            AreaAlignmentResult with region mappings and coordinate transform
        """
        # Group positions by region
        regions: Dict[int, Dict[str, Any]] = {}
        for pos in poslog_positions:
            r = pos["region"]
            if r not in regions:
                regions[r] = {
                    "positions": [],
                    "raster_xs": [],
                    "raster_ys": [],
                }
            regions[r]["positions"].append(pos)
            regions[r]["raster_xs"].append(pos["raster_x"])
            regions[r]["raster_ys"].append(pos["raster_y"])

        # Compute bounds for each region
        for r in regions:
            rx = regions[r]["raster_xs"]
            ry = regions[r]["raster_ys"]
            regions[r]["min_x"] = min(rx)
            regions[r]["max_x"] = max(rx)
            regions[r]["min_y"] = min(ry)
            regions[r]["max_y"] = max(ry)

        # Build position to region lookup
        pos_to_region: Dict[Tuple[int, int], int] = {}
        for r, data in regions.items():
            for pos in data["positions"]:
                pos_to_region[(pos["raster_x"], pos["raster_y"])] = r

        # Match areas to regions by region number.
        # Region number N maps to the Nth Area (0-indexed) from the .mis file.
        # For single-region .d files, the SpotName-derived region number
        # determines which Area the data aligns to.
        region_mappings: List[RegionMapping] = []
        for i, area in enumerate(areas):
            if i not in regions:
                logger.debug(
                    f"Area '{area['name']}' (index {i}) has no matching "
                    f"region in data (available regions: "
                    f"{sorted(regions.keys())})"
                )
                continue

            # Area bounds in image pixels - use min/max of p1 and p2
            # p1 and p2 are two corners of the bounding box, their order
            # indicates scan direction but we use min/max for mapping
            p1_x, p1_y = area["p1"]
            p2_x, p2_y = area["p2"]

            mapping = RegionMapping(
                region_id=i,
                name=area["name"],
                raster_min_x=regions[i]["min_x"],
                raster_max_x=regions[i]["max_x"],
                raster_min_y=regions[i]["min_y"],
                raster_max_y=regions[i]["max_y"],
                image_min_x=min(p1_x, p2_x),
                image_max_x=max(p1_x, p2_x),
                image_min_y=min(p1_y, p2_y),
                image_max_y=max(p1_y, p2_y),
            )
            region_mappings.append(mapping)

            logger.info(
                f"Region {i} -> Area '{area['name']}': "
                f"raster ({mapping.raster_min_x}, {mapping.raster_min_y}) to "
                f"({mapping.raster_max_x}, {mapping.raster_max_y}), "
                f"image ({mapping.image_min_x}, {mapping.image_min_y}) to "
                f"({mapping.image_max_x}, {mapping.image_max_y})"
            )

        lattice = self._fit_lattice(
            areas,
            poslog_positions,
            regions,
            first_raster_x,
            first_raster_y,
            teaching_points,
            reference_point,
            raster_step,
        )

        return AreaAlignmentResult(
            region_mappings=region_mappings,
            first_raster_x=first_raster_x,
            first_raster_y=first_raster_y,
            pos_to_region=pos_to_region,
            lattice=lattice,
        )

    @staticmethod
    def _fit_lattice(
        areas: List[Dict[str, Any]],
        positions: List[Dict[str, Any]],
        regions: Dict[int, Dict[str, Any]],
        first_raster_x: int,
        first_raster_y: int,
        teaching_points: Optional[List[Dict[str, Any]]],
        reference_point: Optional[Tuple[float, float]],
        raster_step: Optional[Tuple[float, float]],
    ) -> Optional[LatticeFit]:
        """The lattice fit, or ``None`` (with the reason logged) to fall back."""
        if not teaching_points or len(teaching_points) < 3:
            logger.warning(
                "The .mis has fewer than three teaching points; optical "
                "alignment stretches each region over its Area outline, which "
                "can be off by up to a raster step."
            )
            return None
        steps = {
            tuple(float(v) for v in areas[r]["raster"])
            for r in regions
            if 0 <= r < len(areas) and areas[r].get("raster")
        }
        if len(steps) > 1:
            logger.warning(
                f"The measured Areas use different raster steps {sorted(steps)}; "
                "one lattice cannot place them, so optical alignment stretches "
                "each region over its Area outline instead."
            )
            return None
        step = next(iter(steps)) if steps else raster_step
        if step is None:
            logger.warning(
                "No raster step for the optical alignment lattice; it stretches "
                "each region over its Area outline instead."
            )
            return None
        return fit_raster_lattice(
            teaching_points=teaching_points,
            reference_point=reference_point,
            raster_step=(float(step[0]), float(step[1])),
            areas=areas,
            positions=positions,
            first_raster_x=first_raster_x,
            first_raster_y=first_raster_y,
        )
