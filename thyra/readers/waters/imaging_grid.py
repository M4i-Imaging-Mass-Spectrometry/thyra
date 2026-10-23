# thyra/readers/waters/imaging_grid.py
"""Imaging grid reconstruction from Waters laser position metadata.

Translates the logic from mzmine's ImagingMetadata.java to reconstruct
a pixel coordinate grid from the laser X/Y positions stored in each scan's
ScanInfo struct. Waters .raw imaging files do not store grid dimensions
explicitly -- they must be derived from the set of unique laser positions.

The derivation fits a *lattice* to those positions rather than ranking the
distinct ones. Ranking is exact only for a raster with every row and column
present: one missing interior row shifts every row below it up by one and
inflates the pitch, because the extent still spans the missing interval
while the interval count does not (issue #227). Fitting origin, pitch and
count, then indexing by ``round((pos - origin) / pitch)``, is right in both
cases and is what lets a function be tested against the raster it claims to
belong to (issue #232).
"""

import logging
import statistics
from dataclasses import dataclass, field
from functools import cached_property
from math import isfinite
from typing import Dict, Iterable, List, Optional, Set, Tuple

from ...errors import ConversionRefused
from .instrument import DeclaredRegion, StagePattern
from .masslynx_lib import NO_POSITION, FunctionType, MassLynxLib, ScanInfoData

logger = logging.getLogger(__name__)

#: How far off its lattice point a stage reading may sit, as a fraction of
#: the pitch, before the raster is not a raster. A quarter pitch is far
#: looser than any real stage -- the four registry acquisitions land on
#: exact multiples of their step -- and still refuses the wobble that used
#: to explode a 20x10 grid into 187x171 (issue #231).
_LATTICE_TOLERANCE = 0.25

#: The fraction of a fitted axis's lattice points that must carry a stage
#: reading. A raster missing a row or two still scores far above this; a
#: pitch mis-estimated from sub-pitch jitter scores a few percent, because
#: the fitted pitch is then the wobble and the lattice has orders of
#: magnitude more points than the acquisition has columns. Held at one half
#: rather than something tighter because a legitimately ragged acquisition
#: -- an aborted raster, a non-rectangular region -- still projects onto
#: nearly every value of each axis, so per-axis occupancy stays high even
#: when the image itself is sparse.
_MIN_LATTICE_OCCUPANCY = 0.5

#: How far a declared region's box is widened on each side when a scan is
#: assigned to it, as a fraction of the axis's declared step. A region's
#: length can be smaller than one step (a line scan declares a 10 um length
#: at a 20 um step), so the stage reading would otherwise sit outside.
_BOX_MARGIN_STEPS = 0.5

#: The widening, in um, for an axis whose step the method does not declare.
_BOX_MARGIN_FALLBACK_UM = 0.5

#: Regions whose own pitches differ by more than this fraction of their
#: median were not rastered on one pitch, and are not placed on one grid.
_REGION_PITCH_AGREEMENT = 0.01

#: A region moved by less than this, in um, is on the lattice for every
#: practical purpose and is not reported.
_REPORTABLE_SHIFT_UM = 0.01

_STAGE_POSITIONS = "stage positions"
_METHOD_STEP = "method step"
_NO_SOURCE = "none"


@dataclass(frozen=True)
class AxisLattice:
    """The regular grid one stage axis was rastered on, in um.

    ``pitch`` is 0.0 and ``count`` is 1 for an axis with a single position:
    a single-line scan has no interval to measure on its short axis. The
    grid then takes that axis's pixel size from the method's declared step
    when it has one (:attr:`ImagingGrid.pixel_size_source`); without one the
    pixel size stays 0.0, which waters_extractor turns into
    ``pixel_size=None`` and so into the actionable "use --pixel-size"
    refusal.
    """

    origin: float
    pitch: float
    count: int

    def index(self, position: float) -> int:
        """The 0-based lattice index a stage reading belongs to."""
        if self.pitch <= 0.0:
            return 0
        return int(round((position - self.origin) / self.pitch))

    def offset(self, position: float) -> float:
        """How far a stage reading sits from its lattice point, in um."""
        if self.pitch <= 0.0:
            return abs(position - self.origin)
        return abs(position - (self.origin + self.index(position) * self.pitch))

    def holds(self, position: float) -> bool:
        """Whether a stage reading lies on this lattice.

        Deliberately unbounded: a raster MassLynx split across functions
        continues past the extent of any one chunk, so the tail of a split
        raster lies on the lattice at an index beyond ``count``. What is
        being asked is whether the position shares the raster's phase, not
        whether it falls inside the part of it seen so far.
        """
        if not isfinite(position):
            return False
        if self.pitch <= 0.0:
            return position == self.origin
        return self.offset(position) <= _LATTICE_TOLERANCE * self.pitch


@dataclass(frozen=True)
class RegionPlacement:
    """Where one declared region sits on a grid that holds several.

    Attributes:
        number: The region's number in the method's table.
        x_lattice: The lattice fitted to this region's readings alone.
        y_lattice: Likewise on y.
        x_offset: Grid column of the region's first column.
        y_offset: Grid row of the region's first row.
        x_shift_um: The region's stage origin minus the stage position of
            the grid column it was placed on, in um, signed (+9.2 um when the
            origin lies 9.2 um beyond that column). Regions are fitted
            separately and placed on the grid's pitch, so a region whose
            origin is not a whole number of pitches from the grid origin
            lies off its column by this much.
        y_shift_um: Likewise on y, against the grid row.
        n_scans: Positioned scans of the fitted functions in this region.
    """

    number: int
    x_lattice: AxisLattice
    y_lattice: AxisLattice
    x_offset: int
    y_offset: int
    x_shift_um: float
    y_shift_um: float
    n_scans: int


@dataclass
class ImagingGrid:
    """Reconstructed imaging grid from Waters laser coordinates.

    Built by scanning all functions/scans and collecting unique laser
    positions. Provides O(1) coordinate lookup during spectrum iteration.
    """

    x_index_map: Dict[float, int]  # position (um) -> 0-based x index
    y_index_map: Dict[float, int]  # position (um) -> 0-based y index
    pixel_count_x: int
    pixel_count_y: int
    pixel_size_x: float  # raster pitch in um; 0.0 when the axis has one position
    pixel_size_y: float  # raster pitch in um; 0.0 when the axis has one position
    lateral_width: float  # max_x - min_x in um, so centre-to-centre, not edge-to-edge
    lateral_height: float  # max_y - min_y in um, so centre-to-centre, not edge-to-edge
    scan_map: Dict[Tuple[int, int], ScanInfoData] = field(repr=False)
    #: The fitted raster, when this grid was built by :func:`build_imaging_grid`.
    #: ``None`` on a hand-built grid, where every position is taken on trust.
    x_lattice: Optional[AxisLattice] = None
    y_lattice: Optional[AxisLattice] = None
    #: Why this grid is no raster, when every positioned scan reports the
    #: same stage position. Only a grid built with ``allow_single_position``
    #: carries one; reading its spectra is refused with this message.
    no_raster_reason: Optional[str] = None
    #: Where each axis's ``pixel_size`` came from: ``"stage positions"`` (the
    #: fitted pitch), ``"method step"`` (a single-position axis whose pitch
    #: could not be measured, taken from the method's declared step) or
    #: ``"none"`` (still 0.0). Derived from the sizes on a hand-built grid.
    pixel_size_source: Optional[Tuple[str, str]] = None
    #: The declared regions that hold scans, one placement each. Empty when
    #: the method declares no regions or they could not be used.
    regions: Tuple[RegionPlacement, ...] = ()
    #: The region a ``region=`` request kept; scans of every other region
    #: have no pixel. ``None`` when the whole run is converted.
    selected_region: Optional[int] = None
    #: The method's declared regions, kept only while they are in use, so
    #: scans of functions outside the fitted set can be assigned too.
    declared: Tuple[DeclaredRegion, ...] = ()
    #: The method's stage pattern, kept so a re-fit reuses it.
    pattern: Optional[StagePattern] = None
    #: Position (um) to ``(x, y, region)`` for every placed scan, on a grid
    #: that holds several regions. A position alone does not say which pixel
    #: it is there: two regions can share an x reading and sit at different
    #: columns, so ``x_index_map`` and ``y_index_map`` keep the first
    #: region's index only.
    position_map: Dict[Tuple[float, float], Tuple[int, int, int]] = field(
        default_factory=dict, repr=False
    )

    def __post_init__(self) -> None:
        """Derive the pixel size sources of a grid built without them."""
        if self.pixel_size_source is None:
            self.pixel_size_source = (
                _STAGE_POSITIONS if self.pixel_size_x > 0 else _NO_SOURCE,
                _STAGE_POSITIONS if self.pixel_size_y > 0 else _NO_SOURCE,
            )

    @property
    def dimensions(self) -> Tuple[int, int, int]:
        """Grid dimensions as (n_x, n_y, n_z). z is always 1 for 2D imaging."""
        return (self.pixel_count_x, self.pixel_count_y, 1)

    @cached_property
    def scans_by_function(self) -> Dict[int, List[ScanInfoData]]:
        """Positioned scan records grouped by function, in scan order.

        Built once and cached. Function selection asks five different
        questions of one function's scans and used to re-sort the whole
        scan map for each: O(F x N log N) over the file, measured at 87 s
        of selection against a 29 s grid build on a 25 M scan raster
        (issue #235).

        Cached, so ``scan_map`` must not be mutated once coordinates have
        been queried.
        """
        grouped: Dict[int, List[ScanInfoData]] = {}
        for (func, _scan), info in sorted(self.scan_map.items()):
            if info.has_position:
                grouped.setdefault(func, []).append(info)
        return grouped

    def region_of(self, scan_info: ScanInfoData) -> Optional[int]:
        """The declared region a positioned scan belongs to.

        ``None`` when the grid uses no regions, or when the scan lies in no
        declared box (or in two). With regions in use a scan in no box is
        not on the raster, whether the whole run or one region is converted.
        """
        if not self.declared or not scan_info.has_position:
            return None
        key = (
            _mm_to_um_key(scan_info.laser_x_pos),
            _mm_to_um_key(scan_info.laser_y_pos),
        )
        placed = self.position_map.get(key)
        if placed is not None:
            return placed[2]
        return _region_at(key[0], key[1], self.declared, self.pattern)

    def _placement(self, number: Optional[int]) -> Optional[RegionPlacement]:
        """The placement of region ``number``, if the grid holds it."""
        for placement in self.regions:
            if placement.number == number:
                return placement
        return None

    def _wanted_region(self, scan_info: ScanInfoData) -> Optional[int]:
        """The scan's region, unless a ``region=`` request leaves it out."""
        region = self.region_of(scan_info)
        if region is None:
            return None
        if self.selected_region is not None and region != self.selected_region:
            return None
        return region

    def _out_of_scope(self, scan_info: ScanInfoData) -> bool:
        """Whether a scan is left out on purpose, so a miss needs no warning."""
        if not self.declared:
            return False
        region = self._wanted_region(scan_info)
        return region is None or self._placement(region) is None

    def lies_on_raster(self, scan_info: ScanInfoData) -> bool:
        """Whether a scan's stage reading shares the fitted raster's phase.

        True when the grid carries no fitted lattice, so a hand-built grid
        behaves as it always did. With regions in use, the declared box is
        the source's statement: a scan outside every box is not on the
        raster, one in a region the grid does not hold yet (a function not
        fitted so far) is, and one in a placed region must lie on that
        region's lattice. A scan of a region left out by ``region=`` is not.
        """
        if self.x_lattice is None or self.y_lattice is None:
            return True
        x_um = _mm_to_um_key(scan_info.laser_x_pos)
        y_um = _mm_to_um_key(scan_info.laser_y_pos)
        if self.declared:
            region = self._wanted_region(scan_info)
            if region is None:
                return False
            placement = self._placement(region)
            if placement is None:
                return True
            return placement.x_lattice.holds(x_um) and placement.y_lattice.holds(y_um)
        return self.x_lattice.holds(x_um) and self.y_lattice.holds(y_um)

    def get_coordinates(
        self, scan_info: ScanInfoData
    ) -> Optional[Tuple[int, int, int]]:
        """Map laser position (mm from DLL) to 0-based pixel coordinates.

        The DLL returns positions in mm. We multiply by 1000 to get um,
        matching the keys stored in x_index_map/y_index_map (which were
        also built from mm * 1000). A scan of a region the grid does not
        hold has no pixel and is not warned about: leaving it out is what
        ``region=`` asks for.
        """
        if not scan_info.has_position:
            return None

        coords = self.locate(scan_info)
        if coords is not None or self._out_of_scope(scan_info):
            return coords

        x_um = _mm_to_um_key(scan_info.laser_x_pos)
        y_um = _mm_to_um_key(scan_info.laser_y_pos)
        logger.warning(
            f"No index found for laser position x={scan_info.laser_x_pos:.4f}mm "
            f"({x_um:.1f}um) -> {self.x_index_map.get(x_um)}, "
            f"y={scan_info.laser_y_pos:.4f}mm "
            f"({y_um:.1f}um) -> {self.y_index_map.get(y_um)}"
        )
        return None

    def locate(self, scan_info: ScanInfoData) -> Optional[Tuple[int, int, int]]:
        """The pixel of a positioned scan, silently ``None`` when it has none."""
        x_um = _mm_to_um_key(scan_info.laser_x_pos)
        y_um = _mm_to_um_key(scan_info.laser_y_pos)

        if self.declared:
            placement = self._placement(self._wanted_region(scan_info))
            if placement is None:
                return None
            if len(self.regions) > 1:
                placed = self.position_map.get((x_um, y_um))
                if placed is not None:
                    return (placed[0], placed[1], 0)
                if not (
                    placement.x_lattice.holds(x_um) and placement.y_lattice.holds(y_um)
                ):
                    return None
                return (
                    placement.x_offset + placement.x_lattice.index(x_um),
                    placement.y_offset + placement.y_lattice.index(y_um),
                    0,
                )

        x_idx = self.x_index_map.get(x_um)
        y_idx = self.y_index_map.get(y_um)

        if x_idx is not None and y_idx is not None:
            return (x_idx, y_idx, 0)

        # A miss is not automatically a failure. The grid is fitted to the
        # functions that define the raster, so a function selection has not
        # accepted yet contributes no key to the maps -- including the tail
        # of a raster MassLynx split across functions, whose pixels are the
        # whole reason it gets rescued. A reading that lies on the raster
        # has a pixel even when the extent fitted so far does not reach it;
        # regrid_for_functions() grows the grid once the tail is kept, so
        # an index past the current count never survives into iteration.
        if (
            self.x_lattice is not None
            and self.y_lattice is not None
            and self.lies_on_raster(scan_info)
        ):
            return (self.x_lattice.index(x_um), self.y_lattice.index(y_um), 0)
        return None


def _margin(step_um: Optional[float]) -> float:
    """How far a region's box is widened along an axis, in um."""
    if step_um is None:
        return _BOX_MARGIN_FALLBACK_UM
    return _BOX_MARGIN_STEPS * step_um + 1e-6


def _boxes_holding(
    x_um: float,
    y_um: float,
    declared: Iterable[DeclaredRegion],
    pattern: Optional[StagePattern],
) -> List[int]:
    """Numbers of the declared boxes a stage position lies in.

    Each box is widened by half a step on each side, since a region's
    length can be smaller than one step.
    """
    margin_x = _margin(pattern.x_step_um if pattern else None)
    margin_y = _margin(pattern.y_step_um if pattern else None)
    holding = []
    for region in declared:
        x_lo, x_hi = sorted((region.x_start_um, region.x_start_um + region.x_length_um))
        y_lo, y_hi = sorted((region.y_start_um, region.y_start_um + region.y_length_um))
        if (
            x_lo - margin_x <= x_um <= x_hi + margin_x
            and y_lo - margin_y <= y_um <= y_hi + margin_y
        ):
            holding.append(region.number)
    return holding


def _region_at(
    x_um: float,
    y_um: float,
    declared: Iterable[DeclaredRegion],
    pattern: Optional[StagePattern],
) -> Optional[int]:
    """The one declared region whose box holds a stage position, else ``None``."""
    holding = _boxes_holding(x_um, y_um, declared, pattern)
    return holding[0] if len(holding) == 1 else None


def _mm_to_um_key(mm_value: float) -> float:
    """Convert mm to um and round for use as dict key.

    Rounding to 2 decimal places (0.01 um precision) avoids floating-point
    comparison issues when using floats as dictionary keys. The mzmine Java
    code uses Float.compare() which does exact bitwise comparison, but
    Python float arithmetic from c_float -> Python float conversion can
    introduce tiny differences.
    """
    return round(mm_value * 1000.0, 2)


def _fit_axis(positions: Iterable[float], axis: str) -> AxisLattice:
    """Fit origin, pitch and count to one axis's distinct stage readings.

    The pitch is seeded from the *median* neighbour interval, which is
    exact when every column is present and survives a minority of gaps: one
    missing row among ten turns a single interval into a double one and
    leaves the median at the true pitch, where ``extent / (n_distinct - 1)``
    reports it 12.5 % high (issue #227).

    It is then refined to ``extent / (count - 1)`` so the two end positions
    land exactly on lattice points, which keeps a complete raster on the
    integers it was acquired at: the four registry acquisitions still fit
    30.0000, 50.0000, 100.0000 and 500.0000 um exactly.

    Raises:
        ValueError: If the readings do not lie on a regular raster.
    """
    unique = sorted(set(positions))
    if len(unique) == 1:
        return AxisLattice(origin=unique[0], pitch=0.0, count=1)

    intervals = [b - a for a, b in zip(unique, unique[1:])]
    seed = statistics.median(intervals)
    span = unique[-1] - unique[0]
    if seed <= 0.0 or span <= 0.0:
        raise ConversionRefused(
            f"Waters stage readings on the {axis} axis have no measurable "
            f"raster pitch ({len(unique)} distinct positions spanning "
            f"{span:.4f} um)."
        )

    count = int(round(span / seed)) + 1
    pitch = span / (count - 1)
    lattice = AxisLattice(origin=unique[0], pitch=pitch, count=count)

    worst = max(unique, key=lattice.offset)
    worst_offset = lattice.offset(worst)
    if worst_offset > _LATTICE_TOLERANCE * pitch:
        raise ConversionRefused(
            f"Waters stage readings on the {axis} axis do not lie on a "
            f"regular raster: with a fitted pitch of {pitch:.4f} um the "
            f"reading at {worst:.2f} um sits {worst_offset:.4f} um "
            f"({100.0 * worst_offset / pitch:.1f} % of a pitch) off the "
            f"nearest raster line. Thyra will not guess which pixel it "
            f"belongs to."
        )

    occupancy = len(unique) / count
    if occupancy < _MIN_LATTICE_OCCUPANCY:
        raise ConversionRefused(
            f"Waters stage readings on the {axis} axis fit a raster of "
            f"{count} lines at {pitch:.4f} um but only {len(unique)} of "
            f"them carry a reading ({100.0 * occupancy:.1f} %). That is "
            f"sub-pitch jitter being read as the pitch, not a raster with "
            f"gaps: the {len(unique)} distinct readings span "
            f"{span:.2f} um. Thyra will not build a grid from it."
        )

    if len(unique) < count:
        logger.warning(
            "The Waters %s axis is missing %d of its %d raster lines "
            "(pitch %.4f um). The grid keeps every line, so the gaps become "
            "empty pixels rather than shifting the rows past them.",
            axis,
            count - len(unique),
            count,
            pitch,
        )

    return lattice


@dataclass
class _PositionCensus:
    """What the stage readings of a set of functions looked like."""

    x_positions: List[float] = field(default_factory=list)
    y_positions: List[float] = field(default_factory=list)
    positioned: int = 0
    unpositioned: int = 0
    non_finite: int = 0
    partial_sentinel: int = 0


def _census(
    scan_map: Dict[Tuple[int, int], ScanInfoData],
    functions: Optional[Iterable[int]] = None,
) -> _PositionCensus:
    """Collect the stage readings of the given functions (all, when None)."""
    wanted = None if functions is None else set(functions)
    census = _PositionCensus()
    for (func, _scan), info in scan_map.items():
        if wanted is not None and func not in wanted:
            continue
        x_mm, y_mm = info.laser_x_pos, info.laser_y_pos
        if not (isfinite(x_mm) and isfinite(y_mm)):
            census.non_finite += 1
            continue
        if not info.has_position:
            census.unpositioned += 1
            continue
        if x_mm == NO_POSITION or y_mm == NO_POSITION:
            census.partial_sentinel += 1
        census.x_positions.append(_mm_to_um_key(x_mm))
        census.y_positions.append(_mm_to_um_key(y_mm))
        census.positioned += 1
    return census


@dataclass
class _Scope:
    """The stage readings one grid is fitted to.

    ``by_region`` and ``counts`` are filled only while the method's declared
    regions are in use, and then hold just the regions the grid is kept to.
    """

    x_positions: Set[float]
    y_positions: Set[float]
    positioned: int
    declared: Tuple[DeclaredRegion, ...] = ()
    selected: Optional[int] = None
    by_region: Dict[int, Set[Tuple[float, float]]] = field(default_factory=dict)
    counts: Dict[int, int] = field(default_factory=dict)


@dataclass
class _Fit:
    """What was fitted to a scope: lattices, placements and pixel sizes."""

    x_lattice: AxisLattice
    y_lattice: AxisLattice
    x_index_map: Dict[float, int]
    y_index_map: Dict[float, int]
    pixel_sizes: Tuple[float, float]
    sources: Tuple[str, str]
    placements: Tuple[RegionPlacement, ...] = ()
    position_map: Dict[Tuple[float, float], Tuple[int, int, int]] = field(
        default_factory=dict
    )


def _assign_regions(
    census: _PositionCensus,
    declared: Tuple[DeclaredRegion, ...],
    pattern: Optional[StagePattern],
) -> Tuple[Optional[List[int]], str]:
    """The declared region of each census position, or ``None`` to use none.

    The second value says why the regions are not used (empty when they are).

    A scan belongs to a region when its position lies inside the region's
    box. If any scan lies in no box or in two, the declared regions say
    nothing reliable about this run: that is warned about once and the
    regions are not used, so the grid is fitted as it was before regions
    were read.
    """
    cache: Dict[Tuple[float, float], List[int]] = {}
    assigned: List[int] = []
    first: Optional[Tuple[Tuple[float, float], List[int]]] = None
    offending = 0
    for position in zip(census.x_positions, census.y_positions):
        holding = cache.get(position)
        if holding is None:
            holding = _boxes_holding(position[0], position[1], declared, pattern)
            cache[position] = holding
        if len(holding) == 1:
            assigned.append(holding[0])
            continue
        offending += 1
        assigned.append(0)
        if first is None:
            first = (position, holding)
    if first is None:
        return assigned, ""
    (x_um, y_um), holding = first
    where = (
        "in the boxes of regions " + ", ".join(str(n) for n in holding)
        if holding
        else "in no region's box"
    )
    reason = (
        f"{offending} scan(s) lie in no region's box or in more than one "
        f"(the first, at x={x_um:.2f} um, y={y_um:.2f} um, lies {where})"
    )
    logger.warning(
        "The method declares %d region(s), but %s. The declared regions are "
        "not used; the run is gridded as one raster.",
        len(declared),
        reason,
    )
    return None, reason


def _check_region_request(
    region: Optional[int], declared: Tuple[DeclaredRegion, ...]
) -> None:
    """Refuse a ``region=`` request that names no region of this run.

    A run with no usable declared regions is one region, numbered 1, as
    :func:`thyra.list_regions` reports it.
    """
    if region is None:
        return
    if declared:
        numbers = [r.number for r in declared]
        if region not in numbers:
            raise ConversionRefused(
                f"region={region} is not a region of this Waters run. The "
                f"method declares region(s) {', '.join(str(n) for n in numbers)}."
            )
    elif region != 1:
        raise ConversionRefused(
            f"region={region} is not a region of this Waters run. It has no "
            f"declared regions, so the whole run is region 1."
        )


def _resolve_scope(
    census: _PositionCensus,
    pattern: Optional[StagePattern],
    region: Optional[int],
) -> _Scope:
    """Decide which stage readings the grid is fitted to.

    Raises:
        ConversionRefused: If ``region`` names no region of the run, or a
            declared region none of the fitted scans lie in.
    """
    declared = pattern.regions if pattern is not None else ()
    assigned, reason = (
        _assign_regions(census, declared, pattern) if declared else (None, "")
    )
    if assigned is None:
        if declared and region is not None:
            raise ConversionRefused(
                f"region={region} cannot be used: the method declares "
                f"{len(declared)} region(s), but they could not be used "
                f"({reason}). This run can only be converted whole."
            )
        _check_region_request(region, ())
        return _Scope(
            set(census.x_positions), set(census.y_positions), census.positioned
        )
    _check_region_request(region, declared)

    by_region: Dict[int, Set[Tuple[float, float]]] = {}
    counts: Dict[int, int] = {}
    for position, number in zip(zip(census.x_positions, census.y_positions), assigned):
        by_region.setdefault(number, set()).add(position)
        counts[number] = counts.get(number, 0) + 1
    if region is not None:
        if region not in by_region:
            raise ConversionRefused(
                f"Region {region} of this Waters run holds no scan of the "
                f"functions being converted."
            )
        by_region = {region: by_region[region]}
        counts = {region: counts[region]}
    pairs = [p for members in by_region.values() for p in members]
    return _Scope(
        {p[0] for p in pairs},
        {p[1] for p in pairs},
        sum(counts.values()),
        declared=declared,
        selected=region,
        by_region=by_region,
        counts=counts,
    )


def _single_pixel_size(
    lattice: AxisLattice, step_um: Optional[float]
) -> Tuple[float, str]:
    """One axis's pixel size and where it came from."""
    if lattice.pitch > 0.0:
        return lattice.pitch, _STAGE_POSITIONS
    if step_um is not None:
        return step_um, _METHOD_STEP
    return 0.0, _NO_SOURCE


def _shared_pitch(
    axis: str, lattices: Dict[int, AxisLattice], step_um: Optional[float]
) -> Tuple[float, str]:
    """The one pitch the regions of a grid are placed on, along an axis.

    The regions that measured a pitch here (more than one position) must
    agree within :data:`_REGION_PITCH_AGREEMENT`; the median of their pitches
    is the grid's. With none, the method's declared step stands in. With
    neither, the regions can be placed only if they sit at one position.

    Raises:
        ConversionRefused: If the measured pitches disagree, or the regions
            cannot be placed on the axis without a step.
    """
    measured = {n: lat.pitch for n, lat in lattices.items() if lat.count > 1}
    if measured:
        median = statistics.median(measured.values())
        if any(
            abs(p - median) > _REGION_PITCH_AGREEMENT * median
            for p in measured.values()
        ):
            listing = ", ".join(f"region {n}: {p:.4f} um" for n, p in measured.items())
            raise ConversionRefused(
                f"The Waters regions were not rastered on one pitch along "
                f"{axis} ({listing}), so they cannot be placed on one grid. "
                f"Convert one region at a time (region=)."
            )
        return float(median), _STAGE_POSITIONS
    if step_um is not None:
        return step_um, _METHOD_STEP
    if len({round(lat.origin, 2) for lat in lattices.values()}) == 1:
        return 0.0, _NO_SOURCE
    raise ConversionRefused(
        f"The Waters regions each hold one position along {axis} and the "
        f"method declares no step there, so they cannot be placed on one "
        f"grid. Convert one region at a time (region=)."
    )


def _place(origin: float, base: float, pitch: float) -> Tuple[int, float]:
    """Grid index of a region's first line and the signed shift that costs, in um."""
    if pitch <= 0.0:
        return 0, 0.0
    offset = int(round((origin - base) / pitch))
    return offset, (origin - base) - offset * pitch


def _warn_on_shifts(
    placements: List[RegionPlacement], pitches: Tuple[float, float]
) -> None:
    """Say which regions lie off the grid line they were placed on."""
    moved = []
    for placement in placements:
        for axis, shift, pitch in (
            ("x", placement.x_shift_um, pitches[0]),
            ("y", placement.y_shift_um, pitches[1]),
        ):
            if abs(shift) > _REPORTABLE_SHIFT_UM and pitch > 0.0:
                line = "pixel column" if axis == "x" else "pixel row"
                moved.append(
                    f"region {placement.number}'s stage origin lies {shift:+.2f} "
                    f"um in {axis} from the {line} it was placed on "
                    f"({100.0 * abs(shift) / pitch:.1f} % of a pitch)"
                )
    if moved:
        logger.warning(
            "The Waters regions were fitted separately and placed on one "
            "grid. Their origins are not whole pitches apart, so each is "
            "placed on the nearest pixel column or row: %s.",
            "; ".join(moved),
        )


def _fit_regions(scope: _Scope, pattern: Optional[StagePattern]) -> _Fit:
    """Fit each region's own lattice and place them all on one grid.

    Every region is fitted alone, so the off-raster refusal only ever
    compares readings that were taken in one raster. They share one pitch
    per axis, and the grid's origin is the smallest region origin.
    """
    numbers = sorted(scope.by_region)
    x_fits = {
        n: _fit_axis({p[0] for p in scope.by_region[n]}, f"x (region {n})")
        for n in numbers
    }
    y_fits = {
        n: _fit_axis({p[1] for p in scope.by_region[n]}, f"y (region {n})")
        for n in numbers
    }
    x_pitch, x_source = _shared_pitch(
        "x", x_fits, pattern.x_step_um if pattern else None
    )
    y_pitch, y_source = _shared_pitch(
        "y", y_fits, pattern.y_step_um if pattern else None
    )
    x_origin = min(lat.origin for lat in x_fits.values())
    y_origin = min(lat.origin for lat in y_fits.values())

    placements: List[RegionPlacement] = []
    position_map: Dict[Tuple[float, float], Tuple[int, int, int]] = {}
    for n in numbers:
        x_off, x_shift = _place(x_fits[n].origin, x_origin, x_pitch)
        y_off, y_shift = _place(y_fits[n].origin, y_origin, y_pitch)
        placements.append(
            RegionPlacement(
                n, x_fits[n], y_fits[n], x_off, y_off, x_shift, y_shift, scope.counts[n]
            )
        )
        for x_um, y_um in scope.by_region[n]:
            position_map.setdefault(
                (x_um, y_um),
                (x_off + x_fits[n].index(x_um), y_off + y_fits[n].index(y_um), n),
            )
    _warn_on_shifts(placements, (x_pitch, y_pitch))

    x_lattice = AxisLattice(
        x_origin, x_pitch, max(p.x_offset + p.x_lattice.count for p in placements)
    )
    y_lattice = AxisLattice(
        y_origin, y_pitch, max(p.y_offset + p.y_lattice.count for p in placements)
    )
    x_index_map: Dict[float, int] = {}
    y_index_map: Dict[float, int] = {}
    for (x_um, y_um), (x_idx, y_idx, _n) in position_map.items():
        x_index_map.setdefault(x_um, x_idx)
        y_index_map.setdefault(y_um, y_idx)
    return _Fit(
        x_lattice,
        y_lattice,
        x_index_map,
        y_index_map,
        (x_pitch, y_pitch),
        (x_source, y_source),
        placements=tuple(placements),
        position_map=position_map,
    )


def _fit_scope(scope: _Scope, pattern: Optional[StagePattern]) -> _Fit:
    """Fit the grid to a scope: one lattice, or one per region placed together."""
    if len(scope.by_region) > 1:
        return _fit_regions(scope, pattern)
    x_lattice = _fit_axis(scope.x_positions, "x")
    y_lattice = _fit_axis(scope.y_positions, "y")
    x_size, x_source = _single_pixel_size(
        x_lattice, pattern.x_step_um if pattern else None
    )
    y_size, y_source = _single_pixel_size(
        y_lattice, pattern.y_step_um if pattern else None
    )
    placements: Tuple[RegionPlacement, ...] = ()
    if scope.by_region:
        (number,) = scope.by_region
        placements = (
            RegionPlacement(
                number, x_lattice, y_lattice, 0, 0, 0.0, 0.0, scope.counts[number]
            ),
        )
    return _Fit(
        x_lattice,
        y_lattice,
        {pos: x_lattice.index(pos) for pos in scope.x_positions},
        {pos: y_lattice.index(pos) for pos in scope.y_positions},
        (x_size, y_size),
        (x_source, y_source),
        placements=placements,
    )


def _warn_on_census(census: _PositionCensus) -> None:
    """Say what the census had to leave out or take at face value."""
    if census.non_finite:
        logger.warning(
            "%d Waters scan(s) report a non-finite stage position (NaN or "
            "infinity) and are treated as unpositioned. A non-finite "
            "reading used to be admitted as its own column, making the "
            "reported pixel size NaN or infinite.",
            census.non_finite,
        )
    if census.partial_sentinel:
        logger.warning(
            "%d Waters scan(s) report the no-position sentinel (%.1f mm) on "
            "one axis and a real reading on the other. They are kept at "
            "face value, so a genuine stage position at %.1f mm on both "
            "axes would be dropped instead; check the raster origin if the "
            "image looks shifted.",
            census.partial_sentinel,
            NO_POSITION,
            NO_POSITION,
        )


def _warn_on_sentinel_pixel(
    census: _PositionCensus, x_lattice: AxisLattice, y_lattice: AxisLattice
) -> None:
    """Say so when the raster passes through the no-position sentinel.

    The sentinel is a coordinate a stage can really visit, and a raster
    that passes through it loses that pixel with no message. Only worth
    saying when the raster actually reaches -1.0 mm on both axes, which is
    why this is checked against the fitted lattice rather than warned about
    on every file: the registry's DESI runs do use negative stage
    millimetres (x from -11.4 mm, y from -6.0 mm), so this is not
    hypothetical for them (issue #231).
    """
    sentinel_um = _mm_to_um_key(NO_POSITION)
    if (
        census.unpositioned
        and x_lattice.holds(sentinel_um)
        and y_lattice.holds(sentinel_um)
        and 0 <= x_lattice.index(sentinel_um) < x_lattice.count
        and 0 <= y_lattice.index(sentinel_um) < y_lattice.count
    ):
        logger.warning(
            "This raster passes through (%.1f, %.1f) mm, which is also the "
            "no-position sentinel, so the %d scan(s) reporting it were "
            "dropped as unpositioned. Pixel (%d, %d) may be empty because "
            "of that rather than because nothing was acquired there.",
            NO_POSITION,
            NO_POSITION,
            census.unpositioned,
            x_lattice.index(sentinel_um),
            y_lattice.index(sentinel_um),
        )


def _grid_from_scan_map(
    scan_map: Dict[Tuple[int, int], ScanInfoData],
    functions: Optional[Iterable[int]] = None,
    allow_single_position: bool = False,
    pattern: Optional[StagePattern] = None,
    region: Optional[int] = None,
) -> ImagingGrid:
    """Fit a grid to the stage readings already collected in ``scan_map``.

    Split out of :func:`build_imaging_grid` so the reader can re-fit the
    raster to the functions it actually converts without touching the DLL
    again. Deriving the pitch from every function -- including the ones
    function selection then excludes -- let one off-raster reference spot
    turn a 10x5 image at 30 um into an 11x6 one at 1097 x 440 um
    (issue #232).

    ``allow_single_position`` returns a one-pixel grid, carrying the
    refusal as ``no_raster_reason``, where a conversion refuses the run:
    a metadata document describes the acquisition whether or not it is
    an image (issue #494).

    ``pattern`` is the stage pattern the run's method declares. Its steps
    stand in for a pitch that cannot be measured on a single-position
    axis. Its regions, when every scan lies in exactly one of them, are
    fitted one by one and placed on a shared grid (see
    :class:`RegionPlacement`); ``region`` keeps one of them.
    """
    census = _census(scan_map, functions)
    _warn_on_census(census)

    if not census.x_positions or not census.y_positions:
        raise ConversionRefused("No valid laser positions found in Waters imaging data")

    scope = _resolve_scope(census, pattern, region)
    fit = _fit_scope(scope, pattern)
    x_lattice, y_lattice = fit.x_lattice, fit.y_lattice

    if len(fit.placements) < 2:
        _warn_on_sentinel_pixel(census, x_lattice, y_lattice)

    # A stage that never moved is not a raster. The check above only catches
    # the total absence of positions; a *constant* position passes it, and
    # everything downstream then agrees the run is a legitimate 1x1 image
    # with a 0.0 um pitch (lateral extent is 0, so both pixel sizes below
    # come out 0.0) and converts the whole acquisition onto one pixel
    # without a word. Measured on the Xevo DESI share, 76 of 105 .raw dirs
    # land here -- every one of them a calibration, a tuning run, a lysis
    # test or a single-spot DDA acquisition, and not one of them an image.
    # Real DESI images on that same share do carry stage coordinates in
    # these fields and grid normally (401x401, 247x140, ...), so this
    # refuses the non-images without costing anything real. See issue #213.
    no_raster_reason: Optional[str] = None
    pixel_sizes, sources = fit.pixel_sizes, fit.sources
    if x_lattice.count == 1 and y_lattice.count == 1:
        no_raster_reason = (
            f"Every positioned scan reports the same stage position "
            f"(x={x_lattice.origin:.1f} um, y={y_lattice.origin:.1f} um), so "
            f"this acquisition has a single pixel: {scope.positioned} "
            f"positioned scans in {len(scan_map)} total. A one-pixel raster "
            f"with a 0.0 um pitch is a spot, calibration or tuning "
            f"acquisition, not an image -- Thyra has no raster to build "
            f"from it."
        )
        if not allow_single_position:
            raise ConversionRefused(no_raster_reason)
        pixel_sizes, sources = (0.0, 0.0), (_NO_SOURCE, _NO_SOURCE)

    # Laser positions are pixel *centres*, so the span from the first to the
    # last is N - 1 pitches, not N. This used to divide by the count, which
    # mixed a centre-to-centre extent with an edge-to-edge one and reported
    # a pitch low by exactly 1/N -- independently per axis, so a square
    # raster came out anisotropic: measured, a 30 um acquisition gave
    # 29.66 x 28.93 and a 100 um one gave 99.40 x 97.83. See issue #217.
    # (mzmine, which this file translates, computes the extent the other way
    # round -- lateralWidth = count * pixelWidth -- and has the
    # divide-by-count form commented out in ImagingParameters.java.)
    lateral_width = x_lattice.pitch * (x_lattice.count - 1)
    lateral_height = y_lattice.pitch * (y_lattice.count - 1)

    logger.info(
        f"Built imaging grid: {x_lattice.count}x{y_lattice.count} pixels, "
        f"pixel size: {pixel_sizes[0]:.1f}x{pixel_sizes[1]:.1f} um, "
        f"lateral: {lateral_width:.1f}x{lateral_height:.1f} um, "
        f"{scope.positioned} positioned scans in {len(scan_map)} total scans"
    )

    grid = ImagingGrid(
        x_index_map=fit.x_index_map,
        y_index_map=fit.y_index_map,
        pixel_count_x=x_lattice.count,
        pixel_count_y=y_lattice.count,
        pixel_size_x=pixel_sizes[0],
        pixel_size_y=pixel_sizes[1],
        lateral_width=lateral_width,
        lateral_height=lateral_height,
        scan_map=scan_map,
        x_lattice=x_lattice,
        y_lattice=y_lattice,
        no_raster_reason=no_raster_reason,
        pixel_size_source=sources,
        regions=fit.placements,
        selected_region=scope.selected,
        declared=scope.declared,
        pattern=pattern,
        position_map=fit.position_map,
    )
    if no_raster_reason is None:
        # On a single position every scan shares the one pixel; saying so
        # adds nothing to the reason the grid carries.
        _warn_on_shared_pixels(grid, functions)
    return grid


def _warn_on_shared_pixels(
    grid: ImagingGrid, functions: Optional[Iterable[int]]
) -> None:
    """Say so when scans of different functions land on the same pixel.

    Those are summed into one pixel by the converter: two mass ranges or
    channels of one acquisition. Scans *one* function took twice on a
    pixel -- the stage reporting the same position again, as on the
    registry's 100 um MALDI set, whose function 2 holds 1275 positioned
    scans on 1274 distinct pixels (issue #233) -- are the reader's to
    combine, and it reports what it did with them (D85). The pixel of each
    scan is the one the grid gives it, so on a grid of several regions two
    regions can never be mistaken for one pixel.
    """
    wanted = None if functions is None else set(functions)
    functions_on: Dict[Tuple[int, int], Set[int]] = {}
    for (func, _scan), info in grid.scan_map.items():
        if wanted is not None and func not in wanted:
            continue
        if not info.has_position:
            continue
        coords = grid.locate(info)
        if coords is None:
            continue
        functions_on.setdefault((coords[0], coords[1]), set()).add(func)

    shared = sorted(p for p, on in functions_on.items() if len(on) > 1)
    if not shared:
        return

    sample = ", ".join(f"({x}, {y})" for x, y in shared[:5])
    if len(shared) > 5:
        sample += ", ..."
    logger.warning(
        "%d pixel(s) hold scans of more than one Waters function, so they "
        "carry the SUM of those functions' scans: %s.",
        len(shared),
        sample,
    )


def build_imaging_grid(
    ml: MassLynxLib,
    handle,
    function_types: Dict[int, FunctionType],
    include_functions: Optional[Iterable[int]] = None,
    allow_single_position: bool = False,
    pattern: Optional[StagePattern] = None,
    region: Optional[int] = None,
) -> ImagingGrid:
    """Build imaging grid by scanning all functions/scans for laser coordinates.

    Translates ImagingMetadata.java constructor logic (lines 56-149).

    Every function's scan records are cached in ``scan_map`` whatever
    ``include_functions`` says, so the raster can be re-fitted later
    (:func:`regrid_for_functions`) without touching the DLL again.

    Args:
        ml: MassLynxLib instance with open handle.
        handle: The opaque file handle from ml.open_file().
        function_types: Pre-classified function types for all functions.
        include_functions: Fit the raster to these functions' stage
            readings only. ``None`` fits it to every function, which is
            what a caller that has not classified them yet wants.
        allow_single_position: Return a one-pixel grid that carries the
            refusal as ``no_raster_reason`` instead of raising it, for a
            caller that reads metadata only.
        pattern: The stage pattern the run's method declares, if any: its
            steps give a single-position axis its pixel size, and its
            regions are placed one by one on a shared grid.
        region: Keep one declared region. A run that declares none is
            region 1 as a whole.

    Returns:
        ImagingGrid with coordinate maps and pixel metadata.

    Raises:
        ValueError: If no scan carries a laser position, if the positions do
            not lie on a regular raster, or if every positioned scan carries
            the *same* one -- a single-pixel "raster" is a spot acquisition,
            not an image -- unless ``allow_single_position`` is set. Also
            if ``region`` names no region of the run.
    """
    scan_map: Dict[Tuple[int, int], ScanInfoData] = {}

    n_functions = ml.get_number_of_functions(handle)
    for func in range(n_functions):
        n_scans = ml.get_number_of_scans_in_function(handle, func)
        for scan in range(n_scans):
            scan_map[(func, scan)] = ml.get_scan_info(handle, func, scan)

    return _grid_from_scan_map(
        scan_map, include_functions, allow_single_position, pattern, region
    )


def regrid_for_functions(
    grid: ImagingGrid,
    functions: Iterable[int],
    allow_single_position: bool = False,
    pattern: Optional[StagePattern] = None,
    region: Optional[int] = None,
) -> ImagingGrid:
    """Re-fit the raster to the stage readings of ``functions`` alone.

    Uses the scan records already cached on ``grid``, so this costs no
    further DLL calls. Returns ``grid`` unchanged when the restriction
    leaves the geometry identical, which is every acquisition whose
    excluded functions sat on the raster with the rest, and when ``grid``
    carries no fitted lattice -- a hand-built grid is taken on trust, the
    same way :meth:`ImagingGrid.lies_on_raster` takes it.
    ``allow_single_position`` is passed on as :func:`build_imaging_grid`
    takes it. ``pattern`` and ``region`` default to the ones ``grid`` was
    built with.
    """
    if grid.x_lattice is None or grid.y_lattice is None:
        return grid
    refitted = _grid_from_scan_map(
        grid.scan_map,
        functions,
        allow_single_position,
        pattern if pattern is not None else grid.pattern,
        region if region is not None else grid.selected_region,
    )
    same = (
        refitted.pixel_count_x == grid.pixel_count_x
        and refitted.pixel_count_y == grid.pixel_count_y
        and refitted.pixel_size_x == grid.pixel_size_x
        and refitted.pixel_size_y == grid.pixel_size_y
        and refitted.x_index_map == grid.x_index_map
        and refitted.y_index_map == grid.y_index_map
        and refitted.regions == grid.regions
        and refitted.position_map == grid.position_map
    )
    if same:
        return grid
    logger.warning(
        "The imaging grid was re-fitted to the converted functions alone: "
        "%dx%d at %.4f x %.4f um, was %dx%d at %.4f x %.4f um. The excluded "
        "functions' stage readings were part of the raster Thyra measured.",
        refitted.pixel_count_x,
        refitted.pixel_count_y,
        refitted.pixel_size_x,
        refitted.pixel_size_y,
        grid.pixel_count_x,
        grid.pixel_count_y,
        grid.pixel_size_x,
        grid.pixel_size_y,
    )
    return refitted
