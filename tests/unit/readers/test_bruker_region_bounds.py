"""Per-region bounding boxes on ``BrukerReader.get_region_info``.

A caller splitting a multi-region ``.d`` into one sample per tissue
section has to show the user which region is which.  Neither field that
was already there can do it: FlexImaging names the areas "01", "02",
"03", and two serial sections of the same tissue have nearly equal
spectrum counts.  Where a region sits on the slide is what tells them
apart, so ``bounds`` is the field that makes the rest usable.

The SQL is exercised against a real in-memory ``MaldiFrameInfo`` rather
than mocked out: the aggregate and the offset subtraction are the entire
substance of the method, and a mock of either would only assert that the
test author and the implementation agree.
"""

from __future__ import annotations

import sqlite3
from typing import Any, Dict, List, Optional, Tuple

from thyra.readers.bruker.timstof.timstof_reader import BrukerReader


class _BoundsHarness:
    """A stand-in exposing only what the two methods under test touch.

    ``get_region_info`` reads ``_region_info`` and ``_mis_metadata`` and
    calls ``_region_bounds``; ``_region_bounds`` reads ``conn`` and
    ``_get_coordinate_offsets``.  Nothing else on the real reader is
    reachable from here, which is why the heavy SDK initialisation can
    be skipped.
    """

    _mis_metadata: Dict[str, Any]
    _region_info: List[Tuple[int, int]]
    conn: sqlite3.Connection
    _offsets: Optional[Tuple[int, int, int]]

    # Bound off the real class so a change in behaviour reaches the test.
    _region_bounds = BrukerReader._region_bounds
    get_region_info = BrukerReader.get_region_info

    def _get_coordinate_offsets(self) -> Optional[Tuple[int, int, int]]:
        return self._offsets


#: One row per acquired frame: ``(RegionNumber, XIndexPos, YIndexPos)``.
#: Region 0 is a 3x2 block at x 10..12, region 1 a 2x2 block at x 20..21,
#: so the two are separable by position alone -- which is the point.
_FRAMES = [
    (0, 10, 5),
    (0, 11, 5),
    (0, 12, 5),
    (0, 10, 6),
    (0, 11, 6),
    (0, 12, 6),
    (1, 20, 7),
    (1, 21, 7),
    (1, 20, 8),
    (1, 21, 8),
]


def _harness(
    *,
    areas: Optional[List[str]] = None,
    region_info: Optional[List[Tuple[int, int]]] = None,
    offsets: Optional[Tuple[int, int, int]] = (10, 5, 0),
    frames: Optional[List[Tuple[int, int, int]]] = None,
    with_table: bool = True,
) -> _BoundsHarness:
    h = _BoundsHarness()
    h._mis_metadata = {"areas": [{"name": n} for n in (areas or ["01", "02"])]}
    h._region_info = region_info if region_info is not None else [(0, 6), (1, 4)]
    h._offsets = offsets
    h.conn = sqlite3.connect(":memory:")
    if with_table:
        h.conn.execute(
            "CREATE TABLE MaldiFrameInfo "
            "(RegionNumber INTEGER, XIndexPos INTEGER, YIndexPos INTEGER)"
        )
        h.conn.executemany(
            "INSERT INTO MaldiFrameInfo VALUES (?, ?, ?)",
            frames if frames is not None else _FRAMES,
        )
    return h


class TestRegionBounds:
    """``_region_bounds``: the aggregate and the normalization."""

    def test_a_box_is_the_min_and_max_of_the_region(self) -> None:
        # Offsets of zero so this test sees the raw extents only.
        h = _harness(offsets=(0, 0, 0))
        assert h._region_bounds() == {0: (10, 5, 12, 6), 1: (20, 7, 21, 8)}

    def test_boxes_are_normalized_by_the_coordinate_offsets(self) -> None:
        """Boxes must land in the frame of ``grid_dims``, not raw index space.

        ``get_region_map`` subtracts the same offsets, and the offsets
        themselves are not public -- so a raw box could not be drawn
        against the grid the caller already has.
        """
        h = _harness(offsets=(10, 5, 0))
        assert h._region_bounds() == {0: (0, 0, 2, 1), 1: (10, 2, 11, 3)}

    def test_absent_offsets_are_treated_as_zero(self) -> None:
        """A reader with no offsets reports raw positions, not a crash."""
        h = _harness(offsets=None)
        assert h._region_bounds() == {0: (10, 5, 12, 6), 1: (20, 7, 21, 8)}

    def test_a_missing_table_is_empty_not_an_exception(self) -> None:
        """No MaldiFrameInfo means no bounds; the rest of the summary stands.

        Bounds are the newest of the three fields and the only one that
        needs this table, so failing to read it must not take
        ``region_number`` and ``n_spectra`` down with it.
        """
        h = _harness(with_table=False)
        assert h._region_bounds() == {}


class TestRegionInfoCarriesBounds:
    """``get_region_info``: what the preview and the converter see."""

    def test_each_region_carries_its_box(self) -> None:
        h = _harness()
        assert h.get_region_info() == [
            {
                "region_number": 0,
                "n_spectra": 6,
                "bounds": (0, 0, 2, 1),
                "name": "01",
            },
            {
                "region_number": 1,
                "n_spectra": 4,
                "bounds": (10, 2, 11, 3),
                "name": "02",
            },
        ]

    def test_a_region_absent_from_the_table_simply_has_no_box(self) -> None:
        """``bounds`` is optional per the base contract, and stays optional.

        A region counted in ``_region_info`` but missing from the bounds
        query is a disagreement between two reads of the same database.
        Omitting the key says so; inventing a zero box would not.
        """
        h = _harness(
            areas=["01", "02", "03"],
            region_info=[(0, 6), (1, 4), (2, 9)],
        )
        summary = h.get_region_info()
        assert summary is not None
        assert "bounds" not in summary[2]
        assert summary[2]["region_number"] == 2
        assert summary[2]["n_spectra"] == 9

    def test_a_single_region_set_is_still_None(self) -> None:
        """Unchanged: one region is not a region list.

        Pinned here because ``_region_bounds`` runs a query and the
        cheapest way to break this would be to run it before the guard.
        """
        h = _harness(areas=["01"], region_info=[(0, 6)])
        assert h.get_region_info() is None

    def test_the_order_is_still_the_readers_own(self) -> None:
        """Frame-count order, NOT region number -- normalized by the caller.

        ``uns["regions"]`` in every already-converted store is in this
        order, so changing it here would silently change what the
        converter writes.  ``preview_msi`` sorts instead.
        """
        h = _harness(region_info=[(1, 4), (0, 6)])
        summary = h.get_region_info()
        assert summary is not None
        assert [r["region_number"] for r in summary] == [1, 0]
