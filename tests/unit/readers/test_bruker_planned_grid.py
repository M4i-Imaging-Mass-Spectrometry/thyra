"""The raster a timsTOF run planned, for the FlexImaging lattice fit (D59).

GlobalMetadata's ImagingArea bounds are the planned raster. On two timsTOF
runs that stopped early they spanned every row of the Area while the
frames covered only the first few, and width times height was the spot
count the FlexImaging info file planned (51 x 48 = 2448). The lattice fit
uses them to place such a run (#497).

The SQL runs against a real in-memory GlobalMetadata table, not a mock.
"""

from __future__ import annotations

import sqlite3
from typing import List, Optional, Tuple

import pytest

from thyra.readers.bruker.timstof.timstof_reader import BrukerReader


def _cursor(rows: List[Tuple[str, str]]) -> sqlite3.Cursor:
    conn = sqlite3.connect(":memory:")
    conn.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
    conn.executemany("INSERT INTO GlobalMetadata VALUES (?, ?)", rows)
    return conn.cursor()


BOUNDS = [
    ("ImagingAreaMinXIndexPos", "8424"),
    ("ImagingAreaMaxXIndexPos", "8474"),
    ("ImagingAreaMinYIndexPos", "4159"),
    ("ImagingAreaMaxYIndexPos", "4206"),
]


def test_the_imaging_area_is_the_planned_grid() -> None:
    grid = BrukerReader._planned_grid(_cursor(BOUNDS))

    assert grid == (8424, 4159, 51, 48)


@pytest.mark.parametrize(
    "rows",
    [
        BOUNDS[:3],
        [*BOUNDS[:3], ("ImagingAreaMaxYIndexPos", "n/a")],
        [*BOUNDS[:3], ("ImagingAreaMaxYIndexPos", "4000")],
        [],
    ],
)
def test_no_planned_grid_without_four_ordered_bounds(
    rows: List[Tuple[str, str]],
) -> None:
    grid: Optional[Tuple[int, int, int, int]] = BrukerReader._planned_grid(
        _cursor(rows)
    )

    assert grid is None


# --- Values read exactly, and empty ones treated as absent (#520) ------------

BIG = 2**53 + 1


def test_bounds_past_two_to_the_53_are_read_exactly() -> None:
    """F133: ``int(float(...))`` read 2**53 + 1 as 2**53, a step off."""
    rows = [
        ("ImagingAreaMinXIndexPos", str(BIG)),
        ("ImagingAreaMaxXIndexPos", str(BIG + 2)),
        ("ImagingAreaMinYIndexPos", "10"),
        ("ImagingAreaMaxYIndexPos", "11"),
    ]
    assert BrukerReader._planned_grid(_cursor(rows)) == (BIG, 10, 3, 2)


@pytest.mark.parametrize("value", ["", "  ", None])
def test_an_empty_bound_is_no_bound(value: Optional[str]) -> None:
    """F144: an empty or NULL value ended the conversion in a traceback."""
    rows = [*BOUNDS[:3], ("ImagingAreaMaxYIndexPos", value)]
    assert BrukerReader._planned_grid(_cursor(rows)) is None  # type: ignore[arg-type]


class _HeaderHarness:
    """Only what ``_build_header_alignment`` reads off a reader."""

    _build_header_alignment = BrukerReader._build_header_alignment
    _planned_grid = staticmethod(BrukerReader._planned_grid)

    def __init__(self, conn: sqlite3.Connection) -> None:
        self._conn = conn
        self._mis_metadata = {"areas": [{"name": "a", "p1": [0, 0], "p2": [1, 1]}]}
        self._selected_region = None

    def _database(self) -> sqlite3.Connection:
        return self._conn


def _frames_db(rows: List[Tuple[str, Optional[str]]]) -> sqlite3.Connection:
    conn = sqlite3.connect(":memory:")
    conn.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
    conn.executemany("INSERT INTO GlobalMetadata VALUES (?, ?)", rows)
    conn.execute(
        "CREATE TABLE MaldiFrameInfo (Frame INTEGER, RegionNumber INTEGER, "
        "XIndexPos INTEGER, YIndexPos INTEGER)"
    )
    conn.executemany(
        "INSERT INTO MaldiFrameInfo VALUES (?, 0, ?, ?)",
        [(i, 100 + i % 3, 10 + i // 3) for i in range(6)],
    )
    return conn


@pytest.mark.parametrize(
    "rows, expected",
    [
        (BOUNDS, (8424, 4159)),
        ([], (100, 10)),
        ([BOUNDS[0]], (8424, 10)),
        ([("ImagingAreaMinXIndexPos", ""), BOUNDS[2]], (100, 4159)),
        (
            [("ImagingAreaMinXIndexPos", None), ("ImagingAreaMinYIndexPos", "")],
            (100, 10),
        ),
    ],
    ids=["stated", "absent", "only x stated", "x empty", "both unreadable"],
)
def test_the_alignment_origin_falls_back_to_the_frames_key_by_key(
    rows: List[Tuple[str, Optional[str]]], expected: Tuple[int, int]
) -> None:
    """F172: a missing key put the origin at (0, 0), or raised IndexError."""
    header = _HeaderHarness(_frames_db(rows))._build_header_alignment()
    assert (header["first_raster_x"], header["first_raster_y"]) == expected
