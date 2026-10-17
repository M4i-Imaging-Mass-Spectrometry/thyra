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
