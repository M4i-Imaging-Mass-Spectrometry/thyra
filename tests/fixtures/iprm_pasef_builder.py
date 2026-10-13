"""Turn a copy of ``synthetic_tims.d`` into a type-9 (DIA-layout) acquisition.

Type-9 (DIA) fragment frames -- the layout MALDI iprm-PASEF is reported to
use; no such MALDI run has been measured -- keep their isolation windows in
the DIA tables: ``DiaFrameMsMsInfo`` maps each
fragment frame to a window group and ``DiaFrameMsMsWindows`` lists the
group's windows, one row per ``[ScanNumBegin, ScanNumEnd)`` slice of the
mobility ramp. A survey frame (type 0) has no row in ``DiaFrameMsMsInfo``.
That is the layout a measured LC diaPASEF run uses: 2,755 survey frames,
22,039 fragment frames, eight window groups of three windows taking turns
frame by frame, adjacent windows of a group sharing their boundary scan.

Three shapes, all from the committed fixture's six frames and six raster
spots, changing only the database:

``"mixed"`` (the default)
    Frames 1-3 (row ``y = 10``) are survey frames, frames 4-6 (row
    ``y = 11``) fragment frames. Frame 7 is added: a fragment frame on
    frame 1's spot ``(100, 10)``, so one spot holds a survey frame and a
    fragment frame. It reuses frame 4's block in ``analysis.tdf_bin``
    (same ``TimsId``, ``NumScans``, ``NumPeaks``, ``SummedIntensities``),
    which the vendor library reads like any other.
``"one_group"``
    All six frames are fragment frames of window group 1. No survey frame,
    one frame per spot.
``"groups_take_turns"``
    All six frames are fragment frames; odd frames use group 1, even
    frames group 2, whose windows are shifted by 20 scans and so overlap
    group 1's, as on the LC run.

Group 1 holds three windows centred on the planted ions' apex scans (80,
120, 160 of 240); the first two share scan 100 as their boundary.

Nothing about a real MALDI iprm-PASEF run is measured here: how its survey
and fragment frames share a spot, and whether it uses one window group or
several, are still open.
"""

from __future__ import annotations

import shutil
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any, Dict, List, Tuple

SOURCE = Path(__file__).resolve().parents[1] / "data" / "fixtures" / "synthetic_tims.d"

SURVEY = 0
DIA = 9

# (ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy)
WINDOWS: List[Tuple[int, int, float, float, float]] = [
    (60, 100, 300.0, 1.0, 30.0),
    (100, 140, 500.0, 1.0, 35.0),
    (140, 180, 700.0, 1.0, 40.0),
]
SHIFTED_WINDOWS: List[Tuple[int, int, float, float, float]] = [
    (b + 20, e + 20, mz + 50.0, w, ce) for b, e, mz, w, ce in WINDOWS
]

SHARED_SPOT = (100, 10)
SHARED_SPOT_FRAMES = (1, 7)

SHAPES = ("mixed", "one_group", "groups_take_turns")


def make_iprm_pasef_copy(target: Path, shape: str = "mixed") -> Dict[str, Any]:
    """Copy the committed fixture to ``target`` and give it a type-9 shape.

    Args:
        target: Where the copy goes; must not exist.
        shape: One of :data:`SHAPES`.

    Returns:
        Facts the tests assert against: frame ids by type, each fragment
        frame's window group, and each frame's ``SummedIntensities``.
    """
    if shape not in SHAPES:
        raise ValueError(f"shape must be one of {SHAPES}, got {shape!r}")
    shutil.copytree(SOURCE, target)
    with closing(sqlite3.connect(target / "analysis.tdf")) as conn:
        if shape == "mixed":
            survey, fragment = (1, 2, 3), (4, 5, 6, 7)
            conn.execute(f"UPDATE Frames SET MsMsType = {SURVEY} WHERE Id <= 3")
            conn.execute(f"UPDATE Frames SET MsMsType = {DIA} WHERE Id > 3")
            _add_frame_on_shared_spot(conn)
        else:
            survey, fragment = (), (1, 2, 3, 4, 5, 6)
            conn.execute(f"UPDATE Frames SET MsMsType = {DIA}")

        groups = {
            frame: (2 if shape == "groups_take_turns" and frame % 2 == 0 else 1)
            for frame in fragment
        }
        layouts = {1: WINDOWS, 2: SHIFTED_WINDOWS}
        for group in sorted(set(groups.values())):
            conn.execute(
                "INSERT INTO DiaFrameMsMsWindowGroups (Id) VALUES (?)", (group,)
            )
            conn.executemany(
                "INSERT INTO DiaFrameMsMsWindows (WindowGroup, ScanNumBegin, "
                "ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                [(group, *window) for window in layouts[group]],
            )
        conn.executemany(
            "INSERT INTO DiaFrameMsMsInfo (Frame, WindowGroup) VALUES (?, ?)",
            sorted(groups.items()),
        )
        conn.commit()
        summed = dict(
            conn.execute("SELECT Id, SummedIntensities FROM Frames").fetchall()
        )

    return {
        "survey_frames": survey,
        "fragment_frames": fragment,
        "window_groups": groups,
        "summed_intensities": summed,
    }


def _add_frame_on_shared_spot(conn: sqlite3.Connection) -> None:
    """Frame 7: frame 4's block, one second after frame 6, on frame 1's spot."""
    cols = [row[1] for row in conn.execute("PRAGMA table_info(Frames)")]
    row = dict(zip(cols, conn.execute("SELECT * FROM Frames WHERE Id = 4").fetchone()))
    row.update(Id=7, Time=7.0, MsMsType=DIA)
    conn.execute(
        f"INSERT INTO Frames ({','.join(cols)}) VALUES ({','.join('?' * len(cols))})",
        [row[c] for c in cols],
    )
    mcols = [row[1] for row in conn.execute("PRAGMA table_info(MaldiFrameInfo)")]
    mrow = dict(
        zip(
            mcols,
            conn.execute("SELECT * FROM MaldiFrameInfo WHERE Frame = 1").fetchone(),
        )
    )
    mrow.update(Frame=7)
    conn.execute(
        f"INSERT INTO MaldiFrameInfo ({','.join(mcols)}) "
        f"VALUES ({','.join('?' * len(mcols))})",
        [mrow[c] for c in mcols],
    )
    conn.execute("UPDATE Segments SET LastFrame = 7")
