"""timsTOF frame positions that are not whole numbers (issue #527, F132).

``int()`` cut a REAL just below a whole number down, so frame 2 of the
synthetic run, moved to ``x = 100.9999999``, was summed into the pixel of
frame 1. A frame whose position is text, NULL or not whole is now skipped
and counted with the frames that have no coordinates, as text and NULL
already were, and the region map, region boxes and metadata bounds leave
it out too.
"""

import logging
import shutil
import sqlite3
from pathlib import Path

import pytest
import spatialdata

from thyra.convert import convert_msi
from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

_TDF = Path(__file__).resolve().parents[2] / "data" / "fixtures" / "synthetic_tims.d"


def _run_with_x(tmp_path: Path, frame: int, x) -> Path:
    run = tmp_path / "run.d"
    shutil.copytree(_TDF, run)
    con = sqlite3.connect(str(run / "analysis.tdf"))
    con.execute("UPDATE MaldiFrameInfo SET XIndexPos = ? WHERE Frame = ?", (x, frame))
    con.commit()
    con.close()
    return run


@pytest.mark.parametrize("x", [100.9999999, "101a", None, 1.0e19])
def test_a_frame_whose_position_is_not_whole_is_skipped_and_counted(
    tmp_path, thyra_logs, x
):
    run = _run_with_x(tmp_path, 2, x)
    out = tmp_path / "out.zarr"
    with thyra_logs("thyra.readers", logging.WARNING) as records:
        assert convert_msi(
            str(run),
            str(out),
            pixel_size_um=20.0,
            include_optical=False,
            reader_options={"tdf_spectrum": "scan_sum"},
        )

    table = next(iter(spatialdata.read_zarr(out).tables.values()))
    assert table.n_obs == 5
    # Frame 2 sat at x 101; no other frame took its pixel.
    pixels = sorted(zip(table.obs["x"].tolist(), table.obs["y"].tolist()))
    assert len(set(pixels)) == 5
    assert "Dropped 1 of 6 frames with no usable coordinates" in records.text


def test_the_region_boxes_leave_the_frame_out(tmp_path):
    run = _run_with_x(tmp_path, 3, 1.0e19)
    reader = BrukerReader(run, metadata_only=True)
    try:
        boxes = [r["bounds"] for r in reader.get_region_info() or []]
        essential = reader.get_essential_metadata()
    finally:
        reader.close()
    assert all(max(box) < 10 for box in boxes)
    assert essential.dimensions[0] <= 3


def test_a_frame_whose_region_number_is_not_whole_is_skipped_too(tmp_path, thyra_logs):
    """One rule for the frame loop and the region map and list."""
    run = tmp_path / "run.d"
    shutil.copytree(_TDF, run)
    con = sqlite3.connect(str(run / "analysis.tdf"))
    con.execute("UPDATE MaldiFrameInfo SET RegionNumber = 0.5 WHERE Frame = 2")
    con.commit()
    con.close()
    out = tmp_path / "out.zarr"
    with thyra_logs("thyra.readers", logging.WARNING) as records:
        assert convert_msi(
            str(run),
            str(out),
            pixel_size_um=20.0,
            include_optical=False,
            reader_options={"tdf_spectrum": "scan_sum"},
        )

    table = next(iter(spatialdata.read_zarr(out).tables.values()))
    assert table.n_obs == 5
    assert -1 not in set(table.obs["region_number"])
    assert "Dropped 1 of 6 frames with no usable coordinates" in records.text
