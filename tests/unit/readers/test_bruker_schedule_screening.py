"""A damaged isolation window is screened out where it is read (issue #498).

The timsTOF reader takes each window's target, width and collision energy
from ``analysis.tdf``. A target must be a finite m/z above 0, a width or
an energy finite and not negative. A row that breaks this is left out of
the schedule with one warning, so it reaches neither the MS/MS table nor
the metadata document. The input is the synthetic six-frame ``.d`` with a
three-window PASEF schedule written into it.
"""

from __future__ import annotations

import json
import logging
import shutil
import sqlite3
from pathlib import Path

import pytest

from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

_TDF = Path(__file__).resolve().parents[2] / "data" / "fixtures" / "synthetic_tims.d"

#: (IsolationMz, IsolationWidth, CollisionEnergy, ScanNumBegin, ScanNumEnd)
_WINDOWS = [
    (936.578, 1.0, 49.6, 20, 60),
    (353.320, 1.0, 35.2, 80, 130),
    (313.275, 1.0, 34.3, 150, 200),
]
_INF = float("inf")


def _damaged(tmp_path: Path, column: int, value: float) -> Path:
    """The fixture as a PASEF run whose second window has one bad value."""
    run = tmp_path / "run.d"
    shutil.copytree(_TDF, run)
    rows = [list(window) for window in _WINDOWS]
    rows[1][column] = value
    with sqlite3.connect(run / "analysis.tdf") as db:
        db.execute("UPDATE Frames SET MsMsType = 8")
        db.execute(
            "CREATE TABLE PasefFrameMsMsInfo (Frame INTEGER, ScanNumBegin INTEGER, "
            "ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, "
            "CollisionEnergy REAL, Precursor INTEGER)"
        )
        for (frame,) in db.execute("SELECT Id FROM Frames").fetchall():
            db.executemany(
                "INSERT INTO PasefFrameMsMsInfo VALUES (?, ?, ?, ?, ?, ?, NULL)",
                [(frame, r[3], r[4], r[0], r[1], r[2]) for r in rows],
            )
    return run


_CASES = [
    pytest.param(0, _INF, id="infinite-target"),
    pytest.param(0, -300.0, id="negative-target"),
    pytest.param(1, -1.0, id="negative-width"),
    pytest.param(1, _INF, id="infinite-width"),
    pytest.param(2, _INF, id="infinite-energy"),
    pytest.param(0, float("nan"), id="nan-target"),
]


def _schedule(run: Path):
    reader = BrukerReader(run, metadata_only=True)
    try:
        return reader.get_fragmentation(), reader.has_precursor_spectra
    finally:
        reader.close()


class TestTheReaderScreensEachWindow:
    @pytest.mark.parametrize("column, value", _CASES)
    def test_the_bad_window_is_left_out_with_a_warning(
        self, tmp_path, thyra_logs, column, value
    ):
        reader = BrukerReader(_damaged(tmp_path, column, value), metadata_only=True)
        try:
            with thyra_logs(
                "thyra.readers.bruker.timstof.timstof_reader", logging.WARNING
            ) as logs:
                schedule = reader.get_fragmentation()
        finally:
            reader.close()

        assert schedule is not None
        assert sorted(w.target for w in schedule.windows) == [313.275, 936.578]
        assert schedule.constant_across_pixels
        screened = [r for r in logs if "isolation windows hold" in r.getMessage()]
        assert len(screened) == 1
        assert screened[0].getMessage().startswith("1 of 3")

    @pytest.mark.parametrize("column, value", _CASES)
    def test_a_window_left_out_with_its_scans_stops_the_split(
        self, tmp_path, column, value
    ):
        """Its scans still hold its ions, which no window would claim (D2)."""
        _, separable = _schedule(_damaged(tmp_path, column, value))

        assert separable is False

    def test_a_negative_energy_is_kept(self, tmp_path, thyra_logs):
        """The schema allows either sign, as signed energies are (#498 review)."""
        run = _damaged(tmp_path, 2, -35.2)
        with thyra_logs(
            "thyra.readers.bruker.timstof.timstof_reader", logging.WARNING
        ) as logs:
            schedule, separable = _schedule(run)

        assert sorted(w.collision_energy for w in schedule.windows) == [
            -35.2,
            34.3,
            49.6,
        ]
        assert separable is True
        assert not [r for r in logs if "isolation windows hold" in r.getMessage()]

    def test_a_table_of_null_targets_is_still_the_one_used(self, tmp_path):
        """A NULL target is screened and counted; the next table is not asked."""
        run = _damaged(tmp_path, 1, 1.0)
        with sqlite3.connect(run / "analysis.tdf") as db:
            db.execute("UPDATE PasefFrameMsMsInfo SET IsolationMz = NULL")
            db.executemany(
                "INSERT INTO FrameMsMsInfo (Frame, Parent, TriggerMass, "
                "IsolationWidth, PrecursorCharge, CollisionEnergy) "
                "VALUES (?, NULL, 500.0, 1.0, NULL, 30.0)",
                [(f,) for f in range(1, 7)],
            )
        schedule, _ = _schedule(run)

        assert schedule is not None and schedule.windows == ()

    def test_isomers_are_listed_by_mobility(self, tmp_path):
        """By first scan descending, so 1/K0 ascending, as the table ranks them."""
        run = _damaged(tmp_path, 0, 313.275)
        schedule, _ = _schedule(run)

        isomers = [w for w in schedule.windows if w.target == 313.275]
        assert [w.scan_begin for w in isomers] == [150, 80]

    def test_a_good_schedule_is_untouched(self, tmp_path, thyra_logs):
        reader = BrukerReader(_damaged(tmp_path, 1, 1.0), metadata_only=True)
        try:
            with thyra_logs(
                "thyra.readers.bruker.timstof.timstof_reader", logging.WARNING
            ) as logs:
                schedule = reader.get_fragmentation()
        finally:
            reader.close()

        assert schedule is not None and len(schedule.windows) == 3
        assert not [r for r in logs if "isolation windows hold" in r.getMessage()]

    @pytest.mark.parametrize(
        "begin, end",
        [(-20, 60), (60, 20), (20, float("inf"))],
        ids=["negative", "reversed", "infinite"],
    )
    def test_a_bad_scan_range_is_left_out(self, tmp_path, begin, end):
        """A window with no usable range occupies no known scans: the split stays."""
        run = _damaged(tmp_path, 3, begin)
        with sqlite3.connect(run / "analysis.tdf") as db:
            db.execute(
                "UPDATE PasefFrameMsMsInfo SET ScanNumEnd = ? WHERE IsolationMz = ?",
                (end, 353.320),
            )
        schedule, separable = _schedule(run)

        assert schedule is not None
        assert sorted(w.target for w in schedule.windows) == [313.275, 936.578]
        assert separable is True

    def test_a_table_with_no_usable_row_is_still_the_one_used(
        self, tmp_path, thyra_logs
    ):
        """The frames use the first table with rows; the next is not asked."""
        run = _damaged(tmp_path, 1, 1.0)
        with sqlite3.connect(run / "analysis.tdf") as db:
            db.execute("UPDATE PasefFrameMsMsInfo SET IsolationWidth = -1")
        reader = BrukerReader(run, metadata_only=True)
        try:
            with thyra_logs(
                "thyra.readers.bruker.timstof.timstof_reader", logging.INFO
            ) as logs:
                schedule = reader.get_fragmentation()
        finally:
            reader.close()

        assert schedule is not None and schedule.ms_level == 2
        assert schedule.windows == ()
        assert not [r for r in logs if "carries a precursor" in r.getMessage()]

    def test_a_single_precursor_table_is_screened(self, tmp_path):
        run = tmp_path / "run.d"
        shutil.copytree(_TDF, run)
        with sqlite3.connect(run / "analysis.tdf") as db:
            db.execute("UPDATE Frames SET MsMsType = 2")
            db.executemany(
                "INSERT INTO FrameMsMsInfo (Frame, Parent, TriggerMass, "
                "IsolationWidth, PrecursorCharge, CollisionEnergy) "
                "VALUES (?, NULL, ?, 1.0, NULL, ?)",
                [(f, 500.0 if f % 2 else float("inf"), 30.0) for f in range(1, 7)],
            )
        reader = BrukerReader(run, metadata_only=True)
        try:
            schedule = reader.get_fragmentation()
        finally:
            reader.close()

        assert schedule is not None
        assert [w.target for w in schedule.windows] == [500.0]


class TestTheStoreKeepsItsMetadata:
    def test_a_negative_width_no_longer_drops_msi_metadata(self, tmp_path):
        """The finding's worst case: exit 0 and no metadata block at all."""
        import spatialdata

        from thyra.convert import convert_msi

        out = tmp_path / "out.zarr"
        assert convert_msi(
            str(_damaged(tmp_path, 1, -1.0)),
            str(out),
            include_optical=False,
            reader_options={"tdf_spectrum": "scan_sum"},
        )
        tables = spatialdata.read_zarr(out).tables

        summed = tables["msi_dataset_z0"]
        fragmentation = summed.uns["msi_metadata"]["ms_analysis"]["fragmentation"]
        windows = fragmentation["windows"]
        if isinstance(windows, str):  # stored as a JSON string
            windows = json.loads(windows)
        assert [w["target"] for w in windows] == [313.275, 936.578]
        # The left-out window still occupies scans 80-130, so the fragments
        # are not split by precursor (D2), rather than split with a gap.
        assert "msi_dataset_z0_msms" not in tables
