"""Thyra reads a vendor ``.d`` without writing to it (issue #290).

Every sqlite open into a Bruker ``.d`` except the reader's own shared
connection was ``sqlite3.connect(path)``, which opens **read-write**.
Three consequences, all measured:

- a ``.tdf`` in WAL mode gained a ``-wal`` and a ``-shm`` next to it on
  the first read, inside the user's acquisition directory;
- a database another program held open lost to the lock, and
  ``_get_frame_count`` reported the failure as **zero frames** -- which
  is indistinguishable from a genuinely empty acquisition, so the caller
  got the empty-conversion refusal blaming their data. The count is
  cached, so the zero was never retried;
- ``with sqlite3.connect(...)`` does not close the handle. It commits a
  transaction and leaves the connection open, so four of the sites
  leaked one per call.

The locked-database tests hold ``BEGIN EXCLUSIVE`` on a second
connection, which is what DataAnalysis does to an open acquisition.
"""

from __future__ import annotations

import sqlite3
from contextlib import closing
from pathlib import Path

import pytest

from thyra.errors import ConversionRefused
from thyra.readers.bruker.timstof.timstof_reader import (
    _get_frame_coordinates,
    _get_frame_count,
)
from thyra.readers.bruker.vendor_db import open_read_only, read_only_uri

N_FRAMES = 6


def _make_tdf(tmp_path: Path, *, wal: bool = False) -> Path:
    """A minimal ``analysis.tdf`` with a Frames table and MALDI positions."""
    db = tmp_path / "analysis.tdf"
    with closing(sqlite3.connect(db)) as conn:
        if wal:
            conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, NumPeaks INTEGER)")
        conn.execute(
            "CREATE TABLE MaldiFrameInfo "
            "(Frame INTEGER PRIMARY KEY, XIndexPos INTEGER, YIndexPos INTEGER)"
        )
        conn.executemany(
            "INSERT INTO Frames VALUES (?, ?)",
            [(i, 10 * i) for i in range(1, N_FRAMES + 1)],
        )
        conn.executemany(
            "INSERT INTO MaldiFrameInfo VALUES (?, ?, ?)",
            [(i, i, 0) for i in range(1, N_FRAMES + 1)],
        )
        conn.commit()
    for side in ("-wal", "-shm"):
        db.with_name(db.name + side).unlink(missing_ok=True)
    return db


class TestTheUri:
    def test_a_drive_letter_path_is_read_only(self, tmp_path):
        uri = read_only_uri(tmp_path / "analysis.tdf")
        assert uri.startswith("file:")
        assert "mode=ro" in uri and "immutable=1" in uri

    def test_immutable_can_be_declined(self, tmp_path):
        # Only the query string: the tmp_path carries the test's own name.
        query = read_only_uri(tmp_path / "analysis.tdf", immutable=False).rsplit("?", 1)
        assert query[1] == "mode=ro"

    def test_a_unc_path_keeps_an_empty_authority(self):
        """``file://server/share`` is rejected as an invalid authority.

        A mapped network drive resolves to a UNC path on Windows, which
        is where this lab's data lives.
        """
        uri = read_only_uri(r"\\server\share\run.d\analysis.tdf")
        assert uri.startswith("file:////server/share/")


class TestTheVendorFileIsNotWrittenTo:
    def test_no_side_files_appear_beside_a_wal_database(self, tmp_path):
        db = _make_tdf(tmp_path, wal=True)

        assert _get_frame_count(db) == N_FRAMES

        assert sorted(p.name for p in tmp_path.iterdir()) == ["analysis.tdf"]

    def test_the_handle_is_closed(self, tmp_path):
        """``with sqlite3.connect(...)`` commits but does not close."""
        db = _make_tdf(tmp_path)
        opened = []
        real = sqlite3.connect

        def spy(*args, **kwargs):
            conn = real(*args, **kwargs)
            opened.append(conn)
            return conn

        import thyra.readers.bruker.vendor_db as vendor_db

        original = vendor_db.sqlite3.connect
        vendor_db.sqlite3.connect = spy
        try:
            _get_frame_count(db)
        finally:
            vendor_db.sqlite3.connect = original

        assert len(opened) == 1
        with pytest.raises(sqlite3.ProgrammingError):
            opened[0].execute("SELECT 1")


class TestALockedDatabase:
    """What DataAnalysis holding the acquisition open looks like."""

    def test_the_frame_count_is_still_read(self, tmp_path):
        db = _make_tdf(tmp_path)
        with closing(sqlite3.connect(db)) as writer:
            writer.execute("BEGIN EXCLUSIVE")
            assert _get_frame_count(db) == N_FRAMES
            writer.rollback()

    def test_coordinates_are_still_read(self, tmp_path):
        db = _make_tdf(tmp_path)
        with closing(sqlite3.connect(db)) as writer:
            writer.execute("BEGIN EXCLUSIVE")
            assert _get_frame_coordinates(db, 3) == (3, 0, 0)
            writer.rollback()

    def test_a_plain_read_only_open_would_not_have_been_enough(self, tmp_path):
        """``mode=ro`` still waits on the lock; ``immutable=1`` is what
        steps past it. This pins *why* the helper defaults that way."""
        db = _make_tdf(tmp_path)
        with closing(sqlite3.connect(db)) as writer:
            writer.execute("BEGIN EXCLUSIVE")
            with pytest.raises(sqlite3.OperationalError, match="locked"):
                with closing(open_read_only(db, immutable=False, timeout=0.1)) as conn:
                    conn.execute("SELECT COUNT(*) FROM Frames").fetchone()
            writer.rollback()


class TestAnUnreadableDatabaseIsRefusedNotCountedAsZero:
    def test_a_missing_file_refuses(self, tmp_path):
        with pytest.raises(ConversionRefused):
            _get_frame_count(tmp_path / "absent.tdf")

    def test_a_file_with_no_frames_table_refuses(self, tmp_path):
        db = tmp_path / "analysis.tdf"
        with closing(sqlite3.connect(db)) as conn:
            conn.execute("CREATE TABLE Other (Id INTEGER)")
            conn.commit()
        with pytest.raises(ConversionRefused, match="frame table"):
            _get_frame_count(db)

    def test_the_message_names_the_likely_cause(self, tmp_path):
        with pytest.raises(ConversionRefused) as excinfo:
            _get_frame_count(tmp_path / "absent.tdf")
        assert "DataAnalysis" in str(excinfo.value)

    def test_an_empty_acquisition_still_counts_zero(self, tmp_path):
        """Zero must stay a real answer, distinguishable from a failure."""
        db = tmp_path / "analysis.tdf"
        with closing(sqlite3.connect(db)) as conn:
            conn.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY)")
            conn.commit()
        assert _get_frame_count(db) == 0
