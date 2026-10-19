"""What the timsTOF reader checks about a ``.d`` before Bruker's library opens it.

Two things, both from issue #518. Which analysis database is read (TSF or
TDF) is decided in one place, and the metadata follows it. And a ``.d``
that the library would have to write into, because a database it opens
is in WAL mode, is refused in one line when the folder cannot be
written. None of this needs the vendor library.
"""

import shutil
import sqlite3
from contextlib import closing
from pathlib import Path
from unittest.mock import patch

import pytest

from thyra.errors import ConversionRefused
from thyra.metadata.extractors.bruker_extractor import analysis_database
from thyra.readers.bruker.timstof import timstof_reader
from thyra.readers.bruker.timstof.timstof_reader import (
    BrukerReader,
    folder_is_writable,
    is_wal_database,
)

FIXTURE = Path(__file__).resolve().parents[2] / "data" / "fixtures" / "synthetic_tims.d"


def _sqlite(path: Path, wal: bool) -> Path:
    with closing(sqlite3.connect(path)) as conn:
        conn.execute("CREATE TABLE t (x INTEGER)")
        if wal:
            conn.execute("PRAGMA journal_mode=WAL")
        conn.commit()
    return path


class TestTheWalHeaderCheck:
    def test_a_wal_database_reads_as_wal(self, tmp_path):
        assert is_wal_database(_sqlite(tmp_path / "a.sqlite", wal=True))

    def test_a_rollback_database_does_not(self, tmp_path):
        assert not is_wal_database(_sqlite(tmp_path / "a.sqlite", wal=False))

    def test_a_short_file_does_not(self, tmp_path):
        short = tmp_path / "short.sqlite"
        short.write_bytes(b"SQLite format 3\x00" + b"\x02" * 10)
        assert not is_wal_database(short)

    def test_a_file_that_is_not_sqlite_does_not(self, tmp_path):
        other = tmp_path / "other.bin"
        other.write_bytes(b"\x02" * 200)
        assert not is_wal_database(other)

    def test_a_missing_file_does_not(self, tmp_path):
        assert not is_wal_database(tmp_path / "absent.sqlite")

    def test_reading_the_header_leaves_nothing_beside_the_file(self, tmp_path):
        database = _sqlite(tmp_path / "a.sqlite", wal=True)
        before = sorted(p.name for p in tmp_path.iterdir())

        is_wal_database(database)

        assert sorted(p.name for p in tmp_path.iterdir()) == before


class TestTheWriteProbe:
    def test_a_writable_folder_is_writable_and_left_as_it_was(self, tmp_path):
        (tmp_path / "keep.txt").write_text("x")

        assert folder_is_writable(tmp_path)
        assert [p.name for p in tmp_path.iterdir()] == ["keep.txt"]

    def test_a_folder_that_refuses_a_new_entry_is_not(self, tmp_path):
        with patch.object(Path, "mkdir", side_effect=PermissionError("denied")):
            assert not folder_is_writable(tmp_path)

    def test_a_missing_folder_is_not(self, tmp_path):
        assert not folder_is_writable(tmp_path / "absent")


def _wal_copy(tmp_path: Path, which: str = "calibration.sqlite") -> Path:
    data_path = tmp_path / "copy.d"
    shutil.copytree(FIXTURE, data_path)
    with closing(sqlite3.connect(data_path / which)) as conn:
        conn.execute("PRAGMA journal_mode=WAL")
    return data_path


@pytest.fixture
def no_library():
    """The SDK load, recorded instead of run."""
    with patch.object(BrukerReader, "_initialize_sdk", autospec=True) as init:
        init.side_effect = lambda self: setattr(self, "sdk", None) or setattr(
            self, "handle", None
        )
        yield init


class TestAFolderTheLibraryCannotWrite:
    def test_a_wal_calibration_in_a_read_only_folder_is_refused(
        self, tmp_path, no_library
    ):
        data_path = _wal_copy(tmp_path)
        with patch.object(timstof_reader, "folder_is_writable", return_value=False):
            with pytest.raises(ConversionRefused) as refused:
                BrukerReader(data_path)

        message = str(refused.value)
        assert "calibration.sqlite" in message
        assert "writable folder" in message
        assert "--no-recalibrated" in message
        assert "\n" not in message
        no_library.assert_not_called()

    def test_a_wal_analysis_database_names_no_option(self, tmp_path, no_library):
        data_path = _wal_copy(tmp_path, "analysis.tdf")
        with patch.object(timstof_reader, "folder_is_writable", return_value=False):
            with pytest.raises(ConversionRefused) as refused:
                BrukerReader(data_path, use_recalibrated_state=False)

        assert "analysis.tdf" in str(refused.value)
        assert "--no-recalibrated" not in str(refused.value)
        no_library.assert_not_called()

    def test_without_the_recalibrated_state_calibration_sqlite_is_not_asked(
        self, tmp_path, no_library
    ):
        data_path = _wal_copy(tmp_path)
        with patch.object(
            timstof_reader, "folder_is_writable", return_value=False
        ) as probe:
            BrukerReader(data_path, use_recalibrated_state=False)

        probe.assert_not_called()
        no_library.assert_called_once()

    def test_databases_in_rollback_mode_ask_nothing(self, no_library):
        with patch.object(
            timstof_reader, "folder_is_writable", return_value=False
        ) as probe:
            BrukerReader(FIXTURE)

        probe.assert_not_called()
        no_library.assert_called_once()

    def test_a_writable_folder_converts(self, tmp_path, no_library):
        data_path = _wal_copy(tmp_path)
        before = sorted(p.name for p in data_path.iterdir())

        BrukerReader(data_path)

        no_library.assert_called_once()
        assert sorted(p.name for p in data_path.iterdir()) == before

    def test_a_metadata_only_reader_is_not_checked(self, tmp_path):
        data_path = _wal_copy(tmp_path)
        with patch.object(
            timstof_reader, "folder_is_writable", return_value=False
        ) as probe:
            BrukerReader(data_path, metadata_only=True)

        probe.assert_not_called()


def _tsf_run(tmp_path: Path) -> Path:
    """A minimal TSF ``.d``: no mobility, no spectra decoded."""
    data_path = tmp_path / "tsf_run.d"
    data_path.mkdir()
    with closing(sqlite3.connect(data_path / "analysis.tsf")) as conn:
        conn.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
        conn.executemany(
            "INSERT INTO GlobalMetadata VALUES (?, ?)",
            [
                ("SchemaType", "TSF"),
                ("MzAcqRangeLower", "100.0"),
                ("MzAcqRangeUpper", "1000.0"),
            ],
        )
        conn.execute(
            "CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, "
            "Polarity CHAR(1), ScanMode INTEGER, MsMsType INTEGER, "
            "NumPeaks INTEGER, MzCalibration INTEGER)"
        )
        conn.executemany(
            "INSERT INTO Frames VALUES (?, ?, '+', 9, 0, 10, 1)",
            [(i, i) for i in (1, 2, 3)],
        )
        conn.commit()
    (data_path / "analysis.tsf_bin").write_bytes(b"\x00" * 64)
    return data_path


class TestWhichDatabaseIsRead:
    def test_tsf_wins_when_both_hold_data(self, tmp_path):
        (tmp_path / "analysis.tsf").write_bytes(b"x")
        (tmp_path / "analysis.tdf").write_bytes(b"x")
        assert analysis_database(tmp_path) == ("tsf", tmp_path / "analysis.tsf")

    def test_an_empty_database_gives_way_to_one_with_data(self, tmp_path):
        (tmp_path / "analysis.tsf").write_bytes(b"")
        (tmp_path / "analysis.tdf").write_bytes(b"x")
        assert analysis_database(tmp_path) == ("tdf", tmp_path / "analysis.tdf")

    def test_a_lone_empty_database_is_still_named(self, tmp_path):
        (tmp_path / "analysis.tdf").write_bytes(b"")
        assert analysis_database(tmp_path) == ("tdf", tmp_path / "analysis.tdf")

    def test_neither_is_none(self, tmp_path):
        assert analysis_database(tmp_path) is None

    def test_a_tsf_with_a_stray_empty_tdf_is_stored_as_tsf(self, tmp_path):
        """The F27 case: real TSF folders carry a 0-byte analysis.tdf."""
        from thyra.metadata.document import read_metadata_document

        data_path = _tsf_run(tmp_path)
        (data_path / "analysis.tdf").write_bytes(b"")

        with BrukerReader(data_path, metadata_only=True) as reader:
            assert reader.file_type == "tsf"
            format_specific = reader.get_comprehensive_metadata().format_specific
        document = read_metadata_document(data_path)

        assert format_specific["bruker_format"] == "bruker_tsf"
        assert Path(format_specific["database_path"]).name == "analysis.tsf"
        assert Path(format_specific["binary_file"]).name == "analysis.tsf_bin"
        assert format_specific["ion_mobility"] == {"present": False}
        assert document["ms_analysis"]["ion_mobility"] == {"present": False}
