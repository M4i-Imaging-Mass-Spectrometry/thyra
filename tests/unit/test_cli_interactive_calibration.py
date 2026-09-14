"""``--interactive-calibration`` on a well-formed Bruker ``.d`` (issue #297).

The flag printed nothing at all. Its query was

    SELECT cs.Id, ci.DateTime
    FROM CalibrationState cs
    LEFT JOIN CalibrationInfo ci ON cs.Id = ci.StateId

and ``CalibrationInfo`` has neither ``DateTime`` nor ``StateId`` -- its
columns are ``(Id, CalibrationState, KeyName, Value)``. So every
well-formed file raised ``no such column`` into a bare
``except Exception: return []``, the display returned on the empty list,
and the run continued at exit 0 with no output and no warning. The reader
five hundred lines away had the schema right the whole time.

Repairing the join alone would not have been enough: the per-row
"recalibrated N times" was ``state["id"] - 1``, a dataset property read
off each row, so a three-state file would have claimed state 2 was
recalibrated once and state 3 twice.

The schema here is the one ``tests/unit/readers/test_bruker_calibration.py``
builds and ``BrukerReader._read_calibration_metadata`` queries, which is
this repository's only record of what a real ``calibration.sqlite``
looks like -- no real Bruker calibration file is committed.
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

from thyra.__main__ import _display_calibration_info
from thyra.readers.bruker.timstof.timstof_reader import read_calibration_states

STATES = [
    (1, "cal-uuid-1", "2025-01-01T12:00:00.000+00:00", "timsTOF"),
    (2, "cal-uuid-2", "2025-02-14T09:30:00.000+00:00", "DataAnalysis"),
    (3, "cal-uuid-3", "2025-03-01T09:00:00.000+00:00", "DataAnalysis"),
]


def _make_d(tmp_path: Path, states=STATES) -> Path:
    """A ``.d`` carrying a calibration.sqlite with the vendor's schema."""
    data_path = tmp_path / "sample.d"
    data_path.mkdir()
    conn = sqlite3.connect(data_path / "calibration.sqlite")
    conn.execute(
        "CREATE TABLE CalibrationState "
        "(Id INTEGER PRIMARY KEY, Key TEXT, DateTime TEXT, Source TEXT)"
    )
    conn.execute(
        "CREATE TABLE CalibrationInfo "
        "(Id INTEGER PRIMARY KEY, CalibrationState INTEGER, "
        "KeyName TEXT, Value TEXT)"
    )
    conn.executemany("INSERT INTO CalibrationState VALUES (?, ?, ?, ?)", states)
    conn.commit()
    conn.close()
    return data_path


class TestTheStatesAreRead:
    def test_every_state_comes_back_oldest_first(self, tmp_path):
        states = read_calibration_states(_make_d(tmp_path))
        assert [s["id"] for s in states] == [1, 2, 3]

    def test_the_datetime_is_the_one_in_the_file(self, tmp_path):
        states = read_calibration_states(_make_d(tmp_path))
        assert states[0]["datetime"] == "2025-01-01T12:00:00.000+00:00"
        assert states[-1]["datetime"] == "2025-03-01T09:00:00.000+00:00"

    def test_the_source_comes_back_too(self, tmp_path):
        states = read_calibration_states(_make_d(tmp_path))
        assert [s["source"] for s in states] == [
            "timsTOF",
            "DataAnalysis",
            "DataAnalysis",
        ]

    def test_a_d_without_a_calibration_file_is_simply_empty(self, tmp_path):
        data_path = tmp_path / "bare.d"
        data_path.mkdir()
        assert read_calibration_states(data_path) == []

    def test_the_vendor_file_is_not_written_to(self, tmp_path):
        """No ``-wal``/``-shm`` may appear beside a vendor database."""
        data_path = _make_d(tmp_path)
        cal = data_path / "calibration.sqlite"
        conn = sqlite3.connect(cal)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.commit()
        conn.close()
        for side in (
            cal.with_name(cal.name + "-wal"),
            cal.with_name(cal.name + "-shm"),
        ):
            side.unlink(missing_ok=True)

        read_calibration_states(data_path)

        assert sorted(p.name for p in data_path.iterdir()) == ["calibration.sqlite"]


class TestTheDisplay:
    def _output(self, data_path: Path, use_recalibrated: bool = True) -> str:
        """Everything the display writes to stdout, read inside the block.

        Read *after* the ``with`` instead and this helper only works on
        click >= 8.3.2. Before that, ``isolation()``'s ``_NamedTextIOWrapper``
        inherited ``TextIOWrapper.__del__``, which closes the buffer it wraps;
        restoring ``sys.stdout`` on the way out of the block drops the last
        reference to the wrapper, CPython collects it there and then, and the
        ``BytesIO`` is closed before ``getvalue()`` is reached. 8.3.2 added a
        no-op ``close()`` to stop that. This project declares ``click>=8.1,
        <9``, so every release from 8.1 through 8.3.1 is supported and closes
        the buffer -- measured on 8.1.8, 8.2.0, 8.2.1, 8.3.0 and 8.3.1, against
        8.3.2 and 8.4.2 (what the lockfile resolves, which is why CI never saw
        this) where it stays open.

        The five tests whose subject prints got away with it: ``click.echo``
        files the wrapper in click's ``_default_text_stdout`` cache, a
        ``WeakKeyDictionary`` that stores each stream as its own value and so
        strongly references its own key. Only the prints-nothing case, which
        never calls ``echo``, was left holding a closed buffer.

        The flush is for a future writer that is not ``click.echo`` -- ``echo``
        flushes itself, a bare ``print`` does not, and unflushed text sits in
        the wrapper rather than in the ``BytesIO`` being read.
        """
        runner = CliRunner()
        with runner.isolation() as streams:
            _display_calibration_info(data_path, use_recalibrated)
            sys.stdout.flush()
            return streams[0].getvalue().decode()

    def test_it_prints_something_at_all(self, tmp_path):
        assert "Calibration Information" in self._output(_make_d(tmp_path))

    def test_every_state_is_listed_with_its_datetime(self, tmp_path):
        out = self._output(_make_d(tmp_path))
        for state_id, _, when, source in STATES:
            assert f"State {state_id}: {when} [{source}]" in out

    def test_the_recalibration_count_is_stated_once_for_the_dataset(self, tmp_path):
        out = self._output(_make_d(tmp_path))
        assert out.count("Recalibrated") == 1
        assert "Recalibrated 2 times since acquisition" in out

    def test_a_never_recalibrated_dataset_says_nothing_about_it(self, tmp_path):
        out = self._output(_make_d(tmp_path, states=STATES[:1]))
        assert "Recalibrated" not in out
        assert "State 1:" in out

    def test_the_highest_id_is_the_active_one(self, tmp_path):
        out = self._output(_make_d(tmp_path))
        active = [line for line in out.splitlines() if "active/will be used" in line]
        assert len(active) == 1
        assert "State 3" in active[0]

    def test_nothing_is_printed_without_a_calibration_file(self, tmp_path):
        data_path = tmp_path / "bare.d"
        data_path.mkdir()
        assert self._output(data_path) == ""


@pytest.mark.parametrize("column", ["DateTime", "StateId"])
def test_calibration_info_has_neither_column_the_old_query_named(tmp_path, column):
    """The schema claim the fix rests on, pinned against the fixture."""
    conn = sqlite3.connect(_make_d(tmp_path) / "calibration.sqlite")
    names = {row[1] for row in conn.execute("PRAGMA table_info(CalibrationInfo)")}
    conn.close()
    assert column not in names
    assert names == {"Id", "CalibrationState", "KeyName", "Value"}
