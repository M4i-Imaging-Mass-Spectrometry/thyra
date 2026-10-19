"""Which calibration Bruker's library applies, checked through the library itself.

Each test edits a copy of ``tests/data/fixtures/synthetic_tims.d`` (the
fixture stays as committed) and opens it through the bundled Bruker
library, so what the reader records is compared with what the library
actually did to the m/z and the 1/K0. Anywhere the library cannot be
loaded these tests skip rather than fail. Issue #518.
"""

from __future__ import annotations

import shutil
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Iterable

import numpy as np
import pytest

from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
from thyra.utils.bruker_exceptions import SDKError

pytestmark = pytest.mark.integration

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "fixtures" / "synthetic_tims.d"

#: Moves the m/z of the fixture's only state by about +20 ppm, so whether
#: the library applied it shows in the m/z.
SHIFT_STATE_1 = "UPDATE FrameMzCalibration SET C1 = C1 * 1.00002"

#: Gives frame 6 a TIMS calibration row of its own in state 1, with the
#: kind of drift a real state shows between frames.
OWN_ROW_FOR_FRAME_6 = (
    "INSERT INTO TimsCalibration SELECT CalibrationState, 2, ModelType, C0, C1, "
    "C2, C3, C4, C5, C6 * 1.007, C7 * 1.007, C8, C9 FROM TimsCalibration "
    "WHERE Id = 1",
    "UPDATE TimsCalibrationMapping SET TimsCalibrationId = 2 WHERE FrameId = 6",
)

#: A second state, a copy of the first that moves m/z by about +20 ppm.
SECOND_STATE = (
    "INSERT INTO CalibrationState VALUES "
    "(2, 'k2', '2026-03-01T00:00:00+00:00', 'RecalTool', NULL, NULL, 1)",
    "INSERT INTO DigitizerConstants SELECT 2, DigitizerTimebase, DigitizerDelay "
    "FROM DigitizerConstants WHERE CalibrationState = 1",
    "INSERT INTO FrameMzCalibration SELECT 2, Id, ModelType, C0, C1 * 1.00002, "
    "C2, C3, C4 FROM FrameMzCalibration WHERE CalibrationState = 1",
    "INSERT INTO FrameMzCalibrationMapping SELECT 2, FrameId, FrameMzCalibrationId "
    "FROM FrameMzCalibrationMapping WHERE CalibrationState = 1",
    "INSERT INTO TimsCalibration SELECT 2, Id, ModelType, C0, C1, C2, C3, C4, C5, "
    "C6, C7, C8, C9 FROM TimsCalibration WHERE CalibrationState = 1",
    "INSERT INTO TimsCalibrationMapping SELECT 2, FrameId, TimsCalibrationId "
    "FROM TimsCalibrationMapping WHERE CalibrationState = 1",
)


def _copy(tmp_path: Path, statements: Iterable[str] = ()) -> Path:
    """A copy of the fixture with ``statements`` run on its calibration.sqlite."""
    data_path = tmp_path / "copy.d"
    shutil.copytree(FIXTURE, data_path)
    with closing(sqlite3.connect(data_path / "calibration.sqlite")) as conn:
        for statement in statements:
            conn.execute(statement)
        conn.commit()
    return data_path


def _open(data_path: Path, use_recalibrated_state: bool = True) -> BrukerReader:
    try:
        return BrukerReader(
            data_path,
            tdf_spectrum="scan_sum",
            use_recalibrated_state=use_recalibrated_state,
        )
    except (SDKError, OSError) as exc:  # the vendor library is not loadable here
        pytest.skip(f"Bruker library not loadable on this platform: {exc}")


def _first_mz(data_path: Path, use_recalibrated_state: bool) -> np.ndarray:
    with _open(data_path, use_recalibrated_state) as reader:
        return reader._read_frame_spectrum(1)[0]


def _ppm(on: np.ndarray, off: np.ndarray) -> float:
    return float(np.median(on / off - 1.0) * 1e6)


class TestTheOptionReachesTheLibrary:
    def test_the_m_z_differ_with_the_option_on_and_off(self, tmp_path):
        """Fails when the option stops reaching the library's open call (F95)."""
        data_path = _copy(tmp_path, [SHIFT_STATE_1])

        on, off = _first_mz(data_path, True), _first_mz(data_path, False)

        assert _ppm(on, off) == pytest.approx(20.0, abs=0.5)


class TestTheStateThatIsApplied:
    def test_a_state_mapping_one_frame_too_few_is_not_named(self, tmp_path):
        """The library skips it, so the store names the analysis database (F35)."""
        data_path = _copy(
            tmp_path,
            [SHIFT_STATE_1, "DELETE FROM FrameMzCalibrationMapping WHERE FrameId = 6"],
        )

        with _open(data_path) as reader:
            on = reader._read_frame_spectrum(1)[0]
            applied = reader.get_applied_mz_calibration()
            metadata = reader._calibration_metadata
        off = _first_mz(data_path, False)

        assert _ppm(on, off) == pytest.approx(0.0, abs=0.01)
        assert applied == {
            "use_recalibrated_state": True,
            "calibration": "analysis.tdf",
        }
        assert metadata["calibration_id"] is None
        assert metadata["latest_calibration_id"] == 1

    def test_a_newest_state_that_does_not_map_falls_back_to_the_one_before(
        self, tmp_path
    ):
        data_path = _copy(
            tmp_path,
            [
                *SECOND_STATE,
                "DELETE FROM FrameMzCalibrationMapping "
                "WHERE CalibrationState = 2 AND FrameId = 6",
                "UPDATE FrameMzCalibration SET C1 = C1 * 1.00001 "
                "WHERE CalibrationState = 1",
            ],
        )

        with _open(data_path) as reader:
            on = reader._read_frame_spectrum(1)[0]
            applied = reader.get_applied_mz_calibration()
            metadata = reader._calibration_metadata
        off = _first_mz(data_path, False)

        # State 1's +10 ppm, not state 2's +20 ppm.
        assert _ppm(on, off) == pytest.approx(10.0, abs=0.5)
        assert applied["calibration_state_id"] == 1
        assert metadata["calibration_id"] == 1
        assert metadata["latest_calibration_id"] == 2

    def test_a_missing_calibration_info_table_costs_only_the_version(self, tmp_path):
        """F37: the state is still named, and its m/z are applied."""
        data_path = _copy(tmp_path, [*SECOND_STATE, "DROP TABLE CalibrationInfo"])

        with _open(data_path) as reader:
            on = reader._read_frame_spectrum(1)[0]
            applied = reader.get_applied_mz_calibration()
            metadata = reader._calibration_metadata
        off = _first_mz(data_path, False)

        assert _ppm(on, off) == pytest.approx(20.0, abs=0.5)
        assert applied["calibration_state_id"] == 2
        assert metadata["calibration_id"] == 2
        assert metadata["calibration_software_version"] is None
        assert metadata["recalibrated"] is True
        assert metadata["num_calibration_versions"] == 2

    def test_a_file_with_no_state_names_the_analysis_database(self, tmp_path):
        """F61: the library applies the analysis database's, which is known."""
        data_path = _copy(tmp_path, [SHIFT_STATE_1, "DELETE FROM CalibrationState"])

        with _open(data_path) as reader:
            on = reader._read_frame_spectrum(1)[0]
            applied = reader.get_applied_mz_calibration()
        off = _first_mz(data_path, False)

        assert _ppm(on, off) == pytest.approx(0.0, abs=0.01)
        assert applied == {
            "use_recalibrated_state": True,
            "calibration": "analysis.tdf",
        }

    def test_the_library_agrees_with_the_state_named(self, tmp_path):
        data_path = _copy(tmp_path, [*SECOND_STATE])

        with _open(data_path) as reader:
            assert reader.sdk.has_recalibrated_state(reader.handle) is True
            assert reader._calibration_metadata["calibration_id"] == 2
        with _open(data_path, use_recalibrated_state=False) as reader:
            assert reader.sdk.has_recalibrated_state(reader.handle) is False
            assert reader._calibration_metadata["calibration_id"] is None


class TestEachFrameTakesItsOwnMobility:
    """F14: a state that calibrates frames apart gives each frame its own 1/K0."""

    def test_a_frame_with_its_own_row_takes_its_own_1_over_k0(self, tmp_path):
        data_path = _copy(tmp_path, OWN_ROW_FOR_FRAME_6)

        with _open(data_path) as reader:
            axis = reader.get_mobility_axis().values
            frames = {frame.frame_id: frame for frame in reader.iter_frame_scans()}
            own = reader.sdk.scannum_to_oneoverk0(
                reader.handle, 6, np.arange(axis.size, dtype=np.float64)
            )
            first = reader.sdk.scannum_to_oneoverk0(
                reader.handle, 1, np.arange(axis.size, dtype=np.float64)
            )
            sixth = frames[6].mobility_points()
            one = frames[1].mobility_points()
            span = reader.get_mobility_value_range()

        assert not np.allclose(own, axis)
        np.testing.assert_array_equal(sixth[1], own[frames[6].scans])
        np.testing.assert_array_equal(one[1], first[frames[1].scans])
        np.testing.assert_array_equal(first, axis)
        # The range covers frame 6's ramp too, so no point is binned past an edge.
        assert span[0] <= own.min() and span[1] >= own.max()
        assert span[0] <= axis.min() and span[1] >= axis.max()

    def test_the_iterator_and_the_frame_records_agree(self, tmp_path):
        data_path = _copy(tmp_path, OWN_ROW_FOR_FRAME_6)

        with _open(data_path) as reader:
            from_records = {
                frame.coords: frame.mobility_points()[1]
                for frame in reader.iter_frame_scans()
            }
            from_iterator = {
                coords: mobility
                for coords, _mz, mobility, _i in reader.iter_mobility_spectra()
            }

        assert from_records.keys() == from_iterator.keys()
        for coords, mobility in from_records.items():
            np.testing.assert_array_equal(mobility, from_iterator[coords])

    def test_with_the_option_off_the_analysis_row_applies_to_every_frame(
        self, tmp_path
    ):
        data_path = _copy(tmp_path, OWN_ROW_FOR_FRAME_6)

        with _open(data_path, use_recalibrated_state=False) as reader:
            axis = reader.get_mobility_axis().values
            sixth = next(f for f in reader.iter_frame_scans() if f.frame_id == 6)
            mobility = sixth.mobility_points()[1]
            assert not reader._mobility_calibrated_per_frame()

        np.testing.assert_array_equal(mobility, axis[sixth.scans])

    def test_the_fixture_as_committed_keeps_one_axis(self):
        """One row for every frame: nothing is converted per frame."""
        with _open(FIXTURE) as reader:
            assert not reader._mobility_calibrated_per_frame()
            assert reader.get_mobility_value_range() == (
                float(reader.get_mobility_axis().values.min()),
                float(reader.get_mobility_axis().values.max()),
            )


class TestTheRecordedTimsCalibration:
    """F36: the block names the row the library applied to the axis frame."""

    def test_it_follows_the_option(self, tmp_path):
        data_path = _copy(
            tmp_path, ["UPDATE TimsCalibration SET C6 = C6 * 1.02, C7 = C7 + 3.0"]
        )
        with closing(sqlite3.connect(data_path / "calibration.sqlite")) as conn:
            state_c6_c7 = conn.execute("SELECT C6, C7 FROM TimsCalibration").fetchone()
        with closing(sqlite3.connect(data_path / "analysis.tdf")) as conn:
            analysis_c6_c7 = conn.execute(
                "SELECT C6, C7 FROM TimsCalibration"
            ).fetchone()

        with _open(data_path) as reader:
            on = reader.get_mobility_axis()
        with _open(data_path, use_recalibrated_state=False) as reader:
            off = reader.get_mobility_axis()

        assert tuple(on.calibration["coefficients"][6:8]) == pytest.approx(
            state_c6_c7, rel=1e-15
        )
        assert on.calibration["stored_in"] == "calibration.sqlite"
        assert on.calibration["calibration_state_id"] == 1
        assert on.calibration["frame_id"] == 1
        assert tuple(off.calibration["coefficients"][6:8]) == pytest.approx(
            analysis_c6_c7, rel=1e-15
        )
        assert off.calibration["stored_in"] == "analysis.tdf"
        assert "calibration_state_id" not in off.calibration
        assert on.values[0] != off.values[0]

    def test_it_counts_the_rows_the_frames_take(self, tmp_path):
        data_path = _copy(tmp_path, OWN_ROW_FOR_FRAME_6)

        with _open(data_path) as reader:
            assert reader.get_mobility_axis().calibration["n_frame_calibrations"] == 2
        with _open(FIXTURE) as reader:
            assert reader.get_mobility_axis().calibration["n_frame_calibrations"] == 1
