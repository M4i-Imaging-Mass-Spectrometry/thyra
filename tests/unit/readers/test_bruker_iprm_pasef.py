"""Type-9 fragment frames: the DIA window tables, and survey frames on a spot.

PASEF frames (``MsMsType = 8``) list their windows per frame in
``PasefFrameMsMsInfo``. Type-9 frames -- diaPASEF, and reportedly MALDI
iprm-PASEF -- name a window group in ``DiaFrameMsMsInfo``, and the group's windows live
in ``DiaFrameMsMsWindows``. Such a run can also hold survey frames
(``MsMsType = 0``) on the same raster spots as its fragment frames.

The schedule is read from the DIA tables (D37). Survey frames among the
fragment frames still make it count as varying, on purpose: the precursor
split cuts every frame by the windows' scan ranges, and a survey frame
isolated nothing. Keeping a spot's survey and fragment frames apart is a
later step; the tests that pin it are marked ``xfail(strict=True)``, so
they turn into failures the day it lands and must then be unmarked.

These need only the database, or a stand-in for the vendor library, so
they run everywhere. The fixture shapes are described in
``tests/fixtures/iprm_pasef_builder.py``.
"""

import logging
import sqlite3
from collections import defaultdict
from contextlib import closing

import pytest

from tests.fixtures.fake_tims_sdk import FakeSDK
from tests.fixtures.iprm_pasef_builder import (
    DIA,
    SHARED_SPOT_FRAMES,
    SHIFTED_WINDOWS,
    SURVEY,
    WINDOWS,
    make_iprm_pasef_copy,
)
from thyra.core.msms import COLLISION_INDUCED_DISSOCIATION_ACCESSION
from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
from thyra.readers.bruker.vendor_db import open_read_only

LATER_STEP = (
    "survey and fragment frames on one spot are not kept apart yet; this "
    "pins the behaviour of that later change"
)


def _copy(tmp_path, shape="mixed"):
    target = tmp_path / f"{shape}.d"
    make_iprm_pasef_copy(target, shape)
    return target


def _open(path):
    return BrukerReader(path, metadata_only=True)


def _schedule(path):
    reader = _open(path)
    try:
        return reader.get_fragmentation(), reader.has_precursor_spectra
    finally:
        reader.close()


def _msms_types(path):
    with closing(open_read_only(path / "analysis.tdf")) as conn:
        return dict(conn.execute("SELECT Id, MsMsType FROM Frames").fetchall())


def _expected(windows):
    return [(b, e, mz, ce) for b, e, mz, _, ce in windows]


def _got(schedule):
    return [
        (w.scan_begin, w.scan_end, w.target, w.collision_energy)
        for w in schedule.windows
    ]


class TestTheFixture:
    """Guards against the tests below passing for the wrong reason."""

    def test_one_spot_holds_a_survey_and_a_fragment_frame(self, tmp_path):
        types = _msms_types(_copy(tmp_path))
        first, second = SHARED_SPOT_FRAMES
        assert (types[first], types[second]) == (SURVEY, DIA)

    def test_both_frames_of_that_spot_reach_the_frame_loop(self, tmp_path):
        """Today's behaviour, measured: the reason the pixel row is summed."""
        reader = _open(_copy(tmp_path))
        try:
            by_coords = defaultdict(list)
            for frame_id, coords in reader._iter_frames():
                by_coords[coords].append(frame_id)
        finally:
            reader.close()
        assert sorted(by_coords[(0, 0, 0)]) == list(SHARED_SPOT_FRAMES)


class TestSchedule:
    def test_the_type_9_windows_are_read(self, tmp_path):
        """Target m/z, half-width offsets, energy and scan range of each window."""
        schedule, _ = _schedule(_copy(tmp_path))

        assert schedule.ms_level == 2
        assert _got(schedule) == _expected(WINDOWS)
        for window, (_, _, _, width, _) in zip(schedule.windows, WINDOWS):
            assert window.lower_offset == pytest.approx(width / 2)
            assert window.upper_offset == pytest.approx(width / 2)
        assert schedule.dissociation_accession == (
            COLLISION_INDUCED_DISSOCIATION_ACCESSION
        )

    def test_a_damaged_dia_window_is_left_out(self, tmp_path, thyra_logs):
        """Issue #498: the DIA tables are screened as PASEF's are."""
        path = _copy(tmp_path, "one_group")
        with sqlite3.connect(path / "analysis.tdf") as db:
            db.execute(
                "UPDATE DiaFrameMsMsWindows SET CollisionEnergy = -1 "
                "WHERE IsolationMz = 500.0"
            )
        with thyra_logs(
            "thyra.readers.bruker.timstof.timstof_reader", logging.WARNING
        ) as logs:
            schedule, _ = _schedule(path)

        assert [w.target for w in schedule.windows] == [300.0, 700.0]
        assert [r for r in logs if "isolation windows hold" in r.getMessage()]

    def test_survey_frames_keep_the_schedule_varying(self, tmp_path):
        """The gate that keeps survey frames from being split by scan range."""
        schedule, separable = _schedule(_copy(tmp_path))
        assert schedule.constant_across_pixels is False
        assert separable is False

    def test_one_window_group_without_survey_frames_is_separable(self, tmp_path):
        schedule, separable = _schedule(_copy(tmp_path, "one_group"))

        assert _got(schedule) == _expected(WINDOWS)
        assert schedule.constant_across_pixels is True
        assert schedule.merges_precursors is True
        assert separable is True

    def test_window_groups_that_take_turns_vary(self, tmp_path):
        """Each group's windows are in half the frames: no global schedule."""
        schedule, separable = _schedule(_copy(tmp_path, "groups_take_turns"))

        assert _got(schedule) == sorted(
            _expected(WINDOWS + SHIFTED_WINDOWS), key=lambda w: w[2]
        )
        assert schedule.constant_across_pixels is False
        assert separable is False

    def test_empty_dia_tables_record_the_level_only(self, tmp_path):
        path = _copy(tmp_path, "one_group")
        with closing(sqlite3.connect(path / "analysis.tdf")) as conn:
            conn.execute("DELETE FROM DiaFrameMsMsWindows")
            conn.commit()
        schedule, separable = _schedule(path)

        assert schedule.ms_level == 2
        assert schedule.windows == ()
        assert separable is False

    @pytest.mark.xfail(strict=True, reason=LATER_STEP)
    def test_survey_frames_do_not_make_the_schedule_vary(self, tmp_path):
        """Constancy judged over the fragment frames alone.

        Every type-9 frame of the mixed shape uses the same window group.
        Only safe once survey frames never reach the precursor split.
        """
        schedule, _ = _schedule(_copy(tmp_path))
        assert schedule.constant_across_pixels is True


class TestTheSplit:
    """Splitting type-9 frames by their group's windows, through a stand-in SDK."""

    def _split(self, path, points=None):
        reader = _open(path)
        reader.sdk = FakeSDK() if points is None else FakeSDK(points)
        reader.handle = 1
        try:
            return list(reader.iter_precursor_spectra())
        finally:
            reader.close()

    def test_one_group_splits_by_its_windows(self, tmp_path):
        # FakeSDK's default points by scan: 70 and 90 fall in window 0
        # [60, 100), 125 in window 1, 160 in window 2; 25, 30, 55, 199 and
        # 220 in none.
        yielded = self._split(_copy(tmp_path, "one_group"))

        assert len(yielded) == 6 * 3
        totals = defaultdict(float)
        for _, window_index, _, intensities in yielded:
            totals[window_index] += float(intensities.sum())
        assert dict(totals) == {0: 6 * (999 + 200), 1: 6 * 300, 2: 6 * 400}

    def test_a_boundary_scan_belongs_to_the_window_it_begins(self, tmp_path):
        """``ScanNumEnd`` is exclusive: scan 100 is window 1's, not window 0's."""
        yielded = self._split(
            _copy(tmp_path, "one_group"), points=[(5, 11, 99), (6, 13, 100)]
        )
        owners = {(w, float(i.sum())) for _, w, _, i in yielded}
        assert owners == {(0, 11.0), (1, 13.0)}

    @pytest.mark.parametrize("shape", ["mixed", "groups_take_turns"])
    def test_a_varying_schedule_is_not_split(self, tmp_path, shape):
        """A direct caller gets a refusal, never survey frames cut into blocks."""
        with pytest.raises(NotImplementedError, match="not the same in every frame"):
            self._split(_copy(tmp_path, shape))


class TestMixedPositions:
    def test_a_position_holding_two_levels_is_warned_about_once(
        self, tmp_path, thyra_logs
    ):
        reader = _open(_copy(tmp_path))
        try:
            with thyra_logs("thyra.readers", logging.WARNING) as records:
                list(reader._iter_frames())
                list(reader._iter_frames())
        finally:
            reader.close()
        mixed = [r.getMessage() for r in records if "MsMsType" in r.getMessage()]
        assert len(mixed) == 1
        assert mixed[0].startswith("1 raster position(s)")
        assert "(0, 9)" in mixed[0]

    def test_only_the_selected_region_is_asked(self, tmp_path, thyra_logs):
        """A spot shared across two regions is two pixels, not one."""
        path = _copy(tmp_path)
        with closing(sqlite3.connect(path / "analysis.tdf")) as conn:
            conn.execute("UPDATE MaldiFrameInfo SET RegionNumber = 1 WHERE Frame = 7")
            conn.commit()
        with thyra_logs("thyra.readers", logging.WARNING) as records:
            reader = BrukerReader(path, metadata_only=True, region=0)
            try:
                list(reader._iter_frames())
            finally:
                reader.close()
        assert not [r for r in records if "MsMsType" in r.getMessage()]

    def test_one_level_per_position_is_not_warned_about(self, tmp_path, thyra_logs):
        reader = _open(_copy(tmp_path, "one_group"))
        try:
            with thyra_logs("thyra.readers", logging.WARNING) as records:
                list(reader._iter_frames())
        finally:
            reader.close()
        assert not [r for r in records if "MsMsType" in r.getMessage()]

    @pytest.mark.xfail(strict=True, reason=LATER_STEP)
    def test_no_pixel_combines_a_survey_and_a_fragment_frame(self, tmp_path):
        """The MSI frame loop never feeds one pixel two ``MsMsType`` values."""
        path = _copy(tmp_path)
        types = _msms_types(path)
        reader = _open(path)
        try:
            kinds = defaultdict(set)
            for frame_id, coords in reader._iter_frames():
                kinds[coords].add(types[frame_id])
        finally:
            reader.close()
        assert {c: k for c, k in kinds.items() if len(k) > 1} == {}
