"""Type-9 (DIA-layout) fragment frames end to end, through the vendor library.

Two shapes from ``tests/fixtures/iprm_pasef_builder.py``:

- ``one_group``: every frame a type-9 fragment frame of one window group.
  The schedule is read from ``DiaFrameMsMsWindows`` (D37) and the frames
  are split into the ``{table}_msms`` sibling by scan number.
- ``mixed``: a survey frame (``MsMsType = 0``) and a fragment frame share a
  raster spot. Today the converter sums both into that spot's pixel row;
  the tests that pin keeping them apart are ``xfail(strict=True)`` until
  that later change lands.

Skips where the library cannot be loaded.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pytest

from tests.fixtures.iprm_pasef_builder import (
    SHARED_SPOT_FRAMES,
    WINDOWS,
    make_iprm_pasef_copy,
)
from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
from thyra.utils.bruker_exceptions import SDKError

EXPECTED = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "fixtures"
    / "synthetic_tims_expected.json"
)
LATER_STEP = (
    "survey and fragment frames on one spot are not kept apart yet; this "
    "pins the behaviour of that later change"
)


def _convert(source, out):
    from thyra.convert import convert_msi

    # convert_msi turns a reader that cannot open into ``False``, so ask the
    # library directly first: a missing library is a skip, not a failure.
    try:
        BrukerReader(source).close()
    except (SDKError, OSError) as exc:
        pytest.skip(f"Bruker library not loadable on this platform: {exc}")
    ok = convert_msi(
        str(source),
        str(out),
        dataset_id="iprm",
        pixel_size_um=20.0,
        reader_options={"tdf_spectrum": "scan_sum"},
        include_optical=False,
    )
    assert ok
    return out


def _tables(out):
    spatialdata = pytest.importorskip("spatialdata")
    from thyra.utils.windows_paths import prepare_zarr_read_path

    return spatialdata.read_zarr(prepare_zarr_read_path(out)).tables


def _rows(table) -> np.ndarray:
    X = table.X
    return np.asarray(X.toarray() if hasattr(X, "toarray") else X)


class TestOneWindowGroup:
    def test_the_frames_are_split_by_their_groups_windows(self, tmp_path):
        """Each pixel's precursor blocks hold exactly the points inside a window.

        The fixture's expected file lists every ``(index, scan, intensity)``
        it wrote, so the current inside ``[ScanNumBegin, ScanNumEnd)`` of
        each window is known per frame; the split must reproduce it.
        """
        make_iprm_pasef_copy(tmp_path / "dia.d", "one_group")
        tables = _tables(_convert(tmp_path / "dia.d", tmp_path / "dia.zarr"))

        assert set(tables) == {"iprm_z0", "iprm_z0_msms"}
        summed, msms = tables["iprm_z0"], tables["iprm_z0_msms"]
        assert list(summed.obs.index) == list(msms.obs.index)
        assert sorted(msms.var["precursor_mz"].unique()) == [w[2] for w in WINDOWS]
        assert summed.uns["msms_schedule"]["resolved_table"] == "iprm_z0_msms"

        frames = json.loads(EXPECTED.read_text(encoding="utf-8"))["frames"]
        x0 = min(f["x"] for f in frames)
        y0 = min(f["y"] for f in frames)
        rows = _rows(msms)
        for f in frames:
            at = np.flatnonzero(
                (msms.obs["x"] == f["x"] - x0) & (msms.obs["y"] == f["y"] - y0)
            )
            assert at.size == 1
            for b, e, mz, _, _ in WINDOWS:
                inside = sum(i for _, s, i in f["pairs"] if b <= s < e)
                block = (msms.var["precursor_mz"] == mz).to_numpy()
                assert rows[at[0], block].sum() == pytest.approx(inside), (f, mz)


class TestSurveyAndFragmentOnOneSpot:
    def test_the_mixed_spot_is_warned_about(self, tmp_path, thyra_logs):
        make_iprm_pasef_copy(tmp_path / "iprm.d")
        with thyra_logs("thyra", logging.WARNING) as records:
            tables = _tables(_convert(tmp_path / "iprm.d", tmp_path / "iprm.zarr"))
        messages = [r.getMessage() for r in records]

        assert any(
            "hold frames of more than one MsMsType (0, 9)" in m for m in messages
        )
        # The schedule is recorded, but no table is split from it.
        assert set(tables) == {"iprm_z0"}
        assert tables["iprm_z0"].uns["msms_schedule"]["n_windows"] == 3

    @pytest.mark.xfail(strict=True, reason=LATER_STEP)
    def test_the_spot_row_holds_the_survey_frame_only(self, tmp_path):
        """The MSI row at the shared spot is the survey frame's spectrum.

        ``AccumulationTime`` is 100 ms in the fixture, so the scan sum of a
        frame equals its ``SummedIntensities`` exactly. Today the row holds
        the survey frame plus the fragment frame.
        """
        facts = make_iprm_pasef_copy(tmp_path / "iprm.d")
        table = _tables(_convert(tmp_path / "iprm.d", tmp_path / "iprm.zarr"))[
            "iprm_z0"
        ]

        survey, fragment = SHARED_SPOT_FRAMES
        summed = facts["summed_intensities"]
        at_spot = (table.obs["x"] == 0) & (table.obs["y"] == 0)
        assert int(at_spot.sum()) == 1
        row_tic = float(_rows(table)[np.flatnonzero(at_spot)[0]].sum())

        assert row_tic != pytest.approx(summed[survey] + summed[fragment])
        assert row_tic == pytest.approx(summed[survey])

    @pytest.mark.xfail(strict=True, reason=LATER_STEP)
    def test_no_position_is_summed_across_ms_levels(self, tmp_path, thyra_logs):
        make_iprm_pasef_copy(tmp_path / "iprm.d")
        with thyra_logs("thyra.converters", logging.WARNING) as records:
            _convert(tmp_path / "iprm.d", tmp_path / "iprm.zarr")
        summed = [
            r.getMessage()
            for r in records
            if "carry more than one spectrum" in r.getMessage()
        ]
        assert summed == []
