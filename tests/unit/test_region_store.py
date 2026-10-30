"""What a store records about its regions (issue #496).

A ``region=`` store describes the region it holds, not the whole run (G1,
G2). A position measured in two regions takes the region of its first
frame, and the reader says so (G4). ``region=`` takes any integer type
(G5). A region with no row has no mean (G6). The input is the synthetic
six-frame ``.d``, rewritten per test: frames 1-3 sit at x 0-2 on the first
row, frames 4-6 at x 0-2 on the second.
"""

from __future__ import annotations

import json
import logging
import shutil
import sqlite3
from pathlib import Path
from typing import Iterable

import numpy as np
import pytest
import spatialdata

from thyra.convert import convert_msi
from thyra.core.region_map import build_region_map
from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

_TDF = Path(__file__).resolve().parents[1] / "data" / "fixtures" / "synthetic_tims.d"

# (IsolationMz, CollisionEnergy, ScanNumBegin, ScanNumEnd), as the PASEF
# schedule fixtures elsewhere write them.
_SCHEDULE = [
    (936.578, 49.6, 20, 60),
    (353.32, 35.2, 80, 130),
    (313.275, 34.3, 150, 200),
]
_OTHER_SCHEDULE = [
    (500.0, 30.0, 20, 60),
    (600.0, 30.0, 80, 130),
    (700.0, 30.0, 150, 200),
]


def _run(tmp_path: Path, *statements: str, name: str = "run.d") -> Path:
    """Copy the synthetic ``.d`` and apply ``statements`` to its database."""
    run = tmp_path / name
    shutil.copytree(_TDF, run)
    with sqlite3.connect(run / "analysis.tdf") as db:
        for statement in statements:
            db.execute(statement)
    return run


def _two_regions(tmp_path: Path, *statements: str, name: str = "run.d") -> Path:
    """Frames 1-3 are region 0, frames 4-6 region 1."""
    return _run(
        tmp_path,
        "UPDATE MaldiFrameInfo SET RegionNumber = (Frame > 3)",
        *statements,
        name=name,
    )


def _pasef(run: Path, frames: Iterable[int], schedule=_SCHEDULE) -> None:
    """Mark ``frames`` as PASEF frames isolating ``schedule``."""
    frames = list(frames)
    with sqlite3.connect(run / "analysis.tdf") as db:
        db.execute(
            "UPDATE Frames SET MsMsType = 8 WHERE Id IN (%s)"
            % ",".join(str(f) for f in frames)
        )
        db.execute(
            "CREATE TABLE IF NOT EXISTS PasefFrameMsMsInfo (Frame INTEGER NOT "
            "NULL, ScanNumBegin INTEGER NOT NULL, ScanNumEnd INTEGER NOT NULL, "
            "IsolationMz REAL NOT NULL, IsolationWidth REAL NOT NULL, "
            "CollisionEnergy REAL NOT NULL, Precursor INTEGER)"
        )
        db.executemany(
            "INSERT INTO PasefFrameMsMsInfo VALUES (?, ?, ?, ?, 1.0, ?, NULL)",
            [(f, lo, hi, mz, ce) for f in frames for mz, ce, lo, hi in schedule],
        )


def _reader(run: Path, region=None) -> BrukerReader:
    return BrukerReader(run, region=region, metadata_only=True)


def _convert(run: Path, out: Path, **kwargs) -> dict:
    assert convert_msi(
        str(run),
        str(out),
        pixel_size_um=20.0,
        include_optical=False,
        reader_options={"tdf_spectrum": "scan_sum"},
        **kwargs,
    )
    return dict(spatialdata.read_zarr(out).tables)


class TestARegionStoreListsOnlyItsRegion:
    """G1: ``uns["regions"]`` of a ``region=`` store."""

    def test_the_reader_summary_holds_the_selected_region(self, tmp_path):
        run = _two_regions(
            tmp_path,
            "UPDATE MaldiFrameInfo SET XIndexPos = XIndexPos + 10, "
            "YIndexPos = YIndexPos + 10 WHERE Frame > 3",
        )
        reader = _reader(run, region=0)
        try:
            info = reader.get_region_info()
        finally:
            reader.close()

        assert [entry["region_number"] for entry in info] == [0]
        assert info[0]["n_spectra"] == 3
        assert info[0]["bounds"] == (0, 0, 2, 0)

    def test_every_region_is_listed_without_a_selection(self, tmp_path):
        reader = _reader(_two_regions(tmp_path))
        try:
            numbers = sorted(e["region_number"] for e in reader.get_region_info())
        finally:
            reader.close()

        assert numbers == [0, 1]

    def test_the_store_lists_one_region_inside_its_raster(self, tmp_path):
        run = _two_regions(
            tmp_path,
            "UPDATE MaldiFrameInfo SET XIndexPos = XIndexPos + 10, "
            "YIndexPos = YIndexPos + 10 WHERE Frame > 3",
        )
        table = next(iter(_convert(run, tmp_path / "out.zarr", region=0).values()))

        regions = json.loads(table.uns["regions"])
        assert [entry["region_number"] for entry in regions] == [0]
        assert regions[0]["n_spectra"] == len(table) == 3
        x_min, y_min, x_max, y_max = regions[0]["bounds"]
        assert (x_min, y_min) == (table.obs["x"].min(), table.obs["y"].min())
        assert (x_max, y_max) == (table.obs["x"].max(), table.obs["y"].max())


class TestARegionScheduleIsJudgedOnItsFrames:
    """G2: frame types and isolation windows under ``region=``."""

    def test_a_fragmenting_region_beside_a_survey_one_is_constant(self, tmp_path):
        run = _two_regions(tmp_path)
        _pasef(run, range(1, 4))
        reader = _reader(run, region=0)
        try:
            schedule = reader.get_fragmentation()
        finally:
            reader.close()

        assert schedule.ms_level == 2
        assert schedule.constant_across_pixels
        assert len(schedule.windows) == 3

    def test_a_survey_region_beside_a_fragmenting_one_is_ms1(self, tmp_path):
        run = _two_regions(tmp_path)
        _pasef(run, range(4, 7))
        reader = _reader(run, region=0)
        try:
            schedule = reader.get_fragmentation()
        finally:
            reader.close()

        assert schedule.ms_level == 1
        assert not schedule.windows

    def test_only_the_regions_own_windows_are_listed(self, tmp_path):
        run = _two_regions(tmp_path)
        _pasef(run, range(1, 4))
        _pasef(run, range(4, 7), schedule=_OTHER_SCHEDULE)
        reader = _reader(run, region=0)
        try:
            schedule = reader.get_fragmentation()
        finally:
            reader.close()

        assert sorted(round(w.target, 3) for w in schedule.windows) == sorted(
            mz for mz, _, _, _ in _SCHEDULE
        )
        assert schedule.constant_across_pixels

    def test_the_whole_run_still_varies(self, tmp_path):
        run = _two_regions(tmp_path)
        _pasef(run, range(1, 4))
        reader = _reader(run)
        try:
            schedule = reader.get_fragmentation()
        finally:
            reader.close()

        assert not schedule.constant_across_pixels

    def test_the_region_store_gets_its_msms_table(self, tmp_path):
        run = _two_regions(tmp_path)
        _pasef(run, range(1, 4))

        tables = _convert(run, tmp_path / "out.zarr", region=0)

        assert sorted(tables) == ["msi_dataset_z0", "msi_dataset_z0_msms"]


class TestAPositionInTwoRegions:
    """G4: a position two regions measured."""

    def test_the_first_scan_names_the_region(self):
        region_map = build_region_map([(0, 0, 0), (1, 0, 0), (0, 0, 1), (2, 0, 1)])

        assert region_map == {(0, 0): 0, (1, 0): 0, (2, 0): 1}

    def test_the_warning_names_the_regions(self, thyra_logs):
        with thyra_logs("thyra.core.region_map", logging.WARNING) as logs:
            build_region_map([(0, 0, 4), (5, 3, 4), (5, 3, 7)])

        assert len(logs) == 1
        message = logs[0].getMessage()
        assert "1 position(s)" in message
        assert "regions 4, 7" in message
        assert "x=5, y=3" in message

    def test_one_region_per_position_is_quiet(self, thyra_logs):
        with thyra_logs("thyra.core.region_map", logging.WARNING) as logs:
            build_region_map([(0, 0, 0), (0, 0, 0), (1, 0, 1)])

        assert not logs

    def test_the_timstof_map_follows_frame_order(self, tmp_path, thyra_logs):
        # Frame 6 (region 1) is measured where frame 1 (region 0) was.
        run = _two_regions(
            tmp_path,
            "UPDATE MaldiFrameInfo SET XIndexPos = "
            "(SELECT XIndexPos FROM MaldiFrameInfo WHERE Frame = 1), "
            "YIndexPos = (SELECT YIndexPos FROM MaldiFrameInfo WHERE Frame = 1) "
            "WHERE Frame = 6",
        )
        reader = _reader(run)
        try:
            with thyra_logs("thyra.core.region_map", logging.WARNING) as logs:
                region_map = reader.get_region_map()
        finally:
            reader.close()

        assert region_map[(0, 0)] == 0
        assert len(logs) == 1
        assert "regions 0, 1" in logs[0].getMessage()

    def test_the_row_takes_the_region_of_its_acquisition_order(self, tmp_path):
        run = _two_regions(
            tmp_path,
            "UPDATE MaldiFrameInfo SET XIndexPos = "
            "(SELECT XIndexPos FROM MaldiFrameInfo WHERE Frame = 1), "
            "YIndexPos = (SELECT YIndexPos FROM MaldiFrameInfo WHERE Frame = 1) "
            "WHERE Frame = 6",
        )
        table = next(iter(_convert(run, tmp_path / "out.zarr").values()))

        row = table.obs[(table.obs["x"] == 0) & (table.obs["y"] == 0)].iloc[0]
        assert row["acquisition_order"] == 1
        assert row["region_number"] == 0

    def test_an_empty_first_frame_does_not_name_the_region(self, tmp_path):
        """Frame 1 holds nothing, so frame 6 gives the row its order and region."""
        run = _two_regions(
            tmp_path,
            "UPDATE MaldiFrameInfo SET XIndexPos = "
            "(SELECT XIndexPos FROM MaldiFrameInfo WHERE Frame = 1), "
            "YIndexPos = (SELECT YIndexPos FROM MaldiFrameInfo WHERE Frame = 1) "
            "WHERE Frame = 6",
            "UPDATE Frames SET NumScans = 0, NumPeaks = 0 WHERE Id = 1",
        )
        table = next(iter(_convert(run, tmp_path / "out.zarr").values()))

        row = table.obs[(table.obs["x"] == 0) & (table.obs["y"] == 0)].iloc[0]
        assert row["acquisition_order"] == 6
        assert row["region_number"] == 1


class TestTheRegionOptionTakesAnyInteger:
    """G5: a numpy integer, as read back from ``obs["region_number"]``."""

    @pytest.mark.parametrize("kind", [np.int32, np.int64, np.uint8])
    def test_a_numpy_integer_selects_the_region(self, tmp_path, kind):
        reader = _reader(_two_regions(tmp_path), region=kind(1))
        try:
            assert reader._selected_region == 1
            assert type(reader._requested_region) is int
        finally:
            reader.close()

    def test_a_bool_is_still_refused(self, tmp_path):
        from thyra.errors import ConversionRefused

        with pytest.raises(ConversionRefused, match="got bool"):
            _reader(_two_regions(tmp_path), region=True)

    def test_a_region_read_back_from_a_store_converts(self, tmp_path):
        run = _two_regions(tmp_path)
        full = next(iter(_convert(run, tmp_path / "full.zarr").values()))
        region = np.unique(full.obs["region_number"])[0]
        assert not isinstance(region, int)

        table = next(iter(_convert(run, tmp_path / "one.zarr", region=region).values()))

        assert set(table.obs["region_number"]) == {int(region)}


class TestARegionWithNoRowHasNoMean:
    """G6: ``uns["average_spectrum_per_region"]``."""

    def test_the_empty_region_has_no_entry(self, tmp_path):
        # Frames 4 and 6 are region 1; frame 5 is alone in region 2 and
        # holds no signal.
        run = _run(
            tmp_path,
            "UPDATE MaldiFrameInfo SET RegionNumber = 1 WHERE Frame IN (4, 6)",
            "UPDATE MaldiFrameInfo SET RegionNumber = 2 WHERE Frame = 5",
            "UPDATE Frames SET NumScans = 0, NumPeaks = 0 WHERE Id = 5",
        )
        table = next(iter(_convert(run, tmp_path / "out.zarr").values()))

        assert set(table.obs["region_number"]) == {0, 1}
        assert sorted(table.uns["average_spectrum_per_region"]) == ["0", "1"]

    def test_two_regions_keep_the_key_when_one_holds_rows(self, tmp_path):
        """The store's shape follows its regions, not which hold rows (D57)."""
        run = _run(
            tmp_path,
            "UPDATE MaldiFrameInfo SET RegionNumber = 2 WHERE Frame = 5",
            "UPDATE Frames SET NumScans = 0, NumPeaks = 0 WHERE Id = 5",
        )
        table = next(iter(_convert(run, tmp_path / "out.zarr").values()))

        assert set(table.obs["region_number"]) == {0}
        assert sorted(table.uns["average_spectrum_per_region"]) == ["0"]


class TestTheIntegerColumnsD26Names:
    """G7: D26 states the types of ``x``, ``y`` and ``region_number``."""

    def test_a_2d_store_writes_them_as_int32_and_int64(self, tmp_path):
        """``region_number`` is int64: a source may number past 2**31 (#527)."""
        table = next(iter(_convert(_TDF, tmp_path / "out.zarr").values()))

        dtypes = {c: str(table.obs[c].dtype) for c in ("x", "y", "region_number")}
        assert dtypes == {"x": "int32", "y": "int32", "region_number": "int64"}
        assert str(table.obs["acquisition_order"].dtype) == "int64"


class TestASingleRegionKeepsItsNumber:
    """G3: one region is stored with the number and name the source states."""

    def test_a_timstof_region_numbered_3(self, tmp_path):
        run = _run(tmp_path, "UPDATE MaldiFrameInfo SET RegionNumber = 3")
        table = next(iter(_convert(run, tmp_path / "out.zarr").values()))

        assert set(table.obs["region_number"]) == {3}
        (entry,) = json.loads(table.uns["regions"])
        assert entry["region_number"] == 3
        assert entry["n_spectra"] == 6
        assert "bounds" in entry
        assert "average_spectrum_per_region" not in table.uns

    def test_two_regions_keep_their_numbers_and_means(self, tmp_path):
        run = _run(
            tmp_path,
            "UPDATE MaldiFrameInfo SET RegionNumber = 0 WHERE Frame <= 3",
            "UPDATE MaldiFrameInfo SET RegionNumber = 3 WHERE Frame > 3",
        )
        table = next(iter(_convert(run, tmp_path / "out.zarr").values()))

        assert set(table.obs["region_number"]) == {0, 3}
        assert sorted(table.uns["average_spectrum_per_region"]) == ["0", "3"]

    def test_an_mzpeak_region_numbered_3_and_named(self, tmp_path):
        from tests.fixtures.mzpeak_builder import (
            Spectrum,
            bruker_region,
            build_mzpeak,
            region_parameter,
        )

        spectra = [
            Spectrum(x, y, [100.0, 101.0], [1.0, 2.0])
            for y in (1, 2)
            for x in (1, 2, 3)
        ]
        archive = build_mzpeak(
            tmp_path / "run.mzpeak",
            spectra,
            scan_columns={"parameters": [region_parameter(3)] * 6},
            file_metadata={
                "bruker_maldi": {
                    "regions": [bruker_region(3, 6, (1, 3), (1, 2), name="tumour")]
                }
            },
        )
        assert convert_msi(
            str(archive),
            str(tmp_path / "out.zarr"),
            pixel_size_um=20.0,
            include_optical=False,
        )
        table = next(iter(spatialdata.read_zarr(tmp_path / "out.zarr").tables.values()))

        assert set(table.obs["region_number"]) == {3}
        (entry,) = json.loads(table.uns["regions"])
        assert (entry["region_number"], entry["name"]) == (3, "tumour")
