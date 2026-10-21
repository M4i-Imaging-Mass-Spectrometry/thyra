"""A conversion that will not fit is refused in one line, through the converter (issue #528).

The 4.3.3 stress test found every one of these ending in a traceback:

* an axis past the count ceiling (F59), built first and then refused by a
  bare ``MemoryError``;
* a source with many regions (F148) or many planes on the 2D route (F116),
  whose extra per-bin cost the guard did not count;
* a matrix bigger than the free disk (F149), which finished the pre-scan
  and then failed on ``[Errno 28]``.

The per-bin costs and free disk are patched rather than reached: nothing
here allocates more than a few kilobytes.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest
from click.testing import CliRunner

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.converters.spatialdata import csc_assembly
from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)
from thyra.resampling import memory_guard

#: The logger the base workflow reports a refusal on.
_REFUSALS = "thyra.core.base_converter"


def _config(n_z: int = 1) -> MockMSIConfig:
    return MockMSIConfig(
        n_x=3, n_y=2, n_z=n_z, n_mz_bins=64, peaks_per_spectrum=(3, 5), seed=5
    )


class _TwoRegions(MockMSIReader):
    """The left column is region 1, the rest region 2."""

    def get_region_map(self):
        return {
            (x, y): (1 if x == 0 else 2)
            for y in range(self.config.n_y)
            for x in range(self.config.n_x)
        }


def _converter(reader, output: Path) -> StreamingSpatialDataConverter:
    return StreamingSpatialDataConverter(
        reader=reader,
        output_path=output,
        dataset_id="m",
        pixel_size_um=10.0,
        include_optical=False,
    )


def _refusals(records) -> list:
    return [r.getMessage() for r in records if r.levelname == "ERROR"]


def _leftovers(folder: Path) -> list:
    return sorted(p.name for p in folder.iterdir() if p.name.startswith(".thyra_"))


@pytest.fixture
def plenty_free(monkeypatch):
    """8 GB free, whatever this machine has, so only the patched costs decide."""
    monkeypatch.setattr(memory_guard, "available_memory_gb", lambda: 8.0)


class TestRegionsAddToTheProjection:
    def test_two_regions_are_counted_before_their_arrays_exist(
        self, tmp_path, monkeypatch, plenty_free, thyra_logs
    ):
        monkeypatch.setattr(memory_guard, "REGION_BYTES_PER_BIN", 10**9)
        output = tmp_path / "regions.zarr"

        with thyra_logs(_REFUSALS, logging.ERROR) as records:
            assert _converter(_TwoRegions(_config()), output).convert() is False

        messages = _refusals(records)
        assert any("for each of 2 regions" in m for m in messages), messages
        assert not any("Traceback" in m for m in messages)
        assert not output.exists()

    def test_one_region_builds_no_region_arrays_and_is_not_charged(
        self, tmp_path, monkeypatch, plenty_free
    ):
        monkeypatch.setattr(memory_guard, "REGION_BYTES_PER_BIN", 10**9)
        output = tmp_path / "one.zarr"

        assert _converter(MockMSIReader(_config()), output).convert() is True


class TestPlaneTablesAddToTheProjection:
    def test_the_tables_the_pre_scan_finds_are_counted(
        self, tmp_path, monkeypatch, plenty_free, thyra_logs
    ):
        """The 2D route writes one table per plane, known only after pass 1."""
        monkeypatch.setattr(memory_guard, "TABLE_BYTES_PER_BIN", 10**9)
        output = tmp_path / "planes.zarr"

        with thyra_logs(_REFUSALS, logging.ERROR) as records:
            assert _converter(MockMSIReader(_config(n_z=3)), output).convert() is False

        messages = _refusals(records)
        assert any("for each of 2 more tables" in m for m in messages), messages
        assert not output.exists()
        assert _leftovers(tmp_path) == []

    def test_one_plane_is_one_table(self, tmp_path, monkeypatch, plenty_free):
        monkeypatch.setattr(memory_guard, "TABLE_BYTES_PER_BIN", 10**9)

        converter = _converter(MockMSIReader(_config()), tmp_path / "one.zarr")
        assert converter.convert() is True

    def test_free_memory_is_read_once_before_the_axis(self, tmp_path, monkeypatch):
        """Later checks compare with that figure, not with what is left by then."""
        reads = []

        def free() -> float:
            reads.append(1)
            return 8.0

        monkeypatch.setattr(memory_guard, "available_memory_gb", free)
        converter = _converter(MockMSIReader(_config(n_z=3)), tmp_path / "x.zarr")

        assert converter.convert() is True
        assert len(reads) == 1


class TestTheDiskIsChecked:
    def test_a_matrix_bigger_than_the_drive_is_refused_before_its_files(
        self, tmp_path, monkeypatch, thyra_logs
    ):
        monkeypatch.setattr(csc_assembly, "free_disk_bytes", lambda folder: 100)
        output = tmp_path / "full.zarr"

        with thyra_logs(_REFUSALS, logging.ERROR) as records:
            assert _converter(MockMSIReader(_config()), output).convert() is False

        messages = _refusals(records)
        assert any(
            "matrix files" in m and "GB free" in m and str(tmp_path) in m
            for m in messages
        ), messages
        assert not output.exists()
        assert _leftovers(tmp_path) == []

    def test_past_half_of_the_drive_it_warns_and_converts(
        self, tmp_path, monkeypatch, thyra_logs
    ):
        monkeypatch.setattr(csc_assembly, "DISK_WARN_FRACTION", 0.0)
        output = tmp_path / "tight.zarr"

        with thyra_logs("thyra", logging.WARNING) as records:
            assert _converter(MockMSIReader(_config()), output).convert() is True

        assert any("the store is written there too" in r.getMessage() for r in records)
        assert output.is_dir()


class TestTheCommandLine:
    def test_an_axis_past_the_ceiling_is_refused_unbuilt_and_without_traceback(
        self, create_minimal_imzml, temp_dir, monkeypatch, thyra_logs
    ):
        """F59: refused at once, whatever memory is free, before the axis exists."""
        from thyra.__main__ import main

        monkeypatch.setattr(memory_guard, "available_memory_gb", lambda: 1e6)
        imzml_path, _, _, _ = create_minimal_imzml
        output = temp_dir / "out.zarr"

        with thyra_logs(_REFUSALS, logging.ERROR) as records:
            result = CliRunner().invoke(
                main,
                [
                    str(imzml_path),
                    str(output),
                    "--pixel-size",
                    "1.0",
                    "--resample-bins",
                    str(memory_guard.MAX_COUNT_KEYS + 1),
                ],
            )

        assert result.exit_code == 1, result.output
        assert "Traceback" not in result.output
        messages = _refusals(records)
        assert any(
            "268,435,457 bins" in m and "--resample-bins" in m for m in messages
        ), messages
        assert not output.exists()
