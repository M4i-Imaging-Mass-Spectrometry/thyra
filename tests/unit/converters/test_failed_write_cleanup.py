"""What a failed or interrupted store write leaves behind, and says (#525).

Three gaps the scratch release (#245, #249, #500) still had:

* A second Ctrl-C landing in the cleanup of a first cut the release loop
  short, and the ``finally`` it ran in was the last call there was, so a
  ``.thyra_summed_*`` folder stayed beside the output for good. The log
  then said the run stopped "before it started reading".
* When one chunk write failed while others still ran, those still held
  views of the scratch memmaps; Windows will not delete a mapped file,
  the single removal attempt failed, and the folder stayed.
* A bare ``MemoryError`` -- what a failed allocation raises -- was logged
  as ``Error saving SpatialData: `` with nothing after the colon.
"""

from __future__ import annotations

import asyncio
import logging
import shutil
import time
from pathlib import Path

import pytest

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra import convert as convert_module
from thyra.converters.spatialdata import csc_assembly
from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)


def _converter(output: Path) -> StreamingSpatialDataConverter:
    reader = MockMSIReader(
        MockMSIConfig(
            n_x=3, n_y=2, n_z=1, n_mz_bins=200, peaks_per_spectrum=(5, 9), seed=3
        )
    )
    return StreamingSpatialDataConverter(
        reader=reader,
        output_path=output,
        dataset_id="m",
        pixel_size_um=10.0,
        include_optical=False,
    )


def _left_beside(folder: Path) -> list[str]:
    return sorted(p.name for p in folder.iterdir())


@pytest.fixture
def imzml(create_minimal_imzml) -> Path:
    path, _, _, _ = create_minimal_imzml
    return path


class TestASecondInterruptInTheCleanup:
    def test_every_scratch_folder_goes_and_the_message_fits(
        self, imzml, tmp_path, monkeypatch, thyra_logs
    ):
        from thyra.converters.spatialdata.base_spatialdata_converter import (
            BaseSpatialDataConverter,
        )

        def interrupted_save(self, state):
            raise KeyboardInterrupt  # the first Ctrl-C, before the write

        real_remove = csc_assembly.remove_scratch
        fired = []

        def remove(path):  # the second, as the cleanup reaches its first folder
            if not fired:
                fired.append(path)
                raise KeyboardInterrupt
            return real_remove(path)

        monkeypatch.setattr(BaseSpatialDataConverter, "_save_output", interrupted_save)
        monkeypatch.setattr(csc_assembly, "remove_scratch", remove)
        out_dir = tmp_path / "out"
        out_dir.mkdir()

        with thyra_logs("thyra", logging.ERROR) as records:
            ok = convert_module.convert_msi(
                str(imzml), str(out_dir / "store.zarr"), pixel_size_um=25.0
            )

        assert ok is False
        assert fired
        assert _left_beside(out_dir) == []
        messages = [r.getMessage() for r in records]
        assert not any("before it started reading" in m for m in messages)


class TestAFailedWrite:
    def test_a_bare_memory_error_names_its_type(
        self, tmp_path, monkeypatch, thyra_logs
    ):
        import spatialdata

        def oom(self, *args, **kwargs):
            raise MemoryError()

        monkeypatch.setattr(spatialdata.SpatialData, "write", oom)

        with thyra_logs("thyra", logging.ERROR) as records:
            assert _converter(tmp_path / "out.zarr").convert() is False

        assert "Error saving SpatialData: MemoryError" in [
            r.getMessage() for r in records
        ]

    def test_the_partial_store_is_set_aside_by_the_converter_itself(
        self, tmp_path, monkeypatch
    ):
        """A caller that builds a converter gets the ``.failed`` move too."""
        import zarr

        def fail_after_writing(self, *args, **kwargs):
            raise OSError(28, "No space left on device")

        monkeypatch.setattr(zarr, "consolidate_metadata", fail_after_writing)

        assert _converter(tmp_path / "out.zarr").convert() is False

        assert _left_beside(tmp_path) == ["out.zarr.failed"]

    def test_the_chunk_writes_still_running_are_waited_for(self):
        from zarr.core.sync import _get_loop

        from thyra.converters.spatialdata.base_spatialdata_converter import (
            _wait_for_pending_store_writes,
        )

        loop = _get_loop()
        finished = []

        async def a_chunk_write_still_running():
            await asyncio.sleep(0.3)
            finished.append(True)

        asyncio.run_coroutine_threadsafe(a_chunk_write_still_running(), loop)
        _wait_for_pending_store_writes()

        assert finished == [True]


class TestTheScratchRemoval:
    def test_a_folder_held_for_a_moment_is_still_removed(
        self, tmp_path, monkeypatch, thyra_logs
    ):
        scratch = tmp_path / ".thyra_summed_x"
        scratch.mkdir()
        (scratch / "data.bin").write_bytes(b"0" * 16)
        real_rmtree = shutil.rmtree
        attempts = []

        def held_twice(path, ignore_errors=False):
            attempts.append(time.monotonic())
            if len(attempts) > 2:
                real_rmtree(path, ignore_errors=ignore_errors)

        monkeypatch.setattr(csc_assembly.shutil, "rmtree", held_twice)

        with thyra_logs("thyra", logging.WARNING) as records:
            csc_assembly.remove_scratch(scratch)

        assert not scratch.exists()
        assert len(attempts) == 3
        assert not records

    def test_a_folder_that_survives_stays_listed_for_the_next_call(
        self, tmp_path, monkeypatch
    ):
        from thyra.converters.spatialdata.sibling_tables import SiblingTables

        siblings = SiblingTables.__new__(SiblingTables)
        scratch = tmp_path / ".thyra_summed_y"
        scratch.mkdir()
        siblings._table_scratch = [(None, scratch)]

        monkeypatch.setattr(csc_assembly, "remove_scratch", lambda path: None)
        siblings.release_scratch()
        assert siblings._table_scratch == [(None, scratch)]

        monkeypatch.undo()
        siblings.release_scratch()
        assert siblings._table_scratch == []
        assert not scratch.exists()
