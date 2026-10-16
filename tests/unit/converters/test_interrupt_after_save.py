"""A Ctrl-C after the store is complete keeps the store (issue #500).

An interrupt is a failed conversion (#245): the partial store is moved
to ``.failed`` and the scratch goes. But ``_save_output`` still has work
to do after the last write -- a log line, dropping the SpatialData, the
scratch release -- and an interrupt there turned a store that held every
row and opened cleanly into ``<name>.zarr.failed``, with a log line that
said nothing usable was written.

The boundary is the metadata consolidation. Before it, the store is not
finished: ``SpatialData.write`` has already consolidated once, but the
optical pixels have not streamed in, so a store on disk that *looks*
consolidated is not proof of anything. After it, the store is complete.
Both sides are pinned here.
"""

from __future__ import annotations

import gc
import logging
from pathlib import Path

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)

_KEPT = "is complete and was kept"
_LOST = "Nothing usable was written"


def _converter(output: Path) -> StreamingSpatialDataConverter:
    config = MockMSIConfig(
        n_x=3, n_y=2, n_z=1, n_mz_bins=200, peaks_per_spectrum=(5, 9), seed=3
    )
    return StreamingSpatialDataConverter(
        reader=MockMSIReader(config),
        output_path=output,
        dataset_id="m",
        pixel_size_um=10.0,
        include_optical=False,
    )


def _interrupt_once(original):
    """Raise ``KeyboardInterrupt`` on the first call, then behave."""
    calls = []

    def replacement(*args, **kwargs):
        calls.append(True)
        if len(calls) == 1:
            raise KeyboardInterrupt
        return original(*args, **kwargs)

    return replacement


def _interrupted_before_the_last_write(output: Path):
    """Between ``SpatialData.write`` and the final consolidation."""
    converter = _converter(output)
    converter.optical.stream_pending_pixels = _interrupt_once(
        converter.optical.stream_pending_pixels
    )
    return converter


def _interrupted_after_the_last_write(output: Path):
    """In the scratch release: the window issue #500 measured."""
    converter = _converter(output)
    converter.siblings.release_scratch = _interrupt_once(
        converter.siblings.release_scratch
    )
    return converter


def _messages(records) -> str:
    return "\n".join(record.getMessage() for record in records)


def _scratch(parent: Path) -> list[str]:
    return sorted(p.name for p in parent.iterdir() if p.name.startswith(".thyra_"))


class TestBeforeTheStoreIsComplete:
    def test_it_is_a_failure(self, tmp_path):
        converter = _interrupted_before_the_last_write(tmp_path / "b.zarr")

        assert converter.convert() is False
        assert converter.store_complete is False

    def test_the_log_says_nothing_usable_was_written(self, tmp_path, thyra_logs):
        converter = _interrupted_before_the_last_write(tmp_path / "b.zarr")

        with thyra_logs("thyra", logging.WARNING) as records:
            converter.convert()

        assert _LOST in _messages(records)
        assert _KEPT not in _messages(records)


class TestAfterTheStoreIsComplete:
    def test_it_is_a_success(self, tmp_path):
        converter = _interrupted_after_the_last_write(tmp_path / "a.zarr")

        assert converter.convert() is True
        assert converter.store_complete is True

    def test_the_store_opens_and_holds_every_row(self, tmp_path):
        import spatialdata

        out = tmp_path / "a.zarr"
        _interrupted_after_the_last_write(out).convert()

        sdata = spatialdata.read_zarr(str(out))
        assert sdata.tables["m_z0"].n_obs == 3 * 2

    def test_the_log_says_it_was_kept(self, tmp_path, thyra_logs):
        converter = _interrupted_after_the_last_write(tmp_path / "a.zarr")

        with thyra_logs("thyra", logging.WARNING) as records:
            converter.convert()

        assert _KEPT in _messages(records)
        assert _LOST not in _messages(records)

    def test_the_scratch_still_goes(self, tmp_path):
        """The release the interrupt cut short runs again in ``convert``'s
        ``finally``, so keeping the store does not leave the scratch."""
        converter = _interrupted_after_the_last_write(tmp_path / "a.zarr")
        converter.convert()

        del converter
        gc.collect()
        assert _scratch(tmp_path) == []


class TestTheMarkFollowsTheSave:
    def test_a_finished_conversion_is_marked(self, tmp_path):
        converter = _converter(tmp_path / "p.zarr")

        assert converter.convert() is True
        assert converter.store_complete is True

    def test_a_failed_write_is_not(self, tmp_path, monkeypatch):
        import spatialdata

        def boom(self, *_args, **_kwargs):
            raise OSError("simulated write failure")

        monkeypatch.setattr(spatialdata.SpatialData, "write", boom)
        converter = _converter(tmp_path / "p.zarr")

        assert converter.convert() is False
        assert converter.store_complete is False
