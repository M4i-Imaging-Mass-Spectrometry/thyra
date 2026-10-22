"""A folder that holds several runs is refused, never converted by guess (#525).

A run is a ``.d`` folder Thyra reads, timsTOF or solariX. ``convert_msi``
on a folder of two runs converted the first in name order and returned
True with one INFO line; ``preview_msi`` and ``thyra metadata`` did the
same, and the CLI with nothing on standard input ended in a bare
``Aborted!``. A folder holding one run is read as that run everywhere,
so ``--interactive-calibration`` shows its display and a folder holding
one solariX run converts, where it was refused.
"""

from __future__ import annotations

import logging
import shutil
import sqlite3
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from thyra.errors import ConversionRefused

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "fixtures" / "synthetic_tims.d"


def _copy_run(target: Path, frames: int = 0) -> Path:
    """A copy of the timsTOF fixture, trimmed to ``frames`` frames if given."""
    shutil.copytree(FIXTURE, target)
    if frames:
        with sqlite3.connect(target / "analysis.tdf") as conn:
            conn.execute("DELETE FROM MaldiFrameInfo WHERE Frame > ?", (frames,))
            conn.execute("DELETE FROM Frames WHERE Id > ?", (frames,))
        conn.close()
    for side in ("-wal", "-shm"):
        (target / f"analysis.tdf{side}").unlink(missing_ok=True)
    return target


@pytest.fixture
def slide(tmp_path: Path) -> Path:
    """A folder holding two runs; ``b_run.d`` keeps 3 of the 6 frames."""
    folder = tmp_path / "slide"
    _copy_run(folder / "a_run.d")
    _copy_run(folder / "b_run.d", frames=3)
    return folder


def _bruker_library() -> None:
    from thyra.readers.bruker.timstof.timstof_reader import BrukerReader
    from thyra.utils.bruker_exceptions import SDKError

    try:
        with BrukerReader(FIXTURE):
            pass
    except (SDKError, OSError) as exc:  # pragma: no cover - platform
        pytest.skip(f"Bruker library not loadable on this platform: {exc}")


class TestTheLibraryRefuses:
    def test_convert_msi_writes_nothing_and_names_both_runs(
        self, slide, tmp_path, thyra_logs
    ):
        from thyra import convert_msi

        out = tmp_path / "out.zarr"
        with thyra_logs("thyra", logging.INFO) as records:
            ok = convert_msi(str(slide), str(out), pixel_size_um=10.0)

        assert ok is False
        assert not out.exists()
        assert not (tmp_path / "out.zarr.failed").exists()
        errors = [r.getMessage() for r in records if r.levelno == logging.ERROR]
        assert len(errors) == 1
        assert "a_run.d" in errors[0] and "b_run.d" in errors[0]
        assert "pass one run's path" in errors[0]
        assert not any("Found .d folder" in r.getMessage() for r in records)

    def test_the_timstof_reader_refuses(self, slide):
        from thyra.readers.bruker.timstof.timstof_reader import BrukerReader

        with pytest.raises(ConversionRefused, match="a_run.d, b_run.d"):
            BrukerReader(slide)

    def test_preview_is_unreadable_and_names_both_runs(self, slide):
        from thyra import preview_msi

        preview = preview_msi(slide)

        assert preview.readable is False
        assert "a_run.d" in preview.error and "b_run.d" in preview.error

    def test_thyra_metadata_exits_1_and_prints_no_document(self, slide):
        from thyra.metadata.schema.cli import metadata_command

        result = CliRunner().invoke(metadata_command, [str(slide)])

        assert result.exit_code == 1
        assert "Could not read" in result.output
        assert "a_run.d" in result.output and "b_run.d" in result.output
        assert "msi_metadata" not in result.output.split("Could not read")[0]

    def test_detection_still_names_the_first_runs_format(self, slide):
        """Detection does not refuse: only the entry points do."""
        from thyra.core.registry import detect_format

        assert detect_format(slide) == "bruker"


class TestTheCliAsks:
    def test_no_answer_ends_with_the_refusal_not_aborted(self, slide, tmp_path):
        from thyra.__main__ import main

        out = tmp_path / "out.zarr"
        result = CliRunner().invoke(main, [str(slide), str(out)], input="")

        assert result.exit_code == 1
        assert "Aborted!" not in result.output
        assert "a_run.d" in result.output and "b_run.d" in result.output
        assert "pipe its number in" in result.output
        assert "Initialized BrukerReader" not in result.output
        assert not out.exists()

    def test_a_piped_number_converts_that_run(self, slide, tmp_path):
        import spatialdata as sd

        from thyra.__main__ import main

        _bruker_library()
        out = tmp_path / "out.zarr"
        result = CliRunner().invoke(
            main, [str(slide), str(out), "--resample-bins", "500"], input="2\n"
        )

        assert result.exit_code == 0, result.output
        assert sd.read_zarr(str(out)).tables["msi_dataset_z0"].n_obs == 3


class TestAFolderHoldingOneRun:
    def test_it_is_read_as_that_run_without_a_question(self, tmp_path):
        from thyra.__main__ import _select_bruker_dataset

        folder = tmp_path / "holder"
        run = _copy_run(folder / "my_run.d")
        (folder / "empty.d").mkdir()

        assert _select_bruker_dataset(folder) == run

    def test_the_calibration_display_is_shown(self, tmp_path):
        from thyra.__main__ import main

        _bruker_library()
        folder = tmp_path / "holder"
        _copy_run(folder / "my_run.d")
        (folder / "empty.d").mkdir()

        result = CliRunner().invoke(
            main,
            [str(folder), str(tmp_path / "out.zarr"), "--interactive-calibration"],
            input="",
        )

        assert result.exit_code == 0, result.output
        assert "Calibration Information" in result.output

    def test_a_solarix_run_converts_like_the_d_itself(self, tmp_path):
        import spatialdata as sd

        from tests.unit.readers.test_solarix_reader import make_solarix_d
        from thyra import convert_msi

        folder = tmp_path / "holder"
        folder.mkdir()
        run = make_solarix_d(folder)
        config = {"method": "nearest_neighbor", "target_bins": 200}

        assert convert_msi(str(run), str(tmp_path / "d.zarr"), resampling_config=config)
        assert convert_msi(
            str(folder), str(tmp_path / "folder.zarr"), resampling_config=config
        )

        direct = sd.read_zarr(str(tmp_path / "d.zarr")).tables["msi_dataset_z0"].X
        via = sd.read_zarr(str(tmp_path / "folder.zarr")).tables["msi_dataset_z0"].X
        assert np.array_equal(direct.toarray(), via.toarray())
