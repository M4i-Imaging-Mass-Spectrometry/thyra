"""A store appears under its output name only when it is complete (#525).

The store used to be written in place. A process killed during the write
(an out-of-memory kill, a native abort) gives Python no chance to move the
partial store to ``.failed``, so it stayed under the output name: the
retry was refused because the output existed, and ``thyra validate``
blamed an old Thyra version. Now the store is written as
``<name>.incomplete-<8 hex>`` beside the output and renamed last.
"""

from __future__ import annotations

import logging
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
from click.testing import CliRunner

import thyra
from thyra import convert_msi
from thyra.utils.store_lifecycle import is_incomplete_store_name

_REPO_ROOT = Path(thyra.__file__).resolve().parent.parent

#: The real CLI in a child process, killed outright at the first write into
#: the store's tables -- what an out-of-memory kill does: no Python cleanup.
_KILLED_CHILD = textwrap.dedent("""
    import os, runpy, sys
    import zarr.storage

    original = zarr.storage.LocalStore.set

    async def set_then_die(self, key, value, *args, **kwargs):
        if key.startswith("tables/"):
            os._exit(137)
        return await original(self, key, value, *args, **kwargs)

    zarr.storage.LocalStore.set = set_then_die
    sys.argv = ["thyra"] + sys.argv[1:]
    runpy.run_module("thyra", run_name="__main__")
    """)


@pytest.fixture
def imzml(create_minimal_imzml) -> Path:
    path, _, _, _ = create_minimal_imzml
    return path


class TestAKilledRun:
    def test_leaves_the_output_name_free_for_the_same_command(self, imzml, tmp_path):
        out = tmp_path / "out.zarr"
        args = [str(imzml), str(out), "--pixel-size", "25"]

        killed = subprocess.run(
            [sys.executable, "-c", _KILLED_CHILD, *args],
            capture_output=True,
            text=True,
            cwd=_REPO_ROOT,
        )
        assert killed.returncode == 137, killed.stderr
        assert not out.exists()
        left = [p.name for p in tmp_path.iterdir() if p.name.startswith("out.zarr")]
        assert len(left) == 1 and is_incomplete_store_name(left[0])

        again = subprocess.run(
            [sys.executable, "-m", "thyra", *args],
            capture_output=True,
            text=True,
            cwd=_REPO_ROOT,
        )
        assert again.returncode == 0, again.stderr
        assert out.is_dir()


class TestTheWritePath:
    def test_no_write_lands_on_the_output_name(self, imzml, tmp_path, monkeypatch):
        import spatialdata
        import zarr
        import zarr.storage

        roots = []
        original = zarr.storage.LocalStore.set

        async def recording(self, key, value, *args, **kwargs):
            roots.append(Path(str(self.root)).name)
            return await original(self, key, value, *args, **kwargs)

        monkeypatch.setattr(zarr.storage.LocalStore, "set", recording)
        out = tmp_path / "out.zarr"

        assert convert_msi(str(imzml), str(out), pixel_size_um=25.0) is True

        assert roots
        assert "out.zarr" not in roots
        assert all(is_incomplete_store_name(root) for root in roots)
        assert [p.name for p in tmp_path.iterdir() if "out.zarr" in p.name] == [
            "out.zarr"
        ]
        spatialdata.read_zarr(str(out))
        zarr.open_consolidated(str(out), mode="r")

    def test_the_windows_budget_is_measured_on_the_temporary_name(
        self, imzml, tmp_path, monkeypatch
    ):
        from thyra.converters.spatialdata import base_spatialdata_converter as module

        measured = []
        original = module.prepare_zarr_output_path

        def spy(path, dataset_id):
            measured.append(Path(path).name)
            return original(path, dataset_id)

        monkeypatch.setattr(module, "prepare_zarr_output_path", spy)

        assert convert_msi(str(imzml), str(tmp_path / "out.zarr"), pixel_size_um=25.0)
        assert len(measured) == 1 and is_incomplete_store_name(measured[0])


class TestARenameThatFails:
    def test_fails_the_run_and_keeps_the_complete_store(
        self, imzml, tmp_path, monkeypatch, thyra_logs
    ):
        import spatialdata

        from thyra.converters.spatialdata import base_spatialdata_converter as module

        def refuse(action):
            raise PermissionError(13, "Access is denied")

        monkeypatch.setattr(module, "retry_on_os_error", refuse)
        out = tmp_path / "out.zarr"

        with thyra_logs("thyra", logging.ERROR) as records:
            assert convert_msi(str(imzml), str(out), pixel_size_um=25.0) is False

        assert not out.exists()
        assert not list(tmp_path.glob("out.zarr.failed*"))
        kept = [p for p in tmp_path.iterdir() if is_incomplete_store_name(p.name)]
        assert len(kept) == 1
        spatialdata.read_zarr(str(kept[0]))
        named = [r.getMessage() for r in records if str(kept[0]) in r.getMessage()]
        assert len(named) == 1 and "complete" in named[0]


class TestValidateNamesAnUnfinishedWrite:
    def _validate(self, path: Path):
        from thyra.metadata.schema.cli import validate_command

        return CliRunner().invoke(validate_command, [str(path)])

    def test_a_root_group_alone_names_both_causes(self, tmp_path):
        store = tmp_path / "out.zarr"
        store.mkdir()
        (store / "zarr.json").write_text('{"zarr_format": 3, "node_type": "group"}')

        result = self._validate(store)

        assert result.exit_code == 1
        assert "stopped before it finished writing" in result.output
        assert "without the metadata schema" in result.output

    def test_a_temporary_folder_is_named_as_one(self, tmp_path):
        store = tmp_path / "out.zarr.incomplete-0123abcd"
        store.mkdir()
        (store / "zarr.json").write_text('{"zarr_format": 3, "node_type": "group"}')

        result = self._validate(store)

        assert result.exit_code == 1
        assert "temporary folder of a conversion that did not finish" in result.output
