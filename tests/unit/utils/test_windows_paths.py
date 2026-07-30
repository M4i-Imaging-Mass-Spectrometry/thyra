# tests/unit/utils/test_windows_paths.py
"""Tests for extended-length output paths on Windows.

A Zarr store's deepest key sits far below the path the caller named, so an
output path well inside the 260 character Windows limit can still push an
internal key past it. The path rewriting is pure string work and is tested
on every platform; the end-to-end conversion is Windows-only.
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import pytest

from thyra.utils import windows_paths
from thyra.utils.windows_paths import (
    DEEPEST_KEY_RESERVE,
    WINDOWS_MAX_PATH,
    prepare_zarr_output_path,
    projected_deepest_key_length,
    to_extended_length_path,
)

EXTENDED_PREFIX = "\\\\?\\"


class TestToExtendedLengthPath:
    """Rewriting must produce syntax Windows actually accepts."""

    def test_drive_path_gets_the_prefix(self):
        result = to_extended_length_path(Path("C:\\data\\out.zarr"))

        assert str(result) == "\\\\?\\C:\\data\\out.zarr"

    def test_already_extended_path_is_unchanged(self):
        original = Path("\\\\?\\C:\\data\\out.zarr")

        assert str(to_extended_length_path(original)) == str(original)

    def test_applying_twice_is_idempotent(self):
        once = to_extended_length_path(Path("C:\\data\\out.zarr"))
        twice = to_extended_length_path(once)

        assert str(twice) == str(once)

    def test_unc_share_uses_the_unc_spelling(self):
        """A bare prefix on a UNC path is invalid syntax."""
        result = to_extended_length_path(Path("\\\\server\\share\\out.zarr"))

        assert str(result) == "\\\\?\\UNC\\server\\share\\out.zarr"
        assert not str(result).startswith("\\\\?\\\\\\")


class TestProjectedLength:
    """The projection must account for the store's internal depth."""

    def test_includes_the_dataset_id_and_reserve(self):
        out = Path("C:\\data\\out.zarr")

        projected = projected_deepest_key_length(out, "msi_dataset")

        assert projected == len(str(out)) + 1 + len("msi_dataset") + DEEPEST_KEY_RESERVE

    def test_longer_dataset_id_projects_further(self):
        out = Path("C:\\data\\out.zarr")

        short = projected_deepest_key_length(out, "ds")
        long = projected_deepest_key_length(out, "a_much_longer_dataset_id")

        assert long > short

    def test_reserve_covers_the_observed_deepest_key(self):
        """The deepest key measured in a real conversion was 95 characters
        below the store root, excluding the dataset id."""
        assert DEEPEST_KEY_RESERVE >= 95


class TestPrepareZarrOutputPath:
    """The prefix must appear only when it is needed."""

    @pytest.fixture
    def on_windows(self, monkeypatch):
        monkeypatch.setattr(sys, "platform", "win32")
        monkeypatch.setattr(windows_paths, "_long_paths_enabled", lambda: False)

    def _long_path(self, dataset_id="msi_dataset"):
        """A path just long enough to need the prefix."""
        needed = WINDOWS_MAX_PATH - DEEPEST_KEY_RESERVE - len(dataset_id)
        return Path("C:\\" + "d" * (needed + 10))

    def _short_path(self):
        return Path("C:\\data\\out.zarr")

    def test_short_path_is_untouched_on_windows(self, on_windows):
        out = self._short_path()

        assert prepare_zarr_output_path(out, "msi_dataset") == out

    def test_long_path_is_extended_on_windows(self, on_windows):
        out = self._long_path()

        result = prepare_zarr_output_path(out, "msi_dataset")

        assert str(result).startswith(EXTENDED_PREFIX)
        assert str(result).endswith(str(out))

    def test_no_op_off_windows(self, monkeypatch):
        monkeypatch.setattr(sys, "platform", "linux")
        out = self._long_path()

        assert prepare_zarr_output_path(out, "msi_dataset") == out

    def test_no_op_when_long_paths_are_enabled(self, monkeypatch):
        monkeypatch.setattr(sys, "platform", "win32")
        monkeypatch.setattr(windows_paths, "_long_paths_enabled", lambda: True)
        out = self._long_path()

        assert prepare_zarr_output_path(out, "msi_dataset") == out

    def test_the_registry_is_only_consulted_when_needed(self, monkeypatch):
        """A short path must not pay for a registry read."""
        monkeypatch.setattr(sys, "platform", "win32")
        calls = {"n": 0}

        def counting():
            calls["n"] += 1
            return False

        monkeypatch.setattr(windows_paths, "_long_paths_enabled", counting)

        prepare_zarr_output_path(self._short_path(), "msi_dataset")

        assert calls["n"] == 0

    def test_a_longer_dataset_id_can_tip_a_path_over(self, on_windows):
        """The dataset id is part of the deepest key, so it matters."""
        dataset_id = "a" * 60
        needed = WINDOWS_MAX_PATH - DEEPEST_KEY_RESERVE
        out = Path("C:\\" + "d" * (needed - 30))

        with_short = prepare_zarr_output_path(out, "ds")
        with_long = prepare_zarr_output_path(out, dataset_id)

        assert not str(with_short).startswith(EXTENDED_PREFIX)
        assert str(with_long).startswith(EXTENDED_PREFIX)


@pytest.mark.skipif(sys.platform != "win32", reason="Windows path limit only")
class TestLongPathConversion:
    """End to end: a store deeper than MAX_PATH must still convert."""

    def test_conversion_succeeds_past_the_limit(self):
        from thyra.convert import convert_msi

        base = Path(tempfile.mkdtemp())
        try:
            deep = base
            while len(str(deep)) < 175:
                deep = deep / "nested_directory_segment"
                deep.mkdir(exist_ok=True)

            src = base / "src"
            src.mkdir()
            imzml = self._write_imzml(src)

            out = deep / "out.zarr"
            assert len(str(out)) > 165, "path must be long enough to matter"

            assert convert_msi(
                str(imzml), str(out), dataset_id="msi_dataset", pixel_size_um=2.5
            ), "a long output path must not fail the conversion"

            import spatialdata as sd

            sdata = sd.read_zarr(EXTENDED_PREFIX + str(out))
            assert "msi_dataset_z0" in sdata.tables
            assert len(sdata.images) >= 1
        finally:
            # rmtree itself trips the limit without the prefix.
            shutil.rmtree(EXTENDED_PREFIX + str(base), ignore_errors=True)

    @staticmethod
    def _write_imzml(directory: Path) -> Path:
        import numpy as np
        from pyimzml.ImzMLWriter import ImzMLWriter

        path = directory / "minimal.imzML"
        mzs = np.linspace(100, 1000, 50)
        with ImzMLWriter(str(path), mode="processed") as writer:
            for x, y, z in [(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)]:
                intensities = np.zeros_like(mzs)
                intensities[10] = 100.0 * x
                intensities[30] = 150.0 * y
                writer.addSpectrum(mzs, intensities, (x, y, z))
        return path
