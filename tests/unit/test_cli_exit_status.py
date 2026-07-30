# tests/unit/test_cli_exit_status.py
"""Tests for the exit status the ``thyra`` command reports to the shell.

A conversion that fails must not look like a success to a calling script
or CI job, and it must not leave an unreadable store sitting at the
destination path where it can be mistaken for a finished conversion.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from thyra.__main__ import main


@pytest.fixture
def runner():
    return CliRunner()


def _invoke(runner, imzml_path, output_path, *extra):
    return runner.invoke(
        main,
        [str(imzml_path), str(output_path), "--pixel-size", "1.0", *extra],
    )


class TestExitStatus:
    """The process exit code must track conversion success."""

    def test_failed_conversion_exits_nonzero(
        self, create_minimal_imzml, temp_dir, monkeypatch, runner
    ):
        imzml_path, _, _, _ = create_minimal_imzml
        output_path = temp_dir / "out.zarr"

        monkeypatch.setattr("thyra.__main__.convert_msi", lambda *a, **k: False)

        result = _invoke(runner, imzml_path, output_path)

        assert result.exit_code != 0, (
            "a failed conversion must report a non-zero exit status; "
            f"got {result.exit_code}"
        )

    def test_successful_conversion_exits_zero(
        self, create_minimal_imzml, temp_dir, monkeypatch, runner
    ):
        imzml_path, _, _, _ = create_minimal_imzml
        output_path = temp_dir / "out.zarr"

        monkeypatch.setattr("thyra.__main__.convert_msi", lambda *a, **k: True)

        result = _invoke(runner, imzml_path, output_path)

        assert result.exit_code == 0, result.output

    def test_real_failing_conversion_exits_one(self, temp_dir, runner):
        """End-to-end: an undetectable input format must exit 1.

        ``convert_msi`` catches the format-detection error and returns
        False, so this exercises the real CLI path rather than a
        monkeypatched stub. Exit code 1 rather than click's 2 confirms the
        failure came from the conversion, not from argument parsing.
        """
        unknown_input = temp_dir / "not_an_msi_file.txt"
        unknown_input.write_text("not MSI data")
        output_path = temp_dir / "out.zarr"

        result = _invoke(runner, unknown_input, output_path)

        assert result.exit_code == 1, result.output


class TestPartialOutputQuarantine:
    """A failed conversion must not leave an unreadable store in place."""

    def test_partial_output_is_moved_aside(
        self, create_minimal_imzml, temp_dir, monkeypatch, runner
    ):
        imzml_path, _, _, _ = create_minimal_imzml
        output_path = temp_dir / "out.zarr"

        def fake_convert(*_args, **_kwargs):
            # Stand in for a conversion that dies part-way through the
            # zarr write, leaving an incomplete store behind.
            output_path.mkdir(parents=True)
            (output_path / "zarr.json").write_text("{}")
            return False

        monkeypatch.setattr("thyra.__main__.convert_msi", fake_convert)

        result = _invoke(runner, imzml_path, output_path)

        assert result.exit_code != 0
        assert not output_path.exists(), (
            "the destination path must be cleared so the partial store "
            "cannot be mistaken for a finished conversion"
        )
        assert (temp_dir / "out.zarr.failed").is_dir()

    def test_quarantine_does_not_overwrite_an_earlier_failure(
        self, create_minimal_imzml, temp_dir, monkeypatch, runner
    ):
        imzml_path, _, _, _ = create_minimal_imzml
        output_path = temp_dir / "out.zarr"
        (temp_dir / "out.zarr.failed").mkdir()

        def fake_convert(*_args, **_kwargs):
            output_path.mkdir(parents=True)
            return False

        monkeypatch.setattr("thyra.__main__.convert_msi", fake_convert)

        result = _invoke(runner, imzml_path, output_path)

        assert result.exit_code != 0
        assert not output_path.exists()
        assert (temp_dir / "out.zarr.failed2").is_dir()

    def test_nothing_to_quarantine_is_not_an_error(
        self, create_minimal_imzml, temp_dir, monkeypatch, runner
    ):
        """Failing before anything is written must still exit non-zero."""
        imzml_path, _, _, _ = create_minimal_imzml
        output_path = temp_dir / "out.zarr"

        monkeypatch.setattr("thyra.__main__.convert_msi", lambda *a, **k: False)

        result = _invoke(runner, imzml_path, output_path)

        assert result.exit_code != 0
        assert not output_path.exists()
        assert not (temp_dir / "out.zarr.failed").exists()
        assert result.exception is None or isinstance(result.exception, SystemExit)


class TestVersionOption:
    """``thyra --version`` must report the installed package version."""

    def test_version_flag_reports_package_version(self, runner):
        from thyra import __version__

        result = runner.invoke(main, ["--version"])

        assert result.exit_code == 0
        assert __version__ in result.output

    def test_version_flag_does_not_require_arguments(self, runner):
        """--version must work without INPUT and OUTPUT."""
        result = runner.invoke(main, ["--version"])

        assert result.exit_code == 0
        assert "Missing argument" not in result.output
