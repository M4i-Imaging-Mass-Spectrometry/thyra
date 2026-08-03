# tests/unit/test_cli_z_spacing.py
"""``--z-spacing`` on the command line.

The flag is only a transport: it becomes ``convert_msi(z_spacing_um=...)``,
and everything that matters about the value happens further down (see
``tests/unit/converters/test_3d_tic_z_spacing.py``). What is worth pinning here
is the transport itself, plus the two ways a user can hold it wrong:

* **Passing it without ``--handle-3d``.** Without 3D handling each slice is
  written as its own 2D image and there is no z axis, so the value is inert. A
  silently ignored calibration flag is the same class of defect as
  ``--optimize-chunks``, which accepted input and did nothing for two releases.
* **Passing a non-positive value.** Rejected before any file is opened, so a
  multi-hour conversion cannot fail on it at the end.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

import thyra.__main__ as cli
from thyra.__main__ import main


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def captured_call(monkeypatch):
    """Record convert_msi's kwargs instead of running a conversion."""
    seen: dict = {}

    def _fake_convert_msi(*args, **kwargs):
        seen.update(kwargs)
        return True

    monkeypatch.setattr(cli, "convert_msi", _fake_convert_msi)
    return seen


def _invoke(runner, tmp_path: Path, *extra: str):
    """Run the CLI against an input that clears its path validation.

    The ``.ibd`` sidecar has to exist: the CLI checks for it before any of the
    conversion runs, and without it every case below would fail on the path
    rather than on the thing it is testing.
    """
    source = tmp_path / "in.imzML"
    source.write_text("not parsed -- convert_msi is stubbed", encoding="utf-8")
    (tmp_path / "in.ibd").write_bytes(b"")
    return runner.invoke(main, [str(source), str(tmp_path / "out.zarr"), *extra])


class TestTransport:
    def test_the_value_reaches_convert_msi(self, runner, tmp_path, captured_call):
        result = _invoke(runner, tmp_path, "--handle-3d", "--z-spacing", "45")

        assert result.exit_code == 0, result.output
        assert captured_call["z_spacing_um"] == 45.0
        assert captured_call["handle_3d"] is True

    def test_omitting_it_forwards_none(self, runner, tmp_path, captured_call):
        """``None`` is what tells the converter to fall back and say so.

        Forwarding a number here instead would erase the distinction between
        an assumed spacing and one somebody chose.
        """
        result = _invoke(runner, tmp_path, "--handle-3d")

        assert result.exit_code == 0, result.output
        assert captured_call["z_spacing_um"] is None

    def test_a_fractional_spacing_survives(self, runner, tmp_path, captured_call):
        """Section thicknesses are not integers."""
        _invoke(runner, tmp_path, "--handle-3d", "--z-spacing", "12.5")

        assert captured_call["z_spacing_um"] == 12.5


class TestRejectedInput:
    @pytest.mark.parametrize("bad", ["0", "-5", "-0.1"])
    def test_a_non_positive_spacing_is_rejected(self, runner, tmp_path, bad):
        result = _invoke(runner, tmp_path, "--handle-3d", "--z-spacing", bad)

        assert result.exit_code != 0

    def test_it_is_rejected_before_the_input_is_opened(
        self, runner, tmp_path, captured_call
    ):
        """Validation runs ahead of the conversion, not after it."""
        _invoke(runner, tmp_path, "--handle-3d", "--z-spacing", "-1")

        assert captured_call == {}, "convert_msi ran despite an invalid z spacing"

    def test_a_non_numeric_spacing_is_rejected_by_click(self, runner, tmp_path):
        result = _invoke(runner, tmp_path, "--z-spacing", "thick")

        assert result.exit_code != 0
        assert "thick" in result.output


class TestFlagSurface:
    def test_help_advertises_it(self, runner):
        result = runner.invoke(main, ["--help"])

        assert result.exit_code == 0
        assert "--z-spacing" in result.output

    def test_help_points_from_handle_3d_to_it(self, runner):
        """Someone reaching for ``--handle-3d`` has to be able to find it."""
        result = runner.invoke(main, ["--help"])

        assert "--z-spacing" in result.output
        handle_3d_line = next(
            line for line in result.output.splitlines() if "--handle-3d" in line
        )
        assert "z-spacing" in handle_3d_line
