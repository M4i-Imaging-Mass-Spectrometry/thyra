# tests/unit/test_cli_spectrum_type.py
"""``--spectrum-type`` on the command line, the way SCiLS exposes ``--rep_type``.

The flag is only a transport: it builds ``reader_options["spectrum_type"]``,
which the imzML reader forwards to the metadata extractor. What is worth
pinning here is the transport, and in particular that ``auto`` -- the CLI's
spelling of "no override" -- is *omitted* rather than forwarded, because
``"auto"`` is not a spectrum representation and the reader would reject it.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from thyra.__main__ import _build_reader_options, main


@pytest.fixture
def runner():
    return CliRunner()


class TestBuildReaderOptions:
    def test_auto_is_not_forwarded(self):
        """Passing "auto" through would raise ValueError in the reader."""
        options = _build_reader_options(True, None, "auto")
        assert "spectrum_type" not in options

    def test_default_is_auto(self):
        """Callers that predate the flag must keep getting no override."""
        assert "spectrum_type" not in _build_reader_options(True, None)

    @pytest.mark.parametrize("value", ["profile", "centroid"])
    def test_an_explicit_choice_is_forwarded_verbatim(self, value):
        assert _build_reader_options(True, None, value)["spectrum_type"] == value

    def test_it_does_not_disturb_the_other_options(self):
        options = _build_reader_options(False, 12.5, "profile")
        assert options["use_recalibrated_state"] is False
        assert options["intensity_threshold"] == 12.5
        assert options["spectrum_type"] == "profile"


class TestFlagSurface:
    def test_help_advertises_the_choices(self, runner):
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "--spectrum-type" in result.output

    def test_an_invalid_choice_is_rejected_by_click(self, runner, tmp_path):
        """click.Choice rejects it before any file is touched."""
        result = runner.invoke(
            main,
            [
                str(tmp_path / "in.imzML"),
                str(tmp_path / "out.zarr"),
                "--spectrum-type",
                "profil",
            ],
        )
        assert result.exit_code != 0
        assert "profil" in result.output
