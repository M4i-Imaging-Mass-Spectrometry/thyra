"""``--sample-id`` (issue #447): the user names the sample, the CLI hands it on."""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from thyra.__main__ import main


@pytest.fixture
def run(create_minimal_imzml, temp_dir, monkeypatch):
    """Invoke the CLI; returns ``(result, kwargs convert_msi was called with)``."""
    imzml_path, _, _, _ = create_minimal_imzml
    seen = {}

    def fake_convert(*_args, **kwargs):
        seen.update(kwargs)
        return True

    monkeypatch.setattr("thyra.convert.convert_msi", fake_convert)

    def invoke(*extra):
        result = CliRunner().invoke(
            main, [str(imzml_path), str(temp_dir / "out.zarr"), *extra]
        )
        return result, seen

    return invoke


def test_the_id_reaches_convert_msi(run):
    result, seen = run("--sample-id", "S1")

    assert result.exit_code == 0, result.output
    assert seen["sample_id"] == "S1"
    assert "sample_name" not in seen


def test_no_id_is_passed_as_none(run):
    result, seen = run()

    assert result.exit_code == 0, result.output
    assert seen["sample_id"] is None


def test_the_option_is_listed_under_other(run):
    result = CliRunner().invoke(main, ["--help"])

    assert "--sample-id" in result.output
