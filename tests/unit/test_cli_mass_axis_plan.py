"""``--mass-axis-plan FILE`` (issue #449, D43).

The plan fixes the mass axis, so the file becomes the whole resampling config
and every option that would also decide the axis is refused beside it. What is
pinned here is that transport and those refusals; building the plan and
converting with it are tested with the planner and the converter.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from thyra.__main__ import main
from thyra.resampling.plan import MassAxisPlan
from thyra.resampling.types import AxisType

_SHA = "ab" * 32


def _plan(method=None) -> MassAxisPlan:
    return MassAxisPlan(
        axis_type=AxisType.CONSTANT,
        min_mz=200.0,
        max_mz=900.0,
        reference_mz=1000.0,
        width_at_mz=0.01,
        n_bins=70000,
        sha256=_SHA,
        method=method,
    )


@pytest.fixture
def run(create_minimal_imzml, temp_dir, monkeypatch):
    """Invoke the CLI with a saved plan; returns ``(result, seen_config)``."""
    imzml_path, _, _, _ = create_minimal_imzml
    seen = {}

    def fake_convert(*_args, **kwargs):
        seen["config"] = kwargs.get("resampling_config")
        return True

    monkeypatch.setattr("thyra.convert.convert_msi", fake_convert)
    plan_file = temp_dir / "axis.json"

    def invoke(*extra, plan=None):
        (plan or _plan()).save(plan_file)
        result = CliRunner().invoke(
            main,
            [
                str(imzml_path),
                str(temp_dir / "out.zarr"),
                "--mass-axis-plan",
                str(plan_file),
                *extra,
            ],
        )
        return result, seen.get("config")

    return invoke


def test_the_plan_becomes_the_resampling_config(run):
    result, config = run()

    assert result.exit_code == 0, result.output
    assert config == _plan().to_resampling_config()
    assert config["mass_axis_sha256"] == _SHA


@pytest.mark.parametrize(
    "flags",
    [
        ["--no-resample"],
        ["--mass-axis-type", "orbitrap"],
        ["--resample-bins", "5000"],
        ["--resample-min-mz", "150"],
        ["--resample-max-mz", "950"],
        ["--resample-width-at-mz", "0.02"],
        ["--resample-reference-mz", "500"],
        ["--resample-reference-mz", "1000"],
        ["--tof-law", "0.02", "0.00001"],
    ],
    ids=lambda flags: flags[0] + (" " + flags[1] if len(flags) > 1 else ""),
)
def test_an_option_that_decides_the_axis_is_refused(run, flags):
    result, config = run(*flags)

    assert result.exit_code != 0
    assert flags[0] in result.output
    assert "mass-axis-plan" in result.output
    assert config is None


def test_a_method_is_allowed_when_the_plan_forces_none(run):
    result, config = run("--resample-method", "tic_preserving")

    assert result.exit_code == 0, result.output
    assert config["method"] == "tic_preserving"
    assert config["mass_axis_sha256"] == _SHA


def test_a_method_is_refused_when_the_plan_forces_one(run):
    from thyra.resampling.types import ResamplingMethod

    result, config = run(
        "--resample-method",
        "tic_preserving",
        plan=_plan(method=ResamplingMethod.NEAREST_NEIGHBOR),
    )

    assert result.exit_code != 0
    assert "nearest_neighbor" in result.output
    assert config is None


def test_the_method_the_plan_forces_may_be_repeated(run):
    from thyra.resampling.types import ResamplingMethod

    result, config = run(
        "--resample-method",
        "nearest_neighbor",
        plan=_plan(method=ResamplingMethod.NEAREST_NEIGHBOR),
    )

    assert result.exit_code == 0, result.output
    assert config["method"] == "nearest_neighbor"


def test_a_hand_edited_plan_is_a_usage_error_not_a_traceback(
    create_minimal_imzml, temp_dir
):
    import json

    imzml_path, _, _, _ = create_minimal_imzml
    data = _plan().to_dict()
    data["min_mz"] = "abc"
    bad = temp_dir / "bad.json"
    bad.write_text(json.dumps(data), encoding="utf-8")

    result = CliRunner().invoke(
        main,
        [str(imzml_path), str(temp_dir / "out.zarr"), "--mass-axis-plan", str(bad)],
    )

    assert result.exit_code == 2, result.output
    assert "min_mz" in result.output
    assert not isinstance(result.exception, ValueError)


def test_a_gap_tolerance_travels_with_the_plan(run):
    result, config = run("--resample-gap-tolerance", "0.5")

    assert result.exit_code == 0, result.output
    assert config["gap_tolerance_da"] == 0.5


def test_a_file_that_is_not_a_plan_is_refused(create_minimal_imzml, temp_dir):
    imzml_path, _, _, _ = create_minimal_imzml
    bad = temp_dir / "bad.json"
    bad.write_text("{}", encoding="utf-8")

    result = CliRunner().invoke(
        main,
        [str(imzml_path), str(temp_dir / "out.zarr"), "--mass-axis-plan", str(bad)],
    )

    assert result.exit_code != 0
    assert "Not a mass axis plan" in result.output
