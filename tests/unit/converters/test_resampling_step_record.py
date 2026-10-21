"""The stored resampling step records what the axis was made with (D56).

The step said what was asked for: a width the 100-bin floor overrode, the
default width beside a bin count that fixed the axis, and a width law on an
axis that does not use one.
"""

import json
from pathlib import Path

import anndata
import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.convert import convert_msi


@pytest.fixture
def source(tmp_path) -> Path:
    path = tmp_path / "src.imzML"
    with ImzMLWriter(str(path), mode="processed", spec_type="centroid") as writer:
        for x in (1, 2):
            writer.addSpectrum(
                np.array([100.0, 500.0, 1000.0]), np.array([1.0, 2.0, 3.0]), (x, 1)
            )
    return path


def _step(source: Path, out: Path, config) -> dict:
    assert convert_msi(
        str(source),
        str(out),
        dataset_id="ds",
        pixel_size_um=10,
        resampling_config=config,
    )
    table = anndata.read_zarr(str(out / "tables" / "ds_z0"))
    steps = json.loads(table.uns["msi_metadata"]["processing"])
    return next(s for s in steps if s["name"] == "mass axis resampling")["parameters"]


def test_a_floored_width_records_the_width_laid(source, tmp_path):
    config = {
        "axis_type": "constant",
        "width_at_mz": 1.0,
        "min_mz": 500.0,
        "max_mz": 510.0,
    }

    step = _step(source, tmp_path / "floor.zarr", config)

    assert step["target_bins"] == 100
    assert step["mass_width_da"] == pytest.approx(10.0 / 99)


def test_a_bin_count_records_the_width_it_lays(source, tmp_path):
    step = _step(
        source, tmp_path / "count.zarr", {"axis_type": "constant", "target_bins": 100}
    )

    # The source spans 100-1000 m/z: 100 points, 99 steps.
    assert step["mass_width_da"] == pytest.approx(900.0 / 99)


def test_a_width_law_off_a_tof_axis_is_not_recorded(source, tmp_path):
    """It is ignored there, with a warning, so it made nothing."""
    config = {"axis_type": "constant", "tof_a": 0.5, "tof_b": 1e-8}

    step = _step(source, tmp_path / "law.zarr", config)

    assert step["axis_type"] == "constant"
    assert not {"tof_a", "tof_b", "bins_per_fwhm"} & set(step)


@pytest.mark.parametrize(
    "method, recorded",
    [("nearest_neighbor", False), ("tic_preserving", True)],
)
def test_a_gap_tolerance_is_recorded_only_where_it_acts(
    source, tmp_path, method, recorded
):
    """Only interpolation reads it (issue #519)."""
    config = {
        "method": method,
        "axis_type": "constant",
        "target_bins": 100,
        "gap_tolerance_da": 0.5,
    }

    step = _step(source, tmp_path / f"{method}.zarr", config)

    assert ("gap_tolerance_da" in step) is recorded
