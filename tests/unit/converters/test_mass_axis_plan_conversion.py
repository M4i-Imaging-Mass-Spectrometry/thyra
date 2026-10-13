"""A conversion given a mass axis plan builds the planned axis, or refuses (issue #449, D43).

Two small imzML sources cover different m/z ranges. Planned together, both
stores must carry the same ``var["mz"]``, and each must record the plan's
fingerprint (D42). A plan the conversion's axis does not match is refused
while the converter initialises, before any spectrum is read or anything is
written.
"""

from pathlib import Path
from typing import Tuple

import numpy as np
import pytest

spatialdata = pytest.importorskip("spatialdata")

from pyimzml.ImzMLWriter import ImzMLWriter  # noqa: E402

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader  # noqa: E402
from thyra.convert import convert_msi  # noqa: E402
from thyra.converters.spatialdata import SpatialDataConverter  # noqa: E402
from thyra.errors import ConversionRefused  # noqa: E402
from thyra.mass_axis_plan import plan_mass_axis  # noqa: E402
from thyra.metadata.schema import read_msi_metadata_blocks  # noqa: E402
from thyra.utils.windows_paths import prepare_zarr_read_path  # noqa: E402


def _write_imzml(path: Path, mz_range: Tuple[float, float]) -> Path:
    """A 2x2 processed imzML whose spectra span ``mz_range``."""
    mzs = np.linspace(mz_range[0], mz_range[1], 40)
    with ImzMLWriter(str(path), mode="processed") as writer:
        for i, (x, y) in enumerate([(1, 1), (1, 2), (2, 1), (2, 2)]):
            intensity = np.zeros_like(mzs)
            intensity[5 + i] = 100.0 + i
            intensity[25] = 50.0
            writer.addSpectrum(mzs, intensity, (x, y, 1))
    return path


def _var_mz(store: Path) -> np.ndarray:
    import zarr

    root = zarr.open_group(str(prepare_zarr_read_path(store)), mode="r")
    (table,) = list(root["tables"])
    return np.asarray(root["tables"][table]["var"]["mz"])


def _recorded(store: Path) -> dict:
    (block,) = read_msi_metadata_blocks(store).values()
    return block["provenance"]["mass_axis"]


@pytest.fixture
def sources(tmp_path):
    return (
        _write_imzml(tmp_path / "a.imzML", (100.0, 600.0)),
        _write_imzml(tmp_path / "b.imzML", (300.0, 1000.0)),
    )


def _convert(source: Path, out: Path, config) -> bool:
    return convert_msi(
        str(source),
        str(out),
        dataset_id="d",
        pixel_size_um=10.0,
        resampling_config=config,
    )


@pytest.mark.parametrize("as_dict", [False, True], ids=["plan", "dict"])
def test_both_sources_get_the_planned_axis(tmp_path, sources, as_dict):
    plan = plan_mass_axis(list(sources), width_at_mz=0.5)
    config = plan.to_resampling_config() if as_dict else plan
    stores = []

    for source in sources:
        out = tmp_path / f"{source.stem}.zarr"
        assert _convert(source, out, config) is True
        stores.append(out)

    first, second = (_var_mz(s) for s in stores)
    assert np.array_equal(first, second)
    assert first.size == plan.n_bins
    for store in stores:
        recorded = _recorded(store)
        assert recorded["sha256"] == plan.sha256
        assert int(recorded["n_bins"]) == plan.n_bins
    # The intersection, each imzML header bound moved out by one unit of its
    # last stated decimal: "300.0" and "600.0" state one decimal (#473).
    assert (plan.min_mz, plan.max_mz) == pytest.approx((299.9, 600.1))


def test_the_plan_leaves_the_resampling_step_as_without_it(tmp_path, sources):
    """D43: a plan changes nothing a store contains.

    The plan's hash is a check, not a setting; the axis's own hash is
    ``provenance.mass_axis`` (D42). So the "mass axis resampling" step of a
    planned conversion is the step of the same explicit request without it.
    """
    plan = plan_mass_axis(list(sources), width_at_mz=0.5)
    plain = plan.to_resampling_config()
    del plain["mass_axis_sha256"]
    source = sources[0]

    steps = []
    for name, config in (("planned", plan), ("plain", plain)):
        out = tmp_path / f"{name}.zarr"
        assert _convert(source, out, config) is True
        block = next(iter(read_msi_metadata_blocks(out).values()))
        step = [s for s in block["processing"] if s["name"] == "mass axis resampling"]
        steps.append(step)

    assert steps[0] == steps[1]
    assert "mass_axis_sha256" not in steps[0][0]["parameters"]


def test_a_plan_whose_fingerprint_is_altered_is_refused_before_anything_is_written(
    tmp_path,
):
    config = MockMSIConfig(n_x=3, n_y=3, n_mz_bins=64, peaks_per_spectrum=(4, 8))
    plan_config = {
        "axis_type": "constant",
        "min_mz": 100.0,
        "max_mz": 2000.0,
        "width_at_mz": 5.0,
        "mass_axis_sha256": "0" * 64,
    }
    out = tmp_path / "out.zarr"
    converter = SpatialDataConverter(
        MockMSIReader(config),
        out,
        dataset_id="d",
        pixel_size_um=10.0,
        resampling_config=plan_config,
    )

    with pytest.raises(ConversionRefused, match="does not describe") as err:
        converter._initialize_conversion()

    message = str(err.value)
    assert "000000000000" in message and "bins" in message
    assert not out.exists()


def test_convert_msi_returns_false_and_leaves_nothing_behind(
    tmp_path, sources, thyra_logs
):
    plan = plan_mass_axis(list(sources), width_at_mz=0.5)
    config = plan.to_resampling_config()
    config["mass_axis_sha256"] = "f" * 64
    out = tmp_path / "out.zarr"

    with thyra_logs() as logs:
        assert _convert(sources[0], out, config) is False

    assert "does not describe this conversion's axis" in logs.text
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("out")) == []


def test_a_conversion_given_a_different_range_than_the_plan_is_refused(
    tmp_path, sources
):
    plan = plan_mass_axis(list(sources), width_at_mz=0.5)
    config = plan.to_resampling_config()
    config["min_mz"] = plan.min_mz - 10.0

    assert _convert(sources[0], tmp_path / "out.zarr", config) is False
    assert not (tmp_path / "out.zarr").exists()


def test_a_plan_from_a_rounded_imzml_header_keeps_every_count(tmp_path):
    """The header states the observed range to 6 decimals; no edge sample is lost (#473)."""
    import re

    import zarr

    mzs = np.linspace(100.12345678, 999.87654321, 200)
    source = tmp_path / "src.imzML"
    with ImzMLWriter(str(source), mode="continuous", mz_dtype=np.float64) as writer:
        for x, y in [(1, 1), (1, 2), (2, 1), (2, 2)]:
            writer.addSpectrum(mzs, np.full(mzs.size, 5.0), (x, y, 1))
    text = source.read_text(encoding="utf-8")
    text = re.sub(
        r'(accession="MS:100052[78]"[^>]*?value=")([^"]*)"',
        lambda m: f'{m.group(1)}{float(m.group(2)):.6f}"',
        text,
    )
    source.write_text(text, encoding="utf-8")

    plan = plan_mass_axis(source)
    out = tmp_path / "out.zarr"
    assert _convert(source, out, plan) is True

    root = zarr.open_group(str(prepare_zarr_read_path(out)), mode="r")
    (table,) = list(root["tables"])
    stored = float(np.asarray(root["tables"][table]["X"]["data"]).sum(dtype=np.float64))
    assert stored == pytest.approx(4 * mzs.size * 5.0)


def test_a_fingerprint_that_is_not_hex_is_refused_by_the_constructor(tmp_path):
    with pytest.raises(ConversionRefused, match="mass_axis_sha256"):
        SpatialDataConverter(
            MockMSIReader(MockMSIConfig(n_x=2, n_y=2, n_mz_bins=64)),
            tmp_path / "out.zarr",
            dataset_id="d",
            pixel_size_um=10.0,
            resampling_config={"mass_axis_sha256": "not-a-hash"},
        )
