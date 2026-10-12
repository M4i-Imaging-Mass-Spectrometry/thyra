"""The SCiLS Lab parity harness, checked on synthetic input.

``tests/parity/`` needs real reference data and skips without it, so its own
arithmetic would otherwise go unchecked in every run that lacks the data,
which is every CI run. These tests feed it arrays whose answer is known.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from tests.parity import metrics, report, thyra_runner
from tests.parity.datasets import PARITY_ENV, discover, parity_root
from tests.parity.scils_reference import load_reference


def _law_axis(exponent: float, lo: float = 100.0, hi: float = 1000.0, n: int = 5000):
    """An axis whose bin width grows as m/z ** exponent."""
    if exponent == 1.0:
        return np.geomspace(lo, hi, n)
    p = 1.0 - exponent
    return np.linspace(lo**p, hi**p, n) ** (1.0 / p)


@pytest.mark.parametrize("law, exponent", sorted(metrics.LAW_EXPONENTS.items()))
def test_spacing_exponent_recovers_each_law(law: str, exponent: float) -> None:
    fitted = metrics.spacing_exponent(_law_axis(exponent))
    assert fitted == pytest.approx(exponent, abs=0.01)
    assert metrics.nearest_law(fitted) == law


def test_axis_agreement_identical_axes() -> None:
    mz = _law_axis(1.0)
    result = metrics.axis_agreement(mz, mz.copy())
    assert result["bins_equal"] == 1.0
    assert result["max_ppm"] == 0.0
    assert result["max_ppm_in_bins"] == 0.0


def test_axis_agreement_half_bin_shift() -> None:
    """Bin centres against bin edges: every value half a bin away."""
    edges = np.geomspace(540.0, 600.0, 1001)
    centres = 0.5 * (edges[1:] + edges[:-1])
    result = metrics.axis_agreement(edges[:-1], centres)
    assert result["max_ppm_in_bins"] == pytest.approx(0.5, abs=1e-3)
    assert result["first_ppm"] > 0


def test_axis_agreement_with_different_bin_counts_pairs_nearest() -> None:
    ref = np.linspace(100.0, 200.0, 101)
    cand = np.linspace(100.0, 200.0, 201)  # every reference value is on it
    result = metrics.axis_agreement(ref, cand)
    assert result["bins_equal"] == 0.0
    assert result["max_abs_da"] == pytest.approx(0.0, abs=1e-9)


def test_intensity_agreement_identical_and_scaled() -> None:
    rng = np.random.default_rng(0)
    ref = rng.random(1000)
    same = metrics.intensity_agreement(ref, ref.copy())
    assert same["rel_l1"] == 0.0
    assert same["pearson_r"] == pytest.approx(1.0)
    assert same["sum_ratio"] == pytest.approx(1.0)
    doubled = metrics.intensity_agreement(ref, 2.0 * ref)
    assert doubled["sum_ratio"] == pytest.approx(2.0)
    assert doubled["rel_l1"] == pytest.approx(1.0)
    assert doubled["pearson_r"] == pytest.approx(1.0)


def test_intensity_agreement_refuses_different_shapes() -> None:
    with pytest.raises(ValueError):
        metrics.intensity_agreement(np.ones(3), np.ones(4))


def test_rebin_sum_keeps_counts_inside_and_drops_outside() -> None:
    mz = np.array([0.5, 1.5, 1.6, 2.5, 9.0])
    intensity = np.array([1.0, 2.0, 3.0, 4.0, 100.0])
    out = metrics.rebin_sum(mz, intensity, np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(out, [5.0, 4.0])


def test_ratio_summary_counts_zero_reference_pixels() -> None:
    out = metrics.ratio_summary(np.array([0.0, 1.0, 2.0]), np.array([5.0, 2.0, 4.0]))
    assert out["zero_reference"] == 1.0
    assert out["median_ratio"] == pytest.approx(2.0)


def test_grid_index_reduces_micrometres_to_steps() -> None:
    np.testing.assert_array_equal(
        metrics.grid_index(np.array([-6831.0, -6811.0, -6771.0])), [0, 1, 3]
    )
    np.testing.assert_array_equal(metrics.grid_index(np.array([4.0, 4.0])), [0, 0])


def _irregular_raster():
    """An L-shaped raster, so the four flips pair different numbers of pixels."""
    xs, ys = np.meshgrid(np.arange(6), np.arange(4))
    keep = ~((xs >= 3) & (ys >= 2))
    return xs[keep].astype(float), ys[keep].astype(float)


def test_match_pixels_finds_a_flipped_y_axis() -> None:
    x, y = _irregular_raster()
    # Reference in micrometres with y counted the other way; candidate as
    # grid indices with an offset.
    ref_x, ref_y = 20.0 * x - 500.0, 20.0 * (y.max() - y)
    cand_x, cand_y = x + 7, y + 3
    match = metrics.match_pixels(ref_x, ref_y, cand_x, cand_y)
    assert match.orientation == "flip_y"
    assert match.reference_index.size == x.size
    assert not match.tie_broken_by_tic
    np.testing.assert_array_equal(match.reference_index, match.candidate_index)


def test_match_pixels_breaks_a_full_rectangle_tie_by_tic() -> None:
    xs, ys = np.meshgrid(np.arange(5), np.arange(3))
    x, y = xs.ravel().astype(float), ys.ravel().astype(float)
    tic = np.arange(x.size, dtype=float) ** 2
    flipped = (y.max() - y) * 5 + x  # candidate row of each reference spot
    cand_tic = np.empty_like(tic)
    cand_tic[flipped.astype(int)] = tic
    match = metrics.match_pixels(x, y, x, y, tic, cand_tic)
    assert match.tie_broken_by_tic
    assert match.orientation == "flip_y"


def _write_reference(scils_dir: Path) -> None:
    scils_dir.mkdir(parents=True)
    mz = np.geomspace(540.0, 600.0, 50)
    spots = np.arange(4)
    rebinned = np.ones((2, mz.size))
    np.savez_compressed(
        scils_dir / "reference.npz",
        mz=mz,
        mean=np.ones(mz.size),
        mean_tic=np.ones(mz.size),
        spot_id=spots,
        spot_x=spots * 20.0,
        spot_y=np.zeros(4),
        spot_z=np.zeros(4),
        spot_sum=np.full(4, float(mz.size)),
        norm_tic=np.ones(4),
        norm_rms=np.ones(4),
        sample_spot_id=np.array([0, 3]),
        sample_rebinned=rebinned,
        sample_raw_mz_0=mz[:3],
        sample_raw_int_0=np.ones(3),
        sample_raw_mz_1=mz[:3],
        sample_raw_int_1=np.ones(3),
    )
    (scils_dir / "reference.json").write_text(
        json.dumps(
            {"import": {"massaxis": {"type": "OTOF"}, "rebinning": "Tic Preserving"}}
        ),
        encoding="utf-8",
    )


def test_discover_and_load_reference(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.delenv(PARITY_ENV, raising=False)
    assert parity_root() is None
    assert discover() == []

    root = tmp_path / "data"
    (root / "one").mkdir(parents=True)
    (root / "one" / "parity.json").write_text(
        json.dumps(
            {
                "source": "source/run.d",
                "spectrum": "centroid",
                "variants": {"other": {"reader_options": {"a": 1}}},
            }
        ),
        encoding="utf-8",
    )
    _write_reference(root / "one" / "scils")
    monkeypatch.setenv(PARITY_ENV, str(root))
    (dataset,) = discover()
    assert dataset.name == "one"
    assert dataset.has_reference
    assert dataset.source == root / "one" / "source" / "run.d"
    assert dataset.variants()["other"].reader_options == {"a": 1}
    assert dataset.reader_options == {}

    reference = load_reference(dataset.scils_dir)
    assert [s.spot_id for s in reference.samples] == [0, 3]
    assert set(reference.normalizations) == {"tic", "rms"}
    assert thyra_runner.import_axis_type(reference) == "reflector_tof"

    config = thyra_runner.matched_resampling_config(dataset, reference)
    assert config == {
        "method": "nearest_neighbor",
        "axis_type": "reflector_tof",
        "min_mz": 540.0,
        "max_mz": 600.0,
        "target_bins": 50,
    }


def _write_store(store: Path, dense: np.ndarray, mz: np.ndarray) -> None:
    """A minimal store holding one table, laid out as Thyra writes it."""
    import anndata as ad
    import pandas as pd
    import zarr
    from scipy import sparse

    n_obs = dense.shape[0]
    table = ad.AnnData(
        X=sparse.csc_matrix(dense),
        obs=pd.DataFrame(
            {
                "x": np.arange(n_obs, dtype=np.int32),
                "y": np.zeros(n_obs, dtype=np.int32),
            },
            index=[str(i) for i in range(n_obs)],
        ),
        var=pd.DataFrame({"mz": mz}, index=[f"mz_{i}" for i in range(mz.size)]),
        uns={"average_spectrum": dense.mean(axis=0)},
    )
    zarr.open_group(str(store / "tables"), mode="w")
    table.write_zarr(str(store / "tables" / "msi_dataset_z0"))


def test_read_output_streams_row_sums_and_rows(tmp_path: Path, monkeypatch) -> None:
    rng = np.random.default_rng(1)
    dense = rng.random((5, 40)) * (rng.random((5, 40)) > 0.6)
    mz = np.linspace(100.0, 200.0, 40)
    store = tmp_path / "s.zarr"
    _write_store(store, dense, mz)
    # A chunk far smaller than the table, so the streaming path is used.
    monkeypatch.setattr(thyra_runner, "_CHUNK", 7)
    out = thyra_runner.read_output(store, rows=[1, 4])
    np.testing.assert_allclose(out.tic, dense.sum(axis=1))
    np.testing.assert_allclose(out.sum_squares, (dense**2).sum(axis=1))
    np.testing.assert_allclose(out.rows[1], dense[1])
    np.testing.assert_allclose(out.rows[4], dense[4])
    np.testing.assert_allclose(out.mz, mz)
    np.testing.assert_allclose(out.mean, dense.mean(axis=0))


def test_verdicts_fail_on_unmeasurable_values() -> None:
    """A NaN that the JSON round trip turned into None fails, it does not raise."""
    axis = {
        "bins_candidate": 10.0,
        "bins_reference": 10.0,
        "bins_equal": 1.0,
        "max_ppm": 0.0,
        "max_ppm_in_bins": 0.0,
        "exponent_candidate": 1.0,
        "exponent_reference": 1.0,
    }
    empty = {"rel_l1": None, "pearson_r": None, "sum_ratio": None}
    result = {
        "comparison": {
            "axis": axis,
            "mean": {"binwise": empty, "coarse": empty},
            "pixels": {"matched": 0, "scils_spots": 4, "matched_fraction": 0.0},
            "tic": None,
            "spots": [],
        }
    }
    checks = {name: ok for name, _, ok in report.verdicts(result)}
    assert checks["axis bin count"]
    assert not any(
        ok for name, ok in checks.items() if name.startswith(("mean", "coarse"))
    )
    assert not any(ok for name, ok in checks.items() if name.startswith("matched"))


@pytest.mark.parametrize(
    "axis_type", ["linear_tof", "reflector_tof", "orbitrap", "fticr"]
)
def test_half_bin_widening_puts_thyra_midpoints_on_the_law_grid(axis_type: str) -> None:
    """With the range widened, Thyra's bin midpoints land on SCiLS's values.

    Thyra's generator for each law is the one under test here, so a change to
    how it places its values shows up as this offset moving off zero.
    """
    from thyra.resampling.mass_axis import (
        FTICRAxisGenerator,
        LinearTOFAxisGenerator,
        OrbitrapAxisGenerator,
        ReflectorTOFAxisGenerator,
    )

    generators = {
        "linear_tof": LinearTOFAxisGenerator,
        "reflector_tof": ReflectorTOFAxisGenerator,
        "orbitrap": OrbitrapAxisGenerator,
        "fticr": FTICRAxisGenerator,
    }
    forward, inverse = thyra_runner._LAW_COORDINATES[axis_type]
    lo, hi, bins = 540.0, 600.0, 2000
    u = np.linspace(forward(lo), forward(hi), bins)
    scils_values = np.array([inverse(v) for v in u])

    plain = generators[axis_type]().generate_axis(lo, hi, bins).mz_values
    wide_lo, wide_hi = thyra_runner._widen_by_half_a_bin(axis_type, lo, hi, bins)
    aligned = generators[axis_type]().generate_axis(wide_lo, wide_hi, bins).mz_values

    plain_offset = metrics.axis_agreement(scils_values, plain)["max_ppm_in_bins"]
    aligned_offset = metrics.axis_agreement(scils_values, aligned)["max_ppm_in_bins"]
    assert plain_offset == pytest.approx(0.5, abs=0.01)
    assert aligned_offset < 1e-3


def test_ratio_summary_drops_non_finite_and_counts_zero_candidates() -> None:
    """A NaN factor (an empty SCiLS spot) must not turn the summary into NaN."""
    ref = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    cand = np.array([1.0, 2.0, 3.0, 0.0, 5.0])
    out = metrics.ratio_summary(ref, cand)
    assert out["non_finite"] == 1.0
    assert out["zero_candidate"] == 1.0
    assert np.isfinite(out["pearson_r"])
    assert out["median_ratio"] == pytest.approx(1.0)


def test_grid_index_ignores_one_irregular_position() -> None:
    values = np.arange(0.0, 2000.0, 20.0)
    values[3] = 59.5  # one off-pitch position must not become the step
    idx = metrics.grid_index(values)
    assert idx[-1] == 99


def test_match_pixels_recovers_a_missing_edge_column() -> None:
    """The candidate lacks the reference's first column: pairs shift by one, not none."""
    xs, ys = np.meshgrid(np.arange(10), np.arange(4))
    x, y = xs.ravel().astype(float), ys.ravel().astype(float)
    keep = x > 0
    tic = np.arange(x.size, dtype=float) ** 1.5
    match = metrics.match_pixels(x, y, x[keep], y[keep], tic, tic[keep])
    assert match.shift == (-1, 0)
    assert match.reference_index.size == int(keep.sum())
    # Every pair joins the same physical pixel.
    np.testing.assert_array_equal(
        x[match.reference_index], x[keep][match.candidate_index]
    )
    np.testing.assert_array_equal(
        y[match.reference_index], y[keep][match.candidate_index]
    )


def test_compare_end_to_end_with_an_empty_sample_spectrum(tmp_path: Path) -> None:
    """compare() on a synthetic reference and store; an empty SCiLS spot is not an error."""
    scils_dir = tmp_path / "scils"
    _write_reference(scils_dir)
    with np.load(scils_dir / "reference.npz") as npz:
        arrays = {k: npz[k] for k in npz.files}
    arrays["sample_raw_mz_1"] = np.zeros(0)
    arrays["sample_raw_int_1"] = np.zeros(0)
    arrays["norm_tic"] = np.array([1.0, 1.0, np.nan, 1.0])
    np.savez_compressed(scils_dir / "reference.npz", **arrays)
    reference = load_reference(scils_dir)

    dense = np.ones((4, reference.mz.size))
    store = tmp_path / "s.zarr"
    _write_store(store, dense, reference.mz)
    output = thyra_runner.read_output(store, rows=[0, 3])
    result = report.compare(reference, output)

    assert result["pixels"]["matched"] == 4
    assert result["mean"]["binwise"]["rel_l1"] == pytest.approx(0.0)
    assert result["tic"]["median_ratio"] == pytest.approx(1.0)
    assert result["tic_factor"]["non_finite"] == 1.0
    spot = {s["spot_id"]: s for s in result["spots"]}
    assert spot[0]["coarse"]["sum_ratio"] == pytest.approx(1.0)
    assert np.isnan(spot[3]["coarse_vs_raw"]["pearson_r"])
    names = [name for name, _, _ in report.verdicts({"comparison": result})]
    assert "axis bin count" in names


def test_summed_table_is_chosen_by_name(tmp_path: Path) -> None:
    """A wider sibling table is not taken for the summed one; two summed tables are refused."""
    import anndata as ad
    import zarr
    from scipy import sparse

    store = tmp_path / "s.zarr"
    _write_store(store, np.ones((3, 10)), np.linspace(100.0, 200.0, 10))
    ad.AnnData(X=sparse.csc_matrix(np.ones((3, 50)))).write_zarr(
        str(store / "tables" / "msi_dataset_z0_mobility")
    )
    assert thyra_runner._summed_table(store) == "msi_dataset_z0"
    zarr.open_group(str(store / "tables"), mode="a").create_group("msi_dataset_z1")
    with pytest.raises(RuntimeError, match="one summed table"):
        thyra_runner._summed_table(store)


def test_read_output_handles_a_csr_table(tmp_path: Path, monkeypatch) -> None:
    import anndata as ad
    from scipy import sparse

    rng = np.random.default_rng(2)
    dense = rng.random((4, 30)) * (rng.random((4, 30)) > 0.5)
    mz = np.linspace(100.0, 200.0, 30)
    store = tmp_path / "s.zarr"
    _write_store(store, dense, mz)
    table = ad.read_zarr(str(store / "tables" / "msi_dataset_z0"))
    table.X = sparse.csr_matrix(dense)
    table.write_zarr(str(store / "tables" / "msi_dataset_z0"))
    monkeypatch.setattr(thyra_runner, "_CHUNK", 5)
    out = thyra_runner.read_output(store, rows=[2])
    np.testing.assert_allclose(out.tic, dense.sum(axis=1))
    np.testing.assert_allclose(out.rows[2], dense[2])


def test_matched_config_widens_only_when_asked(tmp_path: Path) -> None:
    from tests.parity.datasets import ParityDataset

    scils_dir = tmp_path / "scils"
    _write_reference(scils_dir)
    reference = load_reference(scils_dir)
    plain = ParityDataset("d", tmp_path, {"source": "x", "spectrum": "profile"})
    aligned = ParityDataset(
        "d",
        tmp_path,
        {"source": "x", "spectrum": "profile", "axis_alignment": "scils_values"},
    )
    a = thyra_runner.matched_resampling_config(plain, reference)
    b = thyra_runner.matched_resampling_config(aligned, reference)
    assert a["method"] == "tic_preserving"
    assert (a["min_mz"], a["max_mz"]) == (540.0, 600.0)
    assert b["min_mz"] < 540.0 < 600.0 < b["max_mz"]
    assert b["target_bins"] == a["target_bins"]
