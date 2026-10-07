"""A 32-bit imzML converts as exactly as its 64-bit twin (D36).

The imzML reader used to hand pyimzml's 32-bit arrays on unchanged, and
the converter, the threshold and the resamplers computed in that type:
the TIC image missed its own rows, ``tic_preserving`` kept a rounded
total, a threshold kept values just below it, and a typed m/z bound let
peaks just outside it in. Each test below converts a small 32-bit file
and checks the result in 64-bit arithmetic.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import anndata
import numpy as np
import zarr
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.convert import convert_msi

DATASET_ID = "ds"
TABLE = f"tables/{DATASET_ID}_z0"
TIC = f"images/{DATASET_ID}_z0_tic/s0"

Spectrum = Tuple[Sequence[float], Sequence[float]]


def _write(
    path: Path,
    spectra: List[Spectrum],
    mode: str = "processed",
    spec_type: str = "centroid",
    mz_dtype: Any = np.float64,
    intensity_dtype: Any = np.float64,
) -> Path:
    with ImzMLWriter(
        str(path),
        mode=mode,
        spec_type=spec_type,
        mz_dtype=mz_dtype,
        intensity_dtype=intensity_dtype,
    ) as writer:
        for x, (mzs, intensities) in enumerate(spectra, start=1):
            writer.addSpectrum(np.asarray(mzs), np.asarray(intensities), (x, 1, 1))
    return path


def _convert(
    src: Path,
    resampling: Optional[Dict[str, Any]],
    threshold: Optional[float] = None,
) -> anndata.AnnData:
    out = src.with_suffix(".zarr")
    ok = convert_msi(
        str(src),
        str(out),
        dataset_id=DATASET_ID,
        pixel_size_um=10,
        resampling_config=resampling,
        reader_options=(
            {"intensity_threshold": threshold} if threshold is not None else None
        ),
    )
    assert ok is True
    return anndata.read_zarr(str(out / TABLE))


def _row_sums(table: anndata.AnnData) -> np.ndarray:
    return np.asarray(table.X.sum(axis=1), dtype=np.float64).ravel()


def test_no_resample_tic_image_is_each_rows_sum(tmp_path):
    # 2**24 + 3 is not a float32: summed in 32-bit it came out as 2**24.
    rng = np.random.default_rng(7)
    spectra: List[Spectrum] = [([150.0, 250.0, 350.0, 450.0], [2.0**24, 1, 1, 1])]
    spectra += [
        (np.linspace(100, 900, 500), rng.lognormal(6, 2, 500)) for _ in range(3)
    ]
    src = _write(tmp_path / "f32.imzML", spectra, intensity_dtype=np.float32)
    table = _convert(src, resampling=None)

    tic = np.asarray(zarr.open_group(str(src.with_suffix(".zarr")), mode="r")[TIC])
    image = tic[0][table.obs["y"].to_numpy(int), table.obs["x"].to_numpy(int)]
    rows = _row_sums(table)
    assert rows.size == 4
    np.testing.assert_array_equal(image, rows)
    assert rows[table.obs["x"].to_numpy(int) == 0][0] == 2.0**24 + 3


def test_tic_preserving_keeps_the_measured_total(tmp_path):
    rng = np.random.default_rng(7)
    mz = np.linspace(150.0, 900.0, 3000)
    spectra = [rng.lognormal(5.0, 2.0, mz.size).astype(np.float32) for _ in range(3)]
    src = _write(
        tmp_path / "f32.imzML",
        [(mz, s) for s in spectra],
        mode="continuous",
        spec_type="profile",
        intensity_dtype=np.float32,
    )
    table = _convert(src, resampling={"method": "tic_preserving"})

    measured = np.array([spectra[x].astype(np.float64).sum() for x in table.obs["x"]])
    np.testing.assert_allclose(_row_sums(table), measured, rtol=1e-13, atol=0)


def test_threshold_drops_a_32_bit_value_just_below_it(tmp_path):
    # float32(0.7) is 0.699999988: below 0.7, so dropped as a 64-bit source
    # (an mzPeak archive of the same values) drops it.
    spectra: List[Spectrum] = [([150.0, 250.0, 350.0], [0.7, 4.0, 2.0])] * 3
    src = _write(tmp_path / "f32.imzML", spectra, intensity_dtype=np.float32)
    table = _convert(src, resampling=None, threshold=0.7)

    assert table.X.nnz == 6
    assert np.all(table.X.data >= 0.7)


class TestTypedMzBound:
    """float32(400.00001) is 400.0: a peak at 400.0 lies outside the bound."""

    RANGE = {"min_mz": 400.00001, "max_mz": 600.0, "target_bins": 50}
    CLEAN: Spectrum = ([450.0, 500.0], [1.0, 1.0])
    LONE: Spectrum = ([400.0], [1000.0])

    @staticmethod
    def _edge(x: int) -> Spectrum:
        # 1000 counts at m/z 400.0, below the typed minimum.
        return [400.0, 450.0 + x, 500.0], [1000.0, 1.0, 1.0]

    def _convert(self, tmp_path: Path, spectra: List[Spectrum], method: str):
        src = _write(tmp_path / "f32.imzML", spectra, mz_dtype=np.float32)
        return _convert(src, resampling={"method": method, **self.RANGE})

    def test_nearest_neighbor_drops_it_from_every_pixel(self, tmp_path):
        spectra = [self.CLEAN] + [self._edge(x) for x in range(2, 9)]
        table = self._convert(tmp_path, spectra, "nearest_neighbor")
        np.testing.assert_array_equal(_row_sums(table), np.full(8, 2.0))

    def test_nearest_neighbor_converts_when_the_first_pixel_has_it(self, tmp_path):
        # The shared-axis cache and the generic path used to disagree, and
        # the two passes refused the file.
        spectra = [self._edge(x) for x in range(1, 9)]
        table = self._convert(tmp_path, spectra, "nearest_neighbor")
        np.testing.assert_array_equal(_row_sums(table), np.full(8, 2.0))

    def test_tic_preserving_drops_a_pixel_whose_only_peak_is_outside(self, tmp_path):
        spectra = [self.CLEAN] + [self.LONE] * 4
        table = self._convert(tmp_path, spectra, "tic_preserving")
        assert table.n_obs == 1
        np.testing.assert_allclose(_row_sums(table), [2.0])
