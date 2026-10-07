# tests/unit/converters/test_float32_imzml_stores.py
"""Stores of 32-bit imzML files are computed in float64.

The imzML reader widens 32-bit arrays to float64 as it yields them (D36),
so a store of a 32-bit file must behave like a store of its 64-bit twin:
the TIC image equals each pixel's row sum on the native axis, a
``tic_preserving`` row totals its pixel's measured total, an intensity
threshold drops values the typed number is below, and a typed m/z bound
drops the peaks outside it.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import spatialdata
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)
from thyra.readers.imzml.imzml_reader import ImzMLReader


def _convert(
    imzml_path: Path,
    output_path: Path,
    resampling_config,
    intensity_threshold=None,
) -> None:
    reader = ImzMLReader(imzml_path, intensity_threshold=intensity_threshold)
    try:
        assert (
            StreamingSpatialDataConverter(
                reader=reader,
                output_path=output_path,
                dataset_id="f32",
                pixel_size_um=10.0,
                resampling_config=resampling_config,
            ).convert()
            is True
        )
    finally:
        reader.close()


def _row_sums(output_path: Path) -> np.ndarray:
    table = ad.read_zarr(str(output_path / "tables" / "f32_z0"))
    return np.asarray(table.X.sum(axis=1), dtype=np.float64).ravel()


class TestTicImageEqualsRowSums:
    """``--no-resample`` on a 32-bit file: each TIC cell is its row's sum."""

    def test_every_pixel_matches_its_own_row(self, tmp_path):
        src = tmp_path / "tic32.imzML"
        with ImzMLWriter(src, mode="processed", intensity_dtype=np.float32) as w:
            # [2^24, 1, 1, 1]: float32 adds this to 2^24, float64 to 2^24 + 3
            w.addSpectrum(
                np.array([150.0, 250.0, 350.0, 450.0]),
                np.array([2.0**24, 1, 1, 1], np.float32),
                (1, 1, 1),
            )
            rng = np.random.default_rng(7)
            for x in (2, 3):
                w.addSpectrum(
                    np.linspace(100, 900, 1000),
                    rng.lognormal(6, 2, 1000).astype(np.float32),
                    (x, 1, 1),
                )

        out = tmp_path / "tic32.zarr"
        _convert(src, out, resampling_config=None)

        sdata = spatialdata.read_zarr(str(out))
        tic = np.asarray(sdata.images["f32_z0_tic"].data)[0]
        table = ad.read_zarr(str(out / "tables" / "f32_z0"))
        rows = np.asarray(table.X.sum(axis=1), dtype=np.float64).ravel()
        img = tic[table.obs["y"].to_numpy(int), table.obs["x"].to_numpy(int)]

        np.testing.assert_array_equal(img, rows)
        assert img[0] == 2.0**24 + 3


class TestTicPreservingTotals:
    def test_a_row_totals_its_pixels_measured_total(self, tmp_path):
        rng = np.random.default_rng(7)
        mz = np.linspace(150.0, 900.0, 3000)
        spectra = [
            rng.lognormal(5.0, 2.0, mz.size).astype(np.float32) for _ in range(3)
        ]
        src = tmp_path / "tp32.imzML"
        with ImzMLWriter(
            src, mode="continuous", spec_type="profile", intensity_dtype=np.float32
        ) as w:
            for x, s in enumerate(spectra):
                w.addSpectrum(mz, s, (x + 1, 1, 1))

        out = tmp_path / "tp32.zarr"
        _convert(src, out, resampling_config={"method": "tic_preserving"})

        table = ad.read_zarr(str(out / "tables" / "f32_z0"))
        rows = np.asarray(table.X.sum(axis=1), dtype=np.float64).ravel()
        # The measured total: the float32 values, added in float64
        exact = [s.astype(np.float64).sum() for s in spectra]
        want = np.array([exact[int(x)] for x in table.obs["x"]])

        np.testing.assert_allclose(rows, want, rtol=1e-12)


class TestIntensityThreshold:
    def test_no_stored_value_lies_below_the_threshold(self, tmp_path):
        src = tmp_path / "th32.imzML"
        mz = np.array([150.0, 250.0, 350.0])
        it = np.array([0.7, 4.0, 2.0], np.float32)  # float32(0.7) < 0.7
        with ImzMLWriter(src, mode="processed", intensity_dtype=np.float32) as w:
            for x in (1, 2):
                w.addSpectrum(mz, it, (x, 1, 1))

        out = tmp_path / "th32.zarr"
        _convert(
            src,
            out,
            resampling_config=None,
            intensity_threshold=0.7,
        )

        table = ad.read_zarr(str(out / "tables" / "f32_z0"))
        stored = table.X.data
        assert table.X.nnz == 4  # the 4.0 and 2.0 of two pixels
        assert stored[stored < 0.7].size == 0


class TestTypedMzBound:
    """A bound the file stores at float32 precision applies in float64."""

    ROW = np.array([1000.0, 1.0, 1.0])

    def _write(self, tmp_path: Path, dtype) -> Path:
        src = tmp_path / f"bound_{np.dtype(dtype).name}.imzML"
        # The middle m/z differs per pixel so the second spectrum misses the
        # shared-axis cache and takes the generic in-range mask.
        per_pixel = [np.array([400.0, 450.0, 500.0]), np.array([400.0, 450.5, 500.0])]
        with ImzMLWriter(
            src,
            mode="processed",
            spec_type="centroid",
            mz_dtype=dtype,
            intensity_dtype=np.float64,
        ) as w:
            for x, mzs in enumerate(per_pixel, start=1):
                w.addSpectrum(mzs.astype(dtype), self.ROW, (x, 1, 1))
        return src

    def _convert_with_bound(self, src: Path, tmp_path: Path, method: str) -> np.ndarray:
        out = tmp_path / f"{src.stem}_{method}.zarr"
        _convert(
            src,
            out,
            resampling_config={
                "method": method,
                "min_mz": 400.00001,  # float32(400.00001) is 400.0
                "max_mz": 600.0,
                "target_bins": 50,
            },
        )
        return _row_sums(out)

    def test_nearest_neighbor_drops_the_out_of_range_peak(self, tmp_path):
        for dtype in (np.float32, np.float64):
            rows = self._convert_with_bound(
                self._write(tmp_path, dtype), tmp_path, "nearest_neighbor"
            )
            np.testing.assert_array_equal(rows, [2.0, 2.0])

    def test_tic_preserving_drops_a_pixel_whose_only_peak_is_outside(self, tmp_path):
        # Pixel 1's peaks are inside the window; pixel 2 holds one peak at
        # m/z 400.0, outside it: only pixel 1 reaches the store.
        src = tmp_path / "lone32.imzML"
        with ImzMLWriter(
            src,
            mode="processed",
            spec_type="centroid",
            mz_dtype=np.float32,
            intensity_dtype=np.float64,
        ) as w:
            w.addSpectrum(np.array([450.0, 500.0]), np.array([1.0, 1.0]), (1, 1, 1))
            w.addSpectrum(np.array([400.0]), np.array([1000.0]), (2, 1, 1))

        out = tmp_path / "lone32.zarr"
        _convert(
            src,
            out,
            resampling_config={
                "method": "tic_preserving",
                "min_mz": 400.00001,
                "max_mz": 600.0,
                "target_bins": 50,
            },
        )

        table = ad.read_zarr(str(out / "tables" / "f32_z0"))
        assert table.n_obs == 1
        assert float(np.asarray(table.X.sum(axis=1)).ravel()[0]) == pytest.approx(2.0)
