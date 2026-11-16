"""Per-spectrum statistics in every summed table (design decision D90, issue #564).

Each summed MSI table gains five ``obs`` columns -- ``tic``, ``sum_sq``,
``nnz``, ``max_intensity`` and ``median_all`` -- and a
``uns["spectrum_statistics"]`` block. They describe the row as stored in
``X`` over every bin, the implicit zeros of a sparse row included, so each
is checked here against a dense recomputation from ``X`` read back off
disk, not against the numbers the converter held. The cases that can go
wrong are the ones where the stored row is not simply the spectrum the
source gave: a resampled row, a position the source measured twice (the
row is the sum, and exists only after the assembly merged them), a position
left empty, a volume, a plane of a multi-plane source and a store of
several regions. The sibling tables must carry none of it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest
import scipy.sparse as sp

pytest.importorskip("spatialdata")

import spatialdata  # noqa: E402

import thyra  # noqa: E402
from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader  # noqa: E402
from tests.unit.converters.test_fused_passes import (  # noqa: E402
    FusedStubReader,
    UnfusedStubReader,
)
from tests.unit.converters.test_fused_passes import (  # noqa: E402
    _convert as _convert_siblings,
)
from tests.unit.converters.test_fused_passes import (  # noqa: E402
    _read as _read_siblings,
)
from thyra.converters.spatialdata import csc_assembly  # noqa: E402
from thyra.converters.spatialdata.spectrum_statistics import (  # noqa: E402
    RowStatistics,
    median_over_all_bins,
    statistics_block,
)
from thyra.converters.spatialdata.streaming_converter import (  # noqa: E402
    StreamingSpatialDataConverter,
    _column_means,
    _column_sums_and_means,
)
from thyra.core.base_extractor import MetadataExtractor  # noqa: E402
from thyra.core.base_reader import BaseMSIReader  # noqa: E402
from thyra.metadata.schema import MSI_OBS_SPECTRUM_STATISTICS_COLUMNS  # noqa: E402
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata  # noqa: E402

NAMES = ["tic", "sum_sq", "nnz", "max_intensity", "median_all"]
RESAMPLED = {
    "method": "nearest_neighbor",
    "axis_type": "constant",
    "target_bins": 50,
    "min_mz": 100.0,
    "max_mz": 2000.0,
}


# ---------------------------------------------------------------------------
# The arithmetic, against numpy on the dense row
# ---------------------------------------------------------------------------


def _dense_row(values: np.ndarray, n_bins: int, rng) -> np.ndarray:
    """``values`` placed at random bins of an otherwise-zero row."""
    row = np.zeros(n_bins)
    row[rng.choice(n_bins, size=values.size, replace=False)] = values
    return row


class TestMedianOverAllBins:
    @pytest.mark.parametrize("n_bins", [1, 2, 3, 6, 7, 10, 11])
    def test_it_is_numpy_median_of_the_dense_row(self, n_bins):
        rng = np.random.default_rng(n_bins)
        for k in range(n_bins + 1):
            for _ in range(8):
                values = rng.normal(size=k) * 10.0
                row = _dense_row(values, n_bins, rng)
                assert median_over_all_bins(values, n_bins) == np.median(row)

    def test_negative_values_and_an_explicit_zero(self):
        n_bins = 6
        rng = np.random.default_rng(3)
        for values in (
            np.array([-5.0, -1.0, 0.0, 2.0]),
            np.array([-5.0, -1.0, -3.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 4.0, 9.0, 1.0, 2.0]),
            np.array([-2.0, -7.0, -3.0, -1.0]),
        ):
            row = _dense_row(values, n_bins, rng)
            assert median_over_all_bins(values, n_bins) == np.median(row)

    def test_a_sparse_row_has_median_zero_and_no_bins_gives_zero(self):
        assert median_over_all_bins(np.array([5.0, 6.0]), 10) == 0.0
        assert median_over_all_bins(np.array([]), 0) == 0.0

    def test_two_huge_middle_values_do_not_overflow(self):
        values = np.full(4, 1.7e308)
        assert median_over_all_bins(values, 4) == 1.7e308


class TestRowStatisticsRecord:
    def _check(self, values: np.ndarray, n_bins: int, rng) -> None:
        row = _dense_row(values, n_bins, rng)
        stats = RowStatistics(1, n_bins)
        stats.record(0, values)
        assert stats.sum_sq[0] == pytest.approx(float(np.sum(row * row)), rel=1e-12)
        assert stats.nnz[0] == np.count_nonzero(row)
        assert stats.max_intensity[0] == row.max()
        assert stats.median_all[0] == np.median(row)

    def test_a_row_against_the_dense_recomputation(self):
        rng = np.random.default_rng(5)
        for n_bins in (4, 9, 16):
            for k in range(n_bins + 1):
                self._check(rng.uniform(1.0, 50.0, size=k), n_bins, rng)

    def test_an_all_bins_negative_row_has_a_negative_maximum(self):
        n_bins = 5
        values = np.array([-4.0, -1.0, -9.0, -2.0, -3.0])
        stats = RowStatistics(1, n_bins)
        stats.record(0, values)
        assert stats.max_intensity[0] == -1.0
        assert stats.median_all[0] == -3.0

    def test_a_sparse_negative_row_has_maximum_zero(self):
        stats = RowStatistics(1, 8)
        stats.record(0, np.array([-4.0, -1.0]))
        assert stats.max_intensity[0] == 0.0
        assert stats.nnz[0] == 2

    def test_an_empty_row_is_all_zero_and_dtypes_are_fixed(self):
        stats = RowStatistics(2, 8)
        stats.record(1, np.array([]))
        assert stats.sum_sq.dtype == np.float64
        assert stats.max_intensity.dtype == np.float64
        assert stats.median_all.dtype == np.float64
        assert stats.nnz.dtype == np.uint32
        assert (stats.sum_sq == 0).all() and (stats.nnz == 0).all()

    def test_the_scratch_buffer_grows_with_the_largest_row_only(self):
        stats = RowStatistics(3, 1_000_000)
        assert stats._scratch.size == 0
        stats.record(0, np.ones(5))
        stats.record(1, np.ones(3))
        assert stats._scratch.size == 5


class TestRecomputeRows:
    def test_chunked_groups_equal_record_on_the_same_rows(self):
        rng = np.random.default_rng(9)
        n_rows, n_bins = 14, 9
        dense = rng.normal(size=(n_rows, n_bins)) * 5.0
        dense[rng.random(dense.shape) < 0.5] = 0.0
        dense[3] = 0.0  # a stored-empty row
        csc = sp.csc_matrix(dense)
        csr = sp.csr_matrix(dense)
        rows = np.array([1, 3, 4, 8, 11, 13], dtype=np.int64)

        expected = RowStatistics(n_rows, n_bins)
        for row in rows:
            expected.record(int(row), csr.data[csr.indptr[row] : csr.indptr[row + 1]])
        # Seven entries per group: several groups, one row alone.
        got = RowStatistics(n_rows, n_bins)
        got.recompute_rows(csc, rows, chunk_entries=7)

        np.testing.assert_array_equal(got.sum_sq, expected.sum_sq)
        np.testing.assert_array_equal(got.nnz, expected.nnz)
        np.testing.assert_array_equal(got.max_intensity, expected.max_intensity)
        np.testing.assert_array_equal(got.median_all, expected.median_all)
        for row in rows:
            assert got.median_all[row] == np.median(dense[row])

    def test_the_budget_is_read_from_the_assembly_at_call_time(self, monkeypatch):
        """Rows are taken in groups that fit the assembly's sort budget."""

        class _Counting:
            """A matrix that counts how often rows are selected from it."""

            def __init__(self, matrix):
                self._matrix = matrix
                self.indices = matrix.indices
                self.indptr = matrix.indptr
                self.selections = 0

            def __getitem__(self, key):
                self.selections += 1
                return self._matrix[key]

        dense = np.arange(1.0, 13.0).reshape(4, 3)
        rows = np.array([0, 2, 3])

        whole = _Counting(sp.csc_matrix(dense))
        RowStatistics(4, 3).recompute_rows(whole, rows)
        assert whole.selections == 1

        # Three entries a row and a budget of four: every row on its own.
        monkeypatch.setattr(csc_assembly, "SORT_CHUNK_ENTRIES", 4)
        split = _Counting(sp.csc_matrix(dense))
        stats = RowStatistics(4, 3)
        stats.recompute_rows(split, rows)
        assert split.selections == 3
        assert stats.sum_sq[2] == float(np.sum(dense[2] ** 2))
        assert stats.sum_sq[1] == 0.0  # not asked for

    def test_no_rows_is_a_no_op(self):
        stats = RowStatistics(2, 3)
        stats.recompute_rows(sp.csc_matrix(np.ones((2, 3))), np.array([], np.int64))
        assert (stats.sum_sq == 0).all()


class TestColumnSums:
    def test_the_mean_is_the_sum_divided(self):
        matrix = sp.csc_matrix(np.array([[1.0, 0.0], [2.5, 4.0], [0.1, 0.3]]))
        sums, means = _column_sums_and_means(matrix, 3)
        assert np.array_equal(means, sums / 3)
        assert np.array_equal(_column_means(matrix, 3), means)

    def test_an_overflowing_sum_keeps_its_infinity_and_the_mean_stays_finite(self):
        matrix = sp.csc_matrix(np.array([[1.7e308, 1.0], [1.7e308, 2.0]]))
        with np.errstate(over="ignore"):
            sums, means = _column_sums_and_means(matrix, 2)
        assert np.isinf(sums[0]) and sums[1] == 3.0
        assert np.isfinite(means).all() and means[0] == 1.7e308

    def test_the_block_names_its_columns_as_one_json_string(self):
        block = statistics_block(np.array([1.0, 2.0]))
        assert list(block) == ["version", "columns", "column_sum", "computed_by"]
        assert block["version"] == 1
        assert json.loads(block["columns"]) == NAMES
        assert block["computed_by"] == f"thyra/{thyra.__version__}"


# ---------------------------------------------------------------------------
# Whole conversions, read back from disk
# ---------------------------------------------------------------------------

AXIS = np.linspace(100.0, 800.0, 8)

Spectra = List[Tuple[Tuple[int, int, int], List[int], List[float]]]


class _Extractor(MetadataExtractor):
    def __init__(self, dims: Tuple[int, int, int]):
        super().__init__(None)
        self._dims = dims

    def _extract_essential_impl(self) -> EssentialMetadata:
        n = self._dims[0] * self._dims[1] * self._dims[2]
        return EssentialMetadata(
            dimensions=self._dims,
            coordinate_bounds=(
                0.0,
                float(self._dims[0] - 1),
                0.0,
                float(self._dims[1] - 1),
            ),
            mass_range=(100.0, 800.0),
            pixel_size=None,
            n_spectra=n,
            total_peaks=8 * n,
            source_path="/mock/statistics",
        )

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        return ComprehensiveMetadata(
            essential=self._extract_essential_impl(),
            format_specific={},
            acquisition_params={},
            instrument_info={},
            raw_metadata={},
        )


class ListReader(BaseMSIReader):
    """Yields exactly the spectra it is given, on an eight-bin axis."""

    def __init__(self, dims: Tuple[int, int, int], spectra: Spectra):
        super().__init__(Path("/mock/statistics"))
        self._dims = dims
        self._spectra = list(spectra)

    def _create_metadata_extractor(self) -> MetadataExtractor:
        return _Extractor(self._dims)

    def get_common_mass_axis(self) -> np.ndarray:
        return AXIS

    def iter_spectra(self):
        for coords, bins, values in self._spectra:
            yield coords, AXIS[np.asarray(bins, dtype=int)], np.asarray(values, float)

    def close(self) -> None:
        pass


class Float32Reader(ListReader):
    """The same spectra, with intensities left in the 32-bit type of a file."""

    def iter_spectra(self):
        for coords, mzs, values in super().iter_spectra():
            yield coords, mzs, values.astype(np.float32)


class OrderedListReader(ListReader):
    """The same spectra from a reader that knows their acquisition order."""

    @property
    def has_acquisition_order(self) -> bool:
        return True

    def iter_spectra_with_acquisition_order(self):
        for order, (coords, mzs, values) in enumerate(self.iter_spectra()):
            yield coords, 100 + order, mzs, values


def _convert(reader, out: Path, **kwargs) -> Path:
    converter = StreamingSpatialDataConverter(
        reader=reader,
        output_path=out,
        dataset_id="m",
        pixel_size_um=10.0,
        include_optical=False,
        **kwargs,
    )
    assert converter.convert() is True
    return out


def _read(store: Path):
    return spatialdata.read_zarr(str(store))


def _dense(table) -> np.ndarray:
    X = table.X
    return np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=np.float64)


def _tic_at_rows(sdata, key: str, table) -> np.ndarray:
    """The TIC image's value at each row's own position."""
    image = np.asarray(sdata.images[f"{key}_tic"].data)[0]
    obs = table.obs
    if image.ndim == 3:
        return image[obs["z"].to_numpy(), obs["y"].to_numpy(), obs["x"].to_numpy()]
    return image[obs["y"].to_numpy(), obs["x"].to_numpy()]


def _check_table(sdata, key: str) -> None:
    """Every statistic of one table against the dense matrix on disk."""
    table = sdata.tables[key]
    X = _dense(table)
    n_rows, n_bins = X.shape
    obs = table.obs

    columns = list(obs.columns)
    assert columns[-6:-1] == NAMES and columns[-1] == "instance_key"
    assert obs["tic"].dtype == np.float64
    assert obs["sum_sq"].dtype == np.float64
    assert obs["nnz"].dtype == np.uint32
    assert obs["max_intensity"].dtype == np.float64
    assert obs["median_all"].dtype == np.float64

    np.testing.assert_allclose(obs["tic"], X.sum(axis=1), rtol=1e-12)
    np.testing.assert_allclose(obs["sum_sq"], (X * X).sum(axis=1), rtol=1e-12)
    assert (obs["nnz"].to_numpy() == np.count_nonzero(X, axis=1)).all()
    np.testing.assert_allclose(obs["max_intensity"], X.max(axis=1), rtol=1e-12)
    np.testing.assert_allclose(obs["median_all"], np.median(X, axis=1), rtol=1e-12)
    # The same floats as the image, not merely close to them.
    assert (obs["tic"].to_numpy() == _tic_at_rows(sdata, key, table)).all()

    block = table.uns["spectrum_statistics"]
    # Zarr lists a group's members sorted, so the order is pinned on the
    # block itself (TestColumnSums), not on what is read back.
    assert set(block) == {"version", "columns", "column_sum", "computed_by"}
    assert int(block["version"]) == 1
    assert json.loads(block["columns"]) == NAMES
    assert block["computed_by"] == f"thyra/{thyra.__version__}"
    column_sum = np.asarray(block["column_sum"])
    assert column_sum.dtype == np.float64 and column_sum.shape == (n_bins,)
    np.testing.assert_allclose(column_sum, X.sum(axis=0), rtol=1e-12)
    assert np.array_equal(
        np.asarray(table.uns["average_spectrum"]), column_sum / n_rows
    )


def _mock(n_z: int = 1, **kwargs) -> MockMSIReader:
    return MockMSIReader(
        MockMSIConfig(
            n_x=4,
            n_y=3,
            n_z=n_z,
            n_mz_bins=64,
            peaks_per_spectrum=(3, 9),
            seed=21,
            **kwargs,
        )
    )


class TestEndToEnd:
    def test_a_resampled_plane(self, tmp_path):
        store = _convert(_mock(), tmp_path / "a.zarr", resampling_config=RESAMPLED)
        sdata = _read(store)
        assert set(sdata.tables) == {"m_z0"}
        _check_table(sdata, "m_z0")

    def test_a_repeated_position_is_the_sum_of_its_spectra(self, tmp_path):
        first = ([0, 1, 2, 3, 4, 5, 6, 7], [5.0, 1.0, 2.0, 8.0, 3.0, 1.0, 4.0, 2.0])
        second = ([2, 3, 4], [10.0, 20.0, 30.0])
        spectra: Spectra = [
            ((0, 0, 0), [0, 1], [2.0, 3.0]),
            ((1, 0, 0), *first),
            ((0, 1, 0), [3, 4, 5], [1.0, 2.0, 3.0]),
            ((1, 1, 0), [6, 7], [9.0, 1.0]),
            ((1, 0, 0), *second),
        ]
        store = _convert(ListReader((2, 2, 1), spectra), tmp_path / "b.zarr")
        sdata = _read(store)
        _check_table(sdata, "m_z0")

        table = sdata.tables["m_z0"]
        assert table.n_obs == 4
        row = table.obs.index.get_loc("1")
        summed = _dense(table)[row]
        # The sum, bin by bin: 5 1 12 28 33 1 4 2.
        np.testing.assert_array_equal(
            summed, [5.0, 1.0, 12.0, 28.0, 33.0, 1.0, 4.0, 2.0]
        )
        assert table.obs["max_intensity"].iloc[row] == 33.0
        assert table.obs["sum_sq"].iloc[row] == float(np.sum(summed**2))
        for single in (first, second):
            values = np.asarray(single[1])
            assert table.obs["max_intensity"].iloc[row] != values.max()
            assert table.obs["sum_sq"].iloc[row] != float(np.sum(values**2))

    def test_empty_positions_have_no_row_and_a_zero_image(self, tmp_path):
        spectra: Spectra = [
            ((0, 0, 0), [0, 1, 2, 3, 4, 5, 6, 7], [1.0] * 8),
            ((1, 0, 0), [0, 1, 2], [4.0, 5.0, 6.0]),
            # measured, but left empty
            ((0, 1, 0), [0, 1], [0.0, 0.0]),
            ((2, 1, 0), [5, 6, 7], [2.0, 2.0, 3.0]),
            # (2, 0), (1, 1) are never given a spectrum
        ]
        store = _convert(ListReader((3, 2, 1), spectra), tmp_path / "c.zarr")
        sdata = _read(store)
        _check_table(sdata, "m_z0")

        table = sdata.tables["m_z0"]
        assert table.n_obs == 3
        image = np.asarray(sdata.images["m_z0_tic"].data)[0]
        for x, y in ((0, 1), (2, 0), (1, 1)):
            assert image[y, x] == 0.0
            assert not ((table.obs["x"] == x) & (table.obs["y"] == y)).any()

    def test_without_resampling_a_row_over_half_full_has_a_median(self, tmp_path):
        spectra: Spectra = [
            (
                (0, 0, 0),
                [0, 1, 2, 3, 4, 5, 6, 7],
                [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0],
            ),
            ((1, 0, 0), [0, 1, 2, 3, 4, 5], [2.0, 7.0, 1.0, 8.0, 2.0, 8.0]),
            ((0, 1, 0), [3], [5.0]),
            ((1, 1, 0), [1, 6], [4.0, 4.0]),
        ]
        store = _convert(ListReader((2, 2, 1), spectra), tmp_path / "d.zarr")
        sdata = _read(store)
        _check_table(sdata, "m_z0")

        median = sdata.tables["m_z0"].obs["median_all"].to_numpy()
        assert median[0] == 3.5 and median[1] == 2.0
        assert median[2] == 0.0 and median[3] == 0.0

    def test_a_volume(self, tmp_path):
        store = _convert(
            _mock(n_z=3),
            tmp_path / "e.zarr",
            handle_3d=True,
            resampling_config=RESAMPLED,
        )
        sdata = _read(store)
        assert set(sdata.tables) == {"m"}
        _check_table(sdata, "m")
        assert sdata.tables["m"].n_obs == 4 * 3 * 3

    def test_every_plane_of_a_multi_plane_source_has_its_own(self, tmp_path):
        store = _convert(_mock(n_z=3), tmp_path / "f.zarr", resampling_config=RESAMPLED)
        sdata = _read(store)
        assert set(sdata.tables) == {"m_z0", "m_z1", "m_z2"}
        sums = []
        for key in sdata.tables:
            _check_table(sdata, key)
            sums.append(
                np.asarray(sdata.tables[key].uns["spectrum_statistics"]["column_sum"])
            )
        assert not np.array_equal(sums[0], sums[1])

    def test_a_repeated_position_in_a_volume(self, tmp_path):
        """The repeat sits on the second plane, met after the first is full."""
        spectra: Spectra = [
            ((0, 0, 0), [0, 1, 2], [2.0, 3.0, 4.0]),
            ((1, 0, 0), [1, 5], [6.0, 1.0]),
            ((0, 0, 1), [0, 1, 2, 3, 4, 5], [5.0, 1.0, 2.0, 8.0, 3.0, 1.0]),
            ((1, 0, 1), [6, 7], [9.0, 1.0]),
            ((0, 0, 1), [2, 3, 7], [10.0, 20.0, 30.0]),
        ]
        store = _convert(
            ListReader((2, 1, 2), spectra), tmp_path / "v.zarr", handle_3d=True
        )
        sdata = _read(store)
        _check_table(sdata, "m")

        table = sdata.tables["m"]
        assert table.n_obs == 4
        obs = table.obs
        row = int(np.flatnonzero((obs["x"] == 0) & (obs["z"] == 1))[0])
        np.testing.assert_array_equal(
            _dense(table)[row], [5.0, 1.0, 12.0, 28.0, 3.0, 1.0, 0.0, 30.0]
        )
        assert obs["max_intensity"].iloc[row] == 30.0
        assert obs["nnz"].iloc[row] == 7
        assert obs["median_all"].iloc[row] == 4.0

    def test_tic_preserving_resampling(self, tmp_path):
        store = _convert(
            _mock(),
            tmp_path / "t.zarr",
            resampling_config={**RESAMPLED, "method": "tic_preserving"},
        )
        _check_table(_read(store), "m_z0")

    def test_the_columns_follow_the_acquisition_order(self, tmp_path):
        spectra: Spectra = [
            ((0, 0, 0), [0, 1], [2.0, 3.0]),
            ((1, 0, 0), [3, 4, 5], [1.0, 2.0, 3.0]),
        ]
        store = _convert(OrderedListReader((2, 1, 1), spectra), tmp_path / "o.zarr")
        sdata = _read(store)
        _check_table(sdata, "m_z0")
        assert list(sdata.tables["m_z0"].obs.columns) == [
            "y",
            "x",
            "region",
            "spatial_x",
            "spatial_y",
            "region_number",
            "acquisition_order",
            *NAMES,
            "instance_key",
        ]

    def test_a_32_bit_source_without_resampling_keeps_the_image_number(self, tmp_path):
        """``tic`` is the TIC image's number, also where that is a 32-bit sum.

        Without resampling a spectrum reaches the converter in the type the
        reader yields, and pass 1 totals it in that type. A reader that
        yields 32-bit intensities (the Bruker TSF one) so gets an image
        holding a 32-bit sum, while ``X`` holds the values as float64. The
        column follows the image; the other four are taken from the stored
        values.
        """
        values = [16_000_000.25, 3_000_000.5, 7.125, 123_456.75, 9_999_999.5]
        spectra: Spectra = [
            ((0, 0, 0), [0, 1, 2, 3, 4], values),
            ((1, 0, 0), [2, 3, 4, 5, 6], values[::-1]),
        ]
        store = _convert(Float32Reader((2, 1, 1), spectra), tmp_path / "f.zarr")
        sdata = _read(store)

        table = sdata.tables["m_z0"]
        X = _dense(table)
        obs = table.obs
        assert (obs["tic"].to_numpy() == _tic_at_rows(sdata, "m_z0", table)).all()
        np.testing.assert_allclose(obs["tic"], X.sum(axis=1), rtol=1e-6)
        np.testing.assert_allclose(obs["sum_sq"], (X * X).sum(axis=1), rtol=1e-12)
        assert (obs["nnz"].to_numpy() == np.count_nonzero(X, axis=1)).all()
        np.testing.assert_array_equal(obs["max_intensity"], X.max(axis=1))
        np.testing.assert_array_equal(obs["median_all"], np.median(X, axis=1))

    def test_a_multi_region_store(self, tmp_path):
        class _TwoRegions(MockMSIReader):
            def get_region_map(self):
                return {
                    (x, y): (1 if x == 0 else 2) for y in range(3) for x in range(4)
                }

        reader = _TwoRegions(
            MockMSIConfig(n_x=4, n_y=3, n_mz_bins=64, peaks_per_spectrum=(3, 9), seed=2)
        )
        store = _convert(reader, tmp_path / "g.zarr", resampling_config=RESAMPLED)
        sdata = _read(store)
        _check_table(sdata, "m_z0")
        assert "average_spectrum_per_region" in sdata.tables["m_z0"].uns


class TestSiblingsCarryNone:
    @pytest.mark.parametrize("reader_type", [FusedStubReader, UnfusedStubReader])
    def test_mobility_and_msms_tables_have_none_of_it(self, tmp_path, reader_type):
        sdata = _read_siblings(_convert_siblings(reader_type(), tmp_path / "s.zarr"))

        assert set(sdata.tables) == {"stub_z0", "stub_z0_mobility", "stub_z0_msms"}
        summed = sdata.tables["stub_z0"]
        assert set(NAMES) <= set(summed.obs.columns)
        assert "spectrum_statistics" in summed.uns
        for key in ("stub_z0_mobility", "stub_z0_msms"):
            table = sdata.tables[key]
            for name in NAMES:
                assert name not in table.obs.columns
            assert "spectrum_statistics" not in table.uns
            assert list(table.obs.columns)[-1] == "instance_key"

    def test_the_constant_is_the_five_names_in_order(self):
        assert list(MSI_OBS_SPECTRUM_STATISTICS_COLUMNS) == NAMES
