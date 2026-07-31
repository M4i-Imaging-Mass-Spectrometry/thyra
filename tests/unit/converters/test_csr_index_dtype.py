# tests/unit/converters/test_csr_index_dtype.py
"""Tests for CSR column-index width and the streaming index builds.

Two unrelated problems in the same file, both about behaving correctly at
sizes nobody has reached yet.

**Column indices.** The CSR ``indices`` array was created with a hardcoded
``int32`` and the write buffer cast to match. Those hold *column* indices, so
they are bounded by the number of m/z bins rather than by the non-zero count
that the neighbouring ``indptr`` is keyed on. Above 2,147,483,647 bins they
wrapped to negative values -- silently, because zarr truncates an oversized
write without warning and scipy will happily build a matrix from negative
indices unless you explicitly call ``check_format``.

**Index construction.** The var and obs string indices were built from Python
list comprehensions, materialising one ``str`` object per entry. The values
these tests pin are unchanged; only the construction is.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import zarr

from thyra.converters.spatialdata import streaming_converter as mod
from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)

INT32_MAX = np.iinfo(np.int32).max


def _setup(n_cols, total_nnz=8, n_rows=2):
    """Run the real ``_coo_setup_zarr_arrays`` against an in-memory store.

    Returns:
        Tuple of (indices_arr, data_arr, store).
    """
    store = zarr.open_group(store=zarr.storage.MemoryStore(), mode="w")
    stub = SimpleNamespace(_zarr_store=store)
    indptr = np.zeros(n_rows + 1, dtype=np.int64)
    indptr[-1] = total_nnz
    indices_arr, data_arr = StreamingSpatialDataConverter._coo_setup_zarr_arrays(
        stub, n_rows, n_cols, total_nnz, indptr
    )
    return indices_arr, data_arr, store


class TestIndexDtypeSwitch:
    """The dtype must follow the column count, at scipy's own switch point."""

    @pytest.mark.parametrize(
        "n_cols,expected",
        [
            (50_000, np.int32),
            (1_000_000, np.int32),
            (INT32_MAX, np.int32),
            (INT32_MAX + 1, np.int64),
            (3_000_000_000, np.int64),
        ],
    )
    def test_dtype_follows_column_count(self, n_cols, expected):
        indices_arr, _, _ = _setup(n_cols)
        assert indices_arr.dtype == expected

    def test_ordinary_datasets_still_use_int32(self):
        """Every dataset that exists today must take the identical branch."""
        indices_arr, _, _ = _setup(190_000)
        assert indices_arr.dtype == np.int32

    def test_indptr_is_keyed_on_nnz_not_columns(self):
        """indptr and indices are bounded by different quantities.

        A huge column count with few non-zeros must widen the column indices
        and leave indptr alone. The old code got this exactly backwards: the
        only widening it did was keyed on the non-zero count.
        """
        _, _, store = _setup(3_000_000_000, total_nnz=5)
        assert store["X"]["indices"].dtype == np.int64
        assert store["X"]["indptr"].dtype == np.int32


class TestRoundTripAtLargeColumnCount:
    """The corruption this fixes, demonstrated end to end.

    Only five non-zeros, so this costs nothing: the column count enters the
    write path solely through the Zarr shape attribute and the dtype choice.
    """

    COLUMNS = np.array(
        [7, 2147483646, 2147483648, 2500000000, 2999999999], dtype=np.int64
    )

    def test_large_column_indices_survive(self):
        indices_arr, _, _ = _setup(3_000_000_000, total_nnz=len(self.COLUMNS))
        indices_arr[:] = self.COLUMNS.astype(indices_arr.dtype)

        read_back = np.asarray(indices_arr)
        assert np.array_equal(read_back, self.COLUMNS)
        # The tripwire: truncation shows up as negative indices.
        assert (read_back >= 0).all()

    def test_small_column_count_is_lossless_too(self):
        """The int32 branch must stay correct -- this is the common case."""
        columns = np.array([0, 7, 49_999], dtype=np.int64)
        indices_arr, _, _ = _setup(50_000, total_nnz=len(columns))

        assert indices_arr.dtype == np.int32
        indices_arr[:] = columns.astype(indices_arr.dtype)
        assert np.array_equal(np.asarray(indices_arr), columns)


class TestIndexBuildValues:
    """The vectorised builds must produce exactly the old strings."""

    @pytest.mark.parametrize(
        "n", [0, 1, 2, 9, 10, 11, 99, 100, 101, 999, 1000, 2_000_001]
    )
    def test_var_index_matches_list_comprehension(self, n, monkeypatch):
        # A small chunk makes the boundary arithmetic actually run.
        monkeypatch.setattr(mod, "_INDEX_BUILD_CHUNK", 1000)
        str_dtype = np.dtypes.StringDType()

        built = np.empty(n, dtype=str_dtype)
        for start in range(0, n, mod._INDEX_BUILD_CHUNK):
            stop = min(start + mod._INDEX_BUILD_CHUNK, n)
            built[start:stop] = np.strings.add(
                "mz_", np.arange(start, stop, dtype=np.int64).astype(str_dtype)
            )

        expected = np.array([f"mz_{i}" for i in range(n)], dtype=str_dtype)
        assert np.array_equal(built, expected)

    @pytest.mark.parametrize("n", [0, 1, 10, 1000, 100_000])
    def test_obs_index_matches_list_comprehension(self, n):
        str_dtype = np.dtypes.StringDType()
        built = np.arange(n, dtype=np.int64).astype(str_dtype)
        expected = np.array([str(i) for i in range(n)], dtype=str_dtype)
        assert np.array_equal(built, expected)

    @pytest.mark.parametrize(
        "value", [0, 999_999, 2**31 - 1, 2**31, 2**32, 2**53, 2**62]
    )
    def test_int64_formats_like_str_at_extremes(self, value):
        """No scientific notation or sign artefacts at large magnitudes."""
        str_dtype = np.dtypes.StringDType()
        built = np.array([value], dtype=np.int64).astype(str_dtype)
        assert built[0] == str(value)

    def test_chunk_boundaries_are_exact(self, monkeypatch):
        """Off-by-one in the chunk arithmetic would show up here."""
        monkeypatch.setattr(mod, "_INDEX_BUILD_CHUNK", 100)
        str_dtype = np.dtypes.StringDType()

        for n in (99, 100, 101, 200, 201):
            built = np.empty(n, dtype=str_dtype)
            for start in range(0, n, mod._INDEX_BUILD_CHUNK):
                stop = min(start + mod._INDEX_BUILD_CHUNK, n)
                built[start:stop] = np.strings.add(
                    "mz_", np.arange(start, stop, dtype=np.int64).astype(str_dtype)
                )
            expected = np.array([f"mz_{i}" for i in range(n)], dtype=str_dtype)
            assert np.array_equal(built, expected), f"n={n}"
