# thyra/converters/spatialdata/spectrum_statistics.py

"""Per-spectrum statistics of the summed table, taken without a third pass (D90).

What this module owns: the five ``obs`` columns every summed MSI table
carries -- ``tic``, ``sum_sq``, ``nnz``, ``max_intensity`` and
``median_all`` -- and the ``uns["spectrum_statistics"]`` block that names
them and keeps the column sums they were taken against. Each column
describes the row **as stored in X**, over all of the table's bins: a
sparse row's implicit zeros count, so ``median_all`` of a spectrum that is
more than half empty is 0 and ``max_intensity`` of a row of negative
values that leaves a bin unstored is 0 -- exactly what ``numpy`` gives for
the dense row.

Why it lives here and not in ``streaming_converter.py``: the arithmetic is
independent of how the table is assembled, and a test can hold it to a
dense recomputation without a conversion. The converter only decides when
a row is in hand. Pass 2 has each row's final values at the moment it
scatters them, so :meth:`RowStatistics.record` is called there and no
further pass over the data is needed. The one exception is a grid position
the source measured more than once: its stored row is the sum of several
spectra and exists only after the assembly merged them, so those rows are
taken from the finished matrix (:meth:`RowStatistics.recompute_rows`), a
bounded slice at a time.

What it deliberately does not do: touch ``X``, the TIC image or the
sibling tables, which carry none of these columns.
"""

import json
from typing import Any, Dict, Optional

import numpy as np
from numpy.typing import NDArray

from ...metadata.schema import (
    MSI_OBS_MAX_INTENSITY_COLUMN,
    MSI_OBS_MEDIAN_ALL_COLUMN,
    MSI_OBS_NNZ_COLUMN,
    MSI_OBS_SPECTRUM_STATISTICS_COLUMNS,
    MSI_OBS_SUM_SQ_COLUMN,
    MSI_OBS_TIC_COLUMN,
    MSI_SPECTRUM_STATISTICS_VERSION,
)
from . import csc_assembly


def median_over_all_bins(values: NDArray[np.float64], n_bins: int) -> float:
    """The median of a row of ``n_bins`` values, given only its stored ones.

    Equal to ``numpy.median`` of the dense row, without building it. The
    row's ``n_bins - values.size`` other bins are zeros, so in sorted order
    the dense row is the stored values below zero, then every zero (stored
    or not), then the stored values above zero. A middle position is
    therefore a zero, or one known position among the sorted stored
    values, and only those positions are selected (``numpy.partition``):
    nothing is sorted in full.

    Args:
        values: The row's stored values, in any order, any of them zero.
        n_bins: How many bins the row has; at least ``values.size``.

    Returns:
        The median; 0.0 for a row without bins.
    """
    if n_bins <= 0:
        return 0.0
    lo = (n_bins - 1) // 2
    hi = n_bins // 2
    if values.size <= lo:
        # Both middle positions fall among the implicit zeros whatever the
        # signs of the stored values: nothing to look at.
        return 0.0
    implicit = n_bins - int(values.size)
    below = int(np.count_nonzero(values < 0))
    # Dense positions from ``below`` up to here hold a zero.
    through_zero = int(np.count_nonzero(values <= 0)) + implicit

    def slot(position: int) -> int:
        """Where the dense row's element sits among the sorted stored values."""
        if position < below:
            return position
        if position < through_zero:
            return -1
        return position - implicit

    slots = (slot(lo), slot(hi))
    wanted = sorted({s for s in slots if s >= 0})
    if not wanted:
        return 0.0
    selected = np.partition(values, wanted)
    first, second = (0.0 if s < 0 else float(selected[s]) for s in slots)
    if lo == hi:
        return first
    # Halved before the sum, so two large values cannot overflow; equal to
    # ``(a + b) / 2``.
    return 0.5 * first + 0.5 * second


class RowStatistics:
    """The four row statistics that need the row's values, one entry per row.

    ``tic`` is not held here: it is the converter's own pass-1 total, the
    number the TIC image carries (see :meth:`columns`).

    Args:
        n_rows: Rows of the table.
        n_bins: Columns of the table (the length of the mass axis).
    """

    def __init__(self, n_rows: int, n_bins: int) -> None:
        """Allocate the per-row results; the scratch buffer grows on demand."""
        self.n_rows = int(n_rows)
        self.n_bins = int(n_bins)
        self.sum_sq = np.zeros(self.n_rows, dtype=np.float64)
        self.nnz = np.zeros(self.n_rows, dtype=np.uint32)
        self.max_intensity = np.zeros(self.n_rows, dtype=np.float64)
        self.median_all = np.zeros(self.n_rows, dtype=np.float64)
        # The squares of one row. Grown to the largest row seen, never
        # sized by the bin count up front: a row holds a few hundred or a
        # few thousand entries, the axis can hold millions of bins.
        self._scratch = np.empty(0, dtype=np.float64)

    def record(self, row: int, values: Any) -> None:
        """Take the statistics of one row from its stored values.

        Args:
            row: Row of the table.
            values: The row's stored values (zeros allowed, any order).
        """
        values = np.asarray(values, dtype=np.float64)
        k = int(values.size)
        if k == 0:
            self.sum_sq[row] = 0.0
            self.nnz[row] = 0
            self.max_intensity[row] = 0.0
            self.median_all[row] = 0.0
            return
        if self._scratch.size < k:
            self._scratch = np.empty(k, dtype=np.float64)
        squares = self._scratch[:k]
        # Not ``np.dot``: measured 100x slower at 20,000 entries on a
        # laptop, because the BLAS call spins up its threads for each row.
        # The buffer keeps numpy's pairwise sum and allocates nothing.
        np.multiply(values, values, out=squares)
        self.sum_sq[row] = squares.sum()
        self.nnz[row] = np.count_nonzero(values)
        top = float(values.max())
        # A bin the row does not store is a zero, so it is the maximum
        # unless every bin is stored.
        self.max_intensity[row] = top if (top > 0.0 or k >= self.n_bins) else 0.0
        self.median_all[row] = median_over_all_bins(values, self.n_bins)

    def recompute_rows(
        self, matrix: Any, rows: NDArray[np.int64], chunk_entries: Optional[int] = None
    ) -> None:
        """Take the statistics of finished rows out of the matrix.

        For the rows whose final values exist only after the assembly
        summed the spectra of a repeated position. The matrix is a CSC
        matrix over memmaps, so it is never turned into a CSR or dense
        copy as a whole: one pass over its row indices counts each row's
        entries, then the rows are taken in groups of at most
        ``chunk_entries`` entries.

        Args:
            matrix: The finished table matrix (``scipy.sparse.csc_matrix``).
            rows: Row indices to take, ascending.
            chunk_entries: Entries a group may hold; defaults to the
                assembly's sort budget, read at call time so a test can
                shrink it. A row larger than this is taken on its own.
        """
        rows = np.asarray(rows, dtype=np.int64)
        if rows.size == 0:
            return
        if chunk_entries is None:
            chunk_entries = csc_assembly.SORT_CHUNK_ENTRIES
        chunk_entries = max(int(chunk_entries), 1)

        counts = np.zeros(self.n_rows, dtype=np.int64)
        indices = matrix.indices
        n_stored = int(matrix.indptr[-1])
        for start in range(0, n_stored, chunk_entries):
            chunk = np.asarray(indices[start : min(start + chunk_entries, n_stored)])
            counts += np.bincount(chunk, minlength=self.n_rows)

        wanted = counts[rows]
        start = 0
        while start < rows.size:
            end = start + 1
            total = int(wanted[start])
            while end < rows.size and total + int(wanted[end]) <= chunk_entries:
                total += int(wanted[end])
                end += 1
            group = rows[start:end]
            block = matrix[group, :].tocsr()
            for i, row in enumerate(group):
                self.record(int(row), block.data[block.indptr[i] : block.indptr[i + 1]])
            start = end

    def columns(self, tic: NDArray[np.float64]) -> Dict[str, NDArray[Any]]:
        """The five ``obs`` columns, in the order the schema names them.

        Args:
            tic: The rows' totals, taken from the TIC image's numbers.
        """
        return {
            MSI_OBS_TIC_COLUMN: np.asarray(tic, dtype=np.float64),
            MSI_OBS_SUM_SQ_COLUMN: self.sum_sq,
            MSI_OBS_NNZ_COLUMN: self.nnz,
            MSI_OBS_MAX_INTENSITY_COLUMN: self.max_intensity,
            MSI_OBS_MEDIAN_ALL_COLUMN: self.median_all,
        }


def statistics_block(column_sum: NDArray[np.float64]) -> Dict[str, Any]:
    """The ``uns["spectrum_statistics"]`` block of one table.

    Args:
        column_sum: The sum of each column over the rows written, one per
            bin. A column whose sum overflowed keeps the non-finite value.

    Returns:
        ``version``, ``columns`` (one JSON string, like every list of
        strings in ``uns``), ``column_sum`` and ``computed_by``.
    """
    from thyra import __version__

    return {
        "version": MSI_SPECTRUM_STATISTICS_VERSION,
        "columns": json.dumps(list(MSI_OBS_SPECTRUM_STATISTICS_COLUMNS)),
        "column_sum": np.asarray(column_sum, dtype=np.float64),
        "computed_by": f"thyra/{__version__}",
    }
