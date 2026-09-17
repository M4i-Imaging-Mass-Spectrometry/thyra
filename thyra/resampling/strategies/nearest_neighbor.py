# thyra/resampling/strategies/nearest_neighbor.py

"""The operator ``ResamplingMethod.NEAREST_NEIGHBOR`` selects.

Each source peak is accumulated into the one target bin nearest it, so the
summed intensity of a spectrum is preserved exactly and a peak stays where
it was measured. Peaks outside the declared range are dropped rather than
folded into the edge bins, and the mapping itself comes from
:mod:`thyra.resampling.binning`, shared with the sibling tables so every
placement onto one axis is made by one function.

A note on the name, because this class carried it before and did the
opposite with it. Until PR #366 ``NearestNeighborStrategy`` INTERPOLATED:
each target point took the intensity of the nearest source point, so one
source peak was replicated across every target point closer to it than to
any other and the summed intensity scaled with how densely the target axis
was laid. That class had no callers, and the conversion has always BINNED,
by way of the private converter methods this class is built from. The
class now bins. The behaviour of a conversion is unchanged; only the
spelling of where it lives has moved.
"""

from typing import Any, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..binning import SharedAxisNNCache, nn_accumulate, nn_map_to_bins
from ..types import AxisLinearisation
from .base import ResamplingStrategy


class NearestNeighborStrategy(ResamplingStrategy):
    """Bin each peak into its nearest bin on the target axis.

    Args:
        axis: The target mass axis, ascending.
        axis_range: The declared ``(min_mz, max_mz)``, or ``None``.
        linearisation: The axis's own coordinate, already checked by
            :func:`thyra.resampling.binning.usable_linearisation`, or
            ``None`` to place every peak by binary search. The check is
            the caller's because it belongs to the axis and is made once,
            when the axis is built (design decision D21).
    """

    def __init__(
        self,
        axis: NDArray[np.float64],
        axis_range: Optional[Tuple[float, float]],
        linearisation: Optional[AxisLinearisation] = None,
    ) -> None:
        """Bind the operator to one axis, its range and its linearisation."""
        super().__init__(axis, axis_range)
        self.linearisation = linearisation
        # Shared-axis nearest-neighbor cache. Continuous imzML, Rapiflex,
        # Waters and PHI all hand every spectrum the same m/z array, so the
        # peak-to-bin mapping is computed once and verified per spectrum by
        # array equality instead of being re-derived by searchsorted. None
        # means "not built yet"; False means "tried and the data does not
        # share an axis, stop checking". See resample().
        self._shared_cache: Any = None
        self._cache_misses: int = 0

    def map_to_bins(self, mzs: NDArray[np.float64]) -> NDArray[np.int_]:
        """The bin of this axis nearest each m/z value, in the spectrum's order."""
        return nn_map_to_bins(self.axis, mzs, self.linearisation)

    def resample(
        self, mzs: NDArray[np.float64], intensities: NDArray[np.float64]
    ) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
        """Resample spectrum using nearest neighbor binning.

        Maps each m/z value to its nearest bin in the common mass axis and
        accumulates intensities. Returns only non-zero bins for efficiency.

        Peaks outside the axis are **dropped**, not folded into the edge
        bins. They used to be clipped to bin 0 or the last bin and then
        accumulated there, so narrowing the mass range -- the most ordinary
        thing ``--resample-min-mz`` / ``--resample-max-mz`` are for --
        piled everything below the floor onto the first bin and everything
        above the ceiling onto the last. On real ``pea.imzML`` resampled to
        400-800 m/z, bin 0 held 654,158 counts where a real peak there is
        around 80. The total was conserved exactly, so a TIC check could
        not see it; the peak was simply in the wrong place, 1,634x the
        median interior bin.

        "In range" is the **declared** ``[min_mz, max_mz]``, which is the
        outer bin edges and so at most half a bin beyond the first and last
        centre -- see :func:`thyra.resampling.binning.kept_mz_range` for
        why the axis points themselves are the wrong test and why the rule
        stops there. A peak further out is dropped, not clamped.
        :class:`~thyra.resampling.strategies.tic_preserving.TICPreservingStrategy`
        follows the same range, so the two methods agree on what the axis
        covers.

        Args:
            mzs: Original m/z values from spectrum
            intensities: Corresponding intensity values

        Returns:
            Tuple of (bin_indices, accumulated_intensities) containing only non-zero bins
        """
        if mzs.size == 0:
            return np.array([], dtype=np.int_), np.array([], dtype=np.float64)

        # Shared-axis fast path: when every spectrum carries the same m/z
        # array -- continuous imzML; every other reader, Waters and PHI
        # included, reports has_shared_mass_axis = False -- the peak-to-bin
        # mapping is a property of the axis pair, not of the spectrum. It is
        # computed once; each later spectrum only proves its m/z array is
        # the same one -- an exact array comparison, which is far cheaper
        # than re-deriving the mapping by binary search per spectrum.
        # ``False`` means the cache disabled itself (processed-mode data).
        if self._shared_cache is not False:
            result = self._resample_via_cache(mzs, intensities)
            if result is not None:
                return result

        lo_mz, hi_mz = self.kept_range
        in_range = (mzs >= lo_mz) & (mzs <= hi_mz)
        if not in_range.all():
            self._count_out_of_range(int(mzs.size - in_range.sum()), int(mzs.size))
            mzs = mzs[in_range]
            intensities = intensities[in_range]
            if mzs.size == 0:
                return np.array([], dtype=np.int_), np.array([], dtype=np.float64)

        return nn_accumulate(self.map_to_bins(mzs), intensities)

    def _build_shared_cache(
        self, mzs: NDArray[np.float64]
    ) -> Optional[SharedAxisNNCache]:
        """Precompute the nearest-neighbor mapping for one m/z array.

        Returns None when the array cannot be cached -- empty, or not
        ascending, which the contiguous in-range slice below relies on.
        The mapping itself comes from the same :meth:`map_to_bins` the
        generic path uses, so a cache hit and a fresh computation agree
        bin for bin.
        """
        if mzs.size == 0 or not bool(np.all(mzs[1:] >= mzs[:-1])):
            return None

        # Ascending m/z makes the in-range subset one contiguous slice,
        # with the same inclusive endpoints as the generic path's mask --
        # the declared range, not the axis's own span (kept_mz_range).
        lo_mz, hi_mz = self.kept_range
        lo = int(np.searchsorted(mzs, lo_mz, side="left"))
        hi = int(np.searchsorted(mzs, hi_mz, side="right"))

        cache = SharedAxisNNCache()
        cache.key = mzs.copy()
        cache.key_ref = mzs
        cache.n_total = int(mzs.size)
        cache.n_dropped = int(mzs.size - (hi - lo))
        cache.lo = lo
        cache.hi = hi
        if hi <= lo:
            cache.starts = None
            cache.bins = np.array([], dtype=np.int_)
            return cache

        idx = self.map_to_bins(mzs[lo:hi])
        # Ascending m/z onto an ascending axis gives non-decreasing bins,
        # so equal bins form contiguous runs.
        starts = np.concatenate(([0], np.flatnonzero(np.diff(idx)) + 1))
        if starts.size == idx.size:
            # Every peak in its own bin: per-spectrum work reduces to a
            # zero-filter over the raw intensities.
            cache.starts = None
            cache.bins = idx.astype(np.int_, copy=False)
        else:
            cache.starts = starts
            cache.bins = idx[starts].astype(np.int_, copy=False)
        return cache

    def _resample_via_cache(
        self,
        mzs: NDArray[np.float64],
        intensities: NDArray[np.float64],
    ) -> Optional[Tuple[NDArray[np.int_], NDArray[np.float64]]]:
        """Resample through the shared-axis cache, or return None to decline.

        The cache is built from the first spectrum seen. A hit requires the
        spectrum's m/z array to be the cached one -- same object, or equal
        element for element -- so a lying reader cannot get a stale mapping;
        it can only miss. After five consecutive misses the cache disables
        itself so processed-mode data stops paying for the comparison
        (its size check is O(1) in the common case anyway).
        """
        cache = self._shared_cache
        if cache is False:
            return None
        if cache is None:
            cache = self._build_shared_cache(mzs)
            if cache is None:
                self._shared_cache = False
                return None
            self._shared_cache = cache

        if not (
            mzs is cache.key_ref
            or (mzs.size == cache.n_total and bool(np.array_equal(mzs, cache.key)))
        ):
            self._cache_misses += 1
            if self._cache_misses > 4:
                self._shared_cache = False
            return None

        if cache.n_dropped:
            self._count_out_of_range(cache.n_dropped, cache.n_total)
        if cache.hi <= cache.lo:
            return np.array([], dtype=np.int_), np.array([], dtype=np.float64)

        # Match bincount's float64 promotion so cached and generic results
        # carry the same rounding.
        vals = intensities[cache.lo : cache.hi].astype(np.float64, copy=False)
        if cache.starts is None:
            sums = vals
        else:
            sums = np.add.reduceat(vals, cache.starts)
        keep = sums != 0
        return cache.bins[keep], sums[keep]
