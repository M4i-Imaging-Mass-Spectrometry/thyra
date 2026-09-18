# thyra/resampling/binning.py

"""Nearest-neighbour binning onto a common mass axis.

This is the operator that both ``BaseSpatialDataConverter`` and the sibling
tables (the mobility heatmap, the mobility grid, and the MS/MS demultiplexer)
share: map each m/z value to the bin of a target axis nearest it, then
accumulate intensities per bin. It lives in ``thyra.resampling`` rather than
in ``thyra.converters`` because the sibling-table sinks need it too, and a
converter-side home would make importing it from those sinks a cycle -- they
used to reach it with a function-local import for exactly that reason.
"""

import logging
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from .types import AxisLinearisation

logger = logging.getLogger(__name__)


def kept_mz_range(
    axis: NDArray[np.float64],
    axis_range: Optional[Tuple[float, float]],
) -> Tuple[float, float]:
    """The m/z range a peak has to be inside to survive resampling.

    Every physics generator lays ``target_bins + 1`` bin *edges* across the
    requested ``[min_mz, max_mz]`` and returns the midpoints, so the first
    and last axis point sit half a bin *inside* the range that was asked
    for -- ``[50.0001, 999.9975]`` for a source declaring 50-1000. Testing
    membership against the axis points therefore discarded peaks the caller
    had asked to keep, and did it worst on the sources that declare their
    range *as* their first and last sample: PHI TOF-SIMS takes
    ``mass_range`` from the first and last detector channel, so both were
    dropped in every pixel (issue #239).

    The rule is the **declared** range, which is the outer bin edges, which
    is at most half a bin beyond the first and last centre. That bound is
    what keeps this from becoming the clamp
    :meth:`~thyra.resampling.strategies.nearest_neighbor.NearestNeighborStrategy.resample`
    documents: a peak further out than the range is still dropped, so
    narrowing the range with ``--resample-min-mz`` cannot pile the
    discarded part of the spectrum onto bin 0.

    A uniform axis is ``np.linspace(min_mz, max_mz, n)``, whose end points
    *are* the declared bounds, so nothing changes there -- nor for
    ``--no-resample``, where no axis was built and ``axis_range`` is None.

    A module-level function rather than a method for the reason
    :func:`nn_map_to_bins` is: the unbound-call test harnesses drive the
    resampling surface on a ``SimpleNamespace`` that has no methods.

    Args:
        axis: The target mass axis, ascending.
        axis_range: The declared ``(min_mz, max_mz)``, or None when no
            resampled axis was built.

    Returns:
        ``(min_mz, max_mz)``, always covering ``axis`` itself.
    """
    if axis_range is None:
        return float(axis[0]), float(axis[-1])
    # Never narrower than the axis: a generator that returned points
    # outside the range it was handed must not cost anyone a peak that has
    # a bin waiting for it.
    return (
        min(float(axis_range[0]), float(axis[0])),
        max(float(axis_range[1]), float(axis[-1])),
    )


#: How far, in bins, a resampled axis may deviate from its own linearisation
#: before the converter stops computing bin indices and searches for them
#: instead. The +-1 repair in :func:`nn_map_to_bins` is exact below 0.5
#: (see :class:`~thyra.resampling.types.AxisLinearisation`); a quarter of a
#: bin leaves float rounding in the position nowhere to matter. Real axes sit
#: below 0.005: the deviation is second order in one bin's relative width.
NN_LINEARISATION_MARGIN = 0.25


def usable_linearisation(
    linearisation: Optional[AxisLinearisation],
    axis: NDArray[np.float64],
    min_gap: float,
) -> Optional[AxisLinearisation]:
    """The linearisation the nearest-neighbour mapping may compute indices from, or None.

    Checked once, when the axis is built. The closed form is exact only on
    a strictly ascending axis (a duplicated value must still map to its
    first occurrence, which is what ``np.searchsorted`` does and the
    repair does not) whose every centre lies within
    :data:`NN_LINEARISATION_MARGIN` bins of ``u0 + i * du``. Any axis a
    generator lays satisfies both by construction; the check is what turns
    that argument into something the conversion has verified rather than
    trusted, and it is what keeps ``--no-resample`` -- a raw union axis
    with no law at all -- on the search.
    """
    if linearisation is None:
        return None
    if not (min_gap > 0):
        logger.info(
            "Nearest-neighbour bin index: axis is not strictly ascending, "
            "keeping the binary search"
        )
        return None
    deviation = linearisation.deviation(axis)
    if not (deviation < NN_LINEARISATION_MARGIN):
        logger.warning(
            "Nearest-neighbour bin index: axis deviates from its linearisation "
            "by %.3f bins (limit %.2f), keeping the binary search",
            deviation,
            NN_LINEARISATION_MARGIN,
        )
        return None
    logger.info(
        "Nearest-neighbour bin index: closed form, worst deviation %.2e bins",
        deviation,
    )
    return linearisation


def nn_map_to_bins(
    axis: NDArray[np.float64],
    mzs: NDArray[np.float64],
    linearisation: Optional[AxisLinearisation] = None,
) -> NDArray[np.int_]:
    """Map in-range m/z values to their nearest bin on ``axis``.

    Ties (a peak exactly between two bins) resolve to the right bin,
    matching the strict ``<`` comparison this code has always used.

    With a ``linearisation`` -- the coordinate the generator laid the axis
    uniformly in, checked by :func:`usable_linearisation` -- the index is
    computed: round the value's position in that coordinate, then take
    the nearest of that bin and its two neighbours, comparing in m/z with
    the same tie rule. The rounded position is within one bin of the
    answer whenever the axis deviates from its linearisation by under
    half a bin, which every generated axis does (design decision D21), so
    the two routes return the same index for every value. They differ
    only in cost, which is the point: ``np.searchsorted`` on an axis of
    10^5 to 10^6 bins is a cache-missing binary search per peak, where
    this is three linear reads.

    Without one, the bin is found by binary search. That route stays for
    any axis without a law: the raw union axis under ``--no-resample``,
    an axis a caller passes bare, and the sibling tables' sinks.

    A module-level function rather than a method so the unbound-call test
    harnesses (a ``SimpleNamespace`` posing as the converter) keep working.

    Args:
        axis: The target mass axis, ascending.
        mzs: m/z values, all within the range :func:`kept_mz_range`
            reports -- so within half a bin of ``axis``, not necessarily
            within ``[axis[0], axis[-1]]``.
        linearisation: The axis's own, or None to search.

    Returns:
        The nearest-bin index of each m/z value, same length as ``mzs``.
    """
    if linearisation is not None:
        return _nn_computed_bins(axis, mzs, linearisation)

    # Find insertion points using vectorized binary search
    indices = np.searchsorted(axis, mzs)

    # Clip to valid range. Everything reaching here is inside the declared
    # range, so this pins searchsorted's one-past-the-end result and sends
    # a peak in the half-bin skirt of either end into the edge bin it
    # belongs to; it can no longer pull in a peak from outside the range.
    indices_clipped = np.clip(indices, 0, len(axis) - 1)

    # For non-boundary points, check if left is closer
    # Only check where we're not at the left edge
    check_left = indices > 0
    if np.any(check_left):
        # Get distances only for points that need checking
        mz_values = axis[indices_clipped[check_left]]
        mz_values_left = axis[indices_clipped[check_left] - 1]
        mz_query = mzs[check_left]

        # Use left if it's closer
        use_left = np.abs(mz_values_left - mz_query) < np.abs(mz_values - mz_query)
        indices_clipped[check_left] = np.where(
            use_left, indices_clipped[check_left] - 1, indices_clipped[check_left]
        )
    return indices_clipped


def _nn_computed_bins(
    axis: NDArray[np.float64],
    mzs: NDArray[np.float64],
    linearisation: AxisLinearisation,
) -> NDArray[np.int_]:
    """The closed-form route of :func:`nn_map_to_bins`.

    Why the rounded position is never more than one bin off. A value's
    true nearest bin ``j`` has the value between the centres ``j - 1``
    and ``j + 1``; ``forward`` is monotone, so its position lies between
    those two centres' positions, each within ``d < 1/2`` of its own
    index; so the position lies in ``(j - 3/2, j + 3/2)`` and rounds to
    ``j - 1``, ``j`` or ``j + 1``. The comparison among those three is
    then the reference's own comparison, so ties resolve as it resolves
    them: to the right.

    A value in the half-bin skirt beyond either end rounds to ``-1`` or
    ``n``; the clip sends it to the edge bin, whose neighbour is then
    compared like any other -- the same answer the search's clip gives.
    """
    n = axis.size
    k = np.rint(linearisation.positions(mzs)).astype(np.intp)
    np.clip(k, 0, n - 1, out=k)
    left = np.maximum(k - 1, 0)
    right = np.minimum(k + 1, n - 1)
    d_left = np.abs(axis[left] - mzs)
    d_k = np.abs(axis[k] - mzs)
    d_right = np.abs(axis[right] - mzs)
    # Start from the right neighbour and move left only on a strict
    # improvement, so an exact tie keeps the right-hand bin.
    best = right
    best_d = d_right
    closer = d_k < best_d
    best = np.where(closer, k, best)
    best_d = np.where(closer, d_k, best_d)
    closer = d_left < best_d
    return np.where(closer, left, best)


def nn_accumulate(
    idx: NDArray[np.int_], intensities: NDArray[np.float64]
) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
    """Sum intensities per bin and return the non-zero bins, ascending.

    Two equivalent routes. When ``idx`` is non-decreasing -- true whenever
    the spectrum's m/z values are ascending, which every format seen so
    far produces -- equal bins form contiguous runs, so the per-bin sums
    are ``np.add.reduceat`` over the run starts. That is O(n_peaks),
    where the ``np.bincount`` fallback is O(n_bins): for a centroid
    spectrum of hundreds of peaks against an axis of 10^5 bins the
    difference is roughly 7x per spectrum. Both sum left-to-right over
    the same float64 values, so the results are bit-identical; the
    fallback keeps unsorted input correct.

    The kept-bin test is ``!= 0`` rather than ``> 0`` because a bin whose
    accumulated value is negative is still a measurement: dropping it
    silently raises the stored TIC above the input's. Baseline-subtracted
    data can carry negative intensities, though every export seen so far
    filters them out upstream.
    """
    # bincount always promotes its weights to float64; match it so both
    # accumulation routes return the same dtype and the same rounding.
    vals = intensities.astype(np.float64, copy=False)

    if idx.size and bool(np.all(idx[1:] >= idx[:-1])):
        starts = np.concatenate(([0], np.flatnonzero(np.diff(idx)) + 1))
        if starts.size == idx.size:
            # Every peak already sits in its own bin (the common case
            # for both centroid and profile data): nothing to sum.
            bins = idx
            sums = vals
        else:
            bins = idx[starts]
            sums = np.add.reduceat(vals, starts)
        keep = sums != 0
        return bins[keep].astype(np.int_, copy=False), sums[keep]

    accumulated = np.bincount(idx, weights=vals)
    nonzero_mask = accumulated != 0
    nonzero_indices = np.where(nonzero_mask)[0].astype(np.int_)
    nonzero_values = accumulated[nonzero_mask]
    return nonzero_indices, nonzero_values.astype(np.float64)


class SharedAxisNNCache:
    """Precomputed nearest-neighbor mapping for one recurring m/z array.

    Built the first time
    :class:`~thyra.resampling.strategies.nearest_neighbor.NearestNeighborStrategy`
    sees a spectrum, and reused for every later spectrum that carries the
    same m/z array.
    ``key`` is a private copy so a caller-side mutation of the original
    array cannot fool the equality check; ``key_ref`` keeps the original
    object for the O(1) identity test that readers yielding one shared
    array (Rapiflex, Waters, PHI) hit every time.
    """

    __slots__ = (
        "key",
        "key_ref",
        "n_total",
        "n_dropped",
        "lo",
        "hi",
        "starts",
        "bins",
    )

    key: NDArray[np.float64]
    key_ref: NDArray[np.float64]
    n_total: int
    n_dropped: int
    lo: int
    hi: int
    starts: Optional[NDArray[np.int_]]
    bins: NDArray[np.int_]
