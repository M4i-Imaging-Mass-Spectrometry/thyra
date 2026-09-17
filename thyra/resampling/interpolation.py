# thyra/resampling/interpolation.py

"""Sparse TIC-preserving interpolation.

Evaluates linear interpolation onto a common mass axis only where it can be
non-zero, then rescales to the preserved total ion current -- see
``thyra/resampling/tic.py`` for the rule that decides what "preserved" means
and why the rescale is necessary at all.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from .binning import kept_mz_range
from .gaps import zero_across_gaps
from .tic import preserved_tic, rescale_to_preserved_tic


def tic_support_bins(
    axis: NDArray[np.float64],
    mzs: NDArray[np.float64],
    intensities: NDArray[np.float64],
) -> NDArray[np.int_]:
    """Axis indices at which the linear interpolant of a spectrum can be non-zero.

    ``np.interp`` draws straight lines between consecutive source points, so
    the interpolant is non-zero only on segments with a non-zero endpoint.
    Consecutive non-zero samples form runs; each run's support is the open
    interval from the sample before it to the sample after it (those two are
    where the trace touches zero), closed at the spectrum's own ends where
    there is no such neighbour. Runs are separated by at least one zero
    sample, so their supports never overlap and the indices come out sorted.

    A zero-suppressed profile -- Waters MassLynx stores samples in clusters
    around each peak with an explicit zero at either edge -- has supports
    covering a small fraction of a fine axis, which is what makes the sparse
    evaluation cheap. A spectrum with no zeros in it is one run, and its
    support is every axis point the dense evaluation would have populated.

    Args:
        axis: Target mass axis, ascending.
        mzs: Source m/z values, ascending.
        intensities: Source intensities, parallel to ``mzs``.

    Returns:
        Sorted, unique axis indices. Empty when nothing is non-zero.
    """
    nonzero = intensities != 0
    if not nonzero.any():
        return np.array([], dtype=np.int_)

    last = mzs.size - 1
    edges = np.diff(np.concatenate(([False], nonzero, [False])).astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1) - 1

    # Open at a neighbouring zero sample (the interpolant is exactly zero
    # there), closed at the spectrum's own first or last sample.
    lo = np.where(
        starts == 0,
        np.searchsorted(axis, mzs[0], side="left"),
        np.searchsorted(axis, mzs[np.maximum(starts - 1, 0)], side="right"),
    )
    hi = np.where(
        ends == last,
        np.searchsorted(axis, mzs[last], side="right"),
        np.searchsorted(axis, mzs[np.minimum(ends + 1, last)], side="left"),
    )

    lengths = hi - lo
    keep = lengths > 0
    if not keep.any():
        return np.array([], dtype=np.int_)
    lo = lo[keep]
    lengths = lengths[keep]
    offsets = np.cumsum(lengths) - lengths
    total = int(lengths.sum())
    return np.repeat(lo - offsets, lengths) + np.arange(total, dtype=np.int_)


def tic_preserving_sparse(
    axis: NDArray[np.float64],
    mzs: NDArray[np.float64],
    intensities: NDArray[np.float64],
    gap_tolerance_da: Optional[float],
    axis_range: Optional[Tuple[float, float]] = None,
) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
    """TIC-preserving resampling, evaluated only where it can be non-zero.

    This is the operator ``BaseSpatialDataConverter._tic_preserving_resample``
    documents -- interpolate onto the axis, zero unsupported bins, rescale to
    the preserved TIC -- restricted to the axis points
    :func:`tic_support_bins` reports. Every other axis point interpolates
    to exactly zero, so scattering the result into a zero array reproduces
    the dense evaluation bin for bin; the sole difference is the order in
    which the rescale sums its terms.

    Measured on a Waters SELECT SERIES MRT run (13,398 pixels, ~15,000
    stored samples per strong pixel, 1.05M-bin axis): the dense form cost
    570 s against 16 s for nearest-neighbour binning, almost all of it in
    interpolating onto and then scanning a million bins per pixel of which
    ~13,000 were ever non-zero.

    Args:
        axis: Target mass axis, ascending.
        mzs: Source m/z values, any order.
        intensities: Source intensities, parallel to ``mzs``.
        gap_tolerance_da: See :func:`thyra.resampling.gaps.zero_across_gaps`.
        axis_range: The declared ``(min_mz, max_mz)`` the axis was built
            across, which is the share of the spectrum the result is
            entitled to carry. ``None`` falls back to the axis's own span.
            See :func:`thyra.resampling.binning.kept_mz_range`: the two
            resampling methods have to agree on what the axis covers, so
            this is the same range ``_nearest_neighbor_resample`` keeps
            peaks inside.

    Returns:
        ``(bin_indices, intensities)`` holding only the non-zero bins,
        indices ascending.
    """
    empty = (np.array([], dtype=np.int_), np.array([], dtype=np.float64))
    if mzs.size == 0:
        return empty

    kept_range = kept_mz_range(axis, axis_range)

    if np.all(mzs[:-1] <= mzs[1:]):
        mzs_sorted = mzs
        intensities_sorted = intensities
    else:
        order = np.argsort(mzs)
        mzs_sorted = mzs[order]
        intensities_sorted = intensities[order]

    if mzs_sorted.size == 1:
        # np.interp cannot interpolate a lone point onto a grid that does
        # not contain it -- it would return all zeros and lose the peak.
        # Place it in its nearest bin, as the nearest-neighbour path does.
        target_tic = preserved_tic(
            mzs_sorted, intensities_sorted, kept_range[0], kept_range[1]
        )
        if target_tic <= 0.0:
            return empty
        nearest = int(np.argmin(np.abs(axis - mzs_sorted[0])))
        return (
            np.array([nearest], dtype=np.int_),
            np.array([target_tic], dtype=np.float64),
        )

    indices = tic_support_bins(axis, mzs_sorted, intensities_sorted)
    if indices.size == 0:
        return empty

    targets = axis[indices]
    values = np.interp(targets, mzs_sorted, intensities_sorted, left=0.0, right=0.0)

    # Discard bins no source point vouches for, before the rescale so the
    # intensity returns to the bins that were measured rather than being
    # deleted. No-op when no tolerance was configured.
    zero_across_gaps(values, targets, mzs_sorted, gap_tolerance_da)

    # Rescale to the required TIC -- the step that makes the method live up
    # to its name. Reads only the kept range and the sum of ``values``, so
    # the subset evaluation rescales exactly as the dense one would.
    rescale_to_preserved_tic(values, axis, mzs_sorted, intensities_sorted, kept_range)

    keep = values != 0
    return indices[keep], values[keep]
