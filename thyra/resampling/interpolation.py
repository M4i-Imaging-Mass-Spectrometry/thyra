# thyra/resampling/interpolation.py

"""Sparse TIC-preserving resampling.

Evaluates linear interpolation onto a common mass axis only where it can be
non-zero, then rescales to the preserved total ion current -- see
``thyra/resampling/tic.py`` for the rule that decides what "preserved" means
and why the rescale is necessary at all.

Interpolating at a bin centre is sampling, not binning. It is right while a
bin is no wider than the spacing of the source samples, and wrong once a bin
holds several of them: a peak narrower than the bin that falls between two
centres is mostly lost, and the rescale hands its ion current to whatever
the centres did land on (issue #442, design decision D44). So a spectrum
in which some bin holds a sample its centre does not read is rebinned
wherever a bin holds two or more samples -- :func:`wide_bins`,
:func:`rebin_wide_bins` -- and any other spectrum is resampled exactly as
it always was.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from .binning import BinLayout, kept_mz_range
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


def wide_bins(
    layout: BinLayout,
    mzs: NDArray[np.float64],
    intensities: NDArray[np.float64],
) -> NDArray[np.int_]:
    """The bins to rebin rather than interpolate, ascending; often none.

    Interpolating at a bin's centre reads the two samples either side of
    it. A bin that holds any other sample with intensity -- three samples,
    or two on the same side of its centre -- loses that sample's ion
    current to the rescale: that is the defect of issue #442, and a
    spectrum in which it never happens is left to the interpolating
    operator, unchanged. A sample exactly on a centre is read alone, so
    every other sample in its bin counts as skipped. A skipped zero loses
    nothing, so a repeated or
    closely spaced explicit zero does not count; two samples at the same
    m/z count when either carries intensity, since interpolation reads
    only one of them. On an
    evenly sampled source it cannot happen until bins are twice the sample
    spacing, so the default Waters profile route, at 1.14 bins per sample,
    and the default rapifleX one, at 0.2, never come here.

    A spectrum in which it does happen is rebinned wherever a bin holds two
    or more samples, not only in the bins that skip one: between one and
    two samples per bin, interpolating a bin and rebinning its neighbour
    puts the two on different footings, and rebinning both was measured
    the more accurate of the two (design decision D44).

    A sample on an inner edge belongs to the bin on its right; one on the
    top edge belongs to the last bin (:meth:`BinLayout.locate`). Samples
    outside the outer edges are nobody's. Zero-intensity samples count
    towards the bins rebinned: an explicit zero is a measurement.

    Two samples on the same side of a centre are closer together than half
    the gap between that centre and its neighbour, so a spectrum whose
    every spacing is wider than half the gaps near it is settled without
    locating a single sample: one pass over its spacings and a bound read
    once per 64 of them, against the per-sample search the rest costs.

    Args:
        layout: The axis's bins.
        mzs: Source m/z values, ascending.
        intensities: Source intensities, parallel to ``mzs``.

    Returns:
        Sorted, unique bin indices: every bin holding two or more samples,
        or none at all.
    """
    none = np.array([], dtype=np.int_)
    if mzs.size < 2:
        return none
    spacing = np.diff(mzs)
    # A hair of margin: the half-gap of a bin is a midpoint minus a centre,
    # so it can come out an ulp off half the gap.
    margin = 1.0 + 1e-9
    if not float(spacing.min()) < 0.5 * layout.widest * margin:
        return none
    # The bound never falls with m/z, so the one at the top sample of a run
    # of pairs holds for every pair in the run.
    starts = np.arange(0, spacing.size, _RUN)
    tops = np.minimum(starts + _RUN, spacing.size)
    closest = np.minimum.reduceat(spacing, starts)
    if not (closest < 0.5 * layout.gap_near(mzs[tops]) * margin).any():
        return none

    first = int(np.searchsorted(mzs, layout.lower, side="left"))
    stop = int(np.searchsorted(mzs, layout.upper, side="right"))
    if stop - first < 2:
        return none
    inside = mzs[first:stop]
    held = intensities[first:stop] != 0
    bins = layout.locate(inside)
    # ``mzs`` is ascending, so ``bins`` and the (bin, side of its centre)
    # pairs are too: a repeat is a neighbour. Of two neighbours on one
    # side, interpolation reads the one nearer the centre and skips the
    # other: the lower below the centre, the upper at or above it. A
    # sample exactly on the centre is read alone, so the one just below
    # it in the same bin is skipped too.
    centres = layout.axis[bins]
    above = inside >= centres
    side = 2 * bins + above
    same = side[1:] == side[:-1]
    skipped_held = np.where(above[1:], held[1:], held[:-1])
    tied = (inside[1:] == inside[:-1]) & (held[1:] | held[:-1])
    on_centre = (bins[1:] == bins[:-1]) & (inside[1:] == centres[1:]) & held[:-1]
    if not (same & (skipped_held | tied) | on_centre).any():
        return none
    repeated = bins[1:] == bins[:-1]
    return np.unique(bins[1:][repeated]).astype(np.int_)


#: Pairs of neighbouring samples :func:`wide_bins` bounds with one value.
_RUN = 64


def sample_cells(
    mzs: NDArray[np.float64],
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The m/z interval each source sample stands for.

    A profile sample is the signal its detector cell collected, so it
    stands for the stretch of m/z around it: from the midpoint with the
    sample before to the midpoint with the sample after. The cells tile the
    trace, so sharing every sample out by cell conserves its total exactly.

    Where the spacing jumps by more than a factor of two -- the edge of a
    gap in a thresholded or processed spectrum -- the side facing the gap
    is cut back to the width of the other side. Without that, a sample
    beside a 10 Da gap would spread its intensity 5 Da into a region
    nothing was measured in. A smooth grid never trips it, so on a
    detector grid the cells still tile.

    Args:
        mzs: Source m/z values, ascending, at least two.

    Returns:
        ``(lo, hi)``, the cell bounds, parallel to ``mzs``.
    """
    mid = 0.5 * (mzs[1:] + mzs[:-1])
    lo = np.empty_like(mzs)
    hi = np.empty_like(mzs)
    lo[1:] = mid
    hi[:-1] = mid
    # The end samples mirror their one neighbour.
    lo[0] = mzs[0] - (mid[0] - mzs[0])
    hi[-1] = mzs[-1] + (mzs[-1] - mid[-1])
    left = mzs - lo
    right = hi - mzs
    lo = np.where(left > 2.0 * right, mzs - right, lo)
    hi = np.where(right > 2.0 * left, mzs + left, hi)
    return lo, hi


def rebin_wide_bins(
    layout: BinLayout,
    wide: NDArray[np.int_],
    mzs: NDArray[np.float64],
    intensities: NDArray[np.float64],
) -> NDArray[np.float64]:
    """What each wide bin holds: its share of every sample cell it overlaps.

    Each sample's intensity is divided among the bins its cell
    (:func:`sample_cells`) overlaps, in proportion to the overlap. Over a
    bin many samples wide this is summing the samples inside it; at the
    bin's edges, a sample whose cell straddles the edge is split rather
    than given whole to one side. Nothing is sampled, so a peak narrower
    than the bin keeps its ion current wherever it falls.

    Computed as the difference, across each bin, of the running total of
    the cells' intensity from the left, so no (sample, bin) pairs are ever
    listed. The running total only rises, so no bin comes out negative,
    and a bin no cell with intensity overlaps comes out exactly zero.

    Args:
        layout: The axis's bins.
        wide: The bins to fill, ascending, as :func:`wide_bins` returns.
        mzs: Source m/z values, ascending, at least two.
        intensities: Source intensities, parallel to ``mzs``.

    Returns:
        The value of each bin in ``wide``, parallel to it; zeros when no
        sample lies inside the axis.
    """
    lo, hi = sample_cells(mzs)
    # A sample outside the axis is nobody's; one inside keeps only the part
    # of its cell that is.
    inside = slice(
        int(np.searchsorted(mzs, layout.lower, side="left")),
        int(np.searchsorted(mzs, layout.upper, side="right")),
    )
    lo = np.clip(lo[inside], layout.lower, layout.upper)
    hi = np.clip(hi[inside], layout.lower, layout.upper)
    weights = intensities[inside]
    if weights.size == 0:
        return np.zeros(wide.size)
    before = np.concatenate(([0.0], np.cumsum(weights)))

    def total_left_of(edges: NDArray[np.float64]) -> NDArray[np.float64]:
        # Cells starting left of the edge; a zero-width cell (a repeated
        # m/z) on an edge is the right-hand bin's, as a sample on an edge is.
        cell = np.searchsorted(lo, edges, side="left") - 1
        held = np.maximum(cell, 0)
        width = hi[held] - lo[held]
        part = np.divide(
            edges - lo[held], width, out=np.ones_like(width), where=width > 0
        )
        total = before[held] + weights[held] * np.clip(part, 0.0, 1.0)
        total[cell < 0] = 0.0
        return total

    upper = layout.edges_of(wide + 1)
    right = total_left_of(upper)
    # The top edge closes the last bin: a sample on it is that bin's.
    right[wide == layout.axis.size - 1] = before[-1]
    shares: NDArray[np.float64] = right - total_left_of(layout.edges_of(wide))
    return shares


def _merge_wide_bins(
    layout: BinLayout,
    wide: NDArray[np.int_],
    indices: NDArray[np.int_],
    values: NDArray[np.float64],
    mzs: NDArray[np.float64],
    wide_values: NDArray[np.float64],
) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
    """Put the interpolated bins on the wide bins' footing and merge the two.

    A rebinned bin holds ion current: the samples it covers, summed. An
    interpolated value is a height, and the bin it stands for is a fraction
    or a multiple of one sample wide. Scaling it by the bin's width over
    the spacing of the samples around its centre turns it into the ion
    current the bin covers, so a spectrum with both kinds of bin is on one
    scale before the rescale sets its total. Inside a gap the spacing is
    the gap, and the interpolated bridge across it is worth no more than
    the one sample it stands for.

    A spectrum with no wide bin never comes here: its bins are all
    heights, the rescale alone puts them on the TIC, and it is resampled
    exactly as it was before design decision D44.

    Args:
        layout: The axis's bins.
        wide: The rebinned bins, ascending.
        indices: The interpolated bins, ascending, disjoint from ``wide``.
        values: Their interpolated heights.
        mzs: Source m/z values, ascending.
        wide_values: What :func:`rebin_wide_bins` gave each of ``wide``.

    Returns:
        ``(indices, values)`` over both kinds of bin, indices ascending.
    """
    segment = np.clip(
        np.searchsorted(mzs, layout.axis[indices], side="right") - 1,
        0,
        mzs.size - 2,
    )
    spacing = mzs[segment + 1] - mzs[segment]
    width = layout.edges_of(indices + 1) - layout.edges_of(indices)
    values = values * np.divide(
        width, spacing, out=np.zeros_like(width), where=spacing > 0
    )

    merged = np.concatenate((indices, wide))
    order = np.argsort(merged, kind="stable")
    return merged[order], np.concatenate((values, wide_values))[order]


def tic_preserving_sparse(
    axis: NDArray[np.float64],
    mzs: NDArray[np.float64],
    intensities: NDArray[np.float64],
    gap_tolerance_da: Optional[float],
    axis_range: Optional[Tuple[float, float]] = None,
    layout: Optional[BinLayout] = None,
) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
    """TIC-preserving resampling, evaluated only where it can be non-zero.

    This is the operator
    :class:`~thyra.resampling.strategies.tic_preserving.TICPreservingStrategy`
    documents -- interpolate onto the axis, zero unsupported bins, rescale to
    the preserved TIC -- restricted to the axis points
    :func:`tic_support_bins` reports. Every other axis point interpolates
    to exactly zero, so scattering the result into a zero array reproduces
    the dense evaluation bin for bin; the sole difference is the order in
    which the rescale sums its terms.

    Where :func:`wide_bins` finds bins to rebin, they are rebinned
    (:func:`rebin_wide_bins`), the other bins are scaled onto the same
    ion-current footing (:func:`_merge_wide_bins`), and the rescale then
    sets the total as before. A spectrum with none takes exactly the path
    described above. See design decision D44.

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
            this is the same range
            :class:`~thyra.resampling.strategies.nearest_neighbor.NearestNeighborStrategy`
            keeps peaks inside.
        layout: The axis's bins for this range. ``None`` lays them here; a
            caller resampling many spectra onto one axis lays them once.

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

    wide = np.array([], dtype=np.int_)
    if layout is None and axis.size >= 2:
        layout = BinLayout(axis, kept_range)
    if layout is not None:
        wide = wide_bins(layout, mzs_sorted, intensities_sorted)

    indices = tic_support_bins(axis, mzs_sorted, intensities_sorted)
    if wide.size:
        # Interpolate only the bins that are not rebinned.
        indices = indices[~np.isin(indices, wide, assume_unique=True)]
    elif indices.size == 0:
        return empty

    targets = axis[indices]
    values = np.interp(targets, mzs_sorted, intensities_sorted, left=0.0, right=0.0)

    # Discard bins no source point vouches for, before the rescale so the
    # intensity returns to the bins that were measured rather than being
    # deleted. No-op when no tolerance was configured. A wide bin holds
    # samples, so it is vouched for and never reaches this.
    zero_across_gaps(values, targets, mzs_sorted, gap_tolerance_da)

    if wide.size:
        assert layout is not None
        indices, values = _merge_wide_bins(
            layout,
            wide,
            indices,
            values,
            mzs_sorted,
            rebin_wide_bins(layout, wide, mzs_sorted, intensities_sorted),
        )

    # Rescale to the required TIC -- the step that makes the method live up
    # to its name. Reads only the kept range and the sum of ``values``, so
    # the subset evaluation rescales exactly as the dense one would.
    rescale_to_preserved_tic(values, axis, mzs_sorted, intensities_sorted, kept_range)

    keep = values != 0
    return indices[keep], values[keep]
