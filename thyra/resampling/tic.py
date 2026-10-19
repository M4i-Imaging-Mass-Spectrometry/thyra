"""Total ion current preservation: the rule every resampling path shares.

Thyra interpolates spectra onto a common mass axis in one place,
``tic_preserving_sparse`` in ``thyra/resampling/interpolation.py``, which
``TICPreservingStrategy`` calls and the conversion reaches only through
that class; it used to be two (an inline copy on the converter and a
public strategy that shared nothing with it), they drifted, the converter
lost its rescaling step, and keeping the rule in this module is what
stopped that happening again and what the one remaining strategy is built
on (issue #277). On a spectrum whose bins are wide enough to skip a
sample, that function rebins instead of interpolating (design decision
D50); the rule below then sets the total exactly as it does for
interpolated bins.

Why rescaling is needed at all
------------------------------

Linear interpolation *samples* the spectrum at every target point, so the
sum of the interpolated values scales with how densely the target axis is
laid out. Interpolating 4,000 source points onto a 190,000-bin axis returns
roughly 47x the input TIC, and a 150-peak centroid spectrum over 1000x. The
sum only comes out right when the target axis happens to match the source
point density. Rescaling to a TIC computed from the *source* removes that
dependence on axis width.

Why cropping has to reduce the preserved TIC
--------------------------------------------

Interpolation drops everything outside the axis range. Rescaling to the full
input TIC would push that dropped intensity back into the bins that survive,
so a window cropped with ``--resample-min-mz`` / ``--resample-max-mz`` would
claim every ion from the parts that were cut away. The preserved total is
therefore the share of the spectrum that lies inside the axis range.

Why that share is the samples inside the range
----------------------------------------------

The share is the sum of the source samples inside the range, bounds
included: the samples ``nearest_neighbor`` keeps, so the two methods store
the same total for one pixel. It used to be measured by area, the integral
over the range against the integral over the whole spectrum. A sum and a
total are in one unit; an area share is not, once the samples are unevenly
spaced. A spectrum with gaps -- left by ``--intensity-threshold``, or by a
profile file that stores only its peaks -- had straight lines drawn across
each gap, so the share followed the spacing of the samples rather than the
ions: on public files the stored in-window current came out 0.26 to 4.9
times what was measured there, and 126 times on a constructed one (design
decision D71, extending D13).

The area share stays for one case: a window that holds no sample at all.
A 0.6 Da window over a spectrum sampled every 1 Da contains no source
point, and a sum would report zero and blank a spectrum interpolation
represents perfectly well. There the integral keeps a proportionally small
share rather than falling off a cliff.

The common case, where the axis spans the whole spectrum, is short-circuited
so it returns the input TIC exactly rather than to within rounding.
"""

from typing import Any, Optional, Tuple, TypeVar

import numpy as np
import numpy.typing as npt

# ``rescale_to_preserved_tic`` returns the very array it was handed, so the
# element type it was called with is the element type it gives back.
_FloatT = TypeVar("_FloatT", bound=np.floating[Any])


def preserved_tic(
    mzs: npt.NDArray[np.floating[Any]],
    intensities: npt.NDArray[np.floating[Any]],
    axis_min: float,
    axis_max: float,
) -> float:
    """Return the share of a spectrum's TIC that survives a target axis range.

    Args:
        mzs: Source m/z values, ascending. Callers sort before calling;
            the result is meaningless on unsorted input.
        intensities: Source intensities, parallel to ``mzs``.
        axis_min: Lowest m/z on the target axis.
        axis_max: Highest m/z on the target axis.

    Returns:
        The total ion current the resampled spectrum should carry. Equal to
        ``intensities.sum()`` when the axis spans the spectrum, the sum of
        the samples inside ``[axis_min, axis_max]`` when the range cuts it,
        the area share when the range holds no sample, and ``0.0`` when the
        two ranges do not overlap.
    """
    if mzs.size == 0:
        return 0.0

    total_tic = float(np.sum(intensities))
    if total_tic <= 0.0:
        return 0.0

    # The axis covers the whole spectrum: nothing is dropped. Short-circuit
    # so the overwhelmingly common case is exact and cheap.
    if axis_min <= mzs[0] and axis_max >= mzs[-1]:
        return total_tic

    # The range cuts the spectrum: the current measured inside it, bounds
    # included, as nearest_neighbor keeps it (D71).
    inside = (mzs >= axis_min) & (mzs <= axis_max)
    if inside.any():
        return float(np.sum(intensities[inside]))

    # The range holds no sample. Only an area can say how much of the
    # spectrum lies between two samples.
    lo = max(axis_min, float(mzs[0]))
    hi = min(axis_max, float(mzs[-1]))
    if hi <= lo:
        # Disjoint, or touching at a single point of zero width.
        return 0.0

    total_integral = float(np.trapezoid(intensities, mzs))
    if total_integral <= 0.0:
        # Points that carry intensity but enclose no area say nothing
        # about the stretch between two of them.
        return 0.0

    # Interpolated values at the two cut points, so a window narrower than
    # the source spacing still gets a sensible share rather than zero. No
    # sample lies strictly between them.
    kept_mz = np.array([lo, hi])
    kept_integral = float(np.trapezoid(np.interp(kept_mz, mzs, intensities), kept_mz))

    fraction = min(max(kept_integral / total_integral, 0.0), 1.0)
    return total_tic * fraction


def rescale_to_preserved_tic(
    resampled: npt.NDArray[_FloatT],
    axis: npt.NDArray[np.floating[Any]],
    mzs: npt.NDArray[np.floating[Any]],
    intensities: npt.NDArray[np.floating[Any]],
    axis_range: Optional[Tuple[float, float]] = None,
) -> npt.NDArray[_FloatT]:
    """Scale an interpolated spectrum onto the TIC it is meant to carry.

    Modifies ``resampled`` in place and returns it.

    Args:
        resampled: Interpolated intensities on ``axis``, as returned by
            ``np.interp`` with out-of-range values set to zero.
        axis: The target mass axis, ascending.
        mzs: Source m/z values, ascending.
        intensities: Source intensities, parallel to ``mzs``.
        axis_range: The m/z range the axis was built to cover, when that is
            wider than the axis points themselves. A physics axis holds bin
            *centres*, so it stops half a bin short of the range it was
            asked for at either end, and measuring the preserved share
            against the end points alone would forfeit a peak sitting
            exactly on a declared bound (issue #239). ``None`` -- what a
            caller with no declared range passes, handed a bare axis --
            uses the axis's own span, the behaviour this function has
            always had.

    Returns:
        ``resampled``, scaled so that its sum is the preserved TIC.
    """
    if axis_range is None:
        axis_min, axis_max = float(axis[0]), float(axis[-1])
    else:
        axis_min, axis_max = float(axis_range[0]), float(axis_range[1])
    target = preserved_tic(mzs, intensities, axis_min, axis_max)
    if target <= 0.0:
        resampled.fill(0.0)
        return resampled

    current = float(np.sum(resampled))
    if current <= 0.0:
        # Nothing interpolated survived: every source peak landed strictly
        # between two axis points, or the gap tolerance vetoed every bin.
        # There is nothing to scale; tic_preserving_sparse places the
        # current in the nearest bins instead.
        return resampled

    resampled *= target / current
    return resampled
