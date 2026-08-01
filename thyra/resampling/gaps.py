"""Refusing to interpolate across regions that were never measured.

``np.interp`` has no notion of a gap. Handed a spectrum whose m/z values are
sparse or thresholded, it draws a straight line from the last point before the
gap to the first point after it, and every target bin in between comes back
with a fabricated intensity. Nothing was measured there.

On real data this is not a rounding-scale effect. Forcing ``tic_preserving``
onto ``bellini.imzML`` -- 2,222 points per spectrum, median source spacing
0.0072 Da, largest gap 1,269 Da -- puts 80.6% of the output total ion current
into m/z regions more than 0.5 Da from any measured point.

The fix is the parameter Cardinal (``tolerance``) and matter
(``approx1(..., tol=)``) both expose: how far a target m/z may sit from the
nearest source m/z before its interpolated value is discarded rather than
trusted.

Order matters
-------------

Masking happens **before** the TIC rescale, not after. The intensity was
measured; it belongs somewhere. Interpolation had spread it across a gap where
nothing was recorded, so returning it to the bins that do have measurements
behind them keeps ``ResamplingMethod.TIC_PRESERVING``'s promise intact. Masking
after the rescale would instead delete that intensity and quietly break the
method's defining property.
"""

from typing import Any, Optional

import numpy as np
import numpy.typing as npt


def zero_across_gaps(
    resampled: npt.NDArray[Any],
    axis: npt.NDArray[np.floating[Any]],
    mzs: npt.NDArray[np.floating[Any]],
    tolerance: Optional[float],
) -> npt.NDArray[Any]:
    """Zero the bins of ``resampled`` that no source point vouches for.

    Modifies ``resampled`` in place and returns it.

    A bin survives when some source m/z lies within ``tolerance`` of it. On a
    source grid of uniform step ``s`` the farthest any target point can be
    from a source point is ``s / 2``, so any tolerance of at least half the
    source spacing leaves such a spectrum untouched -- which is what makes it
    safe to switch on for continuous profile data. Half the *widest* source
    gap is the rule to size it by; exactly half is a floating-point knife
    edge, so give it a hair of margin.

    Args:
        resampled: Interpolated intensities on ``axis``, modified in place.
        axis: The target mass axis, ascending.
        mzs: Source m/z values, ascending.
        tolerance: Maximum distance in Da from a target bin to the nearest
            source m/z. ``None`` disables the check and returns ``resampled``
            untouched, as does a non-positive tolerance.

    Returns:
        ``resampled``, with unsupported bins set to zero.
    """
    if tolerance is None or tolerance <= 0.0 or mzs.size == 0:
        return resampled

    # Distance from each target bin to the nearest source m/z. searchsorted
    # gives the insertion point; the nearest source point is on one side or
    # the other of it.
    right = np.searchsorted(mzs, axis)
    lo = np.clip(right - 1, 0, mzs.size - 1)
    hi = np.clip(right, 0, mzs.size - 1)
    nearest = np.minimum(np.abs(axis - mzs[lo]), np.abs(axis - mzs[hi]))

    resampled[nearest > tolerance] = 0.0
    return resampled
