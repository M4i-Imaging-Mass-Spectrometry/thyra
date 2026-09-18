# thyra/resampling/strategies/base.py

"""What a resampling strategy is: one spectrum in, the bins it fills out.

A strategy is the per-spectrum operator of a conversion. It is built once,
from the target axis and the range that axis was laid across, and then
called once per spectrum per pass. It holds everything that is a property
of the axis pair rather than of the spectrum -- the kept range, a
linearisation, a shared-axis cache -- which is why it is an object at all
rather than a function the converter calls with six arguments.

The return is sparse on purpose: unique ascending bin indices and their
non-zero values, never a dense array of the axis's length. The dense form
of the TIC-preserving method cost 570 s on a Waters MRT run whose sparse
form is bounded by reading the file, because a zero-suppressed profile
pixel stores about 15,000 samples against an axis of 1.05M bins and every
other bin interpolates to exactly zero. ``to_dense`` is here for the
callers that want the array anyway, so no strategy has to build one.
"""

import logging
from abc import ABC, abstractmethod
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..binning import kept_mz_range

logger = logging.getLogger(__name__)


class ResamplingStrategy(ABC):
    """One resampling method, bound to one target mass axis.

    Args:
        axis: The target mass axis, ascending.
        axis_range: The declared ``(min_mz, max_mz)`` the axis was laid
            across, or ``None`` when no resampled axis was built. A
            physics generator reports bin *centres* and so stops half a
            bin short of this at either end, which is why the two are
            not the same thing -- see
            :func:`thyra.resampling.binning.kept_mz_range`.

    Attributes:
        axis: The target mass axis.
        axis_range: The declared range, as passed in.
        kept_range: The m/z range a peak has to be inside to survive,
            computed once here because it is a property of the axis and
            the range, not of any spectrum.
        out_of_range_peaks: Peaks discarded so far for lying outside it.
    """

    def __init__(
        self,
        axis: NDArray[np.float64],
        axis_range: Optional[Tuple[float, float]],
    ) -> None:
        """Bind the strategy to one axis and the range it was laid across."""
        self.axis: NDArray[np.float64] = np.asarray(axis, dtype=np.float64)
        self.axis_range = axis_range
        self.kept_range = kept_mz_range(self.axis, axis_range)
        # Peaks discarded for falling outside the target mass range, and
        # whether the one-line summary has been emitted yet. See
        # _count_out_of_range().
        self.out_of_range_peaks: int = 0
        self._out_of_range_warned: bool = False

    @abstractmethod
    def resample(
        self, mzs: NDArray[np.float64], intensities: NDArray[np.float64]
    ) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
        """Place one spectrum onto the axis.

        Args:
            mzs: The spectrum's m/z values.
            intensities: Its intensities, parallel to ``mzs``.

        Returns:
            ``(bin_indices, values)``: the bins this spectrum fills,
            ascending and unique, and what they hold. Bins that come out
            exactly zero are not returned -- the store carries no
            explicit zeros -- and neither array is of the axis's length.
        """

    def to_dense(
        self, indices: NDArray[np.int_], values: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Scatter a :meth:`resample` result into an array of the axis's length."""
        dense = np.zeros(len(self.axis))
        dense[indices] = values
        return dense

    def _count_out_of_range(self, n_dropped: int, n_total: int) -> None:
        """Record peaks discarded for lying outside the target mass range.

        Narrowing the mass range is deliberate, so dropping the peaks
        outside it is the correct answer and not an error -- but it is not
        something a user should have to infer either, since the previous
        behaviour conserved the total exactly and so left no trace a TIC
        check could find.

        Warns once per conversion rather than once per spectrum: a
        narrowed range typically excludes peaks in every spectrum, and on
        xenium that is 918,855 identical lines. ``out_of_range_peaks``
        keeps the running total; note it counts resample *calls*, and the
        converter resamples every spectrum twice (once per pass), so it
        is a lower bound on nothing and an upper bound on nothing -- read
        the warning, not the counter, if you want a per-spectrum figure.

        Args:
            n_dropped: Peaks outside the axis in this spectrum.
            n_total: Peaks in this spectrum before filtering.
        """
        self.out_of_range_peaks += n_dropped

        if self._out_of_range_warned:
            return
        self._out_of_range_warned = True
        lo_mz, hi_mz = self.kept_range
        logger.warning(
            "Dropping peaks that fall outside the target mass range "
            "[%.4f, %.4f] m/z -- %d of %d in the first spectrum affected. "
            "They are discarded, not folded into the edge bins. Widen the "
            "resampling range to keep them.",
            lo_mz,
            hi_mz,
            n_dropped,
            n_total,
        )
