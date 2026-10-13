# thyra/resampling/strategies/tic_preserving.py

"""The operator ``ResamplingMethod.TIC_PRESERVING`` selects.

Linear interpolation onto the target axis, followed by rescaling so the
resampled spectrum carries the same total ion current as the input --
evaluated only at the axis points the interpolant can be non-zero at. A
spectrum whose bins are wide enough to skip a sample is rebinned instead
(design decision D48). See :mod:`thyra.resampling.interpolation` for the
evaluation and :mod:`thyra.resampling.tic` for what "preserved" means when
the axis does not span the whole spectrum.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..binning import BinLayout
from ..interpolation import tic_preserving_sparse
from .base import ResamplingStrategy


class TICPreservingStrategy(ResamplingStrategy):
    """Interpolate onto the target axis, preserving total ion current.

    Args:
        axis: The target mass axis, ascending.
        axis_range: The declared ``(min_mz, max_mz)``, or ``None``.
        gap_tolerance_da: Target bins farther than this from any source
            m/z are zeroed before the rescale, so interpolation cannot
            claim regions nothing was measured in. ``None`` disables the
            check. See :mod:`thyra.resampling.gaps`.
    """

    def __init__(
        self,
        axis: NDArray[np.float64],
        axis_range: Optional[Tuple[float, float]],
        gap_tolerance_da: Optional[float] = None,
    ) -> None:
        """Bind the operator to one axis, its range and the gap tolerance."""
        super().__init__(axis, axis_range)
        self.gap_tolerance_da = gap_tolerance_da
        # A property of the axis, so laid once rather than per spectrum.
        self.layout: Optional[BinLayout] = (
            BinLayout(self.axis, self.kept_range) if self.axis.size >= 2 else None
        )

    def resample(
        self, mzs: NDArray[np.float64], intensities: NDArray[np.float64]
    ) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
        """Resample onto the common axis, preserving total ion current.

        Linear interpolation onto the target axis, followed by rescaling so
        the resampled spectrum carries the same total ion current as the
        input -- the behaviour ``ResamplingMethod.TIC_PRESERVING`` and
        docs/resampling.md both describe.

        The rescaling is not cosmetic. Interpolation samples the spectrum at
        every target point, so the raw interpolated sum scales with the
        density of the target axis rather than staying fixed. Onto the
        default 190,000-bin axis, 4,000 source points came back with 47x the
        input TIC, and a 150-peak centroid spectrum with over 1000x.

        The TIC preserved is the share of the spectrum lying inside the
        target axis range -- see ``thyra.resampling.tic``, which holds the
        rule and the reasoning. When the axis spans the spectrum, which is
        the default, that share is the whole input TIC. The range is the
        **declared** one, the same
        :class:`~thyra.resampling.strategies.nearest_neighbor.NearestNeighborStrategy`
        keeps peaks inside, so the two methods cannot disagree about what
        the axis covers.

        Only the axis points the interpolant can be non-zero at are
        evaluated, and only the non-zero results are returned. On a
        zero-suppressed profile source -- a Waters MRT pixel stores
        ~15,000 samples in clusters around its peaks, on a 1.05M-bin
        axis -- this is what turns a 570 s conversion into one that is
        bounded by reading the file. ``to_dense`` reproduces the dense
        evaluation bin for bin.

        Interpolating at a bin centre is sampling. Once a bin is wider than
        the source spacing, a peak narrower than the bin that falls between
        two centres is mostly lost, and the rescale gives its ion current
        to whatever the centres landed on: on 0.73 Da bins over 0.01 Da
        sampling, two narrow peaks kept 2% of their area while the TIC
        still balanced (issue #442). So when some bin holds a sample that
        neither of the two around its centre is, every bin holding two or
        more samples is filled with its share of each sample instead, and
        the spectrum's other bins are put on the same ion-current footing.
        A spectrum in which no centre skips a sample is resampled exactly
        as before. See design decision D48.

        Args:
            mzs: Original m/z values from the spectrum.
            intensities: Corresponding intensity values.

        Returns:
            Tuple of (bin_indices, intensities) containing only non-zero bins.
        """
        return tic_preserving_sparse(
            self.axis,
            mzs,
            intensities,
            self.gap_tolerance_da,
            self.axis_range,
            self.layout,
        )
