"""TIC-preserving linear interpolation strategy for profile data.

This strategy uses linear interpolation while preserving the Total Ion
Current (TIC) of the original spectrum, making it ideal for profile data
from most MS instruments.

The rule for *which* total is preserved lives in ``thyra.resampling.tic``
and is shared with the converters' per-pixel path, so the two
implementations of ``tic_preserving`` cannot drift apart again. In short:
the preserved total is the share of the spectrum lying inside the target
axis range, so resampling onto an axis that crops the spectrum drops the
cropped intensity rather than redistributing it over the bins that remain.
When the axis spans the spectrum, the whole TIC is preserved exactly.
"""

from typing import Any, Optional

import numpy as np
import numpy.typing as npt

from ..gaps import zero_across_gaps
from ..tic import preserved_tic, rescale_to_preserved_tic
from .base import ResamplingStrategy, Spectrum


class TICPreservingStrategy(ResamplingStrategy):
    """TIC-preserving linear interpolation strategy for profile data."""

    def __init__(self, gap_tolerance_da: Optional[float] = None) -> None:
        """Initialize the strategy.

        Args:
            gap_tolerance_da: How far, in Daltons, a target bin may sit from
                the nearest source m/z before its interpolated value is
                discarded. ``None`` (the default) keeps ``np.interp``'s own
                behaviour of drawing straight lines across unmeasured
                regions. See :mod:`thyra.resampling.gaps`.
        """
        self.gap_tolerance_da = gap_tolerance_da

    def resample(
        self, spectrum: Spectrum, target_axis: npt.NDArray[np.floating[Any]]
    ) -> Spectrum:
        """Resample spectrum using TIC-preserving linear interpolation.

        Linearly interpolates onto ``target_axis`` and then scales the
        result so it carries the Total Ion Current the axis is entitled to
        -- the whole input TIC when the axis spans the spectrum, and the
        share lying inside the axis otherwise. See ``thyra.resampling.tic``
        for why the share is measured by integration rather than by counting
        the source points that fall in range.

        Parameters
        ----------
        spectrum : Spectrum
            Input spectrum to resample
        target_axis : npt.NDArray[np.floating[Any]]
            Target mass axis values

        Returns
        -------
        Spectrum
            Resampled spectrum with target_axis as mz values
        """
        if len(spectrum.mz) == 0:
            # Handle empty spectrum
            return Spectrum(
                mz=target_axis.copy(),
                intensity=np.zeros_like(target_axis),
                coordinates=spectrum.coordinates,
                metadata=spectrum.metadata,
            )

        if len(spectrum.mz) == 1:
            # A lone point cannot be interpolated onto a grid that does not
            # contain it, so place it in its nearest bin -- unless it falls
            # outside the axis, in which case it is cropped away like any
            # other out-of-range peak and preserved_tic returns 0.
            resampled_intensity = np.zeros_like(target_axis)
            kept = preserved_tic(
                np.asarray(spectrum.mz, dtype=np.float64),
                np.asarray(spectrum.intensity, dtype=np.float64),
                float(target_axis[0]),
                float(target_axis[-1]),
            )
            if kept > 0.0:
                closest_idx = np.argmin(np.abs(target_axis - spectrum.mz[0]))
                resampled_intensity[closest_idx] = kept

            return Spectrum(
                mz=target_axis.copy(),
                intensity=resampled_intensity,
                coordinates=spectrum.coordinates,
                metadata=spectrum.metadata,
            )

        if np.sum(spectrum.intensity) == 0:
            # Handle zero intensity spectrum
            return Spectrum(
                mz=target_axis.copy(),
                intensity=np.zeros_like(target_axis),
                coordinates=spectrum.coordinates,
                metadata=spectrum.metadata,
            )

        # np.interp and preserved_tic both require ascending m/z.
        order = np.argsort(spectrum.mz)
        mzs_sorted = np.asarray(spectrum.mz, dtype=np.float64)[order]
        intensities_sorted = np.asarray(spectrum.intensity, dtype=np.float64)[order]

        # Interpolate to target axis, dropping anything outside the source
        # range. Equivalent to the previous interp1d(bounds_error=False,
        # fill_value=0.0) for linear interpolation, without the scipy call.
        interpolated_intensity = np.interp(
            target_axis,
            mzs_sorted,
            intensities_sorted,
            left=0.0,
            right=0.0,
        )

        # Ensure no negative values (can happen with extrapolation)
        interpolated_intensity = np.maximum(interpolated_intensity, 0.0)

        # Discard bins no source point vouches for, before the rescale so the
        # intensity returns to the measured bins rather than being deleted.
        zero_across_gaps(
            interpolated_intensity,
            np.asarray(target_axis, dtype=np.float64),
            mzs_sorted,
            self.gap_tolerance_da,
        )

        # Scale onto the TIC the target axis is entitled to carry. Shared
        # with the converters' hot path so the two cannot drift apart again.
        rescale_to_preserved_tic(
            interpolated_intensity,
            np.asarray(target_axis, dtype=np.float64),
            mzs_sorted,
            intensities_sorted,
        )

        return Spectrum(
            mz=target_axis.copy(),
            intensity=np.asarray(interpolated_intensity, dtype=np.float64),
            coordinates=spectrum.coordinates,
            metadata=spectrum.metadata,
        )
