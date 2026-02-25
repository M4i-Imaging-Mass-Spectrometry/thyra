"""TIC-preserving linear interpolation strategy for profile data.

This strategy uses linear interpolation while preserving the Total Ion
Current (TIC) of the original spectrum, making it ideal for profile data
from most MS instruments.
"""

from typing import Any

import numpy as np
import numpy.typing as npt
from scipy import interpolate

from .base import ResamplingStrategy, Spectrum


class TICPreservingStrategy(ResamplingStrategy):
    """TIC-preserving linear interpolation strategy for profile data."""

    def resample(
        self, spectrum: Spectrum, target_axis: npt.NDArray[np.floating[Any]]
    ) -> Spectrum:
        """Resample spectrum using TIC-preserving linear interpolation.

        Uses scipy's linear interpolation and then scales the result
        to preserve the original Total Ion Current (TIC).

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
            # Handle single point spectrum with nearest neighbor
            distances = np.abs(target_axis - spectrum.mz[0])
            closest_idx = np.argmin(distances)
            resampled_intensity = np.zeros_like(target_axis)
            resampled_intensity[closest_idx] = spectrum.intensity[0]

            return Spectrum(
                mz=target_axis.copy(),
                intensity=resampled_intensity,
                coordinates=spectrum.coordinates,
                metadata=spectrum.metadata,
            )

        # Calculate original TIC
        original_tic = np.sum(spectrum.intensity)

        if original_tic == 0:
            # Handle zero intensity spectrum
            return Spectrum(
                mz=target_axis.copy(),
                intensity=np.zeros_like(target_axis),
                coordinates=spectrum.coordinates,
                metadata=spectrum.metadata,
            )

        # Perform linear interpolation
        # Use bounds_error=False and fill_value=0 for extrapolation
        interp_func = interpolate.interp1d(
            spectrum.mz,
            spectrum.intensity,
            kind="linear",
            bounds_error=False,
            fill_value=0.0,
            assume_sorted=False,  # Don't assume sorted for safety
        )

        # Interpolate to target axis
        interpolated_intensity = interp_func(target_axis)

        # Ensure no negative values (can happen with extrapolation)
        interpolated_intensity = np.maximum(interpolated_intensity, 0.0)

        # Calculate new TIC and preserve original TIC
        new_tic = np.sum(interpolated_intensity)

        if new_tic > 0:
            # Scale to preserve TIC
            scaling_factor = original_tic / new_tic
            interpolated_intensity *= scaling_factor
        else:
            # If interpolation resulted in zero TIC, handle gracefully
            # This shouldn't happen with proper data, but be defensive
            interpolated_intensity = np.zeros_like(target_axis)

        return Spectrum(
            mz=target_axis.copy(),
            intensity=interpolated_intensity,
            coordinates=spectrum.coordinates,
            metadata=spectrum.metadata,
        )
