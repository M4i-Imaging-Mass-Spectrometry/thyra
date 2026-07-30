"""Linear TOF mass axis generator with bin size ∝ √m/z spacing."""

import numpy as np

from ..types import AxisType, MassAxis
from .base_generator import BaseAxisGenerator


class LinearTOFAxisGenerator(BaseAxisGenerator):
    """Mass axis generator for Linear Time-of-Flight analyzers.

    Linear TOF has bin size ∝ √m/z, meaning spacing increases with the square root of mass.
    This reflects the fundamental TOF equation: t ∝ √(m/z), so equal time bins
    translate to √m/z spacing in mass.
    """

    def generate_axis(
        self,
        min_mz: float,
        max_mz: float,
        target_bins: int,
        reference_mz: float = 1000.0,
        reference_width: float = 0.005,
    ) -> MassAxis:
        """Generate Linear TOF mass axis with √m/z spacing.

        Parameters
        ----------
        min_mz : float
            Minimum m/z value
        max_mz : float
            Maximum m/z value
        target_bins : int
            Target number of bins
        reference_mz : float
            Reference m/z for width specification (default: 500.0)
        reference_width : float
            Mass width at reference m/z (default: 0.1)

        Returns
        -------
        MassAxis
            Generated mass axis with Linear TOF spacing
        """
        # For Linear TOF: bin_width = k * sqrt(mz)
        # where k is determined by reference_width at reference_mz
        # k = reference_width / np.sqrt(reference_mz)  # Used for scaling

        # Generate axis by solving: integral of 1/sqrt(mz) from min_mz to mz = target_position
        # Integral: 2*sqrt(mz) - 2*sqrt(min_mz) =
        #   target_position * (2*sqrt(max_mz) - 2*sqrt(min_mz)) / target_bins

        sqrt_min = np.sqrt(min_mz)
        sqrt_max = np.sqrt(max_mz)
        # sqrt_range = sqrt_max - sqrt_min  # Used for scaling

        # Create uniform grid in sqrt(mz) space
        sqrt_values = np.linspace(sqrt_min, sqrt_max, target_bins + 1)
        mz_values = sqrt_values**2

        # Use bin centers
        mz_centers = (mz_values[:-1] + mz_values[1:]) / 2

        return MassAxis(
            mz_values=mz_centers,
            min_mz=float(mz_centers[0]),
            max_mz=float(mz_centers[-1]),
            num_bins=len(mz_centers),
            axis_type=AxisType.LINEAR_TOF,
        )

    def calculate_width_at_mz(
        self,
        mz: float,
        reference_mz: float = 500.0,
        reference_width: float = 0.1,
    ) -> float:
        """Calculate expected bin width at given m/z for Linear TOF.

        Parameters
        ----------
        mz : float
            Target m/z value
        reference_mz : float
            Reference m/z for width specification
        reference_width : float
            Width at reference m/z

        Returns
        -------
        float
            Expected bin width at target m/z
        """
        return float(reference_width * np.sqrt(mz / reference_mz))

    def get_axis_type(self) -> AxisType:
        """Return the axis type for Linear TOF."""
        return AxisType.LINEAR_TOF
