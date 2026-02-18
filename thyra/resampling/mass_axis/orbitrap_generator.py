"""Orbitrap mass axis generator with bin size ∝ m/z^1.5 spacing."""

from typing import Any

import numpy as np
import numpy.typing as npt

from ..types import AxisType, MassAxis
from .base_generator import BaseAxisGenerator


class OrbitrapAxisGenerator(BaseAxisGenerator):
    """Mass axis generator for Orbitrap analyzers.

    Orbitrap has bin size ∝ m/z^1.5, meaning spacing increases faster at
    high mass. This reflects the Orbitrap's frequency-based detection
    where f ∝ 1/√m/z, so equal frequency bins translate to m/z^1.5
    spacing.
    """

    def generate_axis(
        self,
        min_mz: float,
        max_mz: float,
        target_bins: int,
        reference_mz: float = 500.0,
        reference_width: float = 0.1,
    ) -> MassAxis:
        """Generate Orbitrap mass axis with m/z^1.5 spacing.

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
            Generated mass axis with Orbitrap spacing
        """
        # For Orbitrap: bin_width = k * mz^1.5
        # where k is determined by reference_width at reference_mz
        # k = reference_width / (reference_mz**1.5)  # Used for scaling

        # Generate axis by solving: integral of 1/mz^1.5 from min_mz to mz = target_position
        # Integral: -2/sqrt(mz) + 2/sqrt(min_mz) =
        #   target_position * (-2/sqrt(max_mz) + 2/sqrt(min_mz)) / target_bins

        inv_sqrt_min = 1.0 / np.sqrt(min_mz)
        inv_sqrt_max = 1.0 / np.sqrt(max_mz)
        # inv_sqrt_range = inv_sqrt_min - inv_sqrt_max  # Note: min > max in 1/sqrt space

        # Create uniform grid in 1/sqrt(mz) space
        inv_sqrt_values = np.linspace(inv_sqrt_max, inv_sqrt_min, target_bins + 1)
        mz_values = 1.0 / (inv_sqrt_values**2)

        # Use bin centers
        mz_centers = (mz_values[:-1] + mz_values[1:]) / 2

        return MassAxis(
            mz_values=mz_centers,
            min_mz=float(mz_centers[0]),
            max_mz=float(mz_centers[-1]),
            num_bins=len(mz_centers),
            axis_type=AxisType.ORBITRAP,
        )

    def calculate_width_at_mz(
        self,
        mz: float,
        reference_mz: float = 500.0,
        reference_width: float = 0.1,
    ) -> float:
        """Calculate expected bin width at given m/z for Orbitrap.

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
        return float(reference_width * ((mz / reference_mz) ** 1.5))

    def get_axis_type(self) -> AxisType:
        """Return the axis type for Orbitrap."""
        return AxisType.ORBITRAP

    def generate_axis_bins(
        self, min_mz: float, max_mz: float, num_bins: int
    ) -> npt.NDArray[np.floating[Any]]:
        """Generate Orbitrap axis with fixed number of bins."""
        axis = self.generate_axis(min_mz, max_mz, num_bins)
        return axis.mz_values

    def generate_axis_width(
        self,
        min_mz: float,
        max_mz: float,
        width_da: float,
        reference_mz: float = 500.0,
    ) -> npt.NDArray[np.floating[Any]]:
        """Generate Orbitrap axis based on mass width at reference m/z."""
        # For Orbitrap: bin_width = k * mz^1.5
        k = width_da / (reference_mz**1.5)

        # Generate axis in 1/sqrt(mz) space for proportional m/z^1.5 spacing
        inv_sqrt_min = 1.0 / np.sqrt(min_mz)
        inv_sqrt_max = 1.0 / np.sqrt(max_mz)

        # Calculate number of bins needed for desired resolution
        inv_sqrt_step = k / (2 * width_da)  # Based on integral derivative
        num_bins = int((inv_sqrt_min - inv_sqrt_max) / inv_sqrt_step) + 1

        inv_sqrt_values = np.linspace(inv_sqrt_max, inv_sqrt_min, num_bins)
        return 1.0 / (inv_sqrt_values**2)
