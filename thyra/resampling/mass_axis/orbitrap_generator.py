"""Orbitrap mass axis generator with bin size ∝ m/z^1.5 spacing."""

import numpy as np

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
        reference_mz: float = 1000.0,
        reference_width: float = 0.005,
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

        # Create uniform grid in 1/sqrt(mz) space. Walk 1/sqrt(mz) downwards,
        # from 1/sqrt(min_mz) to 1/sqrt(max_mz), so the resulting m/z axis
        # comes out ascending: it is assigned straight to the converter's
        # common mass axis, and everything downstream (np.searchsorted
        # binning, the stored var["mz"] column) requires increasing m/z.
        inv_sqrt_values = np.linspace(inv_sqrt_min, inv_sqrt_max, target_bins + 1)
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
