"""FT-ICR mass axis generator with bin size ∝ m/z^2 spacing."""

import numpy as np

from ..types import AxisType, MassAxis
from .base_generator import BaseAxisGenerator


class FTICRAxisGenerator(BaseAxisGenerator):
    """Mass axis generator for FT-ICR (Fourier Transform Ion Cyclotron Resonance) analyzers.

    FT-ICR has bin size ∝ m/z^2, meaning spacing increases quadratically
    with mass. This reflects FT-ICR's cyclotron frequency detection
    where f ∝ 1/m/z, so equal frequency bins translate to m/z^2 spacing.
    """

    def generate_axis(
        self,
        min_mz: float,
        max_mz: float,
        target_bins: int,
        reference_mz: float = 1000.0,
        reference_width: float = 0.005,
    ) -> MassAxis:
        """Generate FT-ICR mass axis with m/z^2 spacing.

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
            Generated mass axis with FT-ICR spacing
        """
        # For FT-ICR: bin_width = k * mz^2
        # where k is determined by reference_width at reference_mz
        # k = reference_width / (reference_mz**2)  # Used for scaling

        # Generate axis by solving: integral of 1/mz^2 from min_mz to mz = target_position
        # Integral: -1/mz + 1/min_mz = target_position * (-1/max_mz + 1/min_mz) / target_bins

        inv_min = 1.0 / min_mz
        inv_max = 1.0 / max_mz
        # inv_range = inv_min - inv_max  # Note: min > max in 1/mz space

        # Create uniform grid in 1/mz space. Walk 1/mz downwards, from
        # 1/min_mz to 1/max_mz, so the resulting m/z axis comes out
        # ascending: it is assigned straight to the converter's common mass
        # axis, and everything downstream (np.searchsorted binning, the
        # stored var["mz"] column) requires increasing m/z.
        inv_values = np.linspace(inv_min, inv_max, target_bins + 1)
        mz_values = 1.0 / inv_values

        # Use bin centers
        mz_centers = (mz_values[:-1] + mz_values[1:]) / 2

        return MassAxis(
            mz_values=mz_centers,
            min_mz=float(mz_centers[0]),
            max_mz=float(mz_centers[-1]),
            num_bins=len(mz_centers),
            axis_type=AxisType.FTICR,
        )

    def calculate_width_at_mz(
        self,
        mz: float,
        reference_mz: float = 500.0,
        reference_width: float = 0.1,
    ) -> float:
        """Calculate expected bin width at given m/z for FT-ICR.

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
        return reference_width * ((mz / reference_mz) ** 2)

    def get_axis_type(self) -> AxisType:
        """Return the axis type for FT-ICR."""
        return AxisType.FTICR
