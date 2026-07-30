"""Linear/uniform mass axis generator with constant spacing."""

import numpy as np

from ..types import AxisType, MassAxis
from .base_generator import BaseAxisGenerator


class LinearAxisGenerator(BaseAxisGenerator):
    """Mass axis generator for uniform/constant spacing.

    Creates equidistant m/z bins with constant spacing across the entire
    range. This is the simplest and most common approach for mass axis
    generation.
    """

    def generate_axis(
        self,
        min_mz: float,
        max_mz: float,
        target_bins: int,
        reference_mz: float = 1000.0,
        reference_width: float = 0.005,
    ) -> MassAxis:
        """Generate uniform mass axis with constant spacing.

        Parameters
        ----------
        min_mz : float
            Minimum m/z value
        max_mz : float
            Maximum m/z value
        target_bins : int
            Number of bins
        reference_mz : float
            Unused. Constant spacing has no reference-m/z dependence; the
            parameter exists so every generator shares one signature.
        reference_width : float
            Unused, for the same reason. The realized width is
            ``(max_mz - min_mz) / (target_bins - 1)``, and the caller sizes
            ``target_bins`` to make that the width it wants.

        Returns
        -------
        MassAxis
            Generated mass axis with uniform spacing
        """
        mz_values = np.linspace(min_mz, max_mz, target_bins)

        return MassAxis(
            mz_values=mz_values,
            min_mz=float(mz_values[0]),
            max_mz=float(mz_values[-1]),
            num_bins=len(mz_values),
            axis_type=AxisType.CONSTANT,
        )

    def get_axis_type(self) -> AxisType:
        """Return the axis type for uniform spacing."""
        return AxisType.CONSTANT
