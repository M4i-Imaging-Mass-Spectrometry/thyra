"""Data types and enums for the resampling module."""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional

import numpy as np
import numpy.typing as npt


class ResamplingMethod(Enum):
    """Available resampling methods.

    Attributes:
        NONE: No resampling -- keep the original mass axis.
        NEAREST_NEIGHBOR: Snap each peak to the nearest target bin.
        TIC_PRESERVING: Redistribute intensity so the total ion count
            is preserved after rebinning (recommended for quantitative
            work).
        LINEAR_INTERPOLATION: Linear interpolation between neighbouring
            bins.
    """

    NONE = "none"
    NEAREST_NEIGHBOR = "nearest_neighbor"
    TIC_PRESERVING = "tic_preserving"
    LINEAR_INTERPOLATION = "linear_interpolation"


class AxisType(Enum):
    """Mass axis spacing model, determined by the analyser physics.

    The axis type controls how target bins are distributed across the
    mass range.  When set to ``None`` in :class:`ResamplingConfig`, the
    type is auto-detected from instrument metadata.

    Attributes:
        CONSTANT: Equidistant spacing (constant Da per bin).
        LINEAR_TOF: Linear TOF -- spacing proportional to
            ``sqrt(m/z)``.
        REFLECTOR_TOF: Reflector TOF -- spacing proportional to ``m/z``.
        ORBITRAP: Orbitrap -- spacing proportional to ``m/z^(3/2)``.
        FTICR: FTICR -- spacing proportional to ``m/z^2``.
        UNKNOWN: Unknown analyser; falls back to constant spacing.
    """

    CONSTANT = "constant"
    LINEAR_TOF = "linear_tof"
    REFLECTOR_TOF = "reflector_tof"
    ORBITRAP = "orbitrap"
    FTICR = "fticr"
    UNKNOWN = "unknown"


@dataclass
class MassAxis:
    """Represents a mass axis with metadata."""

    mz_values: npt.NDArray[np.floating[Any]]
    min_mz: float
    max_mz: float
    num_bins: int
    axis_type: AxisType

    @property
    def spacing(self) -> npt.NDArray[np.floating[Any]]:
        """Calculate spacing between consecutive m/z values."""
        return np.diff(self.mz_values)

    def resolution_at(self, mz: float) -> float:
        """Calculate resolution at given m/z."""
        idx = int(np.searchsorted(self.mz_values, mz))
        if idx > 0 and idx < len(self.mz_values):
            delta_mz = float(self.mz_values[idx] - self.mz_values[idx - 1])
            return mz / delta_mz
        return 0.0


@dataclass
class ResamplingConfig:
    """Configuration for resampling operations.

    All fields default to ``None`` (auto-detect from instrument metadata).
    You can override individual fields while leaving the rest automatic.

    Attributes:
        method: Resampling algorithm.  ``None`` auto-selects based on
            the instrument type.
        axis_type: Mass axis spacing model.  ``None`` auto-detects from
            the instrument metadata.
        target_bins: Number of bins in the resampled axis.  ``None``
            lets the resampler choose a bin count that preserves the
            native resolution.
        mass_width_da: Bin width in Daltons at ``reference_mz``.
            Alternative to ``target_bins`` -- specify one or the other.
        reference_mz: Reference m/z for ``mass_width_da``.  Default
            500.0 Da.
        min_mz: Override the lower bound of the mass range.
        max_mz: Override the upper bound of the mass range.
    """

    method: Optional[ResamplingMethod] = None
    axis_type: Optional[AxisType] = None
    target_bins: Optional[int] = None
    mass_width_da: Optional[float] = None
    reference_mz: float = 500.0
    min_mz: Optional[float] = None
    max_mz: Optional[float] = None
