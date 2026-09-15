"""Data types and enums for the resampling module."""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, Optional

import numpy as np
import numpy.typing as npt

#: Reference m/z that ``mass_width_da`` is anchored to when the caller
#: does not pick one. Matches the ``--resample-reference-mz`` CLI default
#: so the dataclass default and the CLI default cannot drift apart.
DEFAULT_REFERENCE_MZ = 1000.0


class ResamplingMethod(Enum):
    """Available resampling methods.

    Attributes:
        NONE: No resampling -- keep the original mass axis.
        NEAREST_NEIGHBOR: Snap each peak to the nearest target bin.
        TIC_PRESERVING: Redistribute intensity so the total ion count
            is preserved after rebinning (recommended for quantitative
            work).
    """

    NONE = "none"
    NEAREST_NEIGHBOR = "nearest_neighbor"
    TIC_PRESERVING = "tic_preserving"


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
        TOF: General time-of-flight -- spacing proportional to a measured
            peak width ``sqrt(A m + B m^2)``, of which ``LINEAR_TOF``
            (``B = 0``) and ``REFLECTOR_TOF`` (``A = 0``) are the limits.
            Needs the two coefficients; see
            :mod:`thyra.resampling.mass_axis.tof_generator`.
        ORBITRAP: Orbitrap -- spacing proportional to ``m/z^(3/2)``.
        FTICR: FTICR -- spacing proportional to ``m/z^2``.
        UNKNOWN: No analyser identified. Nothing in Thyra ever produces this
            member -- an undetected analyser is left as ``None``, which means
            "auto-detect" and does fall back to constant spacing. Reaching it
            therefore takes a caller who names it, and what happens then
            depends on how it is named. As a string in a ``dict``,
            ``{"axis_type": "unknown"}`` is refused during normalisation. As
            the enum member on a ``ResamplingConfig``, it is NOT refused
            there: that class is a plain dataclass with no validation, and
            normalisation returns a ``ResamplingConfig`` unchanged. It is
            stopped later, by
            :meth:`~thyra.resampling.common_axis.CommonAxisBuilder.build_physics_axis`,
            because there is no spacing model to build an axis from.
    """

    CONSTANT = "constant"
    LINEAR_TOF = "linear_tof"
    REFLECTOR_TOF = "reflector_tof"
    TOF = "tof"
    ORBITRAP = "orbitrap"
    FTICR = "fticr"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class AxisLinearisation:
    """The coordinate a generated axis is uniform in: ``forward(axis[i]) ~= u0 + i * du``.

    Every generator lays its bins uniformly in some analytic coordinate
    ``u = forward(m/z)`` -- m/z itself for constant spacing; ``sqrt``,
    ``ln``, ``1/sqrt`` and ``1/m`` for the single-term laws; the two-term
    law's own cumulative for ``TOF``. The bin nearest a value is therefore
    (nearly) ``(forward(mz) - u0) / du`` rounded, and the converter's
    nearest-neighbour mapping computes it that way instead of binary
    searching the axis (design decision D21).

    "Nearly" because the physics generators lay a uniform grid of bin
    *edges* and report the arithmetic midpoints in m/z, so
    ``forward(centre_i)`` is not exactly ``u0 + i * du``. It is strictly
    within half a step of it for any strictly monotone ``forward``: the
    midpoint lies strictly between its two edges, so its image lies
    strictly between the two edge coordinates. That bound is what makes a
    +-1 repair around the rounded index exact, and :meth:`deviation`
    measures the actual margin so the converter can check it before
    relying on it.
    """

    forward: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]]
    u0: float
    du: float

    def positions(self, mzs: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Where each m/z falls on the axis, in bins, before rounding."""
        u = np.asarray(
            self.forward(np.asarray(mzs, dtype=np.float64)), dtype=np.float64
        )
        return (u - self.u0) / self.du

    def deviation(self, axis: npt.NDArray[np.float64], chunk: int = 1 << 22) -> float:
        """The worst ``|positions(axis[i]) - i|`` over the axis, in bins.

        Chunked, because the axis can run to hundreds of millions of
        entries and this is called once per conversion for one number.
        """
        axis = np.asarray(axis, dtype=np.float64)
        worst = 0.0
        for start in range(0, axis.size, chunk):
            stop = min(start + chunk, axis.size)
            off = self.positions(axis[start:stop])
            off -= np.arange(start, stop, dtype=np.float64)
            worst = max(worst, float(np.max(np.abs(off))))
        return worst


@dataclass
class MassAxis:
    """Represents a mass axis with metadata."""

    mz_values: npt.NDArray[np.floating[Any]]
    min_mz: float
    max_mz: float
    num_bins: int
    axis_type: AxisType
    #: How the generator laid the bins, when it laid them in a coordinate
    #: with a closed form; ``None`` for an axis that was not generated.
    linearisation: Optional[AxisLinearisation] = None

    @property
    def spacing(self) -> npt.NDArray[np.floating[Any]]:
        """Calculate spacing between consecutive m/z values."""
        return np.diff(self.mz_values)

    def resolution_at(self, mz: float) -> float:
        """Return ``mz`` divided by the local bin width.

        This is a property of the AXIS, not of the data on it. A finely
        spaced axis reports a large number here whether or not the spectra
        stored on it resolve anything at that scale, so it must not be
        reported to users as the acquisition's resolving power.
        """
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
            derives the count from ``mass_width_da`` at ``reference_mz``
            across the mass range.  Note that this is a fixed target
            width, not a measurement: the source spectra's own point
            spacing is never consulted, so the resulting axis can be
            finer or coarser than the data it will hold.
        mass_width_da: Bin width in Daltons at ``reference_mz``.
            Alternative to ``target_bins`` -- specify one or the other.
        reference_mz: Reference m/z for ``mass_width_da``.  Default
            1000.0 Da, matching the ``--resample-reference-mz`` CLI
            default.
        min_mz: Override the lower bound of the mass range.
        max_mz: Override the upper bound of the mass range.
        gap_tolerance_da: How far, in Daltons, a target bin may sit from the
            nearest source m/z before ``tic_preserving`` discards its
            interpolated value instead of trusting it.  ``None`` -- the
            default -- means no check, which is how ``np.interp`` behaves:
            straight lines are drawn across regions where nothing was
            measured.  Only affects ``tic_preserving``; ``nearest_neighbor``
            never invents a bin.  See :mod:`thyra.resampling.gaps`.
        tof_a: ``A`` of the two-term TOF width law ``FWHM(m) = sqrt(A m +
            B m^2)`` (mDa^2/Da), for ``AxisType.TOF``.  ``None`` takes the
            pair the detected instrument declares, if any.
        tof_b: ``B`` of the same law (dimensionless).
        bins_per_fwhm: Bins per peak width for ``AxisType.TOF``.  ``None``
            means 3, or -- when ``mass_width_da`` is set -- whatever puts
            a bin of that width at ``reference_mz``.
    """

    method: Optional[ResamplingMethod] = None
    axis_type: Optional[AxisType] = None
    target_bins: Optional[int] = None
    mass_width_da: Optional[float] = None
    reference_mz: float = DEFAULT_REFERENCE_MZ
    min_mz: Optional[float] = None
    max_mz: Optional[float] = None
    gap_tolerance_da: Optional[float] = None
    tof_a: Optional[float] = None
    tof_b: Optional[float] = None
    bins_per_fwhm: Optional[float] = None
