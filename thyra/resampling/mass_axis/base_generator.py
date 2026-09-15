"""Abstract base class for mass axis generators."""

from abc import ABC, abstractmethod
from typing import Any

import numpy as np
import numpy.typing as npt

from ..types import AxisLinearisation, AxisType, MassAxis


class BaseAxisGenerator(ABC):
    """Abstract base class for mass axis generators.

    A generator distributes ``target_bins`` bins across a mass range
    according to one analyser's spacing law. It does not decide *how many*
    bins to use: that comes from
    ``BaseSpatialDataConverter._calculate_bins_from_width``, which
    integrates ``1 / width(m)`` over the range so the width realized at
    ``reference_mz`` is the width the caller asked for. Generators are
    reached through :meth:`~thyra.resampling.common_axis.CommonAxisBuilder.build_physics_axis`.
    """

    @abstractmethod
    def generate_axis(
        self,
        min_mz: float,
        max_mz: float,
        target_bins: int,
        reference_mz: float = 1000.0,
        reference_width: float = 0.005,
    ) -> MassAxis:
        """Distribute ``target_bins`` bins across the mass range.

        The returned axis must be in ascending m/z order: it is assigned
        directly to the converter's common mass axis, where downstream
        binning and the stored ``var["mz"]`` column both require increasing
        m/z.

        Parameters
        ----------
        min_mz : float
            Minimum m/z value.
        max_mz : float
            Maximum m/z value.
        target_bins : int
            Number of bins to distribute.
        reference_mz : float
            Reference m/z the bin count was anchored to. Implementations
            that need only the bin count to realize their spacing law may
            ignore this.
        reference_width : float
            Bin width requested at ``reference_mz``. May likewise be
            ignored when the bin count already determines the spacing.

        Returns
        -------
        MassAxis
            The generated axis, ascending in m/z.
        """

    @abstractmethod
    def get_axis_type(self) -> AxisType:
        """Return the axis type this generator produces."""

    @abstractmethod
    def forward(self, mz: Any) -> npt.NDArray[np.float64]:
        """The coordinate this generator's bins are uniform in, at ``mz``.

        Strictly monotone in m/z. The axis is a uniform grid of bin
        *edges* in this coordinate (a uniform grid of centres for
        constant spacing), which is what lets a bin index be computed
        rather than searched: see :class:`~thyra.resampling.types.AxisLinearisation`.
        Vectorised over an array of m/z values.
        """

    def _midpoint_linearisation(self, u_edges: Any) -> AxisLinearisation:
        """The linearisation of an axis whose centres are m/z midpoints of these edges.

        ``u_edges`` is the uniform grid of ``target_bins + 1`` edge
        coordinates the generator laid with ``np.linspace``. Centre ``i``
        sits between edges ``i`` and ``i + 1``, so its linearised
        coordinate is half a step past edge ``i``.
        """
        u_edges = np.asarray(u_edges, dtype=np.float64)
        du = (float(u_edges[-1]) - float(u_edges[0])) / (u_edges.size - 1)
        return AxisLinearisation(self.forward, float(u_edges[0]) + 0.5 * du, du)
