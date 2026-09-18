# thyra/resampling/strategies/__init__.py

"""The per-spectrum resampling operators, and the factory that picks one.

:func:`build_strategy` is the whole construction contract. Its middle three
arguments -- the axis, the range it was laid across, and the linearisation
that axis may be indexed through -- are exactly the triple the axis
planning step settles once per conversion (issue #352); the method and the
gap tolerance come from the ``ResamplingConfig``. Nothing else reaches a
strategy, which is what keeps the operator independent of how the axis was
decided.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..types import AxisLinearisation, ResamplingMethod
from .base import ResamplingStrategy
from .nearest_neighbor import NearestNeighborStrategy
from .tic_preserving import TICPreservingStrategy

__all__ = [
    "ResamplingStrategy",
    "NearestNeighborStrategy",
    "TICPreservingStrategy",
    "build_strategy",
]


def build_strategy(
    method: ResamplingMethod,
    axis: NDArray[np.float64],
    axis_range: Optional[Tuple[float, float]],
    linearisation: Optional[AxisLinearisation],
    gap_tolerance_da: Optional[float],
) -> ResamplingStrategy:
    """The strategy a resampled conversion runs, one per conversion.

    Args:
        method: The resolved ``ResamplingMethod``. ``NONE`` is not one of
            these: a conversion that does not resample maps its peaks
            onto the reader's own axis and never builds a strategy.
        axis: The common mass axis, ascending.
        axis_range: The declared ``(min_mz, max_mz)`` the axis was laid
            across, or ``None``.
        linearisation: The axis's own coordinate if it may be used, as
            :func:`thyra.resampling.binning.usable_linearisation` returns
            it, else ``None``. Ignored by the methods that do not bin.
        gap_tolerance_da: From the ``ResamplingConfig``; ignored by the
            methods that do not interpolate.

    Returns:
        The strategy, ready to be called once per spectrum.

    Raises:
        ValueError: If ``method`` has no operator.
    """
    if method is ResamplingMethod.NEAREST_NEIGHBOR:
        return NearestNeighborStrategy(axis, axis_range, linearisation)
    if method is ResamplingMethod.TIC_PRESERVING:
        return TICPreservingStrategy(axis, axis_range, gap_tolerance_da)
    raise ValueError(f"No resampling strategy for method {method}")
