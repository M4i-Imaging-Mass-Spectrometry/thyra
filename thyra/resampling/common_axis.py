"""Common axis builder for creating unified mass axes."""

from typing import Dict, Optional, Tuple

from ..errors import ConversionRefused
from .mass_axis import (
    BaseAxisGenerator,
    FTICRAxisGenerator,
    LinearAxisGenerator,
    LinearTOFAxisGenerator,
    OrbitrapAxisGenerator,
    ReflectorTOFAxisGenerator,
    TOFAxisGenerator,
)
from .types import AxisType, MassAxis


class CommonAxisBuilder:
    """Creates optimized common mass axis for datasets."""

    def build_uniform_axis(
        self, min_mz: float, max_mz: float, num_bins: int
    ) -> MassAxis:
        """Create uniform (equidistant) mass axis.

        Delegates to :class:`~thyra.resampling.mass_axis.LinearAxisGenerator`,
        the generator for this axis type, rather than laying the bins a
        second time. Both are ``np.linspace(min_mz, max_mz, num_bins)``
        and always were, so the axis is unchanged -- but only the
        generator reports the coordinate it laid them in, and without
        that every placement onto a constant axis fell back to binary
        search, including the summed table's (issue #348). It is the one
        axis type whose linearisation is the identity, so it was also the
        cheapest one to be missing.

        Args:
            min_mz: Minimum m/z value
            max_mz: Maximum m/z value
            num_bins: Number of bins

        Returns:
            Generated mass axis
        """
        return LinearAxisGenerator().generate_axis(min_mz, max_mz, num_bins)

    def build_physics_axis(
        self,
        min_mz: float,
        max_mz: float,
        num_bins: int,
        axis_type: AxisType,
        reference_mz: float = 500.0,
        reference_width: float = 0.1,
        tof_law: Optional[Tuple[float, float]] = None,
    ) -> MassAxis:
        """Create physics-based mass axis for specific analyzer types.

        Args:
            min_mz: Minimum m/z value
            max_mz: Maximum m/z value
            num_bins: Number of bins
            axis_type: Type of mass analyzer
            reference_mz: Reference m/z for width specification (default:
                500.0)
            reference_width: Mass width at reference m/z (default: 0.1)
            tof_law: ``(A, B)`` of the two-term width law; required for
                ``AxisType.TOF`` and ignored otherwise.

        Returns:
            Generated mass axis with analyzer-specific spacing

        Raises:
            ConversionRefused: If axis_type is not supported, or is ``TOF``
                without a law
        """
        if axis_type is AxisType.TOF:
            if tof_law is None:
                raise ConversionRefused("AxisType.TOF needs tof_law=(A, B)")
            generator: BaseAxisGenerator = TOFAxisGenerator(*tof_law)
        else:
            generator_map: Dict[AxisType, BaseAxisGenerator] = {
                AxisType.CONSTANT: LinearAxisGenerator(),
                AxisType.LINEAR_TOF: LinearTOFAxisGenerator(),
                AxisType.REFLECTOR_TOF: ReflectorTOFAxisGenerator(),
                AxisType.ORBITRAP: OrbitrapAxisGenerator(),
                AxisType.FTICR: FTICRAxisGenerator(),
            }

            if axis_type not in generator_map:
                raise ConversionRefused(f"Unsupported axis type: {axis_type}")

            generator = generator_map[axis_type]

        return generator.generate_axis(
            min_mz, max_mz, num_bins, reference_mz, reference_width
        )
