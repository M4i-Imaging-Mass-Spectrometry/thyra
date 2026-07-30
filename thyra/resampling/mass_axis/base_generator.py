"""Abstract base class for mass axis generators."""

from abc import ABC, abstractmethod

from ..types import AxisType, MassAxis


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
