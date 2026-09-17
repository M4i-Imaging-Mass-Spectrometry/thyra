"""Tests for the interpolation gap tolerance.

``np.interp`` draws a straight line across any gap in the source m/z values,
so every target bin between two distant source points comes back with a
fabricated intensity. ``gap_tolerance_da`` is the parameter Cardinal
(``tolerance``) and matter (``approx1(..., tol=)``) both expose to stop that.

The property that matters most here is the one that makes it safe to switch
on: on a source grid of uniform step ``s``, no target point is farther than
``s / 2`` from a source point, so any tolerance of at least ``s / 2`` is a
no-op. That is the flexImaging/Rapiflex case, which is the only route
auto-selection sends to ``tic_preserving``.
"""

import numpy as np
import pytest

from thyra.resampling.gaps import zero_across_gaps
from thyra.resampling.interpolation import tic_preserving_sparse


class TestZeroAcrossGaps:
    """The masking primitive on its own."""

    def test_none_is_a_no_op(self):
        values = np.ones(5)
        axis = np.linspace(0.0, 4.0, 5)
        mzs = np.array([0.0])
        assert np.array_equal(zero_across_gaps(values, axis, mzs, None), np.ones(5))

    @pytest.mark.parametrize("tolerance", [0.0, -1.0])
    def test_non_positive_tolerance_is_a_no_op(self, tolerance):
        """Zero would otherwise mask everything not exactly on a source point."""
        values = np.ones(5)
        axis = np.linspace(0.0, 4.0, 5)
        mzs = np.array([0.0])
        assert np.array_equal(
            zero_across_gaps(values, axis, mzs, tolerance), np.ones(5)
        )

    def test_empty_source_is_a_no_op(self):
        values = np.ones(3)
        axis = np.linspace(0.0, 2.0, 3)
        assert np.array_equal(
            zero_across_gaps(values, axis, np.array([]), 0.1), np.ones(3)
        )

    def test_bins_beyond_the_tolerance_are_zeroed(self):
        axis = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        mzs = np.array([0.0, 4.0])
        values = np.ones(5)
        zero_across_gaps(values, axis, mzs, 1.0)
        # 0.0 and 4.0 are on source points; 1.0 and 3.0 are exactly one Da
        # away, which is within tolerance; 2.0 is two Da from either.
        assert np.array_equal(values, np.array([1.0, 1.0, 0.0, 1.0, 1.0]))

    def test_boundary_is_inclusive(self):
        axis = np.array([0.0, 1.0])
        mzs = np.array([0.0])
        values = np.ones(2)
        zero_across_gaps(values, axis, mzs, 1.0)
        assert np.array_equal(values, np.ones(2))

    def test_it_masks_in_place_and_returns_the_same_array(self):
        values = np.ones(3)
        out = zero_across_gaps(values, np.array([0.0, 5.0, 10.0]), np.array([0.0]), 1.0)
        assert out is values

    @pytest.mark.parametrize("step", [0.001, 0.005, 0.017, 1.0])
    def test_half_the_source_step_leaves_a_uniform_grid_untouched(self, step):
        """The guarantee that makes this safe on flexImaging data.

        RapiflexReader lays every spectrum out with np.linspace, so this is
        the shape of the only source grid auto-selection interpolates.

        The tolerance carries a relative hair above half the step. Exactly
        half is a floating-point knife edge: the worst-case target bin sits
        at the midpoint, and neither the accumulated grid nor the midpoint is
        exactly representable, so the computed distance lands either side of
        the bound by an ulp or two. Half the *widest* source gap is the
        honest rule, and anything above it is safely a no-op.
        """
        mzs = np.arange(100.0, 200.0, step)
        # A target axis deliberately offset from the source grid, so every
        # bin sits as far from a source point as the geometry allows.
        axis = mzs[:-1] + step / 2.0
        values = np.ones(axis.size)
        tolerance = float(np.diff(mzs).max()) / 2.0 * (1.0 + 1e-9)
        zero_across_gaps(values, axis, mzs, tolerance)
        assert np.count_nonzero(values) == values.size

    def test_exactly_half_the_step_on_an_exact_grid(self):
        """With representable arithmetic the bound holds exactly."""
        mzs = np.linspace(0.0, 1024.0, 1025)  # step 1.0, all exact
        axis = mzs[:-1] + 0.5
        values = np.ones(axis.size)
        zero_across_gaps(values, axis, mzs, 0.5)
        assert np.count_nonzero(values) == values.size


class TestTICPreservingWithGapTolerance:
    """The strategy end to end, and what masking does to the total."""

    def _resample(self, axis, mzs, intensities, gap_tolerance_da=None):
        indices, values = tic_preserving_sparse(
            axis,
            np.asarray(mzs, dtype=np.float64),
            np.asarray(intensities, dtype=np.float64),
            gap_tolerance_da,
            None,
        )
        out = np.zeros(axis.size)
        out[indices] = values
        return out

    def test_default_is_off(self):
        """Unset, the strategy behaves exactly as before."""
        axis = np.linspace(0.0, 100.0, 1001)
        without = self._resample(axis, [0.0, 100.0], [1.0, 1.0])
        assert np.count_nonzero(without) == axis.size

    def test_tolerance_confines_intensity_to_the_measured_regions(self):
        axis = np.linspace(0.0, 100.0, 1001)
        with_tol = self._resample(axis, [0.0, 100.0], [1.0, 1.0], gap_tolerance_da=1.0)
        # Only bins within 1 Da of 0.0 or 100.0 survive: 0.0-1.0 and
        # 99.0-100.0 at 0.1 Da spacing.
        assert np.count_nonzero(with_tol) == 22
        assert np.count_nonzero(with_tol[220:780]) == 0

    def test_tic_is_still_preserved(self):
        """Masking happens before the rescale, so the total is not lost."""
        axis = np.linspace(0.0, 100.0, 1001)
        out = self._resample(axis, [0.0, 100.0], [3.0, 5.0], gap_tolerance_da=1.0)
        assert out.sum() == pytest.approx(8.0)

    def test_uniform_source_is_unaffected(self):
        """A dense uniform source: same answer with and without the tolerance."""
        mzs = np.arange(100.0, 200.0, 0.01)
        rng = np.random.default_rng(0)
        intensities = rng.random(mzs.size)
        axis = np.linspace(100.0, 199.99, 20_000)

        without = self._resample(axis, mzs, intensities)
        with_tol = self._resample(axis, mzs, intensities, gap_tolerance_da=0.005)
        np.testing.assert_allclose(with_tol, without)

    def test_a_tolerance_that_masks_everything_leaves_zeros(self):
        """No surviving bin means nothing to rescale into -- not a crash."""
        axis = np.array([50.0, 51.0])
        out = self._resample(axis, [0.0, 100.0], [1.0, 1.0], gap_tolerance_da=1.0)
        assert np.array_equal(out, np.zeros(2))
