"""The sparse form of ``tic_preserving`` is the dense form, bin for bin.

``tic_preserving`` used to interpolate onto every bin of the axis
and hand a dense array downstream, where it was scanned for non-zeros. On a
zero-suppressed profile source -- a Waters MRT pixel stores ~15,000 samples
in clusters around its peaks and the default axis has 1.05M bins -- that
cost 570 s for a conversion nearest-neighbour binning finished in 16 s. The
sparse form evaluates only the axis points the interpolant can be non-zero
at. These tests pin that it is the same operator: same values, same TIC,
same handling of gaps, cropping, unsorted input and degenerate spectra.
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.resampling.interpolation import (
    _overflow_safe_scale,
    sum_repeated_mz,
    tic_preserving_sparse,
    tic_support_bins,
)
from thyra.resampling.strategies import TICPreservingStrategy


def _sqrt_axis(lo, hi, width_at_1000):
    """A linear_tof axis: uniform in sqrt(m/z), like the digitiser grid."""
    k = width_at_1000 / np.sqrt(1000.0)
    n = int((2.0 / k) * (np.sqrt(hi) - np.sqrt(lo)))
    edges = np.linspace(np.sqrt(lo), np.sqrt(hi), n + 1) ** 2
    return (edges[:-1] + edges[1:]) / 2


def _zero_suppressed_profile(seed=0, n_peaks=40, spacing_at_1000=1.14e-3):
    """Clusters of samples on a sqrt-spaced grid, an explicit zero at each edge.

    This is the shape MassLynx returns for a profile pixel: only samples
    around peaks are stored, the sample positions lie on one global grid,
    and every cluster starts and ends with a stored zero.
    """
    rng = np.random.default_rng(seed)
    grid = _sqrt_axis(300.0, 1000.0, spacing_at_1000)
    centres = np.sort(rng.uniform(320.0, 980.0, n_peaks))
    mzs, its = [], []
    for c in centres:
        sigma = 1.9e-3 * np.sqrt(c / 1000.0)  # ~4 samples per FWHM
        i = np.searchsorted(grid, c)
        lo, hi = max(i - 12, 0), min(i + 13, grid.size)
        m = grid[lo:hi]
        y = rng.uniform(50, 5000) * np.exp(-((m - c) ** 2) / (2 * sigma**2))
        y[y < 20.0] = 0.0
        if y[0] != 0.0 or y[-1] != 0.0:
            continue
        # Trim to one zero on either side, as MassLynx stores it.
        nz = np.flatnonzero(y)
        m, y = m[nz[0] - 1 : nz[-1] + 2], y[nz[0] - 1 : nz[-1] + 2]
        if mzs and m[0] <= mzs[-1][-1]:
            continue
        mzs.append(m)
        its.append(y)
    return np.concatenate(mzs), np.concatenate(its)


AXIS = _sqrt_axis(300.0, 1000.0, 1.3e-3)


def _dense(axis, mzs, its, tol=None):
    strategy = TICPreservingStrategy(axis, None, tol)
    return strategy.to_dense(*strategy.resample(mzs, its))


def _sparse(axis, mzs, its, tol=None):
    return TICPreservingStrategy(axis, None, tol).resample(mzs, its)


def _scatter(axis, idx, vals):
    out = np.zeros(axis.size)
    out[idx] = vals
    return out


class TestSupportBins:
    def test_open_at_the_zero_neighbours_closed_at_the_ends(self):
        axis = np.arange(0.0, 10.5, 0.5)
        mzs = np.array([1.0, 2.0, 3.0, 4.0, 6.0, 7.0, 8.0, 9.0])
        its = np.array([5.0, 3.0, 0.0, 0.0, 0.0, 2.0, 0.0, 4.0])
        idx = tic_support_bins(axis, mzs, its)
        # First run [1, 2] starts at the spectrum's own first sample, so
        # axis == 1.0 is included; it ends before the zero at 3.0.
        # Run [7] sits between zeros at 6.0 and 8.0, open at both.
        # Run [9] ends at the spectrum's last sample: axis == 9.0 included.
        expected = np.flatnonzero(
            ((axis >= 1.0) & (axis < 3.0))
            | ((axis > 6.0) & (axis < 8.0))
            | ((axis > 8.0) & (axis <= 9.0))
        )
        np.testing.assert_array_equal(idx, expected)

    def test_no_zeros_means_the_whole_span(self):
        axis = np.linspace(0.0, 10.0, 41)
        mzs = np.array([2.0, 3.0, 5.0])
        idx = tic_support_bins(axis, mzs, np.array([1.0, 1.0, 1.0]))
        np.testing.assert_array_equal(
            idx, np.flatnonzero((axis >= 2.0) & (axis <= 5.0))
        )

    def test_all_zero_is_empty(self):
        idx = tic_support_bins(np.linspace(0, 1, 5), np.array([0.2, 0.4]), np.zeros(2))
        assert idx.size == 0

    def test_indices_are_sorted_and_unique(self):
        mzs, its = _zero_suppressed_profile()
        idx = tic_support_bins(AXIS, mzs, its)
        assert np.all(np.diff(idx) > 0)


class TestAgreesWithTheDenseForm:
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_zero_suppressed_profile(self, seed):
        mzs, its = _zero_suppressed_profile(seed)
        dense = _dense(AXIS, mzs, its)
        idx, vals = _sparse(AXIS, mzs, its)
        np.testing.assert_allclose(_scatter(AXIS, idx, vals), dense, rtol=1e-12, atol=0)
        assert np.all(vals != 0)
        assert vals.sum() == pytest.approx(its.sum(), rel=1e-12)

    def test_only_a_small_fraction_of_the_axis_is_touched(self):
        """The point of the exercise: work proportional to the samples."""
        mzs, its = _zero_suppressed_profile()
        idx, _ = _sparse(AXIS, mzs, its)
        assert idx.size < 2 * np.count_nonzero(its)
        assert idx.size < 0.05 * AXIS.size

    def test_centroid_like_input_without_zeros(self):
        rng = np.random.default_rng(3)
        mzs = np.sort(rng.uniform(300.0, 1000.0, 150))
        its = rng.uniform(10.0, 1000.0, 150)
        dense = _dense(AXIS, mzs, its)
        idx, vals = _sparse(AXIS, mzs, its)
        np.testing.assert_allclose(_scatter(AXIS, idx, vals), dense, rtol=1e-12, atol=0)

    def test_unsorted_input(self):
        mzs, its = _zero_suppressed_profile(4)
        order = np.random.default_rng(4).permutation(mzs.size)
        dense = _dense(AXIS, mzs, its)
        idx, vals = _sparse(AXIS, mzs[order], its[order])
        np.testing.assert_allclose(_scatter(AXIS, idx, vals), dense, rtol=1e-12, atol=0)

    def test_gap_tolerance_is_applied_the_same_way(self):
        rng = np.random.default_rng(5)
        mzs = np.sort(rng.uniform(300.0, 1000.0, 60))
        its = rng.uniform(10.0, 1000.0, 60)
        dense = _dense(AXIS, mzs, its, tol=0.05)
        idx, vals = _sparse(AXIS, mzs, its, tol=0.05)
        np.testing.assert_allclose(_scatter(AXIS, idx, vals), dense, rtol=1e-12, atol=0)
        assert vals.sum() == pytest.approx(its.sum(), rel=1e-12)

    def test_cropped_axis_keeps_the_same_share(self):
        mzs, its = _zero_suppressed_profile(6)
        axis = _sqrt_axis(500.0, 800.0, 1.3e-3)
        dense = _dense(axis, mzs, its)
        idx, vals = _sparse(axis, mzs, its)
        np.testing.assert_allclose(_scatter(axis, idx, vals), dense, rtol=1e-12, atol=0)
        assert 0 < vals.sum() < its.sum()

    def test_single_point(self):
        idx, vals = _sparse(AXIS, np.array([650.0]), np.array([42.0]))
        assert idx.size == 1 and vals[0] == 42.0
        assert AXIS[idx[0]] == pytest.approx(650.0, abs=1e-3)
        idx, vals = _sparse(AXIS, np.array([1200.0]), np.array([42.0]))
        assert idx.size == 0

    def test_empty_and_all_zero(self):
        idx, vals = _sparse(AXIS, np.array([]), np.array([]))
        assert idx.size == 0 and vals.size == 0
        idx, vals = _sparse(AXIS, np.array([400.0, 401.0]), np.zeros(2))
        assert idx.size == 0
        assert _dense(AXIS, np.array([400.0, 401.0]), np.zeros(2)).sum() == 0.0

    def test_module_function_is_what_the_strategy_calls(self):
        mzs, its = _zero_suppressed_profile(7)
        idx, vals = tic_preserving_sparse(AXIS, mzs, its, None)
        idx2, vals2 = _sparse(AXIS, mzs, its)
        np.testing.assert_array_equal(idx, idx2)
        np.testing.assert_array_equal(vals, vals2)


class TestOnTheDefaultMRTAxis:
    """Bins at 1.14x the sample spacing: nothing doubled, nothing skipped."""

    def test_no_holes_inside_clusters(self):
        mzs, its = _zero_suppressed_profile(8)
        idx, vals = _sparse(AXIS, mzs, its)
        stored = np.zeros(AXIS.size, dtype=bool)
        stored[idx] = True
        nz = its != 0
        for k in np.flatnonzero(nz[1:-1]) + 1:
            inside = (AXIS > mzs[k - 1]) & (AXIS < mzs[k + 1])
            assert stored[inside].all()

    def test_no_bin_exceeds_its_neighbouring_samples(self):
        """Interpolation cannot put two samples' worth into one bin."""
        mzs, its = _zero_suppressed_profile(9)
        idx, vals = _sparse(AXIS, mzs, its)
        scale = vals.sum() / np.interp(AXIS[idx], mzs, its).sum()
        right = np.searchsorted(mzs, AXIS[idx])
        left = np.clip(right - 1, 0, mzs.size - 1)
        right = np.clip(right, 0, mzs.size - 1)
        ceiling = np.maximum(its[left], its[right]) * scale
        assert np.all(vals <= ceiling * (1 + 1e-12))


_TIED_MZ = np.repeat([300.0, 300.5, 301.0, 301.5, 302.0], 3)
_TIED_A = np.array(
    [5, 80, 400, 300, 20, 1, 9, 9, 200, 60, 600, 3, 40, 7, 90], dtype=float
)


def _reverse_each_group(values):
    return values.reshape(5, 3)[:, ::-1].ravel()


class TestSumRepeatedMz:
    def test_ties_are_added_and_the_result_ascends(self):
        mzs, its = sum_repeated_mz(
            np.array([2.0, 1.0, 2.0, 3.0]), np.array([5.0, 7.0, 6.0, 1.0])
        )
        np.testing.assert_array_equal(mzs, [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(its, [7.0, 11.0, 1.0])

    def test_an_ascending_spectrum_is_returned_as_it_came(self):
        mzs = np.array([1.0, 2.0, 3.0])
        its = np.array([4.0, 5.0, 6.0])
        out_mzs, out_its = sum_repeated_mz(mzs, its)
        assert out_mzs is mzs
        assert out_its is its

    def test_the_sum_does_not_depend_on_the_order_of_the_ties(self):
        values = np.array([1e16, 1.0, -1e16, 1.0, 0.5])
        mzs = np.full(values.size, 10.0)
        forward = sum_repeated_mz(mzs, values)[1]
        backward = sum_repeated_mz(mzs, values[::-1].copy())[1]
        assert forward.tolist() == backward.tolist()


class TestTiedMzDoNotDependOnOrder:
    """A source listing one point per scan repeats each m/z (issue #516)."""

    # A kept range never narrows below the axis itself, so the window case
    # lays its axis inside the window.
    CASES = {
        "whole": (np.linspace(299.0, 303.0, 4001), (299.0, 303.0)),
        "window": (np.linspace(300.6, 301.4, 801), (300.6, 301.4)),
    }

    @pytest.mark.parametrize("case", ["whole", "window"])
    def test_reversed_ties_and_summed_ties_give_the_same_spectrum(self, case):
        axis, axis_range = self.CASES[case]
        forward = tic_preserving_sparse(axis, _TIED_MZ, _TIED_A, None, axis_range)
        backward = tic_preserving_sparse(
            axis, _TIED_MZ, _reverse_each_group(_TIED_A), None, axis_range
        )
        summed_mz, summed_a = sum_repeated_mz(_TIED_MZ, _TIED_A)
        summed = tic_preserving_sparse(axis, summed_mz, summed_a, None, axis_range)

        for other in (backward, summed):
            np.testing.assert_array_equal(forward[0], other[0])
            np.testing.assert_array_equal(forward[1], other[1])

    def test_the_whole_range_keeps_every_sample(self):
        axis, axis_range = self.CASES["whole"]
        _, values = tic_preserving_sparse(axis, _TIED_MZ, _TIED_A, None, axis_range)
        assert values.sum() == pytest.approx(_TIED_A.sum(), rel=1e-12)

    def test_a_window_keeps_the_samples_inside_it(self):
        axis, axis_range = self.CASES["window"]
        _, values = tic_preserving_sparse(axis, _TIED_MZ, _TIED_A, None, axis_range)
        # Only the tied group at 301.0 lies inside: 9 + 9 + 200.
        assert values.sum() == pytest.approx(218.0, rel=1e-12)


class TestNothingInterpolatedStillKeepsTheCurrent:
    """A pixel the interpolation leaves empty is placed, not dropped (#516)."""

    AXIS = np.linspace(100.0, 900.0, 800)

    def test_peaks_the_gap_tolerance_vetoes_land_in_their_nearest_bins(self):
        indices, values = tic_preserving_sparse(
            self.AXIS,
            np.array([200.0, 400.0, 600.0]),
            np.array([7.0, 7.0, 7.0]),
            0.05,
            (100.0, 900.0),
        )
        assert values.sum() == pytest.approx(21.0)
        assert indices.size == 3
        for index, mz in zip(indices, (200.0, 400.0, 600.0)):
            assert index == np.argmin(np.abs(self.AXIS - mz))

    def test_two_peaks_between_two_axis_points_are_kept(self):
        _, values = tic_preserving_sparse(
            self.AXIS,
            np.array([400.6, 400.7]),
            np.array([2.0, 3.0]),
            0.05,
            (100.0, 900.0),
        )
        assert values.sum() == pytest.approx(5.0)

    def test_two_samples_either_side_of_a_bin_edge_are_kept_without_a_tolerance(
        self,
    ):
        _, values = tic_preserving_sparse(
            self.AXIS,
            np.array([400.8, 400.95]),
            np.array([2.0, 3.0]),
            None,
            (100.0, 900.0),
        )
        assert values.sum() == pytest.approx(5.0)

    def test_a_spectrum_wholly_outside_the_range_stays_empty(self):
        indices, values = tic_preserving_sparse(
            self.AXIS,
            np.array([50.0, 950.0]),
            np.array([1.0, 1.0]),
            0.05,
            (100.0, 900.0),
        )
        assert indices.size == 0
        assert values.size == 0

    def test_a_spectrum_with_no_positive_current_stays_empty(self):
        indices, values = tic_preserving_sparse(
            self.AXIS,
            np.array([400.6, 400.7]),
            np.array([-1.0, -2.0]),
            None,
            (100.0, 900.0),
        )
        assert indices.size == 0
        assert values.size == 0


class TestHugeIntensitiesStayFinite:
    def test_a_value_near_the_float_maximum_does_not_overflow(self):
        axis = np.linspace(100.0, 300.0, 1001)
        _, values = tic_preserving_sparse(
            axis,
            np.array([100.0, 150.0, 200.0, 250.0, 300.0]),
            np.array([1e308, 1.0, 1.0, 1.0, 1.0]),
            None,
            (100.0, 300.0),
        )
        assert np.isfinite(values).all()
        assert values.sum() == pytest.approx(1e308 + 4.0, rel=1e-12)

    def test_an_ordinary_spectrum_is_not_scaled(self):
        assert _overflow_safe_scale(np.array([1e200, 1.0])) == 1.0
        assert _overflow_safe_scale(np.array([])) == 1.0

    def test_a_huge_spectrum_is_scaled_by_a_power_of_two(self):
        scale = _overflow_safe_scale(np.array([1e308, 1.0]))
        mantissa, _ = np.frexp(scale)
        assert mantissa == 0.5
        assert 1e308 * scale < 1.0


class TestSpectrumOutsideTheRangeIsCounted:
    def test_every_sample_is_counted_and_nothing_is_stored(self):
        strategy = TICPreservingStrategy(np.linspace(100.0, 110.0, 11), (100.0, 110.0))

        indices, values = strategy.resample(
            np.array([50.0, 60.0, 70.0]), np.array([1.0, 2.0, 3.0])
        )

        assert indices.size == 0
        assert values.size == 0
        assert strategy.out_of_range_peaks == 3
