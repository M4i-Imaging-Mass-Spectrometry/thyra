"""tic_preserving onto bins wider than the source sampling (issue #442, D44).

Interpolating at the bin centres is sampling: once a bin holds a sample that
neither of the two around its centre is, a peak narrower than the bin that
falls between two centres is mostly lost and the rescale hands its ion
current to whatever the centres did land on. The total still balances, so
nothing downstream can see it. A spectrum where that happens is rebinned
wherever a bin holds two or more samples; any other is resampled as before.

The reference here is the one the issue used: every source sample summed
into the bin it falls in, which is exact binning when bins are many samples
wide. It is computed independently of the strategy under test.
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.resampling.binning import BinLayout, bin_edges, kept_mz_range
from thyra.resampling.interpolation import sample_cells, wide_bins
from thyra.resampling.strategies import TICPreservingStrategy
from thyra.resampling.tic import rescale_to_preserved_tic


def _dense(strategy, mzs, intensities):
    idx, vals = strategy.resample(mzs, intensities)
    assert np.all(np.diff(idx) > 0)
    assert np.all(vals != 0)
    return strategy.to_dense(idx, vals)


def _summed_into_bins(axis, mzs, intensities):
    """Exact binning: each sample, whole, into the bin whose centre is nearest."""
    edges = 0.5 * (axis[1:] + axis[:-1])
    return np.bincount(
        np.searchsorted(edges, mzs, side="right"),
        weights=intensities,
        minlength=axis.size,
    )


def _two_narrow_peaks():
    """The issue's spectrum: 0.01 Da sampling, FWHM 0.05 Da, flat baseline."""
    mz = np.arange(100.0, 400.0, 0.01)

    def peak(centre, height):
        return height * np.exp(-0.5 * ((mz - centre) / (0.05 / 2.355)) ** 2)

    return mz, peak(167.05, 3865.0) + peak(333.10, 1187.0) + 1.0


class TestCoarseBinsKeepNarrowPeaks:
    """The issue's reproduction, on 0.73 Da bins."""

    AXIS = np.linspace(100.0, 400.0, 412)

    def _near(self, spectrum, centre):
        return spectrum[np.abs(self.AXIS - centre) <= 1.0].sum()

    def test_both_peaks_keep_their_ion_current(self):
        mz, intensity = _two_narrow_peaks()
        out = _dense(TICPreservingStrategy(self.AXIS, (100.0, 400.0)), mz, intensity)
        reference = _summed_into_bins(self.AXIS, mz, intensity)

        # Point sampling kept 421 and 415 of these.
        for centre in (167.05, 333.10):
            assert self._near(out, centre) == pytest.approx(
                self._near(reference, centre), rel=1e-3
            )

    def test_the_total_is_still_the_input_tic(self):
        mz, intensity = _two_narrow_peaks()
        out = _dense(TICPreservingStrategy(self.AXIS, (100.0, 400.0)), mz, intensity)
        assert out.sum() == pytest.approx(intensity.sum(), rel=1e-12)

    def test_the_spectrum_matches_exact_binning(self):
        mz, intensity = _two_narrow_peaks()
        out = _dense(TICPreservingStrategy(self.AXIS, (100.0, 400.0)), mz, intensity)
        reference = _summed_into_bins(self.AXIS, mz, intensity)

        assert np.corrcoef(out, reference)[0, 1] > 0.9999
        assert np.abs(out - reference).sum() / reference.sum() < 0.01

    def test_a_zero_suppressed_peak_between_two_centres_survives(self):
        """Every sample of the cluster lies between two bin centres.

        The interpolant is zero at both centres, so point sampling returned
        an empty spectrum and lost the pixel's whole ion current.
        """
        axis = np.arange(100.0, 110.0, 0.5)
        mz = np.array([104.15, 104.16, 104.17, 104.18, 104.19, 104.20])
        intensity = np.array([0.0, 30.0, 90.0, 100.0, 40.0, 0.0])

        out = _dense(TICPreservingStrategy(axis, None), mz, intensity)

        assert out.sum() == pytest.approx(intensity.sum())
        assert out[np.argmin(np.abs(axis - 104.18))] == pytest.approx(intensity.sum())


class TestSpectraThatSkipNoSampleAreResampledAsBefore:
    """No centre skips a sample: the operator is the old one, bin for bin."""

    @staticmethod
    def _old_operator(axis, mzs, intensities, axis_range):
        values = np.interp(axis, mzs, intensities, left=0.0, right=0.0)
        kept = kept_mz_range(axis, axis_range)
        return rescale_to_preserved_tic(values, axis, mzs, intensities, kept)

    def test_an_axis_finer_than_the_source(self):
        # A rapifleX-like grid (24.9 mDa) onto the default 5 mDa axis.
        rng = np.random.default_rng(0)
        mz = np.arange(100.0, 1000.0, 0.0249237)
        intensity = rng.gamma(0.5, 100.0, mz.size)
        axis = np.arange(100.3, 999.7, 0.005)

        out = _dense(TICPreservingStrategy(axis, None), mz, intensity)

        np.testing.assert_allclose(
            out, self._old_operator(axis, mz, intensity, None), rtol=1e-12, atol=0
        )

    @staticmethod
    def _waters_like_trace():
        """A zero-suppressed sqrt(m/z) trace and the default 1.14x axis."""
        k_axis = 1.3e-3 / np.sqrt(1000.0)
        k_grid = 1.143e-3 / np.sqrt(1000.0)
        n = int((2.0 / k_axis) * (np.sqrt(1000.0) - np.sqrt(300.0)))
        edges = np.linspace(np.sqrt(300.0), np.sqrt(1000.0), n + 1) ** 2
        axis = (edges[:-1] + edges[1:]) / 2
        grid = np.arange(np.sqrt(300.0), np.sqrt(1000.0), k_grid / 2) ** 2
        rng = np.random.default_rng(2)
        intensity = np.zeros(grid.size)
        for c in rng.choice(np.arange(20, grid.size - 20), 200, replace=False):
            intensity[c - 4 : c + 5] = rng.uniform(10, 5000) * np.exp(
                -0.5 * (np.arange(-4, 5) / 1.7) ** 2
            )
        keep = (intensity > 0) | np.roll(intensity > 0, 1) | np.roll(intensity > 0, -1)
        return axis, grid[keep], intensity[keep]

    def test_the_default_waters_profile_axis(self):
        """1.14 bins per sample: some bins hold two, never on one side.

        A zero-suppressed trace on a sqrt(m/z) grid, onto the default
        linear_tof axis of the Waters profile route, is left exactly as it
        was -- the route this fix must not move.
        """
        axis, mz, intensity = self._waters_like_trace()
        layout = BinLayout(axis, kept_mz_range(axis, None))
        bins = layout.locate(mz[(mz >= layout.lower) & (mz <= layout.upper)])
        assert (bins[1:] == bins[:-1]).any(), "the case needs bins holding two"

        out = _dense(TICPreservingStrategy(axis, None), mz, intensity)

        np.testing.assert_allclose(
            out, self._old_operator(axis, mz, intensity, None), rtol=1e-12, atol=0
        )

    def test_a_repeated_explicit_zero_skips_nothing(self):
        """A grid point stored twice with no intensity loses no ion current."""
        axis, mz, intensity = self._waters_like_trace()
        zero = int(np.flatnonzero(intensity == 0)[50])
        mz = np.insert(mz, zero, mz[zero])
        intensity = np.insert(intensity, zero, 0.0)

        out = _dense(TICPreservingStrategy(axis, None), mz, intensity)

        np.testing.assert_allclose(
            out, self._old_operator(axis, mz, intensity, None), rtol=1e-12, atol=0
        )

    def test_a_repeated_mz_with_intensity_is_rebinned(self):
        """Interpolation reads one of two samples at one m/z: that one skips."""
        axis, mz, intensity = self._waters_like_trace()
        peak = int(np.argmax(intensity))
        mz = np.insert(mz, peak, mz[peak])
        intensity = np.insert(intensity, peak, intensity[peak])
        layout = BinLayout(axis, kept_mz_range(axis, None))

        assert wide_bins(layout, mz, intensity).size > 0

    def test_an_axis_on_the_source_grid_returns_the_source(self):
        rng = np.random.default_rng(1)
        mz = np.linspace(500.0, 3000.0, 20001)
        intensity = rng.gamma(0.3, 100.0, mz.size)

        out = _dense(TICPreservingStrategy(mz, (500.0, 3000.0)), mz, intensity)

        np.testing.assert_array_equal(out, intensity)


class TestMixedSpectra:
    def test_both_kinds_of_bin_are_on_one_scale(self):
        """A flat trace gives the same ion current per Da in both regimes.

        The axis is fine on the left (a quarter of a sample per bin) and
        coarse on the right (six samples per bin). Rebinned bins hold ion
        current and interpolated ones a height, so without the width over
        spacing factor the right half would come out at 24 times its share.
        """
        mz = np.arange(100.0, 200.0, 0.01)
        intensity = np.full(mz.size, 7.0)
        axis = np.concatenate(
            (np.arange(100.00125, 150.0, 0.0025), np.arange(150.03, 199.95, 0.06))
        )

        out = _dense(TICPreservingStrategy(axis, None), mz, intensity)

        edges = bin_edges(axis, kept_mz_range(axis, None))
        per_da = out / np.diff(edges)
        interior = slice(5, -5)
        np.testing.assert_allclose(
            per_da[interior], intensity.sum() / (edges[-1] - edges[0]), rtol=0.02
        )

    def test_gap_tolerance_never_zeroes_a_bin_that_holds_samples(self):
        """A wide bin's centre can sit farther than the tolerance from a sample.

        Its samples are inside it, so it is measured, not bridged.
        """
        axis = np.arange(100.0, 110.0, 1.0)
        mz = np.array([104.6, 104.62, 104.64, 104.66])
        intensity = np.array([10.0, 50.0, 60.0, 20.0])

        out = _dense(
            TICPreservingStrategy(axis, None, gap_tolerance_da=0.05), mz, intensity
        )

        assert out[5] == pytest.approx(intensity.sum())

    def test_unsorted_input_and_repeated_mz(self):
        axis = np.linspace(100.0, 110.0, 11)
        mz = np.array([105.2, 101.1, 105.2, 105.25, 101.1, 101.15, 108.0])
        intensity = np.array([3.0, 1.0, 4.0, 2.0, 5.0, 6.0, 0.0])

        strategy = TICPreservingStrategy(axis, None, gap_tolerance_da=0.1)
        out = _dense(strategy, mz, intensity)
        order = np.argsort(mz, kind="stable")
        again = _dense(strategy, mz[order], intensity[order])

        np.testing.assert_allclose(out, again, rtol=1e-15)
        # The tolerance stops the bridge between the clusters, so each
        # cluster's bin holds exactly its own samples, repeats included.
        assert out[1] == pytest.approx(12.0) and out[5] == pytest.approx(9.0)
        assert out.sum() == pytest.approx(intensity.sum())


class TestEdges:
    def test_samples_on_the_outer_edges_are_kept_once(self):
        """On the lower edge and the top edge: in the end bins, whole.

        The tolerance stops the bridge across the middle, so the two end
        bins hold exactly the samples inside them.
        """
        axis = np.linspace(100.0, 110.0, 11)
        mz = np.array([100.0, 100.01, 100.02, 109.98, 109.99, 110.0])
        intensity = np.array([1.0, 2.0, 3.0, 5.0, 6.0, 7.0])

        out = _dense(
            TICPreservingStrategy(axis, (100.0, 110.0), gap_tolerance_da=0.1),
            mz,
            intensity,
        )

        assert out[0] == pytest.approx(6.0) and out[-1] == pytest.approx(18.0)
        assert out.sum() == pytest.approx(intensity.sum())


class TestBuildingBlocks:
    def test_bin_edges_are_the_nearest_centre_boundaries(self):
        axis = np.array([1.0, 2.0, 4.0])
        np.testing.assert_array_equal(
            bin_edges(axis, (0.0, 10.0)), [0.5, 1.5, 3.0, 5.0]
        )
        np.testing.assert_array_equal(bin_edges(axis, (1.0, 4.0)), [1.0, 1.5, 3.0, 4.0])

    def test_two_samples_either_side_of_a_centre_skip_nothing(self):
        layout = BinLayout(np.array([0.5, 1.5, 2.5]), (0.0, 3.0))
        # Bin 1 holds 1.0 (on its lower edge: ties go right) and 1.7.
        assert wide_bins(layout, np.array([0.5, 1.0, 1.7, 2.5]), np.ones(4)).size == 0

    def test_a_sample_on_the_centre_skips_its_neighbour_below(self):
        """Interpolation at 2.0 reads the sample on it and nothing else."""
        axis = np.array([1.0, 2.0, 3.0])
        layout = BinLayout(axis, (0.5, 3.5))
        mz = np.array([1.6, 2.0, 3.0])
        intensity = np.array([10.0, 1.0, 5.0])

        assert wide_bins(layout, mz, intensity).tolist() == [1]
        out = _dense(TICPreservingStrategy(axis, (0.5, 3.5)), mz, intensity)
        # Sampled at the centres, bin 1 kept 1 of the 11 it holds.
        assert out[1] > out[2]
        assert out.sum() == pytest.approx(intensity.sum())

    def test_two_samples_on_one_side_rebin_every_bin_holding_two(self):
        layout = BinLayout(np.array([0.5, 1.5, 2.5]), (0.0, 3.0))
        # Bin 2 holds 2.5 and 3.0 (the top edge closes the last bin), both
        # at or past its centre: the spectrum is rebinned, bin 1 included.
        wide = wide_bins(layout, np.array([0.5, 1.0, 1.5, 2.5, 3.0]), np.ones(5))
        assert wide.tolist() == [1, 2]

    def test_the_layout_draws_the_same_edges_as_bin_edges(self):
        axis = np.sqrt(np.linspace(100.0, 200.0, 1000)) ** 3
        kept = (axis[0] - 0.01, axis[-1] + 0.01)
        layout = BinLayout(axis, kept)
        np.testing.assert_array_equal(
            layout.edges_of(np.arange(axis.size + 1)), bin_edges(axis, kept)
        )

    def test_cells_tile_a_smooth_grid(self):
        mz = np.sqrt(np.linspace(100.0, 200.0, 50)) ** 3
        lo, hi = sample_cells(mz)
        np.testing.assert_array_equal(hi[:-1], lo[1:])

    def test_a_cell_does_not_reach_across_a_gap(self):
        mz = np.array([100.0, 100.01, 100.02, 110.0, 110.01])
        lo, hi = sample_cells(mz)
        assert hi[2] == pytest.approx(100.025)
        assert lo[3] == pytest.approx(109.995)
