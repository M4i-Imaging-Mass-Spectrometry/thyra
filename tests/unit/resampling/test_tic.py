# tests/unit/resampling/test_tic.py
"""Tests for the shared TIC-preservation rule.

``thyra.resampling.tic`` holds the single definition of *which* total
``tic_preserving`` preserves, used by both the converters' per-pixel path
and ``TICPreservingStrategy``. The two used to implement it separately and
drifted: the converter lost its rescaling step entirely, inflating TIC by
47x on the default axis. Pinning the rule here, rather than only through
its two callers, is what keeps them honest.

The rule: preserve the share of the spectrum lying inside the target axis
range, measured by area. Measuring by area rather than by counting the
source points inside the range is what keeps a window narrower than the
source spacing from collapsing to zero.
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.resampling.tic import preserved_tic, rescale_to_preserved_tic


def _flat(lo=400.0, hi=600.0, step=1.0, value=10.0):
    mzs = np.arange(lo, hi + step / 2, step)
    return mzs, np.full_like(mzs, value)


class TestUncroppedIsExact:
    """When the axis spans the spectrum, the whole TIC survives."""

    def test_axis_wider_than_spectrum(self):
        mzs, intensities = _flat()
        assert preserved_tic(mzs, intensities, 100.0, 900.0) == intensities.sum()

    def test_axis_exactly_matching_spectrum(self):
        mzs, intensities = _flat()
        # Exact equality, not approx: the common case is short-circuited
        # rather than computed through the integral ratio.
        assert preserved_tic(mzs, intensities, 400.0, 600.0) == intensities.sum()

    def test_result_is_a_plain_float(self):
        mzs, intensities = _flat()
        assert isinstance(preserved_tic(mzs, intensities, 100.0, 900.0), float)


class TestCropping:
    """A narrower axis keeps proportionally less."""

    @pytest.mark.parametrize(
        "lo,hi,expected_fraction",
        [
            (500.0, 600.0, 0.5),
            (450.0, 550.0, 0.5),
            (400.0, 500.0, 0.5),
            (475.0, 525.0, 0.25),
            (400.0, 600.0, 1.0),
        ],
    )
    def test_flat_spectrum_keeps_width_fraction(self, lo, hi, expected_fraction):
        mzs, intensities = _flat()
        kept = preserved_tic(mzs, intensities, lo, hi)
        assert kept == pytest.approx(intensities.sum() * expected_fraction, rel=1e-6)

    def test_kept_share_is_monotonic_in_width(self):
        mzs, intensities = _flat()
        widths = [0.2, 1.0, 5.0, 25.0, 100.0]
        kept = [preserved_tic(mzs, intensities, 500.0 - w, 500.0 + w) for w in widths]
        assert kept == sorted(kept)

    def test_never_exceeds_the_input_tic(self):
        mzs, intensities = _flat()
        for lo, hi in [(0.0, 1e6), (400.0, 600.0), (399.0, 601.0)]:
            assert preserved_tic(mzs, intensities, lo, hi) <= intensities.sum()


class TestNarrowWindow:
    """The case that a sample-counting rule gets wrong.

    A window between two adjacent source points contains no source point at
    all, so counting samples reports zero and blanks the spectrum.
    """

    def test_window_between_two_samples_is_not_zero(self):
        mzs, intensities = _flat(step=1.0)
        # Source points at 500.0 and 501.0; this window contains neither.
        kept = preserved_tic(mzs, intensities, 500.2, 500.8)
        assert kept > 0.0
        assert kept == pytest.approx(intensities.sum() * 0.6 / 200.0, rel=1e-6)

    def test_no_source_point_inside_still_scales_with_width(self):
        mzs, intensities = _flat(step=1.0)
        narrow = preserved_tic(mzs, intensities, 500.2, 500.4)
        wider = preserved_tic(mzs, intensities, 500.2, 500.8)
        assert 0.0 < narrow < wider


class TestNoOverlap:
    """Disjoint ranges keep nothing."""

    @pytest.mark.parametrize("lo,hi", [(0.0, 100.0), (900.0, 1000.0), (600.0, 700.0)])
    def test_disjoint_or_touching_returns_zero(self, lo, hi):
        mzs, intensities = _flat()
        assert preserved_tic(mzs, intensities, lo, hi) == 0.0


class TestDegenerateInput:
    """Empty, zero and single-point spectra must not raise or divide by zero."""

    def test_empty_spectrum(self):
        assert preserved_tic(np.array([]), np.array([]), 100.0, 200.0) == 0.0

    def test_zero_intensity_spectrum(self):
        mzs, _ = _flat()
        assert preserved_tic(mzs, np.zeros_like(mzs), 100.0, 900.0) == 0.0

    def test_single_point_inside_range(self):
        mzs = np.array([500.0])
        intensities = np.array([42.0])
        assert preserved_tic(mzs, intensities, 400.0, 600.0) == 42.0

    def test_single_point_outside_range(self):
        mzs = np.array([50.0])
        intensities = np.array([42.0])
        assert preserved_tic(mzs, intensities, 400.0, 600.0) == 0.0

    def test_repeated_mz_values_do_not_divide_by_zero(self):
        """Zero-width spectra integrate to nothing; fall back to samples."""
        mzs = np.array([500.0, 500.0, 500.0])
        intensities = np.array([1.0, 2.0, 3.0])
        kept = preserved_tic(mzs, intensities, 400.0, 600.0)
        assert np.isfinite(kept)
        assert kept == pytest.approx(6.0)


class TestRescale:
    """The in-place scaling wrapper."""

    def test_scales_to_the_preserved_total(self):
        mzs, intensities = _flat()
        axis = np.linspace(400.0, 600.0, 5_000)
        resampled = np.interp(axis, mzs, intensities, left=0.0, right=0.0)

        out = rescale_to_preserved_tic(resampled, axis, mzs, intensities)

        assert out.sum() == pytest.approx(intensities.sum(), rel=1e-9)

    def test_modifies_in_place_and_returns_the_same_array(self):
        mzs, intensities = _flat()
        axis = np.linspace(400.0, 600.0, 1_000)
        resampled = np.interp(axis, mzs, intensities, left=0.0, right=0.0)

        out = rescale_to_preserved_tic(resampled, axis, mzs, intensities)

        assert out is resampled

    def test_zero_preserved_total_blanks_the_output(self):
        mzs, intensities = _flat()
        # An axis disjoint from the spectrum keeps nothing.
        axis = np.linspace(900.0, 1000.0, 100)
        resampled = np.full(100, 5.0)

        out = rescale_to_preserved_tic(resampled, axis, mzs, intensities)

        assert not out.any()

    def test_all_zero_interpolation_is_left_alone(self):
        mzs, intensities = _flat()
        axis = np.linspace(400.0, 600.0, 100)
        resampled = np.zeros(100)

        out = rescale_to_preserved_tic(resampled, axis, mzs, intensities)

        assert np.isfinite(out).all()
        assert not out.any()


class TestBothCallersAgree:
    """The converter and the strategy must return the same totals.

    This is the regression that matters: the two implementations diverging
    is the original bug, in a more general form.
    """

    @pytest.mark.parametrize(
        "lo,hi",
        [(400.0, 600.0), (450.0, 550.0), (500.2, 500.8), (300.0, 700.0)],
    )
    @pytest.mark.parametrize("bins", [101, 5_000])
    def test_converter_and_strategy_match(self, lo, hi, bins):
        from types import SimpleNamespace

        from thyra.converters.spatialdata.base_spatialdata_converter import (
            BaseSpatialDataConverter,
        )
        from thyra.resampling.strategies.base import Spectrum
        from thyra.resampling.strategies.tic_preserving import TICPreservingStrategy

        mzs, intensities = _flat()
        axis = np.linspace(lo, hi, bins)

        from_converter = BaseSpatialDataConverter._tic_preserving_resample(
            SimpleNamespace(_common_mass_axis=axis), mzs, intensities
        )
        from_strategy = TICPreservingStrategy().resample(
            Spectrum(mz=mzs, intensity=intensities, coordinates=(0, 0, 0)), axis
        )

        assert from_converter.sum() == pytest.approx(
            float(np.sum(from_strategy.intensity)), rel=1e-9
        )
        # Not just the totals -- the distributions too.
        np.testing.assert_allclose(
            from_converter, from_strategy.intensity, rtol=1e-9, atol=1e-12
        )
