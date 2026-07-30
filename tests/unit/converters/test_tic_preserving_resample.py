# tests/unit/converters/test_tic_preserving_resample.py
"""Tests that ``tic_preserving`` resampling preserves total ion current.

``_tic_preserving_resample`` used to be a bare ``np.interp`` with no
rescaling step, which meant the one property it is named for did not hold.
Interpolation samples the spectrum at every target point, so the output TIC
scaled with the width of the target axis: onto the default 190,000-bin axis
a 4,000-point profile spectrum came back with 47x its input TIC, and a
150-peak centroid spectrum with over 1000x. Nothing caught it because the
method had no direct test coverage.

The invariant these tests pin is that output TIC equals input TIC
*independently of bin count*. A test at a single bin count would have
passed against the broken implementation at 4,000 bins, where the inflation
factor happens to be 1.0 for this data, so the parametrisation over bin
counts is the point rather than incidental thoroughness.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from thyra.converters.spatialdata.base_spatialdata_converter import (
    BaseSpatialDataConverter,
)

MZ_MIN, MZ_MAX = 250.0, 1200.0

# Spans the realistic range: far below the source density, matched to it,
# and the 190,000 the auto default resolves to on this mass range.
BIN_COUNTS = [1_000, 4_000, 10_000, 50_000, 190_000]


def _resample(axis, mzs, intensities):
    """Call the resampler the way the converter does.

    ``_tic_preserving_resample`` reads only ``_common_mass_axis`` off
    ``self``, so a stub carrying it exercises the method without standing
    up a converter (which would need a reader and a writable output path).
    """
    stub = SimpleNamespace(_common_mass_axis=axis)
    return BaseSpatialDataConverter._tic_preserving_resample(stub, mzs, intensities)


def _axis(bins, lo=MZ_MIN, hi=MZ_MAX):
    return np.linspace(lo, hi, bins)


def _profile_spectrum():
    """A 4,000-point profile spectrum: peaks spread over many points."""
    mzs = np.linspace(MZ_MIN, MZ_MAX, 4000)
    intensities = np.zeros(4000)
    for centre, amplitude in [(760.6, 1.0), (782.6, 0.85), (888.6, 0.6)]:
        intensities += amplitude * np.exp(-((mzs - centre) ** 2) / (2 * 0.5**2))
    return mzs, intensities + 4.0


def _centroid_spectrum():
    """A sparse centroid spectrum: 150 discrete peaks."""
    rng = np.random.default_rng(0)
    mzs = np.sort(rng.uniform(MZ_MIN, MZ_MAX, 150))
    return mzs, rng.uniform(10.0, 1000.0, 150)


class TestTICIsPreserved:
    """The property the method is named for, across axis widths."""

    @pytest.mark.parametrize("bins", BIN_COUNTS)
    def test_profile_spectrum(self, bins):
        mzs, intensities = _profile_spectrum()
        out = _resample(_axis(bins), mzs, intensities)
        assert out.sum() == pytest.approx(intensities.sum(), rel=1e-9)

    @pytest.mark.parametrize("bins", BIN_COUNTS)
    def test_centroid_spectrum(self, bins):
        mzs, intensities = _centroid_spectrum()
        out = _resample(_axis(bins), mzs, intensities)
        assert out.sum() == pytest.approx(intensities.sum(), rel=1e-9)

    def test_tic_does_not_drift_with_bin_count(self):
        """The regression: TIC must not track the target axis width.

        Against the unrescaled implementation the ratio between the widest
        and narrowest axis was ~190x rather than 1.
        """
        mzs, intensities = _profile_spectrum()
        tics = [_resample(_axis(b), mzs, intensities).sum() for b in BIN_COUNTS]
        assert max(tics) / min(tics) == pytest.approx(1.0, rel=1e-9)

    @pytest.mark.parametrize("bins", [1_000, 190_000])
    def test_unsorted_input_is_handled(self, bins):
        mzs, intensities = _centroid_spectrum()
        order = np.random.default_rng(1).permutation(len(mzs))
        out = _resample(_axis(bins), mzs[order], intensities[order])
        assert out.sum() == pytest.approx(intensities.sum(), rel=1e-9)


class TestShapeAndDtype:
    """Callers index the result against the full axis."""

    @pytest.mark.parametrize("bins", [1_000, 10_000])
    def test_output_spans_the_axis(self, bins):
        mzs, intensities = _profile_spectrum()
        out = _resample(_axis(bins), mzs, intensities)
        assert out.shape == (bins,)
        assert out.dtype == np.float64

    def test_output_is_non_negative(self):
        mzs, intensities = _profile_spectrum()
        out = _resample(_axis(10_000), mzs, intensities)
        assert (out >= 0).all()


class TestCroppedAxis:
    """Cropping must drop signal, not redistribute it.

    ``np.interp`` zeroes anything outside the axis. Rescaling against the
    whole input TIC would push that dropped intensity back into the bins
    that survive, inventing signal inside the window the caller asked for.
    """

    def test_out_of_range_intensity_is_not_redistributed(self):
        # Five points; the axis covers a window around the middle one.
        mzs = np.array([300.0, 400.0, 500.0, 600.0, 700.0])
        intensities = np.array([100.0, 100.0, 50.0, 100.0, 100.0])
        axis = _axis(5_000, lo=450.0, hi=550.0)

        out = _resample(axis, mzs, intensities)

        # The kept share is measured by area, not by counting the source
        # points that happen to fall inside. Trapezoid area over the whole
        # spectrum is 35,000; over [450, 550] -- taking the interpolated
        # value 75 at each cut point and the sample 50 at 500 -- it is
        # 6,250. So 6250/35000 of a 450 TIC survives.
        assert out.sum() == pytest.approx(450.0 * 6250.0 / 35000.0, rel=1e-9)
        assert out.sum() < intensities.sum()

    def test_full_range_axis_preserves_whole_tic(self):
        """The default path: axis spans the spectrum, nothing is dropped."""
        mzs, intensities = _centroid_spectrum()
        out = _resample(_axis(10_000), mzs, intensities)
        assert out.sum() == pytest.approx(intensities.sum(), rel=1e-9)

    def test_window_narrower_than_source_spacing_is_not_blanked(self):
        """A window containing no source point must still carry signal.

        Counting the source points inside the range reports zero here --
        the nearest samples sit at 500.0 and 501.0, outside a window of
        [500.2, 500.8] -- and would blank a spectrum that interpolation
        represents perfectly well. Measuring the kept share by area does
        not have that discontinuity.
        """
        mzs = np.arange(400.0, 601.0, 1.0)
        intensities = np.full_like(mzs, 10.0)
        axis = _axis(500, lo=500.2, hi=500.8)

        out = _resample(axis, mzs, intensities)

        assert out.sum() > 0.0
        # Flat spectrum: the window holds 0.6 Da of a 200 Da span.
        assert out.sum() == pytest.approx(intensities.sum() * 0.6 / 200.0, rel=1e-6)

    def test_kept_share_grows_with_window_width(self):
        """Widening the window must not shrink what is preserved."""
        mzs = np.arange(400.0, 601.0, 1.0)
        intensities = np.full_like(mzs, 10.0)

        totals = [
            _resample(_axis(500, lo=500.0 - w, hi=500.0 + w), mzs, intensities).sum()
            for w in (0.3, 1.0, 5.0, 25.0)
        ]
        assert totals == sorted(totals)


class TestDegenerateInput:
    """Empty, single-point and zero-intensity spectra must not blow up."""

    def test_empty_spectrum_returns_zeros(self):
        out = _resample(_axis(1_000), np.array([]), np.array([]))
        assert out.shape == (1_000,)
        assert not out.any()

    def test_zero_intensity_spectrum_returns_zeros(self):
        mzs = np.linspace(MZ_MIN, MZ_MAX, 100)
        out = _resample(_axis(1_000), mzs, np.zeros(100))
        assert not out.any()
        assert np.isfinite(out).all()

    @pytest.mark.parametrize("bins", [1_000, 190_000])
    def test_single_peak_survives_and_keeps_its_intensity(self, bins):
        """A lone peak must not be interpolated away to nothing."""
        mzs = np.array([760.6])
        intensities = np.array([1234.5])
        out = _resample(_axis(bins), mzs, intensities)

        assert out.sum() == pytest.approx(1234.5, rel=1e-9)
        assert np.count_nonzero(out) == 1
        # It lands in the bin nearest its m/z.
        assert np.argmax(out) == np.argmin(np.abs(_axis(bins) - 760.6))

    def test_single_peak_outside_axis_is_dropped(self):
        out = _resample(_axis(1_000), np.array([50.0]), np.array([999.0]))
        assert not out.any()


class TestAgreementWithStrategyImplementation:
    """The converter must agree with the repo's other TIC-preserving code.

    ``thyra/resampling/strategies/tic_preserving.py`` already rescaled. The
    converter's inline copy did not, and the divergence is what the fix
    closes.
    """

    @pytest.mark.parametrize("bins", [1_000, 10_000])
    def test_same_tic_as_strategy(self, bins):
        from thyra.resampling.strategies.base import Spectrum
        from thyra.resampling.strategies.tic_preserving import TICPreservingStrategy

        mzs, intensities = _profile_spectrum()
        axis = _axis(bins)

        converter_out = _resample(axis, mzs, intensities)
        strategy_out = TICPreservingStrategy().resample(
            Spectrum(mz=mzs, intensity=intensities, coordinates=(0, 0, 0)), axis
        )

        assert converter_out.sum() == pytest.approx(
            float(np.sum(strategy_out.intensity)), rel=1e-9
        )
