"""
Tests for resampling strategies - Phase 2.
"""

import numpy as np
import pytest

from thyra.resampling.strategies import NearestNeighborStrategy, TICPreservingStrategy
from thyra.resampling.strategies.base import Spectrum


class TestNearestNeighborStrategy:
    """Test nearest neighbor resampling strategy."""

    def setup_method(self):
        """Setup test fixtures."""
        self.strategy = NearestNeighborStrategy()

    def test_basic_resampling(self):
        """Test basic nearest neighbor resampling."""
        # Create test spectrum with centroid-like data
        original_mz = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
        original_intensity = np.array([1000.0, 2000.0, 1500.0, 3000.0, 1200.0])
        spectrum = Spectrum(original_mz, original_intensity, (0, 0, 0))

        # Target axis with some points close to originals
        target_axis = np.array([99.0, 150.0, 201.0, 305.0, 450.0])

        # Resample
        result = self.strategy.resample(spectrum, target_axis)

        # Check structure
        assert np.array_equal(result.mz, target_axis)
        assert len(result.intensity) == len(target_axis)
        assert result.coordinates == (0, 0, 0)

        # Check nearest neighbor assignments
        # 99.0 -> nearest is 100.0 (intensity 1000)
        assert result.intensity[0] == 1000.0
        # 150.0 -> nearest is 100.0 (distance 50) vs 200.0 (distance 50) -> 100.0 wins (left bias)
        assert result.intensity[1] == 1000.0
        # 201.0 -> nearest is 200.0 (intensity 2000)
        assert result.intensity[2] == 2000.0
        # 305.0 -> nearest is 300.0 (intensity 1500)
        assert result.intensity[3] == 1500.0
        # 450.0 -> nearest is 400.0 (distance 50) vs 500.0 (distance 50) -> 400.0 wins (left bias)
        assert result.intensity[4] == 3000.0

    def test_empty_spectrum(self):
        """Test handling of empty spectrum."""
        spectrum = Spectrum(np.array([]), np.array([]), (0, 0, 0))
        target_axis = np.array([100.0, 200.0, 300.0])

        result = self.strategy.resample(spectrum, target_axis)

        assert np.array_equal(result.mz, target_axis)
        assert np.array_equal(result.intensity, np.zeros(3))

    def test_single_point_spectrum(self):
        """Test handling of single point spectrum."""
        spectrum = Spectrum(np.array([250.0]), np.array([1500.0]), (0, 0, 0))
        target_axis = np.array([100.0, 200.0, 300.0, 400.0])

        result = self.strategy.resample(spectrum, target_axis)

        # All points should get the intensity from the single point
        expected_intensity = np.array([1500.0, 1500.0, 1500.0, 1500.0])
        assert np.array_equal(result.intensity, expected_intensity)

    def test_extrapolation(self):
        """Test extrapolation beyond original range."""
        spectrum = Spectrum(
            np.array([200.0, 300.0]), np.array([2000.0, 1500.0]), (0, 0, 0)
        )
        target_axis = np.array([100.0, 250.0, 400.0])

        result = self.strategy.resample(spectrum, target_axis)

        # 100.0 -> nearest is 200.0 (first point)
        assert result.intensity[0] == 2000.0
        # 250.0 -> nearest is 200.0 vs 300.0, 200.0 is closer
        assert result.intensity[1] == 2000.0
        # 400.0 -> nearest is 300.0 (last point)
        assert result.intensity[2] == 1500.0

    def test_metadata_preservation(self):
        """Test that metadata is preserved."""
        spectrum = Spectrum(
            np.array([100.0, 200.0]),
            np.array([1000.0, 2000.0]),
            (5, 10, 2),
            metadata={"instrument": "timsTOF", "acquisition_mode": "centroid"},
        )
        target_axis = np.array([150.0, 250.0])

        result = self.strategy.resample(spectrum, target_axis)

        assert result.coordinates == (5, 10, 2)
        assert result.metadata == {
            "instrument": "timsTOF",
            "acquisition_mode": "centroid",
        }


class TestTICPreservingStrategy:
    """Test TIC-preserving linear interpolation strategy."""

    def setup_method(self):
        """Setup test fixtures."""
        self.strategy = TICPreservingStrategy()

    def test_basic_interpolation(self):
        """Test basic TIC-preserving interpolation."""
        # Create test spectrum with profile-like data
        original_mz = np.array([100.0, 150.0, 200.0, 250.0, 300.0])
        original_intensity = np.array([1000.0, 1500.0, 2000.0, 1500.0, 1000.0])
        spectrum = Spectrum(original_mz, original_intensity, (0, 0, 0))

        original_tic = np.sum(original_intensity)  # 7000.0

        # Target axis with different spacing. It spans the spectrum: this
        # test is about resampling onto a different grid, and the whole TIC
        # is preserved only when nothing is cropped away. Cropping is
        # covered separately by test_cropped_axis_drops_the_cropped_share.
        target_axis = np.array([100.0, 160.0, 220.0, 300.0])

        result = self.strategy.resample(spectrum, target_axis)

        # Check structure
        assert np.array_equal(result.mz, target_axis)
        assert len(result.intensity) == len(target_axis)

        # Check TIC preservation (within tolerance)
        new_tic = np.sum(result.intensity)
        assert abs(new_tic - original_tic) < 1e-10

        # Check that interpolation makes sense
        # All intensities should be positive
        assert np.all(result.intensity >= 0)

    def test_tic_preservation_exact(self):
        """Test exact TIC preservation with simple case."""
        # Simple linear spectrum
        original_mz = np.array([100.0, 200.0, 300.0])
        original_intensity = np.array([1000.0, 2000.0, 1000.0])
        spectrum = Spectrum(original_mz, original_intensity, (0, 0, 0))

        original_tic = 4000.0

        # Resample to higher resolution
        target_axis = np.linspace(100, 300, 21)  # 21 points

        result = self.strategy.resample(spectrum, target_axis)

        new_tic = np.sum(result.intensity)
        assert abs(new_tic - original_tic) < 1e-10

    def test_empty_spectrum(self):
        """Test handling of empty spectrum."""
        spectrum = Spectrum(np.array([]), np.array([]), (0, 0, 0))
        target_axis = np.array([100.0, 200.0, 300.0])

        result = self.strategy.resample(spectrum, target_axis)

        assert np.array_equal(result.mz, target_axis)
        assert np.array_equal(result.intensity, np.zeros(3))

    def test_single_point_spectrum(self):
        """Test handling of single point spectrum."""
        spectrum = Spectrum(np.array([250.0]), np.array([1500.0]), (0, 0, 0))
        target_axis = np.array([200.0, 250.0, 300.0])

        result = self.strategy.resample(spectrum, target_axis)

        # TIC should be preserved
        assert abs(np.sum(result.intensity) - 1500.0) < 1e-10
        # Only the closest point should have intensity
        closest_idx = np.argmin(np.abs(target_axis - 250.0))
        assert result.intensity[closest_idx] == 1500.0
        # Other points should be zero
        for i, intensity in enumerate(result.intensity):
            if i != closest_idx:
                assert intensity == 0.0

    def test_zero_intensity_spectrum(self):
        """Test handling of spectrum with zero intensities."""
        spectrum = Spectrum(np.array([100.0, 200.0]), np.array([0.0, 0.0]), (0, 0, 0))
        target_axis = np.array([150.0, 250.0])

        result = self.strategy.resample(spectrum, target_axis)

        assert np.array_equal(result.intensity, np.zeros(2))

    def test_extrapolation_handling(self):
        """Test extrapolation with fill_value=0."""
        spectrum = Spectrum(
            np.array([200.0, 300.0]), np.array([2000.0, 1000.0]), (0, 0, 0)
        )
        target_axis = np.array([100.0, 250.0, 400.0])

        result = self.strategy.resample(spectrum, target_axis)

        # Points outside range should have zero intensity
        assert result.intensity[0] == 0.0  # 100.0 is below range
        assert result.intensity[2] == 0.0  # 400.0 is above range
        # Point inside range should be interpolated
        assert result.intensity[1] > 0.0  # 250.0 is interpolated

        # TIC should still be preserved
        original_tic = 3000.0
        assert abs(np.sum(result.intensity) - original_tic) < 1e-10

    def test_metadata_preservation(self):
        """Test that metadata is preserved."""
        spectrum = Spectrum(
            np.array([100.0, 200.0]),
            np.array([1000.0, 2000.0]),
            (5, 10, 2),
            metadata={"instrument": "Orbitrap", "acquisition_mode": "profile"},
        )
        target_axis = np.array([120.0, 180.0])

        result = self.strategy.resample(spectrum, target_axis)

        assert result.coordinates == (5, 10, 2)
        assert result.metadata == {
            "instrument": "Orbitrap",
            "acquisition_mode": "profile",
        }

    def test_cropped_axis_drops_the_cropped_share(self):
        """An axis narrower than the spectrum must not keep the whole TIC.

        Rescaling to the full input TIC would push the intensity that
        interpolation dropped back into the bins that survive, so a cropped
        window would claim ions from the parts that were cut away.
        """
        # Flat spectrum over [400, 600]; keep the middle 50 Da.
        original_mz = np.arange(400.0, 601.0, 1.0)
        original_intensity = np.full_like(original_mz, 10.0)
        spectrum = Spectrum(original_mz, original_intensity, (0, 0, 0))

        target_axis = np.linspace(475.0, 525.0, 51)
        result = self.strategy.resample(spectrum, target_axis)

        kept = np.sum(result.intensity)
        assert kept < np.sum(original_intensity)
        # 50 Da of a 200 Da span, on a flat spectrum.
        assert kept == pytest.approx(
            np.sum(original_intensity) * 50.0 / 200.0, rel=1e-6
        )

    def test_axis_spanning_spectrum_preserves_tic_exactly(self):
        """The uncropped case is exact, not merely close."""
        original_mz = np.linspace(100.0, 300.0, 51)
        original_intensity = np.linspace(500.0, 1500.0, 51)
        spectrum = Spectrum(original_mz, original_intensity, (0, 0, 0))

        result = self.strategy.resample(spectrum, np.linspace(100.0, 300.0, 977))

        assert np.sum(result.intensity) == pytest.approx(
            np.sum(original_intensity), rel=1e-12
        )

    def test_single_point_outside_axis_is_dropped(self):
        """A lone peak outside the axis is cropped like any other."""
        spectrum = Spectrum(np.array([50.0]), np.array([999.0]), (0, 0, 0))
        result = self.strategy.resample(spectrum, np.array([200.0, 250.0, 300.0]))
        assert not result.intensity.any()

    def test_no_negative_intensities(self):
        """Test that no negative intensities are produced."""
        # Create a case that might produce negative values with interpolation
        original_mz = np.array([100.0, 110.0, 120.0])
        original_intensity = np.array([1000.0, 100.0, 1000.0])  # Dip in the middle
        spectrum = Spectrum(original_mz, original_intensity, (0, 0, 0))

        # Target points that might cause issues
        target_axis = np.array([95.0, 105.0, 115.0, 125.0])

        result = self.strategy.resample(spectrum, target_axis)

        # All intensities should be non-negative
        assert np.all(result.intensity >= 0.0)


class TestStrategyComparison:
    """Test comparing strategies on the same data."""

    def test_centroid_vs_profile_behavior(self):
        """Test that strategies behave differently for centroid vs profile data."""
        # Centroid-like data (sparse with zeros)
        centroid_mz = np.array([100.0, 200.0, 300.0])
        centroid_intensity = np.array([1000.0, 2000.0, 1500.0])
        centroid_spectrum = Spectrum(centroid_mz, centroid_intensity, (0, 0, 0))

        # Profile-like data (dense)
        profile_mz = np.linspace(100, 300, 21)
        profile_intensity = np.concatenate(
            [
                np.linspace(500, 1000, 7),  # Rising
                np.linspace(1000, 2000, 7),  # Peak
                np.linspace(2000, 500, 7),  # Falling
            ]
        )
        profile_spectrum = Spectrum(profile_mz, profile_intensity, (0, 0, 0))

        # Common target axis. Spans the profile spectrum so its whole TIC
        # is preserved; a narrower axis would legitimately preserve less.
        target_axis = np.linspace(100, 300, 9)

        # Apply both strategies
        nn_strategy = NearestNeighborStrategy()
        tic_strategy = TICPreservingStrategy()

        nn_result = nn_strategy.resample(centroid_spectrum, target_axis)
        tic_result = tic_strategy.resample(profile_spectrum, target_axis)

        # Both should preserve some form of total intensity
        _centroid_tic = np.sum(centroid_spectrum.intensity)  # noqa: F841
        profile_tic = np.sum(profile_spectrum.intensity)

        nn_tic = np.sum(nn_result.intensity)
        tic_tic = np.sum(tic_result.intensity)

        # Nearest neighbor preserves peak intensities but may change TIC
        assert nn_tic > 0
        # TIC preserving should exactly preserve TIC
        assert abs(tic_tic - profile_tic) < 1e-10

        # Nearest neighbor should have step-like behavior
        # TIC preserving should have smooth interpolation
        assert len(np.unique(nn_result.intensity)) <= len(
            centroid_intensity
        )  # Limited unique values
        # TIC result should be smoother (more unique values for dense profile data)
