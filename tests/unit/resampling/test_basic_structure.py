"""
Tests for basic resampling module structure - Phase 1.
"""

import numpy as np
import pytest

from thyra.resampling import (
    AxisType,
    CommonAxisBuilder,
    MassAxis,
    ResamplingConfig,
    ResamplingDecisionTree,
    ResamplingMethod,
)
from thyra.resampling.mass_axis.base_generator import BaseAxisGenerator
from thyra.resampling.strategies.base import ResamplingStrategy, Spectrum


class TestBasicImports:
    """Test that all basic components can be imported."""

    def test_main_classes_import(self):
        """Test main classes can be imported."""
        tree = ResamplingDecisionTree()
        builder = CommonAxisBuilder()
        assert tree is not None
        assert builder is not None

    def test_enums_import(self):
        """Test enums are accessible."""
        methods = [m.value for m in ResamplingMethod]
        axis_types = [a.value for a in AxisType]

        assert "nearest_neighbor" in methods
        assert "tic_preserving" in methods
        assert "constant" in axis_types
        assert "fticr" in axis_types

    def test_dataclasses_import(self):
        """Test dataclasses can be created."""
        config = ResamplingConfig(target_bins=300000)
        assert config.target_bins == 300000
        assert config.reference_mz == 1000.0  # default, matches the CLI default

    def test_resampling_config_with_axis_type(self):
        """Test ResamplingConfig can be created with axis_type field."""
        # Test with explicit axis_type
        config = ResamplingConfig(
            method=ResamplingMethod.NEAREST_NEIGHBOR,
            axis_type=AxisType.REFLECTOR_TOF,
            target_bins=1000,
        )
        assert config.method == ResamplingMethod.NEAREST_NEIGHBOR
        assert config.axis_type == AxisType.REFLECTOR_TOF
        assert config.target_bins == 1000

        # Test without axis_type (should be None by default)
        config_no_axis = ResamplingConfig(method=ResamplingMethod.NEAREST_NEIGHBOR)
        assert config_no_axis.method == ResamplingMethod.NEAREST_NEIGHBOR
        assert config_no_axis.axis_type is None


class TestBasicFunctionality:
    """Test basic functionality works."""

    def test_decision_tree_basic(self):
        """Test decision tree basic behavior."""
        tree = ResamplingDecisionTree()

        # No metadata should raise NotImplementedError
        with pytest.raises(NotImplementedError):
            tree.select_strategy(None)

        # Unknown instruments fall through to DefaultDetector (returns NEAREST_NEIGHBOR)
        method = tree.select_strategy({"instrument_name": "Orbitrap Fusion"})
        assert method == ResamplingMethod.NEAREST_NEIGHBOR

    def test_axis_builder_uniform(self):
        """Test uniform axis building."""
        builder = CommonAxisBuilder()
        axis = builder.build_uniform_axis(100, 1000, 1000)

        assert isinstance(axis, MassAxis)
        assert len(axis.mz_values) == 1000
        assert axis.min_mz == 100
        assert axis.max_mz == 1000
        assert axis.axis_type == AxisType.CONSTANT

        # Test values are evenly spaced
        assert np.allclose(axis.mz_values[0], 100)
        assert np.allclose(axis.mz_values[-1], 1000)

    def test_mass_axis_properties(self):
        """Test MassAxis properties work."""
        mz_values = np.array([100, 200, 300, 400, 500])
        axis = MassAxis(
            mz_values=mz_values,
            min_mz=100,
            max_mz=500,
            num_bins=5,
            axis_type=AxisType.CONSTANT,
        )

        # Test spacing property
        spacing = axis.spacing
        expected_spacing = np.array([100, 100, 100, 100])
        assert np.array_equal(spacing, expected_spacing)

        # Test resolution calculation
        resolution = axis.resolution_at(250)
        # Resolution = m/z / delta_m/z = 250 / 100 = 2.5
        assert np.isclose(resolution, 2.5)


class TestAbstractClasses:
    """Test abstract base classes."""

    def test_resampling_strategy_abstract(self):
        """Test ResamplingStrategy cannot be instantiated."""
        with pytest.raises(TypeError):
            ResamplingStrategy()

    def test_base_axis_generator_abstract(self):
        """Test BaseAxisGenerator cannot be instantiated."""
        with pytest.raises(TypeError):
            BaseAxisGenerator()


class TestSpectrumDataclass:
    """Test Spectrum dataclass functionality."""

    def test_spectrum_creation(self):
        """Test Spectrum can be created."""
        mz = np.array([100, 200, 300])
        intensity = np.array([1000, 2000, 1500])
        coords = (5, 10, 0)

        spectrum = Spectrum(mz=mz, intensity=intensity, coordinates=coords)

        assert np.array_equal(spectrum.mz, mz)
        assert np.array_equal(spectrum.intensity, intensity)
        assert spectrum.coordinates == coords
        assert spectrum.metadata is None

    def test_spectrum_centroid_detection(self):
        """Test centroid detection heuristic."""
        # Small spectrum should be detected as centroid
        small_spectrum = Spectrum(
            mz=np.array([100, 200]),
            intensity=np.array([1000, 2000]),
            coordinates=(0, 0, 0),
        )
        assert small_spectrum.is_centroid is True

        # Large spectrum with many zeros should be centroid
        mz_large = np.linspace(100, 1000, 200)
        intensity_large = np.zeros(200)
        intensity_large[::20] = 1000  # Only every 20th point has intensity

        centroid_spectrum = Spectrum(
            mz=mz_large, intensity=intensity_large, coordinates=(0, 0, 0)
        )
        assert centroid_spectrum.is_centroid is True

        # Large spectrum with no zeros should be profile
        intensity_profile = np.ones(200) * 100  # All points have intensity

        profile_spectrum = Spectrum(
            mz=mz_large, intensity=intensity_profile, coordinates=(0, 0, 0)
        )
        assert profile_spectrum.is_centroid is False
