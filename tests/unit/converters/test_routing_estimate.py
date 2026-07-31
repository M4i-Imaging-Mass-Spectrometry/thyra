# tests/unit/converters/test_routing_estimate.py
"""Tests that CSR-vs-CSC routing is decided on the real bin count.

``_estimate_output_size_gb`` scores a dataset as
``n_pixels * n_mz_bins * 4`` and ``_should_use_pcs`` compares that against
``PCS_SIZE_THRESHOLD_GB``. The resampling branch already resolved the real
bin count (issue #87), but the raw-axis branch still guessed it from
``(max_mass - min_mass) / 0.01`` -- a 10 mDa spacing the data need not have.

That guess is wrong in both directions. A continuous file carrying 4,000
points over 250-1200 m/z was scored as though it had 95,000, inflating it
24x; a processed file whose spectra share no m/z values was scored far too
low, which routed the very largest datasets to the method that holds the
most in memory.

``convert()`` runs ``_initialize_conversion()`` before it asks which method
to use, so the axis is already built when the gate is evaluated and there is
nothing to estimate. These tests pin that the built axis is preferred, and
that the old heuristics still apply when it is not available -- which is the
case when the estimator is called directly, as the older tests in
``test_streaming_converter.py`` do.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pytest

from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)

THRESHOLD = StreamingSpatialDataConverter.PCS_SIZE_THRESHOLD_GB


class _Meta:
    def __init__(self, dimensions, mass_range):
        self.dimensions = dimensions
        self.mass_range = mass_range
        self.n_spectra = dimensions[0] * dimensions[1] * dimensions[2]
        self.pixel_size = (20.0, 20.0)
        self.estimated_memory_gb = 1.0
        self.coordinate_bounds = (0, dimensions[0], 0, dimensions[1])
        self.is_3d = dimensions[2] > 1
        self.has_mass_axis = True
        self.source_path = "mock"


class _Reader:
    def __init__(self, dimensions, mass_range):
        self._meta = _Meta(dimensions, mass_range)

    def get_essential_metadata(self):
        return self._meta


def _converter(
    dimensions: Tuple[int, int, int],
    mass_range: Tuple[float, float] = (250.0, 1200.0),
    axis: Optional[np.ndarray] = None,
    resampling_config=None,
):
    """Build a converter, optionally with the mass axis already resolved."""
    with tempfile.TemporaryDirectory() as tmpdir:
        conv = StreamingSpatialDataConverter(
            reader=_Reader(dimensions, mass_range),
            output_path=Path(tmpdir) / "out.zarr",
            use_csc="auto",
            resampling_config=resampling_config,
        )
    conv._common_mass_axis = axis
    return conv


class TestPrefersTheBuiltAxis:
    """When the axis exists, its length is what counts."""

    def test_uses_real_axis_length(self):
        n_pixels = 10_000
        axis = np.linspace(250.0, 1200.0, 4_000)
        conv = _converter((100, 100, 1), axis=axis)

        expected = n_pixels * 4_000 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)

    def test_real_axis_beats_the_raw_heuristic(self):
        """The heuristic would say 95,000 bins; the axis says 4,000."""
        axis = np.linspace(250.0, 1200.0, 4_000)
        conv = _converter((100, 100, 1), axis=axis)
        with_axis = conv._estimate_output_size_gb()

        conv._common_mass_axis = None
        without_axis = conv._estimate_output_size_gb()

        # (1200 - 250) / 0.01 = 95,000, i.e. 23.75x the real count.
        assert without_axis == pytest.approx(with_axis * 95_000 / 4_000, rel=1e-6)

    def test_real_axis_beats_the_resampling_plan_too(self):
        """A built axis is authoritative even when resampling is configured."""
        axis = np.linspace(250.0, 1200.0, 1_000)
        conv = _converter(
            (100, 100, 1),
            axis=axis,
            resampling_config={"method": "nearest_neighbor", "target_bins": 500_000},
        )
        expected = 10_000 * 1_000 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)


class TestFallbacksStillApply:
    """Without a built axis, the previous behaviour is unchanged."""

    def test_raw_heuristic_when_no_axis_and_no_resampling(self):
        conv = _converter((100, 100, 1), mass_range=(250.0, 1200.0), axis=None)
        expected = 10_000 * 95_000 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)

    def test_resampling_plan_when_no_axis(self):
        conv = _converter(
            (50, 50, 1),
            mass_range=(100.0, 1000.0),
            axis=None,
            resampling_config={"method": "nearest_neighbor", "target_bins": 5_000},
        )
        expected = 2_500 * 5_000 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)


class TestRoutingDecision:
    """The estimate exists to pick a method, so pin the decision itself."""

    def test_narrow_real_axis_routes_to_coo(self):
        """Was PCS on an estimate inflated 24x by the 10 mDa assumption."""
        axis = np.linspace(250.0, 1200.0, 4_000)
        conv = _converter((1000, 600, 1), axis=axis)

        assert conv._estimate_output_size_gb() < THRESHOLD
        assert conv._should_use_pcs() is False

        # Same dataset under the old rule would have been sent to PCS.
        conv._common_mass_axis = None
        assert conv._estimate_output_size_gb() > THRESHOLD

    def test_wide_real_axis_routes_to_pcs(self):
        """Was COO despite being genuinely large."""
        axis = np.linspace(250.0, 1200.0, 200_000)
        conv = _converter((250, 200, 1), axis=axis)

        assert conv._estimate_output_size_gb() > THRESHOLD
        assert conv._should_use_pcs() is True

        # The old rule scored this below the threshold.
        conv._common_mass_axis = None
        assert conv._estimate_output_size_gb() < THRESHOLD

    def test_explicit_use_csc_still_wins(self):
        """The estimate is only consulted in auto mode."""
        axis = np.linspace(250.0, 1200.0, 10)
        conv = _converter((10, 10, 1), axis=axis)

        conv._use_csc = True
        assert conv._should_use_pcs() is True
        conv._use_csc = False
        assert conv._should_use_pcs() is False


class TestDegenerate:
    """Empty and tiny axes must not raise."""

    def test_empty_axis_is_zero_sized(self):
        conv = _converter((10, 10, 1), axis=np.array([]))
        assert conv._estimate_output_size_gb() == 0.0
        assert conv._should_use_pcs() is False

    def test_single_bin_axis(self):
        conv = _converter((10, 10, 1), axis=np.array([500.0]))
        expected = 100 * 1 * 4 / (1024**3)
        assert conv._estimate_output_size_gb() == pytest.approx(expected, rel=1e-9)
