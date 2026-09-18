"""A detector's bin width reaches the axis, behind the caller's own setting.

``reference_params`` is the one place the width and reference m/z are
decided, for both the bin count and the axis generator. Precedence is the
caller's ``--resample-width-at-mz``, then the width the detected instrument
declared (Waters: 2 mDa for the vendor centroid, 1.3 mDa for the MRT
profile), then the per-axis-type default that applied before any detector
could speak.
"""

from __future__ import annotations

import pytest

from thyra.resampling.axis_planner import (
    axis_name,
    bin_count_for_width,
    reference_params,
)
from thyra.resampling.types import AxisType, ResamplingConfig


def _ask(width=None, ref=1000.0):
    """The caller's side of it: what was asked for, and nothing else."""
    return ResamplingConfig(mass_width_da=width, reference_mz=ref)


class TestPrecedence:
    def test_explicit_width_wins(self):
        assert reference_params(
            _ask(0.01, 500.0), "reflector_tof", detected_width=(0.002, 1000.0)
        ) == (
            0.01,
            500.0,
        )

    def test_detected_width_beats_the_axis_default(self):
        assert reference_params(
            _ask(), "linear_tof", detected_width=(0.0013, 1000.0)
        ) == (
            0.0013,
            1000.0,
        )

    @pytest.mark.parametrize(
        "name, expected",
        [
            ("linear_tof", (0.017, 300.0)),
            ("reflector_tof", (0.005, 1000.0)),
            ("constant", (0.005, 1000.0)),
            ("fticr", (0.005, 1000.0)),
        ],
    )
    def test_axis_defaults_when_nothing_was_declared(self, name, expected):
        assert reference_params(_ask(), name) == expected

    def test_no_detector_spoke_at_all(self):
        """The detected pair is optional; the defaults stand without it."""
        assert reference_params(_ask(), "reflector_tof") == (0.005, 1000.0)


class TestBothConsumersAgree:
    def test_bin_count_and_generator_read_the_same_width(self):
        detected = (0.002, 1000.0)
        assert reference_params(
            _ask(), axis_name(AxisType.REFLECTOR_TOF), detected_width=detected
        ) == (0.002, 1000.0)
        bins = bin_count_for_width(
            _ask(), 100.0, 1000.0, AxisType.REFLECTOR_TOF, detected_width=detected
        )
        # ln(10) * 1000 / 0.002
        assert bins == pytest.approx(1_151_292, rel=1e-3)

    def test_mrt_profile_axis_size_on_the_reference_run(self):
        """1.3 mDa at m/z 1000 over 100-1000 on linear_tof: ~1.05M bins."""
        bins = bin_count_for_width(
            _ask(),
            100.0,
            1000.0,
            AxisType.LINEAR_TOF,
            detected_width=(0.0013, 1000.0),
        )
        assert bins == pytest.approx(1_051_941, rel=1e-3)
