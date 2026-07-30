# tests/unit/resampling/test_mass_axis_generators.py
"""Contract tests for the physics-based mass axis generators.

Every generator feeds ``CommonAxisBuilder.build_physics_axis``, whose
output is assigned directly to the converter's common mass axis. Two
things therefore have to hold for all of them:

  * the axis is in ascending m/z order (``np.searchsorted`` binning and
    the stored ``var["mz"]`` column both assume increasing m/z), and
  * bin width scales with m/z as the analyser's physics dictates.

Anchoring that spacing law so the width at ``reference_mz`` is the width
the caller asked for is the bin-count calculation's job, and is covered
by ``tests/unit/converters/test_bin_count_from_width.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.resampling.common_axis import CommonAxisBuilder
from thyra.resampling.types import AxisType

# Every axis type reachable through build_physics_axis, with the exponent
# p in bin_width ∝ (m/z)^p that defines its spacing law.
AXIS_TYPES = [
    (AxisType.CONSTANT, 0.0),
    (AxisType.LINEAR_TOF, 0.5),
    (AxisType.REFLECTOR_TOF, 1.0),
    (AxisType.ORBITRAP, 1.5),
    (AxisType.FTICR, 2.0),
]

MIN_MZ = 100.0
MAX_MZ = 2000.0
REFERENCE_MZ = 1000.0
WIDTH_AT_REF = 0.005  # 5 mDa, the converter's default for most axis types


def _build(axis_type, num_bins, reference_mz=REFERENCE_MZ, width=WIDTH_AT_REF):
    return CommonAxisBuilder().build_physics_axis(
        min_mz=MIN_MZ,
        max_mz=MAX_MZ,
        num_bins=num_bins,
        axis_type=axis_type,
        reference_mz=reference_mz,
        reference_width=width,
    )


class TestAxisOrdering:
    """Generated axes must be ascending in m/z."""

    @pytest.mark.parametrize("axis_type,_p", AXIS_TYPES)
    def test_axis_is_ascending(self, axis_type, _p):
        mz = _build(axis_type, 5000).mz_values

        assert np.all(np.diff(mz) > 0), (
            f"{axis_type.value} axis is not strictly ascending; "
            f"first={mz[0]}, last={mz[-1]}"
        )

    @pytest.mark.parametrize("axis_type,_p", AXIS_TYPES)
    def test_axis_spans_requested_range(self, axis_type, _p):
        axis = _build(axis_type, 5000)
        mz = axis.mz_values

        # Bin-centre conventions mean the ends can sit up to one bin
        # inside the requested range, but not outside it and not far in.
        span = MAX_MZ - MIN_MZ
        assert MIN_MZ - 0.01 <= mz[0] < MIN_MZ + 0.01 * span
        assert MAX_MZ - 0.01 * span < mz[-1] <= MAX_MZ + 0.01
        assert axis.min_mz == pytest.approx(float(mz[0]))
        assert axis.max_mz == pytest.approx(float(mz[-1]))

    @pytest.mark.parametrize("axis_type,_p", AXIS_TYPES)
    def test_bin_widths_are_positive(self, axis_type, _p):
        widths = np.diff(_build(axis_type, 5000).mz_values)

        assert np.all(widths > 0)


class TestSpacingLaw:
    """Bin width must scale with m/z as the analyser physics dictates."""

    @pytest.mark.parametrize("axis_type,exponent", AXIS_TYPES)
    def test_width_ratio_matches_exponent(self, axis_type, exponent):
        """width(m2) / width(m1) must be (m2/m1)^p."""
        mz = _build(axis_type, 200000).mz_values
        widths = np.diff(mz)
        centres = (mz[:-1] + mz[1:]) / 2

        low = int(np.argmin(np.abs(centres - 300.0)))
        high = int(np.argmin(np.abs(centres - 1500.0)))

        observed = widths[high] / widths[low]
        expected = (centres[high] / centres[low]) ** exponent

        assert observed == pytest.approx(expected, rel=0.02), (
            f"{axis_type.value}: width ratio {observed:.4f} does not match "
            f"the expected (m2/m1)^{exponent} = {expected:.4f}"
        )
