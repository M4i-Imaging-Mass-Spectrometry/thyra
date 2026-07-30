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
from thyra.resampling.mass_axis import (
    BaseAxisGenerator,
    FTICRAxisGenerator,
    LinearAxisGenerator,
    LinearTOFAxisGenerator,
    OrbitrapAxisGenerator,
    ReflectorTOFAxisGenerator,
)
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

GENERATORS = [
    LinearAxisGenerator,
    LinearTOFAxisGenerator,
    ReflectorTOFAxisGenerator,
    OrbitrapAxisGenerator,
    FTICRAxisGenerator,
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


class TestGeneratorInterface:
    """All generators share one interface, and it is the one in use."""

    @pytest.mark.parametrize("generator_cls", GENERATORS, ids=lambda c: c.__name__)
    def test_is_concrete(self, generator_cls):
        """Nothing abstract may be left unimplemented."""
        assert isinstance(generator_cls(), BaseAxisGenerator)

    @pytest.mark.parametrize("generator_cls", GENERATORS, ids=lambda c: c.__name__)
    def test_accepts_the_unified_signature(self, generator_cls):
        """Every generator takes the reference parameters, used or not.

        LinearAxisGenerator used to take three arguments while the physics
        generators took five, forcing build_physics_axis to branch on the
        axis type before calling.
        """
        axis = generator_cls().generate_axis(MIN_MZ, MAX_MZ, 1000, REFERENCE_MZ, 0.01)

        assert len(axis.mz_values) > 0

    @pytest.mark.parametrize("generator_cls", GENERATORS, ids=lambda c: c.__name__)
    def test_reported_axis_type_matches_what_is_produced(self, generator_cls):
        generator = generator_cls()

        axis = generator.generate_axis(MIN_MZ, MAX_MZ, 1000)

        assert axis.axis_type == generator.get_axis_type()

    @pytest.mark.parametrize("removed", ["generate_axis_bins", "generate_axis_width"])
    @pytest.mark.parametrize("generator_cls", GENERATORS, ids=lambda c: c.__name__)
    def test_dead_parallel_api_stays_gone(self, generator_cls, removed):
        """Guard against the unused width/bins pair coming back.

        These were a second, never-called route from a target width to an
        axis. Nothing reached them, so they drifted: the FT-ICR variant
        computed its step as ``k / width_da`` instead of ``k``, two of them
        returned descending axes, and the Linear TOF one assigned its step
        twice. Width-driven sizing goes through
        ``_calculate_bins_from_width`` instead.
        """
        assert not hasattr(generator_cls(), removed)


class TestWidthPrediction:
    """calculate_width_at_mz must agree with the axis actually built.

    This cross-checks the two independent expressions of each spacing law
    against each other: the closed-form width prediction, and the width the
    generated axis realizes once the bin count is derived from the physics.
    """

    PREDICTORS = [
        (LinearTOFAxisGenerator, AxisType.LINEAR_TOF),
        (ReflectorTOFAxisGenerator, AxisType.REFLECTOR_TOF),
        (OrbitrapAxisGenerator, AxisType.ORBITRAP),
        (FTICRAxisGenerator, AxisType.FTICR),
    ]

    @pytest.mark.parametrize(
        "generator_cls,axis_type", PREDICTORS, ids=lambda v: getattr(v, "__name__", v)
    )
    @pytest.mark.parametrize("probe_mz", [300.0, 700.0, 1500.0])
    def test_prediction_matches_realized_width(
        self, generator_cls, axis_type, probe_mz
    ):
        from types import SimpleNamespace

        from thyra.converters.spatialdata.base_spatialdata_converter import (
            BaseSpatialDataConverter,
        )

        stub = SimpleNamespace(_width_at_mz=WIDTH_AT_REF, _reference_mz=REFERENCE_MZ)
        bins = BaseSpatialDataConverter._calculate_bins_from_width(
            stub, MIN_MZ, MAX_MZ, axis_type
        )
        mz = _build(axis_type, bins).mz_values
        widths = np.diff(mz)
        centres = (mz[:-1] + mz[1:]) / 2
        realized = widths[int(np.argmin(np.abs(centres - probe_mz)))]

        predicted = generator_cls().calculate_width_at_mz(
            probe_mz, REFERENCE_MZ, WIDTH_AT_REF
        )

        assert realized == pytest.approx(predicted, rel=0.03), (
            f"{axis_type.value} at m/z {probe_mz}: axis realizes "
            f"{realized*1000:.4f} mDa but calculate_width_at_mz predicts "
            f"{predicted*1000:.4f} mDa"
        )
