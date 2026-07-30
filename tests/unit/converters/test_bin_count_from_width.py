# tests/unit/converters/test_bin_count_from_width.py
"""Tests for deriving a bin count from a target width at a reference m/z.

``_calculate_bins_from_width`` turns "I want bins this wide at this m/z"
into a bin count, and ``CommonAxisBuilder.build_physics_axis`` then
distributes that many bins according to the analyser's spacing law. The
two have to agree: if the count is derived with a different law than the
axis is generated with, the axis is still valid but its realized width at
``reference_mz`` is not the width that was asked for.

These tests close that loop -- derive the count, build the axis, and
measure the width actually realized near ``reference_mz``.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from thyra.converters.spatialdata.base_spatialdata_converter import (
    BaseSpatialDataConverter,
)
from thyra.resampling.common_axis import CommonAxisBuilder
from thyra.resampling.types import AxisType

ALL_AXIS_TYPES = [
    AxisType.CONSTANT,
    AxisType.LINEAR_TOF,
    AxisType.REFLECTOR_TOF,
    AxisType.ORBITRAP,
    AxisType.FTICR,
]


def _bin_count(axis_type, width_at_mz, reference_mz, min_mz, max_mz) -> int:
    """Derive a bin count exactly the way the converter does.

    ``_calculate_bins_from_width`` reads only ``_width_at_mz`` and
    ``_reference_mz`` off ``self``, so a stub carrying those two
    attributes exercises it without standing up a whole converter (which
    would need a reader and a writable output path).
    """
    stub = SimpleNamespace(_width_at_mz=width_at_mz, _reference_mz=reference_mz)
    return BaseSpatialDataConverter._calculate_bins_from_width(
        stub, min_mz, max_mz, axis_type
    )


def _realized_width_at(axis_type, width_at_mz, reference_mz, min_mz, max_mz) -> float:
    """Build the axis the converter would build and measure its width.

    Returns the bin width of the bin whose centre sits closest to
    ``reference_mz``.
    """
    bins = _bin_count(axis_type, width_at_mz, reference_mz, min_mz, max_mz)
    axis = CommonAxisBuilder().build_physics_axis(
        min_mz=min_mz,
        max_mz=max_mz,
        num_bins=bins,
        axis_type=axis_type,
        reference_mz=reference_mz,
        reference_width=width_at_mz,
    )
    mz = axis.mz_values
    widths = np.diff(mz)
    centres = (mz[:-1] + mz[1:]) / 2
    return float(widths[int(np.argmin(np.abs(centres - reference_mz)))])


class TestRealizedWidthMatchesRequest:
    """The width realized at reference_mz must be the width requested."""

    @pytest.mark.parametrize("axis_type", ALL_AXIS_TYPES)
    def test_default_width_and_reference(self, axis_type):
        realized = _realized_width_at(axis_type, 0.005, 1000.0, 100.0, 2000.0)

        assert realized == pytest.approx(0.005, rel=0.03), (
            f"{axis_type.value}: realized width {realized*1000:.4f} mDa is "
            f"not within 3% of the requested 5.0000 mDa"
        )

    @pytest.mark.parametrize("axis_type", ALL_AXIS_TYPES)
    def test_fticr_style_narrow_range(self, axis_type):
        """FT-ICR data is typically a narrow high-mass window."""
        realized = _realized_width_at(axis_type, 0.001, 700.0, 400.0, 1000.0)

        assert realized == pytest.approx(0.001, rel=0.03), (
            f"{axis_type.value}: realized width {realized*1000:.4f} mDa is "
            f"not within 3% of the requested 1.0000 mDa"
        )

    @pytest.mark.parametrize("axis_type", ALL_AXIS_TYPES)
    @pytest.mark.parametrize("reference_mz", [200.0, 500.0, 1500.0])
    def test_anchoring_follows_the_reference_mz(self, axis_type, reference_mz):
        """Moving reference_mz must move where the width is anchored."""
        realized = _realized_width_at(axis_type, 0.01, reference_mz, 150.0, 1800.0)

        assert realized == pytest.approx(0.01, rel=0.03), (
            f"{axis_type.value} at reference m/z {reference_mz}: realized "
            f"width {realized*1000:.4f} mDa is not within 3% of the "
            f"requested 10.0000 mDa"
        )


class TestBinCount:
    """Bin-count behaviour that does not depend on the generated axis."""

    def test_fticr_uses_the_inverse_mz_formula(self):
        """bins = (1/min - 1/max) * (reference_mz^2 / width).

        Guards against the FT-ICR axis silently falling through to the
        uniform ``(max - min) / width`` count again.
        """
        min_mz, max_mz, width, ref = 100.0, 2000.0, 0.005, 1000.0
        expected = int((1 / min_mz - 1 / max_mz) * (ref**2 / width))

        assert _bin_count(AxisType.FTICR, width, ref, min_mz, max_mz) == expected

    def test_fticr_differs_from_the_uniform_count(self):
        """The regression this guards had FT-ICR share the constant count."""
        min_mz, max_mz, width, ref = 100.0, 2000.0, 0.005, 1000.0

        fticr = _bin_count(AxisType.FTICR, width, ref, min_mz, max_mz)
        constant = _bin_count(AxisType.CONSTANT, width, ref, min_mz, max_mz)

        assert fticr != constant

    def test_orbitrap_includes_the_factor_of_two(self):
        """bins = 2 * (1/sqrt(min) - 1/sqrt(max)) * (reference_mz^1.5 / width)."""
        min_mz, max_mz, width, ref = 100.0, 2000.0, 0.005, 1000.0
        expected = int(
            2 * (1 / np.sqrt(min_mz) - 1 / np.sqrt(max_mz)) * (ref**1.5 / width)
        )

        assert _bin_count(AxisType.ORBITRAP, width, ref, min_mz, max_mz) == expected

    @pytest.mark.parametrize("axis_type", ALL_AXIS_TYPES)
    def test_bin_count_is_floored_at_100(self, axis_type):
        """A width wider than the whole range must still give a usable axis."""
        assert _bin_count(axis_type, 5000.0, 1000.0, 100.0, 2000.0) == 100

    @pytest.mark.parametrize("axis_type", ALL_AXIS_TYPES)
    def test_narrower_width_gives_more_bins(self, axis_type):
        wide = _bin_count(axis_type, 0.01, 1000.0, 100.0, 2000.0)
        narrow = _bin_count(axis_type, 0.001, 1000.0, 100.0, 2000.0)

        assert narrow > wide


class TestDefaultWidths:
    """With no explicit width, axis-type-specific defaults apply."""

    def test_linear_tof_defaults_to_17_mda_at_300(self):
        """The SCiLS Lab convention for FlexImaging data."""
        stub = SimpleNamespace(_width_at_mz=None, _reference_mz=1000.0)
        bins = BaseSpatialDataConverter._calculate_bins_from_width(
            stub, 100.0, 2000.0, AxisType.LINEAR_TOF
        )
        expected = _bin_count(AxisType.LINEAR_TOF, 0.017, 300.0, 100.0, 2000.0)

        assert bins == expected

    @pytest.mark.parametrize(
        "axis_type",
        [
            AxisType.CONSTANT,
            AxisType.REFLECTOR_TOF,
            AxisType.ORBITRAP,
            AxisType.FTICR,
        ],
    )
    def test_other_axis_types_default_to_5_mda_at_1000(self, axis_type):
        stub = SimpleNamespace(_width_at_mz=None, _reference_mz=300.0)
        bins = BaseSpatialDataConverter._calculate_bins_from_width(
            stub, 100.0, 2000.0, axis_type
        )
        expected = _bin_count(axis_type, 0.005, 1000.0, 100.0, 2000.0)

        assert bins == expected
