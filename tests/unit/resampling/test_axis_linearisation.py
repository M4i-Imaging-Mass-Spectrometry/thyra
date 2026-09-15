# tests/unit/resampling/test_axis_linearisation.py
"""Every generated axis carries the coordinate it is uniform in, within bounds.

``AxisLinearisation`` is what lets the converter compute a bin index
instead of binary-searching for it (design decision D21). The claim it
rests on: for every axis a generator lays, ``forward(axis[i])`` is within
half a step of ``u0 + i * du``. The physics generators lay a uniform grid
of edges and report arithmetic midpoints in m/z, so the deviation is not
zero -- but a midpoint lies strictly between its two edges and ``forward``
is strictly monotone, so its image lies strictly between the two edge
coordinates, which is the bound. These tests pin

* that ``forward`` is the coordinate the axis was actually laid in,
* the closed form of the deviation where one exists (FT-ICR: exactly
  ``(t - 1) / (2 (t + 1))`` for an edge ratio ``t``), and
* that the bound holds even on axes no instrument would ask for, where
  the deviation approaches its limit.
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.resampling.mass_axis import (
    FTICRAxisGenerator,
    LinearAxisGenerator,
    LinearTOFAxisGenerator,
    OrbitrapAxisGenerator,
    ReflectorTOFAxisGenerator,
    TOFAxisGenerator,
)
from thyra.resampling.mass_axis.tof_generator import (
    MRT_TOF_LAW,
    PHI_NANOTOF_LAW,
    TIMSTOF_TOF_LAW,
)
from thyra.resampling.types import AxisLinearisation

GENERATORS = [
    pytest.param(LinearAxisGenerator(), id="constant"),
    pytest.param(LinearTOFAxisGenerator(), id="linear_tof"),
    pytest.param(ReflectorTOFAxisGenerator(), id="reflector_tof"),
    pytest.param(OrbitrapAxisGenerator(), id="orbitrap"),
    pytest.param(FTICRAxisGenerator(), id="fticr"),
    pytest.param(TOFAxisGenerator(*MRT_TOF_LAW), id="tof-mrt"),
    pytest.param(TOFAxisGenerator(*TIMSTOF_TOF_LAW), id="tof-timstof"),
    pytest.param(TOFAxisGenerator(*PHI_NANOTOF_LAW), id="tof-phi"),
]

#: (min_mz, max_mz, target_bins) an instrument might actually ask for.
REALISTIC = [
    (50.0, 1000.0, 50_000),
    (100.0, 2000.0, 190_000),
    (12.0, 400.0, 100_000),
    (540.0, 600.0, 2_000),
    (300.0, 3000.0, 1_050_000),
]


@pytest.mark.parametrize("generator", GENERATORS)
@pytest.mark.parametrize("min_mz, max_mz, target_bins", REALISTIC)
class TestGeneratedAxesCarryTheirLinearisation:
    def test_forward_is_strictly_increasing_in_mz(
        self, generator, min_mz, max_mz, target_bins
    ):
        axis = generator.generate_axis(min_mz, max_mz, target_bins)
        u = axis.linearisation.forward(axis.mz_values)
        # The coordinate may run either way (1/m runs down); it must be
        # strictly monotone, which is what the deviation bound needs.
        steps = np.diff(u)
        assert np.all(steps > 0) or np.all(steps < 0)

    def test_positions_are_the_indices_within_a_small_margin(
        self, generator, min_mz, max_mz, target_bins
    ):
        axis = generator.generate_axis(min_mz, max_mz, target_bins)
        lin = axis.linearisation
        assert isinstance(lin, AxisLinearisation)
        # The deviation is second order in one bin's relative width, so
        # on a real axis it is far below the converter's 0.25 margin.
        assert lin.deviation(axis.mz_values) < 0.01

    def test_deviation_matches_a_direct_evaluation(
        self, generator, min_mz, max_mz, target_bins
    ):
        axis = generator.generate_axis(min_mz, max_mz, target_bins)
        lin = axis.linearisation
        direct = np.max(
            np.abs(lin.positions(axis.mz_values) - np.arange(axis.num_bins))
        )
        # Chunked evaluation must agree with the one-shot one, including
        # when the chunk is smaller than the axis.
        assert lin.deviation(axis.mz_values, chunk=1000) == pytest.approx(direct)


class TestTheBoundHoldsWhereItIsTightest:
    """The deviation approaches, but never reaches, half a bin."""

    @pytest.mark.parametrize("generator", GENERATORS)
    @pytest.mark.parametrize(
        "min_mz, max_mz, target_bins", [(1.0, 1e6, 2), (1.0, 1e6, 3)]
    )
    def test_below_half_a_bin_even_on_an_absurd_axis(
        self, generator, min_mz, max_mz, target_bins
    ):
        axis = generator.generate_axis(min_mz, max_mz, target_bins)
        assert axis.linearisation.deviation(axis.mz_values) < 0.5

    def test_fticr_deviation_has_its_closed_form(self):
        # For u = 1/m with edges m_i < m_{i+1} and t = m_{i+1} / m_i, the
        # centre's deviation is exactly (t - 1) / (2 (t + 1)); it tends to
        # 1/2 only as t -> infinity.
        gen = FTICRAxisGenerator()
        axis = gen.generate_axis(100.0, 1000.0, 10)
        lin = axis.linearisation
        edges = 1.0 / np.linspace(1 / 100.0, 1 / 1000.0, 11)
        t = edges[1:] / edges[:-1]
        expected = (t - 1) / (2 * (t + 1))
        got = np.abs(lin.positions(axis.mz_values) - np.arange(10))
        np.testing.assert_allclose(got, expected, rtol=1e-9, atol=1e-12)

    def test_linear_tof_deviation_is_bounded_by_rms_minus_mean(self):
        # sqrt of the midpoint is the RMS of the two edge coordinates; RMS
        # exceeds the mean by at most (1/sqrt(2) - 1/2) of the step.
        gen = LinearTOFAxisGenerator()
        axis = gen.generate_axis(1.0, 1e6, 2)
        assert (
            axis.linearisation.deviation(axis.mz_values) < 1 / np.sqrt(2) - 0.5 + 1e-12
        )


class TestConstantSpacing:
    def test_the_linspace_is_its_own_linearisation(self):
        axis = LinearAxisGenerator().generate_axis(100.0, 1000.0, 190_000)
        lin = axis.linearisation
        assert lin.u0 == 100.0
        assert lin.du == pytest.approx((1000.0 - 100.0) / (190_000 - 1))
        # Only float rounding separates the two.
        assert lin.deviation(axis.mz_values) < 1e-6

    def test_a_single_point_has_no_step(self):
        axis = LinearAxisGenerator().generate_axis(100.0, 1000.0, 1)
        assert axis.linearisation is None
