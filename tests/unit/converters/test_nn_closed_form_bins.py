# tests/unit/converters/test_nn_closed_form_bins.py
"""The computed bin index must be invisible: same bins as the search, always.

``_nn_map_to_bins`` finds each peak's nearest bin by ``np.searchsorted``
unless it is handed the axis's ``AxisLinearisation``, in which case it
rounds the peak's position in the generator's own coordinate and repairs
by one bin either side (design decision D21). The two routes must agree
element for element on every input the converter can produce, because
the nearest-neighbour path resamples every spectrum twice and pass 2 has
to reproduce pass 1 exactly, and because a wrong bin is a plausible-looking
spectrum that no total-ion check can see.

What is pinned here:

* bit-identity on axis points, exact midpoints, their float neighbours,
  the half-bin skirt beyond either end, and random values, for every axis
  law and for both sorted and unsorted probes;
* the tie rule -- a peak exactly between two bins goes right -- through
  both routes;
* the build-time guard: an axis that is not strictly ascending, or that
  deviates from its linearisation by a quarter bin or more, falls back to
  the search; and
* the converter's resample path, cached and generic, returns the same
  bins and sums with the linearisation as without it.
"""

from __future__ import annotations

from types import MethodType, SimpleNamespace

import numpy as np
import pytest

from thyra.converters.spatialdata.base_spatialdata_converter import (
    NN_LINEARISATION_MARGIN,
    BaseSpatialDataConverter,
    _nn_map_to_bins,
    _usable_linearisation,
)
from thyra.resampling.mass_axis import (
    FTICRAxisGenerator,
    LinearAxisGenerator,
    LinearTOFAxisGenerator,
    OrbitrapAxisGenerator,
    ReflectorTOFAxisGenerator,
    TOFAxisGenerator,
)
from thyra.resampling.mass_axis.tof_generator import PHI_NANOTOF_LAW, TIMSTOF_TOF_LAW
from thyra.resampling.types import AxisLinearisation

GENERATORS = [
    pytest.param(LinearAxisGenerator(), id="constant"),
    pytest.param(LinearTOFAxisGenerator(), id="linear_tof"),
    pytest.param(ReflectorTOFAxisGenerator(), id="reflector_tof"),
    pytest.param(OrbitrapAxisGenerator(), id="orbitrap"),
    pytest.param(FTICRAxisGenerator(), id="fticr"),
    pytest.param(TOFAxisGenerator(*TIMSTOF_TOF_LAW), id="tof-timstof"),
    pytest.param(TOFAxisGenerator(*PHI_NANOTOF_LAW), id="tof-phi"),
]


def _probes(axis, min_mz, max_mz, rng, n_random=50_000):
    """Every value class the mapping has to get right, sorted."""
    mids = (axis[:-1] + axis[1:]) / 2
    parts = [
        axis,
        mids,
        np.nextafter(mids, np.inf),
        np.nextafter(mids, -np.inf),
        np.nextafter(axis, np.inf),
        np.nextafter(axis, -np.inf),
        rng.uniform(min_mz, max_mz, n_random),
        rng.uniform(min_mz, axis[0], 500),  # the skirt below the first centre
        rng.uniform(axis[-1], max_mz, 500),  # and above the last
        np.array([min_mz, max_mz, axis[0], axis[-1]]),
    ]
    return np.sort(np.concatenate(parts))


@pytest.mark.parametrize("generator", GENERATORS)
@pytest.mark.parametrize(
    "min_mz, max_mz, target_bins",
    [(50.0, 1000.0, 2_000), (100.0, 2000.0, 50_000), (540.0, 600.0, 190_000)],
)
class TestComputedBinsMatchTheSearch:
    def test_bit_identical_on_every_value_class(
        self, generator, min_mz, max_mz, target_bins
    ):
        rng = np.random.default_rng(295)
        mass_axis = generator.generate_axis(min_mz, max_mz, target_bins)
        axis = mass_axis.mz_values.astype(np.float64)
        probes = _probes(axis, min_mz, max_mz, rng)
        searched = _nn_map_to_bins(axis, probes)
        computed = _nn_map_to_bins(axis, probes, mass_axis.linearisation)
        np.testing.assert_array_equal(computed, searched)

    def test_unsorted_probes_too(self, generator, min_mz, max_mz, target_bins):
        rng = np.random.default_rng(1)
        mass_axis = generator.generate_axis(min_mz, max_mz, target_bins)
        axis = mass_axis.mz_values.astype(np.float64)
        probes = _probes(axis, min_mz, max_mz, rng, n_random=5_000)
        rng.shuffle(probes)
        np.testing.assert_array_equal(
            _nn_map_to_bins(axis, probes, mass_axis.linearisation),
            _nn_map_to_bins(axis, probes),
        )


class TestTieRule:
    def test_a_peak_exactly_between_two_bins_goes_right_on_both_routes(self):
        mass_axis = LinearAxisGenerator().generate_axis(0.0, 4.0, 5)
        axis = mass_axis.mz_values
        probes = np.array([0.5, 1.5, 2.5, 3.5])
        np.testing.assert_array_equal(_nn_map_to_bins(axis, probes), [1, 2, 3, 4])
        np.testing.assert_array_equal(
            _nn_map_to_bins(axis, probes, mass_axis.linearisation), [1, 2, 3, 4]
        )

    def test_far_outside_values_take_the_edge_bins(self):
        # Not an input the converter produces (the range filter runs
        # first), but the two routes must still agree if one arrives.
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 1000.0, 100)
        axis = mass_axis.mz_values
        probes = np.array([1.0, 50.0, 2000.0, 1e6])
        np.testing.assert_array_equal(
            _nn_map_to_bins(axis, probes, mass_axis.linearisation),
            _nn_map_to_bins(axis, probes),
        )


class TestBuildTimeGuard:
    def test_a_generated_axis_passes(self):
        mass_axis = OrbitrapAxisGenerator().generate_axis(100.0, 2000.0, 50_000)
        axis = mass_axis.mz_values
        assert (
            _usable_linearisation(
                mass_axis.linearisation, axis, float(np.diff(axis).min())
            )
            is mass_axis.linearisation
        )

    def test_no_linearisation_means_no_closed_form(self):
        axis = np.linspace(100.0, 1000.0, 100)
        assert _usable_linearisation(None, axis, 1.0) is None

    def test_a_duplicated_axis_value_keeps_the_search(self):
        # searchsorted maps a duplicated value to its first occurrence;
        # the repair would not, so a non-ascending axis must not compute.
        axis = np.array([100.0, 200.0, 200.0, 300.0])
        lin = AxisLinearisation(lambda m: np.asarray(m, float), 100.0, 100.0)
        assert _usable_linearisation(lin, axis, 0.0) is None

    def test_a_deviation_at_the_margin_keeps_the_search(self):
        axis = np.linspace(100.0, 1000.0, 100)
        # A linearisation that is wrong by exactly the margin.
        step = axis[1] - axis[0]
        lin = AxisLinearisation(
            lambda m: np.asarray(m, float), 100.0 - NN_LINEARISATION_MARGIN * step, step
        )
        assert lin.deviation(axis) == pytest.approx(NN_LINEARISATION_MARGIN)
        assert _usable_linearisation(lin, axis, step) is None

    def test_an_absurd_axis_is_refused_by_the_margin_not_by_luck(self):
        # Two bins over six decades: the deviation is within 2e-6 of the
        # half-bin cliff. The proof still holds there; the margin is what
        # keeps float rounding from ever being asked to.
        mass_axis = FTICRAxisGenerator().generate_axis(1.0, 1e6, 2)
        axis = mass_axis.mz_values
        assert mass_axis.linearisation.deviation(axis) > NN_LINEARISATION_MARGIN
        assert _usable_linearisation(mass_axis.linearisation, axis, 1.0) is None


def _stub(mass_axis, linearisation, cache_enabled):
    """A converter stand-in with just enough state for the resample path."""
    stub = SimpleNamespace(
        _common_mass_axis=np.asarray(mass_axis.mz_values, float),
        _axis_range=(mass_axis.min_mz, mass_axis.max_mz),
        _nn_linearisation=linearisation,
        _nn_shared_cache=None if cache_enabled else False,
        _nn_cache_misses=0,
        _out_of_range_peaks=0,
        _out_of_range_warned=True,
    )
    for name in (
        "_nearest_neighbor_resample",
        "_nn_resample_via_cache",
        "_build_nn_shared_cache",
        "_count_out_of_range",
    ):
        setattr(stub, name, MethodType(getattr(BaseSpatialDataConverter, name), stub))
    return stub


class TestConverterPathsAgree:
    @pytest.mark.parametrize("cache_enabled", [False, True], ids=["generic", "cached"])
    def test_resample_returns_the_same_bins_and_sums(self, cache_enabled):
        rng = np.random.default_rng(7)
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 2000.0, 100_000)
        searching = _stub(mass_axis, None, cache_enabled)
        computing = _stub(mass_axis, mass_axis.linearisation, cache_enabled)
        mzs = np.sort(rng.uniform(90.0, 2010.0, 3_000))  # some out of range
        for _ in range(3):
            ints = rng.lognormal(np.log(50), 1.0, mzs.size)
            ints[rng.random(mzs.size) < 0.3] = 0.0
            a_bins, a_sums = searching._nearest_neighbor_resample(mzs, ints)
            b_bins, b_sums = computing._nearest_neighbor_resample(mzs, ints)
            np.testing.assert_array_equal(b_bins, a_bins)
            np.testing.assert_array_equal(b_sums, a_sums)
        assert computing._out_of_range_peaks == searching._out_of_range_peaks

    def test_a_stub_without_the_attribute_still_searches(self):
        # The unbound-call harnesses predate the attribute; they must keep
        # working, on the search.
        mass_axis = LinearAxisGenerator().generate_axis(100.0, 1000.0, 1_000)
        stub = _stub(mass_axis, None, False)
        del stub._nn_linearisation
        mzs = np.array([100.0, 500.25, 999.9])
        bins, sums = stub._nearest_neighbor_resample(mzs, np.ones(3))
        np.testing.assert_array_equal(bins, _nn_map_to_bins(mass_axis.mz_values, mzs))
