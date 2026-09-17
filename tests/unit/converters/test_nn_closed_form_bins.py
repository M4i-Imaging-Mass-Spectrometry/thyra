# tests/unit/converters/test_nn_closed_form_bins.py
"""The computed bin index must be invisible: same bins as the search, always.

``nn_map_to_bins`` finds each peak's nearest bin by ``np.searchsorted``
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
* the strategy's resample path, cached and generic, returns the same
  bins and sums with the linearisation as without it.
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.resampling.binning import (
    NN_LINEARISATION_MARGIN,
    nn_map_to_bins,
    usable_linearisation,
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
from thyra.resampling.strategies import NearestNeighborStrategy
from thyra.resampling.types import AxisLinearisation, AxisType

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
        searched = nn_map_to_bins(axis, probes)
        computed = nn_map_to_bins(axis, probes, mass_axis.linearisation)
        np.testing.assert_array_equal(computed, searched)

    def test_unsorted_probes_too(self, generator, min_mz, max_mz, target_bins):
        rng = np.random.default_rng(1)
        mass_axis = generator.generate_axis(min_mz, max_mz, target_bins)
        axis = mass_axis.mz_values.astype(np.float64)
        probes = _probes(axis, min_mz, max_mz, rng, n_random=5_000)
        rng.shuffle(probes)
        np.testing.assert_array_equal(
            nn_map_to_bins(axis, probes, mass_axis.linearisation),
            nn_map_to_bins(axis, probes),
        )


class TestTieRule:
    def test_a_peak_exactly_between_two_bins_goes_right_on_both_routes(self):
        mass_axis = LinearAxisGenerator().generate_axis(0.0, 4.0, 5)
        axis = mass_axis.mz_values
        probes = np.array([0.5, 1.5, 2.5, 3.5])
        np.testing.assert_array_equal(nn_map_to_bins(axis, probes), [1, 2, 3, 4])
        np.testing.assert_array_equal(
            nn_map_to_bins(axis, probes, mass_axis.linearisation), [1, 2, 3, 4]
        )

    def test_far_outside_values_take_the_edge_bins(self):
        # Not an input the converter produces (the range filter runs
        # first), but the two routes must still agree if one arrives.
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 1000.0, 100)
        axis = mass_axis.mz_values
        probes = np.array([1.0, 50.0, 2000.0, 1e6])
        np.testing.assert_array_equal(
            nn_map_to_bins(axis, probes, mass_axis.linearisation),
            nn_map_to_bins(axis, probes),
        )


class TestBuildTimeGuard:
    def test_a_generated_axis_passes(self):
        mass_axis = OrbitrapAxisGenerator().generate_axis(100.0, 2000.0, 50_000)
        axis = mass_axis.mz_values
        assert (
            usable_linearisation(
                mass_axis.linearisation, axis, float(np.diff(axis).min())
            )
            is mass_axis.linearisation
        )

    def test_no_linearisation_means_no_closed_form(self):
        axis = np.linspace(100.0, 1000.0, 100)
        assert usable_linearisation(None, axis, 1.0) is None

    def test_a_duplicated_axis_value_keeps_the_search(self):
        # searchsorted maps a duplicated value to its first occurrence;
        # the repair would not, so a non-ascending axis must not compute.
        axis = np.array([100.0, 200.0, 200.0, 300.0])
        lin = AxisLinearisation(lambda m: np.asarray(m, float), 100.0, 100.0)
        assert usable_linearisation(lin, axis, 0.0) is None

    def test_a_deviation_at_the_margin_keeps_the_search(self):
        axis = np.linspace(100.0, 1000.0, 100)
        # A linearisation that is wrong by exactly the margin.
        step = axis[1] - axis[0]
        lin = AxisLinearisation(
            lambda m: np.asarray(m, float), 100.0 - NN_LINEARISATION_MARGIN * step, step
        )
        assert lin.deviation(axis) == pytest.approx(NN_LINEARISATION_MARGIN)
        assert usable_linearisation(lin, axis, step) is None

    def test_an_absurd_axis_is_refused_by_the_margin_not_by_luck(self):
        # Two bins over six decades: the deviation is within 2e-6 of the
        # half-bin cliff. The proof still holds there; the margin is what
        # keeps float rounding from ever being asked to.
        mass_axis = FTICRAxisGenerator().generate_axis(1.0, 1e6, 2)
        axis = mass_axis.mz_values
        assert mass_axis.linearisation.deviation(axis) > NN_LINEARISATION_MARGIN
        assert usable_linearisation(mass_axis.linearisation, axis, 1.0) is None


def _strategy(mass_axis, linearisation, cache_enabled):
    """The nearest-neighbour operator over a generated axis."""
    strategy = NearestNeighborStrategy(
        np.asarray(mass_axis.mz_values, float),
        (mass_axis.min_mz, mass_axis.max_mz),
        linearisation,
    )
    strategy._out_of_range_warned = True
    if not cache_enabled:
        strategy._shared_cache = False
    return strategy


class TestTheStrategysPathsAgree:
    @pytest.mark.parametrize("cache_enabled", [False, True], ids=["generic", "cached"])
    def test_resample_returns_the_same_bins_and_sums(self, cache_enabled):
        rng = np.random.default_rng(7)
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 2000.0, 100_000)
        searching = _strategy(mass_axis, None, cache_enabled)
        computing = _strategy(mass_axis, mass_axis.linearisation, cache_enabled)
        mzs = np.sort(rng.uniform(90.0, 2010.0, 3_000))  # some out of range
        for _ in range(3):
            ints = rng.lognormal(np.log(50), 1.0, mzs.size)
            ints[rng.random(mzs.size) < 0.3] = 0.0
            a_bins, a_sums = searching.resample(mzs, ints)
            b_bins, b_sums = computing.resample(mzs, ints)
            np.testing.assert_array_equal(b_bins, a_bins)
            np.testing.assert_array_equal(b_sums, a_sums)
        assert computing.out_of_range_peaks == searching.out_of_range_peaks

    def test_no_linearisation_still_searches(self):
        # The default: a strategy handed no linearisation places every
        # peak by binary search, and says so.
        mass_axis = LinearAxisGenerator().generate_axis(100.0, 1000.0, 1_000)
        strategy = _strategy(mass_axis, None, False)
        assert strategy.linearisation is None
        mzs = np.array([100.0, 500.25, 999.9])
        bins, _ = strategy.resample(mzs, np.ones(3))
        np.testing.assert_array_equal(bins, nn_map_to_bins(mass_axis.mz_values, mzs))


class TestTheUniformAxisIsLinearised:
    """A constant axis carries its coordinate like every other (issue #348).

    ``CommonAxisBuilder.build_uniform_axis`` laid the bins itself instead
    of asking ``LinearAxisGenerator`` for them, and so reported no
    linearisation. The axis was right; it simply never said what it was
    uniform in, and every placement onto it -- the summed table's
    included -- fell back to the search.
    """

    def test_the_values_are_what_they_always_were(self):
        from thyra.resampling.common_axis import CommonAxisBuilder

        axis = CommonAxisBuilder().build_uniform_axis(90.0, 510.0, 43)

        np.testing.assert_array_equal(axis.mz_values, np.linspace(90.0, 510.0, 43))
        assert axis.axis_type is AxisType.CONSTANT
        assert axis.num_bins == 43

    def test_and_now_it_reports_its_coordinate(self):
        from thyra.resampling.common_axis import CommonAxisBuilder

        axis = CommonAxisBuilder().build_uniform_axis(90.0, 510.0, 43)

        assert axis.linearisation is not None
        # The identity: a constant axis is uniform in m/z itself.
        np.testing.assert_allclose(
            axis.linearisation.positions(axis.mz_values),
            np.arange(43, dtype=np.float64),
            atol=1e-9,
        )
        usable = usable_linearisation(
            axis.linearisation, axis.mz_values, float(np.diff(axis.mz_values).min())
        )
        assert usable is not None

    def test_a_single_bin_axis_has_no_step_to_report(self):
        from thyra.resampling.common_axis import CommonAxisBuilder

        axis = CommonAxisBuilder().build_uniform_axis(90.0, 510.0, 1)

        assert axis.linearisation is None


@pytest.mark.parametrize("fused", [True, False])
def test_the_sibling_tables_are_what_the_search_builds(tmp_path, fused):
    """Every stored table, bit for bit, with the closed form and without it.

    The sinks of the sibling tables were the last callers still searching
    after D21 (issue #348). They map onto the summed table's own axis, so
    they take the same linearisation -- and the tables they build have to
    come out unchanged, which is the constraint the issue set. Checked on
    the store rather than on the mapping: the fused reader hands its
    points indexed and the unfused one flat, so between them this covers
    both ``map_indexed_points_to_axis`` and ``map_points_to_axis``, plus
    the grid's two passes and the MS/MS split.
    """
    spatialdata = pytest.importorskip("spatialdata")

    from tests.unit.converters.test_fused_passes import (
        RESAMPLED,
        FusedStubReader,
        UnfusedStubReader,
    )
    from thyra.converters.spatialdata import base_spatialdata_converter as bsc
    from thyra.converters.spatialdata.streaming_converter import (
        StreamingSpatialDataConverter,
    )
    from thyra.utils.windows_paths import (
        prepare_zarr_output_path,
        prepare_zarr_read_path,
    )

    reader_cls = FusedStubReader if fused else UnfusedStubReader

    def convert(tag, force_search):
        original = bsc.usable_linearisation
        if force_search:
            bsc.usable_linearisation = lambda *a, **k: None
        try:
            out = prepare_zarr_output_path(tmp_path / f"{tag}.zarr", "stub")
            converter = StreamingSpatialDataConverter(
                reader_cls(),
                out,
                dataset_id="stub",
                pixel_size_um=10.0,
                resampling_config=RESAMPLED,
                mobility_grid=True,
            )
            assert converter.convert()
            used = converter._axis_linearisation is not None
            return spatialdata.read_zarr(prepare_zarr_read_path(out)), used
        finally:
            bsc.usable_linearisation = original

    closed, used_closed = convert("closed", force_search=False)
    searched, used_searched = convert("searched", force_search=True)

    # Without this the test would pass by comparing two identical search
    # runs, which is what it did before build_uniform_axis reported one.
    assert used_closed is True
    assert used_searched is False

    assert set(closed.tables) == set(searched.tables)
    assert set(closed.tables) == {"stub_z0", "stub_z0_mobility", "stub_z0_msms"}
    for key in closed.tables:
        a, b = closed.tables[key], searched.tables[key]
        assert a.shape == b.shape, key
        xa, xb = (
            t.X.toarray() if hasattr(t.X, "toarray") else np.asarray(t.X)
            for t in (a, b)
        )
        np.testing.assert_array_equal(xa, xb, err_msg=f"{key}: X differs")
        assert list(a.var.index) == list(b.var.index), key
        for column in a.var.columns:
            np.testing.assert_array_equal(
                a.var[column].to_numpy(),
                b.var[column].to_numpy(),
                err_msg=f"{key}: var[{column}] differs",
            )


def test_an_interpolated_conversion_still_reports_the_axis_coordinate(tmp_path):
    """The linearisation belongs to the axis, not to the method that bins.

    ``tic_preserving`` does not use it to resample -- it interpolates --
    but the sibling tables written beside it still place peaks onto the
    same axis, and they take the linearisation as an argument. Deriving
    it from the strategy instead of holding it beside the axis made this
    ``None`` under the interpolating method, which costs the sinks the
    closed form (D21) for nothing. Pinned so it stays an axis property.
    """
    pytest.importorskip("spatialdata")

    from tests.unit.converters.test_fused_passes import RESAMPLED, FusedStubReader
    from thyra.converters.spatialdata.streaming_converter import (
        StreamingSpatialDataConverter,
    )
    from thyra.resampling.strategies import TICPreservingStrategy
    from thyra.utils.windows_paths import prepare_zarr_output_path

    out = prepare_zarr_output_path(tmp_path / "interpolated.zarr", "stub")
    converter = StreamingSpatialDataConverter(
        FusedStubReader(),
        out,
        dataset_id="stub",
        pixel_size_um=10.0,
        resampling_config={**RESAMPLED, "method": "tic_preserving"},
        mobility_grid=True,
    )
    assert converter.convert()

    assert isinstance(converter._resampler, TICPreservingStrategy)
    assert converter._axis_linearisation is not None
