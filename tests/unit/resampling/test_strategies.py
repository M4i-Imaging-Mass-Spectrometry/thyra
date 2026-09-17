# tests/unit/resampling/test_strategies.py
"""The contract the converter runs a resampling strategy through.

One spectrum in, the bins it fills out: unique ascending indices, no
explicit zeros, and nothing of the axis's length. Both methods answer the
same call, which is what lets the converter hold one of them and not care
which. Pinned here rather than through a conversion because a conversion
would prove it only for the method its fixture happens to select.

The per-method behaviour each class inherited from the converter's private
methods stays where it was measured -- the cache in
``tests/unit/converters/test_nn_shared_axis_cache.py``, the closed-form bin
index in ``test_nn_closed_form_bins.py``, the sparse evaluation in
``test_tic_preserving_sparse.py``. What is new here is the contract itself,
the factory, and the out-of-range accounting now that it belongs to the
strategy rather than to the converter.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from thyra.resampling.binning import nn_map_to_bins
from thyra.resampling.mass_axis import ReflectorTOFAxisGenerator
from thyra.resampling.strategies import (
    NearestNeighborStrategy,
    TICPreservingStrategy,
    build_strategy,
)
from thyra.resampling.types import ResamplingMethod

LOGGER = "thyra.resampling.strategies.base"


@pytest.fixture
def warnings_from_the_strategy(caplog):
    """``caplog``, with the ``thyra`` logger's propagation put back.

    ``setup_logging`` sets ``propagate = False`` on it and any test that
    has called it leaves that behind, so records reach pytest's handler
    only by accident of test order -- which is why
    ``tests/unit/converters/test_out_of_range_peaks.py`` attaches a
    handler of its own instead.
    """
    logger = logging.getLogger("thyra")
    previous = logger.propagate
    logger.propagate = True
    try:
        with caplog.at_level(logging.WARNING, logger=LOGGER):
            yield caplog
    finally:
        logger.propagate = previous


def _axis(n=5_000, lo=100.0, hi=1000.0):
    return np.linspace(lo, hi, n)


def _spectrum(seed=0, n=400, lo=100.0, hi=1000.0):
    rng = np.random.default_rng(seed)
    mzs = np.sort(rng.uniform(lo, hi, n))
    intensities = rng.lognormal(np.log(50), 1.0, n)
    return mzs, intensities


def _both(axis, axis_range=None):
    return [
        NearestNeighborStrategy(axis, axis_range),
        TICPreservingStrategy(axis, axis_range),
    ]


class TestTheSparseContract:
    """What ``_process_spectrum`` hands the CSC assembly, for both methods."""

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_indices_are_unique_and_ascending(self, seed):
        axis = _axis()
        mzs, intensities = _spectrum(seed)

        for strategy in _both(axis):
            indices, _ = strategy.resample(mzs, intensities)
            assert np.all(np.diff(indices) > 0), type(strategy).__name__

    def test_no_value_is_zero(self):
        axis = _axis()
        mzs, intensities = _spectrum(3)
        # A third of the spectrum is an explicit zero, as a profile
        # source stores it; none of them may reach the store.
        intensities[np.random.default_rng(3).random(mzs.size) < 0.33] = 0.0

        for strategy in _both(axis):
            _, values = strategy.resample(mzs, intensities)
            assert np.all(values != 0), type(strategy).__name__

    def test_indices_are_inside_the_axis(self):
        axis = _axis()
        # Well outside on both sides, so the range filter has work to do.
        mzs = np.array([10.0, 250.0, 500.0, 5_000.0])
        intensities = np.array([9.0, 1.0, 2.0, 9.0])

        for strategy in _both(axis):
            indices, _ = strategy.resample(mzs, intensities)
            assert indices.size == 0 or (
                indices.min() >= 0 and indices.max() < axis.size
            )

    def test_the_dense_form_holds_the_same_total(self):
        axis = _axis()
        mzs, intensities = _spectrum(4)

        for strategy in _both(axis):
            indices, values = strategy.resample(mzs, intensities)
            dense = strategy.to_dense(indices, values)
            assert dense.shape == (axis.size,)
            assert dense.sum() == pytest.approx(values.sum(), rel=1e-12)
            assert np.count_nonzero(dense) == values.size

    def test_an_empty_spectrum_gives_empty_arrays(self):
        axis = _axis()
        empty = np.array([], dtype=np.float64)

        for strategy in _both(axis):
            indices, values = strategy.resample(empty, empty)
            assert indices.size == 0 and values.size == 0


class TestTheNearestNeighbourCache:
    def test_it_disables_itself_after_five_misses(self):
        rng = np.random.default_rng(6)
        strategy = NearestNeighborStrategy(_axis(3_000), None)

        for _ in range(8):
            n = int(rng.integers(50, 200))
            mzs = np.sort(rng.uniform(100.0, 1000.0, n))
            strategy.resample(mzs, np.ones(n))

        assert strategy._shared_cache is False

    def test_a_shared_array_keeps_hitting(self):
        strategy = NearestNeighborStrategy(_axis(3_000), None)
        mzs, intensities = _spectrum(7)

        first = strategy.resample(mzs, intensities)
        again = strategy.resample(mzs, intensities)

        np.testing.assert_array_equal(first[0], again[0])
        assert strategy._cache_misses == 0


class TestTheLinearisationIsInvisible:
    """A linearised axis and the same axis without one give the same bins."""

    def test_same_bins_and_sums(self):
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 2000.0, 100_000)
        axis = np.asarray(mass_axis.mz_values, float)
        axis_range = (mass_axis.min_mz, mass_axis.max_mz)
        computing = NearestNeighborStrategy(axis, axis_range, mass_axis.linearisation)
        searching = NearestNeighborStrategy(axis, axis_range, None)

        # Without this the test could pass by comparing two searches.
        assert computing.linearisation is not None
        assert searching.linearisation is None

        mzs, intensities = _spectrum(8, n=3_000, lo=90.0, hi=2_010.0)
        a_bins, a_sums = searching.resample(mzs, intensities)
        b_bins, b_sums = computing.resample(mzs, intensities)

        np.testing.assert_array_equal(b_bins, a_bins)
        np.testing.assert_array_equal(b_sums, a_sums)
        assert computing.out_of_range_peaks == searching.out_of_range_peaks

    def test_map_to_bins_is_the_public_spelling(self):
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 2000.0, 10_000)
        axis = np.asarray(mass_axis.mz_values, float)
        strategy = NearestNeighborStrategy(axis, None, mass_axis.linearisation)
        mzs = np.sort(np.random.default_rng(9).uniform(100.0, 2000.0, 500))

        np.testing.assert_array_equal(
            strategy.map_to_bins(mzs),
            nn_map_to_bins(axis, mzs, mass_axis.linearisation),
        )


class TestTheDeclaredBound:
    """Issue #239, for the interpolating method: a peak on the bound survives."""

    def test_a_peak_on_a_bound_outside_the_axis_span_is_kept(self):
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 1000.0, 5_000)
        axis = np.asarray(mass_axis.mz_values, float)
        # The generator reports bin centres, so the declared range is
        # wider than the axis's own span at both ends.
        assert axis[0] > 100.0 and axis[-1] < 1000.0

        declared = TICPreservingStrategy(axis, (100.0, 1000.0))
        assert declared.kept_range == (100.0, 1000.0)

        indices, values = declared.resample(
            np.array([100.0, 500.0, 1000.0]), np.array([7.0, 5.0, 11.0])
        )

        assert values.sum() == pytest.approx(23.0)
        assert 0 in indices.tolist()
        assert axis.size - 1 in indices.tolist()

    def test_without_the_declared_range_they_fall_outside(self):
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 1000.0, 5_000)
        axis = np.asarray(mass_axis.mz_values, float)

        _, values = TICPreservingStrategy(axis, None).resample(
            np.array([100.0, 500.0, 1000.0]), np.array([7.0, 5.0, 11.0])
        )

        assert values.sum() < 23.0


class TestOutOfRangeAccounting:
    def test_the_count_accumulates_across_spectra(self):
        strategy = NearestNeighborStrategy(np.linspace(100.0, 110.0, 11), None)

        strategy.resample(np.array([90.0, 105.0, 120.0]), np.ones(3))
        strategy.resample(np.array([95.0, 105.0]), np.ones(2))

        assert strategy.out_of_range_peaks == 3

    def test_it_warns_once_and_names_the_kept_range(self, warnings_from_the_strategy):
        strategy = NearestNeighborStrategy(np.linspace(100.0, 110.0, 11), None)

        strategy.resample(np.array([90.0, 105.0]), np.ones(2))
        strategy.resample(np.array([90.0, 105.0]), np.ones(2))

        warnings = [
            record
            for record in warnings_from_the_strategy.records
            if record.name == LOGGER
        ]
        assert len(warnings) == 1, "one line per conversion, not per spectrum"
        assert "100.0000" in warnings[0].getMessage()
        assert "110.0000" in warnings[0].getMessage()

    def test_nothing_is_said_when_nothing_is_dropped(self, warnings_from_the_strategy):
        strategy = NearestNeighborStrategy(np.linspace(100.0, 110.0, 11), None)

        strategy.resample(np.array([101.0, 105.0]), np.ones(2))

        assert strategy.out_of_range_peaks == 0
        assert not [
            record
            for record in warnings_from_the_strategy.records
            if record.name == LOGGER
        ]


class TestTheFactory:
    """``build_strategy`` is the only way the converter makes one."""

    def test_nearest_neighbor_gets_the_linearisation(self):
        mass_axis = ReflectorTOFAxisGenerator().generate_axis(100.0, 1000.0, 1_000)
        axis = np.asarray(mass_axis.mz_values, float)

        strategy = build_strategy(
            ResamplingMethod.NEAREST_NEIGHBOR,
            axis,
            (100.0, 1000.0),
            mass_axis.linearisation,
            0.05,
        )

        assert isinstance(strategy, NearestNeighborStrategy)
        assert strategy.linearisation is mass_axis.linearisation
        assert strategy.axis_range == (100.0, 1000.0)

    def test_tic_preserving_gets_the_gap_tolerance(self):
        strategy = build_strategy(
            ResamplingMethod.TIC_PRESERVING, _axis(), (100.0, 1000.0), None, 0.05
        )

        assert isinstance(strategy, TICPreservingStrategy)
        assert strategy.gap_tolerance_da == 0.05

    def test_none_has_no_operator(self):
        with pytest.raises(ValueError, match="ResamplingMethod.NONE"):
            build_strategy(ResamplingMethod.NONE, _axis(), None, None, None)
