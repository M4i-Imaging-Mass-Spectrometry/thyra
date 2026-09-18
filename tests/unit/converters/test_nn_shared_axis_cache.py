# tests/unit/converters/test_nn_shared_axis_cache.py
"""The shared-axis nearest-neighbor cache must be an invisible optimization.

``NearestNeighborStrategy`` caches the peak-to-bin mapping the first time
it sees a spectrum and reuses it for every later spectrum carrying the same
m/z array -- which is every spectrum, on a shared-axis reader (continuous
imzML, Rapiflex, Waters, PHI). The cache must never change a result:

* a hit must return exactly what the generic path returns;
* an array with different values must miss and take the generic path,
  never receive the stale mapping;
* repeated misses (processed-mode data) must disable the cache;
* out-of-range peaks must be counted identically through both routes.

The bit-parity claim is load-bearing: the streaming PCS route resamples
every spectrum twice and requires pass 2 to reproduce pass 1 exactly, and
mixed hit/miss sequences within one conversion would otherwise let the two
passes disagree.
"""

from types import MethodType, SimpleNamespace

import numpy as np

from thyra.resampling.strategies import NearestNeighborStrategy


def _strategy(axis, cache_enabled=True):
    """The nearest-neighbour operator over one axis, cache on or off."""
    strategy = NearestNeighborStrategy(np.asarray(axis, float), None)
    strategy._out_of_range_warned = True  # count, but stay quiet
    if not cache_enabled:
        strategy._shared_cache = False
    return strategy


def _spectra(rng, mzs, n):
    """n intensity vectors over one shared m/z array, zeros included."""
    out = []
    for _ in range(n):
        ints = rng.lognormal(np.log(50), 1.0, mzs.size)
        ints[rng.random(mzs.size) < 0.4] = 0.0
        out.append(ints)
    return out


class TestHitsMatchTheGenericPath:
    def test_shared_array_bit_parity(self):
        rng = np.random.default_rng(3)
        axis = np.linspace(100.0, 1000.0, 5_000)
        mzs = np.sort(rng.uniform(100.0, 1000.0, 2_000))

        cached = _strategy(axis)
        generic = _strategy(axis, cache_enabled=False)

        for ints in _spectra(rng, mzs, 6):
            got = cached.resample(mzs, ints)
            want = generic.resample(mzs, ints)
            np.testing.assert_array_equal(got[0], want[0])
            np.testing.assert_array_equal(got[1], want[1])  # exact, not approx
            assert got[0].dtype == want[0].dtype
            assert got[1].dtype == want[1].dtype

        assert cached._shared_cache not in (None, False), "cache never built"

    def test_equal_valued_copy_hits(self):
        rng = np.random.default_rng(4)
        axis = np.linspace(100.0, 1000.0, 3_000)
        mzs = np.sort(rng.uniform(100.0, 1000.0, 500))
        ints = rng.exponential(10.0, mzs.size)

        stub = _strategy(axis)
        first = stub.resample(mzs, ints)
        # A fresh, equal-valued array (what continuous imzML yields when the
        # block is re-decoded) must hit, not rebuild or miss.
        again = stub.resample(mzs.copy(), ints)
        np.testing.assert_array_equal(first[0], again[0])
        np.testing.assert_array_equal(first[1], again[1])
        assert stub._cache_misses == 0

    def test_out_of_range_counting_is_identical(self):
        axis = np.linspace(400.0, 800.0, 1_000)
        # Two below, one inside, one above.
        mzs = np.array([300.0, 350.0, 500.0, 900.0])
        ints = np.array([1.0, 2.0, 3.0, 4.0])

        cached = _strategy(axis)
        generic = _strategy(axis, cache_enabled=False)
        for _ in range(3):
            got = cached.resample(mzs, ints)
            want = generic.resample(mzs, ints)
            np.testing.assert_array_equal(got[0], want[0])
            np.testing.assert_array_equal(got[1], want[1])

        assert cached.out_of_range_peaks == generic.out_of_range_peaks == 9


class TestMissesStayCorrect:
    def test_different_values_take_the_generic_path(self):
        rng = np.random.default_rng(5)
        axis = np.linspace(100.0, 1000.0, 3_000)
        first = np.sort(rng.uniform(100.0, 1000.0, 400))
        other = np.sort(rng.uniform(100.0, 1000.0, 400))  # same size, new values
        ints = np.ones(400)

        stub = _strategy(axis)
        stub.resample(first, ints)  # builds the cache
        got = stub.resample(other, ints)
        want = _strategy(axis, cache_enabled=False).resample(other, ints)
        np.testing.assert_array_equal(got[0], want[0])
        np.testing.assert_array_equal(got[1], want[1])
        assert stub._cache_misses == 1

    def test_processed_mode_disables_the_cache(self):
        rng = np.random.default_rng(6)
        axis = np.linspace(100.0, 1000.0, 3_000)
        stub = _strategy(axis)

        for _ in range(8):
            n = int(rng.integers(50, 200))
            mzs = np.sort(rng.uniform(100.0, 1000.0, n))
            stub.resample(mzs, np.ones(n))

        assert stub._shared_cache is False, "five misses must disable it"

    def test_unsorted_mzs_are_not_cached(self):
        axis = np.linspace(100.0, 1000.0, 1_000)
        mzs = np.array([500.0, 300.0, 700.0])  # descending start: not ascending
        ints = np.array([1.0, 2.0, 3.0])

        stub = _strategy(axis)
        got = stub.resample(mzs, ints)
        want = _strategy(axis, cache_enabled=False).resample(mzs, ints)
        np.testing.assert_array_equal(got[0], want[0])
        np.testing.assert_array_equal(got[1], want[1])
        assert stub._shared_cache is False


class TestIdentityMassIndices:
    """_map_mass_to_indices short-circuits mzs == axis to an arange."""

    def _converter_stub(self, axis):
        stub = SimpleNamespace(
            _common_mass_axis=np.asarray(axis),
            _identity_mass_indices=None,
            _axis_strictly_increasing=None,
        )
        from thyra.core.base_converter import BaseMSIConverter

        stub._map_mass_to_indices = MethodType(
            BaseMSIConverter._map_mass_to_indices, stub
        )
        return stub

    def test_identity_on_equal_array(self):
        axis = np.linspace(100.0, 1000.0, 4_000)
        stub = self._converter_stub(axis)
        got = stub._map_mass_to_indices(axis.copy())
        np.testing.assert_array_equal(got, np.arange(axis.size))

    def test_duplicate_axis_keeps_searchsorted_semantics(self):
        # searchsorted maps a duplicated value to its first occurrence, so
        # the arange shortcut is NOT valid here and must not fire.
        axis = np.array([100.0, 200.0, 200.0, 300.0])
        stub = self._converter_stub(axis)
        got = stub._map_mass_to_indices(axis.copy())
        expected = np.searchsorted(axis, axis)  # [0, 1, 1, 3]
        np.testing.assert_array_equal(got, expected)

    def test_subset_still_maps_exactly(self):
        axis = np.linspace(100.0, 1000.0, 4_000)
        stub = self._converter_stub(axis)
        subset = axis[::7]
        got = stub._map_mass_to_indices(subset)
        np.testing.assert_array_equal(got, np.arange(0, axis.size, 7))
