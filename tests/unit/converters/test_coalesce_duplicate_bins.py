"""What ``_coalesce_duplicate_bins`` does, pinned before it is made faster.

It had no test at all -- ``grep -rn coalesce tests/`` was empty -- while
running once per spectrum on the whole mass axis of every ``--no-resample``
conversion (issue #311). These are characterization tests: every one passes
on the tree as it stood, and their job is to fail if the faster ascending
check changes an answer.

Three behaviours, and the third is the one easy to lose. The function is
named for summing duplicates, but its fallback is ``np.unique``, so it also
**sorts** an unsorted input. Anything that skips the slow path has to leave
that intact.
"""

from __future__ import annotations

import numpy as np
import pytest

from thyra.converters.spatialdata.base_spatialdata_converter import (
    BaseSpatialDataConverter,
)

_coalesce = BaseSpatialDataConverter._coalesce_duplicate_bins


class TestItSumsDuplicates:
    def test_repeated_bins_are_summed(self):
        bins, values = _coalesce(
            np.asarray([0, 2, 2, 5, 5, 5, 9]),
            np.asarray([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
        )
        np.testing.assert_array_equal(bins, [0, 2, 5, 9])
        np.testing.assert_array_equal(values, [1.0, 5.0, 15.0, 7.0])

    def test_the_total_is_conserved(self):
        """The reason it exists: the old COO route summed, so a store
        written either way reads back the same."""
        values = np.asarray([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        _, summed = _coalesce(np.asarray([0, 2, 2, 5, 5, 5, 9]), values)
        assert summed.sum() == pytest.approx(values.sum())

    def test_every_bin_the_same(self):
        bins, values = _coalesce(np.asarray([3, 3, 3]), np.asarray([1.0, 2.0, 3.0]))
        np.testing.assert_array_equal(bins, [3])
        np.testing.assert_array_equal(values, [6.0])


class TestItAlsoSorts:
    """Not in the name, and easy to lose: the fallback is ``np.unique``."""

    def test_unsorted_input_comes_back_sorted(self):
        bins, values = _coalesce(np.asarray([5, 1, 3]), np.asarray([1.0, 2.0, 3.0]))
        np.testing.assert_array_equal(bins, [1, 3, 5])
        np.testing.assert_array_equal(values, [2.0, 3.0, 1.0])

    def test_unsorted_with_duplicates(self):
        bins, values = _coalesce(
            np.asarray([5, 1, 5, 1]), np.asarray([1.0, 2.0, 3.0, 4.0])
        )
        np.testing.assert_array_equal(bins, [1, 5])
        np.testing.assert_array_equal(values, [6.0, 4.0])


class TestTheAscendingFastPath:
    """The identity case, which is every ordinary spectrum."""

    def test_a_strictly_ascending_run_is_returned_unchanged(self):
        indices = np.asarray([0, 1, 2, 3])
        intensities = np.asarray([1.0, 2.0, 3.0, 4.0])
        bins, values = _coalesce(indices, intensities)
        # The same objects, not equal copies: the fast path must not
        # allocate, which is the whole point of checking first.
        assert bins is indices
        assert values is intensities

    def test_the_memoised_arange_is_returned_unchanged(self):
        """What ``_map_mass_to_indices`` hands over on a shared axis."""
        indices = np.arange(331_000)
        intensities = np.ones(331_000)
        bins, values = _coalesce(indices, intensities)
        assert bins is indices
        assert values is intensities

    def test_a_gap_is_still_ascending(self):
        indices = np.asarray([0, 7, 8, 400])
        bins, _ = _coalesce(indices, np.ones(4))
        assert bins is indices

    @pytest.mark.parametrize("size", [0, 1])
    def test_too_short_to_have_a_duplicate(self, size):
        indices = np.arange(size)
        intensities = np.ones(size)
        bins, values = _coalesce(indices, intensities)
        assert bins is indices
        assert values is intensities

    def test_the_dtype_survives(self):
        for dtype in (np.int32, np.int64):
            bins, _ = _coalesce(
                np.asarray([0, 2, 2], dtype=dtype), np.asarray([1.0, 2.0, 3.0])
            )
            assert bins.dtype == dtype
