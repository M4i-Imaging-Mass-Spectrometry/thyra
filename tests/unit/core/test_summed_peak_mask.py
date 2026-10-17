"""The intensity threshold's one unit: a peak summed over its listed points (D34)."""

import numpy as np
import pytest

from thyra.core.base_reader import SAME_PEAK_RELATIVE, summed_peak_mask


def test_distinct_mz_are_tested_point_by_point():
    mzs = np.array([100.0, 200.0, 300.0])
    intensities = np.array([5.0, 50.0, 500.0])

    keep = summed_peak_mask(mzs, intensities, 50.0)

    np.testing.assert_array_equal(keep, [False, True, True])


def test_a_repeated_mz_is_kept_or_dropped_as_one_peak():
    # Two scans' points of each peak, listed scan by scan as an archive does.
    mzs = np.array([100.0, 200.0, 100.0, 200.0])
    intensities = np.array([30.0, 10.0, 30.0, 10.0])

    keep = summed_peak_mask(mzs, intensities, 50.0)

    # 100 sums to 60 and stays whole; 200 sums to 20 and goes whole.
    np.testing.assert_array_equal(keep, [True, False, True, False])


def test_copies_that_differ_in_the_last_bit_are_one_peak():
    """A chunk bound and the decoded value of the same bin."""
    mz = 512.123456789
    twin = np.nextafter(mz, np.inf)
    assert twin - mz < SAME_PEAK_RELATIVE * mz
    mzs = np.array([mz, twin, 512.124])
    intensities = np.array([40.0, 40.0, 40.0])

    keep = summed_peak_mask(mzs, intensities, 50.0)

    np.testing.assert_array_equal(keep, [True, True, False])


def test_empty_and_single_point_spectra():
    assert summed_peak_mask(np.array([]), np.array([]), 1.0).size == 0
    np.testing.assert_array_equal(
        summed_peak_mask(np.array([10.0]), np.array([0.5]), 1.0), [False]
    )


def test_a_nan_mz_is_in_no_peak():
    """Never folded into the largest real peak, where argsort puts it.

    The point itself is left for the converter, which refuses or drops a
    non-finite m/z whether or not a threshold is set (D55).
    """
    mzs = np.array([100.0, 200.0, np.nan])
    intensities = np.array([1.0, 30.0, 30.0])

    keep = summed_peak_mask(mzs, intensities, 50.0)

    np.testing.assert_array_equal(keep, [False, False, True])


class TestOnlyMeasurementsCount:
    """D55: the sum holds what the converter would store, nothing else.

    The converter drops a non-finite or negative intensity before
    resampling (``_drop_unusable_intensities``), so a sum that held one
    tested a peak the store never holds.
    """

    def test_a_nan_copy_does_not_remove_the_peak(self):
        mzs = np.array([400.0, 500.0, 500.0])
        intensities = np.array([100.0, np.nan, 100.0])

        keep = summed_peak_mask(mzs, intensities, 50.0)

        np.testing.assert_array_equal(keep, [True, True, True])

    def test_a_negative_copy_is_not_subtracted(self):
        mzs = np.array([400.0, 500.0, 500.0])
        intensities = np.array([100.0, -30.0, 60.0])

        keep = summed_peak_mask(mzs, intensities, 50.0)

        np.testing.assert_array_equal(keep, [True, True, True])

    def test_an_infinite_copy_does_not_lift_the_peak(self):
        mzs = np.array([400.0, 500.0, 500.0])
        intensities = np.array([100.0, np.inf, 10.0])

        keep = summed_peak_mask(mzs, intensities, 50.0)

        # The 10 goes; the infinity is left for the converter to drop.
        np.testing.assert_array_equal(keep, [True, True, False])

    @pytest.mark.parametrize("bad_mz", [np.inf, -np.inf])
    def test_an_infinite_mz_joins_no_finite_peak(self, bad_mz):
        mzs = np.array([400.0, 500.0, bad_mz])
        intensities = np.array([100.0, 10.0, 100.0])

        keep = summed_peak_mask(mzs, intensities, 50.0)

        np.testing.assert_array_equal(keep, [True, False, True])

    def test_two_infinite_mz_are_not_one_peak(self):
        mzs = np.array([np.inf, np.inf, 500.0])
        intensities = np.array([30.0, 30.0, 100.0])

        keep = summed_peak_mask(mzs, intensities, 50.0)

        np.testing.assert_array_equal(keep, [True, True, True])

    def test_a_lone_unusable_point_is_left_in(self):
        """The single-point path follows the same rule."""
        for intensity in (np.nan, np.inf, -5.0):
            np.testing.assert_array_equal(
                summed_peak_mask(np.array([10.0]), np.array([intensity]), 1.0),
                [True],
            )
