"""The intensity threshold's one unit: a peak summed over its listed points (D34)."""

import numpy as np

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


def test_a_nan_mz_is_its_own_peak():
    """Never folded into the largest real peak, where argsort puts it."""
    mzs = np.array([100.0, 200.0, np.nan])
    intensities = np.array([1.0, 30.0, 30.0])

    keep = summed_peak_mask(mzs, intensities, 50.0)

    np.testing.assert_array_equal(keep, [False, False, False])
