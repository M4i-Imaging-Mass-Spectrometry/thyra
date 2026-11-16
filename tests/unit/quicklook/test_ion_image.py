"""From an m/z value to a raster: the window of bins, then the pixels."""

from __future__ import annotations

import re

import numpy as np
import pytest

from tests.fixtures.hand_built_store import (
    L_SHAPE,
    MZ,
    dense_matrix,
    make_table,
    write_store,
)
from thyra.quicklook.ion_image import (
    make_ion_image,
    place_pixels,
    select_bins,
    tallest_peak_near,
)
from thyra.quicklook.store import Positions, QuicklookError, open_store

AXIS = np.array([999.0, 1000.0, 1001.0, 1002.0])


def _positions(points) -> Positions:
    points = np.asarray(points, dtype=np.int64)
    z = points[:, 2] if points.shape[1] == 3 else None
    return Positions(points[:, 0], points[:, 1], z)


class TestSelectBins:
    def test_a_bin_on_either_edge_is_inside(self):
        """1000 ppm of 1000 is exactly 1: both 999 and 1001 belong."""
        window = select_bins(AXIS, 1000.0, 1000.0)

        assert (window.first, window.stop) == (0, 3)
        assert (window.low, window.high, window.n_bins) == (999.0, 1001.0, 3)

    def test_a_bin_just_past_the_edge_is_outside(self):
        window = select_bins(AXIS, 1000.0, 999.0)

        assert (window.first, window.stop) == (1, 2)

    def test_the_window_scales_with_the_mz_value(self):
        """The same 10 ppm is ten times as wide at ten times the m/z."""
        axis = np.array([99.998, 100.0, 100.002, 999.998, 1000.0, 1000.002])

        assert select_bins(axis, 100.0, 10.0).n_bins == 1
        assert select_bins(axis, 1000.0, 10.0).n_bins == 3

    @pytest.mark.parametrize("target", [500.0, 2000.0])
    def test_an_mz_outside_the_axis_is_refused(self, target):
        with pytest.raises(QuicklookError) as refusal:
            select_bins(AXIS, target, 10.0)

        message = str(refusal.value)
        assert "outside the mass axis" in message
        assert "999.0000 to 1002.0000" in message

    def test_a_window_that_reaches_the_first_bin_from_outside_is_kept(self):
        """Only a window that touches no part of the axis is outside it."""
        window = select_bins(AXIS, 998.5, 1000.0)

        assert (window.first, window.stop) == (0, 1)

    def test_an_empty_window_names_the_nearest_bin_and_a_ppm_that_works(self):
        axis = np.array([100.0, 200.0])

        with pytest.raises(QuicklookError) as refusal:
            select_bins(axis, 160.0, 10.0)

        message = str(refusal.value)
        assert "No m/z bin lies within 10 ppm of 160.0000" in message
        assert "nearest bin is at 200.0000" in message
        suggested = int(re.findall(r"--ppm (\d+)", message)[0])
        assert select_bins(axis, 160.0, suggested).n_bins == 1
        with pytest.raises(QuicklookError):
            select_bins(axis, 160.0, suggested - 1)

    def test_an_axis_with_no_bins_is_refused(self):
        with pytest.raises(QuicklookError, match="no m/z bins"):
            select_bins(np.array([]), 100.0, 10.0)


class TestPlacePixels:
    def test_each_value_lands_on_its_row_and_column(self):
        points = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]

        image = place_pixels(_positions(points), np.arange(6.0))

        np.testing.assert_array_equal(image, [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])

    def test_a_position_never_measured_is_nan_not_zero(self):
        points = [(0, 0), (2, 0), (0, 1), (1, 1)]

        image = place_pixels(_positions(points), np.array([0.0, 5.0, 6.0, 7.0]))

        assert image[0, 0] == 0.0  # measured, no signal
        assert np.isnan(image[0, 1]) and np.isnan(image[1, 2])
        assert np.count_nonzero(np.isnan(image)) == 2

    def test_a_region_that_is_not_a_rectangle_keeps_its_shape(self):
        values = np.arange(1.0, len(L_SHAPE) + 1)

        image = place_pixels(_positions(L_SHAPE), values)

        expected = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, np.nan]])
        np.testing.assert_array_equal(image, expected)

    def test_the_raster_is_the_bounding_box_not_the_distance_from_zero(self):
        """The L sits at columns 10 to 12: the raster is 3 wide, not 13."""
        image = place_pixels(_positions(L_SHAPE), np.ones(len(L_SHAPE)))

        assert image.shape == (3, 3)

    def test_a_volume_is_drawn_as_its_slices_side_by_side(self):
        points = [(0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1), (0, 0, 2)]

        image = place_pixels(_positions(points), np.array([1.0, 2.0, 3.0, 4.0, 5.0]))

        # Two slices to a line, one NaN column and row between them.
        nan = np.nan
        expected = [
            [1.0, 2.0, nan, 3.0, 4.0],
            [nan, nan, nan, nan, nan],
            [5.0, nan, nan, nan, nan],
        ]
        np.testing.assert_array_equal(image, expected)

    def test_no_rows_is_refused(self):
        empty = Positions(
            np.array([], dtype=np.int64), np.array([], dtype=np.int64), None
        )

        with pytest.raises(QuicklookError, match="no pixels"):
            place_pixels(empty, np.array([]))


class TestMakeIonImage:
    def test_it_adds_up_exactly_the_bins_inside_the_window(self, store_path):
        """Compared with the matrix the store was built from, not read back."""
        table = open_store(store_path).open_table("msi_dataset_z0")

        ion = make_ion_image(table, 505.0, ppm=2500.0)  # 503.7375 to 506.2625

        expected = dense_matrix(len(L_SHAPE))[:, 4:7].sum(axis=1)
        assert (ion.window.first, ion.window.stop) == (4, 7)
        np.testing.assert_array_equal(ion.image[~np.isnan(ion.image)], expected)
        assert (ion.lowest, ion.highest) == (expected.min(), expected.max())
        assert ion.n_pixels == len(L_SHAPE)

    @pytest.mark.parametrize("column", [0, MZ.size - 1])
    def test_the_first_and_the_last_bin_can_be_read(self, store_path, column):
        table = open_store(store_path).open_table("msi_dataset_z0")

        ion = make_ion_image(table, float(MZ[column]), ppm=10.0)

        expected = dense_matrix(len(L_SHAPE))[:, column]
        np.testing.assert_array_equal(ion.image[~np.isnan(ion.image)], expected)

    def test_a_window_with_bins_but_no_signal_is_an_image_of_zeros(self, store_path):
        table = open_store(store_path).open_table("msi_dataset_z0")

        ion = make_ion_image(table, 503.0, ppm=10.0)

        assert ion.highest == 0.0
        assert np.count_nonzero(ion.image == 0.0) == len(L_SHAPE)

    def test_a_matrix_not_stored_by_column_is_refused(self, tmp_path):
        path = write_store(tmp_path / "old.zarr", {"t": make_table(csr=True)})
        table = open_store(path).open_table("t")

        with pytest.raises(QuicklookError, match="column by column"):
            make_ion_image(table, 505.0)


class TestTallestPeakNear:
    def test_it_finds_the_tallest_bin_of_the_mean_spectrum_nearby(self, store_path):
        table = open_store(store_path).open_table("msi_dataset_z0")

        # Bins 502, 503 and 504 are within reach; 503 is empty, 504 tallest.
        assert tallest_peak_near(table, 503.0) == 504.0

    def test_nothing_nearby_is_none(self, store_path):
        table = open_store(store_path).open_table("msi_dataset_z0")

        assert tallest_peak_near(table, 503.0, reach=0.5) is None
