"""Tests for Waters imaging grid reconstruction.

These tests cover ImagingGrid, _mm_to_um_key, and build_imaging_grid
without requiring native DLLs -- MassLynxLib is mocked.
"""

import logging
from unittest.mock import MagicMock

import pytest

from thyra.errors import ConversionRefused
from thyra.readers.waters.imaging_grid import (
    ImagingGrid,
    _grid_from_scan_map,
    _mm_to_um_key,
    build_imaging_grid,
    regrid_for_functions,
)
from thyra.readers.waters.instrument import DeclaredRegion, StagePattern
from thyra.readers.waters.masslynx_lib import FunctionType, ScanInfoData


def _make_scan_info(x_mm, y_mm, has_pos=True):
    """Create a ScanInfoData with given laser position."""
    return ScanInfoData(
        ms_level=1,
        polarity=0,
        drift_scan_count=0,
        is_profile=0,
        precursor_mz=0.0,
        rt=1.0,
        laser_x_pos=x_mm if has_pos else -1.0,
        laser_y_pos=y_mm if has_pos else -1.0,
    )


class TestMmToUmKey:
    """Test the mm-to-um conversion helper."""

    def test_basic_conversion(self):
        """Test that mm values are converted to um with rounding."""
        assert _mm_to_um_key(1.0) == 1000.0
        assert _mm_to_um_key(0.1) == 100.0
        assert _mm_to_um_key(0.001) == 1.0

    def test_rounding(self):
        """Test that values are rounded to 2 decimal places."""
        # 0.1234 mm = 123.4 um -> rounded to 123.4
        assert _mm_to_um_key(0.1234) == 123.4

    def test_zero(self):
        """Test zero value."""
        assert _mm_to_um_key(0.0) == 0.0

    def test_negative(self):
        """Test negative values (sentinel positions)."""
        assert _mm_to_um_key(-1.0) == -1000.0


class TestScanInfoData:
    """Test the ScanInfoData dataclass."""

    def test_has_position_true(self):
        """Test has_position for valid coordinates."""
        info = _make_scan_info(0.5, 0.3)
        assert info.has_position is True

    def test_has_position_false(self):
        """Test has_position for sentinel -1.0 coordinates."""
        info = _make_scan_info(0.0, 0.0, has_pos=False)
        assert info.has_position is False

    def test_has_position_partial_sentinel(self):
        """Test has_position when only one coord is sentinel."""
        # Both must be -1.0 for has_position to be False
        info = ScanInfoData(
            ms_level=1,
            polarity=0,
            drift_scan_count=0,
            is_profile=0,
            precursor_mz=0.0,
            rt=1.0,
            laser_x_pos=-1.0,
            laser_y_pos=0.5,
        )
        assert info.has_position is True


class TestImagingGrid:
    """Test ImagingGrid coordinate mapping."""

    def _make_grid(self):
        """Create a 3x2 grid for testing."""
        x_map = {100.0: 0, 200.0: 1, 300.0: 2}
        y_map = {50.0: 0, 150.0: 1}
        return ImagingGrid(
            x_index_map=x_map,
            y_index_map=y_map,
            pixel_count_x=3,
            pixel_count_y=2,
            pixel_size_x=100.0,
            pixel_size_y=100.0,
            lateral_width=200.0,
            lateral_height=100.0,
            scan_map={},
        )

    def test_dimensions(self):
        """Test that dimensions returns (x, y, 1)."""
        grid = self._make_grid()
        assert grid.dimensions == (3, 2, 1)

    def test_get_coordinates_valid(self):
        """Test coordinate mapping for a known position."""
        grid = self._make_grid()
        # 0.1 mm = 100.0 um, 0.05 mm = 50.0 um
        info = _make_scan_info(0.1, 0.05)
        coords = grid.get_coordinates(info)
        assert coords == (0, 0, 0)

    def test_get_coordinates_middle(self):
        """Test coordinate mapping for a middle position."""
        grid = self._make_grid()
        # 0.2 mm = 200.0 um, 0.15 mm = 150.0 um
        info = _make_scan_info(0.2, 0.15)
        coords = grid.get_coordinates(info)
        assert coords == (1, 1, 0)

    def test_get_coordinates_last(self):
        """Test coordinate mapping for the last position."""
        grid = self._make_grid()
        # 0.3 mm = 300.0 um, 0.15 mm = 150.0 um
        info = _make_scan_info(0.3, 0.15)
        coords = grid.get_coordinates(info)
        assert coords == (2, 1, 0)

    def test_get_coordinates_no_position(self):
        """Test that scans without position return None."""
        grid = self._make_grid()
        info = _make_scan_info(0.0, 0.0, has_pos=False)
        assert grid.get_coordinates(info) is None

    def test_get_coordinates_unknown_position(self):
        """Test that unknown laser positions return None."""
        grid = self._make_grid()
        # 0.999 mm = 999.0 um -- not in the grid
        info = _make_scan_info(0.999, 0.999)
        assert grid.get_coordinates(info) is None


class TestBuildImagingGrid:
    """Test the build_imaging_grid function with mocked MassLynxLib."""

    def _setup_mock_ml(self, positions):
        """Create a mock MassLynxLib that returns given laser positions.

        Args:
            positions: List of (func, scan, x_mm, y_mm) tuples.
                       Use x_mm=-1.0, y_mm=-1.0 for scans without position.
        """
        mock_ml = MagicMock()

        # Group by function to determine function count and scan counts
        func_scans = {}
        for func, scan, x_mm, y_mm in positions:
            func_scans.setdefault(func, []).append((scan, x_mm, y_mm))

        n_functions = max(f for f, _, _, _ in positions) + 1
        mock_ml.get_number_of_functions.return_value = n_functions

        def get_n_scans(handle, func):
            return len(func_scans.get(func, []))

        mock_ml.get_number_of_scans_in_function.side_effect = get_n_scans

        def get_scan_info(handle, func, scan):
            for f, s, x, y in positions:
                if f == func and s == scan:
                    return _make_scan_info(x, y, has_pos=(x != -1.0))
            return _make_scan_info(0, 0, has_pos=False)

        mock_ml.get_scan_info.side_effect = get_scan_info

        return mock_ml

    def test_basic_3x2_grid(self):
        """Test building a 3x2 grid from 6 positioned scans."""
        positions = [
            (0, 0, 0.1, 0.05),
            (0, 1, 0.2, 0.05),
            (0, 2, 0.3, 0.05),
            (0, 3, 0.1, 0.15),
            (0, 4, 0.2, 0.15),
            (0, 5, 0.3, 0.15),
        ]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS}

        grid = build_imaging_grid(mock_ml, "handle", func_types)

        assert grid.pixel_count_x == 3
        assert grid.pixel_count_y == 2
        assert grid.lateral_width == 200.0  # 300 - 100 um
        assert grid.lateral_height == 100.0  # 150 - 50 um
        assert len(grid.scan_map) == 6

    def test_skips_scans_without_position(self):
        """Test that scans with sentinel positions are excluded from the grid."""
        positions = [
            (0, 0, 0.1, 0.05),
            (0, 1, -1.0, -1.0),  # no position
            (0, 2, 0.2, 0.05),
        ]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS}

        grid = build_imaging_grid(mock_ml, "handle", func_types)

        assert grid.pixel_count_x == 2  # only 0.1 and 0.2
        assert grid.pixel_count_y == 1
        assert len(grid.scan_map) == 3  # all scans in map (including no-pos)

    def test_single_pixel_grid_is_refused(self):
        """A stage that never moved is a spot acquisition, not an image.

        This used to return a 1x1 grid with a 0.0 um pitch and let the
        conversion report success with the whole acquisition on one pixel
        (issue #213). On the Xevo DESI share that is what 76 of 105 .raw
        dirs look like -- calibrations, tuning runs, single-spot DDA.
        """
        positions = [(0, 0, 0.5, 0.5), (0, 1, 0.5, 0.5), (0, 2, 0.5, 0.5)]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS}

        with pytest.raises(ConversionRefused, match="same stage position"):
            build_imaging_grid(mock_ml, "handle", func_types)

    def test_single_pixel_grid_is_returned_for_metadata(self, thyra_logs):
        """A metadata document describes a spot acquisition too (issue #494).

        The grid carries the refusal instead of raising it, with a 0.0 um
        pitch, which the metadata extractor reports as no pixel size. Every
        scan shares the one pixel, so no shared-pixel warning is logged.
        """
        positions = [(0, 0, 0.5, 0.5), (0, 1, 0.5, 0.5), (1, 0, 0.5, 0.5)]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS, 1: FunctionType.MS}

        with thyra_logs("thyra.readers.waters.imaging_grid") as records:
            grid = build_imaging_grid(
                mock_ml, "handle", func_types, allow_single_position=True
            )

        assert (grid.pixel_count_x, grid.pixel_count_y) == (1, 1)
        assert (grid.pixel_size_x, grid.pixel_size_y) == (0.0, 0.0)
        assert grid.no_raster_reason is not None
        assert "same stage position" in grid.no_raster_reason
        assert not [r for r in records if "SUM" in r.getMessage()]
        refitted = regrid_for_functions(grid, [0], allow_single_position=True)
        assert refitted.no_raster_reason is not None
        with pytest.raises(ConversionRefused, match="same stage position"):
            regrid_for_functions(grid, [0])

    def test_a_raster_carries_no_refusal(self):
        positions = [(0, 0, 0.1, 0.05), (0, 1, 0.2, 0.05)]
        mock_ml = self._setup_mock_ml(positions)

        grid = build_imaging_grid(
            mock_ml, "handle", {0: FunctionType.MS}, allow_single_position=True
        )

        assert grid.no_raster_reason is None

    def test_single_row_grid_is_kept(self):
        """One line of a raster is real data, and must survive the refusal.

        The share holds aborted line scans (138x1, 48x1) alongside the
        completed rasters they were retried into. Only a grid that
        collapses on *both* axes is a spot.
        """
        positions = [
            (0, 0, 0.1, 0.05),
            (0, 1, 0.2, 0.05),
            (0, 2, 0.3, 0.05),
        ]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS}

        grid = build_imaging_grid(mock_ml, "handle", func_types)

        assert grid.pixel_count_x == 3
        assert grid.pixel_count_y == 1

    def test_no_valid_positions_raises(self):
        """Test error when no scans have valid positions."""
        positions = [
            (0, 0, -1.0, -1.0),
            (0, 1, -1.0, -1.0),
        ]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS}

        with pytest.raises(ConversionRefused, match="No valid laser positions"):
            build_imaging_grid(mock_ml, "handle", func_types)

    def test_multiple_functions(self):
        """Test grid built from multiple acquisition functions."""
        positions = [
            (0, 0, 0.1, 0.05),
            (0, 1, 0.2, 0.05),
            (1, 0, 0.1, 0.15),  # same x, different y
            (1, 1, 0.2, 0.15),
        ]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS, 1: FunctionType.MS}

        grid = build_imaging_grid(mock_ml, "handle", func_types)

        assert grid.pixel_count_x == 2
        assert grid.pixel_count_y == 2
        assert len(grid.scan_map) == 4

    def test_pixel_size_calculation(self):
        """The pitch is the extent over the *intervals*, not the positions.

        Four positions 100 um apart are a 100 um raster. Dividing the
        300 um centre-to-centre span by the position count reported 75 um
        (issue #217); it is 300 / 3.
        """
        positions = [
            (0, 0, 0.1, 0.1),
            (0, 1, 0.2, 0.1),
            (0, 2, 0.3, 0.1),
            (0, 3, 0.4, 0.1),
        ]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS}

        grid = build_imaging_grid(mock_ml, "handle", func_types)

        assert grid.lateral_width == 300.0
        assert grid.pixel_size_x == 100.0  # 300 / (4 - 1)

    def test_a_square_raster_is_not_reported_as_anisotropic(self):
        """The old error was 1/N per axis, so it invented anisotropy.

        A 100 um raster three wide and two tall came out 66.67 x 50.0.
        Only x reaches the store, so the disagreement never showed.
        """
        positions = [
            (0, 0, 0.1, 0.1),
            (0, 1, 0.2, 0.1),
            (0, 2, 0.3, 0.1),
            (0, 3, 0.1, 0.2),
            (0, 4, 0.2, 0.2),
            (0, 5, 0.3, 0.2),
        ]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS}

        grid = build_imaging_grid(mock_ml, "handle", func_types)

        assert grid.pixel_size_x == 100.0
        assert grid.pixel_size_y == 100.0

    def test_an_axis_with_one_position_reports_no_pitch(self):
        """One position has no interval, so 0.0 rather than a made-up number.

        The extractor turns that into ``pixel_size=None``, which becomes
        the "use --pixel-size" refusal. Aborted single-line scans take
        this path, so it must not become a division by zero.
        """
        positions = [
            (0, 0, 0.1, 0.05),
            (0, 1, 0.2, 0.05),
        ]
        mock_ml = self._setup_mock_ml(positions)
        func_types = {0: FunctionType.MS}

        grid = build_imaging_grid(mock_ml, "handle", func_types)

        assert grid.pixel_size_x == 100.0
        assert grid.pixel_size_y == 0.0


# ---------------------------------------------------------------------------
# Declared regions: gridded one by one, the way Bruker regions are
# ---------------------------------------------------------------------------

X_STEP, Y_STEP = 20.1, 20.0


def _pattern(*regions, x_step=X_STEP, y_step=Y_STEP):
    return StagePattern(
        x_step_um=x_step,
        y_step_um=y_step,
        regions=tuple(
            DeclaredRegion(i + 1, x0, y0, dx, dy)
            for i, (x0, y0, dx, dy) in enumerate(regions)
        ),
    )


def _line(x0_um, y_um, n, pitch=X_STEP, func=0, first_scan=0):
    """One line scan, in the DLL's mm."""
    return {
        (func, first_scan + k): _make_scan_info(
            (x0_um + k * pitch) / 1000.0, y_um / 1000.0
        )
        for k in range(n)
    }


def _two_lines():
    """Two 25-pixel line scans whose x phases differ by 0.46 of a pitch."""
    scan_map = _line(66500.0, 74100.0, 25)
    scan_map.update(_line(66670.0, 63940.0, 25, first_scan=25))
    return scan_map


LPC_BOXES = ((66500.0, 74100.0, 500.0, 10.0), (66670.0, 63940.0, 500.0, 10.0))


class TestRegionsAreGriddedOneByOne:
    def test_the_two_lines_are_refused_without_the_regions(self):
        with pytest.raises(ConversionRefused, match="x axis"):
            _grid_from_scan_map(_two_lines())

    def test_with_the_regions_they_share_one_grid(self, thyra_logs):
        with thyra_logs("thyra.readers.waters", logging.WARNING) as records:
            grid = _grid_from_scan_map(_two_lines(), pattern=_pattern(*LPC_BOXES))

        assert (grid.pixel_count_x, grid.pixel_count_y) == (33, 509)
        assert grid.pixel_size_x == pytest.approx(20.1)
        assert grid.pixel_size_y == pytest.approx(20.0)
        assert grid.pixel_size_source == ("stage positions", "method step")
        first, second = grid.regions
        assert (first.number, second.number) == (1, 2)
        assert (first.x_offset, first.y_offset) == (0, 508)
        assert (second.x_offset, second.y_offset) == (8, 0)
        assert first.x_shift_um == pytest.approx(0.0)
        assert second.x_shift_um == pytest.approx(9.2, abs=0.01)
        assert grid.x_lattice.origin == pytest.approx(66500.0)
        assert grid.y_lattice.origin == pytest.approx(63940.0)
        moved = [
            r.getMessage() for r in records if "placed on one grid" in r.getMessage()
        ]
        assert len(moved) == 1 and "region 2" in moved[0]

    def test_every_scan_is_placed_once_on_its_own_line(self):
        scan_map = _two_lines()
        grid = _grid_from_scan_map(scan_map, pattern=_pattern(*LPC_BOXES))

        placed = {key: grid.get_coordinates(info) for key, info in scan_map.items()}

        assert all(coords is not None for coords in placed.values())
        assert len(set(placed.values())) == 50
        assert [placed[(0, k)] for k in (0, 24)] == [(0, 508, 0), (24, 508, 0)]
        assert [placed[(0, 25 + k)] for k in (0, 24)] == [(8, 0, 0), (32, 0, 0)]
        assert {grid.region_of(scan_map[(0, k)]) for k in range(25)} == {1}
        assert {grid.region_of(scan_map[(0, 25 + k)]) for k in range(25)} == {2}

    def test_a_scan_of_a_function_outside_the_fit_is_placed_by_its_region(self):
        scan_map = _line(66500.0, 74100.0, 25)
        scan_map.update(_line(66670.0, 63940.0, 20, first_scan=25))
        scan_map.update(_line(66670.0 + 20 * X_STEP, 63940.0, 5, func=1))
        grid = _grid_from_scan_map(scan_map, [0], pattern=_pattern(*LPC_BOXES))

        tail = [grid.get_coordinates(scan_map[(1, k)]) for k in range(5)]

        assert tail == [(8 + 20 + k, 0, 0) for k in range(5)]
        assert all(grid.lies_on_raster(scan_map[(1, k)]) for k in range(5))

    def test_a_reading_off_its_regions_lattice_has_no_pixel(self):
        scan_map = _two_lines()
        grid = _grid_from_scan_map(scan_map, pattern=_pattern(*LPC_BOXES))
        # Region 1's box, but half a pitch off region 1's lattice.
        stray = _make_scan_info((66500.0 + 10.0) / 1000.0, 74100.0 / 1000.0)

        assert grid.lies_on_raster(stray) is False
        assert grid.get_coordinates(stray) is None


class TestOneRegionOfSeveral:
    @pytest.mark.parametrize(
        "region, row_um, x0_um", [(1, 74100.0, 66500.0), (2, 63940.0, 66670.0)]
    )
    def test_a_selected_region_is_its_own_raster(self, region, row_um, x0_um):
        grid = _grid_from_scan_map(
            _two_lines(), pattern=_pattern(*LPC_BOXES), region=region
        )

        assert (grid.pixel_count_x, grid.pixel_count_y) == (25, 1)
        assert grid.pixel_size_x == pytest.approx(20.1)
        assert grid.pixel_size_y == pytest.approx(20.0)
        assert grid.pixel_size_source == ("stage positions", "method step")
        assert grid.x_lattice.origin == pytest.approx(x0_um)
        assert grid.y_lattice.origin == pytest.approx(row_um)
        assert grid.selected_region == region

    def test_scans_of_the_other_region_have_no_pixel_and_no_warning(self, thyra_logs):
        scan_map = _two_lines()
        grid = _grid_from_scan_map(scan_map, pattern=_pattern(*LPC_BOXES), region=1)

        with thyra_logs("thyra.readers.waters", logging.WARNING) as records:
            mine = grid.get_coordinates(scan_map[(0, 3)])
            other = grid.get_coordinates(scan_map[(0, 30)])

        assert mine == (3, 0, 0)
        assert other is None
        assert grid.lies_on_raster(scan_map[(0, 30)]) is False
        assert grid.region_of(scan_map[(0, 3)]) == 1
        assert not records

    def test_an_undeclared_region_is_refused_with_the_declared_ones(self):
        with pytest.raises(ConversionRefused, match=r"region=3.*1, 2"):
            _grid_from_scan_map(_two_lines(), pattern=_pattern(*LPC_BOXES), region=3)

    def test_a_declared_region_without_scans_is_refused(self):
        scan_map = _line(66500.0, 74100.0, 25)
        with pytest.raises(ConversionRefused, match="Region 2.*holds no scan"):
            _grid_from_scan_map(scan_map, pattern=_pattern(*LPC_BOXES), region=2)

    def test_a_regrid_keeps_the_region(self):
        grid = _grid_from_scan_map(_two_lines(), pattern=_pattern(*LPC_BOXES), region=2)

        assert regrid_for_functions(grid, [0]) is grid


class TestRunsWithoutUsableRegions:
    def _raster(self):
        scan_map = {}
        for row in range(3):
            scan_map.update(_line(1000.0, 2000.0 + row * Y_STEP, 6, first_scan=row * 6))
        return scan_map

    def test_region_1_is_the_whole_run_when_none_are_declared(self):
        grid = _grid_from_scan_map(self._raster(), pattern=_pattern(), region=1)

        assert (grid.pixel_count_x, grid.pixel_count_y) == (6, 3)
        assert grid.regions == ()
        assert grid.selected_region is None

    def test_any_other_region_is_refused_naming_region_1(self):
        with pytest.raises(ConversionRefused, match=r"region=2.*region 1"):
            _grid_from_scan_map(self._raster(), pattern=_pattern(), region=2)

    def test_no_pattern_at_all_is_the_same_run(self):
        grid = _grid_from_scan_map(self._raster(), region=1)

        assert (grid.pixel_count_x, grid.pixel_count_y) == (6, 3)
        with pytest.raises(ConversionRefused, match="region=2"):
            _grid_from_scan_map(self._raster(), region=2)

    def test_a_scan_in_no_box_drops_the_regions_and_keeps_the_old_refusal(
        self, thyra_logs
    ):
        scan_map = _two_lines()
        scan_map[(1, 0)] = _make_scan_info(70.0, 70.0)
        with pytest.raises(ConversionRefused) as old:
            _grid_from_scan_map(scan_map)

        with thyra_logs("thyra.readers.waters", logging.WARNING) as records:
            with pytest.raises(ConversionRefused) as refused:
                _grid_from_scan_map(scan_map, pattern=_pattern(*LPC_BOXES))

        assert str(refused.value) == str(old.value)
        notes = [r.getMessage() for r in records if "lie in no" in r.getMessage()]
        assert len(notes) == 1
        assert "70000.00" in notes[0] and "no region's box" in notes[0]

    def test_a_scan_in_two_boxes_drops_the_regions(self, thyra_logs):
        overlapping = _pattern(
            (0.0, 0.0, 9000.0, 9000.0), (500.0, 500.0, 9000.0, 9000.0)
        )

        with thyra_logs("thyra.readers.waters", logging.WARNING) as records:
            grid = _grid_from_scan_map(self._raster(), pattern=overlapping)

        assert grid.regions == () and grid.declared == ()
        assert (grid.pixel_count_x, grid.pixel_count_y) == (6, 3)
        assert any("regions 1, 2" in r.getMessage() for r in records)

    def test_one_region_holding_every_scan_keeps_the_old_grid_and_a_map(self):
        scan_map = self._raster()
        pattern = _pattern((900.0, 1990.0, 200.0, 100.0), (9000.0, 9000.0, 10.0, 10.0))
        plain = _grid_from_scan_map(scan_map)

        grid = _grid_from_scan_map(scan_map, pattern=pattern)

        assert grid.dimensions == plain.dimensions
        assert grid.x_index_map == plain.x_index_map
        assert grid.y_index_map == plain.y_index_map
        assert [p.number for p in grid.regions] == [1]
        assert {grid.region_of(i) for i in scan_map.values()} == {1}
        assert grid.get_coordinates(scan_map[(0, 7)]) == plain.get_coordinates(
            scan_map[(0, 7)]
        )


class TestRegionPlacementRefusals:
    def test_regions_rastered_on_different_pitches_are_refused(self):
        scan_map = _line(0.0, 0.0, 6, pitch=20.0)
        scan_map.update(_line(0.0, 20.0, 6, pitch=20.0, first_scan=6))
        scan_map.update(_line(5000.0, 5000.0, 6, pitch=30.0, first_scan=12))
        scan_map.update(_line(5000.0, 5030.0, 6, pitch=30.0, first_scan=18))
        pattern = _pattern(
            (0.0, 0.0, 100.0, 20.0),
            (5000.0, 5000.0, 150.0, 30.0),
            x_step=None,
            y_step=None,
        )

        with pytest.raises(ConversionRefused, match=r"one pitch.*region 1.*region 2"):
            _grid_from_scan_map(scan_map, pattern=pattern)

    def test_single_row_regions_without_a_step_are_refused(self):
        scan_map = _line(0.0, 0.0, 6, pitch=20.0)
        scan_map.update(_line(0.0, 5000.0, 6, pitch=20.0, first_scan=6))
        pattern = _pattern(
            (0.0, 0.0, 100.0, 1.0), (0.0, 5000.0, 100.0, 1.0), y_step=None
        )

        with pytest.raises(ConversionRefused, match=r"no step.*region="):
            _grid_from_scan_map(scan_map, pattern=pattern)

    def test_overlapping_spot_regions_fall_back_to_the_single_pixel_refusal(self):
        scan_map = {
            (0, 0): _make_scan_info(1.0, 2.0),
            (0, 1): _make_scan_info(1.0, 2.0),
        }
        pattern = _pattern((990.0, 1990.0, 20.0, 20.0), (990.0, 1990.0, 20.0, 20.0))

        with pytest.raises(ConversionRefused, match="single pixel"):
            _grid_from_scan_map(scan_map, pattern=pattern)


class TestPixelSizeOfASinglePositionAxis:
    def _line_scan(self):
        return _line(1000.0, 2000.0, 6, pitch=30.0)

    def test_the_declared_step_stands_in(self):
        grid = _grid_from_scan_map(self._line_scan(), pattern=_pattern())

        assert grid.pixel_size_x == pytest.approx(30.0)
        assert grid.pixel_size_y == pytest.approx(Y_STEP)
        assert grid.pixel_size_source == ("stage positions", "method step")
        assert grid.y_lattice.pitch == 0.0
        assert grid.lateral_height == 0.0

    def test_without_a_pattern_it_stays_zero(self):
        grid = _grid_from_scan_map(self._line_scan())

        assert grid.pixel_size_y == 0.0
        assert grid.pixel_size_source == ("stage positions", "none")

    def test_a_pattern_without_that_step_leaves_it_zero(self):
        grid = _grid_from_scan_map(self._line_scan(), pattern=_pattern(y_step=None))

        assert grid.pixel_size_y == 0.0
        assert grid.pixel_size_source == ("stage positions", "none")

    def test_a_single_spot_keeps_no_pixel_size(self):
        scan_map = {
            (0, 0): _make_scan_info(1.0, 2.0),
            (0, 1): _make_scan_info(1.0, 2.0),
        }

        grid = _grid_from_scan_map(
            scan_map, allow_single_position=True, pattern=_pattern()
        )

        assert (grid.pixel_size_x, grid.pixel_size_y) == (0.0, 0.0)
        assert grid.pixel_size_source == ("none", "none")
        assert grid.no_raster_reason is not None


class TestRegionRequestWhenTheRegionsAreDropped:
    @pytest.mark.parametrize("region", [1, 2])
    def test_any_region_request_is_refused_not_read_as_the_whole_run(self, region):
        scan_map = _two_lines()
        scan_map[(1, 0)] = _make_scan_info(70.0, 70.0)

        with pytest.raises(ConversionRefused) as refused:
            _grid_from_scan_map(scan_map, pattern=_pattern(*LPC_BOXES), region=region)

        message = str(refused.value)
        assert f"region={region}" in message
        assert "declares 2 region(s)" in message
        assert "could not be used" in message
        assert "no region's box" in message
        assert "converted whole" in message
        assert "no declared regions" not in message


class TestARegionOnlyInAnotherFunction:
    def _scans(self):
        scan_map = _two_lines()
        # Region 3 sits elsewhere on the slide and only function 1 holds it.
        scan_map.update(_line(66500.0, 70000.0, 10, func=1))
        return scan_map

    def _pattern3(self):
        return _pattern(*LPC_BOXES, (66500.0, 70000.0, 300.0, 10.0))

    def test_its_scans_count_as_on_the_raster_but_have_no_pixel_yet(self):
        scan_map = self._scans()
        grid = _grid_from_scan_map(scan_map, [0], pattern=self._pattern3())

        scan = scan_map[(1, 0)]

        assert grid.lies_on_raster(scan) is True
        assert grid.locate(scan) is None
        assert grid.region_of(scan) == 3

    def test_the_regrid_places_the_region(self):
        scan_map = self._scans()
        grid = _grid_from_scan_map(scan_map, [0], pattern=self._pattern3())

        grown = regrid_for_functions(grid, [0, 1])

        assert [p.number for p in grown.regions] == [1, 2, 3]
        coords = [grown.get_coordinates(scan_map[(1, k)]) for k in range(10)]
        assert all(c is not None for c in coords)
        assert len(set(coords)) == 10
        assert len(grid.regions) == 2

    def test_a_request_for_another_region_leaves_it_off_the_raster(self):
        scan_map = self._scans()
        grid = _grid_from_scan_map(scan_map, [0], pattern=self._pattern3(), region=1)

        assert grid.lies_on_raster(scan_map[(1, 0)]) is False


class TestAScanInNoBoxIsNotOnTheRaster:
    def _scans(self):
        scan_map = _line(66500.0, 74100.0, 25)
        # Function 1 is outside the fitted census and outside every box, but
        # on region 1's lattice.
        scan_map[(1, 0)] = _make_scan_info((66500.0 + 5 * X_STEP) / 1000.0, 80.0)
        return scan_map

    @pytest.mark.parametrize("region", [None, 1])
    def test_it_has_no_region_and_no_pixel(self, region):
        scan_map = self._scans()
        grid = _grid_from_scan_map(
            scan_map, [0], pattern=_pattern(*LPC_BOXES), region=region
        )

        stray = scan_map[(1, 0)]

        assert grid.region_of(stray) is None
        assert grid.lies_on_raster(stray) is False
        assert grid.locate(stray) is None
        assert grid.get_coordinates(stray) is None

    def test_the_whole_run_and_region_1_give_the_same_pixels(self):
        scan_map = self._scans()
        pattern = _pattern(*LPC_BOXES)
        whole = _grid_from_scan_map(scan_map, [0], pattern=pattern)
        one = _grid_from_scan_map(scan_map, [0], pattern=pattern, region=1)

        for key, info in scan_map.items():
            assert whole.get_coordinates(info) == one.get_coordinates(info)
        assert whole.dimensions == one.dimensions
