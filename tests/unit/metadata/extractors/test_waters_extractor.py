"""Tests for the Waters metadata extractor."""

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from thyra.metadata.extractors.waters_extractor import WatersMetadataExtractor
from thyra.readers.waters.imaging_grid import ImagingGrid
from thyra.readers.waters.masslynx_lib import FunctionType, ScanInfoData


def _make_scan_info(x_mm, y_mm, has_pos=True):
    """Create a ScanInfoData."""
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


def _make_grid_and_ml(n_x=3, n_y=2, n_peaks=5):
    """Create a mock MassLynxLib, handle, and ImagingGrid for testing.

    Returns:
        Tuple of (mock_ml, handle, grid, function_types, ms_functions).
    """
    # Build x/y index maps
    x_positions = [100.0 * (i + 1) for i in range(n_x)]
    y_positions = [100.0 * (i + 1) for i in range(n_y)]
    x_index_map = {v: i for i, v in enumerate(x_positions)}
    y_index_map = {v: i for i, v in enumerate(y_positions)}

    # Build scan map: 1 MS function, n_x * n_y scans
    scan_map = {}
    scan_idx = 0
    for yi, y_um in enumerate(y_positions):
        for xi, x_um in enumerate(x_positions):
            x_mm = x_um / 1000.0
            y_mm = y_um / 1000.0
            scan_map[(0, scan_idx)] = _make_scan_info(x_mm, y_mm)
            scan_idx += 1

    n_scans = n_x * n_y

    lateral_width = x_positions[-1] - x_positions[0] if n_x > 1 else 0.0
    lateral_height = y_positions[-1] - y_positions[0] if n_y > 1 else 0.0
    pixel_size_x = lateral_width / n_x if n_x > 1 else 0.0
    pixel_size_y = lateral_height / n_y if n_y > 1 else 0.0

    grid = ImagingGrid(
        x_index_map=x_index_map,
        y_index_map=y_index_map,
        pixel_count_x=n_x,
        pixel_count_y=n_y,
        pixel_size_x=pixel_size_x,
        pixel_size_y=pixel_size_y,
        lateral_width=lateral_width,
        lateral_height=lateral_height,
        scan_map=scan_map,
    )

    # Mock MassLynxLib
    mock_ml = MagicMock()
    mock_ml.get_number_of_scans_in_function.return_value = n_scans
    mock_ml.is_raw_spectrum_profile.return_value = False  # centroid
    mock_ml.get_acquisition_date.return_value = "2026-01-15"
    mock_ml.is_lockmass_corrected.return_value = False
    mock_ml.get_lockmass_function.return_value = -1
    mock_ml.get_number_of_functions.return_value = 1
    mock_ml.get_acquisition_range.return_value = (100.0, 1000.0)

    # Each spectrum returns n_peaks m/z values
    mzs = np.linspace(100.0, 1000.0, n_peaks)
    intensities = np.random.default_rng(42).uniform(10.0, 500.0, n_peaks)
    mock_ml.read_spectrum.return_value = (mzs, intensities)

    handle = "mock_handle"
    function_types = {0: FunctionType.MS}
    ms_functions = [0]

    return mock_ml, handle, grid, function_types, ms_functions


class TestWatersMetadataExtractorEssential:
    """Test essential metadata extraction."""

    def test_dimensions(self):
        """Test that dimensions match the imaging grid."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml(n_x=3, n_y=2)
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.dimensions == (3, 2, 1)

    def test_coordinate_bounds(self):
        """Test coordinate bounds are 0-based pixel indices."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml(n_x=4, n_y=3)
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.coordinate_bounds == (0.0, 3.0, 0.0, 2.0)

    def test_mass_range(self):
        """Test mass range from scanned spectra."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml(n_x=2, n_y=2, n_peaks=10)
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.mass_range[0] == pytest.approx(100.0)
        assert essential.mass_range[1] == pytest.approx(1000.0)

    def test_spectra_count(self):
        """Test that n_spectra matches the number of positioned scans."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml(n_x=3, n_y=2)
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.n_spectra == 6  # 3 * 2

    def test_total_peaks(self):
        """Test total peak count."""
        n_peaks = 20
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml(n_x=2, n_y=2, n_peaks=n_peaks)
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.total_peaks == 4 * n_peaks  # 4 pixels * 20 peaks each

    def test_pixel_size(self):
        """Test pixel size from imaging grid."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml(n_x=3, n_y=2)
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.pixel_size is not None
        assert essential.pixel_size[0] > 0
        assert essential.pixel_size[1] > 0

    def test_spectrum_type_centroid(self):
        """Test centroid spectrum type detection."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        mock_ml.is_raw_spectrum_profile.return_value = False
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.spectrum_type == "centroid spectrum"

    def test_spectrum_type_profile(self):
        """Test profile spectrum type detection."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        mock_ml.is_raw_spectrum_profile.return_value = True
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.spectrum_type == "profile spectrum"

    def test_source_path(self):
        """Test source path is preserved."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.source_path == str(Path("/test/data.raw"))

    def test_caching(self):
        """Test that essential metadata is cached."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        e1 = extractor.get_essential()
        e2 = extractor.get_essential()

        assert e1 is e2

    def test_peak_counts_per_pixel(self):
        """Test per-pixel peak count array."""
        n_peaks = 15
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml(n_x=2, n_y=2, n_peaks=n_peaks)
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        assert essential.peak_counts_per_pixel is not None
        assert len(essential.peak_counts_per_pixel) == 4
        assert all(c == n_peaks for c in essential.peak_counts_per_pixel)


class TestWatersMetadataExtractorComprehensive:
    """Test comprehensive metadata extraction."""

    def test_comprehensive_contains_essential(self):
        """Test that comprehensive metadata includes essential."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        comprehensive = extractor.get_comprehensive()

        assert comprehensive.essential is not None
        assert comprehensive.essential.dimensions == (3, 2, 1)

    def test_format_specific(self):
        """Test Waters-specific metadata fields."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        comprehensive = extractor.get_comprehensive()

        fs = comprehensive.format_specific
        assert fs["data_format"] == "waters_raw"
        assert fs["is_imaging"] is True
        assert fs["pixel_count_x"] == 3
        assert fs["pixel_count_y"] == 2
        assert fs["ms_functions"] == [0]

    def test_instrument_info(self):
        """Test instrument info reports Waters vendor."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        comprehensive = extractor.get_comprehensive()

        assert comprehensive.instrument_info["vendor"] == "Waters"

    def test_acquisition_params(self):
        """Test acquisition parameters extraction."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        comprehensive = extractor.get_comprehensive()

        params = comprehensive.acquisition_params
        assert "acquisition_date" in params
        assert params["acquisition_date"] == "2026-01-15"
        assert "is_lockmass_corrected" in params

    def test_raw_metadata_imaging_grid(self):
        """Test raw metadata contains imaging grid details."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml()
        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        comprehensive = extractor.get_comprehensive()

        raw = comprehensive.raw_metadata
        assert "imaging_grid" in raw
        assert raw["imaging_grid"]["pixel_count_x"] == 3
        assert raw["imaging_grid"]["pixel_count_y"] == 2


class TestWatersMetadataExtractorEdgeCases:
    """Test edge cases in metadata extraction."""

    def test_read_spectrum_error_handled(self):
        """Test that spectrum read errors during metadata scan are handled."""
        mock_ml, handle, grid, ft, ms = _make_grid_and_ml(n_x=2, n_y=1)

        # First spectrum succeeds, second fails
        mzs = np.array([100.0, 500.0, 1000.0])
        ints = np.array([10.0, 20.0, 30.0])
        mock_ml.read_spectrum.side_effect = [
            (mzs, ints),
            RuntimeError("DLL crash"),
        ]

        extractor = WatersMetadataExtractor(
            mock_ml, handle, Path("/test/data.raw"), grid, ft, ms
        )
        essential = extractor.get_essential()

        # Should still have data from the first spectrum
        assert essential.n_spectra == 1
        assert essential.total_peaks == 3
