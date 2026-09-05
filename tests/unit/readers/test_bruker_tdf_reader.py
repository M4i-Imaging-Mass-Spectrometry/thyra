"""BrukerReader's TDF-specific wiring, with the SDK and the database mocked."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from thyra.readers.bruker.timstof.timstof_reader import BrukerReader


def _only_tdf_exists(self: Path) -> bool:
    """``Path.exists`` that makes a fake ``.d`` look like a TDF acquisition."""
    return not str(self).endswith("analysis.tsf")


def _make_reader(file_type: str, **kwargs) -> tuple:
    exists = _only_tdf_exists if file_type == "tdf" else (lambda self: True)
    with (
        patch("thyra.readers.bruker.timstof.timstof_reader.DLLManager") as dll_manager,
        patch(
            "thyra.readers.bruker.timstof.timstof_reader.SDKFunctions"
        ) as sdk_functions,
        patch.object(Path, "exists", new=exists),
        patch.object(Path, "is_dir", return_value=True),
        patch("sqlite3.connect") as connect,
    ):
        dll_manager.return_value = MagicMock()
        sdk = MagicMock()
        sdk_functions.return_value = sdk
        sdk.open_file.return_value = 42
        sdk.read_spectrum.return_value = (np.array([100.0]), np.array([1.0]))
        connect.return_value = MagicMock()
        reader = BrukerReader(Path("/fake/bruker.d"), **kwargs)
    return reader, sdk, sdk_functions


class TestTdfWiring:
    def test_tdf_read_passes_the_frame_id_and_its_num_scans(self):
        reader, sdk, _ = _make_reader("tdf")
        assert reader.file_type == "tdf"
        reader._num_peaks_cache = {3: 500}
        reader._num_scans_cache = {3: 2382}

        reader._read_frame_spectrum(3)

        sdk.read_spectrum.assert_called_once_with(
            42, 3, buffer_size_hint=500, num_scans=2382
        )

    def test_num_scans_is_fetched_from_the_database_when_not_cached(self):
        reader, sdk, _ = _make_reader("tdf")
        reader._num_scans_cache = {}
        reader.conn = MagicMock()
        reader.conn.execute.return_value.fetchone.return_value = (9400,)

        reader._read_frame_spectrum(11)

        sdk.read_spectrum.assert_called_once_with(
            42, 11, buffer_size_hint=None, num_scans=9400
        )
        assert reader._num_scans_cache[11] == 9400

    def test_tsf_read_passes_no_num_scans(self):
        reader, sdk, _ = _make_reader("tsf")
        assert reader.file_type == "tsf"
        reader._num_peaks_cache = {3: 500}

        reader._read_frame_spectrum(3)

        sdk.read_spectrum.assert_called_once_with(
            42, 3, buffer_size_hint=500, num_scans=None
        )

    def test_tdf_spectrum_mode_reaches_the_sdk(self):
        reader, _, sdk_functions = _make_reader("tdf", tdf_spectrum="scan_sum")
        assert reader.tdf_spectrum == "scan_sum"
        sdk_functions.assert_called_once_with(
            reader.dll_manager, "tdf", tdf_spectrum="scan_sum"
        )

    def test_default_mode_is_the_vendor_centroid(self):
        reader, _, _ = _make_reader("tdf")
        assert reader.tdf_spectrum == "vendor_centroid"

    def test_unknown_mode_is_rejected_before_touching_the_sdk(self):
        with pytest.raises(ValueError, match="tdf_spectrum"):
            _make_reader("tdf", tdf_spectrum="bogus")

    def test_tdf_reports_no_per_pixel_peak_counts(self):
        reader, _, _ = _make_reader("tdf")
        reader._num_peaks_cache = {1: 70000, 2: 80000}
        assert reader.get_peak_counts_per_pixel() is None

    def test_num_peaks_above_65535_are_kept(self):
        reader, _, _ = _make_reader("tdf")
        conn = MagicMock()
        cursor = conn.__enter__.return_value.cursor.return_value
        cursor.fetchall.return_value = [(1, 70000, 477), (2, 0, 477), (3, 96845, 477)]
        with patch("sqlite3.connect", return_value=conn):
            cache = reader._preload_frame_num_peaks()

        assert cache == {1: 70000, 3: 96845}
        assert reader._num_scans_cache == {1: 477, 2: 477, 3: 477}
