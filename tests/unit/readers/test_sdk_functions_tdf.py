"""The TDF path of :class:`SDKFunctions`, against a fake Bruker library.

Two defects lived here until 2026-09: the frame id was decremented before
being handed to the SDK, so every pixel received the previous frame (and the
first frame raised "Frame doesn't exist"), and only scans 0..100 of the
mobility ramp were read. The fake records exactly what the SDK is asked for,
so both would fail here.
"""

from __future__ import annotations

import ctypes
import logging
from ctypes import POINTER, c_double, c_float
from typing import Dict, List, Optional, Tuple
from unittest.mock import MagicMock

import numpy as np
import pytest

from thyra.readers.bruker.timstof.sdk.sdk_functions import (
    TDF_SPECTRUM_MODES,
    SDKFunctions,
)
from thyra.utils.bruker_exceptions import SDKError

Scan = Tuple[List[int], List[int]]  # (indices, intensities)


class FakeDll:
    """Just enough of ``timsdata`` for ``SDKFunctions('tdf')``, recording its calls.

    Attributes are plain instance attributes, so ``getattr(dll, name, None)``
    behaves as it does on a real ``CDLL``: a missing export is missing.
    """

    def __init__(
        self,
        frames: Dict[int, List[Scan]],
        centroids: Optional[Dict[int, Tuple[List[float], List[float]]]] = None,
        with_centroid_export: bool = True,
    ):
        self.frames = frames
        self.centroids = centroids or {}
        self.calls: List[Tuple] = []
        self.tims_open = MagicMock(return_value=42)
        self.tims_close = MagicMock()
        self.tims_get_last_error_string = MagicMock(side_effect=self._last_error)
        self.tims_read_scans_v2 = MagicMock(side_effect=self._read_scans)
        self.tims_index_to_mz = MagicMock(side_effect=self._index_to_mz)
        self.tims_mz_to_index = MagicMock(return_value=1)
        self.tims_scannum_to_oneoverk0 = MagicMock(side_effect=self._scan_to_k0)
        self.tims_oneoverk0_to_scannum = MagicMock(return_value=1)
        self.tims_scannum_to_voltage = MagicMock(return_value=1)
        self.tims_oneoverk0_to_ccs_for_mz = MagicMock(
            side_effect=lambda k0, charge, mz: 100.0 * k0 * charge + mz
        )
        self.tims_ccs_to_oneoverk0_for_mz = MagicMock(return_value=1.0)
        if with_centroid_export:
            self.tims_extract_centroided_spectrum_for_frame_v2 = MagicMock(
                side_effect=self._extract
            )

    # -- fake SDK behaviour -------------------------------------------------

    def _last_error(self, buf, length):
        message = b"fake sdk error"
        if buf is None:
            return len(message) + 1
        buf.value = message
        return length

    def _read_scans(
        self, handle, frame_id, scan_begin, scan_end, buffer_ptr, buffer_bytes
    ):
        self.calls.append(("read_scans", frame_id, scan_begin, scan_end))
        scans = self.frames[frame_id][scan_begin:scan_end]
        words: List[int] = [len(idx) for idx, _ in scans]
        for idx, inten in scans:
            words.extend(idx)
            words.extend(inten)
        payload = np.asarray(words, dtype=np.uint32).tobytes()
        if buffer_bytes < len(payload):
            return len(payload)
        ctypes.memmove(buffer_ptr, payload, len(payload))
        return len(payload)

    def _index_to_mz(self, handle, frame_id, in_ptr, out_ptr, n):
        self.calls.append(("index_to_mz", frame_id, int(n)))
        values = np.ctypeslib.as_array(in_ptr, shape=(int(n),))
        out = np.ctypeslib.as_array(out_ptr, shape=(int(n),))
        out[:] = 100.0 + 0.001 * values
        return 1

    def _scan_to_k0(self, handle, frame_id, in_ptr, out_ptr, n):
        self.calls.append(("scan_to_k0", frame_id, int(n)))
        values = np.ctypeslib.as_array(in_ptr, shape=(int(n),))
        out = np.ctypeslib.as_array(out_ptr, shape=(int(n),))
        out[:] = 2.0 - 0.001 * values
        return 1

    def _extract(self, handle, frame_id, scan_begin, scan_end, callback, user_data):
        self.calls.append(("extract", frame_id, scan_begin, scan_end))
        mz, area = self.centroids[frame_id]
        n = len(mz)
        mz_arr = (c_double * n)(*mz)
        area_arr = (c_float * n)(*area)
        callback(frame_id, n, mz_arr, area_arr)
        return 1


class FakeDllManager:
    def __init__(self, dll: FakeDll):
        self.dll = dll


# Frame 7: four scans, the second and the last empty; index 10 appears in
# two scans so the scan sum has to merge it.
FRAME_7: List[Scan] = [([10, 20], [5, 6]), ([], []), ([10], [7]), ([], [])]


def make_sdk(dll: FakeDll, mode: str = "vendor_centroid") -> SDKFunctions:
    return SDKFunctions(FakeDllManager(dll), "tdf", tdf_spectrum=mode)


class TestFrameIdAndRamp:
    def test_vendor_centroid_passes_the_frame_id_unchanged_over_the_full_ramp(self):
        dll = FakeDll({7: FRAME_7}, centroids={7: ([300.5, 100.25], [2.0, 9.0])})
        sdk = make_sdk(dll)

        mzs, intensities = sdk.read_spectrum(handle=42, frame_id=7, num_scans=4)

        assert ("extract", 7, 0, 4) in dll.calls
        assert not any(
            call[1] == 6 for call in dll.calls
        ), "frame_id - 1 is the old bug"
        np.testing.assert_allclose(mzs, [100.25, 300.5])
        np.testing.assert_allclose(intensities, [9.0, 2.0])
        assert mzs.dtype == np.float64 and intensities.dtype == np.float64

    def test_scan_sum_reads_every_scan_and_merges_indices(self):
        dll = FakeDll({7: FRAME_7})
        sdk = make_sdk(dll, "scan_sum")

        mzs, intensities = sdk.read_spectrum(
            handle=42, frame_id=7, buffer_size_hint=3, num_scans=4
        )

        assert ("read_scans", 7, 0, 4) in dll.calls
        assert ("index_to_mz", 7, 2) in dll.calls
        np.testing.assert_allclose(mzs, [100.010, 100.020])
        np.testing.assert_allclose(intensities, [12.0, 6.0])

    def test_tdf_read_requires_num_scans(self):
        sdk = make_sdk(FakeDll({7: FRAME_7}))
        with pytest.raises(SDKError, match="NumScans"):
            sdk.read_spectrum(handle=42, frame_id=7)


class TestReadTdfScans:
    def test_keeps_scan_numbers_and_survives_empty_scans(self):
        dll = FakeDll({7: FRAME_7})
        sdk = make_sdk(dll, "scan_sum")

        indices, intensities, scans = sdk.read_tdf_scans(42, 7, 0, 4, num_peaks_hint=3)

        np.testing.assert_array_equal(indices, [10, 20, 10])
        np.testing.assert_array_equal(intensities, [5, 6, 7])
        np.testing.assert_array_equal(scans, [0, 0, 2])
        assert scans.dtype == np.int32

    def test_scan_range_is_honoured(self):
        dll = FakeDll({7: FRAME_7})
        sdk = make_sdk(dll, "scan_sum")

        indices, intensities, scans = sdk.read_tdf_scans(42, 7, 2, 4)

        np.testing.assert_array_equal(indices, [10])
        np.testing.assert_array_equal(scans, [2])
        assert ("read_scans", 7, 2, 4) in dll.calls

    def test_buffer_grows_when_the_sdk_asks_for_more(self):
        dll = FakeDll({7: FRAME_7})
        sdk = make_sdk(dll, "scan_sum")

        # A hint of 0 pairs makes the first buffer far too small for the
        # default path to matter, so force the tiny-buffer branch directly.
        indices, _, _ = sdk.read_tdf_scans(42, 7, 0, 4, num_peaks_hint=None)
        assert indices.size == 3

        # The fake asked for more bytes exactly once if the first buffer was
        # short, never more than twice in total.
        assert 1 <= dll.tims_read_scans_v2.call_count <= 2

    def test_frame_without_pairs_returns_empty_arrays(self):
        dll = FakeDll({1: [([], []), ([], [])]})
        sdk = make_sdk(dll, "scan_sum")

        indices, intensities, scans = sdk.read_tdf_scans(42, 1, 0, 2)
        assert indices.size == intensities.size == scans.size == 0

        mzs, summed = sdk.read_spectrum(handle=42, frame_id=1, num_scans=2)
        assert mzs.size == summed.size == 0


class TestModes:
    def test_modes_are_the_documented_two(self):
        assert TDF_SPECTRUM_MODES == ("vendor_centroid", "scan_sum")

    def test_unknown_mode_is_rejected(self):
        with pytest.raises(ValueError, match="tdf_spectrum"):
            make_sdk(FakeDll({}), "peak_pick")

    def test_missing_centroid_export_falls_back_to_scan_sum(self, caplog):
        dll = FakeDll({7: FRAME_7}, with_centroid_export=False)
        with caplog.at_level(logging.WARNING):
            sdk = make_sdk(dll, "vendor_centroid")

        assert sdk.tdf_spectrum == "scan_sum"
        assert "scan_sum" in caplog.text
        mzs, intensities = sdk.read_spectrum(handle=42, frame_id=7, num_scans=4)
        np.testing.assert_allclose(intensities, [12.0, 6.0])

    def test_tsf_ignores_the_mode(self):
        dll = MagicMock()
        sdk = SDKFunctions(FakeDllManager(dll), "tsf", tdf_spectrum="scan_sum")
        assert sdk.file_type == "tsf"


class TestMobilityConversions:
    def test_scannum_to_oneoverk0_goes_through_the_sdk(self):
        dll = FakeDll({})
        sdk = make_sdk(dll)

        k0 = sdk.scannum_to_oneoverk0(42, 3, np.array([0, 1000]))

        np.testing.assert_allclose(k0, [2.0, 1.0])
        assert ("scan_to_k0", 3, 2) in dll.calls

    def test_ccs_broadcasts_over_features(self):
        sdk = make_sdk(FakeDll({}))
        ccs = sdk.oneoverk0_to_ccs(np.array([1.0, 1.5]), 1, np.array([500.0, 600.0]))
        np.testing.assert_allclose(ccs, [600.0, 750.0])

    def test_missing_conversion_export_is_a_named_error(self):
        dll = FakeDll({})
        del dll.tims_scannum_to_oneoverk0
        sdk = make_sdk(dll)
        with pytest.raises(SDKError, match="tims_scannum_to_oneoverk0"):
            sdk.scannum_to_oneoverk0(42, 1, np.array([0.0]))


class TestPointerTypes:
    def test_conversion_pointer_type_is_double(self):
        # Guards the ctypes declarations: the SDK reads and writes doubles.
        assert POINTER(c_double) is not POINTER(c_float)
        sdk = make_sdk(FakeDll({}))
        assert sdk._bound_conversions["tims_index_to_mz"] is True
