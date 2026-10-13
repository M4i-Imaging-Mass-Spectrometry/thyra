"""The drift-bin bindings of the MassLynx wrapper, without the native library.

The wrapper is built around a stand-in for ``MLReader``: the tests check the
signatures it declares and the calls it makes, which is what decides whether
a real library call is safe.
"""

from ctypes import POINTER, c_double, c_float, c_int32, c_uint32, c_void_p
from unittest.mock import MagicMock

import numpy as np

from thyra.readers.waters.masslynx_lib import MassLynxLib


def _wrapper(api):
    """A MassLynxLib whose API library is ``api``, without loading anything."""
    ml = MassLynxLib.__new__(MassLynxLib)
    ml._base_lib = MagicMock()
    ml._api_lib = api
    return ml


class _Api:
    """Answers the per-bin read the way MLReader does: count, then fill."""

    def __init__(self, n_points):
        self.n_points = n_points
        self.calls = []

    def getRawMobilityScanDataPoints(
        self, handle, func, scan, drift, mz, inten, nbytes
    ):
        self.calls.append((func, scan, drift, nbytes))
        if self.n_points * 8 <= nbytes:
            for i in range(self.n_points):
                mz[i] = 100.0 + i
                inten[i] = float(drift + 1)
        return self.n_points


class TestDriftBindings:
    def test_signatures_are_declared(self):
        api = MagicMock()
        _wrapper(api)._setup_function_signatures()

        assert api.getNumberOfMobilityScans.argtypes == [c_void_p, c_int32]
        assert api.getNumberOfMobilityScans.restype is c_uint32
        assert api.getMobility.argtypes == [c_void_p, c_int32]
        assert api.getMobility.restype is c_float
        assert api.getRawMobilityScanDataPoints.argtypes == [
            c_void_p,
            c_int32,
            c_int32,
            c_int32,
            POINTER(c_double),
            POINTER(c_double),
            c_uint32,
        ]
        assert api.getRawMobilityScanDataPoints.restype is c_uint32

    def test_drift_times_never_ask_past_the_count(self):
        api = MagicMock()
        api.getMobility.side_effect = lambda _h, i: 0.05 * i

        times = _wrapper(api).get_drift_times("handle", 4)

        np.testing.assert_allclose(times, [0.0, 0.05, 0.10, 0.15])
        assert [c.args[1] for c in api.getMobility.call_args_list] == [0, 1, 2, 3]

    def test_a_bin_read_fits_the_buffer_first_time(self):
        api = _Api(3)

        mzs, intensities = _wrapper(api).read_drift_bin("handle", 0, 5, 2)

        np.testing.assert_array_equal(mzs, [100.0, 101.0, 102.0])
        np.testing.assert_array_equal(intensities, [3.0, 3.0, 3.0])
        assert len(api.calls) == 1
        assert api.calls[0][:3] == (0, 5, 2)

    def test_a_bin_larger_than_the_buffer_is_read_again(self):
        api = _Api(10)

        mzs, intensities = _wrapper(api).read_drift_bin(
            "handle", 1, 0, 0, initial_buf_size=4
        )

        assert mzs.size == 10 and intensities.size == 10
        np.testing.assert_array_equal(mzs, 100.0 + np.arange(10))
        assert len(api.calls) == 2
        assert api.calls[1][3] >= 10 * 8

    def test_an_empty_bin_is_two_empty_arrays(self):
        mzs, intensities = _wrapper(_Api(0)).read_drift_bin("handle", 0, 0, 0)

        assert mzs.size == 0 and intensities.size == 0
        assert mzs.dtype == np.float64
