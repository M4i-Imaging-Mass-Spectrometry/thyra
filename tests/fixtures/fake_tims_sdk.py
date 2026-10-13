"""A stand-in for the two vendor calls the precursor split makes.

Hands the TDF reader a fixed ``(index, intensity, scan)`` point cloud whose
answer is known by construction, so the split runs in CI with no SDK.
Shared by the PASEF and the DIA (type-9) split tests.
"""

import numpy as np

# (index, intensity, scan) triples every frame is given, ordered by scan
# the way tims_read_scans_v2 returns them. Two pairs share a window and a
# digitizer index (scans 25 and 30), one sits in the gap between two
# windows (scan 70), and one sits past the last window (scan 220).
POINTS = [
    (10, 100, 25),
    (10, 5, 30),
    (40, 70, 55),
    (99, 999, 70),
    (20, 200, 90),
    (30, 300, 125),
    (40, 400, 160),
    (10, 50, 199),
    (77, 777, 220),
]


class FakeSDK:
    """The two vendor calls the demultiplexer makes, over a fixed cloud."""

    def __init__(self, points=POINTS):
        self.points = points
        self.scan_ranges = []

    def read_tdf_scans(self, handle, frame_id, scan_begin, scan_end, hint=None):
        self.scan_ranges.append((frame_id, scan_begin, scan_end))
        indices, intensities, scans = (
            np.array([p[i] for p in self.points], dtype=dtype)
            for i, dtype in enumerate((np.uint32, np.uint32, np.int32))
        )
        return indices, intensities, scans

    def index_to_mz(self, handle, frame_id, indices):
        # Monotonic in the index, which is all the reader relies on.
        return 100.0 + np.asarray(indices, dtype=np.float64) * 0.5

    def scannum_to_oneoverk0(self, handle, frame_id, scans):
        return 2.0 - 0.001 * np.asarray(scans, dtype=np.float64)

    def close_file(self, handle):
        pass
