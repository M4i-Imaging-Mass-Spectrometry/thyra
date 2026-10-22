"""Waters lock-mass correction and repeated pixels against a real run (D83, D84).

These tests run only when ``THYRA_WATERS_LOCK_MASS_DATASET`` names the
MassIVE run MSV000092638 ``20210426_9AA_WTBrain1-1.raw`` (SYNAPT XS, MALDI,
negative mode, a red phosphorus lock-mass spot), or a run like it: a brain
section in negative mode whose lock-mass scans hold the P17- cluster anion.
With the variable unset the module skips. The vendor library is bundled for
Windows and Linux only; anywhere it cannot be loaded the tests skip.
"""

import hashlib
import os
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pytest

from thyra.readers.waters.waters_reader import WatersReader

DATASET = os.environ.get("THYRA_WATERS_LOCK_MASS_DATASET")

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not DATASET,
        reason="Set THYRA_WATERS_LOCK_MASS_DATASET to a Waters .raw run to run",
    ),
]

#: The P17- red phosphorus cluster anion.
P17 = 17 * 30.97376163 + 0.00054858

#: [M-H]- of lipids abundant in brain, by monoisotopic mass.
LIPIDS = {"PE 38:4": 766.53923, "PI 38:4": 885.54985, "SHexCer 42:2;O2": 888.62401}


def _open(lock_mass) -> WatersReader:
    try:
        reader = WatersReader(Path(DATASET), lock_mass=lock_mass)  # type: ignore[arg-type]
        reader._ensure_initialized()
    except OSError as exc:  # the vendor library is not loadable here
        pytest.skip(f"Waters library not loadable on this platform: {exc}")
    return reader


def _fingerprint(folder: Path) -> Dict[str, Tuple[int, int, str]]:
    prints = {}
    for path in sorted(folder.iterdir()):
        if path.is_file():
            stat = path.stat()
            small = stat.st_size < 50_000_000
            digest = hashlib.sha256(path.read_bytes()).hexdigest() if small else ""
            prints[path.name] = (stat.st_size, stat.st_mtime_ns, digest)
    return prints


def _lipid_error_by_third(reader: WatersReader) -> List[float]:
    """Mean lipid m/z error (ppm) in the first, middle and last third of the run."""
    ml, handle, grid, _types, _functions = reader._require_initialized()
    n = ml.get_number_of_scans_in_function(handle, 0)
    thirds: List[List[float]] = [[], [], []]
    for third, start in enumerate((0, n // 3, 2 * n // 3)):
        mzs, ints = [], []
        for scan in range(start, start + n // 3, 40):
            info = grid.scan_map[(0, scan)]
            if reader._scan_weight(0, scan) == 0.0:
                continue
            mz, it = ml.read_spectrum(handle, 0, scan)
            mzs.append(mz * reader._mz_scale(info.rt))
            ints.append(it)
        a, w = np.concatenate(mzs), np.concatenate(ints)
        for theory in LIPIDS.values():
            near = (a > theory - 0.03) & (a < theory + 0.03)
            measured = (a[near] * w[near]).sum() / w[near].sum()
            thirds[third].append((measured - theory) / theory * 1e6)
    return [float(np.mean(t)) for t in thirds]


class TestLockMassOnARealRun:
    def test_the_run_folder_is_never_written_to(self):
        folder = Path(DATASET)  # type: ignore[arg-type]
        before = _fingerprint(folder)

        reader = _open(P17)
        reader.close()

        assert _fingerprint(folder) == before

    def test_the_correction_removes_the_drift(self):
        uncorrected = _open("off")
        drifting = _lipid_error_by_third(uncorrected)
        uncorrected.close()
        corrected = _open(P17)
        steady = _lipid_error_by_third(corrected)
        report = corrected.lock_mass_report()
        corrected.close()

        # Measured for D83: +12 ppm rising to +21 ppm uncorrected.
        assert max(drifting) - min(drifting) > 5.0
        assert max(steady) - min(steady) < 2.0
        assert all(abs(error) < 4.0 for error in steady)
        assert report["corrected"] is True
        assert report["source"] == "computed"

    def test_an_uncorrected_run_says_what_its_reference_scans_show(self):
        reader = _open(None)
        report = reader.lock_mass_report()
        reader.close()

        # Measured for D83: P17- the strongest ion, about 11 ppm of drift.
        scans = report["reference_scans"]
        assert report["corrected"] is False
        assert scans["looks_like"] == "red phosphorus"
        assert scans["suggested_reference_mz"] == pytest.approx(P17)
        assert abs(scans["reference_peak_drift_ppm"]) > 5.0


class TestRepeatedPixelsOnARealRun:
    def test_every_left_out_scan_was_cut_short_by_a_lock_mass_break(self):
        reader = _open(None)
        report = reader.repeated_pixels_report()
        _ml, _handle, grid, _types, _functions = reader._require_initialized()
        reference = reader._reference_function()
        breaks = sorted(i.rt for (f, _), i in grid.scan_map.items() if f == reference)
        left_out = report["interrupted_scans_left_out"]
        reader.close()

        # Nine on MSV000092638, one per lock-mass break inside the raster.
        assert left_out
        for entry in left_out:
            info = grid.scan_map[(entry["function"], entry["scan"])]
            again = grid.scan_map[(entry["function"], entry["scan"] + 1)]
            coords = grid.get_coordinates(again)
            assert coords is not None and coords[:2] == tuple(entry["pixel"])
            assert any(info.rt < t < again.rt for t in breaks)
