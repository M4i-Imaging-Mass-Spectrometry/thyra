"""Lock-mass correction of Waters m/z without writing to the run (D82)."""

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from thyra.errors import ConversionRefused
from thyra.readers.waters.lockmass import (
    LockMassCorrection,
    find_lmgt,
    measure_lock_mass,
    read_lmgt,
)

#: The shape the Waters library writes: a version line, reference and
#: tolerance, then one retention time and factor per lock-mass scan.
_LMGT = "1\n526.554504395,0.010000000\n0.1,0.99999\n20.3,0.99998\n40.5,0.99997\n"


def _correction(times, factors, source="lmgt.inf"):
    return LockMassCorrection(
        reference_mz=526.5545,
        tolerance_da=0.01,
        times_min=np.asarray(times, dtype=float),
        factors=np.asarray(factors, dtype=float),
        source=source,
    )


class TestReadingTheCorrectionFile:
    def test_the_file_is_read(self, tmp_path):
        path = tmp_path / "lmgt.inf"
        path.write_text(_LMGT)

        correction = read_lmgt(path)

        assert correction.reference_mz == pytest.approx(526.554504395)
        assert correction.tolerance_da == pytest.approx(0.01)
        assert correction.times_min.tolist() == [0.1, 20.3, 40.5]
        assert correction.factors.tolist() == [0.99999, 0.99998, 0.99997]
        assert correction.source == "lmgt.inf"

    def test_a_file_without_measurements_is_refused(self, tmp_path):
        path = tmp_path / "lmgt.inf"
        path.write_text("1\n526.5545,0.01\n")

        with pytest.raises(ConversionRefused, match="no lock-mass measurements"):
            read_lmgt(path)

    def test_a_damaged_file_is_refused_rather_than_ignored(self, tmp_path):
        path = tmp_path / "lmgt.inf"
        path.write_text("1\n526.5545,0.01\n0.1,not-a-number\n")

        with pytest.raises(ConversionRefused, match="not a lock-mass correction"):
            read_lmgt(path)

    def test_the_file_is_found_whatever_its_case(self, tmp_path):
        (tmp_path / "LMGT.INF").write_text(_LMGT)

        assert find_lmgt(tmp_path) == tmp_path / "LMGT.INF"

    def test_no_file_is_none(self, tmp_path):
        assert find_lmgt(tmp_path) is None


class TestTheFactorOfAScan:
    """Linear between measurements, held at both ends, as the library applies it."""

    def test_between_two_measurements_it_is_interpolated(self):
        correction = _correction([10.0, 20.0], [0.99990, 0.99980])

        assert correction.factor_at(15.0) == pytest.approx(0.99985)

    def test_before_the_first_and_after_the_last_it_is_held(self):
        correction = _correction([10.0, 20.0], [0.99990, 0.99980])

        assert correction.factor_at(1.0) == pytest.approx(0.99990)
        assert correction.factor_at(99.0) == pytest.approx(0.99980)

    def test_the_report_gives_the_range_in_ppm(self):
        report = _correction([10.0, 20.0], [0.99990, 0.99980]).report()

        assert report["n_measurements"] == 2
        assert report["correction_ppm_min"] == pytest.approx(-200.0)
        assert report["correction_ppm_max"] == pytest.approx(-100.0)
        assert report["reference_mz"] == pytest.approx(526.5545)


def _run_folder(tmp_path: Path) -> Path:
    raw = tmp_path / "run.raw"
    raw.mkdir()
    (raw / "_HEADER.TXT").write_text("$$ Instrument: test\n")
    (raw / "_FUNCTNS.INF").write_bytes(b"\x01" * 832)
    (raw / "_FUNC001.DAT").write_bytes(b"\x07" * 4096)
    (raw / "_FUNC001.IDX").write_bytes(b"\x02" * 60)
    (raw / "_func001.cdt").write_bytes(b"\x08" * 2048)
    (raw / "_FUNC002.DAT").write_bytes(bytes(range(256)) * 4)
    (raw / "_FUNC002.IDX").write_bytes(b"\x03" * 30)
    (raw / "_func002.cdt").write_bytes(b"\x09" * 512)
    (raw / "lmgt.inf").write_text("1\n774.3446,0.01\n0.1,0.9\n")
    (raw / "imaging").mkdir()
    return raw


def _snapshot(folder: Path):
    return {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file()}


class TestMeasuringOnACopy:
    """The library measures on a throwaway copy; the run is never written to."""

    def _library(self, seen, applied=True):
        ml = MagicMock()
        ml.open_file.side_effect = lambda path: seen.setdefault("path", Path(path))

        def apply(_handle, reference_mz):
            copy = seen["path"]
            seen["files"] = {p.name: p.read_bytes() for p in copy.iterdir()}
            seen["reference_mz"] = reference_mz
            if applied:
                (copy / "lmgt.inf").write_text(_LMGT)
            return applied

        ml.apply_custom_lockmass_correction.side_effect = apply
        return ml

    def test_the_factors_are_the_ones_the_library_wrote(self, tmp_path):
        raw = _run_folder(tmp_path)
        seen: dict = {}

        correction = measure_lock_mass(self._library(seen), raw, 1, 526.5545)

        assert correction.source == "computed"
        assert correction.factors.tolist() == [0.99999, 0.99998, 0.99997]
        assert seen["reference_mz"] == 526.5545

    def test_the_run_folder_is_left_exactly_as_it_was(self, tmp_path):
        raw = _run_folder(tmp_path)
        before = _snapshot(raw)

        measure_lock_mass(self._library({}), raw, 1, 526.5545)

        assert _snapshot(raw) == before

    def test_the_copy_holds_the_reference_scans_and_placeholders(self, tmp_path):
        raw = _run_folder(tmp_path)
        seen: dict = {}

        measure_lock_mass(self._library(seen), raw, 1, 526.5545)

        files = seen["files"]
        assert seen["path"] != raw
        # The lock-mass function's data, byte for byte.
        assert files["_FUNC002.DAT"] == (raw / "_FUNC002.DAT").read_bytes()
        # Every side file, whole.
        assert files["_HEADER.TXT"] == (raw / "_HEADER.TXT").read_bytes()
        assert files["_FUNC001.IDX"] == (raw / "_FUNC001.IDX").read_bytes()
        # The other data files: the right size, and nothing of the run.
        assert files["_FUNC001.DAT"] == b"\x00" * 4096
        assert files["_func001.cdt"] == b"\x00" * 2048
        assert files["_func002.cdt"] == b"\x00" * 512
        # A previous correction is not measured on top of.
        assert "lmgt.inf" not in files

    def test_the_copy_is_removed_afterwards(self, tmp_path):
        raw = _run_folder(tmp_path)
        seen: dict = {}

        measure_lock_mass(self._library(seen), raw, 1, 526.5545)

        assert not seen["path"].exists()

    def test_no_peak_found_is_refused_by_name(self, tmp_path):
        raw = _run_folder(tmp_path)

        with pytest.raises(ConversionRefused, match="no lock-mass peak at m/z 526.554"):
            measure_lock_mass(self._library({}, applied=False), raw, 1, 526.5545)
