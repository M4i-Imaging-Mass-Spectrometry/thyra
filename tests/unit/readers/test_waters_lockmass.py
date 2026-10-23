"""Lock-mass correction of Waters m/z without writing to the run (D84)."""

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
    survey_reference_scans,
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


# ----------------------------------------------------------------------
# What an uncorrected run's reference scans show (D84)
# ----------------------------------------------------------------------

_P = 30.97376163
_E = 0.00054858
#: A shared profile grid, 1 mDa apart, as an uncorrected run has.
_GRID = np.round(np.arange(100.0, 1200.0, 0.001), 5)


def _profile(peaks, offset_ppm):
    """Gaussian peaks (m/z, height) on the grid, every m/z shifted by offset."""
    intensity = np.zeros_like(_GRID)
    for mz, height in peaks:
        centre = mz * (1 + offset_ppm * 1e-6)
        intensity += height * np.exp(-0.5 * ((_GRID - centre) / 0.004) ** 2)
    return _GRID.copy(), intensity


def _red_phosphorus(negative=True, n_values=range(5, 37, 2)):
    """The cluster ladder, with P17 the strongest ion as on the D84 run."""
    sign = 1.0 if negative else -1.0
    return [(n * _P + sign * _E, 20000.0 if n == 17 else 2000.0) for n in n_values]


def _groups(peaks, offsets_ppm, scans_per_group=3, minutes_apart=20.0):
    scans = []
    for g, offset in enumerate(offsets_ppm):
        for s in range(scans_per_group):
            mz, intensity = _profile(peaks, offset)
            scans.append((g * minutes_apart + 0.02 * s, mz, intensity))
    return scans


class TestSurveyingReferenceScans:
    def test_the_drift_of_the_strongest_peak_needs_no_identification(self):
        survey = survey_reference_scans(
            _groups([(600.0, 5000.0), (300.0, 100.0)], [10.0, 15.0, 22.0]),
            negative=True,
        )

        assert survey is not None
        assert survey.n_groups == 3
        assert survey.drift_ppm == pytest.approx(12.0, abs=1.0)
        assert survey.span_ppm == pytest.approx(12.0, abs=1.0)
        assert survey.peak_mz == pytest.approx(600.0 * (1 + 15.7e-6), abs=0.005)

    def test_a_red_phosphorus_ladder_is_named_with_its_strongest_ion(self):
        survey = survey_reference_scans(
            _groups(_red_phosphorus(), [12.0, 16.0, 20.0]), negative=True
        )

        assert survey.looks_like == "red phosphorus"
        assert survey.matched_ions >= 10
        assert survey.suggested_mz == pytest.approx(17 * _P + _E)
        assert survey.offset_ppm == pytest.approx(16.0, abs=1.0)

    def test_positive_mode_uses_the_cations(self):
        survey = survey_reference_scans(
            _groups(_red_phosphorus(negative=False), [12.0, 16.0]), negative=False
        )

        assert survey.looks_like == "red phosphorus"
        assert survey.suggested_mz == pytest.approx(17 * _P - _E)

    def test_leucine_enkephalin_is_named_as_a_single_ion(self):
        survey = survey_reference_scans(
            _groups([(556.2766, 8000.0), (278.1, 900.0)], [10.0, 14.0]),
            negative=False,
        )

        assert survey.looks_like == "leucine enkephalin"
        assert survey.matched_ions == 1
        assert survey.suggested_mz == pytest.approx(556.2766)

    def test_an_unknown_reference_is_named_nothing(self):
        survey = survey_reference_scans(
            _groups([(622.0290, 6000.0), (922.0098, 3000.0)], [10.0, 14.0]),
            negative=False,
        )

        assert survey.looks_like is None
        assert survey.suggested_mz is None
        # The drift is still reported: it needs no identification.
        assert survey.drift_ppm == pytest.approx(4.0, abs=1.0)

    def test_a_ladder_whose_ions_disagree_is_not_named(self):
        # Every other ion 20 ppm further out: not one compound drifting.
        peaks = [
            (mz * (1 + (20e-6 if i % 2 else 0.0)), h)
            for i, (mz, h) in enumerate(_red_phosphorus())
        ]
        survey = survey_reference_scans(_groups(peaks, [12.0, 16.0]), negative=True)

        assert survey.looks_like is None

    def test_scans_without_signal_give_no_survey(self):
        empty = (1.0, np.array([]), np.array([]))
        assert survey_reference_scans([empty], negative=True) is None
