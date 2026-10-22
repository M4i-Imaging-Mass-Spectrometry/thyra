# tests/unit/readers/test_imzml_continuous_fast_read.py
"""The continuous-mode fast read must be invisible next to pyimzml's reads.

``_read_spectrum_arrays`` decodes the shared m/z block once and reads only
intensity bytes per spectrum -- licensed by ``_continuous_uniform``, which
is measured off the offset tables rather than trusted from the declared
mode. Nothing it yields may differ from ``getspectrum``: everything a
conversion stores flows through these arrays, twice (the streaming routes
iterate the reader in two passes).
"""

from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLParser import ImzMLParser
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.readers.imzml.imzml_reader import ImzMLReader


def _write_processed(path: Path, n_spectra: int = 30) -> None:
    with ImzMLWriter(str(path), mode="processed") as writer:
        for i in range(n_spectra):
            rng = np.random.default_rng([11, i])
            n = int(rng.integers(5, 40))
            mzs = np.sort(rng.uniform(100.0, 900.0, n))
            ints = rng.exponential(50.0, n).astype(np.float32)
            writer.addSpectrum(mzs, ints, (i % 6 + 1, i // 6 + 1, 1))


def _write_continuous(path: Path, n_spectra: int = 20, n_points: int = 300) -> None:
    axis = np.linspace(100.0, 900.0, n_points)
    with ImzMLWriter(str(path), mode="continuous") as writer:
        for i in range(n_spectra):
            rng = np.random.default_rng([13, i])
            ints = rng.exponential(50.0, n_points).astype(np.float32)
            ints[rng.random(n_points) < 0.3] = 0.0
            writer.addSpectrum(axis, ints, (i % 5 + 1, i // 5 + 1, 1))


def _reference_spectra(path: Path):
    """(mzs, intensities) per spectrum straight from pyimzml, no Thyra."""
    parser = ImzMLParser(str(path), parse_lib="ElementTree")
    try:
        return [parser.getspectrum(i) for i in range(len(parser.coordinates))]
    finally:
        parser.m.close()


@pytest.mark.parametrize("mode", ["processed", "continuous"])
def test_iter_spectra_matches_getspectrum(tmp_path, mode):
    path = tmp_path / f"{mode}.imzML"
    if mode == "processed":
        _write_processed(path)
    else:
        _write_continuous(path)
    reference = _reference_spectra(path)

    reader = ImzMLReader(path)
    try:
        got = list(reader.iter_spectra())
        if mode == "continuous":
            # The fast read must actually be active, or this test would
            # silently cover the fallback instead.
            assert reader._continuous_uniform
        assert len(got) == len(reference)
        for (_, mzs, ints), (ref_mzs, ref_ints) in zip(got, reference):
            np.testing.assert_array_equal(np.asarray(mzs), ref_mzs)
            np.testing.assert_array_equal(np.asarray(ints), ref_ints)
            # The same values, widened to float64 (D36).
            assert np.asarray(mzs).dtype == np.float64
            assert np.asarray(ints).dtype == np.float64
    finally:
        reader.close()


@pytest.mark.parametrize("axis_first", [True, False])
def test_continuous_yields_one_shared_mz_object(tmp_path, axis_first):
    """Every yielded m/z IS the common-axis object, whichever is asked first.

    Identity (not mere equality) is what makes the converters' shared-axis
    equality checks O(1) per spectrum.
    """
    path = tmp_path / "cont.imzML"
    _write_continuous(path)

    reader = ImzMLReader(path)
    try:
        if axis_first:
            axis = reader.get_common_mass_axis()
            spectra = list(reader.iter_spectra())
        else:
            spectra = list(reader.iter_spectra())
            axis = reader.get_common_mass_axis()
        assert all(mzs is axis for _, mzs, _ in spectra)
    finally:
        reader.close()


def _write_mislabelled(path: Path) -> None:
    """One m/z block per spectrum, under a continuous declaration (F121, F156)."""
    with ImzMLWriter(str(path), mode="processed") as writer:
        writer.addSpectrum(
            np.array([100.0, 200.0, 300.0]), np.array([1.0, 2.0, 3.0]), (1, 1, 1)
        )
        writer.addSpectrum(
            np.array([50.0, 200.0, 2000.0]), np.array([4.0, 5.0, 6.0]), (2, 1, 1)
        )
        writer.addSpectrum(
            np.array([100.0, 200.0, 300.0]), np.array([7.0, 8.0, 9.0]), (3, 1, 1)
        )
    text = path.read_text(encoding="latin-1")
    path.write_text(
        text.replace(
            '"IMS:1000031" name="processed"', '"IMS:1000030" name="continuous"'
        ),
        encoding="latin-1",
    )


class TestAContinuousDeclarationOverOwnBlocks:
    """Measured from the offsets, the file is processed (D82).

    Spectrum 0's block gave the whole file its mass range and its
    ``--no-resample`` axis, so m/z 50 and 2000 of spectrum 1 were dropped.
    Replacing the block test with ``True`` passed every test before these.
    """

    def test_the_reader_reads_every_spectrum_at_its_own_m_z(self, tmp_path):
        path = tmp_path / "mislabelled.imzML"
        _write_mislabelled(path)
        with ImzMLReader(path) as reader:
            assert reader.has_shared_mass_axis is False
            assert reader.is_continuous is True
            assert not reader._continuous_uniform
            spectra = list(reader.iter_spectra())
            axis = reader.get_common_mass_axis()
        np.testing.assert_array_equal(spectra[1][1], [50.0, 200.0, 2000.0])
        np.testing.assert_array_equal(axis, [50.0, 100.0, 200.0, 300.0, 2000.0])

    def test_the_mass_range_covers_every_spectrum(self, tmp_path):
        path = tmp_path / "mislabelled.imzML"
        _write_mislabelled(path)
        with ImzMLReader(path) as reader:
            assert reader.get_essential_metadata().mass_range == (50.0, 2000.0)

    def test_a_true_continuous_file_still_reads_one_spectrum(self, tmp_path):
        path = tmp_path / "cont.imzML"
        _write_continuous(path)
        with ImzMLReader(path) as reader:
            assert reader.has_shared_mass_axis is True
            assert reader.get_essential_metadata().mass_range == (100.0, 900.0)


def test_two_passes_agree(tmp_path):
    """The streaming converters iterate twice; both passes must match."""
    path = tmp_path / "cont.imzML"
    _write_continuous(path)

    reader = ImzMLReader(path)
    try:
        first = [(c, m.copy(), i.copy()) for c, m, i in reader.iter_spectra()]
        second = list(reader.iter_spectra())
        assert len(first) == len(second)
        for (c1, m1, i1), (c2, m2, i2) in zip(first, second):
            assert c1 == c2
            np.testing.assert_array_equal(m1, np.asarray(m2))
            np.testing.assert_array_equal(i1, np.asarray(i2))
    finally:
        reader.close()
