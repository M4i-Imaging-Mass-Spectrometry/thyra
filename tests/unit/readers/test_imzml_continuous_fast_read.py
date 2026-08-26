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
            assert np.asarray(ints).dtype == ref_ints.dtype
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
