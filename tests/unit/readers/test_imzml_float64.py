"""The imzML reader yields float64 on every path, whatever the file stores (D36).

Many imzML files store m/z and intensity as 32-bit floats. pyimzml decodes
them as stored, and everything after the reader -- the TIC sum, the
intensity threshold, typed m/z bounds, ``tic_preserving`` totals --
computes in the array's own type. The other readers yield float64, so a
32-bit imzML was the one source whose sums and comparisons ran in 32-bit.
"""

from pathlib import Path

import numpy as np
import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.readers.imzml.imzml_reader import ImzMLReader

FIXTURES = Path(__file__).resolve().parents[2] / "data" / "fixtures"

# Values 32-bit floats round: written as float32, read back as their float32
# copies widened to float64.
MZS = np.array([150.1, 250.3, 350.7])
INTENSITIES = [np.array([0.7, 4.1, 2.3]), np.array([1.3, 0.9, 5.5])]


def _write(directory: Path, mode: str) -> Path:
    path = directory / f"f32_{mode}.imzML"
    with ImzMLWriter(
        str(path),
        mode=mode,
        mz_dtype=np.float32,
        intensity_dtype=np.float32,
    ) as writer:
        for x, intensities in enumerate(INTENSITIES, start=1):
            writer.addSpectrum(MZS, intensities, (x, 1, 1))
    return path


def _assert_float64_copy(actual: np.ndarray, written: np.ndarray) -> None:
    assert actual.dtype == np.float64
    np.testing.assert_array_equal(actual, written.astype(np.float32).astype(np.float64))


@pytest.mark.parametrize("mode", ["processed", "continuous"])
class TestSpectraAreFloat64:
    def test_iter_spectra(self, tmp_path, mode):
        with ImzMLReader(_write(tmp_path, mode)) as reader:
            spectra = list(reader.iter_spectra())
        assert len(spectra) == 2
        for (_, mzs, intensities), written in zip(spectra, INTENSITIES):
            _assert_float64_copy(mzs, MZS)
            _assert_float64_copy(intensities, written)

    def test_common_mass_axis_and_read(self, tmp_path, mode):
        with ImzMLReader(_write(tmp_path, mode)) as reader:
            _assert_float64_copy(reader.get_common_mass_axis(), MZS)
            data = reader.read()
        _assert_float64_copy(data["mzs"], MZS)


class TestContinuousSharedAxis:
    """The shared m/z block is widened once, so it stays one object."""

    @pytest.mark.parametrize("axis_first", [True, False])
    def test_one_float64_object_for_every_spectrum(self, tmp_path, axis_first):
        with ImzMLReader(_write(tmp_path, "continuous")) as reader:
            axis = reader.get_common_mass_axis() if axis_first else None
            yielded = [mzs for _, mzs, _ in reader.iter_spectra()]
            if axis is None:
                axis = reader.get_common_mass_axis()
        assert axis.dtype == np.float64
        assert all(mzs is axis for mzs in yielded)


class TestIntensityThreshold:
    def test_a_32_bit_value_below_the_threshold_is_dropped(self, tmp_path):
        # float32(0.7) is 0.699999988, below 0.7. Compared in 32-bit it
        # equalled the threshold and was kept.
        path = _write(tmp_path, "processed")
        with ImzMLReader(path, intensity_threshold=0.7) as reader:
            kept = [intensities for _, _, intensities in reader.iter_spectra()]
        assert all(np.all(values >= 0.7) for values in kept)
        np.testing.assert_array_equal(kept[0], np.array([4.1, 2.3], dtype=np.float32))


class TestMobilityPointCloud:
    def test_mobility_spectra_are_float64(self, monkeypatch):
        # The committed mobility fixture stores 64-bit floats; hand its
        # values back as 32-bit, as pyimzml does for a 32-bit file.
        with ImzMLReader(FIXTURES / "mobility_processed.imzML") as reader:
            reader.get_essential_metadata()
            decode = reader.parser.getspectrum
            monkeypatch.setattr(
                reader.parser,
                "getspectrum",
                lambda idx: tuple(a.astype(np.float32) for a in decode(idx)),
            )
            clouds = list(reader.iter_mobility_spectra())
        assert clouds
        for _, mzs, _, intensities in clouds:
            assert mzs.dtype == np.float64
            assert intensities.dtype == np.float64

    def test_shared_features_are_the_yielded_float64_block(self, monkeypatch):
        # Asked for first, the shared features decode and cache the m/z
        # block; every spectrum must then yield that same float64 object.
        with ImzMLReader(FIXTURES / "mobility_continuous.imzML") as reader:
            reader.get_essential_metadata()
            decode = reader.parser.getspectrum
            monkeypatch.setattr(
                reader.parser,
                "getspectrum",
                lambda idx: tuple(a.astype(np.float32) for a in decode(idx)),
            )
            features = reader.get_shared_mobility_features()
            assert features is not None
            yielded = [mzs for _, mzs, _, _ in reader.iter_mobility_spectra()]
        assert features[0].dtype == np.float64
        assert yielded and all(mzs is features[0] for mzs in yielded)
