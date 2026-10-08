"""The imzML reader yields float64 on every path, whatever the file stores (D36).

Many imzML files store m/z and intensity as 32-bit floats. pyimzml decodes
them as stored, and everything after the reader -- the TIC sum, the
intensity threshold, typed m/z bounds, ``tic_preserving`` totals --
computes in the array's own type. The other readers yield float64, so a
32-bit imzML was the one source whose sums and comparisons ran in 32-bit.
"""

import importlib.util
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


MOB32_MZ = np.array([300.0, 300.0, 450.5], np.float32)
MOB32_K0 = np.array([0.95, 1.10, 1.02], np.float32)
MOB32_INTENSITIES = [
    np.array([10.0, 1.0, 5.0], np.float32),
    np.array([11.0, 2.0, 6.0], np.float32),
]


def _load_fixture_builder():
    spec = importlib.util.spec_from_file_location(
        "build_fixtures_f32", FIXTURES / "build_fixtures.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_32bit_mobility_imzml(directory: Path, monkeypatch) -> Path:
    """A two-pixel processed mobility export whose three arrays are 32-bit.

    pyimzml's writer cannot emit the third array, so the document is
    assembled from the committed fixture builder's pieces, with the
    precision declarations swapped and the .ibd packed as float32.
    """
    builder = _load_fixture_builder()
    monkeypatch.setattr(builder, "FIXTURE_DIR", directory)

    uuid_text = "8F7E3B1A-2C4D-4E5F-9A6B-1C3E5D7B9A01"
    blob = bytearray(builder.uuid_module.UUID(uuid_text).bytes)
    placements = []
    for row in MOB32_INTENSITIES:
        mz_placement = (len(blob), int(MOB32_MZ.size))
        blob += MOB32_MZ.tobytes()
        intensity_placement = (len(blob), int(row.size))
        blob += row.tobytes()
        mobility_placement = (len(blob), int(MOB32_K0.size))
        blob += MOB32_K0.tobytes()
        placements.append((mz_placement, intensity_placement, mobility_placement))
    (directory / "mob32.ibd").write_bytes(bytes(blob))

    float32_line = (
        '      <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"'
        ' value="" />'
    )
    mobility_group = [
        '    <referenceableParamGroup id="mobilityArray">',
        '      <cvParam cvRef="MS" accession="MS:1000576" name="no compression"'
        ' value="" />',
        '      <cvParam cvRef="MS" accession="MS:1003006"'
        ' name="mean inverse reduced ion mobility array" unitCvRef="MS"'
        ' unitAccession="MS:1002814" unitName="volt-second per square centimeter" />',
        float32_line,
        '      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data"'
        ' value="true" />',
        "    </referenceableParamGroup>",
    ]
    lines = builder._canonical_header(
        uuid_text,
        builder.PLAIN_SCAN_SETTINGS,
        [float32_line],
        len(placements),
        spectrum_type_line=builder.CENTROID_LINE,
        instrument_lines=builder.TIMSTOF_INSTRUMENT_LINES,
        extra_param_group_lines=mobility_group,
        file_mode_line=builder.PROCESSED_MODE_LINE,
    )
    # The canonical header declares the intensity array 64-bit.
    swapped = False
    for i, line in enumerate(lines):
        if 'accession="MS:1000523" name="64-bit float"' in line:
            lines[i] = float32_line
            swapped = True
            break
    assert swapped, "canonical header no longer declares a 64-bit intensity"
    for i, placement in enumerate(placements):
        lines += [
            "      " + line
            for line in builder._spectrum_xml_with_mobility(
                i, builder.DENSE_COORDINATES[i], *placement, itemsize=4
            )
        ]
    lines += ["    </spectrumList>", "  </run>", "</mzML>"]
    return builder._write("mob32", lines, newline=b"\n", encoding="utf-8")


class TestMobilityPointCloud:
    def test_a_32_bit_mobility_export_yields_float64(self, tmp_path, monkeypatch):
        path = _write_32bit_mobility_imzml(tmp_path, monkeypatch)

        with ImzMLReader(path) as reader:
            clouds = list(reader.iter_mobility_spectra())

        assert len(clouds) == 2
        for (_, mzs, mobility, intensities), written in zip(clouds, MOB32_INTENSITIES):
            assert mzs.dtype == np.float64
            assert mobility.dtype == np.float64
            assert intensities.dtype == np.float64
            np.testing.assert_array_equal(mzs, MOB32_MZ.astype(np.float64))
            np.testing.assert_array_equal(mobility, MOB32_K0.astype(np.float64))
            np.testing.assert_array_equal(intensities, written.astype(np.float64))

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
