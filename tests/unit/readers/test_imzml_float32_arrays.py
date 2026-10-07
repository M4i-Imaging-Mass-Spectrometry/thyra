# tests/unit/readers/test_imzml_float32_arrays.py
"""32-bit imzML arrays arrive widened to float64 on every read path.

Many imzML files store m/z and intensity as 32-bit floats. The reader used
to hand those arrays on in their storage dtype, and every downstream sum and
comparison then ran in 32-bit: the TIC image missed its own row sums at the
1e-7 level, ``tic_preserving`` rescaled rows to a rounded pixel total, a
threshold of 0.7 kept the value float32 holds as 0.7, and a typed m/z bound
was applied at float32 precision. The reader now widens to float64 once, at
the point it yields each array, so every consumer sees the same values the
other readers already yielded.
"""

import importlib.util
from pathlib import Path

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

from thyra.readers.imzml.imzml_reader import ImzMLReader

FIXTURE_BUILDER = (
    Path(__file__).resolve().parents[2] / "data" / "fixtures" / "build_fixtures.py"
)

MOB32_MZ = np.array([300.0, 300.0, 450.5], np.float32)
MOB32_K0 = np.array([0.95, 1.10, 1.02], np.float32)
MOB32_INTENSITIES = [
    np.array([10.0, 1.0, 5.0], np.float32),
    np.array([11.0, 2.0, 6.0], np.float32),
]


def _write_processed_f32(path: Path) -> tuple:
    """A processed file whose m/z and intensity blocks are both 32-bit."""
    spectra = [
        (np.array([150.5, 250.25, 350.0]), np.array([2.0**24, 1.0, 1.0])),
        (np.array([150.5, 250.25, 350.0]), np.array([1.0, 4.0, 0.5])),
    ]
    with ImzMLWriter(
        str(path), mode="processed", mz_dtype=np.float32, intensity_dtype=np.float32
    ) as writer:
        for i, (mzs, intensities) in enumerate(spectra):
            writer.addSpectrum(mzs, intensities.astype(np.float32), (i + 1, 1, 1))
    return spectra


def _write_continuous_f32(path: Path) -> tuple:
    """A continuous file whose shared m/z block and intensities are 32-bit."""
    mzs = np.linspace(100.0, 500.0, 64)
    spectra = [np.full(64, 3.0), np.arange(64, dtype=np.float64) + 0.5]
    with ImzMLWriter(
        str(path), mode="continuous", mz_dtype=np.float32, intensity_dtype=np.float32
    ) as writer:
        for i, intensities in enumerate(spectra):
            writer.addSpectrum(
                mzs.astype(np.float32),
                intensities.astype(np.float32),
                (i + 1, 1, 1),
            )
    return mzs, spectra


class TestProcessedFile:
    def test_iter_spectra_yields_float64(self, tmp_path):
        spectra = _write_processed_f32(tmp_path / "p32.imzML")

        with ImzMLReader(tmp_path / "p32.imzML") as reader:
            items = list(reader.iter_spectra())

        assert len(items) == 2
        for (_, mzs, intensities), (want_mzs, want_intensities) in zip(items, spectra):
            assert mzs.dtype == np.float64
            assert intensities.dtype == np.float64
            # float32 -> float64 is exact; the values must be the stored ones
            np.testing.assert_array_equal(mzs, want_mzs)
            np.testing.assert_array_equal(intensities, want_intensities)

    def test_common_mass_axis_is_float64(self, tmp_path):
        spectra = _write_processed_f32(tmp_path / "p32.imzML")

        with ImzMLReader(tmp_path / "p32.imzML") as reader:
            axis = reader.get_common_mass_axis()

        assert axis.dtype == np.float64
        np.testing.assert_array_equal(axis, spectra[0][0])


class TestContinuousFile:
    def test_iter_spectra_yields_float64(self, tmp_path):
        mzs, spectra = _write_continuous_f32(tmp_path / "c32.imzML")
        stored_mzs = mzs.astype(np.float32).astype(np.float64)

        with ImzMLReader(tmp_path / "c32.imzML") as reader:
            items = list(reader.iter_spectra())

        assert len(items) == 2
        for (_, got_mzs, intensities), want in zip(items, spectra):
            assert got_mzs.dtype == np.float64
            assert intensities.dtype == np.float64
            np.testing.assert_array_equal(got_mzs, stored_mzs)
            np.testing.assert_array_equal(intensities, want)

    def test_shared_axis_stays_identity_equal_and_float64(self, tmp_path):
        """The converter's O(1) axis check relies on the yielded m/z object
        being the very object get_common_mass_axis cached."""
        mzs, _ = _write_continuous_f32(tmp_path / "c32.imzML")

        with ImzMLReader(tmp_path / "c32.imzML") as reader:
            axis = reader.get_common_mass_axis()
            yielded_mzs = [mzs for _, mzs, _ in reader.iter_spectra()]

        assert axis.dtype == np.float64
        for m in yielded_mzs:
            assert m is axis


class TestMobilityPointCloud:
    def test_iter_mobility_spectra_yields_float64(self, tmp_path, monkeypatch):
        """A 32-bit mobility export: all three arrays come out float64."""
        path = _write_32bit_mobility_imzml(tmp_path, monkeypatch)

        with ImzMLReader(path) as reader:
            items = list(reader.iter_mobility_spectra())

        assert len(items) == 2
        for _, mzs, mobility, intensities in items:
            assert mzs.dtype == np.float64
            assert mobility.dtype == np.float64
            assert intensities.dtype == np.float64
        np.testing.assert_array_equal(items[0][1], MOB32_MZ.astype(np.float64))
        np.testing.assert_array_equal(items[0][2], MOB32_K0.astype(np.float64))
        np.testing.assert_array_equal(items[0][3], MOB32_INTENSITIES[0])


def _load_fixture_builder():
    spec = importlib.util.spec_from_file_location("build_fixtures_f32", FIXTURE_BUILDER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_32bit_mobility_imzml(directory: Path, monkeypatch) -> Path:
    """A two-pixel processed mobility export whose three arrays are 32-bit.

    Reuses the committed fixture builder's header so the document is valid
    by construction; swaps the precision declarations and packs the .ibd
    as float32. pyimzml's writer cannot emit the third array, so the XML
    is hand-assembled here the way build_fixtures.py does.
    """
    builder = _load_fixture_builder()
    monkeypatch.setattr(builder, "FIXTURE_DIR", directory)

    mz32 = MOB32_MZ
    k032 = MOB32_K0
    intensities = MOB32_INTENSITIES

    uuid_text = "8F7E3B1A-2C4D-4E5F-9A6B-1C3E5D7B9A01"
    blob = bytearray(builder.uuid_module.UUID(uuid_text).bytes)
    placements = []
    for row in intensities:
        mz_placement = (len(blob), int(mz32.size))
        blob += mz32.tobytes()
        intensity_placement = (len(blob), int(row.size))
        blob += row.tobytes()
        mobility_placement = (len(blob), int(k032.size))
        blob += k032.tobytes()
        placements.append((mz_placement, intensity_placement, mobility_placement))
    (directory / "mob32.ibd").write_bytes(bytes(blob))

    mz_precision = [
        '      <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"'
        ' value="" />',
    ]
    mobility_group = [
        '    <referenceableParamGroup id="mobilityArray">',
        '      <cvParam cvRef="MS" accession="MS:1000576" name="no compression"'
        ' value="" />',
        '      <cvParam cvRef="MS" accession="MS:1003006"'
        ' name="mean inverse reduced ion mobility array" unitCvRef="MS"'
        ' unitAccession="MS:1002814" unitName="volt-second per square centimeter" />',
        '      <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"'
        ' value="" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data"'
        ' value="true" />',
        "    </referenceableParamGroup>",
    ]
    lines = builder._canonical_header(
        uuid_text,
        builder.PLAIN_SCAN_SETTINGS,
        mz_precision,
        len(placements),
        spectrum_type_line=builder.CENTROID_LINE,
        instrument_lines=builder.TIMSTOF_INSTRUMENT_LINES,
        extra_param_group_lines=mobility_group,
        file_mode_line=builder.PROCESSED_MODE_LINE,
    )
    # The canonical header declares the intensity array 64-bit; this file
    # stores it 32-bit, so swap that one declaration.
    swapped = False
    for i, line in enumerate(lines):
        if 'accession="MS:1000523" name="64-bit float"' in line:
            lines[i] = (
                '      <cvParam cvRef="MS" accession="MS:1000521"'
                ' name="32-bit float" value="" />'
            )
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


class TestIntensityThresholdAtFloat32Values:
    def test_a_threshold_between_float32_neighbours_drops_the_lower(self, tmp_path):
        """0.7 keeps only values that are >= 0.7 as numbers, not as float32.

        float32(0.7) is 0.699999988079071, below 0.7: it must be dropped.
        Compared in the array's own dtype, it passed.
        """
        mzs = np.array([150.0, 250.0, 350.0])
        intensities = np.array([0.7, 4.0, 2.0], np.float32)
        with ImzMLWriter(
            str(tmp_path / "t32.imzML"),
            mode="processed",
            intensity_dtype=np.float32,
        ) as writer:
            writer.addSpectrum(mzs, intensities, (1, 1, 1))

        with ImzMLReader(tmp_path / "t32.imzML", intensity_threshold=0.7) as reader:
            kept = [it for _, _, it in reader.iter_spectra()]

        assert len(kept) == 1
        np.testing.assert_array_equal(kept[0], [4.0, 2.0])
