"""imzML with a third binary array: the ion mobility dimension.

Two hand-authored fixtures reproduce the convention TIMSCONVERT and
TIMSImaging share (``build_fixtures.py``): ``mobility_continuous`` is one
shared (m/z, 1/K0) feature list with m/z repeated where mobility splits a
feature; ``mobility_processed`` is a per-pixel point cloud. Upstream pyimzml
ignores the array, so before this the dimension was dropped silently and the
repeated m/z merged.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from tests.fixtures.imzml_mobility_builder import write_mobility_imzml
from thyra.core.mobility import MobilityAxis, classify_mobility_array
from thyra.errors import ConversionRefused
from thyra.readers.imzml.imzml_reader import ImzMLReader
from thyra.readers.imzml.mobility_array import (
    collect_array_offsets,
    detect_mobility_array,
)

FIXTURES = Path(__file__).resolve().parents[2] / "data" / "fixtures"
CONTINUOUS = FIXTURES / "mobility_continuous.imzML"
PROCESSED = FIXTURES / "mobility_processed.imzML"
PLAIN = FIXTURES / "timstof_flex_export.imzML"

FEATURE_MZ = np.array([300.0, 300.0, 450.5, 600.25, 600.25])
FEATURE_K0 = np.array([0.95, 1.10, 1.02, 1.35, 1.20])
INTENSITIES = np.array(
    [
        [10.0, 1.0, 5.0, 2.0, 20.0],
        [11.0, 2.0, 6.0, 3.0, 21.0],
        [12.0, 3.0, 7.0, 4.0, 22.0],
        [13.0, 4.0, 8.0, 0.0, 23.0],
    ]
)
# DENSE_COORDINATES (1,1) (2,1) (1,2) (2,2), 0-based
COORDS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)]


class TestContinuousExport:
    def test_the_array_is_detected_and_described(self):
        with ImzMLReader(CONTINUOUS) as reader:
            assert reader.has_ion_mobility is True
            axis = reader.get_mobility_axis()
        assert isinstance(axis, MobilityAxis)
        assert axis.kind_accession == "MS:1002815"
        assert axis.kind_name == "inverse reduced ion mobility"
        assert axis.unit_accession == "MS:1002814"
        assert axis.array_accession == "MS:1003006"
        assert axis.source == "imzml"
        assert axis.is_shared
        np.testing.assert_array_equal(axis.values, FEATURE_K0)
        assert axis.acq_range is None
        block = axis.to_uns()
        assert (block["range_lower"], block["range_upper"]) == (0.95, 1.35)
        assert block["n_values"] == 5

    def test_the_feature_axis_is_shared_and_keeps_the_repeats(self):
        with ImzMLReader(CONTINUOUS) as reader:
            assert reader.has_shared_mobility_axis is True
            mz, mobility = reader.get_shared_mobility_features()
        np.testing.assert_array_equal(mz, FEATURE_MZ)
        np.testing.assert_array_equal(mobility, FEATURE_K0)

    def test_the_common_mass_axis_collapses_the_repeats(self):
        with ImzMLReader(CONTINUOUS) as reader:
            axis = reader.get_common_mass_axis()
            spectra = list(reader.iter_spectra())
        np.testing.assert_array_equal(axis, [300.0, 450.5, 600.25])
        assert len(spectra) == 4
        # The plain iteration still yields the shared block as written,
        # repeats included; the converter sums them per axis bin.
        for (coords, mzs, intensities), expected_coords, row in zip(
            spectra, COORDS, INTENSITIES
        ):
            assert coords == expected_coords
            np.testing.assert_array_equal(mzs, FEATURE_MZ)
            np.testing.assert_array_equal(intensities, row)

    def test_mobility_spectra_are_parallel_arrays(self):
        with ImzMLReader(CONTINUOUS) as reader:
            items = list(reader.iter_mobility_spectra())
        assert [c for c, *_ in items] == COORDS
        for (_, mzs, mobility, intensities), row in zip(items, INTENSITIES):
            np.testing.assert_array_equal(mzs, FEATURE_MZ)
            np.testing.assert_array_equal(mobility, FEATURE_K0)
            np.testing.assert_array_equal(intensities, row)

    def test_intensity_threshold_masks_all_three_arrays_together(self):
        """A peak's points go together, by their sum over the scans (D34).

        Pixel 0 sums to 11 at m/z 300, 5 at 450.5 and 22 at 600.25. At 5.5
        the 450.5 peak drops out, and the 1 and 2 points of the other two
        stay, though each is below the threshold on its own.
        """
        with ImzMLReader(CONTINUOUS, intensity_threshold=5.5) as reader:
            _, mzs, mobility, intensities = next(reader.iter_mobility_spectra())
        np.testing.assert_array_equal(mzs, [300.0, 300.0, 600.25, 600.25])
        np.testing.assert_array_equal(mobility, [0.95, 1.10, 1.35, 1.20])
        np.testing.assert_array_equal(intensities, [10.0, 1.0, 2.0, 20.0])

    def test_the_reader_reports_its_threshold(self):
        """What the converter records in the conversion step."""
        with ImzMLReader(CONTINUOUS, intensity_threshold=5.5) as reader:
            assert reader.intensity_threshold == 5.5
        with ImzMLReader(CONTINUOUS) as reader:
            assert reader.intensity_threshold is None

    def test_the_summed_spectrum_keeps_the_same_points(self):
        """The summed table tests the same peak sums as the point cloud."""
        with ImzMLReader(CONTINUOUS, intensity_threshold=2.5) as reader:
            _, mzs, intensities = next(reader.iter_spectra())
        # Every peak of pixel 0 sums above 2.5, so nothing is cut short:
        # 300 holds 11 and 600.25 holds 22 once the converter adds them.
        np.testing.assert_array_equal(mzs, FEATURE_MZ)
        np.testing.assert_array_equal(intensities, INTENSITIES[0])

    def test_extractor_reports_the_declaration(self):
        with ImzMLReader(CONTINUOUS) as reader:
            report = reader.get_comprehensive_metadata().format_specific["ion_mobility"]
        assert report["present"] is True
        assert report["array_accession"] == "MS:1003006"
        assert report["separation_accession"] == "MS:1002815"
        assert report["unit_accession"] == "MS:1002814"

    def test_offsets_are_per_spectrum_copies(self):
        with ImzMLReader(CONTINUOUS) as reader:
            reader._ensure_parser_initialized()
            offsets, lengths, encoded = collect_array_offsets(
                CONTINUOUS, "mobilityArray", 4
            )
        assert np.unique(offsets).size == 4  # the writer copies mobility per spectrum
        np.testing.assert_array_equal(lengths, [5, 5, 5, 5])
        np.testing.assert_array_equal(encoded, [40, 40, 40, 40])


class TestProcessedExport:
    def test_mobility_is_per_pixel(self):
        with ImzMLReader(PROCESSED) as reader:
            assert reader.has_ion_mobility is True
            assert reader.has_shared_mobility_axis is False
            assert reader.get_shared_mobility_features() is None
            axis = reader.get_mobility_axis()
            items = list(reader.iter_mobility_spectra())
        assert axis is not None and not axis.is_shared
        assert axis.kind_accession == "MS:1002815"
        assert len(items) == 4
        for i, (coords, mzs, mobility, intensities) in enumerate(items):
            assert coords == COORDS[i]
            np.testing.assert_array_equal(
                mzs, [100.0 + i, 100.0 + i, 250.5 + i, 400.0 + i]
            )
            np.testing.assert_array_equal(mobility, [0.8, 1.1, 0.9, 1.0])
            np.testing.assert_array_equal(
                intensities, np.array([10.0, 20.0, 30.0, 40.0]) * (i + 1)
            )

    def test_common_axis_is_the_union_of_unique_mz(self):
        with ImzMLReader(PROCESSED) as reader:
            axis = reader.get_common_mass_axis()
        expected = sorted(
            {100.0 + i for i in range(4)}
            | {250.5 + i for i in range(4)}
            | {400.0 + i for i in range(4)}
        )
        np.testing.assert_array_equal(axis, expected)


class TestMobilityValueRange:
    """The range the heatmap and the grid bin over (design decision D46)."""

    def test_a_processed_export_reads_it_from_its_pixels(self):
        with ImzMLReader(PROCESSED) as reader:
            assert reader.get_mobility_value_range() == (0.8, 1.1)
            # The axis description is unchanged: a range read from the
            # pixels is not a shared axis and not a declared range.
            axis = reader.get_mobility_axis()
        assert axis is not None and axis.values is None and axis.acq_range is None

    def test_the_pass_runs_once(self, monkeypatch):
        with ImzMLReader(PROCESSED) as reader:
            first = reader.get_mobility_value_range()

            def _no_second_read(idx):
                raise AssertionError("the mobility arrays were read twice")

            monkeypatch.setattr(reader, "_read_mobility_array", _no_second_read)
            assert reader.get_mobility_value_range() == first

    def test_a_continuous_export_takes_it_from_its_shared_axis(self):
        with ImzMLReader(CONTINUOUS) as reader:
            assert reader.get_mobility_value_range() == (0.95, 1.35)

    def test_a_file_without_mobility_has_none(self):
        with ImzMLReader(PLAIN) as reader:
            assert reader.get_mobility_value_range() is None

    def test_non_finite_values_withhold_the_mobility(self, monkeypatch, thyra_logs):
        """A NaN or infinite mobility leaves the file without mobility (D82)."""
        with ImzMLReader(PROCESSED) as reader:
            original = reader._read_mobility_array

            def _with_nan(idx):
                values = original(idx).copy()
                values[0] = np.nan
                values[1] = np.inf
                return values

            monkeypatch.setattr(reader, "_read_mobility_array", _with_nan)
            with thyra_logs("thyra.readers.imzml.imzml_reader") as records:
                assert reader.has_ion_mobility is False
                assert reader.has_shared_mobility_axis is False
                assert reader.get_mobility_value_range() is None
                assert reader.get_mobility_axis() is None
            # The finite values alone still give the span, were it asked.
            assert reader._survey_mobility().span == (0.9, 1.0)
        messages = [r.getMessage() for r in records if r.levelname == "WARNING"]
        assert len(messages) == 1, messages
        assert "8 of 16 ion mobility values are not finite" in messages[0]
        assert "mobility_processed.imzML" in messages[0]

    def test_values_that_span_nothing_give_none(self, monkeypatch):
        with ImzMLReader(PROCESSED) as reader:
            monkeypatch.setattr(
                reader, "_read_mobility_array", lambda idx: np.full(4, 1.0)
            )
            assert reader.get_mobility_value_range() is None

    def test_an_unreadable_file_declines_once_and_says_so(
        self, monkeypatch, thyra_logs
    ):
        calls = []

        def _refuse(idx):
            calls.append(idx)
            raise ConversionRefused(f"Spectrum {idx}: ion mobility array truncated")

        with ImzMLReader(PROCESSED) as reader:
            # Opening the file checks the offsets; a read fails later.
            reader._ensure_parser_initialized()
            monkeypatch.setattr(reader, "_read_mobility_array", _refuse)
            with thyra_logs("thyra.readers.imzml.imzml_reader") as records:
                assert reader.get_mobility_value_range() is None
                assert reader.get_mobility_value_range() is None
                # Unread is not shown to be non-finite: the file keeps its
                # declared mobility, with no range and no shared axis.
                assert reader.has_ion_mobility
                assert not reader.has_shared_mobility_axis
        assert len(calls) == 1
        messages = [r.getMessage() for r in records]
        assert sum("Could not read the ion mobility values" in m for m in messages) == 1


class TestNonFiniteMobility:
    def test_unplaceable_points_leave_the_cloud_and_are_counted(
        self, monkeypatch, thyra_logs
    ):
        with ImzMLReader(PROCESSED) as reader:
            original = reader._read_mobility_array

            def _with_nan(idx):
                values = original(idx).copy()
                values[0] = np.nan
                return values

            monkeypatch.setattr(reader, "_read_mobility_array", _with_nan)
            with thyra_logs("thyra.readers.imzml.imzml_reader") as records:
                items = list(reader.iter_mobility_spectra())
        assert len(items) == 4
        for i, (_, mzs, mobility, intensities) in enumerate(items):
            np.testing.assert_array_equal(mobility, [1.1, 0.9, 1.0])
            np.testing.assert_array_equal(mzs, [100.0 + i, 250.5 + i, 400.0 + i])
            np.testing.assert_array_equal(
                intensities, np.array([20.0, 30.0, 40.0]) * (i + 1)
            )
        messages = [r.getMessage() for r in records]
        assert any(m.startswith("4 point(s)") for m in messages), messages


def _write_shared_mobility_block(directory, name, term, dtype):
    """Six continuous pixels plus one mobility block they all share (F184).

    The block of four values sits right after the 16-byte UUID, declared
    with ``term`` as its precision; every byte of the ``.ibd`` is declared.
    """
    import re

    from pyimzml.ImzMLWriter import ImzMLWriter

    def cv(ref, accession, label, value="", unit=""):
        return (
            f'<cvParam cvRef="{ref}" accession="{accession}" name="{label}" '
            f'value="{value}"{unit}/>'
        )

    group = (
        '<referenceableParamGroup id="mobilityArray">'
        + cv(
            "MS",
            "MS:1003006",
            "mean inverse reduced ion mobility array",
            unit=' unitCvRef="MS" unitAccession="MS:1002814"',
        )
        + cv("IMS", "IMS:1000101", "external data", "true")
        + term
        + "</referenceableParamGroup>"
    )
    mz = np.array([150.0, 250.0, 350.0, 450.0])
    mobility = np.array([0.70, 0.85, 1.05, 1.30])
    path = directory / f"{name}.imzML"
    with ImzMLWriter(str(path), mode="continuous") as writer:
        for i in range(6):
            writer.addSpectrum(mz, (i + 1) * 10 * mz, (i % 3 + 1, i // 3 + 1, 1))
    block = mobility.astype(dtype).tobytes()
    ibd = path.with_suffix(".ibd")
    raw = ibd.read_bytes()
    ibd.write_bytes(raw[:16] + block + raw[16:])
    xml = re.sub(
        r'(external offset" value=")(\d+)',
        lambda m: m[1] + str(int(m[2]) + len(block)),
        path.read_text(),
    )
    xml = xml.replace(
        "</referenceableParamGroupList>", group + "</referenceableParamGroupList>"
    )
    array = (
        '<binaryDataArray encodedLength="0">'
        '<referenceableParamGroupRef ref="mobilityArray"/>'
        + cv("IMS", "IMS:1000102", "external offset", 16)
        + cv("IMS", "IMS:1000103", "external array length", 4)
        + cv("IMS", "IMS:1000104", "external encoded length", len(block))
        + "</binaryDataArray>"
    )
    path.write_text(
        xml.replace("</binaryDataArrayList>", array + "</binaryDataArrayList>")
    )
    return path


class TestAMobilityArrayOfAnotherPrecision:
    """F184 -- a 16-bit float array was read as 64-bit, from the neighbouring bytes."""

    @pytest.mark.parametrize(
        "label, accession, dtype",
        [
            ("32-bit float", "MS:1000521", np.float32),
            ("16-bit float", "MS:1000520", np.float16),
        ],
    )
    def test_the_values_are_read_at_the_declared_width(
        self, tmp_path, label, accession, dtype
    ):
        term = f'<cvParam cvRef="MS" accession="{accession}" name="{label}" value=""/>'
        path = _write_shared_mobility_block(tmp_path, label[:6], term, dtype)

        with ImzMLReader(path) as reader:
            assert reader.has_ion_mobility is True
            _, _, mobility, _ = next(reader.iter_mobility_spectra())

        assert mobility == pytest.approx([0.70, 0.85, 1.05, 1.30], abs=1e-3)


class TestPlainFilesAreUntouched:
    def test_two_array_file_reports_no_mobility(self):
        with ImzMLReader(PLAIN) as reader:
            assert reader.has_ion_mobility is False
            assert reader.has_shared_mobility_axis is False
            assert reader.get_mobility_axis() is None
            assert reader.get_shared_mobility_features() is None
            with pytest.raises(NotImplementedError):
                next(reader.iter_mobility_spectra())
            report = reader.get_comprehensive_metadata().format_specific["ion_mobility"]
        assert report == {"present": False}


def _metadata(groups):
    """A duck-typed pyimzml Metadata with the given param groups."""
    made = {}
    for group_id, params in groups.items():
        made[group_id] = SimpleNamespace(
            cv_params=[
                (name, accession, None, name, "", unit_name, unit_accession)
                for name, accession, unit_name, unit_accession in params
            ]
        )
    return SimpleNamespace(referenceable_param_groups=made)


class TestDetection:
    def test_mobility_group_with_unit_and_precision(self):
        spec = detect_mobility_array(
            _metadata(
                {
                    "mzArray": [("m/z array", "MS:1000514", None, None)],
                    "mobilityArray": [
                        ("no compression", "MS:1000576", None, None),
                        (
                            "mean inverse reduced ion mobility array",
                            "MS:1003006",
                            "volt-second per square centimeter",
                            "MS:1002814",
                        ),
                        ("32-bit float", "MS:1000521", None, None),
                    ],
                }
            )
        )
        assert spec is not None
        assert spec.group_id == "mobilityArray"
        assert spec.dtype == np.dtype(np.float32)
        assert spec.unit_accession == "MS:1002814"
        assert spec.compressed is False

    def test_zlib_is_flagged(self):
        spec = detect_mobility_array(
            _metadata(
                {
                    "mob": [
                        ("zlib compression", "MS:1000574", None, None),
                        ("ion mobility array", "MS:1002893", None, None),
                    ]
                }
            )
        )
        assert spec is not None and spec.compressed is True
        assert spec.dtype == np.dtype(np.float64)  # undeclared precision

    def test_a_16_bit_float_array_is_read_as_16_bit(self):
        """MS:1000520 was missing from the table, so the array read as 64-bit."""
        spec = detect_mobility_array(
            _metadata(
                {
                    "mobilityArray": [
                        (
                            "mean inverse reduced ion mobility array",
                            "MS:1003006",
                            None,
                            None,
                        ),
                        ("16-bit float", "MS:1000520", None, None),
                    ]
                }
            )
        )
        assert spec is not None
        assert spec.dtype == np.float16
        assert spec.dtype.itemsize == 2

    def test_a_binary_type_thyra_cannot_decode_is_refused(self):
        """MS:1001479 (null-terminated ASCII string) fell through to 64-bit."""
        with pytest.raises(ConversionRefused, match="a precision Thyra cannot decode"):
            detect_mobility_array(
                _metadata(
                    {
                        "mobilityArray": [
                            ("ion mobility array", "MS:1002893", None, None),
                            ("null-terminated ASCII string", "MS:1001479", None, None),
                        ]
                    }
                )
            )

    def test_two_precision_terms_are_refused(self):
        with pytest.raises(ConversionRefused, match="exactly one is required"):
            detect_mobility_array(
                _metadata(
                    {
                        "mobilityArray": [
                            ("ion mobility array", "MS:1002893", None, None),
                            ("32-bit float", "MS:1000521", None, None),
                            ("64-bit float", "MS:1000523", None, None),
                        ]
                    }
                )
            )

    def test_no_mobility_group(self):
        assert (
            detect_mobility_array(
                _metadata({"mzArray": [("m/z array", "MS:1000514", None, None)]})
            )
            is None
        )
        assert detect_mobility_array(None) is None

    @pytest.mark.parametrize(
        "array, unit, expected",
        [
            ("MS:1003006", None, "MS:1002815"),
            ("MS:1002893", "MS:1002814", "MS:1002815"),
            ("MS:1002893", "UO:0000028", "MS:1002476"),
            ("MS:1003153", None, "MS:1002476"),
            ("MS:1002893", None, None),
        ],
    )
    def test_classification(self, array, unit, expected):
        kind, name = classify_mobility_array(array, unit)
        assert kind == expected
        assert (name is None) == (expected is None)

    def test_offsets_refuse_a_spectrum_count_mismatch(self):
        with pytest.raises(ConversionRefused, match="pyimzml reported"):
            collect_array_offsets(CONTINUOUS, "mobilityArray", 3)


FEATURES_MZ = np.array([210.1, 210.1, 305.5, 401.25, 512.0, 512.0, 640.3, 777.7])
FEATURES_K0 = np.array([0.80, 1.10, 0.95, 1.02, 0.70, 1.25, 1.05, 1.30])
THREE_MZ = np.array([200.0, 300.0, 400.0])


class TestTheSharedAxisIsProvenOnEverySpectrum:
    """Q8 (D82): one m/z block, and every mobility copy equal to the first."""

    @pytest.mark.parametrize("odd", [71, 72, 143])
    def test_a_late_copy_that_differs_makes_the_file_per_pixel(self, tmp_path, odd):
        """Spectrum 71 of 144 was never compared by the 64-spectrum sample (F185)."""
        mobility = [FEATURES_K0 + (0.15 if i == odd else 0.0) for i in range(144)]
        path = write_mobility_imzml(
            tmp_path / "odd.imzML", [FEATURES_MZ] * 144, mobility
        )
        with ImzMLReader(path) as reader:
            assert reader.has_ion_mobility is True
            assert reader.has_shared_mobility_axis is False
            assert reader.get_shared_mobility_features() is None
            axis = reader.get_mobility_axis()
            items = list(reader.iter_mobility_spectra())
        assert axis is not None and axis.values is None
        np.testing.assert_allclose(items[odd][2], FEATURES_K0 + 0.15)

    def test_equal_copies_share_the_axis(self, tmp_path):
        path = write_mobility_imzml(
            tmp_path / "even.imzML", [FEATURES_MZ] * 80, [FEATURES_K0] * 80
        )
        with ImzMLReader(path) as reader:
            assert reader.has_shared_mobility_axis is True
            axis = reader.get_mobility_axis()
        assert axis is not None
        np.testing.assert_array_equal(axis.values, FEATURES_K0)

    def test_each_mobility_block_is_read_once(self, tmp_path, monkeypatch):
        path = write_mobility_imzml(
            tmp_path / "copies.imzML", [FEATURES_MZ] * 20, [FEATURES_K0] * 20
        )
        reads = []
        with ImzMLReader(path) as reader:
            original = reader._read_mobility_array

            def _counted(idx):
                reads.append(idx)
                return original(idx)

            monkeypatch.setattr(reader, "_read_mobility_array", _counted)
            assert reader.has_shared_mobility_axis is True
            assert reader.get_mobility_value_range() == (0.70, 1.30)
        assert sorted(reads) == list(range(20))

    def test_a_continuous_declaration_over_own_m_z_blocks_is_per_pixel(
        self, tmp_path, thyra_logs
    ):
        """F186: spectrum 9's own m/z was filed under spectrum 0's."""
        mzs = [FEATURES_MZ + (2.0 if i == 9 else 0.0) for i in range(16)]
        path = write_mobility_imzml(
            tmp_path / "own_mz.imzML", mzs, [FEATURES_K0] * 16, one_mz_block=False
        )
        with thyra_logs("thyra.readers.imzml.imzml_reader") as records:
            with ImzMLReader(path) as reader:
                assert reader.has_shared_mass_axis is False
                assert reader.is_continuous is True
                assert reader.has_shared_mobility_axis is False
                items = list(reader.iter_mobility_spectra())
        np.testing.assert_array_equal(items[9][1], FEATURES_MZ + 2.0)
        messages = [r.getMessage() for r in records]
        assert any(
            "declares continuous mode, but its spectra point at 16" in m
            for m in messages
        ), messages


class TestNonFiniteMobilityWithholdsTheMobility:
    """Q13 (D82): a NaN or infinite mobility writes no mobility products."""

    def test_byte_identical_copies_holding_a_nan_are_not_said_to_differ(
        self, tmp_path, thyra_logs
    ):
        """F205: NaN != NaN made equal copies read as per-pixel mobility."""
        k0 = np.array([0.8, 1.0, np.nan])
        path = write_mobility_imzml(tmp_path / "nan.imzML", [THREE_MZ] * 16, [k0] * 16)
        with thyra_logs("thyra.readers.imzml.imzml_reader") as records:
            with ImzMLReader(path) as reader:
                assert reader._survey_mobility().shared is not None
                assert reader.has_ion_mobility is False
                assert reader.get_mobility_axis() is None
        messages = [r.getMessage() for r in records]
        assert not any("differ between spectra" in m for m in messages), messages
        said = [m for m in messages if "not finite numbers" in m]
        assert len(said) == 1, messages
        assert "16 of 48 ion mobility values" in said[0]
        assert "nan.imzML" in said[0]

    @pytest.mark.parametrize("bad", [np.inf, -np.inf, np.nan])
    def test_one_shared_block_with_a_bad_value(self, tmp_path, bad):
        """F168: the block every spectrum points at is never re-read or compared."""
        k0 = np.array([0.8, bad, 1.2])
        path = write_mobility_imzml(
            tmp_path / "block.imzML", [THREE_MZ] * 4, [k0] * 4, one_mobility_block=True
        )
        with ImzMLReader(path) as reader:
            assert reader.has_ion_mobility is False
            assert reader.has_shared_mobility_axis is False
            assert reader.get_mobility_value_range() is None

    def test_finite_values_keep_their_mobility(self, tmp_path, monkeypatch):
        """One block for every spectrum is read once, and is the shared axis."""
        k0 = np.array([0.8, 1.0, 1.2])
        path = write_mobility_imzml(
            tmp_path / "fine.imzML", [THREE_MZ] * 4, [k0] * 4, one_mobility_block=True
        )
        reads = []
        with ImzMLReader(path) as reader:
            original = reader._read_mobility_array

            def _counted(idx):
                reads.append(idx)
                return original(idx)

            monkeypatch.setattr(reader, "_read_mobility_array", _counted)
            assert reader.has_ion_mobility is True
            assert reader.has_shared_mobility_axis is True
            assert reader.get_mobility_value_range() == (0.8, 1.2)
        assert reads == [0]
