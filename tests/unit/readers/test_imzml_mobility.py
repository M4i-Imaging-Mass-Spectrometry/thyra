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

from thyra.core.mobility import MobilityAxis, classify_mobility_array
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
        with ImzMLReader(CONTINUOUS, intensity_threshold=2.5) as reader:
            _, mzs, mobility, intensities = next(reader.iter_mobility_spectra())
        np.testing.assert_array_equal(mzs, [300.0, 450.5, 600.25])
        np.testing.assert_array_equal(mobility, [0.95, 1.02, 1.20])
        np.testing.assert_array_equal(intensities, [10.0, 5.0, 20.0])

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
        with pytest.raises(ValueError, match="pyimzml reported"):
            collect_array_offsets(CONTINUOUS, "mobilityArray", 3)
