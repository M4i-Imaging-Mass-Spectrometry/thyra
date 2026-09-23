"""A store carries no person its source names, and none of its source's folders.

Vendor software records who calibrated the instrument and who ran it, and
paths on the acquisition PC whose folders carry a user's name. None of it
does anything in a store, and the source still holds all of it; see
:mod:`thyra.metadata.personal_data`. These tests pin the rule itself, that
every extractor passes through it, and two real table layouts that carry a
person. The PHI and solariX conversions are checked end to end beside their
own fixtures, in ``tests/unit/readers``.
"""

import json
import sqlite3
from pathlib import Path

import numpy as np
import pytest

from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra.core.base_extractor import MetadataExtractor
from thyra.metadata.extractors.bruker_extractor import BrukerMetadataExtractor
from thyra.metadata.personal_data import (
    file_name,
    is_path,
    strip_personal_data,
    without_personal_data,
)
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata
from thyra.readers.mzpeak import MzPeakReader

FIXTURE = Path(__file__).resolve().parents[2] / "data" / "fixtures" / "synthetic_tims.d"

# A name and a folder that appear nowhere else, so finding either one
# anywhere in a serialised result means it leaked.
PERSON = "Qwerty Operatorperson"
FOLDER = "qwertyfolder"


def _nowhere_in(result) -> bool:
    blob = json.dumps(result, default=str)
    return PERSON not in blob and FOLDER not in blob


class TestPeopleAreNotCopied:
    @pytest.mark.parametrize(
        "key",
        [
            "OperatorName",  # Bruker GlobalMetadata, solariX Properties
            "Operator",  # PHI header
            "User Name",  # PHI [Data Manager]
            "CalibrationUser",  # timsTOF calibration.sqlite
            "MobilityCalibrationUser",  # timsTOF calibration.sqlite
            "operator_name",  # Thyra's own spelling, solariX
            "calibration_user",  # Thyra's own spelling, timsTOF
            "user-name",
        ],
    )
    def test_every_spelling_of_a_person_is_dropped(self, key):
        assert strip_personal_data({key: PERSON, "kept": 1}) == {"kept": 1}

    @pytest.mark.parametrize(
        "key",
        ["User Company", "Operator Mode", "resampling_operator", "SampleName"],
    )
    def test_a_key_that_only_mentions_one_is_kept(self, key):
        # An organisation, a mode, Thyra's own vocabulary, and free text
        # the lab chose -- none of them is a person field.
        assert strip_personal_data({key: "value"}) == {key: "value"}

    def test_a_person_is_dropped_at_any_depth(self):
        result = strip_personal_data(
            {
                "global_metadata": {"OperatorName": PERSON, "InstrumentName": "t"},
                "states": [{"CalibrationUser": PERSON, "Id": 3}],
            }
        )
        assert result == {
            "global_metadata": {"InstrumentName": "t"},
            "states": [{"Id": 3}],
        }

    def test_the_psi_contact_terms_are_dropped_and_the_affiliation_kept(self):
        params = [
            {"accession": "MS:1000586", "name": "contact name", "value": PERSON},
            {"accession": "MS:1000589", "name": "contact email", "value": "q@x.org"},
            {"accession": "MS:1001755", "name": "contact phone number", "value": "1"},
            {"accession": "MS:1000590", "name": "contact affiliation", "value": "Uni"},
        ]
        assert strip_personal_data({"contacts": params}) == {"contacts": [params[3]]}

    def test_an_email_address_is_dropped_whatever_its_key(self):
        result = strip_personal_data({"Notify": "qwerty.person@example.org", "n": 1})
        assert result == {"n": 1}


class TestPathsKeepTheirFileName:
    @pytest.mark.parametrize(
        "value, expected",
        [
            (f"D:\\SmartSoft\\Data\\{FOLDER}\\run42.raw", "run42.raw"),
            (f"D:/Methods/{FOLDER}/lipids.m", "lipids.m"),
            (f"\\\\labserver\\share\\{FOLDER}\\run.d", "run.d"),
            (f"{FOLDER}\\run.raw", "run.raw"),
            (f"file:///home/{FOLDER}/raw/run.mzML", "run.mzML"),
            (f"/home/{FOLDER}/run.imzML", "run.imzML"),
            (f"C:\\Data\\{FOLDER}\\data\\sample.d\\", "sample.d"),
        ],
    )
    def test_a_path_becomes_its_last_component(self, value, expected):
        assert is_path(value)
        assert strip_personal_data({"Recorded": value}) == {"Recorded": expected}

    @pytest.mark.parametrize(
        "value",
        [
            "m/z",
            "1/K0",
            "Vs/cm2",
            "06/23/2026 21:12:35",  # PHI AcqFileDate, month first
            "2025-04-22T08:59:34.395+02:00",
            "https://example.org/a/b",
            "timsTOF fleX",
        ],
    )
    def test_what_only_looks_like_a_path_is_left_alone(self, value):
        assert not is_path(value)
        assert strip_personal_data({"Recorded": value}) == {"Recorded": value}

    @pytest.mark.parametrize(
        "document",
        [
            {"Raster Size (µm)": "8.0"},  # json.dumps escapes the micro sign
            {"Note": 'probe "A"'},  # and a quote
            ["a", "b"],
        ],
    )
    def test_a_json_document_is_not_a_path(self, document):
        # PHI keeps its header as JSON. JSON escapes start with a
        # backslash, and a rule that took any backslash for a path cut
        # the whole document down to what followed the last one.
        text = json.dumps(document)
        assert not is_path(text)
        assert strip_personal_data({"header": text}) == {"header": text}

    def test_the_paths_thyra_records_about_its_source_stay_whole(self):
        # Where the source was read, on the machine that holds the store:
        # it finds the source again, which is why #384 kept it.
        own = {
            "data_path": f"C:\\Data\\{FOLDER}\\x.d",
            "database_path": f"C:\\Data\\{FOLDER}\\x.d\\analysis.tdf",
            "binary_file": f"C:\\Data\\{FOLDER}\\x.d\\analysis.tdf_bin",
            "ibd_file": f"C:\\Data\\{FOLDER}\\x.ibd",
            "mis_file": f"C:\\Data\\{FOLDER}\\x.mis",
        }
        assert strip_personal_data(own) == own

    def test_only_at_the_top_level_of_a_dictionary(self):
        # Nested, the same key is the vendor's, and a vendor path it is.
        nested = {"vendor": {"data_path": f"C:\\Data\\{FOLDER}\\x.d"}}
        assert strip_personal_data(nested) == {"vendor": {"data_path": "x.d"}}


class TestTheCopy:
    def test_the_input_is_not_modified(self):
        source = {"Operator": PERSON, "nested": {"Path": f"C:\\{FOLDER}\\a.raw"}}
        before = json.dumps(source)
        strip_personal_data(source)
        assert json.dumps(source) == before

    def test_containers_keep_their_type_and_other_values_pass_through(self):
        array = np.arange(3)
        result = strip_personal_data(
            {
                "pair": (1.0, f"C:\\{FOLDER}\\a.raw"),
                "list": [2, None, True],
                "array": array,
            }
        )
        assert result["pair"] == (1.0, "a.raw")
        assert result["list"] == [2, None, True]
        assert result["array"] is array


class TestFileName:
    def test_it_is_the_reduction_method_file_always_used(self):
        # The builder's helper moved here; what it returns did not move.
        assert file_name("D:\\Methods\\x.m") == "x.m"
        assert file_name("D:/Methods/x.m") == "x.m"
        assert file_name("x.m") == "x.m"
        assert file_name("D:\\Methods\\") is None


class _Reporting(MetadataExtractor):
    """An extractor that reports a person and a path in every section."""

    def _extract_essential_impl(self) -> EssentialMetadata:
        return EssentialMetadata(
            dimensions=(1, 1, 1),
            coordinate_bounds=(0.0, 0.0, 0.0, 0.0),
            mass_range=(100.0, 200.0),
            pixel_size=None,
            n_spectra=1,
            total_peaks=1,
            source_path=f"C:\\Data\\{FOLDER}\\sample.d",
        )

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        people = {"OperatorName": PERSON, "Method": f"D:\\{FOLDER}\\m.m", "n": 1}
        return ComprehensiveMetadata(
            essential=self.get_essential(),
            format_specific=dict(people),
            acquisition_params=dict(people),
            instrument_info=dict(people),
            raw_metadata={"dump": dict(people)},
        )


class TestEveryExtractorPassesThroughIt:
    def test_get_comprehensive_returns_no_person_and_no_folder(self):
        comprehensive = _Reporting(None).get_comprehensive()
        for section in (
            comprehensive.format_specific,
            comprehensive.acquisition_params,
            comprehensive.instrument_info,
            comprehensive.raw_metadata,
        ):
            assert _nowhere_in(section)
        assert comprehensive.acquisition_params == {"Method": "m.m", "n": 1}

    def test_the_source_path_stays_whole(self):
        # The one path a store keeps: #384.
        comprehensive = _Reporting(None).get_comprehensive()
        assert comprehensive.essential.source_path.endswith(f"{FOLDER}\\sample.d")

    def test_without_personal_data_leaves_a_missing_section_missing(self):
        extractor = _Reporting(None)
        comprehensive = extractor._extract_comprehensive_impl()
        comprehensive.instrument_info = None
        assert without_personal_data(comprehensive).instrument_info is None


class TestRealTableLayouts:
    def test_the_bruker_global_metadata_dump_has_no_operator(self):
        # The hand-written TDF fixture carries OperatorName in
        # GlobalMetadata, the way every timsTOF acquisition does; the raw
        # dump copies that table whole.
        uri = "file:" + (FIXTURE / "analysis.tdf").as_posix() + "?mode=ro"
        conn = sqlite3.connect(uri, uri=True)
        try:
            raw = (
                BrukerMetadataExtractor(conn, FIXTURE).get_comprehensive().raw_metadata
            )
        finally:
            conn.close()
        global_metadata = raw["global_metadata"]
        assert "OperatorName" not in global_metadata
        assert global_metadata["InstrumentName"]  # the rest is still there

    def test_mzpeak_contact_terms_and_source_folders_are_not_carried(self, tmp_path):
        archive = build_mzpeak(
            tmp_path / "people.mzpeak",
            grid_spectra(2, 2),
            index_metadata={
                "contacts": [
                    {
                        "accession": "MS:1000586",
                        "name": "contact name",
                        "value": PERSON,
                    },
                    {
                        "accession": "MS:1000590",
                        "name": "contact affiliation",
                        "value": "Example University",
                    },
                ],
                "source_location": f"file:///C:/Data/{FOLDER}/raw/run.raw",
            },
        )
        with MzPeakReader(archive) as reader:
            raw = reader.get_comprehensive_metadata().raw_metadata
        assert _nowhere_in(raw)
        assert raw["contacts"] == [
            {
                "accession": "MS:1000590",
                "name": "contact affiliation",
                "value": "Example University",
            }
        ]
        assert raw["source_location"] == "run.raw"
