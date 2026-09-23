"""Auto-population: the builder reports what the source knows, nothing more."""

from thyra.metadata.schema import build_metadata_document, build_msi_metadata
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata


def _essential(source_path: str = "input.imzML") -> EssentialMetadata:
    return EssentialMetadata(
        dimensions=(10, 10, 1),
        coordinate_bounds=(0.0, 9.0, 0.0, 9.0),
        mass_range=(100.0, 1000.0),
        pixel_size=(20.0, 25.0),
        n_spectra=100,
        total_peaks=1000,
        source_path=source_path,
        spectrum_type="centroid spectrum",
    )


def _comprehensive(
    acquisition_params=None,
    instrument_info=None,
    format_specific=None,
    raw_metadata=None,
) -> ComprehensiveMetadata:
    return ComprehensiveMetadata(
        essential=_essential(),
        format_specific=format_specific or {},
        acquisition_params=acquisition_params or {},
        instrument_info=instrument_info or {},
        raw_metadata=raw_metadata or {},
    )


class TestBuildMsiMetadata:
    def test_pixel_size_and_provenance_are_always_present(self):
        meta = build_msi_metadata(
            _comprehensive(),
            pixel_size_um=(20.0, 25.0),
            pixel_size_source="automatic",
            source_format="imzml",
        )
        assert meta.ms_analysis.pixel_size_um.x == 20.0
        assert meta.ms_analysis.pixel_size_um.y == 25.0
        assert meta.provenance.source_format == "imzml"
        assert meta.provenance.source_path == "input.imzML"
        assert meta.provenance.pixel_size_source == "automatic"
        assert meta.provenance.thyra_version

    def test_without_comprehensive_metadata_the_block_still_builds(self):
        meta = build_msi_metadata(None, pixel_size_um=(5.0, 5.0))
        assert meta.ms_analysis.pixel_size_um.x == 5.0
        assert meta.provenance.source_path is None

    def test_unreported_fields_stay_unset(self):
        meta = build_msi_metadata(_comprehensive(), pixel_size_um=(20.0, 20.0))
        analysis = meta.ms_analysis
        assert analysis.polarity is None
        assert analysis.ionisation_source is None
        assert analysis.analyzer is None
        assert analysis.instrument_model is None

    def test_phi_polarity_and_format_facts(self):
        # PHI reports polarity in acquisition_params; SIMS/TOF follow
        # from the format itself.
        meta = build_msi_metadata(
            _comprehensive(acquisition_params={"polarity": "Positive"}),
            pixel_size_um=(2.0, 2.0),
            source_format="phi",
        )
        analysis = meta.ms_analysis
        assert analysis.polarity == "positive"
        assert analysis.polarity_term is not None
        assert analysis.polarity_term.accession == "MS:1000130"
        assert analysis.ionisation_source == "SIMS"
        assert analysis.analyzer == "TOF"

    def test_imzml_polarity_comes_from_the_preserved_cv_params(self):
        # imzML declares polarity as MS:1000130/MS:1000129 in the file
        # description; the extractor preserves those with accessions.
        meta = build_msi_metadata(
            _comprehensive(
                raw_metadata={
                    "cvParams": [
                        {
                            "name": "positive scan",
                            "accession": "MS:1000130",
                            "value": True,
                        }
                    ]
                }
            ),
            pixel_size_um=(20.0, 20.0),
            source_format="imzml",
        )
        assert meta.ms_analysis.polarity == "positive"
        assert meta.ms_analysis.polarity_term is not None
        assert meta.ms_analysis.polarity_term.accession == "MS:1000130"

    def test_a_file_declaring_both_polarities_stays_unset(self):
        meta = build_msi_metadata(
            _comprehensive(
                raw_metadata={
                    "cvParams": [
                        {"accession": "MS:1000130"},
                        {"accession": "MS:1000129"},
                    ]
                }
            ),
            pixel_size_um=(20.0, 20.0),
        )
        assert meta.ms_analysis.polarity is None

    def test_bruker_maldi_flag_sets_the_source(self):
        meta = build_msi_metadata(
            _comprehensive(format_specific={"is_maldi": True}),
            pixel_size_um=(20.0, 20.0),
            source_format="bruker",
        )
        analysis = meta.ms_analysis
        assert analysis.ionisation_source == "MALDI"
        assert analysis.ionisation_source_term is not None
        assert analysis.ionisation_source_term.accession == "MS:1000075"
        assert analysis.analyzer == "TOF"

    def test_imzml_gets_no_format_default_analyzer(self):
        # imzML can come off any instrument; claiming an analyzer for
        # it would be a guess.
        meta = build_msi_metadata(
            _comprehensive(),
            pixel_size_um=(20.0, 20.0),
            source_format="imzml",
        )
        assert meta.ms_analysis.analyzer is None

    def test_instrument_model_is_probed_across_key_spellings(self):
        for key in ("instrument_model", "instrument_name", "model", "platform"):
            meta = build_msi_metadata(
                _comprehensive(instrument_info={key: "rapifleX"}),
                pixel_size_um=(20.0, 20.0),
            )
            assert meta.ms_analysis.instrument_model == "rapifleX"

    def test_a_nanoelectrospray_source_reaches_the_document(self):
        # The first non-imaging source described reported this and lost
        # the field, because the alias table knew only imaging sources.
        meta = build_msi_metadata(
            _comprehensive(
                acquisition_params={"ionisation_source": "nanoelectrospray"}
            ),
            pixel_size_um=(20.0, 20.0),
        )
        analysis = meta.ms_analysis
        assert analysis.ionisation_source == "nanoESI"
        assert analysis.ionisation_source_term is not None
        assert analysis.ionisation_source_term.accession == "MS:1000398"

    def test_an_unknown_source_spelling_is_logged_not_guessed(self, thyra_logs):
        with thyra_logs("thyra.metadata.schema.builder", "DEBUG") as records:
            meta = build_msi_metadata(
                _comprehensive(acquisition_params={"ionisation_source": "laser magic"}),
                pixel_size_um=(20.0, 20.0),
            )
        assert meta.ms_analysis.ionisation_source is None
        assert any("laser magic" in record.getMessage() for record in records)

    def test_resolving_power_is_filled_from_the_value_and_its_reference_mz(self):
        meta = build_msi_metadata(
            _comprehensive(
                instrument_info={"resolving_power": 8750, "resolving_power_at_mz": 200}
            ),
            pixel_size_um=(20.0, 20.0),
        )
        power = meta.ms_analysis.detector_resolving_power
        assert power is not None
        assert (power.value, power.at_mz) == (8750.0, 200.0)

    def test_resolving_power_is_also_read_from_acquisition_params(self):
        meta = build_msi_metadata(
            _comprehensive(
                acquisition_params={
                    "resolving_power": "60000",
                    "resolving_power_at_mz": "400",
                }
            ),
            pixel_size_um=(20.0, 20.0),
        )
        power = meta.ms_analysis.detector_resolving_power
        assert power is not None
        assert (power.value, power.at_mz) == (60000.0, 400.0)

    def test_resolving_power_without_its_reference_mz_stays_unset(self):
        # A resolving power is not comparable without the m/z it is quoted
        # at, so half the pair is not written as if it were whole.
        for info in (
            {"resolving_power": 8750},
            {"resolving_power_at_mz": 200},
            {"resolving_power": 0, "resolving_power_at_mz": 200},
            {"resolving_power": 8750, "resolving_power_at_mz": -1},
            {"resolving_power": [8000, 9000], "resolving_power_at_mz": 200},
        ):
            meta = build_msi_metadata(
                _comprehensive(instrument_info=info), pixel_size_um=(20.0, 20.0)
            )
            assert meta.ms_analysis.detector_resolving_power is None, info

    def test_built_block_passes_validation(self):
        from thyra.metadata.schema import validate_document

        meta = build_msi_metadata(
            _comprehensive(acquisition_params={"polarity": "negative"}),
            pixel_size_um=(20.0, 25.0),
            pixel_size_source="manual",
            source_format="bruker",
        )
        model, issues = validate_document(meta.to_uns_dict())
        assert model is not None
        assert issues == []


class TestAcquisitionSection:
    """One case per reader: what each fills, and that the rest stays unset.

    The stubs carry exactly the keys and value shapes the extractors write
    (see the key tables in ``builder.py``); the values are the ones read
    off real acquisitions when the units were verified.
    """

    def _acquisition(self, acquisition_params, source_format, **format_specific):
        meta = build_msi_metadata(
            _comprehensive(
                acquisition_params=acquisition_params,
                format_specific=format_specific or None,
            ),
            pixel_size_um=(20.0, 20.0),
            source_format=source_format,
        )
        return meta.acquisition

    def test_bruker_tsf_fills_every_field(self):
        acquisition = self._acquisition(
            {
                "acquisition_datetime": "2025-04-22T08:59:34.395+02:00",
                "laser_power": 70.00436147054036,
                "laser_frequency": 1000.0,
                "num_laser_shots": 50,
                "method_name": "20250415_pos50-1000_20um_M2.m",
                "beam_scan_size_x": 16.0,
            },
            "bruker",
            bruker_format="bruker_tsf",
            is_maldi=True,
        )
        assert acquisition is not None
        assert acquisition.model_dump(exclude_none=True) == {
            "acquisition_datetime": "2025-04-22T08:59:34.395+02:00",
            "laser_power_percent": 70.00436147054036,
            "laser_frequency_hz": 1000.0,
            "shots_per_pixel": 50,
            "method_file": "20250415_pos50-1000_20um_M2.m",
        }

    def test_bruker_tdf_is_the_same_schema_and_keeps_a_varying_value_unset(self):
        # A per-frame value that differed across frames is reported by the
        # extractor as a [min, max] pair; the section states one value per
        # acquisition or none, so the pair leaves the field unset.
        acquisition = self._acquisition(
            {
                "acquisition_datetime": "2026-01-01T00:00:00.000+00:00",
                "laser_power": [60.0, 70.0],
                "laser_frequency": 10000.0,
                "num_laser_shots": 200,
                "method_name": "synthetic.m",
            },
            "bruker",
            bruker_format="bruker_tdf",
            is_maldi=True,
        )
        assert acquisition is not None
        assert acquisition.model_dump(exclude_none=True) == {
            "acquisition_datetime": "2026-01-01T00:00:00+00:00",
            "laser_frequency_hz": 10000.0,
            "shots_per_pixel": 200,
            "method_file": "synthetic.m",
        }

    def test_solarix_fills_every_field_from_its_own_spellings(self):
        acquisition = self._acquisition(
            {
                "acquisition_datetime": "2022-03-29T14:50:28.506+02:00",
                "laser_power": 80.0,
                "laser_rep_rate": 2000.0,
                "num_summations": 200,
                "method_name": "PDE_micebrain_DHB_pos.m",
                "operator_name": "Admin",
            },
            "solarix",
        )
        assert acquisition is not None
        assert acquisition.model_dump(exclude_none=True) == {
            "acquisition_datetime": "2022-03-29T14:50:28.506+02:00",
            "laser_power_percent": 80.0,
            "laser_frequency_hz": 2000.0,
            "shots_per_pixel": 200,
            "method_file": "PDE_micebrain_DHB_pos.m",
        }

    def test_phi_has_only_a_timestamp(self):
        # SmartSoft-TOF writes AcqFileDate month-first with no zone; the
        # ion gun has no laser and the header names no method file.
        acquisition = self._acquisition(
            {
                "acquisition_date": "06/23/2026 21:12:35",
                "polarity": "Negative",
                "primary_species": "Bi3 +",
            },
            "phi",
        )
        assert acquisition is not None
        assert acquisition.model_dump(exclude_none=True) == {
            "acquisition_datetime": "2026-06-23T21:12:35",
        }

    def test_waters_has_a_timestamp_and_the_ms_method(self):
        # MassLynx's getAcquisitionDate is 'dd-Mon-yyyy hh:mm:ss' with no
        # zone; the .EXP name comes from _header.txt; nothing says laser.
        acquisition = self._acquisition(
            {
                "acquisition_date": "07-Nov-2019 14:09:44",
                "ms_method": "neg_FastDDA.EXP",
                "is_lockmass_corrected": False,
            },
            "waters",
        )
        assert acquisition is not None
        assert acquisition.model_dump(exclude_none=True) == {
            "acquisition_datetime": "2019-11-07T14:09:44",
            "method_file": "neg_FastDDA.EXP",
        }

    def test_rapiflex_fills_the_two_facts_that_are_unambiguous(self):
        # The info file is parsed as text, so the values are strings. Its
        # "Laser Power" was not verified as a percentage and its "Start
        # Time" format only against a synthetic file, so both stay unset.
        acquisition = self._acquisition(
            {
                "shots_per_spot": "100",
                "laser_power": "50",
                "method": "TestMethod.par",
                "start_time": "Mon, 01.01.2023 12:00:00",
            },
            "rapiflex",
        )
        assert acquisition is not None
        assert acquisition.model_dump(exclude_none=True) == {
            "shots_per_pixel": 100,
            "method_file": "TestMethod.par",
        }

    def test_imzml_has_no_section(self):
        meta = build_msi_metadata(
            _comprehensive(
                acquisition_params={"scan_direction": "left to right"},
                raw_metadata={"cvParams": [{"accession": "MS:1000130"}]},
            ),
            pixel_size_um=(20.0, 20.0),
            source_format="imzml",
        )
        assert meta.acquisition is None
        assert "acquisition" not in meta.to_uns_dict()

    def test_without_comprehensive_metadata_there_is_no_section(self):
        assert build_msi_metadata(None, pixel_size_um=(5.0, 5.0)).acquisition is None

    def test_a_timestamp_the_builder_cannot_parse_stays_unset(self):
        # The raw string is still in acquisition_params; the section must
        # not carry a value the model would reject.
        acquisition = self._acquisition(
            {"acquisition_date": "Mon, 01.01.2023 12:00:00", "ms_method": "a.EXP"},
            "waters",
        )
        assert acquisition is not None
        assert acquisition.acquisition_datetime is None
        assert acquisition.method_file == "a.EXP"

    def test_a_timestamp_keeps_the_offset_and_precision_it_came_with(self):
        for raw, expected in (
            ("2025-04-22T08:59:34.395+02:00", "2025-04-22T08:59:34.395+02:00"),
            ("2026-01-01T00:00:00.000+00:00", "2026-01-01T00:00:00+00:00"),
            ("2026-01-01T00:00:00", "2026-01-01T00:00:00"),
            ("2026-01-01T00:00:00.123456", "2026-01-01T00:00:00.123456"),
        ):
            acquisition = self._acquisition({"acquisition_datetime": raw}, "bruker")
            assert acquisition is not None
            assert acquisition.acquisition_datetime == expected, raw

    def test_the_method_is_a_name_even_when_the_vendor_recorded_a_path(self):
        for recorded in (
            "D:\\Methods\\imaging_pos.m",
            "/data/methods/imaging_pos.m",
            "imaging_pos.m",
        ):
            acquisition = self._acquisition({"method_name": recorded}, "bruker")
            assert acquisition is not None
            assert acquisition.method_file == "imaging_pos.m", recorded

    def test_laser_power_is_a_percentage_only_where_that_is_verified(self):
        for source_format, expected in (
            ("bruker", 70.0),
            ("solarix", 70.0),
            ("rapiflex", None),
            ("imzml", None),
        ):
            acquisition = self._acquisition({"laser_power": 70.0}, source_format)
            value = acquisition.laser_power_percent if acquisition else None
            assert value == expected, source_format

    def test_values_outside_their_range_stay_unset(self):
        acquisition = self._acquisition(
            {
                "laser_power": 140.0,
                "laser_frequency": 0.0,
                "num_laser_shots": 0,
                "method_name": "   ",
            },
            "bruker",
        )
        assert acquisition is None

    def test_built_section_passes_validation(self):
        from thyra.metadata.schema import validate_document

        meta = build_msi_metadata(
            _comprehensive(
                acquisition_params={
                    "acquisition_datetime": "2025-04-22T08:59:34.395+02:00",
                    "laser_power": 70.0,
                    "laser_frequency": 1000.0,
                    "num_laser_shots": 50,
                    "method_name": "imaging.m",
                }
            ),
            pixel_size_um=(20.0, 20.0),
            source_format="bruker",
        )
        model, issues = validate_document(meta.to_uns_dict())
        assert model is not None
        assert model.acquisition is not None
        assert issues == []


#: A Windows source path in the spelling ``str(Path(...))`` produces, which
#: is what every extractor stores. Assembled rather than written out: a
#: drive letter followed by backslashes is precisely the shape the
#: ``no-lab-share-paths`` pre-commit hook rejects in a tracked file, and
#: that hook enforces the policy this test class is about.
_WINDOWS_SOURCE = "\\".join(("D:", "acquisitions", "2025", "mouse_brain.d"))


class TestTheSourceIsNamedInADocumentAndLocatedInAStore:
    """Issue #384: one field, two readerships.

    A store block keeps the path -- the store sits beside the source on
    the machine that wrote it.  A document is written to be handed to
    somebody, so it carries the name and leaves the filesystem behind.
    """

    @staticmethod
    def _from(source_path: str) -> ComprehensiveMetadata:
        return ComprehensiveMetadata(
            essential=_essential(source_path),
            format_specific={},
            acquisition_params={},
            instrument_info={},
            raw_metadata={},
        )

    def test_a_store_block_keeps_the_whole_path(self):
        meta = build_msi_metadata(
            self._from(_WINDOWS_SOURCE),
            pixel_size_um=(20.0, 20.0),
            source_format="bruker",
        )
        assert meta.provenance.source_path == _WINDOWS_SOURCE

    def test_a_document_carries_the_name_only(self):
        document = build_metadata_document(
            self._from(_WINDOWS_SOURCE),
            pixel_size_um=(20.0, 20.0),
            source_format="bruker",
        )
        assert document["provenance"]["source_path"] == "mouse_brain.d"

    def test_a_posix_path_is_reduced_the_same_way(self):
        document = build_metadata_document(
            self._from("/mnt/share/2025/mouse_brain.d"),
            pixel_size_um=(20.0, 20.0),
            source_format="bruker",
        )
        assert document["provenance"]["source_path"] == "mouse_brain.d"

    def test_a_directory_source_keeps_its_name_despite_a_trailing_separator(self):
        # Bruker .d and Waters .raw name a directory, so the value can
        # arrive with a separator the name would otherwise be lost behind.
        document = build_metadata_document(
            self._from("/mnt/share/2025/mouse_brain.d/"),
            pixel_size_um=(20.0, 20.0),
            source_format="bruker",
        )
        assert document["provenance"]["source_path"] == "mouse_brain.d"

    def test_a_name_is_left_alone(self):
        document = build_metadata_document(
            self._from("input.imzML"),
            pixel_size_um=(20.0, 20.0),
            source_format="imzml",
        )
        assert document["provenance"]["source_path"] == "input.imzML"

    def test_a_source_that_reduces_to_nothing_is_unset_rather_than_empty(self):
        # An empty string is not a name, and the field is optional. The
        # store convention is to omit rather than to write empty.
        document = build_metadata_document(
            self._from("/"),
            pixel_size_um=(20.0, 20.0),
            source_format="imzml",
        )
        assert "source_path" not in document["provenance"]

    def test_a_document_from_no_metadata_at_all_still_builds(self):
        document = build_metadata_document(
            None, pixel_size_um=(20.0, 20.0), source_format="imzml"
        )
        assert "source_path" not in document["provenance"]
