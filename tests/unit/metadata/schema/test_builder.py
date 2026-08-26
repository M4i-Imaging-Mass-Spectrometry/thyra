"""Auto-population: the builder reports what the source knows, nothing more."""

from thyra.metadata.schema import build_msi_metadata
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata


def _essential(source_path: str = "input.imzML") -> EssentialMetadata:
    return EssentialMetadata(
        dimensions=(10, 10, 1),
        coordinate_bounds=(0.0, 9.0, 0.0, 9.0),
        mass_range=(100.0, 1000.0),
        pixel_size=(20.0, 25.0),
        n_spectra=100,
        total_peaks=1000,
        estimated_memory_gb=0.1,
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
