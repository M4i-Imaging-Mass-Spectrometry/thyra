"""METASPACE export: a straight mapping that never fabricates data."""

from thyra.metadata.schema import MSIMetadata, to_metaspace


def _full_meta() -> MSIMetadata:
    return MSIMetadata.model_validate(
        {
            "schema_version": "0.1.0",
            "sample": {
                "organism": "Mus musculus",
                "organism_part": "liver",
                "condition": "wildtype",
            },
            "preparation": {
                "sample_stabilisation": "fresh frozen",
                "matrix": "2,5-dihydroxybenzoic acid (DHB)",
                "matrix_application": "ImagePrep",
                "solvent": "methanol",
            },
            "ms_analysis": {
                "polarity": "positive",
                "ionisation_source": "MALDI",
                "analyzer": "Orbitrap",
                "detector_resolving_power": {"value": 130000, "at_mz": 400},
                "pixel_size_um": {"x": 20, "y": 40},
            },
            "provenance": {"thyra_version": "1.0.0"},
        }
    )


class TestToMetaspace:
    def test_complete_document_maps_without_warnings(self):
        document, warnings = to_metaspace(_full_meta())
        assert warnings == []
        assert document == {
            "Data_Type": "Imaging MS",
            "Sample_Information": {
                "Organism": "Mus musculus",
                "Organism_Part": "liver",
                "Condition": "wildtype",
                "Sample_Growth_Conditions": "",
            },
            "Sample_Preparation": {
                "Sample_Stabilisation": "fresh frozen",
                "Tissue_Modification": "",
                "MALDI_Matrix": "2,5-dihydroxybenzoic acid (DHB)",
                "MALDI_Matrix_Application": "ImagePrep",
                "Solvent": "methanol",
            },
            "MS_Analysis": {
                "Polarity": "Positive",
                "Ionisation_Source": "MALDI",
                "Analyzer": "Orbitrap",
                "Pixel_Size": {"Xaxis": 20, "Yaxis": 40},
                "Detector_Resolving_Power": {"mz": 400, "Resolving_Power": 130000},
            },
        }

    def test_missing_required_fields_are_empty_and_warned(self):
        meta = _full_meta().model_copy(deep=True)
        meta.sample.organism = None
        meta.ms_analysis.polarity = None
        document, warnings = to_metaspace(meta)
        assert document["Sample_Information"]["Organism"] == ""
        assert document["MS_Analysis"]["Polarity"] == ""
        joined = " ".join(warnings)
        assert "Sample_Information.Organism" in joined
        assert "MS_Analysis.Polarity" in joined

    def test_matrix_free_source_truthfully_reports_no_matrix(self):
        meta = _full_meta().model_copy(deep=True)
        meta.preparation.matrix = None
        meta.preparation.matrix_application = None
        meta.ms_analysis.ionisation_source = "DESI"
        document, warnings = to_metaspace(meta)
        assert document["Sample_Preparation"]["MALDI_Matrix"] == "none"
        assert document["Sample_Preparation"]["MALDI_Matrix_Application"] == "none"
        assert not any("MALDI_Matrix" in w for w in warnings)

    def test_missing_matrix_on_maldi_is_warned_not_invented(self):
        meta = _full_meta().model_copy(deep=True)
        meta.preparation.matrix = None
        document, warnings = to_metaspace(meta)
        assert document["Sample_Preparation"]["MALDI_Matrix"] == ""
        assert any("MALDI_Matrix" in w for w in warnings)

    def test_missing_resolving_power_is_omitted_and_warned(self):
        meta = _full_meta().model_copy(deep=True)
        meta.ms_analysis.detector_resolving_power = None
        document, warnings = to_metaspace(meta)
        assert "Detector_Resolving_Power" not in document["MS_Analysis"]
        assert any("Detector_Resolving_Power" in w for w in warnings)

    def test_fractional_pixel_size_stays_a_float(self):
        meta = _full_meta().model_copy(deep=True)
        meta.ms_analysis.pixel_size_um.x = 12.5
        document, _ = to_metaspace(meta)
        assert document["MS_Analysis"]["Pixel_Size"]["Xaxis"] == 12.5
