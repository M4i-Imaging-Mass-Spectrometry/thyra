"""Structural validation: the pydantic models enforce the schema's shape."""

import pytest
from pydantic import ValidationError

from thyra.metadata.schema import (
    MSI_METADATA_SCHEMA_VERSION,
    MSAnalysis,
    MSIMetadata,
    OntologyTerm,
    PixelSizeUm,
    Provenance,
)


def _minimal() -> MSIMetadata:
    return MSIMetadata(
        ms_analysis=MSAnalysis(pixel_size_um=PixelSizeUm(x=20.0, y=20.0)),
        provenance=Provenance(thyra_version="1.0.0"),
    )


class TestMSIMetadata:
    def test_minimal_document_is_valid(self):
        meta = _minimal()
        assert meta.schema_version == MSI_METADATA_SCHEMA_VERSION

    def test_round_trips_through_dump_and_validate(self):
        meta = _minimal()
        assert MSIMetadata.model_validate(meta.model_dump()) == meta

    def test_pixel_size_is_required(self):
        with pytest.raises(ValidationError):
            MSAnalysis()

    def test_pixel_size_must_be_positive(self):
        with pytest.raises(ValidationError):
            PixelSizeUm(x=0.0, y=20.0)

    def test_unknown_fields_are_rejected(self):
        doc = _minimal().model_dump()
        doc["ms_analysis"]["lazer_power"] = 3.0
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(doc)

    def test_accession_must_be_a_curie(self):
        with pytest.raises(ValidationError):
            OntologyTerm(accession="not a curie", name="x")

    def test_polarity_vocabulary_is_closed(self):
        doc = _minimal().model_dump()
        doc["ms_analysis"]["polarity"] = "both"
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(doc)

    def test_contradicting_polarity_term_is_rejected(self):
        with pytest.raises(ValidationError, match="contradicts"):
            MSAnalysis(
                pixel_size_um=PixelSizeUm(x=20.0, y=20.0),
                polarity="positive",
                polarity_term=OntologyTerm(
                    accession="MS:1000129", name="negative scan"
                ),
            )

    def test_agreeing_polarity_term_is_accepted(self):
        analysis = MSAnalysis(
            pixel_size_um=PixelSizeUm(x=20.0, y=20.0),
            polarity="negative",
            polarity_term=OntologyTerm(accession="MS:1000129", name="negative scan"),
        )
        assert analysis.polarity == "negative"


class TestProcessing:
    def test_processing_steps_round_trip(self):
        from thyra.metadata.schema import ProcessingStep, SoftwareRef

        meta = _minimal().model_copy(deep=True)
        meta.processing = [
            ProcessingStep(
                name="mass axis resampling",
                software=SoftwareRef(name="thyra", version="1.0.0"),
                parameters={"target_bins": 50000, "method": "nearest_neighbor"},
            )
        ]
        restored = MSIMetadata.model_validate(meta.model_dump())
        assert restored.processing[0].parameters["target_bins"] == 50000

    def test_empty_processing_is_omitted_from_uns(self):
        assert "processing" not in _minimal().to_uns_dict()

    def test_step_requires_software(self):
        from thyra.metadata.schema import ProcessingStep

        with pytest.raises(ValidationError):
            ProcessingStep(name="normalisation")


class TestToUnsDict:
    def test_none_fields_are_dropped(self):
        data = _minimal().to_uns_dict()
        assert "polarity" not in data["ms_analysis"]
        assert "source_format" not in data["provenance"]

    def test_empty_user_sections_are_omitted_not_written_empty(self):
        data = _minimal().to_uns_dict()
        assert "sample" not in data
        assert "preparation" not in data

    def test_populated_user_sections_are_kept(self):
        meta = _minimal().model_copy(deep=True)
        meta.sample.organism = "Mus musculus"
        data = meta.to_uns_dict()
        assert data["sample"] == {"organism": "Mus musculus"}

    def test_uns_dict_validates_back(self):
        data = _minimal().to_uns_dict()
        assert MSIMetadata.model_validate(data) is not None
