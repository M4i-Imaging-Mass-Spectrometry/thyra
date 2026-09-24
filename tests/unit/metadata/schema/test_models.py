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


class TestCvBindings:
    def test_bound_fields_are_exactly_the_expected_set(self):
        from thyra.metadata.schema.models import field_cv_bindings

        accessions = {cv["accession"] for cv in field_cv_bindings().values()}
        assert accessions == {
            "IMS:1000046",  # pixel size (x)
            "IMS:1000047",  # pixel size y
            "MS:1000465",  # scan polarity
            "MS:1000008",  # ionization type
            "MS:1000443",  # mass analyzer type
            "MS:1000031",  # instrument model
            "MS:1001269",  # instrument vendor
            "MS:1000529",  # instrument serial number
            "MS:1000800",  # mass resolving power
            "MS:1002892",  # ion mobility attribute
            "MS:1000511",  # ms level
            "MS:1000827",  # isolation window target m/z
            "MS:1000828",  # isolation window lower offset
            "MS:1000829",  # isolation window upper offset
            "MS:1000045",  # collision energy
            "IMS:1006000",  # repetition rate
            "IMS:1006001",  # laser shots per spectrum
            "IMS:1006008",  # optical image location
            "IMS:1006017",  # method used to align optical image
        }

    def test_every_binding_resolves_in_the_local_tables(self):
        from thyra.metadata.ontology.cache import ONTOLOGY
        from thyra.metadata.schema.models import field_cv_bindings

        for path, cv in field_cv_bindings().items():
            entry = ONTOLOGY.terms.get(cv["accession"])
            assert entry is not None, f"{path}: unknown {cv['accession']}"
            assert (
                entry[0] == cv["name"]
            ), f"{path}: bound name {cv['name']!r} vs CV label {entry[0]!r}"

    def test_bindings_land_in_the_committed_json_schema(self):
        # The claim must be checkable from the artifact alone, without
        # importing Thyra.
        import json
        from pathlib import Path

        from thyra.metadata.schema import models

        artifact = (
            Path(models.__file__).parent / models.SCHEMA_JSON_FILENAME
        ).read_text(encoding="utf-8")
        rendered = json.dumps(json.loads(artifact))
        assert '"IMS:1000046"' in rendered
        assert '"MS:1000443"' in rendered

    def test_candidate_concepts_are_declared(self):
        from thyra.metadata.schema import models

        concepts = [c for c, _ in models.CANDIDATE_CV_CONCEPTS]
        assert any("resampling" in c for c in concepts)
        assert any("stage offset" in c for c in concepts)
        # The acquisition section's fields without a CV term are listed
        # with the ones it does not bind, so the gap is on record.
        paths = [p for _, p in models.CANDIDATE_CV_CONCEPTS]
        assert "acquisition.acquisition_datetime" in paths
        assert "acquisition.laser_power_percent" in paths
        assert "acquisition.method_file" in paths
        # PSI-MS has one calibration term, MS:1001485, and it names the
        # processing action; every fact the calibration section states is
        # a gap on record instead of a binding.
        for field in (
            "calibration_datetime",
            "recalibrated",
            "software",
            "mz_standard_deviation_ppm",
            "lock_mass_corrected",
        ):
            assert f"calibration.{field}" in paths
        assert "ms_analysis.n_spectra" in paths
        # IMS names the optical image and the alignment method, and has no
        # term for the teaching points a flexImaging registration rests on.
        assert "alignment.teaching_points" in paths

    def test_the_calibration_step_term_resolves_locally(self):
        from thyra.metadata.ontology.cache import ONTOLOGY

        assert ONTOLOGY.terms["MS:1001485"][0] == "m/z calibration"


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


class TestAcquisition:
    def _document(self, **acquisition) -> dict:
        doc = _minimal().model_dump()
        doc["acquisition"] = acquisition
        return doc

    def test_the_section_validates(self):
        from thyra.metadata.schema import Acquisition

        section = Acquisition(
            acquisition_datetime="2025-04-22T08:59:34.395+02:00",
            laser_power_percent=70.0,
            laser_frequency_hz=1000.0,
            shots_per_pixel=50,
            method_file="imaging.m",
        )
        meta = MSIMetadata.model_validate(self._document(**section.model_dump()))
        assert meta.acquisition == section

    def test_an_absent_section_is_absent_in_the_uns_dict(self):
        assert _minimal().acquisition is None
        assert "acquisition" not in _minimal().to_uns_dict()

    def test_extra_keys_are_forbidden(self):
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(self._document(laser_power=70.0))

    @pytest.mark.parametrize(
        "value",
        [
            "07-Nov-2019 14:09:44",  # MassLynx, unparsed
            "06/23/2026 21:12:35",  # SmartSoft-TOF, unparsed
            "2019-11-07",  # a date is not a datetime
            "2019-11-07T25:00:00",  # not a time
            "not a timestamp",
        ],
    )
    def test_an_unparseable_datetime_is_rejected(self, value):
        with pytest.raises(ValidationError, match="acquisition_datetime"):
            MSIMetadata.model_validate(self._document(acquisition_datetime=value))

    @pytest.mark.parametrize(
        "value",
        [
            "2019-11-07T14:09:44",
            "2019-11-07T14:09:44.123+02:00",
            "2019-11-07T14:09:44Z",
        ],
    )
    def test_iso_8601_with_or_without_an_offset_is_accepted(self, value):
        meta = MSIMetadata.model_validate(self._document(acquisition_datetime=value))
        assert meta.acquisition is not None
        assert meta.acquisition.acquisition_datetime == value

    def test_the_method_must_be_a_name_not_a_path(self):
        for path in ("D:\\Methods\\imaging.m", "/data/imaging.m"):
            with pytest.raises(ValidationError, match="not a path"):
                MSIMetadata.model_validate(self._document(method_file=path))

    def test_laser_power_is_bounded_to_a_percentage(self):
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(self._document(laser_power_percent=140.0))
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(self._document(laser_frequency_hz=0.0))
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(self._document(shots_per_pixel=0))


class TestCalibration:
    def _document(self, **calibration) -> dict:
        doc = _minimal().model_dump()
        doc["calibration"] = calibration
        return doc

    def test_the_section_validates(self):
        from thyra.metadata.schema import Calibration

        section = Calibration(
            calibration_datetime="2026-07-27T17:29:40",
            recalibrated=True,
            original_calibration_datetime="2025-04-22T08:43:30+02:00",
            software="timsTOF",
            software_version="4.1.12",
            n_reference_peaks=14,
            mz_standard_deviation_ppm=0.411709,
            lock_mass_corrected=False,
        )
        meta = MSIMetadata.model_validate(self._document(**section.model_dump()))
        assert meta.calibration == section

    def test_an_absent_section_is_absent_in_the_uns_dict(self):
        assert _minimal().calibration is None
        assert "calibration" not in _minimal().to_uns_dict()

    def test_extra_keys_are_forbidden(self):
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(self._document(calibration_user="someone"))

    @pytest.mark.parametrize(
        "field", ["calibration_datetime", "original_calibration_datetime"]
    )
    @pytest.mark.parametrize(
        "value", ["08/15/19 11:50", "2019-08-15", "2019-08-15T25:00", "not a time"]
    )
    def test_a_timestamp_must_be_an_iso_8601_date_and_time(self, field, value):
        fields = {field: value}
        if field == "original_calibration_datetime":
            fields["recalibrated"] = True
        with pytest.raises(ValidationError, match=field):
            MSIMetadata.model_validate(self._document(**fields))

    def test_minute_precision_is_a_valid_timestamp(self):
        # MassLynx records its calibration time to the minute; the value
        # keeps that precision rather than gaining seconds nobody recorded.
        meta = MSIMetadata.model_validate(
            self._document(calibration_datetime="2019-08-15T11:50")
        )
        assert meta.calibration is not None
        assert meta.calibration.calibration_datetime == "2019-08-15T11:50"

    def test_an_original_time_needs_a_recalibration_to_be_the_original_of(self):
        with pytest.raises(ValidationError, match="recalibrated is not true"):
            MSIMetadata.model_validate(
                self._document(original_calibration_datetime="2025-04-22T08:43:30")
            )
        with pytest.raises(ValidationError, match="recalibrated is not true"):
            MSIMetadata.model_validate(
                self._document(
                    recalibrated=False,
                    original_calibration_datetime="2025-04-22T08:43:30",
                )
            )

    def test_a_standard_deviation_needs_two_peaks_to_rest_on(self):
        for peaks in ({}, {"n_reference_peaks": 1}):
            with pytest.raises(ValidationError, match="at least 2"):
                MSIMetadata.model_validate(
                    self._document(mz_standard_deviation_ppm=0.4, **peaks)
                )
        meta = MSIMetadata.model_validate(
            self._document(mz_standard_deviation_ppm=0.4, n_reference_peaks=2)
        )
        assert meta.calibration is not None

    def test_a_zero_standard_deviation_is_not_a_value(self):
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(
                self._document(mz_standard_deviation_ppm=0.0, n_reference_peaks=4)
            )

    def test_the_peak_count_stands_on_its_own(self):
        meta = MSIMetadata.model_validate(self._document(n_reference_peaks=1))
        assert meta.calibration is not None
        assert meta.calibration.n_reference_peaks == 1
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(self._document(n_reference_peaks=0))

    def test_a_0_7_0_document_still_validates(self):
        # 0.8.0 only adds optional fields, so a document written by the
        # previous minor version is valid as it stands.
        doc = _minimal().to_uns_dict()
        doc["schema_version"] = "0.7.0"
        doc["ms_analysis"]["manufacturer"] = "Bruker"
        assert MSIMetadata.model_validate(doc).calibration is None


# The first teaching point of a real flexImaging sequence file. Its stage x is
# negative: stage positions are signed.
_POINT = {
    "image_x_px": 4780.0,
    "image_y_px": 784.0,
    "stage_x_um": -26352.0,
    "stage_y_um": 26386.0,
}


class TestAlignment:
    def _document(self, **alignment) -> dict:
        doc = _minimal().model_dump()
        doc["alignment"] = alignment
        return doc

    def test_the_section_validates(self):
        from thyra.metadata.schema import Alignment, TeachingPoint

        section = Alignment(
            optical_image_file="slide_0000.tif",
            method="teaching points",
            teaching_points=[
                TeachingPoint(**_POINT),
                TeachingPoint(
                    image_x_px=13648,
                    image_y_px=11296,
                    stage_x_um=-8793,
                    stage_y_um=5388,
                ),
                TeachingPoint(
                    image_x_px=32156,
                    image_y_px=724,
                    stage_x_um=28003,
                    stage_y_um=26505,
                ),
            ],
        )
        meta = MSIMetadata.model_validate(self._document(**section.model_dump()))
        assert meta.alignment == section

    def test_an_absent_section_is_absent_in_the_uns_dict(self):
        assert _minimal().alignment is None
        assert "alignment" not in _minimal().to_uns_dict()

    def test_extra_keys_are_forbidden(self):
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(self._document(areas=[]))
        with pytest.raises(ValidationError):
            MSIMetadata.model_validate(
                self._document(teaching_points=[{**_POINT, "name": "01"}])
            )

    def test_a_teaching_point_is_both_of_its_positions(self):
        for missing in _POINT:
            point = {k: v for k, v in _POINT.items() if k != missing}
            with pytest.raises(ValidationError, match=missing):
                MSIMetadata.model_validate(self._document(teaching_points=[point]))

    def test_the_image_must_be_a_name_not_a_path(self):
        for path in ("D:\\Data\\slide_0000.tif", "../scans/slide_0000.tif"):
            with pytest.raises(ValidationError, match="not a path"):
                MSIMetadata.model_validate(self._document(optical_image_file=path))

    def test_the_teaching_points_are_packed_in_the_uns_dict(self):
        # A list of objects does not survive AnnData/zarr, so it travels as
        # JSON, exactly as the isolation windows do.
        import json

        meta = MSIMetadata.model_validate(
            self._document(optical_image_file="slide.tif", teaching_points=[_POINT])
        )
        stored = meta.to_uns_dict()["alignment"]
        assert json.loads(stored["teaching_points"]) == [_POINT]
        assert stored["optical_image_file"] == "slide.tif"

    def test_a_0_8_0_document_still_validates(self):
        doc = _minimal().to_uns_dict()
        doc["schema_version"] = "0.8.0"
        doc["calibration"] = {"lock_mass_corrected": False}
        assert MSIMetadata.model_validate(doc).alignment is None


class TestSpectrumCount:
    def test_it_is_optional_and_positive(self):
        assert _minimal().ms_analysis.n_spectra is None
        assert "n_spectra" not in _minimal().to_uns_dict()["ms_analysis"]
        assert MSAnalysis(pixel_size_um=PixelSizeUm(x=1, y=1), n_spectra=713)
        with pytest.raises(ValidationError):
            MSAnalysis(pixel_size_um=PixelSizeUm(x=1, y=1), n_spectra=0)


class TestProcessingStepTerm:
    def test_a_step_may_name_its_psi_ms_action(self):
        from thyra.metadata.schema import ProcessingStep, SoftwareRef

        step = ProcessingStep(
            name="m/z calibration",
            action_term=OntologyTerm(accession="MS:1001485", name="m/z calibration"),
            software=SoftwareRef(name="thyra", version="1.0.0"),
            parameters={"use_recalibrated_state": True},
        )
        meta = _minimal().model_copy(deep=True)
        meta.processing = [step]
        restored = MSIMetadata.model_validate(meta.model_dump())
        assert restored.processing[0].action_term == step.action_term

    def test_the_term_is_optional(self):
        from thyra.metadata.schema import ProcessingStep, SoftwareRef

        step = ProcessingStep(
            name="conversion", software=SoftwareRef(name="thyra", version="1")
        )
        assert step.action_term is None
        assert "action_term" not in step.model_dump(exclude_none=True)
