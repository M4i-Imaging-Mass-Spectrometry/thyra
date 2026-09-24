"""The block a converter writes must read back and validate.

The uns parity suite (tests/unit/converters/test_uns_provenance_parity.py)
asserts the block is identical across all four write paths; this file
covers the consumer side: reading it out of a real store, validating
it, and exporting METASPACE metadata from it.
"""

import json

import pytest
from click.testing import CliRunner

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.metadata.schema import (
    MSI_METADATA_SCHEMA_VERSION,
    read_msi_metadata_blocks,
    validate_document,
)


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    from thyra.converters.spatialdata import SpatialDataConverter

    output = tmp_path_factory.mktemp("schema_store") / "out.zarr"
    converter = SpatialDataConverter(
        reader=MockMSIReader(
            MockMSIConfig(n_x=4, n_y=4, n_mz_bins=200, peaks_per_spectrum=(10, 20))
        ),
        output_path=output,
        dataset_id="mock",
        pixel_size_um=10.0,
    )
    assert converter.convert() is True
    return output


class TestStoreRoundTrip:
    def test_block_reads_back_and_validates(self, store):
        blocks = read_msi_metadata_blocks(store)
        assert blocks, "converted store carries no msi_metadata block"
        for block in blocks.values():
            assert block["schema_version"] == MSI_METADATA_SCHEMA_VERSION
            assert block["ms_analysis"]["pixel_size_um"] == {"x": 10.0, "y": 10.0}
            meta, issues = validate_document(block)
            assert meta is not None
            assert not [i for i in issues if i.severity == "error"]

    def test_processing_history_records_the_conversion(self, store):
        from thyra import __version__

        for block in read_msi_metadata_blocks(store).values():
            steps = block["processing"]
            assert steps[0]["name"] == "conversion"
            assert steps[0]["software"] == {
                "name": "thyra",
                "version": __version__,
            }
            # No resampling was configured, so no resampling step.
            assert [s["name"] for s in steps] == ["conversion"]

    def test_root_attrs_carry_the_explicit_affine(self, store):
        import zarr

        cs_global = dict(zarr.open_group(str(store), mode="r").attrs)[
            "coordinate_systems"
        ]["global"]
        assert cs_global["raster_to_global_affine"] == [
            [10.0, 0.0, 0.0],
            [0.0, 10.0, 0.0],
            [0.0, 0.0, 1.0],
        ]

    def test_validate_cli_passes_on_a_real_store(self, store):
        from thyra.metadata.schema.cli import validate_command

        result = CliRunner().invoke(validate_command, [str(store)])
        assert result.exit_code == 0, result.output

    def test_the_mock_reader_reports_no_acquisition_facts_so_the_section_is_absent(
        self, store
    ):
        for block in read_msi_metadata_blocks(store).values():
            assert "acquisition" not in block


class _ReaderWithAcquisitionFacts(MockMSIReader):
    """The mock reader with the acquisition keys a Bruker tsf/tdf extractor writes."""

    def _create_metadata_extractor(self):
        from dataclasses import replace

        extractor = super()._create_metadata_extractor()
        implementation = extractor._extract_comprehensive_impl

        def with_acquisition_params():
            return replace(
                implementation(),
                acquisition_params={
                    "acquisition_datetime": "2025-04-22T08:59:34.395+02:00",
                    "laser_power": 70.0,
                    "laser_frequency": 1000.0,
                    "num_laser_shots": 50,
                    "method_name": "imaging_pos.m",
                },
            )

        extractor._extract_comprehensive_impl = with_acquisition_params
        return extractor


@pytest.fixture(scope="module")
def store_with_acquisition(tmp_path_factory):
    from thyra.converters.spatialdata import SpatialDataConverter

    output = tmp_path_factory.mktemp("schema_store_acq") / "out.zarr"
    converter = SpatialDataConverter(
        reader=_ReaderWithAcquisitionFacts(
            MockMSIConfig(n_x=4, n_y=4, n_mz_bins=200, peaks_per_spectrum=(10, 20))
        ),
        output_path=output,
        dataset_id="mock",
        pixel_size_um=10.0,
    )
    assert converter.convert() is True
    return output


class TestAcquisitionSectionRoundTrip:
    """The section survives the store as plain scalars and validates back."""

    def test_the_section_reads_back_as_written(self, store_with_acquisition):
        blocks = read_msi_metadata_blocks(store_with_acquisition)
        assert blocks
        for block in blocks.values():
            # No laser_power_percent: the converter driven directly reports
            # no source format, and a percentage is only trusted for the
            # formats whose unit was verified (see builder.py).
            assert block["acquisition"] == {
                "acquisition_datetime": "2025-04-22T08:59:34.395+02:00",
                "laser_frequency_hz": 1000.0,
                "shots_per_pixel": 50,
                "method_file": "imaging_pos.m",
            }
            # Plain strings stay plain strings: the JSON encoding the store
            # applies to string *lists* must not have touched them.
            assert isinstance(block["acquisition"]["acquisition_datetime"], str)
            assert isinstance(block["acquisition"]["method_file"], str)
            meta, issues = validate_document(block)
            assert meta is not None and meta.acquisition is not None
            assert not [i for i in issues if i.severity == "error"]

    def test_the_uns_group_holds_scalars_not_encoded_strings(
        self, store_with_acquisition
    ):
        import spatialdata as sd

        from thyra.metadata import sanitize_uns_string_arrays

        sdata = sd.read_zarr(store_with_acquisition)
        for table in sdata.tables.values():
            uns = sanitize_uns_string_arrays(table.uns)
            section = uns["msi_metadata"]["acquisition"]
            assert section["acquisition_datetime"] == "2025-04-22T08:59:34.395+02:00"
            assert section["shots_per_pixel"] == 50
            # Nothing in the section is a JSON string to decode.
            assert not any(
                isinstance(v, str) and v.startswith(("[", "{"))
                for v in section.values()
            )

    def test_validate_cli_passes_with_the_section_present(self, store_with_acquisition):
        from thyra.metadata.schema.cli import validate_command

        result = CliRunner().invoke(validate_command, [str(store_with_acquisition)])
        assert result.exit_code == 0, result.output

    def test_export_metaspace_cli_works_on_a_real_store(self, store, tmp_path):
        from thyra.metadata.schema.cli import export_metaspace_command

        output = tmp_path / "submission.json"
        result = CliRunner().invoke(
            export_metaspace_command, [str(store), "-o", str(output)]
        )
        assert result.exit_code == 0, result.output
        document = json.loads(output.read_text(encoding="utf-8"))
        assert document["MS_Analysis"]["Pixel_Size"] == {"Xaxis": 10, "Yaxis": 10}

    def test_resampling_step_serialises_the_config(self, tmp_path):
        from thyra.converters.spatialdata import SpatialDataConverter

        converter = SpatialDataConverter(
            reader=MockMSIReader(
                MockMSIConfig(n_x=4, n_y=4, n_mz_bins=200, peaks_per_spectrum=(10, 20))
            ),
            output_path=tmp_path / "out.zarr",
            dataset_id="mock",
            pixel_size_um=10.0,
            resampling_config={"method": "nearest_neighbor", "target_bins": 100},
        )
        steps = converter.uns._processing_provenance(converter._uns_context())
        assert [s.name for s in steps] == ["conversion", "mass axis resampling"]
        parameters = steps[1].parameters
        assert parameters["method"] == "nearest_neighbor"
        assert parameters["target_bins"] == 100
        # Unset config fields are dropped, not serialised as nulls.
        assert "min_mz" not in parameters

    def test_resolved_resampling_lands_in_the_stored_step(self, tmp_path):
        from thyra.converters.spatialdata import SpatialDataConverter

        output = tmp_path / "resampled.zarr"
        converter = SpatialDataConverter(
            reader=MockMSIReader(
                MockMSIConfig(n_x=4, n_y=4, n_mz_bins=200, peaks_per_spectrum=(10, 20))
            ),
            output_path=output,
            dataset_id="mock",
            pixel_size_um=10.0,
            resampling_config={"method": "nearest_neighbor", "target_bins": 100},
        )
        assert converter.convert() is True

        block = read_msi_metadata_blocks(output)["mock_z0"]
        step = next(
            s for s in block["processing"] if s["name"] == "mass axis resampling"
        )
        # The step declares what was DONE: the resolved axis, not just
        # the requested config.
        parameters = step["parameters"]
        assert parameters["method"] == "nearest_neighbor"
        assert parameters["target_bins"] == 100
        assert parameters["axis_type"]
        assert parameters["min_mz"] < parameters["max_mz"]

    def test_a_zarr_group_that_is_not_spatialdata_is_reported_cleanly(self, tmp_path):
        import zarr

        zarr.open_group(str(tmp_path / "plain.zarr"), mode="a")
        with pytest.raises(ValueError, match="tables"):
            read_msi_metadata_blocks(tmp_path / "plain.zarr")


class _ReaderWithCalibrationFacts(MockMSIReader):
    """The mock reader with the calibration keys a PHI extractor writes.

    PHI because it exercises everything the section can hold that a store
    must round-trip: a boolean, a timestamp, an integer and a float. The
    reader also reports the calibration it applied, which becomes the
    ``m/z calibration`` processing step.
    """

    def get_applied_mz_calibration(self):
        return {"calibration": "appended"}

    def _create_metadata_extractor(self):
        from dataclasses import replace

        extractor = super()._create_metadata_extractor()
        implementation = extractor._extract_comprehensive_impl

        def with_calibration():
            return replace(
                implementation(),
                raw_metadata={
                    "calibration": {
                        "source": "appended",
                        "recalibrated": True,
                        "recalibration_date": "07/27/2026 17:29:40",
                        "recalibration_calibrants": json.dumps(
                            [
                                {"measured_mz": 26.003016, "theoretical_mz": 26.003099},
                                {"measured_mz": 41.998143, "theoretical_mz": 41.998001},
                                {"measured_mz": 57.975196, "theoretical_mz": 57.975201},
                                {"measured_mz": 117.971046, "theoretical_mz": 117.9711},
                            ]
                        ),
                    }
                },
            )

        extractor._extract_comprehensive_impl = with_calibration
        return extractor


@pytest.fixture(scope="module")
def store_with_calibration(tmp_path_factory):
    from thyra.converters.spatialdata import SpatialDataConverter

    output = tmp_path_factory.mktemp("schema_store_cal") / "out.zarr"
    converter = SpatialDataConverter(
        reader=_ReaderWithCalibrationFacts(
            MockMSIConfig(n_x=4, n_y=4, n_mz_bins=200, peaks_per_spectrum=(10, 20))
        ),
        output_path=output,
        dataset_id="mock",
        pixel_size_um=10.0,
    )
    assert converter.convert() is True
    return output


class TestCalibrationSectionRoundTrip:
    """The section and the step survive the store and validate back."""

    def test_the_section_reads_back_as_written(self, store_with_calibration):
        blocks = read_msi_metadata_blocks(store_with_calibration)
        assert blocks
        for block in blocks.values():
            assert block["calibration"] == {
                "calibration_datetime": "2026-07-27T17:29:40",
                "recalibrated": True,
                "n_reference_peaks": 4,
                "mz_standard_deviation_ppm": 2.69798,
            }
            meta, issues = validate_document(block)
            assert meta is not None and meta.calibration is not None
            assert meta.calibration.recalibrated is True
            assert not [i for i in issues if i.severity == "error"]

    def test_the_applied_calibration_is_a_step_bound_to_its_term(
        self, store_with_calibration
    ):
        for block in read_msi_metadata_blocks(store_with_calibration).values():
            steps = block["processing"]
            assert [s["name"] for s in steps] == ["conversion", "m/z calibration"]
            assert steps[1]["action_term"] == {
                "accession": "MS:1001485",
                "name": "m/z calibration",
            }
            assert steps[1]["parameters"] == {"calibration": "appended"}

    def test_the_mock_reader_counts_its_spectra(self, store_with_calibration):
        for block in read_msi_metadata_blocks(store_with_calibration).values():
            assert block["ms_analysis"]["n_spectra"] == 16

    def test_validate_cli_passes_with_the_section_present(self, store_with_calibration):
        from thyra.metadata.schema.cli import validate_command

        result = CliRunner().invoke(validate_command, [str(store_with_calibration)])
        assert result.exit_code == 0, result.output

    def test_a_store_without_the_facts_has_no_section_and_no_step(self, store):
        for block in read_msi_metadata_blocks(store).values():
            assert "calibration" not in block
            assert "m/z calibration" not in [s["name"] for s in block["processing"]]
