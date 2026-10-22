"""The instrument comes from the configuration the converted spectra use (D81).

An imzML spectrum uses the configuration its scan names, else the run's
default; an mzPeak spectrum the one its scan's ``instrument_configuration_id``
names, else ``run.default_instrument_id``. A configuration nothing uses
describes no stored spectrum, and configurations the spectra use that state
different facts leave those facts unset. An imzML and its archive resolve
alike, as D29 intends.
"""

import logging
from pathlib import Path

import pytest
from pyimzml.ImzMLWriter import ImzMLWriter

import thyra
from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra.readers.imzml.imzml_reader import ImzMLReader
from thyra.readers.mzpeak import MzPeakReader


def _config(cid: str, model: str, extra: str = "") -> str:
    return (
        f'<instrumentConfiguration id="{cid}"><cvParam cvRef="MS" '
        f'accession="MS:1000031" name="instrument model" value="{model}"/>'
        f"{extra}</instrumentConfiguration>"
    )


ORBITRAP = _config("ICorbi", "Q Exactive")
TOF = _config("ICtof", "rapifleX")
ONE_CONFIGURATION = (
    '<instrumentConfigurationList count="1">\n'
    '    <instrumentConfiguration id="IC1">\n'
    "    </instrumentConfiguration>"
)


def _write_imzml(path: Path, configurations: str, scan_refs, run_default="ICtof"):
    """Four spectra; ``scan_refs`` names each scan's configuration."""
    with ImzMLWriter(str(path), mode="processed") as writer:
        for i in range(4):
            writer.addSpectrum(
                [400.0, 500.0, 600.0], [1.0, 2.0, 3.0], (i % 2 + 1, i // 2 + 1, 1)
            )
    text = path.read_text(encoding="ISO-8859-1")
    text = text.replace(
        ONE_CONFIGURATION,
        f'<instrumentConfigurationList count="2">{configurations}',
    )
    text = text.replace(
        'defaultInstrumentConfigurationRef="IC1"',
        f'defaultInstrumentConfigurationRef="{run_default}"',
    )
    parts = text.split('instrumentConfigurationRef="instrumentConfiguration0"')
    assert len(parts) == len(scan_refs) + 1
    text = parts[0] + "".join(
        f'instrumentConfigurationRef="{ref}"' + part
        for ref, part in zip(scan_refs, parts[1:])
    )
    path.write_text(text, encoding="ISO-8859-1")
    return path


def _conversion_instrument(path: Path):
    with ImzMLReader(path) as reader:
        return reader.get_comprehensive_metadata().instrument_info


def _document_model(path: Path):
    return thyra.read_metadata_document(path)["ms_analysis"].get("instrument_model")


class TestImzML:
    @pytest.mark.parametrize("configurations", [ORBITRAP + TOF, TOF + ORBITRAP])
    def test_the_scans_configuration_decides_in_either_order(
        self, tmp_path, configurations
    ):
        """F175: the list's order no longer picks the model."""
        path = _write_imzml(tmp_path / "two.imzML", configurations, ["ICtof"] * 4)

        assert _conversion_instrument(path)["instrument_model"] == "rapifleX"
        assert _document_model(path) == "rapifleX"
        assert thyra.preview_msi(path).instrument_type == (
            thyra.preview_msi(
                _write_imzml(tmp_path / "tof.imzML", TOF + TOF, ["ICtof"] * 4)
            ).instrument_type
        )

    def test_the_run_default_decides_for_a_scan_naming_none(self, tmp_path):
        """pyimzml's writer names a configuration the list does not hold."""
        path = _write_imzml(
            tmp_path / "default.imzML",
            ORBITRAP + TOF,
            ["instrumentConfiguration0"] * 4,
            run_default="ICtof",
        )

        assert _conversion_instrument(path)["instrument_model"] == "rapifleX"
        assert _document_model(path) == "rapifleX"

    def test_spectra_using_configurations_that_disagree_record_no_model(
        self, tmp_path, thyra_logs
    ):
        path = _write_imzml(
            tmp_path / "split.imzML",
            ORBITRAP + TOF,
            ["ICtof", "ICorbi", "ICtof", "ICorbi"],
        )

        with thyra_logs("thyra.metadata", logging.WARNING) as records:
            info = _conversion_instrument(path)

        assert "instrument_model" not in info
        assert any("instrument model" in r.getMessage() for r in records)

    def test_a_pyimzml_file_keeps_its_one_configuration(self, tmp_path):
        """The dangling scan reference falls back to the run default IC1."""
        path = tmp_path / "plain.imzML"
        with ImzMLWriter(str(path), mode="processed") as writer:
            for i in range(2):
                writer.addSpectrum([400.0, 500.0], [1.0, 2.0], (i + 1, 1, 1))
        text = path.read_text(encoding="ISO-8859-1").replace(
            '<instrumentConfiguration id="IC1">',
            '<instrumentConfiguration id="IC1"><cvParam cvRef="MS" '
            'accession="MS:1000031" name="instrument model" value="rapifleX"/>',
        )
        path.write_text(text, encoding="ISO-8859-1")

        assert _conversion_instrument(path)["instrument_model"] == "rapifleX"

    def test_the_solvent_of_the_used_source_is_read(self, tmp_path):
        """F72: IMS:1001211 on the ion source reaches preparation.solvent."""
        source = (
            '<componentList count="1"><source order="1"><cvParam cvRef="IMS" '
            'accession="IMS:1001211" name="solvent" value="methanol"/>'
            "</source></componentList>"
        )
        path = _write_imzml(
            tmp_path / "desi.imzML",
            _config("ICorbi", "Q Exactive") + _config("ICtof", "rapifleX", source),
            ["ICtof"] * 4,
        )

        with ImzMLReader(path) as reader:
            acquisition = reader.get_comprehensive_metadata().acquisition_params
        assert acquisition["solvent"] == "methanol"
        document = thyra.read_metadata_document(path)
        assert document["preparation"] == {"solvent": "methanol"}


def _parameter(name, accession, value=None):
    return {"name": name, "accession": accession, "value": value, "unit": None}


def _mzpeak_config(cid, model, serial=None, analyzer=None):
    parameters = [_parameter("instrument model", "MS:1000031", model)]
    if serial is not None:
        parameters.append(_parameter("instrument serial number", "MS:1000529", serial))
    components = []
    if analyzer is not None:
        components.append(
            {
                "component_type": "analyzer",
                "order": 1,
                "parameters": [_parameter(*analyzer)],
            }
        )
    return {"id": cid, "components": components, "parameters": parameters}


def _archive_info(path: Path):
    with MzPeakReader(path) as reader:
        return reader.get_comprehensive_metadata().instrument_info


ORBITRAP_CONFIG = _mzpeak_config(0, "Q Exactive", analyzer=("orbitrap", "MS:1000484"))
TOF_CONFIG = _mzpeak_config(1, "rapifleX", analyzer=("time-of-flight", "MS:1000084"))


class TestMzPeak:
    def test_the_scans_configuration_id_decides(self, tmp_path):
        """F175: the configuration every scan names, not the first listed."""
        archive = build_mzpeak(
            tmp_path / "scans.mzpeak",
            grid_spectra(2, 2),
            file_metadata={
                "instrument_configuration_list": [ORBITRAP_CONFIG, TOF_CONFIG]
            },
            scan_columns={"instrument_configuration_id": [1, 1, 1, 1]},
        )

        info = _archive_info(archive)

        assert info["instrument_model"] == "rapifleX"
        assert info["analyzer"] == "time-of-flight"
        assert "instrument_type" not in info

    def test_the_run_default_decides_without_a_scan_column(self, tmp_path):
        archive = build_mzpeak(
            tmp_path / "default.mzpeak",
            grid_spectra(2, 2),
            file_metadata={
                "instrument_configuration_list": [ORBITRAP_CONFIG, TOF_CONFIG],
                "run": {"id": "r", "start_time": None, "default_instrument_id": 1},
            },
        )

        assert _archive_info(archive)["instrument_model"] == "rapifleX"

    def test_scans_naming_configurations_that_disagree_record_no_model(self, tmp_path):
        archive = build_mzpeak(
            tmp_path / "split.mzpeak",
            grid_spectra(2, 2),
            file_metadata={
                "instrument_configuration_list": [ORBITRAP_CONFIG, TOF_CONFIG]
            },
            scan_columns={"instrument_configuration_id": [0, 1, 0, 1]},
        )

        info = _archive_info(archive)

        assert "instrument_model" not in info
        assert "instrument_type" not in info

    @pytest.mark.parametrize(
        "serial, expected",
        [
            ("1877407.00228", "1877407.00228"),
            (1877407.00228, "1877407.00228"),
            ("0", "0"),
            (0, "0"),
            ("SN-ABC12", "SN-ABC12"),
        ],
    )
    def test_a_serial_number_is_kept_as_text(self, tmp_path, serial, expected):
        """F216: a serial stated as a JSON number is not dropped."""
        archive = build_mzpeak(
            tmp_path / "serial.mzpeak",
            grid_spectra(2, 2),
            file_metadata={
                "instrument_configuration_list": [
                    _mzpeak_config(0, "timsTOF Maldi 2", serial=serial)
                ]
            },
        )

        assert _archive_info(archive)["instrument_serial_number"] == expected

    @pytest.mark.parametrize(
        "accession, name, model, family",
        [
            # The generic term with no value names no model, as on imzML.
            ("MS:1000031", "Q Exactive HF", None, None),
            # A term the model table lacks: its ontology name, and the
            # family its name identifies.
            ("MS:1003245", "Q Exactive UHMR", "Q Exactive UHMR", "Orbitrap"),
        ],
    )
    def test_valueless_model_terms_resolve_as_on_imzml(
        self, tmp_path, accession, name, model, family
    ):
        """F202: a valueless term reads as True; a known accession takes its name."""
        configuration = {
            "id": 0,
            "components": [],
            "parameters": [_parameter(name, accession)],
        }
        archive = build_mzpeak(
            tmp_path / "model.mzpeak",
            grid_spectra(3, 2),
            file_metadata={"instrument_configuration_list": [configuration]},
        )

        info = _archive_info(archive)

        assert info.get("instrument_model") == model
        assert info.get("instrument_type") == family

    def test_the_solvent_of_the_used_source_is_read(self, tmp_path):
        """F72, on the archive route."""
        configuration = _mzpeak_config(0, "Q Exactive")
        configuration["components"] = [
            {
                "component_type": "ionsource",
                "order": 1,
                "parameters": [_parameter("solvent", "IMS:1001211", "methanol")],
            }
        ]
        archive = build_mzpeak(
            tmp_path / "desi.mzpeak",
            grid_spectra(2, 2),
            file_metadata={"instrument_configuration_list": [configuration]},
        )

        with MzPeakReader(archive) as reader:
            acquisition = reader.get_comprehensive_metadata().acquisition_params
        assert acquisition["solvent"] == "methanol"
