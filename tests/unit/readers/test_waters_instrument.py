"""Which Waters instrument wrote a run, and what the reader does about it.

The three cases the detection has to get right are the ones that exist on
disk: a SELECT SERIES MRT run (``OpticMode = MRT``, ``Resolution`` ~225,000,
``$$ Instrument: MRT#``), a Synapt G2-Si run (no ``OpticMode``, no
instrument line, ``Resolution`` 10,000, tab-run key padding), and a directory
with neither side file at all. The field values and formatting below are
copied from those two runs.
"""

from __future__ import annotations

import logging
from unittest.mock import MagicMock, patch

import pytest

from thyra.readers.waters.imaging_grid import ImagingGrid
from thyra.readers.waters.instrument import (
    MRT_RESOLUTION_THRESHOLD,
    WatersInstrument,
    identify_waters_instrument,
    parse_extern_inf,
    parse_function_parameters,
    parse_header_txt,
    read_cyclic_sequence,
    read_stage_pattern,
)
from thyra.readers.waters.masslynx_lib import FunctionType, ScanInfoData
from thyra.readers.waters.waters_reader import WatersReader

# The MRT pads keys with spaces then one tab; the Synapt with runs of tabs.
MRT_EXTERN = (
    "Instrument Configuration:\n"
    "Lteff                    \t48400.0\n"
    "Veff                    \t7221.664\n"
    "Resolution              \t225909.053\n"
    "Acquisition Device      \tWatersADC\n"
    "ADC Sample Frequency (GHz)\t1.35\n"
    "\n"
    "Function Parameters - Function 1 - TOF MS FUNCTION\n"
    "ADC Sample Frequency (GHz)                   \t1.35\n"
    "ADC Pusher Period (us)                       \t187.0\n"
    "OpticMode                                    \tMRT\n"
    "Polarity                                     \tPositive\n"
)
MRT_HEADER = (
    "$$ Version: 01.00\n"
    "$$ Acquired Name: 20260821_MRT_FTW2\n"
    "$$ Instrument: MRT#\n"
    "$$ Cal Function 1: 0.000000000000000E+00,1.000000000000000E+00,T1\n"
)
SYNAPT_EXTERN = (
    "Parameters for E:\\MassLynx.PRO\\Acqudb\\synapt on MALDI MS.exp\n"
    "Created by 4.1 SCN957\n"
    " \n"
    "Instrument Configuration:\n"
    "Lteff\t\t\t\t\t\t1800.0\n"
    "Veff\t\t\t\t\t\t7200.50\n"
    "Resolution\t\t\t\t\t10000\n"
    "Min Points in Peak\t\t\t\t2\n"
    "Acquisition Device\t\t\t\tWatersADC\n"
    "ADC Sample Frequency (GHz)\t\t\t3.0\n"
    "Pusher Cycle Time (\xb5s)\t\t\t\tAutomatic\n"
    "PusherInterval\t\t\t\t\t54.000000\n"
)
SYNAPT_HEADER = (
    "$$ Version: 01.00\n"
    "$$ Acquired Name: 20201216_MDA468_slide20200630_synapt\n"
    "$$ Job Code: 20201216_MDA468_slide20200630_synapt\n"
)

# A Cyclic IMS full-scan run: the field set and formatting of a real one,
# cut to what the reader reads. No "Acquisition Algorithm" line anywhere.
CYCLIC_EXTERN = (
    "Instrument Configuration:\n"
    "Lteff                    \t2190.0\n"
    "Veff                    \t6925.44\n"
    "Resolution              \t30000.00\n"
    "Acquisition Device\tWatersADC\n"
    "ADC Sample Frequency (GHz)\t7.2\n"
    "\n"
    "Experimental Instrument Parameters\n"
    "Quad.MSFixedMass.Setting                                         \t500.0\n"
    "TofADC.IMSCycleStartDelay.Setting                                \t4682\n"
    "\n"
    "Function Parameters - Function 1 - TOF HDMS FUNCTION\n"
    "ADC Sample Frequency (GHz)                   \t7.2\n"
    "ADC Pusher Period (us)                       \t66.0\n"
    "IMS Cycle Start Delay (ms)                   \t4682\n"
    "Pushes Per IMS Increment                     \t3\n"
    "MSMSMode                                     \tMS\n"
    "OpticMode                                    \tV-Mode\n"
    "TOFMode                                      \tIMS\n"
    "CCS Calibration                              \tNo\n"
    "[COLLISION ENERGY]\n"
    "Using Auto Trap MS Collision Energy (eV)     \t6.0000\n"
)
CYCLIC_HEADER = (
    "$$ Version: 04.20\n"
    "$$ Acquired Name: cyclic_run_01\n"
    "$$ Instrument: Cyclic-IMS\n"
    "#GBA001\n"
    "\n"
)
#: The two function sections of the public Synapt XS image
#: (MSV000092638), cut to a few of their fields.
SYNAPT_XS_EXTERN = (
    "Instrument Configuration:\n"
    "Acquisition Device\t\t\t\tWatersADC\n"
    "Acquisition Algorithm\t\t\t\tADC Mode\n"
    "ADC Pushes Per IMS Increment\t\t\t1\n"
    "Function Parameters - Function 1 - MALDI MOBILITY TOF MS FUNCTION\n"
    "[MS]\n"
    "Start Mass\t\t\t\t\t50.0\n"
    "ADC Pusher Frequency (\xb5s)\t\t\t66.0\n"
    "Function Parameters - Function 2 - MALDI MOBILITY TOF MS FUNCTION\n"
    "[MS]\n"
    "Start Mass\t\t\t\t\t50.0\n"
    "ADC Pusher Frequency (\xb5s)\t\t\t66.0\n"
)
#: What a Synapt G1 run writes in _header.txt: a serial, no model.
G1_HEADER = "$$ Version: 01.00\n$$ Instrument: UAA232\n"
CYCLIC_METHOD_XML = """<?xml version="1.0"?>
<MsMethod InstrumentType="QTof" InstrumentModel="Cyclic IMS" Version="1.0">
    <Settings>
        <Setting Name="PushesPerBin" Value="3"/>
        <Table Name="CyclicSequence" Filename="LPC 20 passes" IMSCycleTime="348.58">
            <Row>
                <Setting Name="CyclicFunction" Value="0"/>
                <Setting Name="CyclicArrayMode" Value="0"/>
                <Setting Name="Time" Value="10"/>
            </Row>
            <Row>
                <Setting Name="CyclicFunction" Value="1"/>
                <Setting Name="CyclicArrayMode" Value="3"/>
                <Setting Name="Time" Value="298.98"/>
            </Row>
            <Row>
                <Setting Name="CyclicFunction" Value="6"/>
                <Setting Name="CyclicArrayMode" Value="7"/>
            </Row>
        </Table>
    </Settings>
</MsMethod>
"""


def _raw_dir(tmp_path, extern=None, header=None):
    raw = tmp_path / "run.raw"
    raw.mkdir(parents=True)
    (raw / "_FUNC001.DAT").write_bytes(b"\x00" * 64)
    if extern is not None:
        (raw / "_extern.inf").write_text(extern, encoding="latin-1")
    if header is not None:
        (raw / "_header.txt").write_text(header, encoding="latin-1")
    return raw


class TestParsing:
    def test_mrt_formatting(self):
        fields = parse_extern_inf(MRT_EXTERN)
        assert fields["OpticMode"] == "MRT"
        assert fields["Resolution"] == "225909.053"
        assert fields["ADC Pusher Period (us)"] == "187.0"

    def test_synapt_tab_runs_and_latin1(self):
        fields = parse_extern_inf(SYNAPT_EXTERN)
        assert fields["Resolution"] == "10000"
        assert fields["Lteff"] == "1800.0"
        assert fields["Pusher Cycle Time (\xb5s)"] == "Automatic"
        assert "OpticMode" not in fields

    def test_first_occurrence_wins(self):
        fields = parse_extern_inf("K\t1\nK\t2\n")
        assert fields["K"] == "1"

    def test_section_headers_are_skipped(self):
        assert "Instrument Configuration:" not in parse_extern_inf(MRT_EXTERN)

    def test_header(self):
        fields = parse_header_txt(MRT_HEADER)
        assert fields["Instrument"] == "MRT#"
        assert fields["Acquired Name"] == "20260821_MRT_FTW2"
        assert "Instrument" not in parse_header_txt(SYNAPT_HEADER)


class TestDecision:
    def test_mrt_run(self, tmp_path):
        inst = identify_waters_instrument(_raw_dir(tmp_path, MRT_EXTERN, MRT_HEADER))
        assert inst.is_mrt
        assert inst.decided_by == "OpticMode=MRT"
        assert inst.name == "SELECT SERIES MRT"
        assert inst.resolution == pytest.approx(225909.053)

    def test_synapt_run(self, tmp_path):
        inst = identify_waters_instrument(
            _raw_dir(tmp_path, SYNAPT_EXTERN, SYNAPT_HEADER)
        )
        assert not inst.is_mrt
        assert inst.decided_by == "Resolution=10000"
        assert inst.optic_mode is None
        assert inst.instrument_label is None

    def test_neither_field(self, tmp_path):
        inst = identify_waters_instrument(_raw_dir(tmp_path))
        assert not inst.is_mrt
        assert inst.decided_by == "no instrument field found"
        assert inst.resolution is None

    def test_optic_mode_outranks_resolution(self, tmp_path):
        # OpticMode says a non-MRT optic even though the resolution is high.
        extern = MRT_EXTERN.replace("\tMRT\n", "\tV\n")
        inst = identify_waters_instrument(_raw_dir(tmp_path, extern))
        assert not inst.is_mrt
        assert inst.decided_by == "OpticMode=V"

    def test_header_label_when_optic_mode_is_absent(self, tmp_path):
        inst = identify_waters_instrument(_raw_dir(tmp_path, header=MRT_HEADER))
        assert inst.is_mrt
        assert inst.decided_by == "$$ Instrument=MRT#"

    def test_resolution_fallback(self, tmp_path):
        extern = "Resolution\t225909.053\n"
        inst = identify_waters_instrument(_raw_dir(tmp_path, extern))
        assert inst.is_mrt
        assert inst.decided_by == "Resolution=225909"
        assert MRT_RESOLUTION_THRESHOLD == 100_000.0

    def test_unparseable_resolution_is_not_a_decision(self, tmp_path):
        inst = identify_waters_instrument(_raw_dir(tmp_path, "Resolution\tn/a\n"))
        assert not inst.is_mrt
        assert inst.decided_by == "no instrument field found"

    def test_decision_is_logged(self, tmp_path, thyra_logs):
        with thyra_logs("thyra.readers.waters.instrument", logging.INFO) as records:
            identify_waters_instrument(_raw_dir(tmp_path, MRT_EXTERN, MRT_HEADER))
        assert "SELECT SERIES MRT (OpticMode=MRT)" in records.text


class TestSampleSpacing:
    """Predicted from Lteff, Veff and the ADC clock; measured values beside."""

    def test_mrt(self, tmp_path):
        inst = identify_waters_instrument(_raw_dir(tmp_path, MRT_EXTERN))
        # 1.16 mDa measured at m/z 990-1000, 1.04 at 795-805.
        assert inst.profile_sample_spacing_da(1000.0) == pytest.approx(
            1.143e-3, rel=0.01
        )
        assert inst.profile_sample_spacing_da(800.0) == pytest.approx(
            1.022e-3, rel=0.01
        )

    def test_synapt(self, tmp_path):
        inst = identify_waters_instrument(_raw_dir(tmp_path, SYNAPT_EXTERN))
        # 13.4 mDa measured at m/z 900-1000.
        assert inst.profile_sample_spacing_da(1000.0) == pytest.approx(
            13.8e-3, rel=0.01
        )

    def test_sqrt_law(self, tmp_path):
        inst = identify_waters_instrument(_raw_dir(tmp_path, MRT_EXTERN))
        ratio = inst.profile_sample_spacing_da(400.0) / inst.profile_sample_spacing_da(
            1600.0
        )
        assert ratio == pytest.approx(0.5)

    def test_missing_constants(self, tmp_path):
        inst = identify_waters_instrument(_raw_dir(tmp_path))
        assert inst.profile_sample_spacing_da(1000.0) is None
        assert WatersInstrument(is_mrt=False, decided_by="x", flight_path_mm=1.0)
        assert (
            WatersInstrument(
                is_mrt=False,
                decided_by="x",
                flight_path_mm=0.0,
                effective_voltage_v=1.0,
                adc_sample_frequency_ghz=1.0,
            ).profile_sample_spacing_da(500.0)
            is None
        )


def _grid_one_pixel():
    scan_map = {
        (0, 0): ScanInfoData(
            ms_level=1,
            polarity=0,
            drift_scan_count=0,
            is_profile=1,
            precursor_mz=0.0,
            rt=1.0,
            laser_x_pos=0.1,
            laser_y_pos=0.1,
        )
    }
    return ImagingGrid(
        x_index_map={100.0: 0},
        y_index_map={100.0: 0},
        pixel_count_x=1,
        pixel_count_y=1,
        pixel_size_x=0.0,
        pixel_size_y=0.0,
        lateral_width=0.0,
        lateral_height=0.0,
        scan_map=scan_map,
    )


class TestReaderDefault:
    """``use_centroid=None`` follows the instrument; an explicit value wins."""

    def test_mrt_defaults_to_profile(self, tmp_path):
        reader = WatersReader(_raw_dir(tmp_path, MRT_EXTERN, MRT_HEADER))
        assert reader.use_centroid is False
        assert reader.instrument.is_mrt

    def test_synapt_defaults_to_centroid(self, tmp_path):
        reader = WatersReader(_raw_dir(tmp_path, SYNAPT_EXTERN, SYNAPT_HEADER))
        assert reader.use_centroid is True

    def test_unknown_defaults_to_centroid(self, tmp_path):
        reader = WatersReader(_raw_dir(tmp_path))
        assert reader.use_centroid is True
        assert reader.instrument.decided_by == "no instrument field found"

    @pytest.mark.parametrize("explicit", [True, False])
    def test_explicit_choice_overrides_either_way(self, tmp_path, explicit):
        mrt = WatersReader(
            _raw_dir(tmp_path / "a", MRT_EXTERN, MRT_HEADER), use_centroid=explicit
        )
        synapt = WatersReader(
            _raw_dir(tmp_path / "b", SYNAPT_EXTERN, SYNAPT_HEADER),
            use_centroid=explicit,
        )
        assert mrt.use_centroid is explicit
        assert synapt.use_centroid is explicit

    def test_default_is_logged_with_its_reason(self, tmp_path, thyra_logs):
        with thyra_logs("thyra.readers.waters.waters_reader", logging.INFO) as records:
            WatersReader(_raw_dir(tmp_path, MRT_EXTERN, MRT_HEADER))
        assert "profile trace by default (OpticMode=MRT" in records.text

    @patch("thyra.readers.waters.waters_reader.MassLynxLib")
    @patch("thyra.readers.waters.waters_reader.build_imaging_grid")
    def test_mrt_default_reaches_the_native_mode_switch(
        self, mock_build_grid, mock_ml_cls, tmp_path
    ):
        """The decision is only real once ``set_centroid`` sees it."""
        mock_ml = MagicMock()
        mock_ml.is_ion_mobility_function.return_value = False
        mock_ml_cls.get_instance.return_value = mock_ml
        mock_ml.open_file.return_value = "handle"
        mock_ml.is_imaging_file.return_value = True
        mock_ml.get_number_of_functions.return_value = 1
        mock_ml.classify_function.return_value = FunctionType.MS
        mock_build_grid.return_value = _grid_one_pixel()

        reader = WatersReader(_raw_dir(tmp_path, MRT_EXTERN, MRT_HEADER))
        reader._ensure_initialized()
        mock_ml.set_centroid.assert_called_once_with("handle", False)

        extractor = reader._create_metadata_extractor()
        assert extractor._use_centroid is False
        assert extractor._instrument is reader.instrument


class TestMobilityReadable:
    """Whether the library can read the run's drift bins, from _extern.inf.

    The two shapes are copied from real runs: a Synapt G2-Si / XS declares
    the acquisition algorithm under the instrument configuration; an older
    Synapt G1 run carries a ``Pusher Interval`` field and no such line.
    """

    def test_adc_mode_is_readable(self, tmp_path):
        extern = SYNAPT_EXTERN + "Acquisition Algorithm\t\t\t\tADC Mode\n"
        decision = identify_waters_instrument(_raw_dir(tmp_path, extern, None))

        assert decision.acquisition_algorithm == "ADC Mode"
        assert decision.mobility_readable is True

    def test_an_old_run_without_the_field_is_not(self, tmp_path):
        extern = "Instrument Configuration:\nPusher Interval (\xb5S)\t90\n"
        decision = identify_waters_instrument(_raw_dir(tmp_path, extern, None))

        assert decision.acquisition_algorithm is None
        assert decision.mobility_readable is False

    def test_no_side_file_is_not(self, tmp_path):
        decision = identify_waters_instrument(_raw_dir(tmp_path, None, None))

        assert decision.mobility_readable is False

    def test_a_cyclic_run_is_readable(self, tmp_path):
        decision = identify_waters_instrument(
            _raw_dir(tmp_path, CYCLIC_EXTERN, CYCLIC_HEADER)
        )

        assert decision.acquisition_algorithm is None
        assert decision.acquisition_device == "WatersADC"
        assert decision.is_cyclic_ims
        assert decision.mobility_readable is True
        assert decision.name == "SELECT SERIES Cyclic IMS"
        assert not decision.is_mrt

    def test_the_header_label_alone_is_not_enough(self, tmp_path):
        extern = CYCLIC_EXTERN.replace("Acquisition Device\tWatersADC\n", "")
        decision = identify_waters_instrument(_raw_dir(tmp_path, extern, CYCLIC_HEADER))

        assert decision.is_cyclic_ims
        assert decision.mobility_readable is False

    @pytest.mark.parametrize(
        "extern, header",
        [
            (SYNAPT_EXTERN, SYNAPT_HEADER),
            (MRT_EXTERN, MRT_HEADER),
            (CYCLIC_EXTERN, G1_HEADER),
        ],
        ids=["synapt", "mrt", "g1-header"],
    )
    def test_the_adc_alone_is_not_enough(self, tmp_path, extern, header):
        # A Synapt G2-Si and an MRT declare WatersADC too, and no algorithm.
        decision = identify_waters_instrument(_raw_dir(tmp_path, extern, header))

        assert decision.acquisition_device == "WatersADC"
        assert decision.mobility_readable is False

    def test_a_g1_run_stays_unreadable(self, tmp_path):
        extern = "Instrument Configuration:\nPusher Interval (\xb5S)\t90\n"
        decision = identify_waters_instrument(_raw_dir(tmp_path, extern, G1_HEADER))

        assert not decision.is_cyclic_ims
        assert decision.mobility_readable is False

    def test_the_header_is_found_in_capitals(self, tmp_path):
        raw = _raw_dir(tmp_path, CYCLIC_EXTERN)
        (raw / "_HEADER.TXT").write_text(CYCLIC_HEADER, encoding="latin-1")

        assert identify_waters_instrument(raw).mobility_readable is True


class TestFunctionParameters:
    """Each function's own section of _extern.inf (D72)."""

    def test_sections_are_keyed_from_zero(self):
        sections = parse_function_parameters(CYCLIC_EXTERN)

        assert list(sections) == [0]
        assert sections[0]["MSMSMode"] == "MS"
        assert sections[0]["Pushes Per IMS Increment"] == "3"
        # The sub-header carries no tab; the line after it is still read.
        assert sections[0]["Using Auto Trap MS Collision Energy (eV)"] == "6.0000"
        # The instrument-wide setting is not a function parameter.
        assert "TofADC.IMSCycleStartDelay.Setting" not in sections[0]

    def test_each_function_keeps_its_own_values(self):
        text = (
            "Function Parameters - Function 1 - TOF MS FUNCTION\n"
            "MSMSMode\tMS\nMSMSMode\tSIM\n"
            "Function Parameters - Function 2 - TOF MSMS FUNCTION\n"
            "MSMSMode\tMSMS\n"
        )
        sections = parse_function_parameters(text)

        assert sections[0]["MSMSMode"] == "MS"
        assert sections[1]["MSMSMode"] == "MSMS"
        # The whole-file parser keeps only the first function's value.
        assert parse_extern_inf(text)["MSMSMode"] == "MS"

    def test_a_synapt_run_has_none(self):
        assert parse_function_parameters(SYNAPT_EXTERN) == {}
        # A real Synapt XS section: its own fields, none of the Cyclic ones.
        sections = parse_function_parameters(SYNAPT_XS_EXTERN)
        assert sections[0]["ADC Pusher Frequency (\xb5s)"] == "66.0"
        assert sections[1]["Start Mass"] == "50.0"

    def test_msms_mode_and_bin_timing(self, tmp_path):
        decision = identify_waters_instrument(
            _raw_dir(tmp_path, CYCLIC_EXTERN, CYCLIC_HEADER)
        )

        assert decision.msms_mode(0) == "MS"
        assert decision.msms_mode(1) is None
        assert decision.drift_bin_timing(0) == (4682, 3)
        assert decision.drift_bin_timing(1) is None

    def test_a_synapt_run_declares_no_bin_timing(self, tmp_path):
        decision = identify_waters_instrument(_raw_dir(tmp_path, SYNAPT_XS_EXTERN))

        assert decision.function_parameters
        for function in (0, 1):
            assert decision.drift_bin_timing(function) is None
            assert decision.msms_mode(function) is None

    @pytest.mark.parametrize(
        "fields",
        [
            "Pushes Per IMS Increment\t2\n",
            "IMS Cycle Start Delay (ms)\t4682\n",
            "IMS Cycle Start Delay (ms)\t309.002\nPushes Per IMS Increment\t3\n",
            "IMS Cycle Start Delay (ms)\t4682\nPushes Per IMS Increment\t0\n",
            "IMS Cycle Start Delay (ms)\t-1\nPushes Per IMS Increment\t3\n",
        ],
        ids=["no-delay", "no-pushes", "delay-in-ms", "zero-pushes", "negative"],
    )
    def test_a_timing_that_places_nothing_is_refused(self, tmp_path, fields):
        extern = "Function Parameters - Function 1 - TOF HDMS FUNCTION\n" + fields
        decision = identify_waters_instrument(_raw_dir(tmp_path, extern))

        with pytest.raises(ValueError, match="do not place its drift bins"):
            decision.drift_bin_timing(0)


class TestCyclicSequence:
    """The separation sequence from the run's method file, as declared."""

    def test_the_sequence_is_read(self, tmp_path):
        raw = _raw_dir(tmp_path, CYCLIC_EXTERN, CYCLIC_HEADER)
        (raw / "cyclic_run_01.xml").write_text(CYCLIC_METHOD_XML, encoding="utf-8")

        sequence = read_cyclic_sequence(raw)

        assert sequence == {
            "sequence_name": "LPC 20 passes",
            "ims_cycle_time_ms": pytest.approx(348.58),
            "steps": [
                {"function": 0, "array_mode": 0, "time_ms": 10.0},
                {"function": 1, "array_mode": 3, "time_ms": pytest.approx(298.98)},
                {"function": 6, "array_mode": 7},
            ],
        }

    def test_no_method_file_gives_none(self, tmp_path):
        assert read_cyclic_sequence(_raw_dir(tmp_path, CYCLIC_EXTERN)) is None

    def test_a_method_without_a_sequence_gives_none(self, tmp_path):
        raw = _raw_dir(tmp_path)
        (raw / "run.xml").write_text(
            '<?xml version="1.0"?><MsMethod><Settings/></MsMethod>', encoding="utf-8"
        )
        (raw / "broken.xml").write_text("<MsMethod", encoding="utf-8")

        assert read_cyclic_sequence(raw) is None

    def test_an_entity_declaration_is_refused_not_expanded(self, tmp_path):
        raw = _raw_dir(tmp_path)
        (raw / "run.xml").write_text(
            '<?xml version="1.0"?><!DOCTYPE m [<!ENTITY n "x">]>'
            '<MsMethod><Settings><Table Name="CyclicSequence" Filename="&n;"/>'
            "</Settings></MsMethod>",
            encoding="utf-8",
        )

        assert read_cyclic_sequence(raw) is None


def _method_xml(settings="", table=""):
    return (
        '<?xml version="1.0"?>\n<MsMethod><Settings>'
        f"{settings}{table}</Settings></MsMethod>"
    )


def _region_row(x, y, dx, dy):
    return (
        "<Row>"
        f'<Setting Name="MaldiXStart" Value="{x}"/>'
        f'<Setting Name="MaldiYStart" Value="{y}"/>'
        f'<Setting Name="MaldiXLength" Value="{dx}"/>'
        f'<Setting Name="MaldiYLength" Value="{dy}"/>'
        "</Row>"
    )


STEPS = (
    '<Setting Name="MaldiXStep" Value="0.0201"/>'
    '<Setting Name="MaldiYStep" Value="0.0200"/>'
)
TWO_REGIONS = (
    '<Table Name="MultiRegions">'
    + _region_row("66.5000", "74.1000", "0.5", "0.01")
    + _region_row("66.6700", "63.9400", "0.5", "0.01")
    + "</Table>"
)


class TestStagePattern:
    """The pixel steps and regions a run's method declares."""

    def test_two_regions_and_the_steps_are_read_in_um(self, tmp_path):
        raw = _raw_dir(tmp_path)
        (raw / "run.xml").write_text(_method_xml(STEPS, TWO_REGIONS), encoding="utf-8")

        pattern = read_stage_pattern(raw)

        assert pattern.x_step_um == pytest.approx(20.1)
        assert pattern.y_step_um == pytest.approx(20.0)
        assert [r.number for r in pattern.regions] == [1, 2]
        first, second = pattern.regions
        assert (first.x_start_um, first.y_start_um) == pytest.approx((66500.0, 74100.0))
        assert (first.x_length_um, first.y_length_um) == pytest.approx((500.0, 10.0))
        assert (second.x_start_um, second.y_start_um) == pytest.approx(
            (66670.0, 63940.0)
        )

    def test_a_single_pattern_run_is_not_a_region(self, tmp_path):
        """Top-level start and length settings describe one pattern."""
        raw = _raw_dir(tmp_path)
        settings = (
            STEPS
            + '<Setting Name="MaldiXStart" Value="66.5" '
            + 'Mapping="Maldi.Stage.Pattern.XStart"/>'
            + '<Setting Name="MaldiXLength" Value="0.5"/>'
        )
        (raw / "run.xml").write_text(_method_xml(settings), encoding="utf-8")

        pattern = read_stage_pattern(raw)

        assert pattern.regions == ()
        assert pattern.x_step_um == pytest.approx(20.1)

    def test_a_method_with_neither_gives_none(self, tmp_path):
        raw = _raw_dir(tmp_path)
        (raw / "run.xml").write_text(_method_xml(), encoding="utf-8")

        assert read_stage_pattern(raw) is None

    def test_no_method_file_gives_none(self, tmp_path):
        assert read_stage_pattern(_raw_dir(tmp_path)) is None

    def test_a_missing_directory_gives_none(self, tmp_path):
        assert read_stage_pattern(tmp_path / "gone.raw") is None

    @pytest.mark.parametrize("value", ["0", "-0.02", "nan", "inf", "abc", ""])
    def test_an_unusable_step_is_none(self, tmp_path, value):
        raw = _raw_dir(tmp_path)
        settings = f'<Setting Name="MaldiXStep" Value="{value}"/>'
        (raw / "run.xml").write_text(
            _method_xml(settings, TWO_REGIONS), encoding="utf-8"
        )

        pattern = read_stage_pattern(raw)

        assert pattern.x_step_um is None
        assert len(pattern.regions) == 2

    def test_one_bad_row_drops_the_whole_table(self, tmp_path, thyra_logs):
        raw = _raw_dir(tmp_path)
        table = (
            '<Table Name="MultiRegions">'
            + _region_row("66.5", "74.1", "0.5", "0.01")
            + _region_row("66.67", "oops", "0.5", "0.01")
            + "</Table>"
        )
        (raw / "run.xml").write_text(_method_xml(STEPS, table), encoding="utf-8")

        with thyra_logs("thyra.readers.waters", logging.WARNING) as records:
            pattern = read_stage_pattern(raw)

        assert pattern.regions == ()
        assert pattern.x_step_um == pytest.approx(20.1)
        assert any("Row 2" in r.getMessage() for r in records)

    def test_the_file_named_after_the_run_comes_first(self, tmp_path):
        raw = _raw_dir(tmp_path)
        (raw / "a_other.xml").write_text(
            _method_xml('<Setting Name="MaldiXStep" Value="0.5"/>'), encoding="utf-8"
        )
        (raw / "run.xml").write_text(_method_xml(STEPS), encoding="utf-8")

        assert read_stage_pattern(raw).x_step_um == pytest.approx(20.1)

    def test_an_entity_declaration_is_refused_not_expanded(self, tmp_path):
        raw = _raw_dir(tmp_path)
        (raw / "run.xml").write_text(
            '<?xml version="1.0"?><!DOCTYPE m [<!ENTITY n "0.02">]>'
            '<MsMethod><Settings><Setting Name="MaldiXStep" Value="&n;"/>'
            "</Settings></MsMethod>",
            encoding="utf-8",
        )

        assert read_stage_pattern(raw) is None
