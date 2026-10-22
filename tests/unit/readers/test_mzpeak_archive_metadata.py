"""What an mzPeak archive states beyond its spectra, and what Thyra does with it.

Each class is one finding of issue #422, made on the public example archives
and on archives from mzpeak-convert 0.16.0: the instrument, where the image
sits, its regions, spectra of other MS levels, ion mobility and the images
an archive embeds. The file-level blocks are spelled as that converter
writes them.
"""

from __future__ import annotations

import json
import logging
import uuid
import zipfile

import numpy as np
import pytest

from tests.fixtures.mzpeak_builder import INDEX_MEMBER, SCANS_MEMBER, Spectrum
from tests.fixtures.mzpeak_builder import bruker_region as _region
from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from tests.fixtures.mzpeak_builder import region_parameter as _region_parameter
from thyra.core.base_reader import OpticalImageLabel
from thyra.errors import ConversionRefused
from thyra.readers.mzpeak import MzPeakReader
from thyra.readers.mzpeak.mzpeak_reader import _on_tdf_scale
from thyra.resampling.data_characteristics import DataCharacteristics
from thyra.resampling.instrument_detectors import InstrumentDetectorChain


def _parameter(name, accession, value=None, unit=None):
    """One CV parameter, as the archive writes it."""
    return {"name": name, "accession": accession, "value": value, "unit": unit}


#: The instrument block of the public bladder example (LTQ Orbitrap).
ORBITRAP_CONFIGURATION = [
    {
        "components": [
            {
                "component_type": "ionsource",
                "order": 1,
                "parameters": [
                    _parameter("electrospray ionization", "MS:1000073"),
                ],
            },
            {
                "component_type": "analyzer",
                "order": 2,
                "parameters": [
                    _parameter("orbitrap", "MS:1000484"),
                    _parameter("accuracy", "MS:1000014", 0, "MS:1000040"),
                ],
            },
        ],
        "parameters": [
            _parameter("LTQ Orbitrap Discovery", "MS:1000555"),
            _parameter("instrument serial number", "MS:1000529", "none"),
        ],
        "software_reference": "Xcalibur",
        "id": 0,
    }
]

#: The instrument block of a Bruker TSF archive from mzpeak-convert 0.16.0.
TIMSTOF_CONFIGURATION = [
    {
        "components": [
            {
                "component_type": "analyzer",
                "order": 1,
                "parameters": [_parameter("time-of-flight", "MS:1000084")],
            }
        ],
        "parameters": [
            _parameter("Bruker Daltonics timsTOF series", "MS:1003123"),
            _parameter("instrument model", "MS:1000031", "timsTOF Maldi 2"),
            _parameter("instrument serial number", "MS:1000529", "1877407.00228"),
        ],
        "software_reference": "timsTOF",
        "id": 0,
    }
]


def _comprehensive(path):
    """Read comprehensive metadata and close the archive."""
    with MzPeakReader(path) as reader:
        return reader.get_comprehensive_metadata()


def _detected(comprehensive):
    """The detector the resampling chain picks for this metadata."""
    metadata = {
        "essential_metadata": {
            "spectrum_type": comprehensive.essential.spectrum_type,
            "total_peaks": comprehensive.essential.total_peaks,
            "n_spectra": comprehensive.essential.n_spectra,
        },
        "instrument_info": comprehensive.instrument_info,
        "format_specific": comprehensive.format_specific,
    }
    return InstrumentDetectorChain().detect(DataCharacteristics.from_metadata(metadata))


class TestInstrument:
    """The instrument is read as the imzML route reads it."""

    def test_an_orbitrap_archive_gets_the_orbitrap_axis(self, tmp_path):
        """Without the family, an Orbitrap archive was binned as unknown."""
        archive = build_mzpeak(
            tmp_path / "orbitrap.mzpeak",
            grid_spectra(2, 2),
            file_metadata={"instrument_configuration_list": ORBITRAP_CONFIGURATION},
        )

        comprehensive = _comprehensive(archive)

        info = comprehensive.instrument_info
        assert info["instrument_model"] == "LTQ Orbitrap Discovery"
        assert info["instrument_serial_number"] == "none"
        assert info["analyzer"] == "orbitrap"
        assert info["instrument_type"] == "Orbitrap"
        assert _detected(comprehensive).name == "Orbitrap"

    def test_a_timstof_archive_is_recognised_by_its_model(self, tmp_path):
        """The free-text model names the timsTOF, as it does for imzML."""
        archive = build_mzpeak(
            tmp_path / "timstof.mzpeak",
            grid_spectra(2, 2),
            signal="centroid",
            file_metadata={"instrument_configuration_list": TIMSTOF_CONFIGURATION},
        )

        comprehensive = _comprehensive(archive)

        assert comprehensive.instrument_info["instrument_model"] == "timsTOF Maldi 2"
        assert comprehensive.instrument_info["analyzer"] == "time-of-flight"
        assert "instrument_type" not in comprehensive.instrument_info
        assert _detected(comprehensive).name == "timsTOF"

    def test_the_raw_lists_are_kept(self, tmp_path):
        """What the archive wrote stays beside what was read from it."""
        archive = build_mzpeak(
            tmp_path / "raw.mzpeak",
            grid_spectra(2, 2),
            file_metadata={"instrument_configuration_list": ORBITRAP_CONFIGURATION},
        )

        info = _comprehensive(archive).instrument_info

        assert info["instrument_configuration_list"] == ORBITRAP_CONFIGURATION

    def test_no_instrument_block_states_no_instrument(self, tmp_path):
        """An empty list resolves to nothing, not to a guess."""
        archive = build_mzpeak(tmp_path / "none.mzpeak", grid_spectra(2, 2))

        info = _comprehensive(archive).instrument_info

        for key in ("instrument_model", "analyzer", "instrument_type"):
            assert key not in info


class TestMissingPixelSize:
    """An archive with no pixel size is refused once, saying what is missing."""

    def test_the_refusal_is_said_once(self, tmp_path, thyra_logs):
        """It was logged three times, the hint on a line of its own."""
        from thyra.convert import convert_msi

        archive = build_mzpeak(
            tmp_path / "nopixel.mzpeak", grid_spectra(2, 2), pixel_size=None
        )

        with thyra_logs("thyra.convert", logging.ERROR) as records:
            assert convert_msi(str(archive), str(tmp_path / "out.zarr")) is False

        messages = [r.getMessage() for r in records]
        assert len(messages) == 1
        assert "Pixel size not found in metadata" in messages[0]
        assert "--pixel-size" in messages[0]

    def test_the_archive_is_said_to_state_none(self, tmp_path, thyra_logs):
        """The extractor names the two terms it looked for."""
        archive = build_mzpeak(
            tmp_path / "nopixel.mzpeak", grid_spectra(2, 2), pixel_size=None
        )

        with thyra_logs("thyra.metadata", logging.INFO) as records:
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

        assert any(
            "states no pixel size" in r.getMessage()
            and "IMS:1000046" in r.getMessage()
            and "IMS:1000047" in r.getMessage()
            for r in records
        )


class TestWhereTheImageSits:
    """The base is the one D14 sets for imzML; a writer's shift is added back."""

    def test_a_declared_base_of_zero_is_used(self, tmp_path):
        """A 0-based archive that says so keeps its first row and column."""
        archive = build_mzpeak(
            tmp_path / "zero.mzpeak",
            grid_spectra(2, 2, base=2),
            file_metadata={"imaging": {"coordinate_base": 0, "is_imaging": True}},
        )

        with MzPeakReader(archive) as reader:
            essential = reader.get_essential_metadata()

        assert essential.coordinate_offsets == (0, 0, 0)
        assert essential.dimensions == (4, 4, 1)

    def test_a_position_below_the_declared_base_folds_it_down(
        self, tmp_path, thyra_logs
    ):
        """Positions from 0 in an archive that declares 1 lose nothing."""
        archive = build_mzpeak(tmp_path / "low.mzpeak", grid_spectra(2, 2, base=0))

        with thyra_logs("thyra.readers.mzpeak", logging.INFO) as records:
            with MzPeakReader(archive) as reader:
                coords = [c for c, _, _ in reader.iter_spectra()]
                essential = reader.get_essential_metadata()

        assert coords == [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)]
        assert essential.coordinate_offsets == (0, 0, 0)
        assert any("smallest is x=0, y=0" in r.getMessage() for r in records)

    def test_the_position_offset_is_added_back(self, tmp_path):
        """mzpeak-convert 0.16.0 shifts Bruker raster indices to start at 1.

        With the shift added back, the offset is the raster index of the
        first pixel, which is what a conversion of the ``.d`` records.
        """
        archive = build_mzpeak(
            tmp_path / "shifted.mzpeak",
            grid_spectra(3, 2),
            file_metadata={
                "imaging": {
                    "coordinate_base": 1,
                    "is_imaging": True,
                    "position_offset": {"x": 668, "y": 108},
                }
            },
        )

        with MzPeakReader(archive) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]
            essential = reader.get_essential_metadata()

        assert coords[0] == (0, 0, 0)
        assert essential.coordinate_offsets == (669, 109, 0)
        assert essential.dimensions == (3, 2, 1)

    def test_an_offset_that_is_not_two_whole_numbers_is_not_applied(
        self, tmp_path, thyra_logs
    ):
        """A malformed shift is named in the log and left alone."""
        archive = build_mzpeak(
            tmp_path / "badshift.mzpeak",
            grid_spectra(2, 2),
            file_metadata={"imaging": {"position_offset": {"x": "668", "y": 108}}},
        )

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                essential = reader.get_essential_metadata()

        assert essential.coordinate_offsets == (1, 1, 0)
        assert any("position_offset" in r.getMessage() for r in records)


class TestPositionZ:
    """The plane an archive states is its z offset; none stated is 0 (D30).

    mzpeak-convert 0.16.0 writes a z column when its source imzML states z,
    and none when it does not. The imzML route records the same.
    """

    def _offsets(self, archive):
        with MzPeakReader(archive) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]
            essential = reader.get_essential_metadata()
        assert {z for _, _, z in coords} == {0}
        assert essential.dimensions[2] == 1
        return essential.coordinate_offsets

    def test_a_stated_plane_is_the_z_offset(self, tmp_path):
        """As the glioma archive from 0.16.0 and its imzML both say."""
        archive = build_mzpeak(
            tmp_path / "z1.mzpeak", grid_spectra(2, 2), position_z=[1, 1, 1, 1]
        )

        assert self._offsets(archive) == (1, 1, 1)

    def test_a_plane_other_than_one_is_kept(self, tmp_path):
        """One plane of a series keeps its place in it."""
        archive = build_mzpeak(
            tmp_path / "z3.mzpeak", grid_spectra(2, 2), position_z=[3, 3, 3, 3]
        )

        assert self._offsets(archive) == (1, 1, 3)

    def test_an_archive_without_z_records_zero(self, tmp_path):
        """As DESI and the examples do: neither they nor their imzML state z."""
        archive = build_mzpeak(tmp_path / "noz.mzpeak", grid_spectra(2, 2))

        assert self._offsets(archive) == (1, 1, 0)

    def test_the_column_is_found_by_its_conventional_name(self, tmp_path):
        """An index without bindings still names it the way 0.12.0 did."""
        archive = build_mzpeak(
            tmp_path / "noname.mzpeak",
            grid_spectra(2, 2),
            column_mapping_key=None,
            position_z=[1, 1, 1, 1],
            position_z_column="opt_IMS_1000052_position_z",
        )

        assert self._offsets(archive) == (1, 1, 1)

    def test_a_column_of_nulls_states_no_plane(self, tmp_path):
        """A column a writer always adds, empty when the source had no z."""
        archive = build_mzpeak(
            tmp_path / "nullz.mzpeak",
            grid_spectra(2, 2),
            position_z=[None, None, None, None],
        )

        assert self._offsets(archive) == (1, 1, 0)

    def test_a_scan_on_no_pixel_needs_no_plane(self, tmp_path):
        """Only the spectra on a pixel are asked for one."""
        spectra = grid_spectra(2, 2) + [Spectrum(None, None, [100.0], [1.0])]
        archive = build_mzpeak(
            tmp_path / "offpixel.mzpeak", spectra, position_z=[1, 1, 1, 1, None]
        )

        assert self._offsets(archive) == (1, 1, 1)

    def test_a_plane_stated_for_some_spectra_only_is_refused(self, tmp_path):
        """A spectrum on a pixel either has a plane or none of them does."""
        archive = build_mzpeak(
            tmp_path / "somez.mzpeak", grid_spectra(2, 2), position_z=[1, 1, None, 1]
        )

        with pytest.raises(ConversionRefused, match="3 of 4 spectra"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_two_planes_are_refused(self, tmp_path):
        """Before, the z column went unread and both planes became one."""
        archive = build_mzpeak(
            tmp_path / "twoz.mzpeak", grid_spectra(2, 2), position_z=[1, 1, 2, 2]
        )

        with pytest.raises(ConversionRefused, match="2 planes"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_scans_listed_out_of_order_keep_their_plane(self, tmp_path):
        """Each spectrum's z comes from its own scan, joined on source_index.

        The scan on no pixel, which states no z, is listed first. Read in row
        order, its missing z would land on a spectrum on a pixel and refuse
        the archive.
        """
        spectra = grid_spectra(2, 2) + [Spectrum(None, None, [100.0], [1.0])]
        archive = build_mzpeak(
            tmp_path / "shuffled.mzpeak", spectra, position_z=[1, 1, 1, 1, None]
        )
        _reverse_scans(archive)

        assert self._offsets(archive) == (1, 1, 1)


def _z_parameter(z, name="position z", accession="IMS:1000052"):
    """A scan's position z parameter, as mzpeak-convert 0.12.0 and 0.14.0 write it."""
    if z is None:
        return []
    return [
        {
            "value": {"integer": z, "float": None, "string": None, "boolean": None},
            "accession": accession,
            "name": name,
            "unit": None,
        }
    ]


def _reverse_scans(archive):
    """List the rows of an archive's scans member in reverse order."""
    import io

    import pyarrow.parquet as pq

    with zipfile.ZipFile(archive) as zf:
        members = {name: zf.read(name) for name in zf.namelist()}
    scans = pq.read_table(io.BytesIO(members[SCANS_MEMBER]))
    scans = scans.take(list(range(scans.num_rows - 1, -1, -1)))
    sink = io.BytesIO()
    pq.write_table(scans, sink)
    members[SCANS_MEMBER] = sink.getvalue()
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_STORED) as zf:
        for name, data in members.items():
            zf.writestr(name, data)


class TestPositionZInScanParameters:
    """mzpeak-convert 0.12.0 and 0.14.0 state z in each scan's parameters list.

    They write no z column. Read only from a column, two planes were summed
    into one and a stated plane was recorded as 0 (D39).
    """

    def _archive(self, tmp_path, zs, spectra=None, **parameter):
        return build_mzpeak(
            tmp_path / "params.mzpeak",
            spectra or grid_spectra(2, 2),
            scan_columns={"parameters": [_z_parameter(z, **parameter) for z in zs]},
        )

    def _offsets(self, archive):
        with MzPeakReader(archive) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]
            essential = reader.get_essential_metadata()
        assert {z for _, _, z in coords} == {0}
        return essential.coordinate_offsets

    def test_a_stated_plane_is_the_z_offset(self, tmp_path):
        """As gbm_Test_P15_r2 and ltpmsi-chilli state z = 1 on every scan."""
        archive = self._archive(tmp_path, [1, 1, 1, 1])

        assert self._offsets(archive) == (1, 1, 1)

    def test_two_planes_are_refused(self, tmp_path):
        """The same refusal as for a z column."""
        archive = self._archive(tmp_path, [1, 1, 2, 2])

        with pytest.raises(ConversionRefused, match="2 planes"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_a_plane_stated_for_some_spectra_only_is_refused(self, tmp_path):
        """A scan whose list lacks z, beside scans that state it."""
        archive = self._archive(tmp_path, [1, 1, None, 1])

        with pytest.raises(ConversionRefused, match="3 of 4 spectra"):
            with MzPeakReader(archive) as reader:
                reader.get_essential_metadata()

    def test_the_parameter_is_found_by_its_accession(self, tmp_path):
        """A writer may spell the name differently; the accession is the contract."""
        archive = self._archive(tmp_path, [2, 2, 2, 2], name="z position")

        assert self._offsets(archive) == (1, 1, 2)

    def test_the_parameter_is_found_by_its_name(self, tmp_path):
        """A writer may leave the accession out, as 0.17 does for its region."""
        archive = self._archive(tmp_path, [2, 2, 2, 2], accession=None)

        assert self._offsets(archive) == (1, 1, 2)

    def test_lists_without_z_record_zero(self, tmp_path):
        """Parameters that name no plane state none."""
        archive = self._archive(tmp_path, [None, None, None, None])

        assert self._offsets(archive) == (1, 1, 0)

    def test_a_scan_on_no_pixel_needs_no_plane(self, tmp_path):
        """Only the spectra on a pixel are asked for one."""
        spectra = grid_spectra(2, 2) + [Spectrum(None, None, [100.0], [1.0])]
        archive = self._archive(tmp_path, [1, 1, 1, 1, None], spectra=spectra)

        assert self._offsets(archive) == (1, 1, 1)

    def test_scans_listed_out_of_order_keep_their_plane(self, tmp_path):
        """The parameters list follows its scan, as the z column does."""
        spectra = grid_spectra(2, 2) + [Spectrum(None, None, [100.0], [1.0])]
        archive = self._archive(tmp_path, [1, 1, 1, 1, None], spectra=spectra)
        _reverse_scans(archive)

        assert self._offsets(archive) == (1, 1, 1)

    def test_the_first_z_of_a_scan_wins(self, tmp_path):
        """A scan that states z twice is read at its first statement (D39)."""
        archive = build_mzpeak(
            tmp_path / "twice.mzpeak",
            grid_spectra(2, 2),
            scan_columns={
                "parameters": [_z_parameter(1) + _z_parameter(2) for _ in range(4)]
            },
        )

        assert self._offsets(archive) == (1, 1, 1)

    @pytest.mark.parametrize(
        "integer",
        [None, "1", 1.0],
        ids=["null-typed", "string-typed", "float-typed"],
    )
    def test_an_integer_slot_of_another_type_states_no_plane(self, tmp_path, integer):
        """Only an integer-typed slot states z (D39), so no plane is stated."""
        parameter = _z_parameter(1)
        parameter[0]["value"] = {
            "integer": integer,
            "float": 1.0,
            "string": None,
            "boolean": None,
        }
        archive = build_mzpeak(
            tmp_path / "typed.mzpeak",
            grid_spectra(2, 2),
            scan_columns={"parameters": [parameter] * 4},
        )

        assert self._offsets(archive) == (1, 1, 0)


def _two_regions(tmp_path, regions, name="regions.mzpeak"):
    """Two 2x2 blocks side by side, raster indices 101-102 and 104-105.

    Written shifted to start at 1, with the shift kept in
    ``imaging.position_offset``, as mzpeak-convert does for a Bruker run.
    """
    spectra = [
        Spectrum(x, y, [100.0, 101.0], [float(x), float(y)])
        for y in (1, 2)
        for x in (1, 2, 4, 5)
    ]
    return build_mzpeak(
        tmp_path / name,
        spectra,
        file_metadata={
            "imaging": {"coordinate_base": 1, "position_offset": {"x": 100, "y": 10}},
            "bruker_maldi": {"regions": regions},
        },
    )


class TestRegions:
    """Regions are read when their boxes place every pixel once."""

    def test_two_regions_are_read(self, tmp_path):
        """Each pixel gets the region whose box holds it."""
        archive = _two_regions(
            tmp_path,
            [
                _region(0, 4, (101, 102), (11, 12), name="left"),
                _region(1, 4, (104, 105), (11, 12)),
            ],
        )

        with MzPeakReader(archive) as reader:
            region_map = reader.get_region_map()
            info = reader.get_region_info()

        assert region_map == {
            (0, 0): 0,
            (1, 0): 0,
            (0, 1): 0,
            (1, 1): 0,
            (3, 0): 1,
            (4, 0): 1,
            (3, 1): 1,
            (4, 1): 1,
        }
        assert info == [
            {
                "region_number": 0,
                "n_spectra": 4,
                "bounds": (0, 0, 1, 1),
                "name": "left",
            },
            {"region_number": 1, "n_spectra": 4, "bounds": (3, 0, 4, 1)},
        ]

    def test_one_region_is_a_map_as_from_the_d(self, tmp_path):
        """One listed region keeps its number and box, as a ``.d`` does (D57)."""
        archive = _two_regions(tmp_path, [_region(3, 8, (101, 105), (11, 12))])

        with MzPeakReader(archive) as reader:
            assert set(reader.get_region_map().values()) == {3}
            (entry,) = reader.get_region_info()

        assert entry["region_number"] == 3
        assert entry["n_spectra"] == 8

    @pytest.mark.parametrize(
        "regions",
        [
            # The boxes overlap on the middle column.
            [_region(0, 4, (101, 104), (11, 12)), _region(1, 4, (104, 105), (11, 12))],
            # A frame count that is not the number of pixels in the box.
            [_region(0, 5, (101, 102), (11, 12)), _region(1, 4, (104, 105), (11, 12))],
            # A pixel in no box.
            [_region(0, 2, (101, 102), (11, 11)), _region(1, 4, (104, 105), (11, 12))],
        ],
        ids=["overlap", "count", "outside"],
    )
    def test_boxes_that_do_not_place_every_pixel_once_are_not_read(
        self, tmp_path, thyra_logs, regions
    ):
        """No scan names its region, so a box that is not exact is no answer."""
        archive = _two_regions(tmp_path, regions)

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                assert reader.get_region_map() is None
                assert reader.get_region_info() is None

        assert any("acquisition regions" in r.getMessage() for r in records)

    def test_selecting_a_region_is_refused(self, tmp_path):
        """It was accepted and ignored, which converted every region."""
        archive = _two_regions(
            tmp_path,
            [
                _region(0, 4, (101, 102), (11, 12)),
                _region(1, 4, (104, 105), (11, 12)),
            ],
        )

        with pytest.raises(ConversionRefused, match="region=1"):
            MzPeakReader(archive, region=1)

    def test_an_incomplete_region_stops_the_reading(self, tmp_path, thyra_logs):
        """A region without its index ranges cannot be checked."""
        archive = _two_regions(
            tmp_path,
            [_region(0, 4, (101, 102), (11, 12)), {"region_number": 1, "frames": 4}],
        )

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                assert reader.get_region_map() is None

        assert any("regions are not read" in r.getMessage() for r in records)


class TestMsLevels:
    """Spectra of other MS levels are not summed into the MS1 pixels."""

    def test_ms2_spectra_beside_ms1_are_left_out(self, tmp_path, thyra_logs):
        """Each pixel keeps its MS1 spectrum, and the log counts the rest."""
        ms1 = grid_spectra(2, 1)
        ms2 = [Spectrum(s.x, s.y, [50.0, 60.0], [1e6, 1e6]) for s in ms1]
        archive = build_mzpeak(
            tmp_path / "levels.mzpeak",
            [ms1[0], ms2[0], ms1[1], ms2[1]],
            ms_levels=[1, 2, 1, 2],
        )

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                emitted = list(reader.iter_spectra())
                n_spectra = reader.n_spectra
                axis = reader.get_common_mass_axis()
                essential = reader.get_essential_metadata()

        assert n_spectra == 2
        assert [coords for coords, _, _ in emitted] == [(0, 0, 0), (1, 0, 0)]
        for (_, mzs, _), expected in zip(emitted, ms1):
            np.testing.assert_array_equal(mzs, expected.mzs)
        assert any("level 2: 2" in r.getMessage() for r in records)
        # Nor do they shape the axis, the range or the count.
        np.testing.assert_array_equal(axis, ms1[0].mzs)
        assert essential.mass_range == (100.0, float(ms1[0].mzs.max()))
        assert essential.total_peaks == sum(s.mzs.size for s in ms1)

    def test_an_ms1_scan_on_no_pixel_does_not_cost_the_placed_ones(self, tmp_path):
        """A survey scan with no position, then an MS2 image."""
        survey = Spectrum(None, None, [100.0, 101.0], [1.0, 1.0])
        archive = build_mzpeak(
            tmp_path / "survey.mzpeak",
            [survey] + grid_spectra(2, 1),
            ms_levels=[1, 2, 2],
        )

        with MzPeakReader(archive) as reader:
            coords = [c for c, _, _ in reader.iter_spectra()]

        assert coords == [(0, 0, 0), (1, 0, 0)]

    def test_a_level_that_is_not_stated_is_kept(self, tmp_path):
        """mzpeak-convert writes level 0 when the source states none."""
        archive = build_mzpeak(
            tmp_path / "unstated.mzpeak", grid_spectra(2, 1), ms_levels=[1, 0]
        )

        with MzPeakReader(archive) as reader:
            assert reader.n_spectra == 2

    def test_an_archive_without_ms1_is_read_whole(self, tmp_path):
        """Leaving out every spectrum would convert nothing."""
        archive = build_mzpeak(
            tmp_path / "ms2only.mzpeak", grid_spectra(2, 1), ms_levels=[2, 2]
        )

        with MzPeakReader(archive) as reader:
            assert reader.n_spectra == 2


class TestSharedPixels:
    """Spectra that share a pixel are summed there, and the log says so."""

    def test_two_spectra_on_one_pixel_are_reported(self, tmp_path, thyra_logs):
        """Three spectra, two pixels."""
        spectra = grid_spectra(2, 1) + [Spectrum(1, 1, [100.0], [5.0])]
        archive = build_mzpeak(tmp_path / "shared.mzpeak", spectra)

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                _ = reader.n_spectra

        assert any("places 3 spectra on 2 pixels" in r.getMessage() for r in records)


class TestIonMobility:
    """Ion mobility is not read, and the log says so."""

    def test_scans_with_a_mobility_value_are_reported(self, tmp_path, thyra_logs):
        """One message, with the count."""
        archive = build_mzpeak(
            tmp_path / "with_values.mzpeak",
            grid_spectra(2, 1),
            ion_mobility=[0.8, None],
        )

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                _ = reader.n_spectra

        assert any(
            "1 of 2 scans" in r.getMessage() and "ion mobility" in r.getMessage()
            for r in records
        )

    def test_a_column_of_nulls_is_not_reported(self, tmp_path, thyra_logs):
        """mzpeak-convert writes the column for every run, empty when unused."""
        archive = build_mzpeak(
            tmp_path / "empty_column.mzpeak",
            grid_spectra(2, 1),
            ion_mobility=[None, None],
        )

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                _ = reader.n_spectra

        assert not any("mobility" in r.getMessage() for r in records)


def _with_images(archive, members, entity_type="image"):
    """Add members to an archive, listed in its index as ``entity_type``."""
    with zipfile.ZipFile(archive) as source:
        contents = {name: source.read(name) for name in source.namelist()}
    index = json.loads(contents[INDEX_MEMBER])
    for name, payload in members.items():
        contents[name] = payload
        index["files"].append(
            {
                "name": name,
                "entity_type": entity_type,
                "data_kind": "proprietary",
                "column_mapping": [],
                "parameters": [],
            }
        )
    contents[INDEX_MEMBER] = json.dumps(index).encode()
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_STORED) as target:
        for name, payload in contents.items():
            target.writestr(name, payload)
    return archive


def _corrupted(archive, payload):
    """Flip one byte of a stored member's data, so its CRC no longer matches."""
    data = bytearray(archive.read_bytes())
    start = bytes(data).index(payload)
    data[start + len(payload) // 2] ^= 0xFF
    archive.write_bytes(bytes(data))
    return archive


def _png_bytes(height=4, width=5):
    """A small, valid PNG."""
    import io

    from PIL import Image

    buffer = io.BytesIO()
    Image.fromarray(np.zeros((height, width, 3), np.uint8)).save(buffer, "PNG")
    return buffer.getvalue()


def _described(members):
    """``imaging.images`` naming each member's media type, as 0.12.0 writes it."""
    return {
        "imaging": {
            "images": [
                {"archive_path": name, "media_type": media_type, "role": "optical"}
                for name, media_type in members.items()
            ]
        }
    }


class TestEmbeddedImages:
    """Images the archive embeds reach the optical image loader."""

    def test_an_image_is_extracted_and_removed_on_close(self, tmp_path):
        """Byte for byte, under its own name, gone once the reader closes."""
        payload = b"\x89PNG\r\n\x1a\n not decoded here"
        archive = _with_images(
            build_mzpeak(tmp_path / "image.mzpeak", grid_spectra(2, 2)),
            {"images/image_0000.png": payload},
        )

        reader = MzPeakReader(archive)
        (path,) = reader.get_optical_image_paths()
        assert path.name == "image_0000.png"
        assert path.read_bytes() == payload
        assert reader.get_optical_image_paths() == [path]
        reader.close()

        assert not path.exists()
        assert not path.parent.exists()
        assert reader.get_optical_image_label(path) is None

    def test_a_tiff_under_another_name_is_read_as_a_tiff(self, tmp_path):
        """The glioma example embeds an Aperio slide scan as ``.svs``.

        The archive declares it ``image/tiff``; the loader picks its decoder
        by suffix, so the copy is a ``.tif``. The store still records the
        member it came from.
        """
        payload = b"II*\x00 stands in for a slide scan"
        member = "images/image_0000.svs"
        archive = _with_images(
            build_mzpeak(
                tmp_path / "svs.mzpeak",
                grid_spectra(2, 2),
                file_metadata=_described({member: "image/tiff"}),
            ),
            {member: payload},
        )

        with MzPeakReader(archive) as reader:
            (path,) = reader.get_optical_image_paths()
            label = reader.get_optical_image_label(path)
            content = path.read_bytes()

        assert path.name == "image_0000.tif"
        assert content == payload
        assert label == OpticalImageLabel("image_0000", member)

    def test_the_declared_media_type_decides(self, tmp_path):
        """A member without an image suffix is read as what it is declared."""
        archive = _with_images(
            build_mzpeak(
                tmp_path / "declared.mzpeak",
                grid_spectra(2, 2),
                file_metadata=_described({"images/scan.bin": "image/png"}),
            ),
            {"images/scan.bin": b"\x89PNG\r\n\x1a\n"},
        )

        with MzPeakReader(archive) as reader:
            names = [path.name for path in reader.get_optical_image_paths()]

        assert names == ["scan.png"]

    def test_a_declared_format_the_loader_does_not_read_is_skipped(
        self, tmp_path, thyra_logs
    ):
        """The log names the member and the type it is declared as."""
        member = "images/overlay.svg"
        archive = _with_images(
            build_mzpeak(
                tmp_path / "svg_declared.mzpeak",
                grid_spectra(2, 2),
                file_metadata=_described({member: "image/svg+xml"}),
            ),
            {member: b"<svg/>"},
        )

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                assert reader.get_optical_image_paths() == []

        assert any(
            member in r.getMessage() and "image/svg+xml" in r.getMessage()
            for r in records
        )

    def test_each_image_is_named_after_its_member(self, tmp_path):
        """Not the high resolution scan and the derived image of a vendor folder."""
        members = {
            "images/image_0000.png": b"first",
            "images/image_0001.png": b"second",
        }
        archive = _with_images(
            build_mzpeak(tmp_path / "two.mzpeak", grid_spectra(2, 2)), members
        )

        with MzPeakReader(archive) as reader:
            labels = [
                reader.get_optical_image_label(path)
                for path in reader.get_optical_image_paths()
            ]

        assert labels == [
            OpticalImageLabel("image_0000", "images/image_0000.png"),
            OpticalImageLabel("image_0001", "images/image_0001.png"),
        ]

    @pytest.mark.parametrize("upper_first", [True, False])
    def test_names_that_differ_only_in_case_stay_apart(self, tmp_path, upper_first):
        """The store lowercases an element name, so two must not meet there.

        Both orders: the first name taken is looked up by the second, and
        either side being case-sensitive lets one order through.
        """
        members = {"images/Scan.png": b"upper", "images/scan.png": b"lower"}
        if not upper_first:
            members = dict(reversed(list(members.items())))
        archive = _with_images(
            build_mzpeak(tmp_path / "case.mzpeak", grid_spectra(2, 2)), members
        )

        with MzPeakReader(archive) as reader:
            paths = reader.get_optical_image_paths()
            names = [reader.get_optical_image_label(path).name for path in paths]
            contents = [path.read_bytes() for path in paths]

        assert len({name.lower() for name in names}) == 2
        assert len({path.name.lower() for path in paths}) == 2
        assert contents == list(members.values())

    def test_images_of_one_name_in_two_folders_are_both_kept(self, tmp_path):
        """A folder inside the archive is dropped, so names can meet."""
        members = {
            "images/a_1.png": b"first",
            "images/a.png": b"second",
            "other/a.png": b"third",
        }
        archive = _with_images(
            build_mzpeak(tmp_path / "names.mzpeak", grid_spectra(2, 2)), members
        )

        with MzPeakReader(archive) as reader:
            paths = reader.get_optical_image_paths()
            contents = [path.read_bytes() for path in paths]

        assert len(set(paths)) == 3
        assert contents == [b"first", b"second", b"third"]

    def test_a_format_the_loader_does_not_read_is_skipped(self, tmp_path, thyra_logs):
        """Named in the log, not handed to the loader."""
        archive = _with_images(
            build_mzpeak(tmp_path / "svg.mzpeak", grid_spectra(2, 2)),
            {"images/overlay.svg": b"<svg/>"},
        )

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                assert reader.get_optical_image_paths() == []

        assert any("images/overlay.svg" in r.getMessage() for r in records)

    def test_an_archive_without_images_gives_none(self, tmp_path):
        """No folder is made for nothing."""
        archive = build_mzpeak(tmp_path / "plain.mzpeak", grid_spectra(2, 2))

        with MzPeakReader(archive) as reader:
            assert reader.get_optical_image_paths() == []
            assert reader._image_dir is None

    def test_a_member_cannot_choose_where_its_copy_goes(self, tmp_path):
        """F80: a name with ``..`` or folders is a name inside the archive."""
        # A name of its own, so a copy left in the shared temporary folder by
        # another run cannot pass or fail this one.
        stray = f"escape_probe_{uuid.uuid4().hex}"
        members = {
            f"../{stray}.png": b"up",
            "deep/sub/x.png": b"down",
            "/abs/rooted.png": b"root",
        }
        archive = _with_images(
            build_mzpeak(tmp_path / "escape.mzpeak", grid_spectra(2, 2)), members
        )

        with MzPeakReader(archive) as reader:
            paths = reader.get_optical_image_paths()
            image_dir = reader._image_dir
            assert len(paths) == 3
            assert all(path.parent == image_dir for path in paths)
            assert sorted(path.name for path in image_dir.iterdir()) == sorted(
                [f"{stray}.png", "rooted.png", "x.png"]
            )
            assert not (image_dir.parent / f"{stray}.png").exists()
        assert not (tmp_path / f"{stray}.png").exists()

    def test_the_declared_media_type_wins_over_a_wrong_suffix(self, tmp_path):
        """F94: a PNG named ``.jpg``, a TIFF named ``.dat`` declared in mixed case."""
        members = {
            "images/photo.jpg": b"\x89PNG\r\n\x1a\n",
            "images/scan.dat": b"II*\x00",
        }
        archive = _with_images(
            build_mzpeak(
                tmp_path / "wrong_suffix.mzpeak",
                grid_spectra(2, 2),
                file_metadata=_described(
                    {"images/photo.jpg": "image/png", "images/scan.dat": " Image/TIFF "}
                ),
            ),
            members,
        )

        with MzPeakReader(archive) as reader:
            names = [path.name for path in reader.get_optical_image_paths()]

        assert names == ["photo.png", "scan.tif"]

    def test_a_member_that_is_not_an_image_entity_is_not_copied(self, tmp_path):
        """F204: an ``annotation`` member is not an image, whatever its suffix."""
        archive = _with_images(
            build_mzpeak(tmp_path / "annotation.mzpeak", grid_spectra(2, 2)),
            {"images/overlay.png": _png_bytes()},
            entity_type="annotation",
        )

        with MzPeakReader(archive) as reader:
            assert reader.get_optical_image_paths() == []
            assert reader._image_dir is None

    def test_a_long_mixed_case_member_makes_a_short_lowercase_element(self, tmp_path):
        """F204: 41 characters, upper and lower case, a 30-character element name."""
        from thyra.converters.spatialdata.optical_image import OpticalImages

        stem = "Scan" + "A" * 37
        assert len(stem) == 41
        archive = _with_images(
            build_mzpeak(tmp_path / "long.mzpeak", grid_spectra(2, 2)),
            {f"images/{stem}.png": _png_bytes()},
        )

        with MzPeakReader(archive) as reader:
            optical = OpticalImages(
                reader,
                lambda: tmp_path / "long.zarr",
                "ds",
                pixel_size_xy=lambda: (10.0, 10.0),
                include=True,
                apply_alignment=True,
            )
            optical.compute_alignment()
            optical.add_images({})

        (name,) = optical.pending
        suffix = name.removeprefix("ds_optical_")
        assert suffix == stem.lower()[:30]
        assert len(suffix) <= 30 and suffix == suffix.lower()
        assert optical.sources == {name: f"images/{stem}.png"}

    def test_a_member_with_a_bad_crc_loses_only_its_image(self, tmp_path, thyra_logs):
        """F41: the damaged member is named in a warning; the rest are copied."""
        bad = b"CORRUPT-PAYLOAD-0123456789"
        archive = _corrupted(
            _with_images(
                build_mzpeak(tmp_path / "crc.mzpeak", grid_spectra(2, 2)),
                {"images/bad.png": bad, "images/good.png": b"good bytes"},
            ),
            bad,
        )

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                paths = reader.get_optical_image_paths()
                names = sorted(path.name for path in reader._image_dir.iterdir())
                content = [path.read_bytes() for path in paths]

        assert [path.name for path in paths] == ["good.png"]
        assert content == [b"good bytes"]
        assert names == ["good.png"]  # no half-written copy is left behind
        messages = [r.getMessage() for r in records if r.levelno == logging.WARNING]
        assert any(
            "images/bad.png" in m and "BadZipFile" in m for m in messages
        ), messages

    def test_a_broken_deflated_member_loses_only_its_image(self, tmp_path, thyra_logs):
        """F41: a deflate stream that breaks off raises zlib.error, not BadZipFile."""
        archive = _with_images(
            build_mzpeak(tmp_path / "deflate.mzpeak", grid_spectra(2, 2)),
            {"images/good.png": b"good bytes"},
        )
        member = "images/bad.png"
        with zipfile.ZipFile(archive) as source:
            contents = {name: source.read(name) for name in source.namelist()}
        index = json.loads(contents[INDEX_MEMBER])
        index["files"].append(
            {"name": member, "entity_type": "image", "data_kind": "proprietary"}
        )
        contents[INDEX_MEMBER] = json.dumps(index).encode()
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_STORED) as target:
            for name, payload in contents.items():
                target.writestr(name, payload)
            target.writestr(member, bytes(range(256)) * 64, zipfile.ZIP_DEFLATED)
        with zipfile.ZipFile(archive) as source:
            info = source.getinfo(member)
        data = bytearray(archive.read_bytes())
        header = info.header_offset
        name_length = int.from_bytes(data[header + 26 : header + 28], "little")
        extra_length = int.from_bytes(data[header + 28 : header + 30], "little")
        start = header + 30 + name_length + extra_length
        for offset in range(2, info.compress_size - 2):
            data[start + offset] ^= 0xA5
        archive.write_bytes(bytes(data))

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                paths = reader.get_optical_image_paths()

        assert [path.name for path in paths] == ["good.png"]
        messages = [r.getMessage() for r in records if r.levelno == logging.WARNING]
        assert any(member in m and "cannot be copied" in m for m in messages), messages

    def test_a_member_that_cannot_be_opened_loses_only_its_image(
        self, tmp_path, thyra_logs, monkeypatch
    ):
        """F41: an OSError from the archive is a skipped image, not a failure."""
        from thyra.readers.mzpeak.mzpeak_reader import MzPeakArchive

        archive = _with_images(
            build_mzpeak(tmp_path / "oserror.mzpeak", grid_spectra(2, 2)),
            {"images/locked.png": b"locked", "images/open.png": b"open"},
        )
        real = MzPeakArchive.open_member

        def refuse(self, name):
            if name == "images/locked.png":
                raise OSError("device not ready")
            return real(self, name)

        monkeypatch.setattr(MzPeakArchive, "open_member", refuse)

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                names = [path.name for path in reader.get_optical_image_paths()]

        assert names == ["open.png"]
        assert any(
            "images/locked.png" in r.getMessage() and "OSError" in r.getMessage()
            for r in records
        )

    def test_a_member_name_a_disk_refuses_is_copied_under_a_safe_name(self, tmp_path):
        """F41: ``?`` and ``:`` are not file-name characters on every platform."""
        members = {"images/what?.png": b"question", "images/a:b.png": b"colon"}
        archive = _with_images(
            build_mzpeak(tmp_path / "unsafe.mzpeak", grid_spectra(2, 2)), members
        )

        with MzPeakReader(archive) as reader:
            paths = reader.get_optical_image_paths()
            contents = [path.read_bytes() for path in paths]
            labels = [reader.get_optical_image_label(path) for path in paths]

        assert [path.name for path in paths] == ["what_.png", "a_b.png"]
        assert contents == [b"question", b"colon"]
        assert [label.source_file for label in labels] == list(members)

    @pytest.mark.parametrize(
        "member, stem",
        [
            ("images/CON.png", "CON_"),
            ("images/nul.png", "nul_"),
            ("images/Com1.png", "Com1_"),
            ("images/con.tar.png", "con.tar_"),
            ("images/.hidden .png", "hidden"),
            ("images/" + "d" * 99 + ".x.png", "d" * 99),
            ("images/....png", "image"),
            ("images/caf\u00e9.png", "caf_"),
            ("images/" + "n" * 150 + ".png", "n" * 100),
            ("images/COM10.png", "COM10"),
        ],
    )
    def test_the_copy_stem_is_safe_on_any_disk(self, member, stem):
        assert MzPeakReader._copy_stem(member) == stem

    def test_a_warning_for_an_undecodable_image_names_the_member(
        self, tmp_path, thyra_logs
    ):
        """F165: the member, never the temporary folder the copy sits in."""
        from thyra.converters.spatialdata.optical_image import OpticalImages

        member = "images/overview.jpeg"
        archive = _with_images(
            build_mzpeak(tmp_path / "broken.mzpeak", grid_spectra(2, 2)),
            {member: b"this is not a jpeg"},
        )

        with MzPeakReader(archive) as reader:
            optical = OpticalImages(
                reader,
                lambda: tmp_path / "broken.zarr",
                "ds",
                pixel_size_xy=lambda: (10.0, 10.0),
                include=True,
                apply_alignment=True,
            )
            optical.compute_alignment()
            with thyra_logs("thyra.converters", logging.WARNING) as records:
                optical.add_images({})
            image_dir = str(reader._image_dir)

        messages = [r.getMessage() for r in records]
        assert optical.pending == {}
        assert any(member in m for m in messages), messages
        assert not any(image_dir in m or "thyra-mzpeak-images-" in m for m in messages)


# ----------------------------------------------------------------------
# Issue #429: what mzpeak-convert 0.17 states for a Bruker run
# ----------------------------------------------------------------------


def _named_regions(tmp_path, numbers, regions, name="named.mzpeak"):
    """The layout of :func:`_two_regions`, each scan naming its region."""
    spectra = [
        Spectrum(x, y, [100.0, 101.0], [float(x), float(y)])
        for y in (1, 2)
        for x in (1, 2, 4, 5)
    ]
    return build_mzpeak(
        tmp_path / name,
        spectra,
        file_metadata={
            "imaging": {"coordinate_base": 1, "position_offset": {"x": 100, "y": 10}},
            "bruker_maldi": {"regions": regions},
        },
        scan_columns={
            "parameters": [
                [] if number is None else _region_parameter(number)
                for number in numbers
            ]
        },
    )


class TestRegionNamedByEachScan:
    """From 0.17 each scan names its region, and that is what is read."""

    #: Left block region 0, right block region 1, both rows.
    NUMBERS = [0, 0, 1, 1, 0, 0, 1, 1]

    def test_the_named_region_wins_where_the_boxes_cannot_tell(self, tmp_path):
        """Overlapping boxes, which alone leave the archive one region."""
        archive = _named_regions(
            tmp_path,
            self.NUMBERS,
            [
                _region(0, 4, (101, 104), (11, 12), name="left"),
                _region(1, 4, (104, 105), (11, 12)),
            ],
        )

        with MzPeakReader(archive) as reader:
            region_map = reader.get_region_map()
            info = reader.get_region_info()

        assert region_map == {
            (0, 0): 0,
            (1, 0): 0,
            (0, 1): 0,
            (1, 1): 0,
            (3, 0): 1,
            (4, 0): 1,
            (3, 1): 1,
            (4, 1): 1,
        }
        assert info == [
            {
                "region_number": 0,
                "n_spectra": 4,
                "bounds": (0, 0, 3, 1),
                "name": "left",
            },
            {"region_number": 1, "n_spectra": 4, "bounds": (3, 0, 4, 1)},
        ]

    def test_a_region_no_listing_describes_is_boxed_by_its_pixels(self, tmp_path):
        """The scans are read even when ``bruker_maldi`` lists nothing."""
        archive = _named_regions(tmp_path, self.NUMBERS, [])

        with MzPeakReader(archive) as reader:
            info = reader.get_region_info()

        assert info == [
            {"region_number": 0, "n_spectra": 4, "bounds": (0, 0, 1, 1)},
            {"region_number": 1, "n_spectra": 4, "bounds": (3, 0, 4, 1)},
        ]

    def test_one_scan_without_a_region_falls_back_to_the_boxes(self, tmp_path):
        """A region read for some pixels and guessed for others is neither."""
        numbers = list(self.NUMBERS)
        numbers[0] = None
        archive = _named_regions(
            tmp_path,
            numbers,
            [_region(0, 4, (101, 102), (11, 12)), _region(1, 4, (104, 105), (11, 12))],
        )

        with MzPeakReader(archive) as reader:
            region_map = reader.get_region_map()

        assert region_map is not None
        assert region_map[(0, 0)] == 0
        assert region_map[(3, 0)] == 1

    def test_one_named_region_is_a_map_as_from_the_d(self, tmp_path):
        """As for a ``.d`` of one region (D57)."""
        archive = _named_regions(tmp_path, [3] * 8, [])

        with MzPeakReader(archive) as reader:
            assert set(reader.get_region_map().values()) == {3}
            assert [e["region_number"] for e in reader.get_region_info()] == [3]

    def test_the_first_region_of_a_scan_wins(self, tmp_path):
        """A scan that names two regions is read at the first one."""
        spectra = [
            Spectrum(x, y, [100.0, 101.0], [float(x), float(y)])
            for y in (1, 2)
            for x in (1, 2, 4, 5)
        ]
        archive = build_mzpeak(
            tmp_path / "twice.mzpeak",
            spectra,
            file_metadata={
                "imaging": {
                    "coordinate_base": 1,
                    "position_offset": {"x": 100, "y": 10},
                }
            },
            scan_columns={
                "parameters": [
                    _region_parameter(number) + _region_parameter(1 - number)
                    for number in self.NUMBERS
                ]
            },
        )

        with MzPeakReader(archive) as reader:
            region_map = reader.get_region_map()

        assert region_map is not None
        assert region_map[(0, 0)] == 0
        assert region_map[(3, 0)] == 1

    def test_a_null_typed_region_falls_back_to_the_boxes(self, tmp_path):
        """No scan names a region as a whole number, so the boxes decide."""
        spectra = [
            Spectrum(x, y, [100.0, 101.0], [float(x), float(y)])
            for y in (1, 2)
            for x in (1, 2, 4, 5)
        ]
        parameter = _region_parameter(0)
        parameter[0]["value"] = {
            "integer": None,
            "float": 0.0,
            "string": None,
            "boolean": None,
        }
        archive = build_mzpeak(
            tmp_path / "typed.mzpeak",
            spectra,
            file_metadata={
                "imaging": {
                    "coordinate_base": 1,
                    "position_offset": {"x": 100, "y": 10},
                },
                "bruker_maldi": {
                    "regions": [
                        _region(0, 4, (101, 102), (11, 12)),
                        _region(1, 4, (104, 105), (11, 12)),
                    ]
                },
            },
            scan_columns={"parameters": [parameter] * 8},
        )

        with MzPeakReader(archive) as reader:
            region_map = reader.get_region_map()

        assert region_map is not None
        assert region_map[(0, 0)] == 0
        assert region_map[(3, 0)] == 1


def _source(accession, name):
    """``file_description`` naming one source file of a format."""
    return {
        "file_description": {
            "contents": [
                {
                    "name": "profile spectrum",
                    "accession": "MS:1000128",
                    "value": None,
                    "unit": None,
                }
            ],
            "source_files": [
                {
                    "id": name,
                    "name": name,
                    "location": "file://",
                    "parameters": [_parameter("Bruker format", accession)],
                }
            ],
        }
    }


TDF_SOURCE = _source("MS:1002817", "analysis.tdf")
TSF_SOURCE = _source("MS:1003282", "analysis.tsf")


def _tdf_archive(tmp_path, times, file_metadata=TDF_SOURCE, name="tdf.mzpeak"):
    """Two pixels of raw counts; the first holds two points on one m/z."""
    spectra = [
        Spectrum(1, 1, [100.0, 100.0, 101.0, 102.0], [1.0, 1.0, 3.0, 5.0]),
        Spectrum(2, 1, [100.0, 103.0], [7.0, 2.0]),
    ]
    columns = {} if times is None else {"ion_injection_time": times}
    return build_mzpeak(
        tmp_path / name,
        spectra,
        file_metadata=file_metadata,
        scan_columns=columns,
    )


class TestTdfIntensityScale:
    """A TDF archive's raw counts are put on the scale the ``.d`` gives."""

    def test_counts_are_scaled_half_up_and_summed_per_m_z(self, tmp_path):
        """At 200 ms every odd count is a tie; Bruker's library rounds it up.

        1, 1, 3 and 5 become 0.5, 0.5, 1.5 and 2.5, which half up gives
        1, 1, 2 and 3, and the two points on m/z 100 sum to 2. Rounding to
        even would give 0, 0, 2 and 2.
        """
        archive = _tdf_archive(tmp_path, [200.0, 50.0])

        with MzPeakReader(archive) as reader:
            spectra = list(reader.iter_spectra())

        (_, mzs0, ints0), (_, mzs1, ints1) = spectra
        np.testing.assert_array_equal(mzs0, [100.0, 101.0, 102.0])
        np.testing.assert_array_equal(ints0, [2.0, 2.0, 3.0])
        np.testing.assert_array_equal(mzs1, [100.0, 103.0])
        np.testing.assert_array_equal(ints1, [14.0, 4.0])

    def test_the_store_records_the_scale(self, tmp_path):
        """Beside the data, so a reader of the store knows what it holds."""
        archive = _tdf_archive(tmp_path, [200.0, 50.0])

        with MzPeakReader(archive) as reader:
            facts = reader.get_comprehensive_metadata().format_specific

        assert "accumulation time" in facts["intensity_scale"]

    def test_another_source_keeps_its_values(self, tmp_path):
        """A TSF run's values are what the ``.d`` gives already."""
        archive = _tdf_archive(tmp_path, [200.0, 50.0], file_metadata=TSF_SOURCE)

        with MzPeakReader(archive) as reader:
            (_, mzs, ints), _ = list(reader.iter_spectra())
            facts = reader.get_comprehensive_metadata().format_specific

        np.testing.assert_array_equal(mzs, [100.0, 100.0, 101.0, 102.0])
        np.testing.assert_array_equal(ints, [1.0, 1.0, 3.0, 5.0])
        assert "intensity_scale" not in facts

    @pytest.mark.parametrize("times", [None, [200.0, None], [200.0, 0.0]])
    def test_without_every_time_the_counts_stay_raw_and_it_says_so(
        self, tmp_path, thyra_logs, times
    ):
        """A scale for some pixels and not others would be no scale."""
        archive = _tdf_archive(tmp_path, times)

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                (_, _, ints), _ = list(reader.iter_spectra())

        np.testing.assert_array_equal(ints, [1.0, 1.0, 3.0, 5.0])
        assert any("accumulation time" in r.getMessage() for r in records)


class TestVendorMassRange:
    """A Bruker archive takes the run's acquisition range, as the ``.d`` does."""

    def test_the_acquisition_range_is_the_mass_range(self, tmp_path):
        """Not the observed one, which is narrower and shifts every bin."""
        archive = build_mzpeak(
            tmp_path / "range.mzpeak",
            grid_spectra(2, 2),
            file_metadata={
                "vendor_metadata": {
                    "MzAcqRangeLower": "50.000000",
                    "MzAcqRangeUpper": "2000.000000",
                }
            },
        )

        with MzPeakReader(archive) as reader:
            assert reader.get_essential_metadata().mass_range == (50.0, 2000.0)

    @pytest.mark.parametrize(
        "vendor",
        [
            None,
            {"MzAcqRangeLower": "50"},
            {"MzAcqRangeLower": "x", "MzAcqRangeUpper": "9"},
            {"MzAcqRangeLower": "900", "MzAcqRangeUpper": "100"},
        ],
        ids=["absent", "half", "not-a-number", "reversed"],
    )
    def test_otherwise_the_observed_range(self, tmp_path, vendor):
        """An archive of any other source, or one without a usable range."""
        spectra = grid_spectra(2, 2)
        archive = build_mzpeak(
            tmp_path / "observed.mzpeak",
            spectra,
            file_metadata={} if vendor is None else {"vendor_metadata": vendor},
        )

        with MzPeakReader(archive) as reader:
            low, high = reader.get_essential_metadata().mass_range

        assert low == min(float(s.mzs.min()) for s in spectra)
        assert high == max(float(s.mzs.max()) for s in spectra)


#: image px -> pixel position, slightly rotated: 100 image pixels a pixel.
TEACH_POINTS_MATRIX = [0.01, 0.0004, -1.5, -0.0003, 0.0098, 2.25]


#: The ``registration`` block mzpeak-convert 0.17 writes beside the affine.
REGISTRATION = {
    "sequence": "run.mis",
    "teach_points": [
        {"image_px": [1204.0, 778.0], "stage_um": [-22965.0, 15855.0]},
        {"image_px": [6706.0, 648.0], "stage_um": [21901.0, 17125.0]},
        {"image_px": [2550.0, 4990.0], "stage_um": [-11848.0, -18292.0]},
    ],
}


def _registered(
    tmp_path,
    quality="teach_points",
    matrix=TEACH_POINTS_MATRIX,
    regions=(0, 0, 1) * 2,
    listed=(),
    media_type="image/png",
    member="images/image_0000.png",
):
    """A 3 x 2 run at raster 101-103, 11-12, its image registered to it."""
    spectra = [
        Spectrum(x, y, [100.0, 101.0], [1.0, 2.0]) for y in (1, 2) for x in (1, 2, 3)
    ]
    archive = build_mzpeak(
        tmp_path / f"{quality}.mzpeak",
        spectra,
        file_metadata={
            "imaging": {
                "coordinate_base": 1,
                "position_offset": {"x": 100, "y": 10},
                "images": [
                    {
                        "archive_path": member,
                        "media_type": media_type,
                        "role": "optical",
                        "source_name": "IMG_0000.jpg",
                        "registration": REGISTRATION,
                        "affine": {
                            "maps": "image_px -> ms_px",
                            "matrix": matrix,
                            "registration_quality": quality,
                            "type": "affine",
                        },
                    }
                ],
            },
            "bruker_maldi": {"regions": list(listed)},
        },
        scan_columns={"parameters": [_region_parameter(n) for n in regions]},
    )
    return _with_images(archive, {member: b"\x89PNG not decoded here"})


class TestTeachPointsImage:
    """An image fitted to the teaching points places the pixels on it."""

    def test_each_pixel_centre_lands_where_the_archive_maps_it(self, tmp_path):
        """The archive's affine takes that image point to the pixel's position."""
        archive = _registered(tmp_path)
        matrix = np.array(TEACH_POINTS_MATRIX).reshape(2, 3)

        with MzPeakReader(archive) as reader:
            alignment = reader.get_image_alignment()
            primary = reader.get_primary_optical_image_path()
            label = reader.get_optical_image_label(primary)

        assert alignment is not None and alignment.lattice is not None
        assert label.source_file == "images/image_0000.png"
        assert (alignment.first_raster_x, alignment.first_raster_y) == (101, 11)
        for x in range(3):
            for y in range(2):
                image = alignment.transform_point(x, y)
                position = matrix @ np.array([image[0], image[1], 1.0])
                np.testing.assert_allclose(position, [x + 1, y + 1], atol=1e-9)

    def test_a_cell_is_one_pixel_wide_on_the_image(self, tmp_path):
        """Its corners are half a pixel each side of its centre."""
        archive = _registered(tmp_path)
        matrix = np.array(TEACH_POINTS_MATRIX).reshape(2, 3)

        with MzPeakReader(archive) as reader:
            alignment = reader.get_image_alignment()

        corners = alignment.cell_corners(1, 0)
        positions = [matrix @ np.array([cx, cy, 1.0]) for cx, cy in corners]
        np.testing.assert_allclose(
            positions, [[1.5, 0.5], [2.5, 0.5], [2.5, 1.5], [1.5, 1.5]], atol=1e-9
        )

    def test_the_regions_are_the_ones_the_scans_name(self, tmp_path):
        """One mapping per region, in raw raster indices."""
        archive = _registered(tmp_path)

        with MzPeakReader(archive) as reader:
            alignment = reader.get_image_alignment()

        bounds = {
            m.region_id: (
                m.raster_min_x,
                m.raster_max_x,
                m.raster_min_y,
                m.raster_max_y,
            )
            for m in alignment.region_mappings
        }
        assert bounds == {0: (101, 102, 11, 12), 1: (103, 103, 11, 12)}
        assert alignment.pos_to_region[(103, 12)] == 1

    @pytest.mark.parametrize(
        "quality, matrix",
        [
            ("assumed_full_extent", TEACH_POINTS_MATRIX),
            ("teach_points", TEACH_POINTS_MATRIX[:4]),
            ("teach_points", [0.0] * 6),
        ],
        ids=["not-registered", "short-matrix", "singular"],
    )
    def test_anything_else_is_carried_unaligned(self, tmp_path, quality, matrix):
        """A stretched image is no registration; a broken matrix is none."""
        archive = _registered(tmp_path, quality=quality, matrix=matrix)

        with MzPeakReader(archive) as reader:
            assert reader.get_image_alignment() is None
            assert reader.get_primary_optical_image_path() is None
            assert len(reader.get_optical_image_paths()) == 1

    def test_one_region_keeps_its_number_and_name(self, tmp_path):
        """The store's crop names the region as the archive does, not ``0``."""
        archive = _registered(
            tmp_path,
            regions=[3] * 6,
            listed=[_region(3, 6, (101, 103), (11, 12), name="vc_rugose_1")],
        )

        with MzPeakReader(archive) as reader:
            assert set(reader.get_region_map().values()) == {3}
            alignment = reader.get_image_alignment()

        (mapping,) = alignment.region_mappings
        assert (mapping.region_id, mapping.name) == (3, "vc_rugose_1")
        assert set(alignment.pos_to_region.values()) == {3}

    @pytest.mark.parametrize(
        "media_type, member",
        [("image/svg+xml", "images/image_0000.svg"), ("image/png", "notes.txt")],
        ids=["unreadable-format", "not-an-image-member"],
    )
    def test_an_image_the_store_cannot_hold_places_nothing(
        self, tmp_path, thyra_logs, media_type, member
    ):
        """The pixels are never placed on an image the store does not carry."""
        archive = _registered(tmp_path, media_type=media_type, member=member)
        index_name = INDEX_MEMBER
        if member == "notes.txt":
            # Listed, but not as an image.
            with zipfile.ZipFile(archive) as source:
                contents = {n: source.read(n) for n in source.namelist()}
            index = json.loads(contents[index_name])
            for entry in index["files"]:
                if entry["name"] == member:
                    entry["entity_type"] = "other"
            contents[index_name] = json.dumps(index).encode()
            with zipfile.ZipFile(archive, "w") as target:
                for name, payload in contents.items():
                    target.writestr(name, payload)

        with thyra_logs("thyra.readers.mzpeak", logging.WARNING) as records:
            with MzPeakReader(archive) as reader:
                assert reader.get_image_alignment() is None
                assert reader.get_primary_optical_image_path() is None
                assert reader.get_image_alignment() is None

        said = [r for r in records if "does not embed it" in r.getMessage()]
        assert len(said) == 1

    def test_the_metadata_document_states_the_teaching_points(self, tmp_path):
        """As the ``.d`` route states them from its ``.mis``."""
        archive = _registered(tmp_path)

        with MzPeakReader(archive) as reader:
            raw = reader.get_comprehensive_metadata().raw_metadata

        assert raw["mis_metadata"] == {
            "ImageFile": "IMG_0000.jpg",
            "teaching_points": [
                {"image": [1204.0, 778.0], "stage": [-22965.0, 15855.0]},
                {"image": [6706.0, 648.0], "stage": [21901.0, 17125.0]},
                {"image": [2550.0, 4990.0], "stage": [-11848.0, -18292.0]},
            ],
        }

    def test_an_unregistered_archive_states_no_teaching_points(self, tmp_path):
        """No sequence is made up for an image stretched over the run."""
        archive = _registered(tmp_path, quality="assumed_full_extent")

        with MzPeakReader(archive) as reader:
            raw = reader.get_comprehensive_metadata().raw_metadata

        assert "mis_metadata" not in raw


class TestOneTofBin:
    """Copies of one TOF bin are one m/z, though a bound replaced one of them."""

    def test_copies_a_bit_apart_are_summed(self):
        """A chunk's end is the converter's bound; the next copy is decoded."""
        bound = 1107.680159057425
        decoded = np.nextafter(bound, 2000.0)
        neighbour = bound * (1 + 2e-6)
        mzs = np.array([neighbour, decoded, bound])
        counts = np.array([4.0, 2.0, 2.0])

        out_mz, out_intensity = _on_tdf_scale(mzs, counts, 100.0)

        np.testing.assert_array_equal(out_mz, [bound, neighbour])
        np.testing.assert_array_equal(out_intensity, [4.0, 4.0])

    def test_numpress_copies_of_one_bin_are_summed(self):
        """F189: each scan's row is decoded on its own MS-Numpress scale (D82).

        A row reaching m/z 1000 stores a value to within 1000 / 2**32, so
        two scans' copies of a bin at m/z 50 differ by some 5e-9 relative.
        Tested apart, a peak whose sum passes the threshold was lost.
        """
        bin_mz = 50.123456789
        quantum = 1000.0 / 2**31
        copies = bin_mz + np.array([-0.5, 0.0, 0.25, 0.5]) * quantum
        neighbour = bin_mz * (1 + 2e-6)
        mzs = np.concatenate((copies, [neighbour]))
        counts = np.array([30.0, 30.0, 30.0, 30.0, 7.0])

        out_mz, out_intensity = _on_tdf_scale(mzs, counts, 100.0)

        assert out_mz.size == 2
        np.testing.assert_array_equal(out_intensity, [120.0, 7.0])
        assert out_mz[1] == neighbour


class TestOnlyMeasurementsAreSummed:
    """A point that is not a measurement joins no TOF bin (D61)."""

    @staticmethod
    def _read(tmp_path, mzs, counts, threshold=None):
        """The one pixel of a TDF archive at 100 ms, so counts stay as they are."""
        archive = build_mzpeak(
            tmp_path / "frame.mzpeak",
            [Spectrum(1, 1, mzs, counts)],
            file_metadata=TDF_SOURCE,
            scan_columns={"ion_injection_time": [100.0]},
        )
        with MzPeakReader(archive, intensity_threshold=threshold) as reader:
            ((_, out_mz, out_intensity),) = list(reader.iter_spectra())
        return out_mz, out_intensity

    def test_an_infinite_m_z_stays_out_of_the_largest_bin(self, tmp_path):
        """Its step from the largest finite m/z is no step, relative to inf."""
        out_mz, out_intensity = self._read(
            tmp_path, [100.0, 101.0, 102.0, np.inf], [2.0, 3.0, 5.0, 7.0]
        )

        np.testing.assert_array_equal(out_mz, [100.0, 101.0, 102.0, np.inf])
        np.testing.assert_array_equal(out_intensity, [2.0, 3.0, 5.0, 7.0])

    def test_minus_infinite_m_z_are_not_summed_together(self, tmp_path):
        """Each is the converter's to drop and count, one point at a time."""
        out_mz, out_intensity = self._read(
            tmp_path, [-np.inf, -np.inf, 100.0], [7.0, 4.0, 2.0]
        )

        np.testing.assert_array_equal(out_mz, [-np.inf, -np.inf, 100.0])
        np.testing.assert_array_equal(out_intensity, [7.0, 4.0, 2.0])

    @pytest.mark.parametrize("bad", [np.nan, np.inf, -4.0])
    def test_a_bad_count_leaves_its_bin_whole(self, tmp_path, bad):
        """A NaN would turn the bin into NaN; a negative would shrink it."""
        out_mz, out_intensity = self._read(
            tmp_path, [100.0, 100.0, 101.0], [9.0, bad, 3.0]
        )

        np.testing.assert_array_equal(out_mz, [100.0, 100.0, 101.0])
        np.testing.assert_array_equal(out_intensity, [9.0, bad, 3.0])

    def test_an_infinite_m_z_does_not_lift_a_bin_over_the_threshold(self, tmp_path):
        """The 7 at inf would carry the 5 at m/z 102 past a threshold of 6."""
        out_mz, out_intensity = self._read(
            tmp_path, [100.0, 101.0, 102.0, np.inf], [2.0, 3.0, 5.0, 7.0], 6.0
        )

        np.testing.assert_array_equal(out_mz, [np.inf])
        np.testing.assert_array_equal(out_intensity, [7.0])

    @pytest.mark.parametrize(
        "mzs, counts, expected, said",
        [
            # Stored inside the bin at m/z 102 before: 2, 3 and 12.
            (
                [100.0, 101.0, 102.0, np.inf],
                [2.0, 3.0, 5.0, 7.0],
                [2.0, 3.0, 5.0],
                "(0 non-finite, 0 negative, 1 at a non-finite m/z)",
            ),
            # Made the bin at m/z 100 NaN before, and the 2 was dropped with it.
            (
                [100.0, 100.0, 101.0, 102.0],
                [2.0, np.nan, 3.0, 5.0],
                [2.0, 3.0, 5.0],
                "(1 non-finite, 0 negative, 0 at a non-finite m/z)",
            ),
            # Shrank the bin at m/z 100 to 5 before, with no warning.
            (
                [100.0, 100.0, 101.0, 102.0],
                [9.0, -4.0, 3.0, 5.0],
                [9.0, 3.0, 5.0],
                "(0 non-finite, 1 negative, 0 at a non-finite m/z)",
            ),
        ],
        ids=["inf-mz", "nan-count", "negative-count"],
    )
    def test_the_store_holds_the_measured_bins(
        self, tmp_path, thyra_logs, mzs, counts, expected, said
    ):
        """The point that is not a measurement is dropped by the converter, said once."""
        import anndata

        from thyra.convert import convert_msi

        archive = build_mzpeak(
            tmp_path / "frame.mzpeak",
            [
                Spectrum(1, 1, mzs, counts),
                Spectrum(2, 1, [100.0, 101.0, 102.0], [1.0, 1.0, 1.0]),
            ],
            file_metadata=TDF_SOURCE,
            scan_columns={"ion_injection_time": [100.0, 100.0]},
        )
        out = tmp_path / "frame.zarr"

        with thyra_logs("thyra.converters", logging.WARNING) as records:
            assert convert_msi(
                str(archive),
                str(out),
                dataset_id="ds",
                resampling_config={
                    "method": "nearest_neighbor",
                    "min_mz": 99.5,
                    "max_mz": 102.5,
                    "target_bins": 7,
                },
            )

        table = anndata.read_zarr(str(out / "tables" / "ds_z0"))
        first = table.X.toarray()[0]
        stored = table.var["mz"].to_numpy()[first > 0]
        np.testing.assert_allclose(stored, [100.0, 101.0, 102.0])
        np.testing.assert_array_equal(first[first > 0], expected)
        dropped = [r.getMessage() for r in records if "Dropping" in r.getMessage()]
        assert len(dropped) == 1
        assert "1 of 4 in the first spectrum affected " + said in dropped[0]
