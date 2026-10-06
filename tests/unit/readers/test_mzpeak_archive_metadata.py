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
import zipfile

import numpy as np
import pytest

from tests.fixtures.mzpeak_builder import (
    INDEX_MEMBER,
    Spectrum,
    build_mzpeak,
    grid_spectra,
)
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


def _region(number, frames, x_index, y_index, name=None):
    """One ``bruker_maldi.regions`` entry, as mzpeak-convert writes it."""
    return {
        "frames": frames,
        "name": name,
        "raster_step_um": None,
        "region_number": number,
        "x_index": list(x_index),
        "y_index": list(y_index),
    }


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

    def test_one_region_is_none_as_from_the_d(self, tmp_path):
        """The timsTOF reader gives no region map for one region either."""
        archive = _two_regions(tmp_path, [_region(0, 8, (101, 105), (11, 12))])

        with MzPeakReader(archive) as reader:
            assert reader.get_region_map() is None
            assert reader.get_region_info() is None

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

        with pytest.raises(ConversionRefused, match="--region 1"):
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


def _with_images(archive, members):
    """Add image members to an archive, listed in its index as images."""
    with zipfile.ZipFile(archive) as source:
        contents = {name: source.read(name) for name in source.namelist()}
    index = json.loads(contents[INDEX_MEMBER])
    for name, payload in members.items():
        contents[name] = payload
        index["files"].append(
            {
                "name": name,
                "entity_type": "image",
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

    def test_names_that_differ_only_in_case_stay_apart(self, tmp_path):
        """The store lowercases an element name, so two must not meet there."""
        members = {"images/Scan.png": b"upper", "images/scan.png": b"lower"}
        archive = _with_images(
            build_mzpeak(tmp_path / "case.mzpeak", grid_spectra(2, 2)), members
        )

        with MzPeakReader(archive) as reader:
            paths = reader.get_optical_image_paths()
            names = [reader.get_optical_image_label(path).name for path in paths]
            contents = [path.read_bytes() for path in paths]

        assert len({name.lower() for name in names}) == 2
        assert contents == [b"upper", b"lower"]

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


# ----------------------------------------------------------------------
# Issue #429: what mzpeak-convert 0.17 states for a Bruker run
# ----------------------------------------------------------------------


def _region_parameter(number):
    """A scan's ``acquisition region`` parameter, as mzpeak-convert 0.17 writes it."""
    return [
        {
            "value": {
                "integer": number,
                "float": None,
                "string": None,
                "boolean": None,
            },
            "accession": None,
            "name": "acquisition region",
            "unit": None,
        }
    ]


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

    def test_one_named_region_is_none_as_from_the_d(self, tmp_path):
        """As for a ``.d`` of one region."""
        archive = _named_regions(tmp_path, [0] * 8, [])

        with MzPeakReader(archive) as reader:
            assert reader.get_region_map() is None


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
            assert reader.get_region_map() is None
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
