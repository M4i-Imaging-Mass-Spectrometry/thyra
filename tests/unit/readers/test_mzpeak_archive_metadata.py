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
from thyra.errors import ConversionRefused
from thyra.readers.mzpeak import MzPeakReader
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
