"""One pixel size in an mzPeak archive: a length, an area, or not believed.

IMS:1000046 was named "pixel size" and gave the AREA of a pixel until
imagingMS.obo commit 421481e of 2017-09-07. Since then it is
"pixel size (x)", a length, and IMS:1000047 -- "image shape" before -- is
"pixel size y". Files written under the old vocabulary are still published
and the reference converter copies the parameter into the archive unchanged,
so an archive can carry either meaning under the same accession.

The archives are built here with ``pyarrow`` and ``zipfile``, through the
fixture builder. The scan settings are spelled as the reference converter
writes them: ``name``, ``accession``, ``value`` as a number, ``unit`` as a
unit accession or null.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import pytest

from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra.readers.mzpeak import MzPeakReader

EXTRACTOR_LOGGER = "thyra.metadata.extractors.mzpeak_extractor"

MICROMETRE = "UO:0000017"
MILLIMETRE = "UO:0000016"
CENTIMETRE = "UO:0000015"


def _parameter(accession: str, name: str, value, unit: Optional[str] = None) -> dict:
    """One scan-settings parameter, as the reference converter spells it."""
    return {"name": name, "accession": accession, "value": value, "unit": unit}


def old_pixel_size(value, unit: Optional[str] = None) -> dict:
    """IMS:1000046 under the name it had when it gave an area."""
    return _parameter("IMS:1000046", "pixel size", value, unit)


def pixel_size_x(value, unit: Optional[str] = None) -> dict:
    """IMS:1000046 under its name since 2017."""
    return _parameter("IMS:1000046", "pixel size (x)", value, unit)


def pixel_size_y(value, unit: Optional[str] = None) -> dict:
    """IMS:1000047 under its name since 2017."""
    return _parameter("IMS:1000047", "pixel size y", value, unit)


def count_x(value) -> dict:
    """Declared pixel count along x."""
    return _parameter("IMS:1000042", "max count of pixels x", value)


def count_y(value) -> dict:
    """Declared pixel count along y."""
    return _parameter("IMS:1000043", "max count of pixels y", value)


def extent_x(value, unit: Optional[str] = MICROMETRE) -> dict:
    """Declared extent along x."""
    return _parameter("IMS:1000044", "max dimension x", value, unit)


def extent_y(value, unit: Optional[str] = MICROMETRE) -> dict:
    """Declared extent along y."""
    return _parameter("IMS:1000045", "max dimension y", value, unit)


def _pixel_size(tmp_path, parameters: List[dict]):
    """Build an archive with these scan settings and read its pixel size."""
    archive = build_mzpeak(
        tmp_path / "pixel_size.mzpeak",
        grid_spectra(5, 3),
        pixel_size=None,
        footer_metadata=False,
        index_metadata={
            "scan_settings_list": [
                {
                    "id": "scanSettings1",
                    "source_file_refs": [],
                    "targets": [],
                    "parameters": parameters,
                }
            ]
        },
    )
    with MzPeakReader(archive) as reader:
        return reader.get_essential_metadata().pixel_size


def _warnings(records) -> str:
    """Every WARNING the extractor logged, as one string."""
    return "\n".join(records.messages)


class TestAreaUnderTheOldName:
    """A value the grid shows to be an area gives its square root."""

    def test_the_reported_archive(self, tmp_path, thyra_logs):
        """Pixel size 10000, 5 pixels on x, extent 500 um: a 100 um pixel.

        This is the archive the reference converter makes from an imzML
        written under the old vocabulary. Thyra 4.2.0 wrote 10000 um.
        """
        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            pixel_size = _pixel_size(
                tmp_path,
                [old_pixel_size(10000), extent_x(500), count_x(5), count_y(3)],
            )

        assert pixel_size == pytest.approx((100.0, 100.0))
        logged = _warnings(records)
        assert "AREA" in logged
        assert "10000" in logged
        assert "100 um" in logged

    @pytest.mark.parametrize(
        ("value", "count", "extent", "side"),
        [
            (10000, 174, 17400, 100.0),
            (100, 2704, 27040, 10.0),
        ],
    )
    def test_values_found_in_public_files(self, tmp_path, value, count, extent, side):
        """Two of the public headers in which the area reading is exact."""
        pixel_size = _pixel_size(
            tmp_path, [old_pixel_size(value), count_x(count), extent_x(extent)]
        )

        assert pixel_size == pytest.approx((side, side))

    def test_both_axes_agree(self, tmp_path):
        """When the file declares both axes, both confirm the side."""
        pixel_size = _pixel_size(
            tmp_path,
            [
                old_pixel_size(10000),
                count_x(5),
                extent_x(500),
                count_y(3),
                extent_y(300),
            ],
        )

        assert pixel_size == pytest.approx((100.0, 100.0))

    def test_the_new_name_does_not_shield_an_area(self, tmp_path):
        """The test is on the numbers. A lone term is tested under any name."""
        pixel_size = _pixel_size(
            tmp_path, [pixel_size_x(10000), count_x(5), extent_x(500)]
        )

        assert pixel_size == pytest.approx((100.0, 100.0))

    def test_image_shape_beside_it_is_not_a_pixel_size(self, tmp_path):
        """Under the old vocabulary IMS:1000047 is "image shape"."""
        shape = {
            "name": "image shape",
            "accession": "IMS:1000047",
            "value": None,
            "unit": None,
        }
        pixel_size = _pixel_size(
            tmp_path, [old_pixel_size(10000), shape, count_x(5), extent_x(500)]
        )

        assert pixel_size == pytest.approx((100.0, 100.0))

    def test_a_unit_scales_the_side(self, tmp_path):
        """0.01 in millimetre, read as an area, is a side of 0.1 mm."""
        pixel_size = _pixel_size(
            tmp_path, [old_pixel_size(0.01, MILLIMETRE), count_x(5), extent_x(500)]
        )

        assert pixel_size == pytest.approx((100.0, 100.0))


class TestLengthConfirmed:
    """A lone value the grid shows to be a length is taken as written."""

    def test_a_lone_length(self, tmp_path, thyra_logs):
        """Value x count is the extent. Nothing to warn about."""
        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            pixel_size = _pixel_size(
                tmp_path, [pixel_size_x(100), count_x(5), extent_x(500)]
            )

        assert pixel_size == pytest.approx((100.0, 100.0))
        assert _warnings(records) == ""

    def test_the_old_name_on_a_length(self, tmp_path):
        """A writer may keep the old name and mean a length."""
        pixel_size = _pixel_size(
            tmp_path, [old_pixel_size(100), count_x(5), extent_x(500)]
        )

        assert pixel_size == pytest.approx((100.0, 100.0))

    def test_the_old_name_on_a_length_keeps_the_declared_y(self, tmp_path):
        """Shown to be a length, the file speaks the new vocabulary."""
        pixel_size = _pixel_size(
            tmp_path,
            [
                old_pixel_size(100, MICROMETRE),
                pixel_size_y(50, MICROMETRE),
                count_x(5),
                extent_x(500),
            ],
        )

        assert pixel_size == pytest.approx((100.0, 50.0))

    def test_units_are_folded_before_the_test(self, tmp_path):
        """0.1 mm against an extent in micrometres."""
        pixel_size = _pixel_size(
            tmp_path, [pixel_size_x(0.1, MILLIMETRE), count_x(5), extent_x(500)]
        )

        assert pixel_size == pytest.approx((100.0, 100.0))

    def test_a_value_of_one_is_the_same_either_way(self, tmp_path, thyra_logs):
        """The square root of 1 is 1, so the two readings agree."""
        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            pixel_size = _pixel_size(
                tmp_path, [old_pixel_size(1), count_x(5), extent_x(5, None)]
            )

        assert pixel_size == pytest.approx((1.0, 1.0))
        assert "AREA" not in _warnings(records)

    def test_a_rounded_extent_passes(self, tmp_path):
        """9.9248 x 100 is 992.48; the file says 992."""
        pixel_size = _pixel_size(
            tmp_path, [pixel_size_x(9.9248), count_x(100), extent_x(992)]
        )

        assert pixel_size == pytest.approx((9.9248, 9.9248))


class TestNotBelieved:
    """No test, or no reading that fits: no pixel size, and a reason."""

    def _refused(self, tmp_path, thyra_logs, parameters) -> str:
        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            pixel_size = _pixel_size(tmp_path, parameters)
        assert pixel_size is None
        return _warnings(records)

    def test_no_count_and_no_extent(self, tmp_path, thyra_logs):
        """The value alone is what Thyra 4.2.0 turned into a square pixel."""
        logged = self._refused(tmp_path, thyra_logs, [old_pixel_size(10000)])

        assert "IMS:1000046 = 10000" in logged
        assert "--pixel-size" in logged

    def test_a_count_without_an_extent(self, tmp_path, thyra_logs):
        """31 public files are like this. Nothing to test the value against."""
        logged = self._refused(
            tmp_path, thyra_logs, [old_pixel_size(100), count_x(256), count_y(128)]
        )

        assert "--pixel-size" in logged

    def test_neither_reading_fits(self, tmp_path, thyra_logs):
        """10000 x 5 is not 700, and neither is 100 x 5."""
        logged = self._refused(
            tmp_path, thyra_logs, [old_pixel_size(10000), count_x(5), extent_x(700)]
        )

        assert "fit neither reading" in logged

    def test_the_axes_disagree(self, tmp_path, thyra_logs):
        """One number describes both axes, so both have to confirm it."""
        self._refused(
            tmp_path,
            thyra_logs,
            [
                old_pixel_size(10000),
                count_x(5),
                extent_x(500),
                count_y(3),
                extent_y(900),
            ],
        )

    def test_an_extent_outside_the_tolerance(self, tmp_path, thyra_logs):
        """2% off is not a rounded extent."""
        self._refused(
            tmp_path, thyra_logs, [pixel_size_x(100), count_x(5), extent_x(510)]
        )

    @pytest.mark.parametrize("value", [0, -100])
    def test_a_value_that_is_no_size(self, tmp_path, thyra_logs, value):
        """Zero and negative values have no reading at all."""
        self._refused(
            tmp_path, thyra_logs, [old_pixel_size(value), count_x(5), extent_x(500)]
        )

    def test_a_lone_y(self, tmp_path, thyra_logs):
        """IMS:1000047 alone was "image shape" before 2017."""
        logged = self._refused(
            tmp_path, thyra_logs, [pixel_size_y(100, MICROMETRE), count_y(3)]
        )

        assert "image shape" in logged

    def test_a_lone_value_in_an_unknown_unit(self, tmp_path, thyra_logs):
        """The unit is refused before any test is made."""
        logged = self._refused(
            tmp_path,
            thyra_logs,
            [old_pixel_size(100, CENTIMETRE), count_x(5), extent_x(500)],
        )

        assert CENTIMETRE in logged


class TestDeclaredPair:
    """Both terms under their current names are taken as written."""

    def test_a_pair_with_a_unit_is_not_tested(self, tmp_path, thyra_logs):
        """The declared extent does not overrule a declared pair.

        A real file declares an extent of 500 beside 128 pixels of
        4.40625 um. Its pixel size is still what it says.
        """
        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            pixel_size = _pixel_size(
                tmp_path,
                [
                    pixel_size_x(30, MICROMETRE),
                    pixel_size_y(40, MICROMETRE),
                    count_x(5),
                    extent_x(9999),
                ],
            )

        assert pixel_size == pytest.approx((30.0, 40.0))
        assert _warnings(records) == ""

    def test_a_pair_with_no_unit_is_micrometres(self, tmp_path):
        """As on the imzML path, so a file and its archive agree."""
        pixel_size = _pixel_size(tmp_path, [pixel_size_x(25), pixel_size_y(25)])

        assert pixel_size == pytest.approx((25.0, 25.0))

    def test_centimetre_is_refused(self, tmp_path, thyra_logs):
        """The public files that declare centimetre name it "micrometer".

        A factor of 10000 would be wrong for every one of them, and the
        number as written would trust a unit the file contradicts.
        """
        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            pixel_size = _pixel_size(
                tmp_path,
                [pixel_size_x(100, CENTIMETRE), pixel_size_y(100, CENTIMETRE)],
            )

        assert pixel_size is None
        assert CENTIMETRE in _warnings(records)

    def test_one_unreadable_axis_is_not_copied_from_the_other(
        self, tmp_path, thyra_logs
    ):
        """Thyra 4.2.0 filled the missing axis in from the readable one."""
        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            pixel_size = _pixel_size(
                tmp_path,
                [pixel_size_x(25, MICROMETRE), pixel_size_y("wide", MICROMETRE)],
            )

        assert pixel_size is None
        assert "--pixel-size" in _warnings(records)
