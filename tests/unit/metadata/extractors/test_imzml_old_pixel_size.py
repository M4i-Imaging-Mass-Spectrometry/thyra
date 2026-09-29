"""An imzML pixel size written under the vocabulary of before 2017.

IMS:1000046 was named "pixel size" and gave the AREA of a pixel until
imagingMS.obo commit 421481e of 2017-09-07, and IMS:1000047 was
"image shape". Public files still use that form: the term alone, no unit,
and a value whose square root times the pixel count is the extent.

The imzML path does not read IMS:1000046 under that name. The file gives no
pixel size, the log says why, and the conversion asks for ``--pixel-size``.
These tests hold that on the fast path, in the preview, and on the full-XML
fallback, and hold that the current names are still read.
"""

from __future__ import annotations

import json
import logging
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace

import pytest

from thyra.convert import convert_msi
from thyra.metadata.extractors.imzml_extractor import ImzMLMetadataExtractor
from thyra.preview import preview_msi
from thyra.readers.imzml.imzml_reader import ImzMLReader

EXTRACTOR_LOGGER = "thyra.metadata.extractors.imzml_extractor"

FIXTURE = (
    Path(__file__).resolve().parents[3] / "data" / "fixtures" / "unit_nanometre.imzML"
)

#: The two pixel-size lines of the fixture, replaced by each variant below.
NEW_X = (
    '<cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)" '
    'value="4406.25" unitCvRef="UO" unitAccession="UO:0000018" '
    'unitName="nanometer" />'
)
NEW_Y = NEW_X.replace("IMS:1000046", "IMS:1000047").replace(
    "pixel size (x)", "pixel size y"
)

#: A 100 um pixel as an area, beside the extent of the fixture's two pixels.
OLD_AREA = (
    '<cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size" '
    'value="10000" />\n      '
    '<cvParam cvRef="IMS" accession="IMS:1000044" name="max dimension x" '
    'value="200" unitCvRef="UO" unitAccession="UO:0000017" '
    'unitName="micrometer" />'
)

#: IMS:1000047 with a number, under its old name and under its current one.
NUMBER_BESIDE_IT = {
    "image shape": (
        '<cvParam cvRef="IMS" accession="IMS:1000047" name="image shape" '
        'value="10000" />'
    ),
    "pixel size y": (
        '<cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y" '
        'value="10000" />'
    ),
}

pytestmark = pytest.mark.filterwarnings("ignore:Accession IMS:UserWarning")


def _variant(tmp_path: Path, x_line: str, y_line: str = "") -> Path:
    """Copy the fixture pair with its pixel-size lines replaced."""
    text = FIXTURE.read_bytes().decode("utf-8")
    assert NEW_X in text and NEW_Y in text
    text = text.replace(NEW_X, x_line).replace(NEW_Y, y_line)

    imzml = tmp_path / "old_pixel_size.imzML"
    imzml.write_bytes(text.encode("utf-8"))
    shutil.copyfile(FIXTURE.with_suffix(".ibd"), imzml.with_suffix(".ibd"))
    return imzml


def _pixel_size(imzml: Path):
    """The pixel size a conversion of this file would start from."""
    reader = ImzMLReader(imzml)
    try:
        return reader.get_essential_metadata().pixel_size
    finally:
        reader.close()


def _fallback(imzml: Path):
    """Run the full-XML fallback on the document itself.

    pyimzml's ``Metadata`` keeps no ``root``, so on a real parser the
    fallback returns ``None`` at its first guard for every file. Handing it
    the document is the only way to reach the part that reads the terms.
    """
    extractor = ImzMLMetadataExtractor.__new__(ImzMLMetadataExtractor)
    extractor.parser = SimpleNamespace(
        metadata=SimpleNamespace(root=ET.parse(imzml).getroot())
    )
    extractor.imzml_path = imzml
    return extractor._extract_pixel_size_from_xml()


class TestLoneOldPixelSize:
    """IMS:1000046 alone, under the old name: no pixel size is taken."""

    def test_the_value_reaches_the_parser(self, tmp_path):
        """pyimzml hands the area over as ``pixel size x``, a bare number."""
        reader = ImzMLReader(_variant(tmp_path, OLD_AREA))
        try:
            reader.get_essential_metadata()
            assert reader.parser.imzmldict["pixel size x"] == pytest.approx(10000.0)
            assert "pixel size y" not in reader.parser.imzmldict
        finally:
            reader.close()

    def test_the_parser_keeps_the_name_the_file_wrote(self, tmp_path):
        """pyimzml renames the term; the raw name is what the guard reads."""
        reader = ImzMLReader(_variant(tmp_path, OLD_AREA))
        try:
            reader.get_essential_metadata()
            params = [
                param
                for group in reader.parser.metadata.scan_settings.values()
                for param in group.cv_params
                if param[1] == "IMS:1000046"
            ]
        finally:
            reader.close()

        assert [(param[0], param[3]) for param in params] == [
            ("pixel size (x)", "pixel size")
        ]

    def test_the_fast_path_takes_no_pixel_size(self, tmp_path):
        """Not 10000 um, and not the 100 um it means either."""
        reader = ImzMLReader(_variant(tmp_path, OLD_AREA))
        try:
            essential = reader.get_essential_metadata()
            assert essential.pixel_size is None
            assert essential.has_pixel_size is False
        finally:
            reader.close()

    def test_the_log_says_why(self, tmp_path, thyra_logs):
        """ "Pixel size not found" alone would hide that the file has one."""
        imzml = _variant(tmp_path, OLD_AREA)

        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            assert _pixel_size(imzml) is None

        logged = "\n".join(records.messages)
        assert "old_pixel_size.imzML" in logged
        assert "area of a pixel" in logged
        assert "10000" in logged
        assert "--pixel-size" in logged

    def test_nothing_is_recorded_beside_the_data(self, tmp_path):
        """The fallback fills ``acquisition_params``; it must stay empty."""
        reader = ImzMLReader(_variant(tmp_path, OLD_AREA))
        try:
            params = reader.get_comprehensive_metadata().acquisition_params
            assert "pixel_size_x_um" not in params
            assert "pixel_size_y_um" not in params
        finally:
            reader.close()

    def test_the_fallback_takes_no_pixel_size(self, tmp_path):
        """Reached by hand, the fallback does not read it either."""
        assert _fallback(_variant(tmp_path, OLD_AREA)) is None


@pytest.mark.parametrize("name_of_y", sorted(NUMBER_BESIDE_IT))
class TestOldNameBesideANumber:
    """The old name on IMS:1000046, with a NUMBER on IMS:1000047.

    Both imzML routes read such a pair as two lengths and gave 10000 um by
    10000 um. Under the old vocabulary IMS:1000047 is "image shape" and has
    no value type, and no public file of this form was found. The name
    decides all the same, whatever IMS:1000047 is called.
    """

    def test_the_fast_path_takes_no_pixel_size(self, tmp_path, name_of_y):
        """It read (10000, 10000)."""
        imzml = _variant(tmp_path, OLD_AREA, NUMBER_BESIDE_IT[name_of_y])

        assert _pixel_size(imzml) is None

    def test_the_log_says_why(self, tmp_path, thyra_logs, name_of_y):
        """The reason is logged where the pixel size is dropped."""
        imzml = _variant(tmp_path, OLD_AREA, NUMBER_BESIDE_IT[name_of_y])

        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            _pixel_size(imzml)

        logged = "\n".join(records.messages)
        assert "area of a pixel" in logged
        assert "--pixel-size" in logged

    def test_the_preview_takes_no_pixel_size(self, tmp_path, name_of_y):
        """The header-only extractor shares the fast path."""
        imzml = _variant(tmp_path, OLD_AREA, NUMBER_BESIDE_IT[name_of_y])

        preview = preview_msi(imzml)

        assert preview.readable
        assert preview.pixel_size_um is None

    def test_the_fallback_takes_no_pixel_size(self, tmp_path, name_of_y):
        """It read (10000, 10000) when handed the document."""
        imzml = _variant(tmp_path, OLD_AREA, NUMBER_BESIDE_IT[name_of_y])

        assert _fallback(imzml) is None

    def test_the_conversion_asks_for_the_pixel_size(self, tmp_path, name_of_y):
        """``convert_msi`` returned ``True`` and wrote 10000 um."""
        imzml = _variant(tmp_path, OLD_AREA, NUMBER_BESIDE_IT[name_of_y])
        store = tmp_path / "refused.zarr"

        assert convert_msi(imzml, store, dataset_id="old") is False

        assert not store.exists()

    def test_a_stated_pixel_size_is_used(self, tmp_path, name_of_y):
        """``--pixel-size`` is the way through, and it is what gets stored."""
        imzml = _variant(tmp_path, OLD_AREA, NUMBER_BESIDE_IT[name_of_y])
        store = tmp_path / "stated.zarr"

        assert convert_msi(imzml, store, dataset_id="old", pixel_size_um=100.0)

        attributes = json.loads((store / "zarr.json").read_text(encoding="utf-8"))
        system = attributes["attributes"]["coordinate_systems"]["global"]
        assert system["pixel_size_um_x"] == pytest.approx(100.0)
        assert system["pixel_size_um_y"] == pytest.approx(100.0)


class TestCurrentNamesAreStillRead:
    """The guard is on one name. Every other spelling reads as before."""

    @pytest.mark.parametrize(
        "name_of_x",
        [
            "pixel size (x)",  # the vocabulary since 2017
            "pixel size x",  # what every real file at hand writes
        ],
    )
    def test_the_fast_path(self, tmp_path, thyra_logs, name_of_x):
        """4406.25 nanometre is 4.40625 um, and nothing is warned about."""
        imzml = _variant(tmp_path, NEW_X.replace("pixel size (x)", name_of_x), NEW_Y)

        with thyra_logs(EXTRACTOR_LOGGER, logging.WARNING) as records:
            pixel_size = _pixel_size(imzml)

        assert pixel_size == pytest.approx((4.40625, 4.40625))
        assert records.messages == []

    def test_the_fallback(self, tmp_path):
        """The guard must not be a fallback that reads nothing."""
        imzml = _variant(tmp_path, NEW_X, NEW_Y)

        assert _fallback(imzml) == pytest.approx((4.40625, 4.40625))
