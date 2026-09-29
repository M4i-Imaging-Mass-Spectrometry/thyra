"""An imzML pixel size written under the vocabulary of before 2017.

IMS:1000046 was named "pixel size" and gave the AREA of a pixel until
imagingMS.obo commit 421481e of 2017-09-07, and IMS:1000047 was
"image shape". Public files still use that form: the term alone, no unit,
and a value whose square root times the pixel count is the extent.

Both imzML routes to a pixel size ask for IMS:1000046 AND IMS:1000047, so a
file of that form gives none and the conversion asks for ``--pixel-size``.
These tests hold that, on the fast path and on the full-XML fallback.
"""

from __future__ import annotations

import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace

import pytest

from thyra.metadata.extractors.imzml_extractor import ImzMLMetadataExtractor
from thyra.readers.imzml.imzml_reader import ImzMLReader

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


def _variant(tmp_path: Path, x_line: str, y_line: str = "") -> Path:
    """Copy the fixture pair with its pixel-size lines replaced."""
    text = FIXTURE.read_bytes().decode("utf-8")
    assert NEW_X in text and NEW_Y in text
    text = text.replace(NEW_X, x_line).replace(NEW_Y, y_line)

    imzml = tmp_path / "old_pixel_size.imzML"
    imzml.write_bytes(text.encode("utf-8"))
    shutil.copyfile(FIXTURE.with_suffix(".ibd"), imzml.with_suffix(".ibd"))
    return imzml


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

    @pytest.mark.filterwarnings("ignore:Accession IMS:UserWarning")
    def test_the_value_reaches_the_parser(self, tmp_path):
        """pyimzml hands the area over as ``pixel size x``, a bare number."""
        reader = ImzMLReader(_variant(tmp_path, OLD_AREA))
        try:
            reader.get_essential_metadata()
            assert reader.parser.imzmldict["pixel size x"] == pytest.approx(10000.0)
            assert "pixel size y" not in reader.parser.imzmldict
        finally:
            reader.close()

    @pytest.mark.filterwarnings("ignore:Accession IMS:UserWarning")
    def test_the_fast_path_takes_no_pixel_size(self, tmp_path):
        """Not 10000 um, and not the 100 um it means either."""
        reader = ImzMLReader(_variant(tmp_path, OLD_AREA))
        try:
            essential = reader.get_essential_metadata()
            assert essential.pixel_size is None
            assert essential.has_pixel_size is False
        finally:
            reader.close()

    @pytest.mark.filterwarnings("ignore:Accession IMS:UserWarning")
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
        """Reached by hand, the fallback still wants both terms."""
        assert _fallback(_variant(tmp_path, OLD_AREA)) is None

    def test_the_fallback_reads_a_declared_pair(self, tmp_path):
        """The guard above must not be a fallback that reads nothing."""
        imzml = _variant(tmp_path, NEW_X, NEW_Y)

        assert _fallback(imzml) == pytest.approx((4.40625, 4.40625))


class TestOldNameBesideANumber:
    """The one form both imzML routes would misread.

    IMS:1000046 under the old name, with a NUMBER on IMS:1000047. Under the
    old vocabulary that term is "image shape" and has no value type, so a
    file of this form is not expected. None was found: in 408 public files
    that give both terms, sampled one file per deposit and writer, every
    IMS:1000046 is named "pixel size x".
    """

    OLD_PAIR_Y = (
        '<cvParam cvRef="IMS" accession="IMS:1000047" name="image shape" '
        'value="10000" />'
    )

    @pytest.mark.xfail(
        strict=True,
        raises=AssertionError,
        reason="the imzML path reads a pair as lengths without looking at "
        "the name of IMS:1000046; no file of this form is known",
    )
    @pytest.mark.filterwarnings("ignore:Accession IMS:UserWarning")
    def test_the_fast_path(self, tmp_path):
        """Reads (10000, 10000) today."""
        reader = ImzMLReader(_variant(tmp_path, OLD_AREA, self.OLD_PAIR_Y))
        try:
            assert reader.get_essential_metadata().pixel_size is None
        finally:
            reader.close()

    @pytest.mark.xfail(
        strict=True,
        raises=AssertionError,
        reason="the full-XML fallback reads a pair as lengths without "
        "looking at the name of IMS:1000046; no real parser reaches it",
    )
    def test_the_fallback(self, tmp_path):
        """Reads (10000, 10000) when handed the document."""
        imzml = _variant(tmp_path, OLD_AREA, self.OLD_PAIR_Y)

        assert _fallback(imzml) is None
