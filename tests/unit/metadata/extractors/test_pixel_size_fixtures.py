"""The five pixel size test files, read on both paths.

``tests/data/fixtures/pixel_size_*.imzML`` hold one synthetic acquisition,
3 by 2 pixels of 50 um, five times. Each file states the pixel size in one
of the ways the headers of public imzML files do. The truth is known, because
the files are made by script, so every reading can be judged.

What is held here: a pixel size taken from one of these files is 50 um, and
where Thyra takes none the log asks for ``--pixel-size``. No path gives
another number.

The imzML path reads the files themselves. The mzPeak path reads an archive
whose scan settings are the geometry terms of the file, copied as a
converter copies them, built with the fixture builder.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pytest

from tests.fixtures.mzpeak_builder import build_mzpeak, grid_spectra
from thyra.preview import preview_msi
from thyra.readers.imzml.imzml_reader import ImzMLReader
from thyra.readers.mzpeak import MzPeakReader

IMZML_LOGGER = "thyra.metadata.extractors.imzml_extractor"
MZPEAK_LOGGER = "thyra.metadata.extractors.mzpeak_extractor"

FIXTURE_DIR = Path(__file__).resolve().parents[3] / "data" / "fixtures"
PREFIX = "pixel_size_"

CASES = [
    "unit_declared",
    "unit_absent",
    "area_old_name",
    "area_old_name_no_extent",
    "unit_contradiction",
]

COUNT_X, COUNT_Y = 3, 2
TRUE_SIZE = (50.0, 50.0)

MICROMETRE = ("UO:0000017", "micrometer")
NO_UNIT = (None, None)

#: accession, name, value, unit accession and unit name of IMS:1000044 to
#: IMS:1000047, as each header gives them.
EXTENT = [
    ("IMS:1000044", "max dimension x", "150", *MICROMETRE),
    ("IMS:1000045", "max dimension y", "100", *MICROMETRE),
]
GEOMETRY = {
    "unit_declared": [
        ("IMS:1000046", "pixel size (x)", "50", *MICROMETRE),
        ("IMS:1000047", "pixel size y", "50", *MICROMETRE),
        *EXTENT,
    ],
    "unit_absent": [
        ("IMS:1000046", "pixel size x", "50", *NO_UNIT),
        ("IMS:1000047", "pixel size y", "50", *NO_UNIT),
        *EXTENT,
    ],
    "area_old_name": [
        ("IMS:1000046", "pixel size", "2500", *NO_UNIT),
        *EXTENT,
    ],
    "area_old_name_no_extent": [
        ("IMS:1000046", "pixel size", "2500", *NO_UNIT),
    ],
    "unit_contradiction": [
        ("IMS:1000046", "pixel size", "50", "UO:0000015", "micrometer"),
    ],
}

#: What the imzML path takes from each file.
IMZML_PIXEL_SIZE = {
    "unit_declared": TRUE_SIZE,
    "unit_absent": TRUE_SIZE,
    "area_old_name": None,
    "area_old_name_no_extent": None,
    "unit_contradiction": None,
}

#: What the mzPeak path takes from an archive with the same terms. The area
#: beside an extent is the one case the two paths answer differently: the
#: mzPeak path tests the value against the grid, the imzML path asks.
MZPEAK_PIXEL_SIZE = {
    "unit_declared": TRUE_SIZE,
    "unit_absent": TRUE_SIZE,
    "area_old_name": TRUE_SIZE,
    "area_old_name_no_extent": None,
    "unit_contradiction": None,
}

TAG = re.compile(r"<cvParam\b[^>]*>")
ATTRIBUTE = re.compile(r'(\w+)="([^"]*)"')

pytestmark = pytest.mark.filterwarnings("ignore:Accession IMS:UserWarning")

Term = Tuple[str, str, str, Optional[str], Optional[str]]


def _path(case: str, suffix: str = "imzML") -> Path:
    return FIXTURE_DIR / f"{PREFIX}{case}.{suffix}"


def _header_terms(case: str) -> Dict[str, Term]:
    """Every cvParam in front of ``<run>``, by accession, first one wins."""
    head = _path(case).read_bytes().decode("ascii").split("<run", 1)[0]
    terms: Dict[str, Term] = {}
    for tag in TAG.findall(head):
        found = dict(ATTRIBUTE.findall(tag))
        terms.setdefault(
            found["accession"],
            (
                found["accession"],
                found["name"],
                found.get("value", ""),
                found.get("unitAccession"),
                found.get("unitName"),
            ),
        )
    return terms


def _imzml_pixel_size(case: str):
    """The pixel size a conversion of this file would start from."""
    reader = ImzMLReader(_path(case))
    try:
        return reader.get_essential_metadata().pixel_size
    finally:
        reader.close()


def _archive_parameters(case: str, unit_from: str = "accession") -> List[dict]:
    """The geometry terms of a file, spelled as an archive holds them.

    Name and value are copied. The value is a number in an archive. The unit
    is one accession or null, and where the accession and the name of a unit
    disagree a converter has to take one of the two.
    """
    names = {"micrometer": "UO:0000017"}
    terms = _header_terms(case)
    parameters = []
    for accession in sorted(terms):
        if not "IMS:1000042" <= accession <= "IMS:1000047":
            continue
        _, name, value, unit_accession, unit_name = terms[accession]
        unit = unit_accession
        if unit_from == "name" and unit_name is not None:
            unit = names[unit_name]
        parameters.append(
            {"name": name, "accession": accession, "value": int(value), "unit": unit}
        )
    return parameters


def _mzpeak_pixel_size(tmp_path, parameters: List[dict]):
    """Build an archive with these scan settings and read its pixel size."""
    archive = build_mzpeak(
        tmp_path / "pixel_size.mzpeak",
        grid_spectra(COUNT_X, COUNT_Y),
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


@pytest.mark.parametrize("case", CASES)
class TestTheFilesAreWhatTheReadmeSays:
    """Each file differs from the others in its geometry terms and no more."""

    def test_the_geometry_terms(self, case):
        """The header states the pixel size in the way the case is named for."""
        terms = _header_terms(case)
        geometry = [
            terms[accession]
            for accession in sorted(terms)
            if "IMS:1000044" <= accession <= "IMS:1000047"
        ]

        assert geometry == sorted(GEOMETRY[case])

    def test_the_grid(self, case):
        """3 by 2 pixels, declared and stored."""
        terms = _header_terms(case)

        assert terms["IMS:1000042"][2] == str(COUNT_X)
        assert terms["IMS:1000043"][2] == str(COUNT_Y)
        reader = ImzMLReader(_path(case))
        try:
            metadata = reader.get_essential_metadata()
            assert metadata.dimensions == (COUNT_X, COUNT_Y, 1)
            assert metadata.n_spectra == COUNT_X * COUNT_Y
        finally:
            reader.close()

    def test_the_scan_pattern(self, case):
        """All five declare a meander, so a dropped scan term shows."""
        assert _header_terms(case)["IMS:1000410"][1] == "meandering"

    def test_the_binary_file_is_the_one_the_header_names(self, case):
        """UUID and SHA-1 of the ``.ibd`` are the ones the header states."""
        terms = _header_terms(case)
        ibd = _path(case, "ibd").read_bytes()

        assert uuid.UUID(terms["IMS:1000080"][2]).bytes == ibd[:16]
        # Not a security use: the header of an imzML file names this digest.
        digest = hashlib.sha1(ibd, usedforsecurity=False).hexdigest()
        assert terms["IMS:1000091"][2] == digest.upper()

    def test_every_line_ends_in_crlf(self, case):
        """The line ends the files had when the converters were run on them."""
        data = _path(case).read_bytes()

        assert data.count(b"\r\n") > 0
        assert data.count(b"\n") == data.count(b"\r\n")

    def test_the_expected_outcome_is_recorded(self, case):
        """``pixel_size_expected.json`` gives the true size or none."""
        expected = json.loads((FIXTURE_DIR / f"{PREFIX}expected.json").read_bytes())

        assert sorted(expected) == sorted(CASES)
        assert expected[case]["pixel_size_um"] in ([50, 50], None)


@pytest.mark.parametrize("case", CASES)
class TestImzMLPath:
    """The file itself: 50 um, or no pixel size and a reason."""

    def test_the_pixel_size(self, case):
        """A declared pair is read, with or without a unit. A lone term is not."""
        pixel_size = _imzml_pixel_size(case)

        if IMZML_PIXEL_SIZE[case] is None:
            assert pixel_size is None
        else:
            assert pixel_size == pytest.approx(IMZML_PIXEL_SIZE[case])

    def test_the_log_asks_for_the_pixel_size(self, case, thyra_logs):
        """Where no pixel size is taken the log says why and what to pass."""
        with thyra_logs(IMZML_LOGGER, logging.WARNING) as records:
            _imzml_pixel_size(case)

        logged = "\n".join(records.messages)
        if IMZML_PIXEL_SIZE[case] is None:
            assert _path(case).name in logged
            assert "--pixel-size" in logged
        else:
            assert "--pixel-size" not in logged

    def test_the_preview_agrees(self, case):
        """The preview reads the head only and gives the same answer."""
        preview = preview_msi(_path(case))

        assert preview.readable
        assert preview.grid_dims == (COUNT_X, COUNT_Y)
        if IMZML_PIXEL_SIZE[case] is None:
            assert preview.pixel_size_um is None
        else:
            assert preview.pixel_size_um == pytest.approx(TRUE_SIZE[0])


@pytest.mark.parametrize("case", CASES)
class TestMzPeakPath:
    """An archive with the same terms: 50 um, or no pixel size and a reason."""

    def test_the_pixel_size(self, tmp_path, case):
        """The area beside an extent is tested against the grid and gives 50."""
        pixel_size = _mzpeak_pixel_size(tmp_path, _archive_parameters(case))

        if MZPEAK_PIXEL_SIZE[case] is None:
            assert pixel_size is None
        else:
            assert pixel_size == pytest.approx(MZPEAK_PIXEL_SIZE[case])

    def test_a_unit_taken_from_its_name(self, tmp_path, case):
        """A converter that takes the unit by its name changes no answer."""
        parameters = _archive_parameters(case, unit_from="name")

        pixel_size = _mzpeak_pixel_size(tmp_path, parameters)

        if MZPEAK_PIXEL_SIZE[case] is None:
            assert pixel_size is None
        else:
            assert pixel_size == pytest.approx(MZPEAK_PIXEL_SIZE[case])


class TestMzPeakPathSaysWhy:
    """Every reading that is not the number as written is logged."""

    @pytest.mark.parametrize(
        ("case", "unit_from", "reason"),
        [
            ("area_old_name_no_extent", "accession", "Pass --pixel-size"),
            ("unit_contradiction", "name", "Pass --pixel-size"),
            ("unit_contradiction", "accession", "unsupported unit 'UO:0000015'"),
        ],
    )
    def test_no_pixel_size_has_a_reason(
        self, tmp_path, thyra_logs, case, unit_from, reason
    ):
        """One value that cannot be tested, and the centimetre, are refused."""
        with thyra_logs(MZPEAK_LOGGER, logging.WARNING) as records:
            pixel_size = _mzpeak_pixel_size(
                tmp_path, _archive_parameters(case, unit_from)
            )

        assert pixel_size is None
        assert reason in "\n".join(records.messages)

    def test_the_log_names_the_area(self, tmp_path, thyra_logs):
        """2500 beside 3 pixels and 150 um: read as an area, side 50 um."""
        with thyra_logs(MZPEAK_LOGGER, logging.WARNING) as records:
            pixel_size = _mzpeak_pixel_size(
                tmp_path, _archive_parameters("area_old_name")
            )

        assert pixel_size == pytest.approx(TRUE_SIZE)
        logged = "\n".join(records.messages)
        assert "AREA" in logged
        assert "2500" in logged
        assert "50 um" in logged
