"""Regenerate the five pixel size test files in this directory.

Every file holds the same synthetic acquisition: 3 by 2 pixels, each 50
micrometre wide, scanned as a meander. Only the way the header states the
pixel size differs. The five ways are the ones found in the headers of public
imzML files, see ``README.md`` beside this script.

Unlike the hand-authored corpus these files come from pyimzml's own
``ImzMLWriter``. That is on purpose: the subject is what a reader or a
converter does with the geometry terms, and the rest of the file should be as
plain as a file can be. The geometry lines are put into the header as text,
because the writer has no argument for them.

Run it from the repository root::

    python tests/data/fixtures/build_pixel_size_fixtures.py

The output is the same byte for byte on every run, so a clean run leaves
``git status`` clean. Three things make it so. The UUID of each file is fixed
below, where the writer would draw a new one. The writer names the run after
the path it is given, so each file is written under its bare name in a
scratch directory. And the writer ends lines as the platform does, so the
line ends are set to CRLF afterwards, which is what the files had when they
were first made.

The committed bytes are the fixture. This script is provenance, and no test
invokes it.
"""

from __future__ import annotations

import json
import os
import tempfile
import uuid as uuid_module
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple
from unittest import mock

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

FIXTURE_DIR = Path(__file__).resolve().parent

#: Every file of the set starts with this, so the set sorts together.
PREFIX = "pixel_size_"

COUNT_X, COUNT_Y = 3, 2
SIDE_UM = 50

#: The geometry lines go in front of the first pixel count line of the header.
MARKER = '<cvParam cvRef="IMS" accession="IMS:1000042"'
INDENT = "\n      "

MICROMETRE = 'unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"'

#: The unit accession of the centimetre under the name of the micrometre.
CONTRADICTION = 'unitCvRef="UO" unitAccession="UO:0000015" unitName="micrometer"'

EXTENT = [
    '<cvParam cvRef="IMS" accession="IMS:1000044" name="max dimension x" '
    f'value="{COUNT_X * SIDE_UM}" {MICROMETRE}/>',
    '<cvParam cvRef="IMS" accession="IMS:1000045" name="max dimension y" '
    f'value="{COUNT_Y * SIDE_UM}" {MICROMETRE}/>',
]


class Case(NamedTuple):
    """One way of stating the pixel size."""

    uuid: str
    lines: List[str]
    pixel_size_um: Optional[Tuple[int, int]]
    pixel_size_source: str


CASES: Dict[str, Case] = {
    "unit_declared": Case(
        uuid="2487DB05-6594-49DA-8B16-D4C1FA24CA57",
        lines=[
            '<cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)" '
            f'value="{SIDE_UM}" {MICROMETRE}/>',
            '<cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y" '
            f'value="{SIDE_UM}" {MICROMETRE}/>',
            *EXTENT,
        ],
        pixel_size_um=(SIDE_UM, SIDE_UM),
        pixel_size_source="declared",
    ),
    "unit_absent": Case(
        uuid="8F539171-3663-4176-93BF-AC5897310C7C",
        lines=[
            '<cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size x" '
            f'value="{SIDE_UM}"/>',
            '<cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y" '
            f'value="{SIDE_UM}"/>',
            *EXTENT,
        ],
        pixel_size_um=(SIDE_UM, SIDE_UM),
        pixel_size_source="unit_assumed",
    ),
    "area_old_name": Case(
        uuid="23217EBF-8C63-4185-B420-05DACE3D3DC1",
        lines=[
            '<cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size" '
            f'value="{SIDE_UM * SIDE_UM}"/>',
            *EXTENT,
        ],
        pixel_size_um=(SIDE_UM, SIDE_UM),
        pixel_size_source="derived_from_area",
    ),
    "area_old_name_no_extent": Case(
        uuid="862CEA61-0E22-464D-A179-5EEE11414B42",
        lines=[
            '<cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size" '
            f'value="{SIDE_UM * SIDE_UM}"/>',
        ],
        pixel_size_um=None,
        pixel_size_source="unknown",
    ),
    "unit_contradiction": Case(
        uuid="FF7158DC-6989-45E7-B947-77D7456A3115",
        lines=[
            '<cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size" '
            f'value="{SIDE_UM}" {CONTRADICTION}/>',
        ],
        pixel_size_um=None,
        pixel_size_source="unknown",
    ),
}


def _spectra() -> List[Tuple[Tuple[int, int, int], np.ndarray, np.ndarray]]:
    """Six profile spectra with one peak at m/z 184.07, on one shared axis."""
    rng = np.random.default_rng(20260929)
    spectra = []
    for y in range(1, COUNT_Y + 1):
        for x in range(1, COUNT_X + 1):
            mzs = np.arange(183.8, 184.4, 0.01)
            peak = np.exp(-0.5 * ((mzs - 184.07) / 0.03) ** 2)
            intensities = np.round(rng.gamma(2.0, 500.0) * peak + 1.0)
            spectra.append(((x, y, 1), mzs, intensities.astype(np.float32)))
    return spectra


def _write_pair(name: str, case: Case, scratch: Path) -> Tuple[bytes, bytes]:
    """Write one file with pyimzml and hand back its imzML and .ibd bytes."""
    fixed = uuid_module.UUID(case.uuid)
    here = os.getcwd()
    os.chdir(scratch)  # the writer names the run after the path it is given
    try:
        with mock.patch.object(uuid_module, "uuid4", return_value=fixed):
            with ImzMLWriter(
                f"{name}.imzML",
                mode="continuous",
                spec_type="profile",
                scan_direction="top_down",
                line_scan_direction="line_left_right",
                scan_pattern="meandering",
                scan_type="horizontal_line",
            ) as writer:
                for position, mzs, intensities in _spectra():
                    writer.addSpectrum(mzs, intensities, position)
    finally:
        os.chdir(here)

    text = (scratch / f"{name}.imzML").read_bytes().decode("ascii")
    text = text.replace("\r\n", "\n")
    assert text.count(MARKER) == 1
    lines = "".join(line + INDENT for line in case.lines)
    text = text.replace(MARKER, lines + MARKER)
    imzml = text.replace("\n", "\r\n").encode("ascii")
    return imzml, (scratch / f"{name}.ibd").read_bytes()


def main() -> None:
    """Write the five pairs and ``pixel_size_expected.json``."""
    expected = {}
    with tempfile.TemporaryDirectory() as scratch:
        for name, case in CASES.items():
            imzml, ibd = _write_pair(name, case, Path(scratch))
            (FIXTURE_DIR / f"{PREFIX}{name}.imzML").write_bytes(imzml)
            (FIXTURE_DIR / f"{PREFIX}{name}.ibd").write_bytes(ibd)
            expected[name] = {
                "pixel_size_um": (
                    None if case.pixel_size_um is None else list(case.pixel_size_um)
                ),
                "pixel_size_source": case.pixel_size_source,
            }
            print(f"wrote {PREFIX}{name}: {len(imzml)} + {len(ibd)} bytes")

    text = json.dumps(expected, indent=2) + "\n"
    (FIXTURE_DIR / f"{PREFIX}expected.json").write_bytes(text.encode("ascii"))


if __name__ == "__main__":
    main()
