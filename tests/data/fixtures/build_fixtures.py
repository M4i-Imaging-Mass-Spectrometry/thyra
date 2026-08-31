"""Regenerate the hand-authored imzML corpus in this directory.

Every imzML the rest of the suite reads was produced by pyimzml's own
``ImzMLWriter``, so the parser and the writer agree on each other's mistakes and
nine structural features of real vendor files are unreachable (audit #15). The
files here exist to break that loop: the XML below is written out literally,
character by character, and the ``.ibd`` is packed by hand, so nothing pyimzml
does can influence what a fixture looks like.

Run it from the repository root with the worktree first on the path::

    PYTHONPATH=. python tests/data/fixtures/build_fixtures.py

Output is byte-for-byte reproducible, so a clean run leaves ``git status``
clean. The committed bytes are the fixture -- this script is provenance, not a
build step, and no test invokes it. If you change it, commit the regenerated
files in the same commit and re-read the blobs (see this directory's README for
why reading the worktree file is not enough).
"""

from __future__ import annotations

import uuid as uuid_module
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np

FIXTURE_DIR = Path(__file__).resolve().parent

# Byte 0..15 of an .ibd is the binary form of the UUID the header declares in
# IMS:1000080. pyimzml never checks it; real files always carry it.
UUID_HEADER_BYTES = 16

Coordinate = Tuple[int, int]


def _pack_ibd(
    uuid_text: str,
    spectra: Sequence[Tuple[np.ndarray, np.ndarray]],
) -> Tuple[bytes, List[Tuple[int, int, int, int]]]:
    """Pack spectra into an .ibd payload and report where each array landed.

    Args:
        uuid_text: The UUID as it appears in IMS:1000080, braces optional.
        spectra: One ``(mzs, intensities)`` pair per spectrum, already in the
            dtype the header will declare.

    Returns:
        The ``.ibd`` bytes, and one ``(mz_offset, mz_length, intensity_offset,
        intensity_length)`` tuple per spectrum. Lengths are element counts, the
        way IMS:1000103 expresses them.
    """
    blob = bytearray(uuid_module.UUID(uuid_text).bytes)
    assert len(blob) == UUID_HEADER_BYTES

    placements: List[Tuple[int, int, int, int]] = []
    for mzs, intensities in spectra:
        mz_offset = len(blob)
        blob += mzs.tobytes()
        intensity_offset = len(blob)
        blob += intensities.tobytes()
        placements.append(
            (mz_offset, int(mzs.size), intensity_offset, int(intensities.size))
        )
    return bytes(blob), placements


def _spectrum_xml(
    index: int,
    coordinate: Coordinate,
    placement: Tuple[int, int, int, int],
    itemsize: int,
    spectrum_group: str = "spectrum1",
) -> List[str]:
    """Emit one ``<spectrum>`` element, with 0-based ``@index`` as vendors write it.

    Args:
        index: The 0-based spectrum index. pyimzml's writer emits 1-based
            indices; all three real files are 0-based.
        coordinate: The 1-based ``(x, y)`` grid position.
        placement: The tuple ``_pack_ibd`` returned for this spectrum.
        itemsize: Bytes per element, used for IMS:1000104.
        spectrum_group: The referenceableParamGroup carrying the spectrum type.

    Returns:
        The element's lines, without line terminators.
    """
    mz_offset, mz_length, intensity_offset, intensity_length = placement
    x, y = coordinate
    return [
        f'<spectrum id="Scan={index + 1}" defaultArrayLength="0" index="{index}">'
        f'<referenceableParamGroupRef ref="{spectrum_group}" />',
        '<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795"'
        ' name="no combination" value="" />',
        '<scan><referenceableParamGroupRef ref="scan1" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000050" name="position x"'
        f' value="{x}" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000051" name="position y"'
        f' value="{y}" />',
        "</scan>",
        "</scanList>",
        '<binaryDataArrayList count="2"><binaryDataArray encodedLength="0">'
        '<referenceableParamGroupRef ref="mzArray" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000103"'
        f' name="external array length" value="{mz_length}" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000102" name="external offset"'
        f' value="{mz_offset}" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000104"'
        f' name="external encoded length" value="{mz_length * itemsize}" />',
        "</binaryDataArray>",
        '<binaryDataArray encodedLength="0">'
        '<referenceableParamGroupRef ref="intensityArray" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000103"'
        f' name="external array length" value="{intensity_length}" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000102" name="external offset"'
        f' value="{intensity_offset}" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000104"'
        f' name="external encoded length"'
        f' value="{intensity_length * itemsize}" />',
        "</binaryDataArray>",
        "</binaryDataArrayList>",
        "</spectrum>",
    ]


def _make_spectra(
    n: int, dtype: np.dtype, n_peaks: int = 5
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Build n small ascending processed-mode spectra.

    Each spectrum gets its own m/z values so the file is genuinely processed
    mode; the offsets differ per spectrum, which is what makes a truncated or
    past-EOF .ibd reachable from these fixtures later.
    """
    spectra = []
    for i in range(n):
        mzs = np.linspace(100.0 + i, 500.0 + i, n_peaks).astype(dtype)
        intensities = (np.arange(1, n_peaks + 1) * (i + 1) * 10.0).astype(dtype)
        spectra.append((mzs, intensities))
    return spectra


def _write(stem: str, lines: Sequence[str], newline: bytes, encoding: str) -> Path:
    """Write the imzML with the exact bytes asked for, bypassing text mode."""
    path = FIXTURE_DIR / f"{stem}.imzML"
    payload = newline.join(line.encode(encoding) for line in lines) + newline
    path.write_bytes(payload)
    return path


# --------------------------------------------------------------------------
# iontof_sparse
# --------------------------------------------------------------------------

IONTOF_UUID = "FC37F303-A9C0-4CD3-A28E-1D18E523C269"

# Six acquired positions out of a declared 4x4 grid. Real IONTOF acquisitions
# are a sparse subset of the declared raster, which is what makes the empty-row
# handling of the different write routes observable at all.
IONTOF_COORDINATES: List[Coordinate] = [
    (1, 1),
    (2, 1),
    (4, 2),
    (1, 3),
    (3, 3),
    (4, 4),
]


def build_iontof_sparse() -> Path:
    """The whole bellini class in one file.

    ISO-8859-1, CRLF, unindented, no ``mzML/@version``, no IMS:1000052, three
    misspelled ontology names, and ``max dimension`` contradicting
    ``pixel size``. Every one of those is a real trait of the shipped
    ``bellini.imzML`` and none of them is reachable through ``ImzMLWriter``.

    The unindented-CRLF layout is also the exact shape that makes lxml raise
    ``XMLSyntaxError: xmlSAX2Characters`` -- see the ``parse_lib`` rationale in
    ``thyra/readers/imzml/imzml_reader.py``.
    """
    spectra = _make_spectra(len(IONTOF_COORDINATES), np.dtype(np.float64))
    ibd, placements = _pack_ibd(IONTOF_UUID, spectra)
    (FIXTURE_DIR / "iontof_sparse.ibd").write_bytes(ibd)

    lines = [
        '<?xml version="1.0" encoding="ISO-8859-1"?>'
        # No version attribute. mzML 1.1.0 requires one; IONTOF omits it and
        # pyimzml never looks.
        '<mzML xmlns="http://psi.hupo.org/ms/mzml">'
        '<cvList count="3"><cv id="MS" fullName="Proteomics Standards Initiative'
        ' Mass Spectrometry Ontology" version="1.3.1"'
        ' URI="http://psidev.info/ms/mzML/psi-ms.obo" />',
        '<cv id="UO" fullName="Unit Ontology" version="1.15"'
        ' URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo" />',
        '<cv id="IMS" fullName="Imaging MS Ontology" version="0.9.1"'
        ' URI="http://www.maldi-msi.org/download/imzml/imagingMS.obo" />',
        "</cvList>",
        '<fileDescription><fileContent><cvParam cvRef="MS" accession="MS:1000128"'
        ' name="profile spectrum" value="" />',
        '<cvParam cvRef="IMS" accession="IMS:1000031" name="processed" value="" />',
        '<cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value="" />',
        f'<cvParam cvRef="IMS" accession="IMS:1000080"'
        f' name="universally unique identifier" value="{{{IONTOF_UUID}}}" />',
        "</fileContent>",
        # sourceFileList: present on all three real files, absent from every
        # ImzMLWriter fixture.
        '<sourceFileList count="1"><sourceFile id="sf1" name="2019_03_0025" />',
        "</sourceFileList>",
        "</fileDescription>",
        '<referenceableParamGroupList count="4">'
        '<referenceableParamGroup id="mzArray">'
        '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value="" />',
        '<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" />',
        '<cvParam cvRef="IMS" accession="IMS:1000101" name="external data"'
        ' value="true" />',
        '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value="" />',
        "</referenceableParamGroup>",
        '<referenceableParamGroup id="intensityArray">'
        '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value="" />',
        '<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" />',
        '<cvParam cvRef="IMS" accession="IMS:1000101" name="external data"'
        ' value="true" />',
        '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value="" />',
        "</referenceableParamGroup>",
        '<referenceableParamGroup id="scan1"><cvParam cvRef="MS"'
        ' accession="MS:1000093" name="increasing m/z scan" value="" />',
        '<cvParam cvRef="MS" accession="MS:1000095" name="linear" value="" />',
        "</referenceableParamGroup>",
        '<referenceableParamGroup id="spectrum1"><cvParam cvRef="MS"'
        ' accession="MS:1000128" name="profile spectrum" value="" />',
        '<cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value="" />',
        "</referenceableParamGroup>",
        "</referenceableParamGroupList>",
        '<sampleList count="1"><sample id="sample1" name="Unnamed Sample" />',
        "</sampleList>",
        '<softwareList count="2"><software id="IONTOF SurfaceLab 7.5" />',
        '<software id="IONTOF imzML exporter"><cvParam cvRef="MS"'
        ' accession="MS:1000799" name="IONTOF SurfaceLab 7.5" value="" />',
        "</software>",
        "</softwareList>",
        # Misspelled names below. The ontology calls IMS:1000042/43
        # "max count of pixels x/y" and IMS:1000046 "pixel size (x)"; IONTOF
        # writes the singular and the unparenthesised form, and pyimzml emits a
        # UserWarning per correction. IMS:1000047's canonical name really is
        # "pixel size y", so only three of these four warn.
        #
        # The contradiction: 4 pixels x 4.40625 um = 17.625 um of raster, but
        # max dimension says 15. Deriving um/pixel from max dimension gives
        # 3.75, reading pixel size gives 4.40625 -- 17.5% apart with no unit on
        # either cvParam to arbitrate. Real bellini is 12.8% apart the same way.
        '<scanSettingsList count="1"><scanSettings id="scanSettings1">'
        '<cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixel x"'
        ' value="4" />',
        '<cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixel y"'
        ' value="4" />',
        '<cvParam cvRef="IMS" accession="IMS:1000044" name="max dimension x"'
        ' value="15" />',
        '<cvParam cvRef="IMS" accession="IMS:1000045" name="max dimension y"'
        ' value="15" />',
        '<cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size x"'
        ' value="4.406250e+00" />',
        '<cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y"'
        ' value="4.406250e+00" />',
        "</scanSettings>",
        "</scanSettingsList>",
        # componentList: present on all three real files, absent from every
        # ImzMLWriter fixture.
        '<instrumentConfigurationList count="1">'
        '<instrumentConfiguration id="IC1"><componentList count="3">'
        '<source order="1"><cvParam cvRef="MS" accession="MS:1000402"'
        ' name="secondary ionization" value="" />',
        "</source>",
        '<analyzer order="2"><cvParam cvRef="MS" accession="MS:1000084"'
        ' name="time-of-flight" value="" />',
        "</analyzer>",
        '<detector order="3"><cvParam cvRef="MS" accession="MS:1000114"'
        ' name="microchannel plate detector" value="" />',
        "</detector>",
        "</componentList>",
        "</instrumentConfiguration>",
        "</instrumentConfigurationList>",
        '<dataProcessingList count="1">'
        '<dataProcessing id="IONTOF SurfaceLab 7.5 Data Export">'
        '<processingMethod order="1" softwareRef="IONTOF SurfaceLab 7.5">'
        '<cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML"'
        ' value="" />',
        "</processingMethod>",
        "</dataProcessing>",
        "</dataProcessingList>",
        '<run defaultInstrumentConfigurationRef="IC1">'
        f'<spectrumList count="{len(IONTOF_COORDINATES)}"'
        ' defaultDataProcessingRef="IONTOF SurfaceLab 7.5 Data Export">',
    ]
    for i, (coordinate, placement) in enumerate(zip(IONTOF_COORDINATES, placements)):
        lines += _spectrum_xml(i, coordinate, placement, itemsize=8)
    lines += ["</spectrumList>", "</run>", "</mzML>"]

    return _write("iontof_sparse", lines, newline=b"\r\n", encoding="iso-8859-1")


# --------------------------------------------------------------------------
# The three LF fixtures
# --------------------------------------------------------------------------


def _canonical_header(
    uuid_text: str,
    scan_settings_lines: Sequence[str],
    mz_precision_lines: Sequence[str],
    n_spectra: int,
) -> List[str]:
    """The parts the three LF fixtures share.

    Deliberately conventional everywhere the fixture is not making a point:
    correct ontology names, ``mzML/@version``, a single 32-bit-float intensity
    declaration. Each fixture varies exactly one thing so a failing assertion
    points at that thing.
    """
    return [
        '<?xml version="1.0" encoding="utf-8"?>',
        '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">',
        '  <cvList count="3">',
        '    <cv id="MS" fullName="Proteomics Standards Initiative Mass'
        ' Spectrometry Ontology" version="4.1.35"'
        ' URI="http://psidev.info/ms/mzML/psi-ms.obo" />',
        '    <cv id="UO" fullName="Unit Ontology" version="releases/2020-03-10"'
        ' URI="http://ontologies.berkeleybop.org/uo.obo" />',
        '    <cv id="IMS" fullName="Imaging MS Ontology" version="0.9.1"'
        ' URI="http://www.maldi-msi.org/download/imzml/imagingMS.obo" />',
        "  </cvList>",
        "  <fileDescription>",
        "    <fileContent>",
        '      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"'
        ' value="" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000031" name="processed"'
        ' value="" />',
        '      <cvParam cvRef="MS" accession="MS:1000130" name="positive scan"'
        ' value="" />',
        f'      <cvParam cvRef="IMS" accession="IMS:1000080"'
        f' name="universally unique identifier" value="{{{uuid_text}}}" />',
        "    </fileContent>",
        "  </fileDescription>",
        '  <referenceableParamGroupList count="3">',
        '    <referenceableParamGroup id="mzArray">',
        '      <cvParam cvRef="MS" accession="MS:1000576" name="no compression"'
        ' value="" />',
        '      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array"'
        ' value="" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data"'
        ' value="true" />',
        *mz_precision_lines,
        "    </referenceableParamGroup>",
        '    <referenceableParamGroup id="intensityArray">',
        '      <cvParam cvRef="MS" accession="MS:1000576" name="no compression"'
        ' value="" />',
        '      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array"'
        ' value="" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data"'
        ' value="true" />',
        '      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"'
        ' value="" />',
        "    </referenceableParamGroup>",
        '    <referenceableParamGroup id="scan1">',
        '      <cvParam cvRef="MS" accession="MS:1000093" name="increasing m/z scan"'
        ' value="" />',
        "    </referenceableParamGroup>",
        '    <referenceableParamGroup id="spectrum1">',
        '      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"'
        ' value="" />',
        '      <cvParam cvRef="MS" accession="MS:1000130" name="positive scan"'
        ' value="" />',
        "    </referenceableParamGroup>",
        "  </referenceableParamGroupList>",
        '  <sampleList count="1">',
        '    <sample id="sample1" name="hand-authored fixture" />',
        "  </sampleList>",
        '  <softwareList count="1">',
        '    <software id="hand-authored" version="1" />',
        "  </softwareList>",
        *scan_settings_lines,
        # pyimzml dereferences instrumentConfigurationList unconditionally in
        # __readimzmlmeta and raises AttributeError on None, so even a fixture
        # with nothing to say about the instrument has to declare one.
        '  <instrumentConfigurationList count="1">',
        '    <instrumentConfiguration id="IC1">',
        '      <componentList count="2">',
        '        <source order="1"><cvParam cvRef="MS" accession="MS:1000075"'
        ' name="matrix-assisted laser desorption ionization" value="" /></source>',
        '        <detector order="2"><cvParam cvRef="MS" accession="MS:1000114"'
        ' name="microchannel plate detector" value="" /></detector>',
        "      </componentList>",
        "    </instrumentConfiguration>",
        "  </instrumentConfigurationList>",
        '  <dataProcessingList count="1">',
        '    <dataProcessing id="export">',
        '      <processingMethod order="1" softwareRef="hand-authored">',
        '        <cvParam cvRef="MS" accession="MS:1000544"'
        ' name="Conversion to mzML" value="" />',
        "      </processingMethod>",
        "    </dataProcessing>",
        "  </dataProcessingList>",
        '  <run defaultInstrumentConfigurationRef="IC1">',
        f'    <spectrumList count="{n_spectra}" defaultDataProcessingRef="export">',
    ]


DEFAULT_MZ_PRECISION = [
    '      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"'
    ' value="" />',
]

# A dense 2x2 acquisition: these three fixtures are not about sparsity.
DENSE_COORDINATES: List[Coordinate] = [(1, 1), (2, 1), (1, 2), (2, 2)]


def _build_lf_fixture(
    stem: str,
    uuid_text: str,
    scan_settings_lines: Sequence[str],
    mz_precision_lines: Sequence[str] = tuple(DEFAULT_MZ_PRECISION),
) -> Path:
    """Pack the .ibd and emit an LF/UTF-8 imzML around it."""
    spectra = _make_spectra(len(DENSE_COORDINATES), np.dtype(np.float64))
    ibd, placements = _pack_ibd(uuid_text, spectra)
    (FIXTURE_DIR / f"{stem}.ibd").write_bytes(ibd)

    lines = _canonical_header(
        uuid_text, scan_settings_lines, mz_precision_lines, len(DENSE_COORDINATES)
    )
    for i, (coordinate, placement) in enumerate(zip(DENSE_COORDINATES, placements)):
        lines += [
            "      " + line
            for line in _spectrum_xml(i, coordinate, placement, itemsize=8)
        ]
    lines += ["    </spectrumList>", "  </run>", "</mzML>"]

    return _write(stem, lines, newline=b"\n", encoding="utf-8")


NANOMETRE_UUID = "1D0C6F1B-9D7A-4E2B-9C33-0A7A5B2E4411"


def build_unit_nanometre() -> Path:
    """IMS:1000046/47 carrying UO:0000018 nanometre (audit #9).

    The physical pixel is 4.40625 um, declared as 4406.25 nm. pyimzml's
    ``convert_cv_param`` takes no unit argument, so ``imzmldict['pixel size x']``
    is the bare number 4406.25 -- which Thyra used to label micrometres, a
    1000x spatial error that ``convert_msi`` reported as success.

    The unit is not lost from the document, only from the dict pyimzml builds:
    it survives on ``metadata.scan_settings[...].cv_params``, which is the path
    ``ImzMLMetadataExtractor`` now reads to convert the value to micrometres.
    """
    scan_settings = [
        '  <scanSettingsList count="1">',
        '    <scanSettings id="scanSettings1">',
        '      <cvParam cvRef="IMS" accession="IMS:1000042"'
        ' name="max count of pixels x" value="2" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000043"'
        ' name="max count of pixels y" value="2" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)"'
        ' value="4406.25" unitCvRef="UO" unitAccession="UO:0000018"'
        ' unitName="nanometer" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y"'
        ' value="4406.25" unitCvRef="UO" unitAccession="UO:0000018"'
        ' unitName="nanometer" />',
        "    </scanSettings>",
        "  </scanSettingsList>",
    ]
    return _build_lf_fixture("unit_nanometre", NANOMETRE_UUID, scan_settings)


TWO_SCANSETTINGS_UUID = "2E1D7A2C-8B6E-4F1A-8D44-1B8B6C3F5522"


def build_two_scansettings() -> Path:
    """Two ``<scanSettings>`` blocks with different pixel sizes (audit #9's neighbour).

    ``__readimzmlmeta`` resolves each accession by first-match-anywhere below
    ``scanSettingsList``, so the dict it returns is a chimera: the pixel counts
    and the 10.0 um pixel size come from block 1, the max dimensions come from
    block 2 (the only block declaring them), and the resulting triple describes
    neither block. Block 2's own 25.0 um pixel size is never seen.
    """
    scan_settings = [
        '  <scanSettingsList count="2">',
        '    <scanSettings id="scanSettings1">',
        '      <cvParam cvRef="IMS" accession="IMS:1000042"'
        ' name="max count of pixels x" value="2" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000043"'
        ' name="max count of pixels y" value="2" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)"'
        ' value="10.0" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y"'
        ' value="10.0" />',
        "    </scanSettings>",
        '    <scanSettings id="scanSettings2">',
        '      <cvParam cvRef="IMS" accession="IMS:1000044" name="max dimension x"'
        ' value="200" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000045" name="max dimension y"'
        ' value="200" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)"'
        ' value="25.0" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y"'
        ' value="25.0" />',
        "    </scanSettings>",
        "  </scanSettingsList>",
    ]
    return _build_lf_fixture("two_scansettings", TWO_SCANSETTINGS_UUID, scan_settings)


TWO_PRECISION_UUID = "3F2E8B3D-7C5F-4A0B-9E55-2C9C7D4A6633"

TWO_PRECISION_MZ_LINES = [
    '      <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"'
    ' value="" />',
    '      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"'
    ' value="" />',
]


def build_two_precision_terms() -> Path:
    """One param group declaring both 32-bit and 64-bit float (audit #7).

    Schema-valid and ambiguous. ``__process_metadata`` keeps the *last* match
    while scanning ``PRECISION_DICT`` in insertion order, so the winner is
    "64-bit float" regardless of document order -- nothing in pyimzml notices
    that the group said two contradictory things.

    The .ibd here is genuinely float64, so today the guess happens to be right
    and the file converts correctly. That is deliberate: the fixture isolates
    the ambiguity itself, so a validator that refuses it is refusing the
    declaration and not reacting to corrupt data.
    """
    scan_settings = [
        '  <scanSettingsList count="1">',
        '    <scanSettings id="scanSettings1">',
        '      <cvParam cvRef="IMS" accession="IMS:1000042"'
        ' name="max count of pixels x" value="2" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000043"'
        ' name="max count of pixels y" value="2" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)"'
        ' value="10.0" />',
        '      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y"'
        ' value="10.0" />',
        "    </scanSettings>",
        "  </scanSettingsList>",
    ]
    return _build_lf_fixture(
        "two_precision_terms",
        TWO_PRECISION_UUID,
        scan_settings,
        mz_precision_lines=TWO_PRECISION_MZ_LINES,
    )


def main() -> None:
    """Rebuild all four pairs and report their sizes."""
    for builder in (
        build_iontof_sparse,
        build_unit_nanometre,
        build_two_scansettings,
        build_two_precision_terms,
    ):
        imzml_path = builder()
        ibd_path = imzml_path.with_suffix(".ibd")
        print(
            f"{imzml_path.name:>26}  {imzml_path.stat().st_size:>6} bytes"
            f"   {ibd_path.name:>24}  {ibd_path.stat().st_size:>6} bytes"
        )


if __name__ == "__main__":
    main()
