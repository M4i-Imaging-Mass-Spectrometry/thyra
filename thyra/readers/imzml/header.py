"""The head of an imzML document, read without the spectra behind it.

An imzML file states its raster, its pixel pitch, its instrument and how
many spectra it holds in the block before ``<run>``.  Everything a preview
card shows, except the mass range, is already there -- so a preview does
not have to parse the spectrum list, and must not decode the ``.ibd``.

It used to do both.  ``ImzMLReader`` swallowed ``metadata_only`` through
``**kwargs``, so ``preview_msi`` built an ordinary parser: pyimzml walked
every ``<spectrum>`` element to collect the offsets, and then
``_get_mass_range_processed`` decoded each spectrum's m/z array out of the
``.ibd`` to find the extrema.  Measured on a 2.0 GiB imzML with 918,855
spectra: **64 s**, against the ``<500 ms`` ``preview.py`` promises and the
"no spectra are read" it promises alongside it.  Reading the head of the
same file takes **0.4 ms**, and returns the same numbers -- 918,855
spectra, a 1007 x 1469 raster, a 5 um pitch (issue #360).  This is the
same defect PHI had in #240, in a different format.

Two things the head cannot answer, and what is done about each:

- **The mass range.**  Writers that follow the spec record
  ``MS:1000528``/``MS:1000527`` -- the lowest and highest observed m/z --
  on every ``<spectrum>``, and those are exactly the extrema the ``.ibd``
  scan recomputes: verified identical on two files, to the six decimals
  the XML stores.  :func:`scan_observed_mz_range` reads them off the XML,
  which costs the document's size and not the acquisition's (53 ms for a
  29 MB file, 3.9 s for a 2.0 GiB one) and never opens the binary.  A
  writer that records neither -- IONTOF SurfaceLab is one -- leaves the
  range genuinely unknown, and it is reported as unknown rather than
  guessed from one spectrum, which on a real file was 28 Da narrow.

- **The raster, when the declaration is missing or ambiguous.**
  ``IMS:1000042`` / ``IMS:1000043`` are what the preview reports as the
  grid.  A file that omits them cannot be previewed from its head; nor
  can a 0-based one, because the number means a count to the spec and a
  largest-coordinate to some writers, and those differ by one exactly
  there.  Both fall back to the full parse rather than invent a number;
  see :func:`~thyra.metadata.extractors.imzml_header_extractor.declared_raster`
  and :meth:`ImzMLReader._create_metadata_extractor`.

Note the declared raster is the acquisition's *intended* extent, while a
conversion sizes its store from the coordinates actually present.  The two
agree on every local file, and disagree exactly when an acquisition was
interrupted -- the preview then describes what was set up, the store what
arrived.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Dict, NamedTuple, Optional, Tuple, Union
from xml.etree import ElementTree as ET  # nosec B405

from ._pyimzml_compat import ensure_lenient_cv_param_values

logger = logging.getLogger(__name__)

MZML_NS = "{http://psi.hupo.org/ms/mzml}"

# The accessions pyimzml's own ``__readimzmlmeta`` collects, under the
# names it files them under, so ``imzmldict`` means the same thing on this
# path as on the full one.  Keyed by accession rather than by the CV name
# the file spells out: real writers disagree on the name ("max count of
# pixels x" against IONTOF's "max count of pixel x") and agree on the
# number.
_SCAN_SETTINGS_PARAMS = (
    ("max count of pixels x", "IMS:1000042"),
    ("max count of pixels y", "IMS:1000043"),
    ("max dimension x", "IMS:1000044"),
    ("max dimension y", "IMS:1000045"),
    ("pixel size x", "IMS:1000046"),
    ("pixel size y", "IMS:1000047"),
    ("matrix solution concentration", "MS:1000835"),
)

_INSTRUMENT_CONFIG_PARAMS = (
    ("wavelength", "MS:1000843"),
    ("focus diameter x", "MS:1000844"),
    ("focus diameter y", "MS:1000845"),
    ("pulse energy", "MS:1000846"),
    ("pulse duration", "MS:1000847"),
    ("attenuation", "MS:1000848"),
)

LOWEST_OBSERVED_MZ = "MS:1000528"
HIGHEST_OBSERVED_MZ = "MS:1000527"

POSITION_X = "IMS:1000050"
POSITION_Y = "IMS:1000051"

EXTERNAL_OFFSET = "IMS:1000102"
EXTERNAL_ENCODED_LENGTH = "IMS:1000104"

_LOWEST = re.compile(
    rb'accession="' + LOWEST_OBSERVED_MZ.encode() + rb'"[^>]*?value="([^"]*)"'
)
_HIGHEST = re.compile(
    rb'accession="' + HIGHEST_OBSERVED_MZ.encode() + rb'"[^>]*?value="([^"]*)"'
)

# Big enough that the read is sequential, small enough that the file is
# never held in memory.
_SCAN_CHUNK_BYTES = 8 << 20
# Any single cvParam element is a few hundred bytes, so this overlap cannot
# cut one in half.  Elements inside it are matched twice, which costs
# nothing: the scan keeps a minimum and a maximum, and both are idempotent.
_SCAN_OVERLAP_BYTES = 4096


class ImzMLHeaderParser:
    """What the extractor needs from a parser, read from the head alone.

    Stands in for :class:`ImzMLParser` on the metadata-only path and
    carries the same three attributes the metadata extractor reads --
    ``metadata``, ``imzmldict`` and (here) the spectrum count.

    It deliberately has **no** ``coordinates`` and no ``getspectrum``.
    Both are the expensive half, both are what this class exists to
    avoid, and a stub returning an empty list would turn "we did not
    read the spectra" into "this file has no spectra" somewhere far from
    here.  Anything that needs them gets an ``AttributeError`` naming the
    attribute, which is the whole diagnosis.
    """

    def __init__(
        self,
        root: ET.Element,
        n_spectra: Optional[int],
        first_position: Optional[Tuple[int, int]] = None,
    ) -> None:
        """Wrap a parsed document head.

        Args:
            root: The ``mzML`` element, built as far as the first
                spectrum by :func:`_parse_head`.
            n_spectra: The count ``<spectrumList>`` declares, or ``None``
                when it declares none.
            first_position: The first spectrum's ``(x, y)``, or ``None``
                when it states none.  Read for what it says about the
                coordinate base, not for the pixel.
        """
        from pyimzml.metadata import Metadata

        self.metadata = Metadata(root)
        self.imzmldict: Dict[str, Any] = _imzmldict(self.metadata)
        self.n_spectra = n_spectra
        self.first_position = first_position


def read_header(path: Union[str, Path]) -> ImzMLHeaderParser:
    """Parse an imzML document up to its spectrum list.

    Args:
        path: The ``.imzML`` file.

    Returns:
        An :class:`ImzMLHeaderParser` over the document's head.

    Raises:
        Exception: Whatever ``ElementTree`` raises on a document it
            cannot parse.  A malformed head is a malformed file, and the
            caller (``preview_msi``) turns it into ``readable=False``.
    """
    ensure_lenient_cv_param_values()
    root, n_spectra, first_position = _parse_head(Path(path))
    return ImzMLHeaderParser(root, n_spectra, first_position)


def _parse_head(
    path: Path,
) -> Tuple[ET.Element, Optional[int], Optional[Tuple[int, int]]]:
    """The document tree, built only as far as the first spectrum.

    ``iterparse`` attaches each element to its parent as it opens it, so
    by the time the spectrum list starts, every block :class:`Metadata`
    reads -- file description, referenceable param groups, samples,
    software, scan settings, instrument configurations, data processing --
    is complete in the tree.  Stopping there leaves ``<run>`` present and
    nearly empty, which ``Metadata`` never looks at.

    Three things are taken on the way:

    - the tree, for ``Metadata``;
    - ``<spectrumList count=...>``, the number of spectra the file says it
      holds, which is what the coordinate list would have been counted
      for;
    - the first spectrum's ``IMS:1000050``/``IMS:1000051`` position, which
      is what says whether the declared raster can be trusted.  See
      :func:`~thyra.metadata.extractors.imzml_header_extractor.declared_raster`.
    """
    n_spectra: Optional[int] = None
    position: Dict[str, int] = {}
    with path.open("rb") as handle:
        context = ET.iterparse(handle, events=("start",))  # nosec B314
        root: Optional[ET.Element] = None
        spectra_seen = 0
        for _, element in context:
            if root is None:
                root = element
            tag = element.tag
            if tag.endswith("}spectrumList"):
                count = element.get("count")
                if count is not None:
                    try:
                        n_spectra = int(count)
                    except ValueError:
                        logger.debug("Unparseable spectrumList count %r", count)
            elif tag.endswith("}spectrum"):
                spectra_seen += 1
                # The head is read for what precedes the spectra plus the
                # first of them; the second is where it stops.
                if spectra_seen > 1:
                    break
            elif tag.endswith("}cvParam") and spectra_seen == 1:
                accession = element.get("accession", "")
                if accession in (POSITION_X, POSITION_Y):
                    try:
                        position[accession] = int(element.get("value", ""))
                    except ValueError:
                        logger.debug("Unparseable %s on the first spectrum", accession)
                    if len(position) == 2:
                        break
        if root is None:
            raise ValueError(f"{path} contains no XML elements")

    first_position = (
        (position[POSITION_X], position[POSITION_Y]) if len(position) == 2 else None
    )
    return root, n_spectra, first_position


def _imzmldict(metadata: Any) -> Dict[str, Any]:
    """Rebuild pyimzml's ``imzmldict`` from the parsed head.

    pyimzml fills the same dict by searching the tree for each accession
    in document order, first match wins; the parsed ``ParamGroup``s hold
    those values already converted, so they are read from there and the
    same first-match rule applied over the groups.
    """
    found: Dict[str, Any] = {}
    _collect(found, getattr(metadata, "scan_settings", None), _SCAN_SETTINGS_PARAMS)
    _collect(
        found,
        getattr(metadata, "instrument_configurations", None),
        _INSTRUMENT_CONFIG_PARAMS,
    )
    return found


def _collect(
    into: Dict[str, Any],
    groups: Any,
    wanted: Tuple[Tuple[str, str], ...],
) -> None:
    """Fill ``into`` from whichever group first declares each accession.

    ``groups`` is typed ``Any`` because it comes off pyimzml's ``Metadata``
    by name: the shape is checked here rather than promised by a
    signature, the way the rest of the imzML metadata code treats that
    object.
    """
    if not isinstance(groups, dict):
        return
    for group in groups.values():
        by_accession = getattr(group, "param_by_accession", None)
        if not isinstance(by_accession, dict):
            continue
        for name, accession in wanted:
            if name in into:
                continue
            value = by_accession.get(accession)
            if value is not None:
                into[name] = value


def scan_observed_mz_range(
    path: Union[str, Path],
) -> Optional[Tuple[float, float]]:
    """The mass range the file itself records, or ``None``.

    Reads every ``MS:1000528`` (lowest observed m/z) and ``MS:1000527``
    (highest observed m/z) in the document and returns their extremes.
    These are the writer's own record of what each spectrum contains, and
    on the files checked they reproduce a full ``.ibd`` scan exactly --
    ``(250.0125002799743, 1200.0)`` from the binary against
    ``(250.0125, 1200.0)`` from the XML, which is the same number to the
    six decimals the XML carries.

    ``None`` means the file records neither term, not that it has no
    spectra: IONTOF SurfaceLab writes no observed-m/z cvParams at all.
    Callers report that as unknown.  Taking the first spectrum's own
    extrema instead was measured 28 Da narrow on a real file, which is a
    range a reader would believe.

    The scan is a byte scan rather than an XML parse: it is the m/z
    terms that are wanted, not the document structure, and the cost is
    then the file's size rather than its element count -- 53 ms for
    29 MB, 3.9 s for 2.0 GiB.
    """
    path = Path(path)
    lowest: Optional[float] = None
    highest: Optional[float] = None

    with path.open("rb") as handle:
        carry = b""
        while True:
            chunk = handle.read(_SCAN_CHUNK_BYTES)
            if not chunk:
                break
            buffer = carry + chunk
            for match in _LOWEST.finditer(buffer):
                value = _as_float(match.group(1))
                if value is not None and (lowest is None or value < lowest):
                    lowest = value
            for match in _HIGHEST.finditer(buffer):
                value = _as_float(match.group(1))
                if value is not None and (highest is None or value > highest):
                    highest = value
            carry = buffer[-_SCAN_OVERLAP_BYTES:]

    if lowest is None or highest is None:
        logger.info(
            "%s records no observed m/z range (MS:1000528/MS:1000527); "
            "a metadata-only read reports it as unknown.",
            path.name,
        )
        return None
    return (lowest, highest)


_ARRAY_BLOCK = re.compile(rb"<binaryDataArray\b.*?</binaryDataArray>", re.DOTALL)
_EXTERNAL_OFFSET = re.compile(
    rb'accession="' + EXTERNAL_OFFSET.encode() + rb'"[^>]*?value="([0-9]+)"'
)
_EXTERNAL_ENCODED_LENGTH = re.compile(
    rb'accession="' + EXTERNAL_ENCODED_LENGTH.encode() + rb'"[^>]*?value="([0-9]+)"'
)
_SPECTRUM_INDEX = re.compile(rb"<spectrum\b[^>]*?index=\"([0-9]+)\"")

# A spectrum's record is on the order of a kilobyte; this holds the last
# several of them whatever the acquisition's size.
_TAIL_BYTES = 64 << 10


class DeclaredExtent(NamedTuple):
    """Where the document says its last binary array ends."""

    end_byte: int
    spectrum_index: Optional[int]


def declared_binary_extent(path: Union[str, Path]) -> Optional[DeclaredExtent]:
    """The end byte of the last array the document declares, from its tail.

    A preview that reads only the head cannot check the ``.ibd`` the way
    :meth:`ImzMLReader._validate_ibd_extent` does -- that check walks every
    spectrum's offsets, which is the parse being avoided.  The last one is
    enough to catch the case that matters, a truncated or partially copied
    binary, and costs a single 64 KB read at the end of the file.

    Returns ``None`` when the tail holds no complete ``<binaryDataArray>``,
    which is not a verdict -- the caller simply has nothing to check.

    One-sided by construction, and deliberately so: the last array in
    document order is the furthest into the binary for every writer met,
    but that is not guaranteed (the full validator warns when offsets are
    not monotonic), so this can miss a bad file.  It cannot report a good
    one, because it only fires when the XML claims bytes the ``.ibd`` does
    not have.
    """
    path = Path(path)
    size = path.stat().st_size
    with path.open("rb") as handle:
        handle.seek(max(0, size - _TAIL_BYTES))
        tail = handle.read()

    blocks = _ARRAY_BLOCK.findall(tail)
    if not blocks:
        return None
    last = blocks[-1]

    offset = _EXTERNAL_OFFSET.search(last)
    length = _EXTERNAL_ENCODED_LENGTH.search(last)
    if offset is None or length is None:
        return None

    indices = _SPECTRUM_INDEX.findall(tail)
    return DeclaredExtent(
        end_byte=int(offset.group(1)) + int(length.group(1)),
        spectrum_index=int(indices[-1]) if indices else None,
    )


def _as_float(raw: bytes) -> Optional[float]:
    """A cvParam value as a finite float, or ``None``.

    A blank or malformed value is skipped rather than raised on: one bad
    spectrum out of a million should not cost the whole range, and the
    extremes are taken over whatever else parsed.
    """
    try:
        value = float(raw)
    except ValueError:
        return None
    if value != value or value in (float("inf"), float("-inf")):
        return None
    return value
