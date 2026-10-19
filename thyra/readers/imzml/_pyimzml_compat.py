"""Make pyimzml survive a cvParam that has a numeric type but no value.

``pyimzml.ontology.ontology.convert_xml_value`` converts a cvParam's value
to the type its ontology entry declares. It catches ``ValueError`` -- an
empty ``value=""`` on a float-typed term becomes ``None`` -- but not
``TypeError``, so a term with **no** ``value`` attribute at all
(``float(None)``) aborts the whole metadata parse, before a single spectrum
is read.

That is exactly what the pyimzML fork behind TIMSCONVERT and TIMSImaging
writes for the ion mobility array it adds: a ``referenceableParamGroup``
with ``MS:1003006 mean inverse reduced ion mobility array`` carrying a unit
and no value, and pyimzml's table types ``MS:1003006`` as ``xsd:float``.
Upstream pyimzml therefore cannot open any such export. Measured on the
hand-authored ``mobility_continuous.imzML``: ``TypeError: float() argument
must be a string or a real number, not 'NoneType'`` from
``ImzMLParser.__init__``.

The patch is as narrow as the defect: a missing value on a typed term
converts to ``None``, the same result an empty value already produces.
Every other conversion is untouched. It is installed once, the first time
Thyra builds a parser, and is idempotent.
"""

import logging
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, cast

from ...utils.imzml_coordinate_base import (
    Z_FROM_GROUP,
    Z_GROUP_ATTRIBUTE,
    Z_STATED,
    Z_STATES_ATTRIBUTE,
    Z_UNSTATED,
)

logger = logging.getLogger(__name__)

_INSTALLED = False


def ensure_lenient_cv_param_values() -> None:
    """Patch pyimzml so a value-less typed cvParam parses as ``None``."""
    global _INSTALLED
    if _INSTALLED:
        return

    from pyimzml import metadata as pyimzml_metadata
    from pyimzml.ontology import ontology as pyimzml_ontology

    original = pyimzml_ontology.convert_xml_value

    def convert_xml_value(dtype, value):
        try:
            return original(dtype, value)
        except TypeError:
            if value is None:
                return None
            raise

    convert_xml_value.__doc__ = original.__doc__
    convert_xml_value.__wrapped__ = original  # type: ignore[attr-defined]

    pyimzml_ontology.convert_xml_value = convert_xml_value
    # metadata.py imports the function by name, so its own binding must be
    # replaced too; lookup_and_convert_cv_param resolves it through the
    # ontology module's globals, which the assignment above already covers.
    if getattr(pyimzml_metadata, "convert_xml_value", None) is original:
        pyimzml_metadata.convert_xml_value = convert_xml_value

    _INSTALLED = True
    logger.debug("Installed the lenient cvParam value conversion for pyimzml")


_SPECTRUM_TERMS_INSTALLED = False


def ensure_spectrum_term_fix() -> None:
    """Patch pyimzml to correct a mislabelled spectrum term as it does elsewhere.

    pyimzml corrects a cvParam an exporter wrote under the wrong accession
    ("positive scan" as ``MS:1000128``, the profile spectrum term) where it
    parses terms itself: ``fileContent`` and the parameter groups. The
    per-spectrum values asked for with ``include_spectra_metadata`` are
    looked up by accession alone, so a file stating that pair on each
    spectrum had no polarity at all, while the same pair in a group gave
    "positive" (D55). After the patch, ``spectrum_metadata_fields`` reads as
    if the file were written correctly. Nothing else the parser records
    changes, and a parser built without ``include_spectra_metadata`` is
    untouched.

    The hook is pyimzml's per-spectrum step, which is private.
    ``test_polarity_statements`` converts a file that states the pair on
    each spectrum, so a pyimzml release that renames the step fails there
    rather than losing the polarity without a word. Idempotent.
    """
    global _SPECTRUM_TERMS_INSTALLED
    if _SPECTRUM_TERMS_INSTALLED:
        return

    from pyimzml import ImzMLParser as parser_module

    from ...metadata.constants import MISLABELLED_POLARITY_TERMS

    parser_class = parser_module.ImzMLParser
    # Name-mangled, as pyimzml's own calls of it are.
    step = "_ImzMLParser__process_spectrum"
    original = getattr(parser_class, step)
    prefix = parser_module.XMLNS_PREFIX
    corrections = [
        (f'.//{prefix}cvParam[@accession="{accession}"][@name="{name}"]', fixed)
        for (accession, name), fixed in MISLABELLED_POLARITY_TERMS.items()
    ]

    def process_spectrum(self, elem, include_spectra_metadata):
        original(self, elem, include_spectra_metadata)
        fields = getattr(self, "spectrum_metadata_fields", None)
        if not isinstance(fields, dict):
            return
        for path, fixed in corrections:
            values = fields.get(fixed)
            # The step above appended this spectrum's value last; it is None
            # when the spectrum states no correct term.
            if values and values[-1] is None and elem.find(path) is not None:
                values[-1] = True

    process_spectrum.__wrapped__ = original  # type: ignore[attr-defined]
    setattr(parser_class, step, process_spectrum)
    _SPECTRUM_TERMS_INSTALLED = True
    logger.debug("Installed the spectrum term correction for pyimzml")


_SPECTRUM_CHECK_INSTALLED = False

_POSITION_ACCESSIONS = (
    ("x", "IMS:1000050"),
    ("y", "IMS:1000051"),
    ("z", "IMS:1000052"),
)
_OFFSET_ACCESSIONS = (
    ("external offset", "IMS:1000102"),
    ("external array length", "IMS:1000103"),
)


def ensure_spectrum_check() -> None:
    """Patch pyimzml to record each spectrum's z and to name what it cannot read.

    Two things happen in pyimzml's one pass over the spectra, so neither
    costs a second pass over the document:

    - Each spectrum's position z is recorded as one byte on the parser
      (``Z_STATES_ATTRIBUTE``): stated on the scan with a value, not
      stated, or stated only in a parameter group the scan references.
      pyimzml itself gives the last two the same stand-in z of 1. The
      reader judges the file from these bytes (D63).
    - A spectrum pyimzml fails on raises :class:`ConversionRefused` that
      names the spectrum and the cause, instead of the ``TypeError``,
      ``ValueError`` or ``AttributeError`` it raised from deep inside its
      constructor: a position held in a parameter group or a later
      ``<scan>``, a position that is not a whole number, a binary array
      read through a parameter group pyimzml did not pick.

    Nothing changes for a spectrum pyimzml reads. The hook is pyimzml's
    per-spectrum step, which is private, as for
    :func:`ensure_spectrum_term_fix`; the z tests convert files through it,
    so a pyimzml release that renames the step fails there. Idempotent.
    """
    global _SPECTRUM_CHECK_INSTALLED
    if _SPECTRUM_CHECK_INSTALLED:
        return

    from pyimzml import ImzMLParser as parser_module

    from ...errors import ConversionRefused

    parser_class = parser_module.ImzMLParser
    step = "_ImzMLParser__process_spectrum"
    original = getattr(parser_class, step)

    # Runs once per spectrum, so it is written for speed: precomputed tags,
    # plain-tag ``find`` (in C) and a loop over the scan's few children. A
    # path or an attribute predicate goes through ElementPath in Python,
    # which cost about 5 s on a 918,855-spectrum file.
    prefix = parser_module.XMLNS_PREFIX
    scan_list_tag = prefix + "scanList"
    scan_tag = prefix + "scan"
    cv_param_tag = prefix + "cvParam"

    def process_spectrum(self, elem, include_spectra_metadata):
        # Read the way pyimzml reads z: the first <scan>, its own cvParam.
        own = self.__dict__
        states = own.get(Z_STATES_ATTRIBUTE)
        if states is None:
            states = own[Z_STATES_ATTRIBUTE] = bytearray()
        state = Z_UNSTATED
        scan_list = elem.find(scan_list_tag)
        scan = scan_list.find(scan_tag) if scan_list is not None else None
        if scan is not None:
            for node in scan:
                if node.tag == cv_param_tag and node.get("accession") == "IMS:1000052":
                    if node.get("value") is not None:
                        state = Z_STATED
                    break
            else:
                group = _group_holding(self, scan, "IMS:1000052")
                if group is not None:
                    state = Z_FROM_GROUP
                    own.setdefault(Z_GROUP_ATTRIBUTE, (len(states), group))
        index = len(states)
        states.append(state)
        try:
            original(self, elem, include_spectra_metadata)
        except ConversionRefused:
            raise
        except (AttributeError, KeyError, TypeError, ValueError) as e:
            raise ConversionRefused(_why_unreadable(self, elem, index, e)) from e

    process_spectrum.__wrapped__ = original  # type: ignore[attr-defined]
    setattr(parser_class, step, process_spectrum)
    _SPECTRUM_CHECK_INSTALLED = True
    logger.debug("Installed the per-spectrum check for pyimzml")


def _group_positions(parser: Any) -> Dict[str, FrozenSet[str]]:
    """The position accessions each referenceable param group carries.

    Built once per parser, on the first spectrum that needs it.
    """
    cached = parser.__dict__.get("_thyra_group_positions")
    if cached is not None:
        return cast(Dict[str, FrozenSet[str]], cached)
    wanted = {accession for _axis, accession in _POSITION_ACCESSIONS}
    metadata = getattr(parser, "metadata", None)
    groups = getattr(metadata, "referenceable_param_groups", None)
    found: Dict[str, FrozenSet[str]] = {}
    if isinstance(groups, dict):
        for group_id, group in groups.items():
            held = frozenset(
                param[1]
                for param in getattr(group, "cv_params", [])
                if len(param) > 1 and param[1] in wanted
            )
            if held:
                found[str(group_id)] = held
    parser.__dict__["_thyra_group_positions"] = found
    return found


def _group_holding(parser: Any, scan: Any, accession: str) -> Optional[str]:
    """The id of a group ``scan`` references that carries ``accession``."""
    positions = _group_positions(parser)
    if not positions:
        # No group in the file holds a position: nothing to look up, and
        # the case of nearly every file, so it costs one dict check.
        return None
    refs = scan.findall(f"{parser.sl}referenceableParamGroupRef")
    for ref in refs:
        group_id = ref.get("ref")
        if group_id is not None and accession in positions.get(group_id, ()):
            return str(group_id)
    return None


def _why_unreadable(parser: Any, elem: Any, index: int, error: Exception) -> str:
    """One sentence on why pyimzml could not read spectrum ``index``."""
    source = parser.filename
    name = Path(str(getattr(source, "name", source))).name
    reason = _array_problem(parser, elem) or _position_problem(parser, elem)
    if reason is None:
        reason = f"the imzML parser failed on it with {type(error).__name__}: {error}"
    return f"imzML {name}, spectrum {index}: {reason}"


def _array_problem(parser: Any, elem: Any) -> Optional[str]:
    """What is wrong with the spectrum's binary arrays, or ``None``."""
    ns = parser.sl
    arrays = elem.find(f"{ns}binaryDataArrayList")
    if arrays is None:
        return "it has no <binaryDataArrayList>, so it holds no m/z or intensity array."
    refs: List[str] = []
    for array in arrays:
        ref = array.find(f"{ns}referenceableParamGroupRef")
        if ref is None or ref.get("ref") is None:
            return (
                "one of its binary arrays references no parameter group, so "
                "the parser cannot tell what the array holds."
            )
        refs.append(str(ref.get("ref")))

    if parser.mzGroupId not in refs or parser.intGroupId not in refs:
        return _group_mismatch(parser, refs)

    for label, group_id in (
        ("m/z", parser.mzGroupId),
        ("intensity", parser.intGroupId),
    ):
        array = arrays[refs.index(group_id)]
        for term, accession in _OFFSET_ACCESSIONS:
            node = array.find(f'{ns}cvParam[@accession="{accession}"]')
            value = None if node is None else node.get("value")
            if value is None:
                return f"its {label} array states no {term} ({accession})."
            try:
                int(value)
            except ValueError:
                return (
                    f"its {label} array states its {term} ({accession}) as "
                    f"{value!r}, which is not a whole number."
                )
    return None


def _group_mismatch(parser: Any, refs: List[str]) -> str:
    """Why a spectrum's arrays miss the groups pyimzml reads them through."""
    groups = getattr(parser.metadata, "referenceable_param_groups", None) or {}
    sentence = (
        f"its binary arrays reference the parameter groups "
        f"{', '.join(map(repr, refs))}, but the parser reads the m/z array "
        f"only through group {parser.mzGroupId!r} and the intensity array "
        f"only through group {parser.intGroupId!r}."
    )
    for label, term in (("m/z", "m/z array"), ("intensity", "intensity array")):
        declared = [str(g) for g, group in groups.items() if term in group]
        if len(declared) > 1:
            sentence += (
                f" The file declares the {label} array in {len(declared)} "
                f"groups ({', '.join(map(repr, declared))}), and the parser "
                f"takes one. Declare each array in one group."
            )
    return sentence


def _position_problem(parser: Any, elem: Any) -> Optional[str]:
    """What is wrong with the spectrum's position, or ``None``."""
    ns = parser.sl
    scans = elem.findall(f"{ns}scanList/{ns}scan")
    if not scans:
        return "it has no <scan>, so it states no position (IMS:1000050, IMS:1000051)."
    scan = scans[0]
    for axis, accession in _POSITION_ACCESSIONS:
        node = scan.find(f'{ns}cvParam[@accession="{accession}"]')
        value = None if node is None else node.get("value")
        if value is None:
            if axis == "z":
                # pyimzml stands in 1 for a missing z; the reader judges
                # that from the recorded states, not here.
                continue
            return _missing_position(parser, scans, axis, accession)
        try:
            int(value)
        except ValueError:
            return (
                f"it states position {axis} ({accession}) as {value!r}. "
                f"Positions are whole numbers (xsd:nonNegativeInteger)."
            )
    return None


def _missing_position(parser: Any, scans: List[Any], axis: str, accession: str) -> str:
    """Where a position the first scan lacks is, if it is anywhere."""
    ns = parser.sl
    if scans[0].find(f'{ns}cvParam[@accession="{accession}"]') is not None:
        return f"its position {axis} ({accession}) has no value."
    group = _group_holding(parser, scans[0], accession)
    if group is not None:
        return (
            f"it takes its position {axis} ({accession}) from parameter group "
            f"{group!r}. The imzML parser reads a position only from the scan "
            f"itself; state it there."
        )
    for number, other in enumerate(scans[1:], start=2):
        if other.find(f'{ns}cvParam[@accession="{accession}"]') is not None:
            return (
                f"it states its position {axis} ({accession}) in <scan> "
                f"{number} of its scan list, and the imzML parser reads only "
                f"the first <scan>."
            )
    return f"it states no position {axis} ({accession})."
