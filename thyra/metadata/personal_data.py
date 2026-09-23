# thyra/metadata/personal_data.py
"""What a store does not copy out of a source's own metadata.

Vendor software records the person at the instrument -- who calibrated it,
who ran it -- and paths on the acquisition PC, whose folders carry a user's
name. The extractors used to hand both on as read, so they reached every
store (a table's ``uns``, and the root attributes as well) and everything
else that reads the vendor dictionaries.

Neither does anything in a store. Nothing reads them, and the source they
came from still holds every one -- ``provenance.source_path`` says where.
A store, meanwhile, is the thing that gets handed on, and whoever hands it
on cannot be expected to know what the vendor dictionaries hold. So the
choice is made once, where every format's metadata passes
(:meth:`MetadataExtractor.get_comprehensive`), rather than left to be
remembered:

* **A field that names a person is not copied.** A name has no shape, so
  it is matched by what the field is: the vendor spellings in
  :data:`PERSON_KEYS`, whether they are a key or the name of a user
  parameter, and, for sources that state their metadata as
  controlled-vocabulary parameters, the PSI-MS contact terms in
  :data:`PERSON_ACCESSIONS`. A value that *is* an e-mail address is
  dropped whatever its key, because an address does have a shape.
* **A path keeps its last component.** The file name the acquisition was
  saved under is the part a reader can use; the folders above it describe
  a machine the reader does not have. Matched by *shape* (see
  :func:`is_path`) rather than by key, as the repository's own
  ``no-lab-share-paths`` and ``no-user-home-paths`` hooks are, so a key a
  vendor adds later is caught too.

The paths Thyra itself records about where it read the source are the
exception, named in :data:`LOCATING_KEYS`: they are on the machine that
holds the store and they find the source again, which is why #384 kept
``provenance.source_path`` whole in a store.

Free text the lab chose -- a sample name, a method name, a project name --
is not searched for names. It identifies the data, and a lab that writes a
person's name into it has made the name part of the data's identity;
nothing here could tell it from any other word.

A dictionary that an extractor turns into a JSON string before handing it
on is opaque here, so that extractor strips it first: PHI does, for its
header and for the blocks appended to its file.
"""

import re
from collections.abc import Mapping
from dataclasses import replace
from functools import lru_cache
from typing import Any, Dict, FrozenSet, Optional

from .types import ComprehensiveMetadata

__all__ = [
    "LOCATING_KEYS",
    "PERSON_ACCESSIONS",
    "PERSON_KEYS",
    "file_name",
    "is_path",
    "strip_personal_data",
    "without_personal_data",
]

#: Vendor keys whose value is a person, compared case-blind and ignoring
#: spaces, underscores and hyphens. Each is the spelling of a source Thyra
#: reads: ``OperatorName`` (Bruker ``GlobalMetadata`` and the solariX
#: ``Properties`` table), ``CalibrationUser`` and ``MobilityCalibrationUser``
#: (timsTOF ``calibration.sqlite``), ``Operator`` and ``User Name`` (the PHI
#: header). Thyra's own spellings of the same facts, ``operator_name`` and
#: ``calibration_user``, fold onto these.
PERSON_KEYS: FrozenSet[str] = frozenset(
    {
        "operatorname",
        "operator",
        "username",
        "calibrationuser",
        "mobilitycalibrationuser",
    }
)

#: The PSI-MS contact attributes (MS:1000585) whose value reaches the person
#: rather than the institution: everything but ``contact affiliation``
#: (MS:1000590), which names an organisation and is kept.
PERSON_ACCESSIONS: FrozenSet[str] = frozenset(
    {
        "MS:1000586",  # contact name
        "MS:1000587",  # contact address
        "MS:1000588",  # contact URL
        "MS:1000589",  # contact email
        "MS:1001755",  # contact phone number
        "MS:1001756",  # contact fax number
        "MS:1001757",  # contact toll-free phone number
    }
)

#: The keys under which Thyra records where it read the source, at the top
#: level of a vendor dictionary: ``data_path`` (Bruker tsf/tdf and Waters),
#: ``database_path`` and ``binary_file`` (Bruker tsf/tdf), ``ibd_file``
#: (imzML) and ``mis_file`` (the ``.mis`` a solariX reader found). Matched
#: exactly: a vendor's own key is never spelled this way.
LOCATING_KEYS: FrozenSet[str] = frozenset(
    {"data_path", "database_path", "binary_file", "ibd_file", "mis_file"}
)

_PATH_SHAPE = re.compile(
    r"""
    ^[A-Za-z]:[\\/]     # a drive letter, with either separator
    | ^\\\\[^\\]        # a UNC share
    | ^file:/           # a file URI
    | ^/[^/\s]+/        # a POSIX absolute path of two components or more
    | ^[^\\/:*?"<>|\r\n]+(?:\\[^\\/:*?"<>|\r\n]+)*
      \\[^\\/:*?"<>|\r\n]+\.[A-Za-z0-9]{1,8}$
                        # a relative path written on Windows, ending in a
                        # file with an extension, with nothing in it that
                        # Windows forbids in a name
    """,
    re.VERBOSE | re.IGNORECASE,
)

_EMAIL_SHAPE = re.compile(
    r"^[^@\s:/<>]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)*\.[A-Za-z]{2,}$"
)

# A value of one of these types is kept as it is, and is the common case:
# a per-frame table is thousands of rows of numbers.
_PLAIN = (bool, int, float, bytes, type(None))

_DROP = object()


def is_path(value: str) -> bool:
    """Whether ``value`` has the shape of a path on some machine.

    Deliberately blind to a forward-slash string with no leading slash,
    drive letter or scheme: ``m/z``, ``1/K0``, ``Vs/cm2`` and the
    month-first date PHI writes (``06/23/2026 21:12:35``) all look like
    that, and none of them is a path. A relative path written with
    backslashes counts only when it ends in a file name with an
    extension, so free text with a backslash in it (``brain slice 3``
    written with one, or a pattern) is left alone.

    A backslash alone is never enough. A vendor dictionary can hold a
    JSON document -- PHI keeps its whole header as one -- and JSON writes
    a micro sign or a quote as an escape that starts with a backslash.
    Every JSON document has quotes in it, which Windows forbids in a name,
    so none has the shape of a relative path; taken for one, it would be
    cut down to whatever followed its last backslash. Quotes around a
    value as a whole are ignored, so a quoted path is still a path.
    """
    return bool(_PATH_SHAPE.match(_unquoted(value)))


def file_name(value: str) -> Optional[str]:
    """The last component of ``value`` on either path separator, or ``None``.

    Vendor software records a file as it was opened, which on a Bruker,
    Waters or PHI acquisition PC is an absolute path; the name is the part
    worth keeping, and the machine it came from is not.
    """
    name = re.split(r"[\\/]", value)[-1].strip()
    return name or None


def strip_personal_data(mapping: Mapping[str, Any]) -> Dict[str, Any]:
    """A copy of one vendor dictionary without the people in it.

    Recurses through nested mappings and lists. A person key, a person
    parameter and an e-mail address are dropped; a path becomes its file
    name. A key in :data:`LOCATING_KEYS` at the top level is kept as it
    is. ``mapping`` is not modified, and a nested mapping or list with
    nothing to take out is shared with it rather than copied.
    """
    result: Dict[str, Any] = {}
    for key, value in mapping.items():
        if key in LOCATING_KEYS:
            result[key] = value
        elif not _is_person_key(key):
            stripped = _strip(value)
            if stripped is not _DROP:
                result[key] = stripped
    return result


def without_personal_data(
    comprehensive: ComprehensiveMetadata,
) -> ComprehensiveMetadata:
    """Apply :func:`strip_personal_data` to each vendor dictionary.

    The essential metadata is left as it is: its ``source_path`` is the
    path Thyra read the source from, and it stays whole in a store.
    """
    return replace(
        comprehensive,
        format_specific=_strip_section(comprehensive.format_specific),
        acquisition_params=_strip_section(comprehensive.acquisition_params),
        instrument_info=_strip_section(comprehensive.instrument_info),
        raw_metadata=_strip_section(comprehensive.raw_metadata),
    )


def _strip_section(section: Any) -> Any:
    if isinstance(section, Mapping):
        return strip_personal_data(section)
    return section


def _unquoted(value: str) -> str:
    return value.strip().strip("\"'")


def _is_email(value: str) -> bool:
    return bool(_EMAIL_SHAPE.match(value.strip()))


@lru_cache(maxsize=4096)
def _folded(key: str) -> str:
    return re.sub(r"[\s_-]", "", key).lower()


def _is_person_key(key: Any) -> bool:
    return isinstance(key, str) and _folded(key) in PERSON_KEYS


def _names_a_person(param: Mapping) -> bool:
    """A contact CV parameter, or a user parameter named like a person field."""
    accession = param.get("accession")
    if isinstance(accession, str):
        return accession in PERSON_ACCESSIONS
    return "value" in param and _is_person_key(param.get("name"))


def _strip_text(value: str) -> Any:
    if _is_email(value):
        return _DROP
    if is_path(value):
        return file_name(_unquoted(value).rstrip("\\/")) or _DROP
    return value


def _strip(value: Any) -> Any:
    """``value`` without the people in it; ``value`` itself when it had none.

    Handing an untouched container back rather than a copy is what keeps
    a per-frame table -- one row per pixel, a million rows on a large
    slide -- from being duplicated in memory for nothing. Nothing edits a
    vendor dictionary after extraction, so sharing the rows is safe.
    """
    if isinstance(value, str):
        return _strip_text(value)
    if isinstance(value, _PLAIN):
        return value
    if isinstance(value, Mapping):
        if _names_a_person(value):
            return _DROP
        kept: Dict[Any, Any] = {}
        changed = False
        for key, item in value.items():
            if _is_person_key(key):
                changed = True
            elif isinstance(item, _PLAIN):
                kept[key] = item
            else:
                stripped = _strip(item)
                if stripped is not _DROP:
                    kept[key] = stripped
                changed = changed or stripped is not item
        return value if not changed and type(value) is dict else kept
    if isinstance(value, (list, tuple)):
        items = []
        changed = False
        for item in value:
            stripped = item if isinstance(item, _PLAIN) else _strip(item)
            if stripped is not _DROP:
                items.append(stripped)
            changed = changed or stripped is not item
        if not changed:
            return value
        return tuple(items) if isinstance(value, tuple) else items
    return value
