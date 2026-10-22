"""What each imzML spectrum states about its representation and instrument.

The spectrum type and the instrument follow the spectra that are converted,
not the file's summary (D78). So two things are recorded per spectrum, by
the hook in :mod:`thyra.readers.imzml._pyimzml_compat` during pyimzml's one
pass over the spectra:

- the representation it states, as one byte (see
  :data:`~thyra.metadata.constants.REPRESENTATION_UNSTATED` and its
  siblings), read on the spectrum and in the parameter groups it
  references;
- the ``instrumentConfigurationRef`` of its first ``<scan>``, as one byte
  coding into a short list of the references met.

The metadata extractor reads both for the spectra that are converted. A
parser built without the hook records neither, and the extractor then
falls back to what the file states as a whole.
"""

from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from ..metadata.constants import REPRESENTATION_CENTROID, REPRESENTATION_PROFILE

#: Parser attribute: one byte per spectrum, the representation it states.
REPRESENTATIONS_ATTRIBUTE = "thyra_spectrum_representations"

#: Parser attribute: one byte per spectrum, 0 where its scan names no
#: configuration, else one more than the reference's place in
#: :data:`CONFIGURATION_REFS_ATTRIBUTE`.
CONFIGURATION_CODES_ATTRIBUTE = "thyra_configuration_codes"

#: Parser attribute: the configuration references the scans name, in the
#: order they were first met.
CONFIGURATION_REFS_ATTRIBUTE = "thyra_configuration_refs"

#: References past this many are recorded as none, so the run default
#: decides for them. No file met names more than a handful.
_MAX_CODED_REFS = 254


def record_configuration(own: Dict[str, Any], ref: Optional[str]) -> None:
    """Append one spectrum's configuration reference to the parser's record.

    Args:
        own: The parser's ``__dict__``.
        ref: The reference its first scan names, or ``None``.
    """
    codes = own.get(CONFIGURATION_CODES_ATTRIBUTE)
    if codes is None:
        codes = own[CONFIGURATION_CODES_ATTRIBUTE] = bytearray()
        own[CONFIGURATION_REFS_ATTRIBUTE] = []
    if ref is None:
        codes.append(0)
        return
    refs: List[Optional[str]] = own[CONFIGURATION_REFS_ATTRIBUTE]
    try:
        codes.append(refs.index(ref) + 1)
    except ValueError:
        if len(refs) < _MAX_CODED_REFS:
            refs.append(ref)
            codes.append(len(refs))
        else:
            # A reference past the coded ones reads as unstated, so the
            # run default decides for it.
            codes.append(0)


def stated_representations(
    parser: Any, kept: Optional[Sequence[int]]
) -> Optional[Tuple[int, int]]:
    """How many of the given spectra state centroid, and how many profile.

    Args:
        parser: A pyimzml parser (or the head parser).
        kept: The spectra that are converted, by list position, or ``None``
            for every recorded spectrum.

    Returns:
        ``(centroid, profile)``; a spectrum stating both counts in both.
        ``None`` when the parser recorded nothing.
    """
    own = getattr(parser, "__dict__", {})
    recorded = own.get(REPRESENTATIONS_ATTRIBUTE)
    if not isinstance(recorded, bytearray):
        return None
    values = (
        recorded if kept is None else (recorded[i] for i in kept if i < len(recorded))
    )
    centroid = profile = 0
    for value in values:
        if value & REPRESENTATION_CENTROID:
            centroid += 1
        if value & REPRESENTATION_PROFILE:
            profile += 1
    return centroid, profile


def configuration_references(
    parser: Any, kept: Optional[Sequence[int]]
) -> Optional[Set[Optional[str]]]:
    """The configuration references the given spectra's scans name.

    Args:
        parser: A pyimzml parser (or the head parser).
        kept: The spectra that are converted, or ``None`` for every one.

    Returns:
        The set of references, ``None`` among them for a scan that names
        none. ``None`` when the parser recorded nothing.
    """
    own = getattr(parser, "__dict__", {})
    codes = own.get(CONFIGURATION_CODES_ATTRIBUTE)
    refs = own.get(CONFIGURATION_REFS_ATTRIBUTE)
    if not isinstance(codes, bytearray) or not isinstance(refs, list):
        return None
    used = set(codes) if kept is None else {codes[i] for i in kept if i < len(codes)}
    return {refs[code - 1] if code else None for code in used}
