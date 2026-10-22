"""The controlled vocabulary the metadata layer speaks.

Accession codes and the CV names behind them: what an imzML file says
about itself, in the spelling ``uns["msi_metadata"]`` stores. These
described data before anything resampled it, and they describe data no
resampler ever touches -- ``ImzMLAccessions.UUID`` is the binary file's
identity, ``BinaryDataType`` is how the ``.ibd`` is laid out.

**They lived in ``thyra/resampling/constants.py`` until issue #381.** That
gave three vendor metadata extractors and the imzML reader an import of
the resampling package, which is backwards -- resampling reads a
description of the data, it does not define one -- and it is how importing
a metadata extractor came to load all 25 resampling modules. The old
module re-exports every name here, so both spellings keep working;
``Thresholds``, the one genuine resampling number in it, stayed there.
"""

from typing import Any, Dict, Iterable, NamedTuple, Optional, Tuple

from ..errors import ConversionRefused


class ImzMLAccessions:
    """PSI-MS and imzML controlled vocabulary accession codes."""

    # Spectrum type (MS ontology)
    CENTROID_SPECTRUM = "MS:1000127"
    PROFILE_SPECTRUM = "MS:1000128"

    # Scan polarity (MS ontology)
    POSITIVE_SCAN = "MS:1000130"
    NEGATIVE_SCAN = "MS:1000129"

    # Binary data type (imzML ontology)
    CONTINUOUS_BINARY = "IMS:1000030"
    PROCESSED_BINARY = "IMS:1000031"

    # Pixel size (imzML ontology)
    PIXEL_SIZE_X = "IMS:1000046"
    PIXEL_SIZE_Y = "IMS:1000047"

    # Software identifiers
    SCILS_LAB = "MS:1002384"

    # Binary-file UUID (imzML ontology). pyimzml keys ``param_by_name`` by
    # resolved CV *name*, not by accession, so the name is what a lookup
    # needs; the accession is kept beside it for error messages and docs.
    UUID = "IMS:1000080"
    UUID_NAME = "universally unique identifier"


class SpectrumType:
    """Spectrum type string constants."""

    CENTROID = "centroid spectrum"
    PROFILE = "profile spectrum"


# What a caller may write for a spectrum representation, mapped to the CV name
# Thyra stores. The bare words match SCiLS Lab's ``--rep_type PROFILE|CENTROID``
# (2026b User Guide p.81); the full CV names are what MS:1000127 / MS:1000128
# are called, and are what comes back out of the store, so both are accepted.
SPECTRUM_TYPE_ALIASES: Dict[str, str] = {
    "centroid": SpectrumType.CENTROID,
    "centroided": SpectrumType.CENTROID,
    "centroid spectrum": SpectrumType.CENTROID,
    "profile": SpectrumType.PROFILE,
    "profile spectrum": SpectrumType.PROFILE,
}


def normalize_spectrum_type(value: Optional[object]) -> Optional[str]:
    """Map a user-supplied spectrum representation onto a :class:`SpectrumType`.

    Args:
        value: ``"profile"``, ``"centroid"``, either CV name, or ``None``.
            Case and surrounding whitespace are ignored.

    Returns:
        The matching :class:`SpectrumType` constant, or ``None`` for ``None``.

    Raises:
        ValueError: If the value is not a recognised representation. A typo
            here would otherwise be indistinguishable from "no override" and
            would silently hand the file back to auto-detection, so this fails
            loudly instead.
    """
    if value is None:
        return None

    if not isinstance(value, str):
        raise ConversionRefused(
            f"spectrum_type must be a string or None, got {type(value).__name__}"
        )

    key = value.strip().lower()
    if key in SPECTRUM_TYPE_ALIASES:
        return SPECTRUM_TYPE_ALIASES[key]

    accepted = ", ".join(repr(k) for k in sorted(SPECTRUM_TYPE_ALIASES))
    raise ConversionRefused(
        f"Unknown spectrum_type {value!r}. Accepted values: {accepted}."
    )


#: The polarity each scan polarity term states.
POLARITY_OF_ACCESSION: Dict[str, str] = {
    ImzMLAccessions.POSITIVE_SCAN: "positive",
    ImzMLAccessions.NEGATIVE_SCAN: "negative",
}

#: A scan polarity term some exporters write under the wrong accession, and
#: the accession it means. "positive scan" goes with MS:1000130, but some
#: files pair it with MS:1000128, the profile spectrum term. pyimzml corrects
#: the pair where it parses terms itself (``fileContent`` and the parameter
#: groups), not on a spectrum, so Thyra corrects it there (D55). Only the
#: name tells it apart from a true profile spectrum term.
MISLABELLED_POLARITY_TERMS: Dict[Tuple[str, str], str] = {
    (ImzMLAccessions.PROFILE_SPECTRUM, "positive scan"): ImzMLAccessions.POSITIVE_SCAN,
}


def polarity_accession(accession: str, name: Optional[str]) -> Optional[str]:
    """The scan polarity accession a cvParam states, ``None`` for none.

    Args:
        accession: The cvParam's ``accession``.
        name: The cvParam's ``name``, which decides a mislabelled term.

    Returns:
        A key of :data:`POLARITY_OF_ACCESSION`, or ``None``.
    """
    if accession in POLARITY_OF_ACCESSION:
        return accession
    # The name as written, as pyimzml matches it; a looser match here would
    # read a term the conversion route does not.
    return MISLABELLED_POLARITY_TERMS.get((accession, name or ""))


#: What an extractor records for a source that states both polarities:
#: pyimzml's word for it. No polarity term matches it, so the document's
#: field stays unset, and the file-level terms cannot fill it instead.
MIXED_POLARITY = "mixed"


def agreed_polarity(stated: Iterable[str]) -> Optional[str]:
    """The polarity every statement of a source agrees on (D35).

    ``stated`` holds ``"positive"`` or ``"negative"`` once for each place
    the source states one: its file-level terms, a parameter group, a
    spectrum. One value is the answer, and None means the source says
    nothing. Both give :data:`MIXED_POLARITY`: the source alternated or
    contradicts itself, which has no single truthful value.
    """
    values = set(stated)
    if not values:
        return None
    if len(values) > 1:
        return MIXED_POLARITY
    return next(iter(values))


#: One byte per spectrum for the representation it states: none, centroid,
#: profile, or both (a spectrum whose own terms or groups name both).
REPRESENTATION_UNSTATED = 0
REPRESENTATION_CENTROID = 1
REPRESENTATION_PROFILE = 2
REPRESENTATION_BOTH = REPRESENTATION_CENTROID | REPRESENTATION_PROFILE


class StatedSpectrumType(NamedTuple):
    """What a source's statements give for its spectrum type.

    Attributes:
        spectrum_type: :attr:`SpectrumType.CENTROID`,
            :attr:`SpectrumType.PROFILE`, or ``None``.
        stated: Whether anything was stated. ``False`` only when neither
            the spectra nor the file name a type, which is the one case a
            reader may still guess.
        warning: The one line to log, or ``None``.
    """

    spectrum_type: Optional[str]
    stated: bool
    warning: Optional[str]


def spectrum_type_from_statements(
    centroid: int, profile: int, declared: Iterable[str], source: str
) -> StatedSpectrumType:
    """The spectrum type of a source, by the one rule both routes apply.

    The type is the one every converted spectrum states. ``centroid`` and
    ``profile`` count the converted spectra stating each; a spectrum that
    states both counts in both. When the spectra state both types, no
    type is recorded. The file-level declaration (``declared``, the CV
    names it lists) counts only when no spectrum states a type, and only
    if it names exactly one. The spectra win over a declaration naming
    only the other type, with a warning.

    Args:
        centroid: Converted spectra stating centroid.
        profile: Converted spectra stating profile.
        declared: The types the file-level declaration names.
        source: The source's name, for the warning.
    """
    named = set(declared) & {SpectrumType.CENTROID, SpectrumType.PROFILE}
    if centroid and profile:
        return StatedSpectrumType(
            None,
            True,
            f"{source}: the converted spectra state both types ({centroid} "
            f"centroid, {profile} profile); no spectrum type is recorded.",
        )
    if centroid or profile:
        spectra = SpectrumType.CENTROID if centroid else SpectrumType.PROFILE
        warning = None
        if named and spectra not in named:
            warning = (
                f"{source}: the file declares {next(iter(named))!r} but the "
                f"converted spectra state {spectra!r}; recording {spectra!r}."
            )
        return StatedSpectrumType(spectra, True, warning)
    if len(named) == 1:
        return StatedSpectrumType(next(iter(named)), True, None)
    if named:
        return StatedSpectrumType(
            None,
            True,
            f"{source}: the file declares both centroid and profile spectra "
            "and no spectrum states its own; no spectrum type is recorded.",
        )
    return StatedSpectrumType(None, False, None)


class BinaryDataType:
    """Binary data type string constants."""

    CONTINUOUS = "continuous"
    PROCESSED = "processed"


#: The longest sample id or name accepted. Both are labels a person types,
#: stored once in the metadata block; a longer value is a pasted document.
MAX_SAMPLE_VALUE_LENGTH = 256


def sample_value_problem(name: str, value: Any) -> Optional[str]:
    """Why ``value`` cannot be a sample id or name, or ``None``.

    ``None`` is fine: both are optional and Thyra never fills them in.
    Anything given must be a printable, non-blank string of at most
    :data:`MAX_SAMPLE_VALUE_LENGTH` characters. ``convert_msi`` checks it
    before the reader opens, and the SpatialData converter checks it again
    for a caller that builds one directly: a value the metadata model
    refuses at write time would drop the whole ``msi_metadata`` block.

    Args:
        name: The argument's name, for the message.
        value: What the caller passed.
    """
    if value is None:
        return None
    if not isinstance(value, str):
        return f"{name} must be a string, not {type(value).__name__}"
    if not value.strip():
        return f"{name} must not be empty or blank"
    if any(not c.isprintable() for c in value):
        return f"{name} must not contain control characters or line breaks"
    if len(value) > MAX_SAMPLE_VALUE_LENGTH:
        return (
            f"{name} must be at most {MAX_SAMPLE_VALUE_LENGTH} characters, "
            f"got {len(value)}"
        )
    return None
