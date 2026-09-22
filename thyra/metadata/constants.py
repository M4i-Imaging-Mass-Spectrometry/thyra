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

from typing import Dict, Optional

from ..errors import ConversionRefused


class ImzMLAccessions:
    """PSI-MS and imzML controlled vocabulary accession codes."""

    # Spectrum type (MS ontology)
    CENTROID_SPECTRUM = "MS:1000127"
    PROFILE_SPECTRUM = "MS:1000128"

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


class BinaryDataType:
    """Binary data type string constants."""

    CONTINUOUS = "continuous"
    PROCESSED = "processed"
