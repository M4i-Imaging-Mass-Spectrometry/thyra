# thyra/metadata/schema/builder.py
"""Build the ``msi_metadata`` block from extracted metadata.

Auto-population is best-effort and honest: a field the source does not
report is left unset rather than guessed.  The only inferences made are
facts that follow from the format itself (a PHI raw file is a TOF-SIMS
acquisition) or from vendor metadata that directly encodes the fact
(a Bruker dataset with a ``MaldiFrameLaserInfo`` table is MALDI).
"""

import logging
from typing import Any, Dict, Literal, Optional, Tuple, cast

from ..types import ComprehensiveMetadata
from .models import MSAnalysis, MSIMetadata, PixelSizeUm, Provenance
from .vocab import normalize_analyzer, normalize_ionisation_source, normalize_polarity

logger = logging.getLogger(__name__)

# Facts that follow from the source format itself.  Kept deliberately
# conservative: only entries where every dataset of that format shares
# the value.  PHI raw files come from TOF-SIMS instruments; the Bruker
# formats Thyra reads (.d with analysis.tsf/.tdf) are all
# timsTOF-family, i.e. TOF analyzers.
_FORMAT_DEFAULTS: Dict[str, Dict[str, str]] = {
    "phi": {"ionisation_source": "SIMS", "analyzer": "TOF"},
    "bruker": {"analyzer": "TOF"},
    "tsf": {"analyzer": "TOF"},
    "tdf": {"analyzer": "TOF"},
}

# Key spellings the extractors use for the instrument model, in
# preference order (imzML: instrument_model; Bruker: instrument_name /
# model; PHI: platform).
_INSTRUMENT_MODEL_KEYS = (
    "instrument_model",
    "instrument_name",
    "model",
    "platform",
)


def _first_string(mapping: Dict[str, Any], keys: Tuple[str, ...]) -> Optional[str]:
    """The first non-empty string value among ``keys``, if any."""
    for key in keys:
        value = mapping.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return None


def _build_ms_analysis(
    acquisition: Dict[str, Any],
    instrument: Dict[str, Any],
    format_specific: Dict[str, Any],
    pixel_size_um: Tuple[float, float],
    source_format: Optional[str],
) -> MSAnalysis:
    """Assemble the acquisition section from what the extractors report."""
    fields: Dict[str, Any] = {}

    polarity = normalize_polarity(acquisition.get("polarity"))
    if polarity is not None:
        fields["polarity"], fields["polarity_term"] = polarity

    source = normalize_ionisation_source(
        _first_string(acquisition, ("ionisation_source", "ion_source", "technique"))
    )
    if source is None and format_specific.get("is_maldi"):
        source = normalize_ionisation_source("maldi")
    fmt_defaults = _FORMAT_DEFAULTS.get((source_format or "").lower(), {})
    if source is None and "ionisation_source" in fmt_defaults:
        source = normalize_ionisation_source(fmt_defaults["ionisation_source"])
    if source is not None:
        fields["ionisation_source"], fields["ionisation_source_term"] = source

    analyzer = normalize_analyzer(
        _first_string(instrument, ("analyzer", "mass_analyzer"))
        or _first_string(acquisition, ("analyzer", "mass_analyzer"))
    )
    if analyzer is None and "analyzer" in fmt_defaults:
        analyzer = normalize_analyzer(fmt_defaults["analyzer"])
    if analyzer is not None:
        fields["analyzer"], fields["analyzer_term"] = analyzer

    instrument_model = _first_string(instrument, _INSTRUMENT_MODEL_KEYS)
    if instrument_model is not None:
        fields["instrument_model"] = instrument_model

    return MSAnalysis(
        pixel_size_um=PixelSizeUm(x=pixel_size_um[0], y=pixel_size_um[1]),
        **fields,
    )


def build_msi_metadata(
    comprehensive: Optional[ComprehensiveMetadata],
    *,
    pixel_size_um: Tuple[float, float],
    pixel_size_source: Optional[str] = None,
    source_format: Optional[str] = None,
) -> MSIMetadata:
    """Build an :class:`MSIMetadata` document from extracted metadata.

    Args:
        comprehensive: The reader's comprehensive metadata, or ``None``
            when unavailable -- the document is still built, carrying
            the pixel size and provenance.
        pixel_size_um: Resolved in-plane pixel pitch ``(x_um, y_um)``.
            Required: conversion refuses to run without one, so a block
            without it describes no store Thyra ever wrote.
        pixel_size_source: How the pixel size was determined
            (``"automatic"`` / ``"manual"`` / ``"default"``).
        source_format: Detected input format name (``"imzml"``,
            ``"bruker"``, ...), when known.

    Returns:
        The populated document.  Fields the source does not report are
        left unset.
    """
    from thyra import __version__

    acquisition: Dict[str, Any] = {}
    instrument: Dict[str, Any] = {}
    format_specific: Dict[str, Any] = {}
    source_path: Optional[str] = None
    if comprehensive is not None:
        acquisition = dict(comprehensive.acquisition_params or {})
        instrument = dict(comprehensive.instrument_info or {})
        format_specific = dict(comprehensive.format_specific or {})
        essential = comprehensive.essential
        if essential is not None:
            source_path = str(essential.source_path)

    return MSIMetadata(
        ms_analysis=_build_ms_analysis(
            acquisition, instrument, format_specific, pixel_size_um, source_format
        ),
        provenance=Provenance(
            thyra_version=__version__,
            source_format=source_format,
            source_path=source_path,
            pixel_size_source=cast(
                Optional[Literal["default", "manual", "automatic"]],
                pixel_size_source,
            ),
        ),
    )
