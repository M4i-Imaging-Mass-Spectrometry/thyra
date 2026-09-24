# thyra/metadata/schema/builder.py
"""Build the ``msi_metadata`` block from extracted metadata.

Auto-population is best-effort and honest: a field the source does not
report is left unset rather than guessed.  The only inferences made are
facts that follow from the format itself (a PHI raw file is a TOF-SIMS
acquisition) or from vendor metadata that directly encodes the fact
(a Bruker dataset with a ``MaldiFrameLaserInfo`` table is MALDI).
"""

import json
import logging
import math
import re
from datetime import datetime
from typing import Any, Dict, List, Literal, Optional, Tuple, cast

from ...core.msms import distinct_precursors
from ..personal_data import file_name as _file_name
from ..types import ComprehensiveMetadata
from .models import (
    MSI_METADATA_UNS_KEY,
    Acquisition,
    Fragmentation,
    IonMobility,
    IsolationWindow,
    MobilityGrid,
    MSAnalysis,
    MSIMetadata,
    PixelSizeUm,
    ProcessingStep,
    Provenance,
    ResolvingPower,
)
from .vocab import (
    normalize_analyzer,
    normalize_ionisation_source,
    normalize_polarity,
    term_from_accession,
)

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
#: Key spellings for who built the instrument, in preference order.  The
#: extractors disagree on the word: solariX, Waters and rapiflex write
#: ``manufacturer``, PHI writes ``vendor``.  ``manufacturer`` is
#: preferred because it is the spelling three of the four already use;
#: the PSI concept it binds to is MS:1001269 "instrument vendor", and
#: the two name the same fact.
_MANUFACTURER_KEYS = ("manufacturer", "vendor")

#: Key spellings for the instrument's serial number, in preference
#: order.  imzML (from MS:1000529) and the Bruker tsf/tdf extractor write
#: ``instrument_serial_number``; solariX and rapiflex write
#: ``serial_number``.
_SERIAL_NUMBER_KEYS = ("instrument_serial_number", "serial_number")

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


def _first_number(mapping: Dict[str, Any], keys: Tuple[str, ...]) -> Optional[float]:
    """The first finite scalar among ``keys``, if any.

    A ``[min, max]`` pair (how the Bruker extractors report a per-frame
    value that varied across the acquisition) is not a scalar and is
    skipped: the section states one value per acquisition or none. A
    numeric string counts, since the FlexImaging info file is parsed as
    text.
    """
    for key in keys:
        value = mapping.get(key)
        if isinstance(value, bool):
            continue
        if isinstance(value, str):
            try:
                value = float(value.strip())
            except ValueError:
                continue
        if isinstance(value, (int, float)) and math.isfinite(value):
            return float(value)
    return None


# --- the ``acquisition`` section -------------------------------------------
#
# Key spellings per field, in preference order.  Each is a name an
# extractor already writes into ``acquisition_params``; the raw dict is
# read here and never rewritten, so the vendor spelling stays available
# beside the normalised one.  The vendor evidence behind each unit is in
# the extractor that writes the key.
_ACQUISITION_DATETIME_KEYS = (
    "acquisition_datetime",  # Bruker tsf/tdf and solariX, ISO 8601 with offset
    "acquisition_date",  # PHI (AcqFileDate) and Waters (MassLynx), no offset
)
_LASER_POWER_KEYS = ("laser_power",)  # Bruker tsf/tdf, solariX, rapiflex
_LASER_FREQUENCY_KEYS = (
    "laser_frequency",  # Bruker tsf/tdf: MaldiFrameInfo.LaserRepRate, Hz
    "laser_rep_rate",  # solariX: Spectra.LaserRepRate, Hz
)
_SHOTS_PER_PIXEL_KEYS = (
    "num_laser_shots",  # Bruker tsf/tdf: MaldiFrameInfo.NumLaserShots
    "num_summations",  # solariX: Spectra.NumSummations, which is the shot count
    "shots_per_spot",  # rapiflex: "Number of Shots" in the info file
)
_METHOD_FILE_KEYS = (
    "method_name",  # Bruker tsf/tdf (GlobalMetadata.MethodName) and solariX (*.m)
    "ms_method",  # Waters: "$$ MS Method" in _header.txt
    "method",  # rapiflex: "Method" in the info file
)

# Formats whose ``laser_power`` is verified to be a percentage of the
# laser's range: Bruker tsf/tdf (``MaldiFrameInfo.LaserPower``, the value
# timsControl shows as "Laser Power %") and solariX (``Spectra.LaserPower``
# equals the method's ``LaserAttn`` on a ``LaserPowerRange`` of 100).  The
# rapiflex info file's "Laser Power" was not checked against an
# acquisition, so it is not trusted as a percentage; imzML never carries
# the key.
_LASER_POWER_PERCENT_FORMATS = frozenset({"bruker", "tsf", "tdf", "solarix"})

_MONTH_ABBREVIATIONS = {
    name: number
    for number, name in enumerate(
        (
            "jan",
            "feb",
            "mar",
            "apr",
            "may",
            "jun",
            "jul",
            "aug",
            "sep",
            "oct",
            "nov",
            "dec",
        ),
        start=1,
    )
}
# Waters MassLynx ``getAcquisitionDate``: '07-Nov-2019 14:09:44'.
_WATERS_DATETIME = re.compile(
    r"(?P<day>\d{1,2})-(?P<month>[A-Za-z]{3})-(?P<year>\d{4})"
    r"\s+(?P<hour>\d{1,2}):(?P<minute>\d{2}):(?P<second>\d{2})"
)
# PHI SmartSoft-TOF ``AcqFileDate``: '06/23/2026 21:12:35', month first.
_PHI_DATETIME = re.compile(
    r"(?P<month>\d{1,2})/(?P<day>\d{1,2})/(?P<year>\d{4})"
    r"\s+(?P<hour>\d{1,2}):(?P<minute>\d{2}):(?P<second>\d{2})"
)


def _parse_vendor_datetime(value: str) -> Optional[datetime]:
    """Parse the timestamp formats the readers report, or ``None``.

    Three formats, each tied to the vendor that writes it: ISO 8601 with
    a UTC offset (Bruker tsf/tdf ``GlobalMetadata.AcquisitionDateTime``
    and solariX ``Properties.AcquisitionDateTime``), MassLynx's
    ``dd-Mon-yyyy hh:mm:ss`` (Waters) and SmartSoft-TOF's
    ``mm/dd/yyyy hh:mm:ss`` (PHI).  They cannot be confused with each
    other, so the order does not matter.  The month names are matched
    against a fixed table rather than ``strptime``'s ``%b``, whose
    reading depends on the process locale.  Anything else stays
    unparsed: the raw string is still in ``acquisition_params``.
    """
    text = value.strip()
    if "T" in text:
        try:
            return datetime.fromisoformat(text)
        except ValueError:
            return None
    for pattern in (_WATERS_DATETIME, _PHI_DATETIME):
        match = pattern.fullmatch(text)
        if match is None:
            continue
        parts = match.groupdict()
        month_text = parts["month"]
        month = (
            _MONTH_ABBREVIATIONS.get(month_text.lower())
            if month_text.isalpha()
            else int(month_text)
        )
        if month is None:
            return None
        try:
            return datetime(
                int(parts["year"]),
                month,
                int(parts["day"]),
                int(parts["hour"]),
                int(parts["minute"]),
                int(parts["second"]),
            )
        except ValueError:
            return None
    return None


def _format_iso_8601(moment: datetime) -> str:
    """ISO 8601 text keeping the source's precision and offset, if any.

    Whole seconds stay whole and milliseconds stay milliseconds, so a
    Bruker ``2025-04-22T08:59:34.395+02:00`` is written back exactly.
    """
    if moment.microsecond == 0:
        return moment.isoformat(timespec="seconds")
    if moment.microsecond % 1000 == 0:
        return moment.isoformat(timespec="milliseconds")
    return moment.isoformat()


def _build_acquisition(
    acquisition: Dict[str, Any], source_format: Optional[str]
) -> Optional[Acquisition]:
    """The ``acquisition`` section from what the extractor reported.

    Every field follows the builder's rule -- best-effort and honest: a
    vendor value is taken only where its meaning and unit are verified
    (see the key tables above), an unparseable timestamp leaves the
    field unset with the raw string untouched in ``acquisition_params``,
    and a per-frame value that varied across the acquisition (reported
    as a pair, not a scalar) is not summarised into one number.

    Returns ``None`` when nothing was reported, so the section is absent
    from the block rather than written empty.
    """
    fields: Dict[str, Any] = {}

    raw_datetime = _first_string(acquisition, _ACQUISITION_DATETIME_KEYS)
    if raw_datetime is not None:
        moment = _parse_vendor_datetime(raw_datetime)
        if moment is not None:
            fields["acquisition_datetime"] = _format_iso_8601(moment)
        else:
            logger.debug(
                "Acquisition timestamp %r is in no format the builder parses; "
                "left unset",
                raw_datetime,
            )

    if (source_format or "").lower() in _LASER_POWER_PERCENT_FORMATS:
        power = _first_number(acquisition, _LASER_POWER_KEYS)
        if power is not None and 0.0 <= power <= 100.0:
            fields["laser_power_percent"] = power

    frequency = _first_number(acquisition, _LASER_FREQUENCY_KEYS)
    if frequency is not None and frequency > 0.0:
        fields["laser_frequency_hz"] = frequency

    shots = _first_number(acquisition, _SHOTS_PER_PIXEL_KEYS)
    if shots is not None and shots >= 1.0 and shots == int(shots):
        fields["shots_per_pixel"] = int(shots)

    method = _first_string(acquisition, _METHOD_FILE_KEYS)
    if method is not None:
        name = _file_name(method)
        if name is not None:
            fields["method_file"] = name

    if not fields:
        return None
    return Acquisition(**fields)


def _polarity_from_cv_params(raw_metadata: Dict[str, Any]) -> Optional[str]:
    """Polarity declared by the raw file's own cvParams, if unambiguous.

    imzML declares polarity as MS:1000130 (positive scan) / MS:1000129
    (negative scan) in the file description; the extractor preserves
    those with their accessions.  A file declaring both (alternating
    polarity) has no single truthful value and returns ``None``.
    """
    cv_params = raw_metadata.get("cvParams")
    if not isinstance(cv_params, list):
        return None
    accessions = {
        entry.get("accession") for entry in cv_params if isinstance(entry, dict)
    }
    positive = "MS:1000130" in accessions
    negative = "MS:1000129" in accessions
    if positive == negative:
        return None
    return "positive" if positive else "negative"


def _build_instrument_fields(
    acquisition: Dict[str, Any],
    instrument: Dict[str, Any],
    format_specific: Dict[str, Any],
    source_format: Optional[str],
) -> Dict[str, Any]:
    """What instrument produced the data: source, analyzer and model.

    Grouped because all three resolve the same way -- what the extractor
    reported, then what the format itself implies (a PHI raw file is a
    TOF-SIMS acquisition), and nothing when neither says. The per-format
    defaults are consulted by all of them, which is what makes this one
    step rather than three.
    """
    fields: Dict[str, Any] = {}
    fmt_defaults = _FORMAT_DEFAULTS.get((source_format or "").lower(), {})

    source = _resolve_source(acquisition, format_specific, fmt_defaults)
    if source is not None:
        fields["ionisation_source"], fields["ionisation_source_term"] = source

    analyzer = _resolve_analyzer(acquisition, instrument, fmt_defaults)
    if analyzer is not None:
        fields["analyzer"], fields["analyzer_term"] = analyzer

    instrument_model = _first_string(instrument, _INSTRUMENT_MODEL_KEYS)
    if instrument_model is not None:
        fields["instrument_model"] = instrument_model

    manufacturer = _first_string(instrument, _MANUFACTURER_KEYS)
    if manufacturer is not None:
        fields["manufacturer"] = manufacturer

    serial_number = _first_string(instrument, _SERIAL_NUMBER_KEYS)
    if serial_number is not None:
        fields["serial_number"] = serial_number

    resolving_power = _build_resolving_power(acquisition, instrument)
    if resolving_power is not None:
        fields["detector_resolving_power"] = resolving_power

    return fields


def _resolve_source(
    acquisition: Dict[str, Any],
    format_specific: Dict[str, Any],
    fmt_defaults: Dict[str, Any],
) -> Optional[Tuple[str, Any]]:
    """The ionisation source: reported, else implied by the format."""
    reported = _first_string(
        acquisition, ("ionisation_source", "ion_source", "technique")
    )
    source = normalize_ionisation_source(reported)
    if reported is not None and source is None:
        # "Unset beats guessed" is the rule, but a spelling the alias table
        # does not know should be findable in a log rather than only as a
        # field that is quietly missing (issue #388).
        logger.debug(
            "Ionisation source %r matches no known spelling; left unset", reported
        )
    if source is None and format_specific.get("is_maldi"):
        source = normalize_ionisation_source("maldi")
    if source is None and "ionisation_source" in fmt_defaults:
        source = normalize_ionisation_source(fmt_defaults["ionisation_source"])
    return source


def _resolve_analyzer(
    acquisition: Dict[str, Any],
    instrument: Dict[str, Any],
    fmt_defaults: Dict[str, Any],
) -> Optional[Tuple[str, Any]]:
    """The mass analyzer: reported, else implied by the format."""
    reported = _first_string(
        instrument, ("analyzer", "mass_analyzer")
    ) or _first_string(acquisition, ("analyzer", "mass_analyzer"))
    analyzer = normalize_analyzer(reported)
    if reported is not None and analyzer is None:
        logger.debug("Analyzer %r matches no known spelling; left unset", reported)
    if analyzer is None and "analyzer" in fmt_defaults:
        analyzer = normalize_analyzer(fmt_defaults["analyzer"])
    return analyzer


# The two keys an extractor writes when the source states its resolving
# power: the value and the m/z it is quoted at.  Both are needed -- a
# resolving power without its reference m/z is not comparable between
# analyzers, which is why the schema field is the pair.  Probed in
# ``instrument_info`` first, then ``acquisition_params``.  No shipped
# extractor writes them yet: the Waters extractor's ``declared_resolution``
# is not mapped here because what MassLynx means by it has not been
# verified against an acquisition, and a Bruker method's resolution setting
# is a mode, not a number.  A Thermo scan trailer states both directly.
_RESOLVING_POWER_KEYS = ("resolving_power",)
_RESOLVING_POWER_AT_MZ_KEYS = ("resolving_power_at_mz",)


def _build_resolving_power(
    acquisition: Dict[str, Any], instrument: Dict[str, Any]
) -> Optional[ResolvingPower]:
    """The stated resolving power and its reference m/z, or ``None``."""
    for mapping in (instrument, acquisition):
        value = _first_number(mapping, _RESOLVING_POWER_KEYS)
        at_mz = _first_number(mapping, _RESOLVING_POWER_AT_MZ_KEYS)
        if value is None and at_mz is None:
            continue
        if value is None or at_mz is None or value <= 0.0 or at_mz <= 0.0:
            logger.debug(
                "Resolving power reported without a usable value and reference "
                "m/z (%r at %r); left unset",
                value,
                at_mz,
            )
            return None
        return ResolvingPower(value=value, at_mz=at_mz)
    return None


def _build_ms_analysis(
    acquisition: Dict[str, Any],
    instrument: Dict[str, Any],
    format_specific: Dict[str, Any],
    raw_metadata: Dict[str, Any],
    pixel_size_um: Tuple[float, float],
    source_format: Optional[str],
    mobility_resolved_table: Optional[str] = None,
    mobility_grid: Optional[Dict[str, Any]] = None,
    fragmentation: Any = None,
    msms_resolved_table: Optional[str] = None,
) -> MSAnalysis:
    """Assemble the acquisition section from what the extractors report."""
    fields: Dict[str, Any] = {}

    polarity = normalize_polarity(
        acquisition.get("polarity") or _polarity_from_cv_params(raw_metadata)
    )
    if polarity is not None:
        fields["polarity"], fields["polarity_term"] = polarity

    fields.update(
        _build_instrument_fields(
            acquisition, instrument, format_specific, source_format
        )
    )

    ion_mobility = _build_ion_mobility(
        format_specific.get("ion_mobility"), mobility_resolved_table, mobility_grid
    )
    if ion_mobility is not None:
        fields["ion_mobility"] = ion_mobility

    fragmentation_block = _build_fragmentation(fragmentation, msms_resolved_table)
    if fragmentation_block is not None:
        fields["fragmentation"] = fragmentation_block

    return MSAnalysis(
        pixel_size_um=PixelSizeUm(x=pixel_size_um[0], y=pixel_size_um[1]),
        **fields,
    )


def _optional_term(accession: Any) -> Optional[Any]:
    """The ontology term for an accession the extractor reported, if resolvable."""
    if not isinstance(accession, str) or not accession:
        return None
    try:
        return term_from_accession(accession)
    except KeyError:
        logger.debug("Mobility accession %s is not in the local ontology", accession)
        return None


def _build_ion_mobility(
    reported: Any,
    resolved_table: Optional[str] = None,
    grid: Optional[Dict[str, Any]] = None,
) -> Optional[IonMobility]:
    """The mobility block from what a reader's extractor reported.

    The Bruker and imzML extractors report one (``present`` True for TDF
    and for an imzML with a mobility array, False for TSF); readers that
    say nothing leave the field unset, which is honest -- "not reported"
    is not the same as "no mobility". A resolved table written beside
    the summed one is named here whatever the extractor said, since its
    existence proves the dimension.

    ``grid`` is present only when that table was *binned* onto a common
    mobility grid rather than read off a shared feature axis. It is the
    one thing in the store that says which of the two mechanisms filled
    the table, and it is a description: the table itself is the same
    shape either way.
    """
    if not isinstance(reported, dict) or "present" not in reported:
        if resolved_table:
            return IonMobility(
                present=True,
                resolved_table=resolved_table,
                grid=_mobility_grid(grid),
            )
        return None
    present = bool(reported["present"])
    if not present and not resolved_table:
        return IonMobility(present=False)

    fields: Dict[str, Any] = {"present": True}
    if resolved_table:
        fields["resolved_table"] = resolved_table
    grid_block = _mobility_grid(grid)
    if grid_block is not None:
        fields["grid"] = grid_block
    fields.update(_mobility_axis_fields(reported))
    return IonMobility(**fields)


def _mobility_grid(reported: Any) -> Optional[MobilityGrid]:
    """The grid block from what the converter resolved, or ``None``.

    Checked for shape rather than trusted, like everything else here: a
    grid that will not validate is dropped, since an invented one would
    let a consumer map a heatmap box onto channels that do not exist.
    """
    if not isinstance(reported, dict):
        return None
    try:
        return MobilityGrid(
            law=str(reported["law"]),
            lower=float(reported["lower"]),
            upper=float(reported["upper"]),
            n_channels=int(reported["n_channels"]),
        )
    except (KeyError, TypeError, ValueError) as e:
        logger.debug("Mobility grid block is not usable and was dropped: %s", e)
        return None


def _mobility_axis_fields(reported: Dict[str, Any]) -> Dict[str, Any]:
    """The axis description an extractor reported, in the block's field names.

    Everything is optional and checked for shape rather than trusted: an
    unresolvable accession is dropped, never invented, and a malformed
    range or scan count is ignored.
    """
    fields: Dict[str, Any] = {}
    separation = reported.get("separation")
    if isinstance(separation, str) and separation.strip():
        fields["separation"] = separation.strip()
    for field, source in (
        ("separation_term", "separation_accession"),
        ("unit_term", "unit_accession"),
    ):
        term = _optional_term(reported.get(source))
        if term is not None:
            fields[field] = term
    mobility_range = reported.get("one_over_k0_range") or reported.get("range")
    if isinstance(mobility_range, (list, tuple)) and len(mobility_range) == 2:
        try:
            fields["range_lower"] = float(mobility_range[0])
            fields["range_upper"] = float(mobility_range[1])
        except (TypeError, ValueError):
            pass
    num_scans = reported.get("num_scans_max", reported.get("num_scans"))
    if isinstance(num_scans, (int, float)) and num_scans >= 1:
        fields["num_scans"] = int(num_scans)
    return fields


def _build_fragmentation(
    reported: Any, resolved_table: Optional[str] = None
) -> Optional[Fragmentation]:
    """The fragmentation block from what a reader reported.

    ``None`` in, ``None`` out: a reader that cannot tell says nothing,
    and an unset block is honest where ``present=False`` would be a
    claim. Everything is checked for shape rather than trusted -- a
    window without a usable target m/z is dropped rather than invented,
    since a precursor list is exactly the thing a consumer would act on.

    A demultiplexed table written beside the summed one is named here,
    the way :func:`_build_ion_mobility` names the mobility sibling, so
    both kinds are discoverable from this versioned block alone.
    """
    if not isinstance(reported, dict) or "ms_level" not in reported:
        return None
    try:
        ms_level = int(reported["ms_level"])
    except (TypeError, ValueError):
        return None
    if ms_level < 1:
        return None

    present = bool(reported.get("present", ms_level > 1))
    if not present:
        return Fragmentation(present=False, ms_level=1)

    windows = [
        window
        for window in (
            _build_isolation_window(entry) for entry in reported.get("windows") or []
        )
        if window is not None
    ]
    term = _optional_term(reported.get("dissociation_accession"))
    return Fragmentation(
        present=True,
        ms_level=max(ms_level, 2),
        constant_across_pixels=bool(reported.get("constant_across_pixels", True)),
        merges_precursors=distinct_precursors(windows) > 1,
        dissociation_term=term if windows else None,
        windows=windows,
        resolved_table=resolved_table,
    )


def _build_isolation_window(entry: Any) -> Optional[IsolationWindow]:
    """One isolation window, or ``None`` when it carries no usable target."""
    if not isinstance(entry, dict):
        return None
    try:
        target = float(entry["isolation_window_target"])
    except (KeyError, TypeError, ValueError):
        return None
    if not target > 0:
        return None

    fields: Dict[str, Any] = {"target": target}
    for field, key in (
        ("lower_offset", "isolation_window_lower_offset"),
        ("upper_offset", "isolation_window_upper_offset"),
        ("collision_energy", "collision_energy"),
    ):
        value = entry.get(key)
        if isinstance(value, (int, float)):
            fields[field] = float(value)
    begin, end = entry.get("scan_begin"), entry.get("scan_end")
    if isinstance(begin, int) and isinstance(end, int) and end > begin >= 0:
        fields["scan_begin"] = begin
        fields["scan_end"] = end
    return IsolationWindow(**fields)


def build_msi_metadata(
    comprehensive: Optional[ComprehensiveMetadata],
    *,
    pixel_size_um: Tuple[float, float],
    pixel_size_source: Optional[str] = None,
    source_format: Optional[str] = None,
    processing: Optional[List[ProcessingStep]] = None,
    mobility_resolved_table: Optional[str] = None,
    mobility_grid: Optional[Dict[str, Any]] = None,
    fragmentation: Any = None,
    msms_resolved_table: Optional[str] = None,
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
        processing: Ordered processing steps performed so far, oldest
            first (see :class:`ProcessingStep`).
        mobility_resolved_table: Element key of the mobility-resolved
            sibling table written beside the summed table, when one was.
        mobility_grid: The common mobility grid that table was binned
            onto, as
            :meth:`thyra.resampling.mobility_grid.MobilityGrid.to_schema_report`
            renders it. ``None`` for a table read off a shared feature
            axis, which was binned onto nothing.
        fragmentation: What the reader reported about fragmentation, as
            :meth:`thyra.core.msms.FragmentationSchedule.to_extractor_report`
            renders it. ``None`` means the reader did not say, which is
            not the same as "MS1" and leaves the block unset.
        msms_resolved_table: Element key of the demultiplexed MS/MS
            sibling table written beside the summed table, when one was.

    Returns:
        The populated document.  Fields the source does not report are
        left unset, and the ``acquisition`` section is absent when the
        reader reported none of its facts.
    """
    from thyra import __version__

    acquisition: Dict[str, Any] = {}
    instrument: Dict[str, Any] = {}
    format_specific: Dict[str, Any] = {}
    raw_metadata: Dict[str, Any] = {}
    source_path: Optional[str] = None
    if comprehensive is not None:
        acquisition = dict(comprehensive.acquisition_params or {})
        instrument = dict(comprehensive.instrument_info or {})
        format_specific = dict(comprehensive.format_specific or {})
        raw_metadata = dict(comprehensive.raw_metadata or {})
        essential = comprehensive.essential
        if essential is not None:
            source_path = str(essential.source_path)

    return MSIMetadata(
        ms_analysis=_build_ms_analysis(
            acquisition,
            instrument,
            format_specific,
            raw_metadata,
            pixel_size_um,
            source_format,
            mobility_resolved_table,
            mobility_grid,
            fragmentation,
            msms_resolved_table,
        ),
        acquisition=_build_acquisition(acquisition, source_format),
        processing=list(processing or []),
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


# A pitch to get ``MSAnalysis`` built when the source states none. It is
# removed from the document in the same breath, inside
# :func:`build_metadata_document`, and never reaches a caller: the field
# is mandatory in schema 0.6.0 and the model cannot be instantiated
# without one, which is precisely the constraint a non-imaging document
# runs into.
_PITCH_STAND_IN = (1.0, 1.0)


def _name_the_source(document: Dict[str, Any]) -> None:
    """Reduce ``provenance.source_path`` to the source's name, in place.

    A block inside a store keeps the whole path.  The store sits on the
    machine that wrote it, the two are usually neighbours, and which
    folder an element was converted from is provenance a reader of that
    store can act on.

    A document is the opposite on every count: it is small, it is
    portable, and ``thyra metadata`` exists to produce one that can be
    handed to somebody else.  An absolute path in it describes a
    filesystem the recipient does not have, while carrying a user
    directory and the acquisition folder's name to a machine that has no
    use for either.  The name is the part a recipient can reconcile
    against the acquisition; the rest is the machine it was read on.

    Same reduction as ``acquisition.method_file`` and for the same
    reason (see :func:`_file_name`).  It runs here rather than in
    :func:`build_msi_metadata` because the difference is not in the
    field, it is in what is being written: only a document is meant to
    travel.

    Trailing separators are stripped first, which :func:`_file_name`
    does not do for a method file and does not need to.  Two of the
    formats Thyra reads name a *directory* -- Bruker ``.d`` and Waters
    ``.raw`` -- so a source path ending in a separator is ordinary here
    where it never is for a file.
    """
    # Indexed, not probed: ``provenance`` is required on the model and
    # ``thyra_version`` is required within it, so the section always
    # survives ``exclude_none`` and ``to_uns_dict`` never drops it.  A
    # guard here would be unreachable, and for a step whose whole job is
    # to take something out it would be the wrong failure mode anyway --
    # a shape this does not recognise should stop the document being
    # written, not silently leave the path in it.
    provenance = document["provenance"]
    source = provenance.get("source_path")
    if not isinstance(source, str):
        # Unset, which is what a source Thyra could not name looks like.
        return
    name = _file_name(source.rstrip("/\\"))
    if name is None:
        provenance.pop("source_path", None)
    else:
        provenance["source_path"] = name


def build_metadata_document(
    comprehensive: Optional[ComprehensiveMetadata],
    *,
    pixel_size_um: Optional[Tuple[float, float]],
    pixel_size_source: Optional[str] = None,
    source_format: Optional[str] = None,
    processing: Optional[List[ProcessingStep]] = None,
    fragmentation: Any = None,
) -> Dict[str, Any]:
    """The ``msi_metadata`` document for a source that was not converted.

    The same shape :func:`~thyra.metadata.schema.store_io.read_msi_metadata_blocks`
    hands back for a converted store, so one dataset's metadata reads
    the same whether it came from the raw file or from the store: real
    lists for ``processing`` and for the isolation windows, where
    :meth:`MSIMetadata.to_uns_dict` packs both into JSON strings because
    AnnData/zarr cannot round-trip a list of objects.

    Args:
        comprehensive: The reader's comprehensive metadata, or ``None``.
        pixel_size_um: The in-plane raster pitch, or ``None`` when the
            source states none.  Nothing that writes a store passes
            ``None``: a conversion refuses without a pitch, so every
            stored block has one.  A metadata-only read has no such
            guarantee -- a Bruker ``.d`` that imaged nothing has no
            raster and no pitch -- and the field is then left out rather
            than filled with a number nobody measured.  Schema 0.6.0
            makes it mandatory, so such a document does not validate;
            that is the finding this command exists to show, not a
            defect in the document.
        pixel_size_source: How the pitch was determined, when it was.
        source_format: Detected input format name, when known.
        processing: Steps performed so far; empty for an unconverted
            source, which has had nothing done to it.
        fragmentation: What the reader reported about fragmentation, as
            :meth:`thyra.core.msms.FragmentationSchedule.to_extractor_report`
            renders it.

    Returns:
        The document as a plain dict.  Sibling-table fields are absent:
        no sibling was written, because nothing was written.
    """
    from .store_io import decode_isolation_windows

    meta = build_msi_metadata(
        comprehensive,
        pixel_size_um=pixel_size_um or _PITCH_STAND_IN,
        pixel_size_source=pixel_size_source,
        source_format=source_format,
        processing=processing,
        fragmentation=fragmentation,
    )
    document = meta.to_uns_dict()
    if pixel_size_um is None:
        document.get("ms_analysis", {}).pop("pixel_size_um", None)
    _name_the_source(document)

    stored = document.get("processing")
    if isinstance(stored, str):
        document["processing"] = json.loads(stored)
    decode_isolation_windows(document, MSI_METADATA_UNS_KEY)
    return document


# The sibling tables named in ``ms_analysis``, and the fields of each
# section that describe the sibling rather than the acquisition.
#
# ``ion_mobility.grid`` is here with ``resolved_table`` because it
# describes the binning *of that table*: it is present only when the
# table was binned onto a common mobility grid rather than read off a
# shared feature axis (see :func:`_build_ion_mobility`). Keeping it after
# the table is unnamed would say a table was binned onto a grid that
# nobody wrote.
_SIBLING_TABLE_FIELDS: Dict[str, Tuple[str, ...]] = {
    "ion_mobility": ("resolved_table", "grid"),
    "fragmentation": ("resolved_table",),
}


def forget_resolved_table(meta_uns: Dict[str, Any], section: str) -> None:
    """Unname the sibling table ``ms_analysis.<section>`` points at.

    The converter names a sibling table in the summed table's metadata
    *before* the sibling is built, so that the two agree; a builder that
    then declines by returning ``None`` -- a decision, not a failure --
    would otherwise leave the block pointing at an element nobody wrote
    (issue #343).

    Takes the serialised block, not the model: by the time a builder can
    decline, :meth:`MSIMetadata.to_uns_dict` has already run and the
    ``uns`` entry is a plain dict.  The path into it lives here, beside
    the builder that wrote it, rather than in the converter.

    The section itself is left in place.  ``ion_mobility.present`` and
    the fragmentation windows describe the *acquisition*, which is no
    less true for the sibling not being written, and the converter only
    ever names a sibling for a source that has the dimension.
    """
    analysis = meta_uns.get("ms_analysis")
    if not isinstance(analysis, dict):
        return
    block = analysis.get(section)
    if not isinstance(block, dict):
        return
    for field in _SIBLING_TABLE_FIELDS[section]:
        block.pop(field, None)
