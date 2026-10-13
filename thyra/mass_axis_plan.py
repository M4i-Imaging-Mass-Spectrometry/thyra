"""Plan one mass axis for several sources (issue #449, design decision D43).

Two stores can be compared bin for bin only if their ``var["mz"]`` are the
same array. :func:`plan_mass_axis` reads each source's header (the same
metadata-only read :func:`thyra.preview_msi` does; no spectrum is decoded),
decides the axis settings the sources share, builds that axis and keeps its
fingerprint (D42) in a :class:`~thyra.resampling.plan.MassAxisPlan`. A
conversion given the plan builds the same axis, and refuses when its own
hashes differently.

What a plan fixes is the axis and nothing else. The axis follows from the
range, the axis type, the bin count and, for ``tof``, the width law; the
resampling method does not enter it. Each source therefore keeps its own
detector's method unless the caller forces one, so a centroid source and a
profile source can share an axis.

The decision is split in two. :func:`_plan_from_metadata` takes what the
detector chain is told about each source and reads no file, so every rule
can be exercised without one. :func:`plan_mass_axis` reads the sources and
calls it.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from .core.registry import detect_format
from .errors import ConversionRefused
from .metadata.schema import mass_axis_sha256
from .preview import _build_reader, _close_quietly, _resampling_metadata_dict
from .resampling.axis_planner import (
    AXIS_TYPE_BY_NAME,
    METHOD_BY_NAME,
    AxisPlanner,
    _resolve_config_enum,
    axis_name,
    normalize_resampling_config,
    reference_params,
)
from .resampling.decision_tree import ResamplingDecisionTree
from .resampling.mass_axis.tof_generator import DEFAULT_BINS_PER_FWHM, TOFAxisGenerator
from .resampling.plan import MassAxisPlan
from .resampling.types import (
    DEFAULT_REFERENCE_MZ,
    AxisType,
    ResamplingConfig,
    ResamplingMethod,
)

logger = logging.getLogger(__name__)

#: How a bin's width grows with m/z on each axis law, as the exponent ``p`` of
#: ``width ~ m^p``. Read off the generators in ``thyra/resampling/mass_axis``:
#: a width ``w`` stated at ``r`` is ``w * (R / r) ** p`` at ``R``.
_WIDTH_EXPONENT = {
    "constant": 0.0,
    "linear_tof": 0.5,
    "reflector_tof": 1.0,
    "orbitrap": 1.5,
    "fticr": 2.0,
}

#: Relative slack when asking whether one width is narrower than another, so
#: two widths that agree to rounding do not warn.
_WIDTH_TOLERANCE = 1e-9


@dataclass(frozen=True)
class _Source:
    """What the planner knows about one source, with no file behind it.

    Attributes:
        name: The source's name in messages (a file name, not a path).
        metadata: The dict the detector chain reads; its
            ``essential_metadata`` entry carries ``mass_range`` and
            ``spectrum_type``.
        range_known: Whether the source states a mass range at all.
        source_format: The registry format name (``"imzml"``, ``"phi"``,
            ...), or ``None`` when unknown.
    """

    name: str
    metadata: Dict[str, Any]
    range_known: bool = True
    source_format: Optional[str] = None

    @property
    def mass_range(self) -> Tuple[float, float]:
        """The source's own ``(min_mz, max_mz)``; ``(0, 0)`` when unknown."""
        raw = (self.metadata.get("essential_metadata") or {}).get("mass_range")
        if raw is None:
            return 0.0, 0.0
        return float(raw[0]), float(raw[1])

    @property
    def spectrum_type(self) -> Optional[str]:
        """The source's spectrum type string, if it states one."""
        value = (self.metadata.get("essential_metadata") or {}).get("spectrum_type")
        return None if value is None else str(value)


def _as_source(item: Any) -> _Source:
    """Accept a ``_Source`` or a ``(name, metadata)`` pair."""
    if isinstance(item, _Source):
        return item
    name, metadata = item
    known = (metadata.get("essential_metadata") or {}).get("mass_range") is not None
    return _Source(str(name), metadata, range_known=known)


class _Notes:
    """Where each setting came from, and what the plan costs a source."""

    def __init__(self) -> None:
        self.origins: Dict[str, str] = {}
        self.warnings: List[str] = []


def _positive(name: str, value: Optional[float]) -> None:
    if value is not None and not (math.isfinite(value) and value > 0):
        raise ConversionRefused(
            f"{name} must be a finite positive number, got {value}."
        )


def _check_overrides(
    min_mz: Optional[float],
    max_mz: Optional[float],
    width_at_mz: Optional[float],
    reference_mz: Optional[float],
    tof_a: Optional[float],
    tof_b: Optional[float],
    bins_per_fwhm: Optional[float],
) -> None:
    """Refuse a forced number that cannot be part of an axis."""
    _positive("min_mz", min_mz)
    _positive("max_mz", max_mz)
    _positive("width_at_mz", width_at_mz)
    _positive("reference_mz", reference_mz)
    _positive("bins_per_fwhm", bins_per_fwhm)
    if (tof_a is None) != (tof_b is None):
        raise ConversionRefused(
            "tof_a and tof_b are one pair: the width law needs both, or neither."
        )
    if tof_a is not None and tof_b is not None:
        if not (math.isfinite(tof_a) and math.isfinite(tof_b)):
            raise ConversionRefused("tof_a and tof_b must be finite numbers.")
        if tof_a < 0 or tof_b < 0 or (tof_a == 0 and tof_b == 0):
            raise ConversionRefused(
                "tof_a and tof_b must not be negative, and not both zero: a law "
                "with no width lays no bins."
            )
    if width_at_mz is not None and bins_per_fwhm is not None:
        raise ConversionRefused(
            "width_at_mz and bins_per_fwhm each fix the bin count of a tof axis; "
            "pass one of them."
        )


def _fmt_each(items: Sequence[Tuple[str, str]]) -> str:
    return "; ".join(f"{name}: {value}" for name, value in items)


# ----------------------------------------------------------------------
# The rules
# ----------------------------------------------------------------------


def _plan_axis_type(sources: Sequence[_Source], forced: Any, notes: _Notes) -> AxisType:
    """The forced axis type, else the one every source's detector agrees on."""
    chosen = _resolve_config_enum(forced, AXIS_TYPE_BY_NAME, "axis_type")
    if chosen is not None:
        notes.origins["axis_type"] = "forced"
        return chosen  # type: ignore[no-any-return]
    tree = ResamplingDecisionTree()
    detected = [(s.name, tree.select_axis_type(s.metadata)) for s in sources]
    kinds = {kind for _, kind in detected}
    if len(kinds) > 1:
        raise ConversionRefused(
            "The sources' detectors want different axis types ("
            + _fmt_each([(n, k.value) for n, k in detected])
            + "). Pass axis_type to choose one."
        )
    kind = detected[0][1]
    notes.origins["axis_type"] = "all sources agree"
    return kind


def _plan_range(
    sources: Sequence[_Source],
    min_mz: Optional[float],
    max_mz: Optional[float],
    notes: _Notes,
) -> Tuple[float, float]:
    """The forced bounds, else the range every source covers."""
    known = [s for s in sources if s.range_known]
    unknown = [s.name for s in sources if not s.range_known]
    if unknown and (min_mz is None or max_mz is None):
        raise ConversionRefused(
            "No mass range is stated by " + ", ".join(unknown) + ", and finding "
            "one would mean decoding spectra. Pass both min_mz and max_mz."
        )
    low = min_mz if min_mz is not None else max(s.mass_range[0] for s in known)
    high = max_mz if max_mz is not None else min(s.mass_range[1] for s in known)
    if not low < high:
        ranges = _fmt_each(
            [(s.name, f"[{s.mass_range[0]:g}, {s.mass_range[1]:g}]") for s in known]
        )
        if min_mz is not None and max_mz is not None:
            raise ConversionRefused(
                f"The range [{low:g}, {high:g}] m/z is empty: the minimum has "
                "to be below the maximum."
            )
        raise ConversionRefused(
            f"The sources have no m/z in common: the planned range "
            f"[{low:g}, {high:g}] is empty ({ranges}). Pass min_mz and max_mz "
            "to plan a range anyway."
        )
    notes.origins["min_mz"] = (
        "forced" if min_mz is not None else "intersection of the sources' ranges"
    )
    notes.origins["max_mz"] = (
        "forced" if max_mz is not None else "intersection of the sources' ranges"
    )
    _warn_range(known, low, high, min_mz is None, max_mz is None, notes)
    return float(low), float(high)


def _warn_range(
    known: Sequence[_Source],
    low: float,
    high: float,
    low_derived: bool,
    high_derived: bool,
    notes: _Notes,
) -> None:
    """Say which sources the planned range does not fit, and which cut others."""
    for s in known:
        lo, hi = s.mass_range
        if lo > low or hi < high:
            notes.warnings.append(
                f"{s.name}: its own range [{lo:g}, {hi:g}] does not cover the "
                f"planned range [{low:g}, {high:g}]; bins outside "
                f"[{lo:g}, {hi:g}] are always zero in this source's store."
            )
    # An imzML range is the span of the data observed, not the acquisition
    # window, so an imzML source that sets a bound may cut signal the other
    # sources hold near the edge. A PHI range is its full channel axis. An
    # mzPeak range is observed data only when the archive states no scan
    # window, which the metadata does not say cheaply, so it is not flagged.
    for s in known:
        if (s.source_format or "").lower() != "imzml":
            continue
        lo, hi = s.mass_range
        others = [o for o in known if o is not s]
        bounds = []
        if low_derived and lo == low and any(o.mass_range[0] < low for o in others):
            bounds.append("minimum")
        if high_derived and hi == high and any(o.mass_range[1] > high for o in others):
            bounds.append("maximum")
        if bounds:
            notes.warnings.append(
                f"{s.name}: the planned {' and '.join(bounds)} m/z comes from "
                "this imzML source's observed range, which is the span of its "
                "data and not its acquisition window; the intersection may cut "
                "signal the other sources hold near the edge. Pass min_mz and "
                "max_mz to widen it."
            )


def _plan_tof(
    sources: Sequence[_Source],
    tof_a: Optional[float],
    tof_b: Optional[float],
    bins_per_fwhm: Optional[float],
    width_at_mz: Optional[float],
    reference_mz: float,
    notes: _Notes,
) -> Tuple[float, float, float]:
    """``(A, B, bins_per_fwhm)`` of a ``tof`` axis."""
    if tof_a is not None and tof_b is not None:
        law = (float(tof_a), float(tof_b))
        notes.origins["tof_law"] = "forced"
    else:
        tree = ResamplingDecisionTree()
        declared = [(s.name, tree.select_tof_law(s.metadata)) for s in sources]
        declared = [(n, d) for n, d in declared if d is not None]
        if not declared:
            raise ConversionRefused(
                "A tof axis needs tof_a and tof_b: none of the sources "
                "declares a width law."
            )
        laws = {(float(d[0]), float(d[1])) for _, d in declared}
        if len(laws) > 1:
            raise ConversionRefused(
                "The sources declare different TOF width laws ("
                + _fmt_each([(n, f"A={d[0]:g} B={d[1]:g}") for n, d in declared])
                + "). Pass tof_a and tof_b to choose one."
            )
        law = laws.pop()
        notes.origins["tof_law"] = (
            "all sources that declare one agree ("
            + ", ".join(n for n, _ in declared)
            + ")"
        )
    generator = TOFAxisGenerator(*law)
    if bins_per_fwhm is not None:
        k = float(bins_per_fwhm)
        notes.origins["bins_per_fwhm"] = "forced"
    elif width_at_mz is not None:
        k = generator.bins_per_fwhm_for(reference_mz, float(width_at_mz))
        notes.origins["bins_per_fwhm"] = (
            f"derived from the forced width {width_at_mz:g} Da at m/z {reference_mz:g}"
        )
    else:
        k = DEFAULT_BINS_PER_FWHM
        notes.origins["bins_per_fwhm"] = "default"
    realised = float(generator.bin_width_at(reference_mz, k))
    notes.origins["width_at_mz"] = (
        f"{realised:.4g} Da at m/z {reference_mz:g}, from the width law "
        "(a result, not a setting of the plan)"
    )
    return law[0], law[1], float(k)


def _own_width_at(source: _Source, kind: AxisType, reference_mz: float) -> float:
    """The width this source's own conversion would use, stated at ``reference_mz``."""
    name = axis_name(kind)
    detected = ResamplingDecisionTree().select_reference_width(source.metadata)
    width, at = reference_params(ResamplingConfig(), name, detected_width=detected)
    scale: float = (reference_mz / float(at)) ** _WIDTH_EXPONENT[name]
    return float(width) * scale


def _source_method(
    source: _Source, forced: Optional[ResamplingMethod]
) -> Optional[ResamplingMethod]:
    """The method this source converts with: the caller's, else its detector's."""
    if forced is not None:
        return forced
    try:
        return ResamplingDecisionTree().select_strategy(source.metadata)
    except Exception as exc:  # advisory only: a source the chain cannot classify
        logger.debug("No detector method for %s: %s", source.name, exc)
        return None


def _plan_width(
    sources: Sequence[_Source],
    kind: AxisType,
    width_at_mz: Optional[float],
    reference_mz: float,
    forced_method: Optional[ResamplingMethod],
    notes: _Notes,
) -> float:
    """The forced width, else the widest width any source declares."""
    own = [(s, _own_width_at(s, kind, reference_mz)) for s in sources]
    if width_at_mz is not None:
        width = float(width_at_mz)
        notes.origins["width_at_mz"] = "forced"
    else:
        width = max(w for _, w in own)
        notes.origins["width_at_mz"] = (
            "widest declared width: "
            + _fmt_each([(s.name, f"{w:.4g} Da") for s, w in own])
            + f" (each stated at m/z {reference_mz:g})"
        )
    for s, w in own:
        if w >= width * (1.0 - _WIDTH_TOLERANCE):
            continue
        text = (
            f"{s.name}: its own bin width at m/z {reference_mz:g} is {w:.4g} Da, "
            f"narrower than the planned {width:.4g} Da; it loses resolution in "
            "this store."
        )
        if _source_method(s, forced_method) is ResamplingMethod.TIC_PRESERVING:
            text += (
                " Under tic_preserving, bins coarser than the source's own can "
                "merge narrow peaks (#442)."
            )
        notes.warnings.append(text)
    return width


def _warn_spectrum_types(sources: Sequence[_Source], notes: _Notes) -> None:
    kinds = {s.spectrum_type for s in sources if s.spectrum_type is not None}
    if len(kinds) > 1:
        notes.warnings.append(
            "The sources hold different spectrum types ("
            + _fmt_each([(s.name, s.spectrum_type or "not stated") for s in sources])
            + "); a centroid and a profile source on one axis are each placed "
            "by their own method."
        )


def _lay_axis(first: _Source, plan: MassAxisPlan) -> Tuple[int, str]:
    """Build the planned axis the way a conversion does, and fingerprint it.

    Every setting is explicit in the config, so the detected values are never
    consulted and any source gives the same axis. The first source's metadata
    stands in for all of them, with its range replaced by the planned one.
    """
    config = plan.to_resampling_config()
    config.pop("mass_axis_sha256")
    metadata = dict(first.metadata)
    essential = dict(metadata.get("essential_metadata") or {})
    essential["mass_range"] = (plan.min_mz, plan.max_mz)
    metadata["essential_metadata"] = essential
    planner = AxisPlanner.from_metadata(metadata, normalize_resampling_config(config))
    axis = planner.settle(_refuse_wide_axis).axis
    return int(len(axis)), mass_axis_sha256(axis)


def _refuse_wide_axis(n_bins: int) -> None:
    """Refuse an axis no conversion here could hold, before it is laid.

    The conversion's own memory check: a plan for an axis that every
    conversion would refuse is no use, and laying it first can take more
    memory than the machine has.
    """
    from .converters.spatialdata.csc_assembly import mass_axis_refusal

    refusal = mass_axis_refusal(
        int(n_bins),
        levers="pass a wider width_at_mz, or fewer bins_per_fwhm on a tof axis",
    )
    if refusal is not None:
        raise ConversionRefused(refusal)


def _plan_from_metadata(
    sources: Sequence[Any],
    *,
    axis_type: Any = None,
    min_mz: Optional[float] = None,
    max_mz: Optional[float] = None,
    width_at_mz: Optional[float] = None,
    reference_mz: Optional[float] = None,
    tof_a: Optional[float] = None,
    tof_b: Optional[float] = None,
    bins_per_fwhm: Optional[float] = None,
    method: Any = None,
) -> MassAxisPlan:
    """Plan one axis from what the detector chain is told about each source.

    Reads no file. ``sources`` holds ``_Source`` items, or ``(name,
    metadata)`` pairs whose range is taken as known.

    Raises:
        ConversionRefused: As :func:`plan_mass_axis` describes.
    """
    srcs = [_as_source(s) for s in sources]
    if not srcs:
        raise ConversionRefused("A mass axis plan needs at least one source.")
    _check_overrides(
        min_mz, max_mz, width_at_mz, reference_mz, tof_a, tof_b, bins_per_fwhm
    )
    forced_method = _resolve_config_enum(method, METHOD_BY_NAME, "method")

    notes = _Notes()
    kind = _plan_axis_type(srcs, axis_type, notes)
    low, high = _plan_range(srcs, min_mz, max_mz, notes)

    is_tof = kind is AxisType.TOF
    if not is_tof and (tof_a is not None or bins_per_fwhm is not None):
        raise ConversionRefused(
            f"tof_a, tof_b and bins_per_fwhm describe a tof axis, but the "
            f"planned axis type is {kind.value}."
        )
    if reference_mz is not None:
        ref = float(reference_mz)
        notes.origins["reference_mz"] = "forced"
    elif is_tof:
        ref = DEFAULT_REFERENCE_MZ
        notes.origins["reference_mz"] = "default"
    else:
        ref = float(reference_params(ResamplingConfig(), axis_name(kind))[1])
        notes.origins["reference_mz"] = f"default for a {kind.value} axis"

    fields: Dict[str, Any] = {}
    if is_tof:
        a, b, k = _plan_tof(srcs, tof_a, tof_b, bins_per_fwhm, width_at_mz, ref, notes)
        fields.update(tof_a=a, tof_b=b, bins_per_fwhm=k)
    else:
        fields["width_at_mz"] = _plan_width(
            srcs, kind, width_at_mz, ref, forced_method, notes
        )
    _warn_spectrum_types(srcs, notes)

    draft = MassAxisPlan(
        axis_type=kind,
        min_mz=low,
        max_mz=high,
        reference_mz=ref,
        n_bins=0,
        sha256="0" * 64,
        method=forced_method,
        **fields,
    )
    try:
        n_bins, digest = _lay_axis(srcs[0], draft)
    except ConversionRefused as exc:
        raise ConversionRefused(f"The planned axis cannot be built: {exc}") from exc

    plan = replace(
        draft,
        n_bins=n_bins,
        sha256=digest,
        origins=dict(notes.origins),
        warnings=tuple(notes.warnings),
        sources=tuple(s.name for s in srcs),
    )
    for setting, origin in plan.origins.items():
        logger.info("Mass axis plan, %s: %s", setting, origin)
    for warning in plan.warnings:
        logger.warning("Mass axis plan: %s", warning)
    return plan


# ----------------------------------------------------------------------
# Reading the sources
# ----------------------------------------------------------------------


def _source_names(paths: Sequence[Path]) -> List[str]:
    """File names, with the parent folder added where two would collide."""
    names = [p.name for p in paths]
    out = []
    for path, name in zip(paths, names):
        if names.count(name) > 1:
            name = f"{path.parent.name}/{name}"
        out.append(name)
    for i, name in enumerate(out):
        if out.count(name) > 1:
            out[i] = f"{name} ({i + 1})"
    return out


def _read_source(path: Path, name: str) -> _Source:
    """Read one source's header, never its spectra."""
    if not path.exists():
        raise ConversionRefused(f"{name}: the path does not exist: {path}")
    reader, error = _build_reader(path)
    if reader is None:
        raise ConversionRefused(f"{name}: could not be read: {error}")
    try:
        essential = reader.get_essential_metadata()
        comprehensive = None
        try:
            comprehensive = reader.get_comprehensive_metadata()
        except Exception as exc:
            logger.debug("Comprehensive metadata unavailable for %s: %s", name, exc)
        metadata = _resampling_metadata_dict(essential, comprehensive)
        known = bool(getattr(essential, "mass_range_known", True))
    except Exception as exc:
        raise ConversionRefused(f"{name}: could not be read: {exc}") from exc
    finally:
        _close_quietly(reader)
    try:
        source_format: Optional[str] = detect_format(path)
    except Exception:
        source_format = None
    if known and (source_format or "").lower() == "imzml":
        metadata = _with_outward_range(path, metadata)
    return _Source(name, metadata, range_known=known, source_format=source_format)


def _with_outward_range(path: Path, metadata: Dict[str, Any]) -> Dict[str, Any]:
    """An imzML source's range moved out by the precision of its header text.

    The header states the lowest and highest observed m/z as rounded text:
    ``100.123457`` for a first sample at ``100.12345678``. A conversion
    reads the exact range from the ``.ibd``, but a plan reads no spectra. A
    plan bound taken from the rounded text can sit inside the data, and a
    conversion given the plan then drops every sample beyond it (#473). So
    the header is read again for bounds moved out by one unit of the last
    digit each states. On a continuous file that read stops after the first
    spectra, as the metadata read does.
    """
    from .metadata.constants import BinaryDataType
    from .readers.imzml.header import scan_observed_mz_range

    mode = (metadata.get("format_specific") or {}).get("file_mode")
    outward = scan_observed_mz_range(
        path, shared_axis=mode == BinaryDataType.CONTINUOUS, outward=True
    )
    if outward is None:
        return metadata
    essential = dict(metadata.get("essential_metadata") or {})
    essential["mass_range"] = outward
    moved = dict(metadata)
    moved["essential_metadata"] = essential
    return moved


def plan_mass_axis(
    sources: Union[str, Path, Sequence[Union[str, Path]]],
    *,
    axis_type: Optional[Union[str, AxisType]] = None,
    min_mz: Optional[float] = None,
    max_mz: Optional[float] = None,
    width_at_mz: Optional[float] = None,
    reference_mz: Optional[float] = None,
    tof_a: Optional[float] = None,
    tof_b: Optional[float] = None,
    bins_per_fwhm: Optional[float] = None,
    method: Optional[Union[str, ResamplingMethod]] = None,
) -> MassAxisPlan:
    """Plan one mass axis for several sources.

    Each source is read the way :func:`thyra.preview_msi` reads it: its
    header only, no spectra. Every setting you pass is used as given. The
    rest are decided from the sources, and the plan records where each came
    from in ``origins`` and what it costs a source in ``warnings``.

    * Axis type: every source's detector agrees, else refused.
    * Range: the intersection of the sources' ranges.
    * Width (non-``tof`` axes): the widest width any source declares, each
      first stated at ``reference_mz`` under the axis law, so that on the
      shared axis no source's peaks scatter over bins finer than they are.
      The narrower sources lose resolution and are named in the warnings.
    * ``tof`` axes: the width law every declaring source agrees on.

    The plan fixes the axis only. Each conversion keeps its own detector's
    resampling method unless you force ``method``. Pass the plan as
    ``resampling_config`` to :func:`thyra.convert_msi`; the conversion
    refuses when the axis it builds hashes differently from the plan's.

    Args:
        sources: Paths of the sources, or one path.
        axis_type: Force the axis type (a name or an ``AxisType``).
        min_mz: Force the lower bound of the range.
        max_mz: Force the upper bound of the range.
        width_at_mz: Force the bin width in Da at ``reference_mz``.
        reference_mz: The m/z the width is stated at. Default: 300 for a
            ``linear_tof`` axis, 1000 otherwise.
        tof_a: Force ``A`` of the TOF width law (needs ``tof_b``).
        tof_b: Force ``B`` of the TOF width law.
        bins_per_fwhm: Force the bins per peak width of a ``tof`` axis.
        method: Force one resampling method for every conversion.

    Returns:
        The plan.

    Raises:
        ConversionRefused: If a source cannot be read or states no m/z range
            (and ``min_mz`` and ``max_mz`` were not both given), if the ranges
            have no m/z in common, if the detectors want different axis types
            or TOF laws and none was forced, or if the axis cannot be built.
    """
    if isinstance(sources, (str, Path)):
        sources = [sources]
    paths = [Path(p) for p in sources]
    names = _source_names(paths)
    read = [_read_source(p, n) for p, n in zip(paths, names)]
    return _plan_from_metadata(
        read,
        axis_type=axis_type,
        min_mz=min_mz,
        max_mz=max_mz,
        width_at_mz=width_at_mz,
        reference_mz=reference_mz,
        tof_a=tof_a,
        tof_b=tof_b,
        bins_per_fwhm=bins_per_fwhm,
        method=method,
    )
