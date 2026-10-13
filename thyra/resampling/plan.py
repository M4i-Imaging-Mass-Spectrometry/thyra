# thyra/resampling/plan.py

"""A mass axis fixed once, for several conversions (issue #449, D43).

The axis a conversion builds depends on five settings and nothing else: the
range, the axis type, the bin count and, for a ``tof`` axis, the width law
(:class:`~thyra.resampling.axis_planner.AxisPlan`). The resampling method
does not enter it. A :class:`MassAxisPlan` therefore holds those settings
and the fingerprint of the axis they lay (D42), and leaves each source's
method to its own detector unless the caller forced one.

This module is pure data: it imports no reader and no converter, so a plan
can be saved, loaded and passed around without either. Deciding the
settings from several sources is :func:`thyra.mass_axis_plan.plan_mass_axis`.
"""

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

from ..errors import ConversionRefused
from .types import AxisType, ResamplingMethod

#: The ``format`` key of a saved plan, and the layout version it is in.
PLAN_FORMAT = "thyra-mass-axis-plan"
PLAN_VERSION = 1

_SHA256 = re.compile(r"[0-9a-f]{64}")


def check_sha256(value: Any, where: str) -> str:
    """Return ``value`` when it is 64 lower-case hex characters.

    Args:
        value: What the caller supplied.
        where: The key it came from, named in the refusal.

    Raises:
        ConversionRefused: If ``value`` is anything else.
    """
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ConversionRefused(
            f"{where} must be a SHA-256 as 64 lower-case hexadecimal "
            f"characters, got {value!r}."
        )
    return value


def _positive_number(
    data: Dict[str, Any], key: str, *, allow_zero: bool = False
) -> float:
    """``data[key]`` as a finite positive float, or a refusal naming the key.

    ``allow_zero`` is for the two TOF law terms, either of which may be 0.

    A plan file is edited by hand as often as it is written by
    :meth:`MassAxisPlan.save`, so a string, a list or a zero is said as a
    refusal rather than raised as a ``TypeError`` from ``float()``.
    """
    value = data.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConversionRefused(
            f"The mass axis plan's {key!r} must be a number, not {value!r}."
        )
    number = float(value)
    if not (math.isfinite(number) and (number >= 0 if allow_zero else number > 0)):
        kind = "non-negative" if allow_zero else "positive"
        raise ConversionRefused(
            f"The mass axis plan's {key!r} must be a finite {kind} number, "
            f"not {value!r}."
        )
    return number


def _optional_positive_number(
    data: Dict[str, Any], key: str, *, allow_zero: bool = False
) -> Optional[float]:
    if data.get(key) is None:
        return None
    return _positive_number(data, key, allow_zero=allow_zero)


def _strings(data: Dict[str, Any], key: str) -> Tuple[str, ...]:
    """``data[key]`` as a tuple of strings; absent means empty."""
    value = data.get(key) or ()
    if not isinstance(value, (list, tuple)) or not all(
        isinstance(item, str) for item in value
    ):
        raise ConversionRefused(
            f"The mass axis plan's {key!r} must be a list of strings."
        )
    return tuple(value)


@dataclass(frozen=True)
class MassAxisPlan:
    """One mass axis, fixed for several sources.

    Attributes:
        axis_type: The spacing law of the bins.
        min_mz: Lower end of the range the bins span.
        max_mz: Upper end of it.
        reference_mz: The m/z at which the bin width is stated.
        n_bins: How many bins the axis has.
        sha256: The fingerprint of the axis (D42). A conversion given
            this plan refuses when the axis it builds hashes differently.
        width_at_mz: Bin width in Da at ``reference_mz``. ``None`` on a
            ``tof`` axis, whose width follows from its law.
        tof_a: ``A`` of the TOF width law (``tof`` axis only).
        tof_b: ``B`` of the TOF width law (``tof`` axis only).
        bins_per_fwhm: Bins per peak width (``tof`` axis only).
        method: The resampling method, only when the caller forced one.
            ``None`` leaves each source to its own detector.
        origins: One plain line per setting saying where the value came
            from.
        warnings: What the plan costs a source, each naming its source.
        sources: The names of the sources the plan was made from.
    """

    axis_type: AxisType
    min_mz: float
    max_mz: float
    reference_mz: float
    n_bins: int
    sha256: str
    width_at_mz: Optional[float] = None
    tof_a: Optional[float] = None
    tof_b: Optional[float] = None
    bins_per_fwhm: Optional[float] = None
    method: Optional[ResamplingMethod] = None
    origins: Dict[str, str] = field(default_factory=dict)
    warnings: Tuple[str, ...] = ()
    sources: Tuple[str, ...] = ()

    def to_resampling_config(self) -> Dict[str, Any]:
        """The explicit request that rebuilds this axis, as a config dict.

        Never holds both ``target_bins`` and ``width_at_mz``: the bin count
        follows from the settings, and a conversion refuses the pair. The
        fingerprint travels as ``mass_axis_sha256``.
        """
        config: Dict[str, Any] = {
            "axis_type": self.axis_type.value,
            "min_mz": self.min_mz,
            "max_mz": self.max_mz,
            "reference_mz": self.reference_mz,
        }
        if self.axis_type is AxisType.TOF:
            config["tof_a"] = self.tof_a
            config["tof_b"] = self.tof_b
            config["bins_per_fwhm"] = self.bins_per_fwhm
        else:
            config["width_at_mz"] = self.width_at_mz
        if self.method is not None:
            config["method"] = self.method.value
        config["mass_axis_sha256"] = self.sha256
        return config

    def to_dict(self) -> Dict[str, Any]:
        """A JSON-safe dict: enums as their values, tuples as lists."""
        return {
            "format": PLAN_FORMAT,
            "version": PLAN_VERSION,
            "axis_type": self.axis_type.value,
            "min_mz": self.min_mz,
            "max_mz": self.max_mz,
            "reference_mz": self.reference_mz,
            "width_at_mz": self.width_at_mz,
            "tof_a": self.tof_a,
            "tof_b": self.tof_b,
            "bins_per_fwhm": self.bins_per_fwhm,
            "method": None if self.method is None else self.method.value,
            "n_bins": self.n_bins,
            "sha256": self.sha256,
            "origins": dict(self.origins),
            "warnings": list(self.warnings),
            "sources": list(self.sources),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MassAxisPlan":
        """Rebuild a plan from :meth:`to_dict` output.

        Raises:
            ConversionRefused: If ``data`` is not a plan of this layout, or
                lacks a setting every plan holds.
        """
        if not isinstance(data, dict) or data.get("format") != PLAN_FORMAT:
            raise ConversionRefused(
                f"Not a mass axis plan: the 'format' key must be {PLAN_FORMAT!r}."
            )
        if data.get("version") != PLAN_VERSION:
            raise ConversionRefused(
                f"Unsupported mass axis plan version {data.get('version')!r}; "
                f"this Thyra reads version {PLAN_VERSION}."
            )
        required = ("axis_type", "min_mz", "max_mz", "reference_mz", "n_bins", "sha256")
        missing = [key for key in required if data.get(key) is None]
        if missing:
            raise ConversionRefused(
                "The mass axis plan lacks " + ", ".join(repr(k) for k in missing) + "."
            )
        try:
            axis_type = AxisType(data["axis_type"])
            raw_method = data.get("method")
            method = None if raw_method is None else ResamplingMethod(raw_method)
        except (ValueError, TypeError) as exc:
            raise ConversionRefused(f"The mass axis plan is not valid: {exc}.") from exc
        if method is ResamplingMethod.NONE:
            raise ConversionRefused(
                "The mass axis plan's 'method' cannot be 'none': a plan resamples."
            )
        n_bins = data["n_bins"]
        if isinstance(n_bins, bool) or not isinstance(n_bins, int) or n_bins < 2:
            raise ConversionRefused(
                f"The mass axis plan's 'n_bins' must be a whole number of at "
                f"least 2, not {n_bins!r}."
            )
        origins = data.get("origins") or {}
        if not isinstance(origins, dict) or not all(
            isinstance(k, str) and isinstance(v, str) for k, v in origins.items()
        ):
            raise ConversionRefused(
                "The mass axis plan's 'origins' must map setting names to text."
            )
        return cls(
            axis_type=axis_type,
            min_mz=_positive_number(data, "min_mz"),
            max_mz=_positive_number(data, "max_mz"),
            reference_mz=_positive_number(data, "reference_mz"),
            width_at_mz=_optional_positive_number(data, "width_at_mz"),
            tof_a=_optional_positive_number(data, "tof_a", allow_zero=True),
            tof_b=_optional_positive_number(data, "tof_b", allow_zero=True),
            bins_per_fwhm=_optional_positive_number(data, "bins_per_fwhm"),
            method=method,
            n_bins=n_bins,
            sha256=check_sha256(data["sha256"], "The plan's 'sha256'"),
            origins=dict(origins),
            warnings=_strings(data, "warnings"),
            sources=_strings(data, "sources"),
        )

    def save(self, path: Union[str, Path]) -> None:
        """Write the plan as UTF-8 JSON."""
        text = json.dumps(self.to_dict(), indent=2, ensure_ascii=False)
        Path(path).write_text(text + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path: Union[str, Path]) -> "MassAxisPlan":
        """Read a plan written by :meth:`save`.

        Raises:
            ConversionRefused: If the file cannot be read, is not JSON, or
                is not a plan.
        """
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise ConversionRefused(
                f"Could not read the mass axis plan {Path(path).name}: {exc}"
            ) from exc
        return cls.from_dict(data)
