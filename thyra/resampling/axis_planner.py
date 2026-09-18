# thyra/resampling/axis_planner.py

"""What decides the common mass axis of one conversion.

Everything here answers one question -- *which* m/z values the store's
``var`` will carry -- and nothing here touches a spectrum. The operator
that places peaks onto that axis is the other half of the same concern and
lives in :mod:`thyra.resampling.strategies`; the two meet at
:func:`~thyra.resampling.strategies.build_strategy`, whose middle three
arguments are exactly the triple :meth:`AxisPlanner.settle` hands back.

The decision is where the design decisions are (D1-D7, D13-D16, D20, D21):
which instrument the detector chain reads out of the source's metadata,
which law its bins follow, how wide one bin is at a reference m/z, and how
many of them span the range. It used to be a dozen converter methods over
sixteen attributes assigned in three places, so the only way to ask what a
given source and a given ``--resample-*`` pair come to was to stand up a
converter with a reader and a writable output path. A planner takes a
reader and a :class:`~thyra.resampling.types.ResamplingConfig` -- or, with
no reader at all, the metadata dict the detectors read and that same config
-- and answers it (issue #352).

A planner exists for every conversion, including one that does not resample:
``config=None`` means the axis is the source's own, which is as much a
decision about the axis as any generator's, and keeping both routes here is
what holds the axis and its linearisation together. The two must never
disagree about which axis they describe, and a raw union axis was laid by no
generator, so it has none.
"""

import logging
from dataclasses import dataclass
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Tuple,
    Union,
    cast,
)

import numpy as np
from numpy.typing import NDArray

from ..errors import MALFORMED_METADATA, ConversionRefused
from .binning import usable_linearisation
from .common_axis import CommonAxisBuilder
from .decision_tree import ResamplingDecisionTree
from .mass_axis.tof_generator import DEFAULT_BINS_PER_FWHM, TOFAxisGenerator
from .types import (
    DEFAULT_REFERENCE_MZ,
    AxisLinearisation,
    AxisType,
    ResamplingConfig,
    ResamplingMethod,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..core.base_reader import BaseMSIReader
    from ..metadata.types import ComprehensiveMetadata, EssentialMetadata

logger = logging.getLogger(__name__)

#: Axis entries differenced at once by :func:`bin_width_range`. 32 MB of
#: float64 at a time; read at call time so a test can shrink it.
BIN_WIDTH_CHUNK = 1 << 22


def bin_width_range(axis: NDArray[np.float64]) -> Tuple[float, float]:
    """The narrowest and widest gap between consecutive axis entries.

    Chunked, because the obvious ``np.diff(axis)`` allocates a second
    array the length of the axis to produce two numbers for a log line --
    1.6 GB of it on the 200M-bin axis issue #251 measures. Chunks overlap
    by one entry so no gap falls between two of them.
    """
    axis = np.asarray(axis, dtype=np.float64)
    if axis.size < 2:
        return 0.0, 0.0
    low = np.inf
    high = -np.inf
    for start in range(0, axis.size - 1, BIN_WIDTH_CHUNK):
        stop = min(start + BIN_WIDTH_CHUNK + 1, axis.size)
        widths = np.diff(axis[start:stop])
        low = min(low, float(widths.min()))
        high = max(high, float(widths.max()))
    return low, high


def _resolve_config_enum(raw: Any, by_name: Dict[str, Any], key: str) -> Any:
    """Resolve one resampling-config value to its enum member.

    ``None``, ``"auto"`` and ``""`` all mean "decide this automatically"
    and resolve to ``None``.

    The CLI is protected by ``click.Choice``, but the Python API takes
    whatever the caller passes. Reject anything unrecognised rather than
    dropping it: silently treating ``"tic_preserving "`` or a typo as
    "auto-detect" hands back a conversion that ignored the request
    without saying so.

    Args:
        raw: The value as supplied by the caller.
        by_name: Accepted string spellings mapped to enum members.
        key: The config key, used in error messages.

    Raises:
        ConversionRefused: If ``raw`` names no accepted value, or is
            neither a string nor one of the accepted enum members.
    """
    if raw is None:
        return None

    valid = ", ".join(repr(name) for name in ["auto", *by_name])

    if isinstance(raw, str):
        if raw in ("auto", ""):
            return None
        if raw not in by_name:
            raise ConversionRefused(
                f"Unknown resampling_config[{key!r}] value {raw!r}. "
                f"Valid values are: {valid}."
            )
        return by_name[raw]

    if raw in by_name.values():
        return raw

    raise ConversionRefused(
        f"Unsupported resampling_config[{key!r}] value {raw!r}. "
        f"Pass one of {valid}, or the matching enum member."
    )


#: Every key :func:`normalize_resampling_config` reads out of a
#: ``resampling_config`` dict. Anything else in the dict is a typo or a
#: leftover, and is warned about rather than ignored.
_RESAMPLING_CONFIG_KEYS = frozenset(
    {
        "method",
        "axis_type",
        "target_bins",
        "width_at_mz",
        "reference_mz",
        "min_mz",
        "max_mz",
        "gap_tolerance_da",
        "tof_a",
        "tof_b",
        "bins_per_fwhm",
    }
)


def normalize_resampling_config(
    config: Union[Dict[str, Any], ResamplingConfig],
) -> ResamplingConfig:
    """Normalise a resampling config dict or dataclass to a ResamplingConfig.

    Accepts either a plain dict (as produced by _build_resampling_config in
    __main__.py) or an already-constructed ResamplingConfig dataclass and
    returns a ResamplingConfig in both cases.

    Raises:
        ConversionRefused: If ``method`` or ``axis_type`` is not a
            recognised value.
    """
    if isinstance(config, ResamplingConfig):
        return config

    # A key this function does not read changes nothing, and used to
    # change nothing silently: ``{"target_bin": 4000}`` was accepted and
    # the axis built from the default (issue #250). Said once, with the
    # keys that would have worked, rather than refused -- a caller
    # passing an extra key is not necessarily wrong, but a caller with a
    # typo always is.
    unknown = sorted(set(config) - _RESAMPLING_CONFIG_KEYS)
    if unknown:
        logger.warning(
            "resampling_config holds %s, which %s not read and %s no effect. "
            "Known keys: %s.",
            ", ".join(repr(key) for key in unknown),
            "is" if len(unknown) == 1 else "are",
            "has" if len(unknown) == 1 else "have",
            ", ".join(sorted(_RESAMPLING_CONFIG_KEYS)),
        )

    # Only the methods the resampling pipeline actually implements, and
    # only the axis types CommonAxisBuilder has a generator for. These
    # deliberately match the CLI's click.Choice lists.
    method_by_name = {
        "nearest_neighbor": ResamplingMethod.NEAREST_NEIGHBOR,
        "tic_preserving": ResamplingMethod.TIC_PRESERVING,
    }
    axis_type_by_name = {
        "constant": AxisType.CONSTANT,
        "linear_tof": AxisType.LINEAR_TOF,
        "reflector_tof": AxisType.REFLECTOR_TOF,
        "tof": AxisType.TOF,
        "orbitrap": AxisType.ORBITRAP,
        "fticr": AxisType.FTICR,
    }

    reference_mz = config.get("reference_mz")

    def _optional_float(key: str) -> Optional[float]:
        value = config.get(key)
        return None if value is None else float(value)

    return ResamplingConfig(
        method=_resolve_config_enum(config.get("method"), method_by_name, "method"),
        axis_type=_resolve_config_enum(
            config.get("axis_type"), axis_type_by_name, "axis_type"
        ),
        target_bins=config.get("target_bins"),
        mass_width_da=config.get("width_at_mz"),
        reference_mz=(
            DEFAULT_REFERENCE_MZ if reference_mz is None else float(reference_mz)
        ),
        min_mz=config.get("min_mz"),
        max_mz=config.get("max_mz"),
        gap_tolerance_da=config.get("gap_tolerance_da"),
        tof_a=_optional_float("tof_a"),
        tof_b=_optional_float("tof_b"),
        bins_per_fwhm=_optional_float("bins_per_fwhm"),
    )


def axis_name(axis_type: Any) -> str:
    """The lower-case name of an axis type, enum member or not.

    ``AxisType`` members carry it as ``.value``; the branch below is for
    the stubs and the strings the Python API accepts, which reach the
    width law the same way.
    """
    if hasattr(axis_type, "value"):
        return str(axis_type.value)
    return str(axis_type).split(".")[-1].lower()


def tof_plan(
    config: ResamplingConfig,
    *,
    detected_tof_law: Optional[Tuple[float, float]] = None,
) -> Tuple[float, float, float]:
    """``(A, B, bins_per_fwhm)`` for an ``AxisType.TOF`` axis.

    The law is the caller's ``tof_a``/``tof_b`` pair, else the pair the
    detected instrument declared. ``bins_per_fwhm`` is derived from
    ``--resample-width-at-mz`` at the reference m/z when that was given, so
    the width flag means the same thing on every axis type; otherwise it is
    the API's ``bins_per_fwhm``, else 3.

    Args:
        config: The caller's resampling request.
        detected_tof_law: The ``(A, B)`` the detected instrument declares,
            when the caller named none.

    Returns:
        ``(A, B, bins_per_fwhm)``.

    Raises:
        ConversionRefused: If neither the caller nor the instrument has a
            law, which leaves the bin width undefined everywhere.
    """
    a = config.tof_a
    b = config.tof_b
    if a is None or b is None:
        if detected_tof_law is None:
            raise ConversionRefused(
                "A 'tof' mass axis needs the width law's coefficients: pass "
                "--tof-law A B, or convert a run whose instrument declares "
                "them (SELECT SERIES MRT centroid, timsTOF, PHI nanoTOF)."
            )
        a, b = detected_tof_law
    generator = TOFAxisGenerator(float(a), float(b))
    if config.mass_width_da is not None:
        k: Optional[float] = generator.bins_per_fwhm_for(
            float(config.reference_mz), float(config.mass_width_da)
        )
    else:
        k = config.bins_per_fwhm
        if k is None:
            k = DEFAULT_BINS_PER_FWHM
    return float(a), float(b), float(k)


def reference_params(
    config: ResamplingConfig,
    name: str,
    *,
    detected_width: Optional[Tuple[float, float]] = None,
    detected_tof_law: Optional[Tuple[float, float]] = None,
) -> Tuple[float, float]:
    """``(width_da, reference_mz)`` for the axis about to be built.

    Precedence: the caller's ``--resample-width-at-mz`` /
    ``--resample-reference-mz``; then the width the detected instrument
    declared (auto path only); then the per-axis-type default -- 17 mDa at
    m/z 300 for ``linear_tof``, chosen to be close to the axis SCiLS Lab
    produces for FlexImaging data, and 5 mDa at m/z 1000 for everything
    else.

    One function, because the bin count and the generator have to be given
    the same width or the axis is valid and the wrong width.

    Args:
        config: The caller's resampling request.
        name: The axis type's lower-case name, as :func:`axis_name` gives it.
        detected_width: ``(width_da, reference_mz)`` the detected
            instrument asked for, or ``None``.
        detected_tof_law: The detected ``(A, B)``, read only by the ``tof``
            branch.

    Returns:
        ``(width_da, reference_mz)``.
    """
    if config.mass_width_da is not None:
        return config.mass_width_da, config.reference_mz
    if name == "tof":
        # The law and the bins-per-FWHM fix the width everywhere; report
        # the one realised at the reference m/z.
        a, b, k = tof_plan(config, detected_tof_law=detected_tof_law)
        ref_mz = float(config.reference_mz)
        return float(TOFAxisGenerator(a, b).bin_width_at(ref_mz, k)), ref_mz
    if detected_width is not None:
        return float(detected_width[0]), float(detected_width[1])
    if name == "linear_tof":
        return 0.017, 300.0
    return 0.005, 1000.0


def bin_count_for_width(
    config: ResamplingConfig,
    min_mz: float,
    max_mz: float,
    axis_type: Any,
    *,
    detected_width: Optional[Tuple[float, float]] = None,
    detected_tof_law: Optional[Tuple[float, float]] = None,
) -> int:
    """How many bins of the requested width span a mass range.

    The count and :meth:`CommonAxisBuilder.build_physics_axis` have to
    apply the same spacing law: an axis whose count was derived under a
    different law is still a valid axis, but the width it realises at
    ``reference_mz`` is not the width that was asked for.

    Args:
        config: The caller's resampling request.
        min_mz: Lower end of the mass range.
        max_mz: Upper end of the mass range.
        axis_type: The axis type, which decides the spacing law.
        detected_width: What the detected instrument asked for, if
            anything; see :func:`reference_params`.
        detected_tof_law: The detected ``(A, B)``, for a ``tof`` axis.

    Returns:
        The bin count, never below 100.
    """
    name = axis_name(axis_type)
    width_at_mz, ref_mz = reference_params(
        config,
        name,
        detected_width=detected_width,
        detected_tof_law=detected_tof_law,
    )

    logger.info(
        f"Calculating bins for {width_at_mz*1000:.1f} mDa width at m/z {ref_mz:.1f}"
    )

    if name == "reflector_tof":
        # REFLECTOR_TOF: constant relative resolution (width ∝ m/z)
        # relative_resolution = reference_mz / width_at_mz
        # For logarithmic spacing: bins ≈ ln(max_mz/min_mz) * (reference_mz / width_at_mz)
        relative_resolution = ref_mz / width_at_mz
        bins = int(np.log(max_mz / min_mz) * relative_resolution)

    elif name == "tof":
        # TOF: bin width = sqrt(A m + B m^2) / k. The generator carries
        # the closed-form integral of 1 / width, so the count is exact.
        a, b, k = tof_plan(config, detected_tof_law=detected_tof_law)
        bins = TOFAxisGenerator(a, b).bin_count(min_mz, max_mz, k)
        logger.info(
            "TOF width law A=%.4g mDa^2/Da, B=%.4g, %.2f bins per FWHM",
            a,
            b,
            k,
        )

    elif name == "linear_tof":
        # LINEAR_TOF: bin width = k * sqrt(m/z), where k = width_at_mz / sqrt(reference_mz)
        # Number of bins: n = (2/k) * (sqrt(max_mz) - sqrt(min_mz))
        # This matches SCiLS Lab's "Linear TOF" mass axis calculation
        k = width_at_mz / np.sqrt(ref_mz)
        bins = int((2.0 / k) * (np.sqrt(max_mz) - np.sqrt(min_mz)))

    elif name == "orbitrap":
        # ORBITRAP: bin_width = k * (m/z)^1.5 with k = width_at_mz /
        # reference_mz^1.5. Integrating dm / w(m) over the range gives
        #   bins = 2 * (1/sqrt(min_mz) - 1/sqrt(max_mz)) *
        #          (reference_mz^1.5 / width_at_mz)
        # The factor 2 comes from d(m^-0.5)/dm = -1/2 * m^-1.5 and was
        # missing, which halved the bin count and made every bin twice
        # the requested width.
        scaling_factor = (ref_mz**1.5) / width_at_mz
        bins = int(2 * (1 / np.sqrt(min_mz) - 1 / np.sqrt(max_mz)) * scaling_factor)

    elif name == "fticr":
        # FTICR: width ∝ m/z^2, i.e. bin_width = k * (m/z)^2 with
        # k = width_at_mz / reference_mz^2. FTICRAxisGenerator lays the
        # axis out uniformly in 1/mz, so the bin count is the 1/mz span
        # divided by the step k:
        #   bins = (1/min_mz - 1/max_mz) * (reference_mz^2 / width_at_mz)
        # Without this branch an FT-ICR axis took its bin count from the
        # uniform formula below, so the realized width at reference_mz was
        # not the width that was asked for.
        scaling_factor = (ref_mz**2) / width_at_mz
        bins = int((1 / min_mz - 1 / max_mz) * scaling_factor)

    else:
        # LINEAR/CONSTANT: uniform spacing
        # bins = (max_mz - min_mz) / width_at_mz
        bins = int((max_mz - min_mz) / width_at_mz)

    # Ensure minimum bin count
    bins = max(100, bins)

    logger.info(f"Calculated {bins} bins for {name} axis type")
    return bins


@dataclass(frozen=True)
class AxisPlan:
    """What a resampled conversion settled on, before any axis is laid.

    Everything expensive follows from these four numbers, which is why
    they are resolved on their own: the per-bin cost of a conversion is
    known from ``target_bins`` alone, and a bin count that cannot fit is
    refused before the axis it describes is allocated (issue #251).

    Attributes:
        min_mz: Lower end of the range the bins span.
        max_mz: Upper end of it.
        axis_type: The law the bins follow -- an
            :class:`~thyra.resampling.types.AxisType` on every path the CLI
            and the Python API reach, and read through :func:`axis_name`
            rather than by identity so a caller's own spelling still lands
            on the right law.
        target_bins: How many of them.
    """

    min_mz: float
    max_mz: float
    axis_type: Any
    target_bins: int


@dataclass(frozen=True)
class SettledAxis:
    """The common mass axis of one conversion, and what is true about it.

    Attributes:
        axis: The axis itself, ascending float64.
        axis_range: The ``(min_mz, max_mz)`` the bins were laid across,
            or ``None`` on a raw axis. It, not the axis's own span,
            decides whether a peak is in range: a physics axis reports bin
            *centres* and so stops half a bin short at either end (issue
            #239).
        linearisation: The coordinate the axis may be indexed through, or
            ``None`` when every placement onto it has to be a search
            (D21). A raw union axis was laid by no generator and has none.
        method: The resampling method the conversion runs, after the
            detector chain and the TIC-preserving gate have both had their
            say; ``None`` when nothing is resampled.
        provenance: The resolved plan as the processing block records it,
            or ``None`` on a raw axis. With "auto" settings the requested
            config says nothing about the method, axis and bin width that
            were actually used.
    """

    axis: NDArray[np.float64]
    axis_range: Optional[Tuple[float, float]]
    linearisation: Optional[AxisLinearisation]
    method: Optional[ResamplingMethod]
    provenance: Optional[Dict[str, Any]]


class AxisPlanner:
    """The common mass axis of one conversion, decided from arguments.

    Two entry points. :meth:`resolve` settles the four numbers a resampled
    axis follows from, and :meth:`settle` lays the axis itself -- or takes
    the source's own when the conversion does not resample. Everything
    else is how those two reach their answer.
    """

    def __init__(
        self,
        reader: Optional["BaseMSIReader"] = None,
        config: Optional[ResamplingConfig] = None,
        *,
        detection_metadata: Optional[Dict[str, Any]] = None,
        essential_metadata: Optional["EssentialMetadata"] = None,
    ) -> None:
        """Build the axis planner of one conversion.

        The method is chosen here, from the detector chain or from the
        caller, because an explicit ``--resample-method`` has to be
        checked against the detector while there is still someone to warn
        (issue #246) -- and because the caches the extraction fills are
        the same ones the plan reads later.

        Args:
            reader: The MSI source. Read for its essential, comprehensive
                and spectrum metadata, and -- when nothing is resampled --
                for its own mass axis. ``None`` is allowed, and then
                ``detection_metadata`` is all the planner knows.
            config: The caller's resampling request, already normalised by
                :func:`normalize_resampling_config`. ``None`` means the
                conversion does not resample.
            detection_metadata: The dict the detector chain reads, for a
                caller that already has one. Supplying it keeps the
                planner off the reader entirely.
            essential_metadata: The source's essential metadata when it
                has already been read, so a second query is not made for
                the mass range.
        """
        self.reader = reader
        self.config = config
        # What the detector chain settled, and whether it was asked at all.
        # The TIC-preserving gate governs auto-selection only -- an explicit
        # method is taken as given and merely warned about (D15) -- so
        # :meth:`resolve` needs to know which this was before it re-asks the
        # gate against the axis it lands on.
        self._method: Optional[ResamplingMethod] = None
        self._method_was_auto = False
        # Filled by :meth:`resolve` when the detected instrument declares a
        # bin width or a TOF width law and the caller named neither.
        self._detected_reference_width: Optional[Tuple[float, float]] = None
        self._detected_tof_law: Optional[Tuple[float, float]] = None
        self._plan: Optional[AxisPlan] = None
        # Metadata caches. Every extractor swallows its failures at DEBUG,
        # so an extraction run before these exist reports nothing and hands
        # the decision tree an empty dict -- which is how *every* reader's
        # method was once chosen by DefaultDetector while the axis type,
        # resolved later, came from the right detector.
        self._essential_metadata = essential_metadata
        self._comprehensive_metadata: Optional["ComprehensiveMetadata"] = None
        self._spectrum_metadata: Optional[Dict[str, Any]] = None
        self._detection_metadata = detection_metadata

        if self.config is not None:
            self._select_method()

    @classmethod
    def from_metadata(
        cls,
        metadata: Dict[str, Any],
        config: Optional[ResamplingConfig] = None,
    ) -> "AxisPlanner":
        """A planner that has already been told everything about the source.

        The dict is the one the detector chain reads, as
        :meth:`detection_metadata` builds it; its ``essential_metadata``
        entry carries the mass range, so no reader is needed to plan an
        axis or to ask what the defaults come to for a given source.

        Args:
            metadata: What the detectors are told about the source.
            config: The caller's resampling request.

        Returns:
            A planner with no reader behind it.
        """
        return cls(None, config, detection_metadata=metadata)

    @property
    def method(self) -> Optional[ResamplingMethod]:
        """The resampling method this conversion runs, or ``None``."""
        return self._method

    @property
    def method_was_auto(self) -> bool:
        """Whether the detector chain chose the method rather than the caller."""
        return self._method_was_auto

    # ------------------------------------------------------------------
    # What the detectors are told
    # ------------------------------------------------------------------

    def adopt_essential_metadata(self, essential: "EssentialMetadata") -> None:
        """Take the source's essential metadata from whoever already read it.

        The conversion reads it once at initialization; handing it over
        here is what stops the plan asking the reader for it again, which
        on a Bruker source is a query rather than an attribute.
        """
        self._essential_metadata = essential

    def detection_metadata(self) -> Dict[str, Any]:
        """Everything the detector chain is told about the source, cached.

        One method, deliberately: there used to be two -- one that always
        re-extracted and one that returned the cache -- and they returned
        the same dict, because the extraction reads the same cached
        metadata objects either way.
        """
        if self._detection_metadata is not None:
            return self._detection_metadata

        metadata: Dict[str, Any] = {}
        self._extract_essential_metadata(metadata)
        self._extract_comprehensive_metadata(metadata)
        self._extract_spectrum_metadata(metadata)

        self._detection_metadata = metadata
        return metadata

    def _extract_essential_metadata(self, metadata: Dict[str, Any]) -> None:
        """Extract essential metadata for resampling decisions.

        The reader call keeps a broad catch: it is the boundary, and what
        surfaces there was raised by a vendor SDK, an XML parser or sqlite.
        Reading the result is narrowed, because that is Thyra's own code
        against an object Thyra built (issue #280).
        """
        try:
            essential = self._read_essential_metadata()
        except Exception as e:
            logger.debug(f"Could not read essential metadata: {e}")
            return
        if essential is None:
            return

        try:
            if hasattr(essential, "source_path"):
                metadata["source_path"] = str(essential.source_path)

            # Add essential metadata for resampling decisions
            metadata["essential_metadata"] = {
                "spectrum_type": getattr(essential, "spectrum_type", None),
                "dimensions": essential.dimensions,
                "mass_range": essential.mass_range,
                "source_path": str(essential.source_path),
                "total_peaks": getattr(essential, "total_peaks", None),
                "n_spectra": getattr(essential, "n_spectra", None),
            }
        except MALFORMED_METADATA as e:
            logger.debug(f"Could not extract essential metadata: {e}")

    def _extract_comprehensive_metadata(self, metadata: Dict[str, Any]) -> None:
        """Extract comprehensive metadata including Bruker GlobalMetadata.

        Split the same way :meth:`_extract_essential_metadata` is: the
        reader call is the boundary and keeps its breadth, the two
        extractions below it are Thyra reading an object Thyra built and
        are narrowed (issue #280).
        """
        try:
            # Use cached comprehensive metadata if available
            if self._comprehensive_metadata is not None:
                comp_meta = self._comprehensive_metadata
            elif self.reader is None:
                return
            else:
                comp_meta = self.reader.get_comprehensive_metadata()
                self._comprehensive_metadata = comp_meta
        except Exception as e:
            logger.debug(f"Could not read comprehensive metadata: {e}")
            return

        try:
            self._extract_bruker_metadata(metadata, comp_meta)
            self._extract_instrument_info(metadata, comp_meta)
        except MALFORMED_METADATA as e:
            logger.debug(f"Could not extract comprehensive metadata: {e}")

    @staticmethod
    def _extract_bruker_metadata(metadata: Dict[str, Any], comp_meta: Any) -> None:
        """Extract Bruker GlobalMetadata from comprehensive metadata."""
        if (
            hasattr(comp_meta, "raw_metadata")
            and "global_metadata" in comp_meta.raw_metadata
        ):
            metadata["GlobalMetadata"] = comp_meta.raw_metadata["global_metadata"]
            logger.debug(
                f"Extracted Bruker GlobalMetadata with keys: "
                f"{list(metadata['GlobalMetadata'].keys())}"
            )

    @staticmethod
    def _extract_instrument_info(metadata: Dict[str, Any], comp_meta: Any) -> None:
        """Extract instrument_info for fallback detection."""
        if hasattr(comp_meta, "instrument_info"):
            metadata["instrument_info"] = comp_meta.instrument_info
            logger.debug(f"Extracted instrument_info: {comp_meta.instrument_info}")

        # Extract format_specific for FlexImaging detection
        if hasattr(comp_meta, "format_specific"):
            metadata["format_specific"] = comp_meta.format_specific
            logger.debug(f"Extracted format_specific: {comp_meta.format_specific}")

        # Extract acquisition_params for additional detection
        if hasattr(comp_meta, "acquisition_params"):
            metadata["acquisition_params"] = comp_meta.acquisition_params
            logger.debug(
                f"Extracted acquisition_params: {comp_meta.acquisition_params}"
            )

    def _extract_spectrum_metadata(self, metadata: Dict[str, Any]) -> None:
        """Extract ImzML-specific spectrum metadata."""
        try:
            if hasattr(self.reader, "get_spectrum_metadata"):
                # Use cached spectrum metadata if available
                if self._spectrum_metadata is not None:
                    spec_meta = self._spectrum_metadata
                else:
                    spec_meta = self.reader.get_spectrum_metadata()
                    self._spectrum_metadata = spec_meta

                if spec_meta:
                    metadata.update(spec_meta)
        except Exception as e:
            logger.debug(f"Could not extract spectrum metadata: {e}")

    def _read_essential_metadata(self) -> Optional["EssentialMetadata"]:
        """The source's essential metadata, read at most once."""
        if self._essential_metadata is None and self.reader is not None:
            self._essential_metadata = self.reader.get_essential_metadata()
        return self._essential_metadata

    # ------------------------------------------------------------------
    # The method
    # ------------------------------------------------------------------

    def _select_method(self) -> None:
        """Pick the resampling method: the caller's, or the detector chain's."""
        config = self.config
        assert config is not None  # _select_method runs only with one
        method = config.method

        # If method is None or "auto", use DecisionTree to determine strategy
        if method is None:
            try:
                # Get metadata from reader for instrument detection
                metadata = self.detection_metadata()
                tree = ResamplingDecisionTree()
                detected_method = tree.select_strategy(metadata)
                logger.info(f"Auto-detected resampling method: {detected_method}")
                self._method = detected_method
            except NotImplementedError as e:
                logger.error(f"Auto-detection failed: {e}")
                logger.info("Falling back to nearest_neighbor for resampling")
                self._method = ResamplingMethod.NEAREST_NEIGHBOR
        else:
            # Use provided method directly (already an enum)
            self._method = method
            self._warn_if_override_contradicts_detector(method, config)

        self._method_was_auto = method is None

        logger.info(f"Using resampling method: {self._method}")

        if config.gap_tolerance_da is not None:
            logger.info(
                f"Interpolation gap tolerance: {config.gap_tolerance_da} Da "
                "(target bins farther than this from any source m/z are zeroed)"
            )

    def _warn_if_override_contradicts_detector(
        self, method: ResamplingMethod, config: Any
    ) -> None:
        """Say so when an explicit ``--resample-method`` overrules the detector.

        The detector has a verdict for every source; until #246 only the
        ``auto`` path ever asked for it, so an explicit method was applied
        with nothing checked and nothing said. The asymmetry is the whole
        defect: ``tic_preserving`` on a Bruker TDF -- for which the
        detector chooses nearest-neighbour -- interpolates across the gaps
        of a sparse centroid list and fills the axis. Measured on a
        713-frame PASEF acquisition: 423,386,757 stored non-zeros against
        302,106, a 583 MB table against 7.6 MB, 5.9 GB of peak RSS against
        0.5. Per-pixel TIC is identical either way, so the TIC identity
        cannot see it; what breaks is the siblings, and quietly (the
        heatmap marginal against the stored mean spectrum came to rel 68).

        This is the same bug class as #168 on PHI ToF-SIMS, which was fixed
        *by* adding a detector -- which is exactly why a detector is not
        enough on its own. A detector only steers ``auto``.

        Why a warning and not a refusal or an automatic gap tolerance is
        design decision D15 in ``docs/design-decisions.md``. In short: the
        remedy already has a flag, and nothing stored changes.

        Args:
            method: The method the caller asked for.
            config: The resampling config, read for a gap tolerance that
                is already in force.
        """
        try:
            detected = ResamplingDecisionTree().select_strategy(
                self.detection_metadata()
            )
        except Exception as exc:
            # Detection is advisory here. A source it cannot classify must
            # still convert with the method that was actually asked for.
            #
            # Broad on purpose, and it stays broad (issue #280). Narrowing
            # it to the exceptions a detector chain can raise looks right
            # and breaks ``test_a_source_the_detector_cannot_classify_
            # still_converts``, which pins the stronger promise this
            # docstring makes: *whatever* goes wrong in the advisory check,
            # the method the caller asked for is still applied. A
            # conversion that failed because an advisory check raised would
            # be a worse defect than any this catch can hide, and the catch
            # hides it at DEBUG rather than swallowing it silently.
            logger.debug(
                "Could not check --resample-method against the detector: %s",
                str(exc),
            )
            return

        if detected is method:
            return

        message = (
            "Resampling method %s was given explicitly, but this source's "
            "detector chose %s for it. "
        )
        args: List[Any] = [method.name, detected.name]

        if method is ResamplingMethod.TIC_PRESERVING:
            tolerance = getattr(config, "gap_tolerance_da", None)
            if tolerance is None:
                message += (
                    "Interpolating a source the detector reads as sparse "
                    "fills the whole axis: every bin between two measured "
                    "points gets a fabricated intensity, the stored matrix "
                    "grows by orders of magnitude, and per-pixel TIC still "
                    "balances so no total reveals it. Pass "
                    "--resample-gap-tolerance to discard bins no measured "
                    "m/z vouches for, or drop the override."
                )
            else:
                message += (
                    "--resample-gap-tolerance %s Da is set, so bins further "
                    "than that from a measured m/z are discarded rather "
                    "than interpolated across."
                )
                args.append(tolerance)

        logger.warning(message, *args)

    def _regate_tic_preserving(
        self, tree: ResamplingDecisionTree, axis_type: Any
    ) -> None:
        """Re-ask the TIC-preserving gate about the axis that will be built.

        The method is picked from the detector chain, which gates
        ``TIC_PRESERVING`` on the detector's *own* axis choice.
        ``--mass-axis-type`` overrides that choice, and it is applied in
        :meth:`resolve`, after the method has been chosen. The two then
        come apart and the gate's premise is gone: TIC-preserving
        resampling is exact only when the source grid law and the target
        axis are the same law.

        Both detectors that reach ``TIC_PRESERVING`` -- ``RapiflexDetector``
        (CONSTANT) and ``WatersProfileDetector`` (LINEAR_TOF) -- declare a
        source law equal to their own axis, so the early gate always cleared
        and the conversion then interpolated onto whatever axis was asked for.
        ``docs/resampling.md`` measures the cost: two ions of equal abundance
        come back with their ratio distorted by up to 13.4x. Nothing downstream
        sees it, because the per-pixel TIC still balances exactly -- preserving
        it is what the operator does.

        Auto-selected methods only. An explicit ``--resample-method`` is the
        caller's decision, honoured with a warning instead (D15); overruling it
        here is what #246 deliberately rejected. On the auto path with no axis
        override this re-derives the same answer, so it is a no-op.

        Args:
            tree: The decision tree to ask.
            axis_type: The axis the conversion will build.
        """
        if not self._method_was_auto:
            return
        if self._method is not ResamplingMethod.TIC_PRESERVING:
            return

        regated = tree.select_strategy_for_axis(self.detection_metadata(), axis_type)
        if regated is not self._method:
            logger.info(
                "Resampling method downgraded to %s: the mass axis resolved "
                "to %s, which is not the source grid law.",
                regated.name,
                axis_type.name,
            )
            self._method = regated

    # ------------------------------------------------------------------
    # The plan, and the axis
    # ------------------------------------------------------------------

    def _source_mass_range(self) -> Tuple[float, float]:
        """The source's own mass range, from the reader or the metadata dict."""
        essential = self._read_essential_metadata()
        if essential is not None:
            return essential.mass_range
        declared = (self._detection_metadata or {}).get("essential_metadata", {})
        mass_range = declared.get("mass_range")
        if mass_range is None:
            raise ValueError(
                "The mass range is unknown: the planner has neither a reader "
                "nor essential metadata to read it from."
            )
        return cast(Tuple[float, float], mass_range)

    def resolve(self) -> AxisPlan:
        """Resolve the mass range, the axis type and the bin count.

        Resolved once and remembered: the memory guard reads the bin count
        before the axis is laid and the build reads all four afterwards,
        and re-deriving them would run the detector chain a second time to
        reach the same answer.

        Returns:
            The four numbers a resampled axis follows from.

        Raises:
            ConversionRefused: If the range is empty or inverted, or fewer
                than two bins span it.
        """
        if self._plan is not None:
            return self._plan
        config = self.config
        if config is None:
            raise ValueError("Nothing is resampled, so there is no plan to resolve.")

        mass_range = self._source_mass_range()
        min_mz = mass_range[0] if config.min_mz is None else config.min_mz
        max_mz = mass_range[1] if config.max_mz is None else config.max_mz

        # An inverted range has no axis to lay. Unchecked, it built a
        # descending one, dropped every peak against it, reported "4 of 3
        # in the first spectrum affected" from a negative count, and
        # succeeded with an empty store (issue #250).
        if not (min_mz < max_mz):
            raise ConversionRefused(
                f"The resampling mass range [{min_mz:g}, {max_mz:g}] m/z is "
                "empty: the minimum has to be below the maximum. "
                + (
                    "Both come from the source's own mass range."
                    if config.min_mz is None and config.max_mz is None
                    else "Check --resample-min-mz / --resample-max-mz."
                )
            )

        tree = ResamplingDecisionTree()
        if config.axis_type is not None:
            axis_type = config.axis_type
        else:
            metadata = self.detection_metadata()
            axis_type = tree.select_axis_type(metadata)
            # A detector's width goes with the axis law it chose, so it is
            # consulted only on the auto path and only when the caller has
            # not set a width of their own.
            if config.mass_width_da is None:
                self._detected_reference_width = tree.select_reference_width(metadata)

        self._regate_tic_preserving(tree, axis_type)

        # A TOF axis without the caller's own coefficients takes the pair
        # the instrument declares -- on the auto path (an MRT centroid
        # conversion) and when --mass-axis-type tof was asked for by name
        # (a timsTOF opting in).
        if axis_type is AxisType.TOF and (config.tof_a is None or config.tof_b is None):
            self._detected_tof_law = tree.select_tof_law(self.detection_metadata())
        elif axis_type is not AxisType.TOF and config.tof_a is not None:
            logger.warning(
                "--tof-law was given but the mass axis resolved to %s, which "
                "does not use a width law; it is ignored. Pass "
                "--mass-axis-type tof to use it.",
                getattr(axis_type, "value", axis_type),
            )

        if config.mass_width_da is not None or config.target_bins is None:
            target_bins = self.bin_count(min_mz, max_mz, axis_type)
        else:
            target_bins = config.target_bins

        # A one-point axis has no bin width, and the log line that
        # reports one reduced an empty array: "zero-size array to
        # reduction operation minimum which has no identity", raised out
        # of initialization rather than said as a refusal (issue #250).
        if target_bins < 2:
            raise ConversionRefused(
                f"A mass axis needs at least 2 bins, got {target_bins}. "
                "One point is a single m/z value, not an axis: it has no bin "
                "width and nothing can be resampled onto it."
            )

        self._plan = AxisPlan(min_mz, max_mz, axis_type, target_bins)
        return self._plan

    def bin_count(self, min_mz: float, max_mz: float, axis_type: Any) -> int:
        """How many bins of the requested width span a range, for this source.

        :func:`bin_count_for_width` with what the detectors declared for
        this source filled in.
        """
        config = self.config
        if config is None:
            raise ValueError("Nothing is resampled, so no bin count is asked for.")
        return bin_count_for_width(
            config,
            min_mz,
            max_mz,
            axis_type,
            detected_width=self._detected_reference_width,
            detected_tof_law=self._detected_tof_law,
        )

    def reference_params(self, axis_type: Any) -> Tuple[float, float]:
        """``(width_da, reference_mz)`` for this source and this axis type."""
        config = self.config
        if config is None:
            raise ValueError("Nothing is resampled, so no bin width is asked for.")
        return reference_params(
            config,
            axis_name(axis_type),
            detected_width=self._detected_reference_width,
            detected_tof_law=self._detected_tof_law,
        )

    def settle(self, refuse_wide: Callable[[int], None]) -> SettledAxis:
        """The axis this conversion writes, built or taken from the source.

        Args:
            refuse_wide: The conversion's per-bin memory budget, given a
                bin count. Asked *before* a resampled axis is materialised
                and *after* a raw one has been handed over, because that is
                when each becomes knowable -- 300M bins were once refused
                by the count array's own ceiling in 2.3 s, but the process
                had reached 7.44 GB building the axis first (issue #251).

        Returns:
            The axis, and everything the conversion needs to know about it.
        """
        config_status = "SET" if self.config else "NOT SET"
        logger.info(f"Mass axis mode: resampling_config={config_status}")
        if self.config is not None:
            return self._build_axis(refuse_wide)
        return self._source_axis(refuse_wide)

    def _build_axis(self, refuse_wide: Callable[[int], None]) -> SettledAxis:
        """Lay the resampled axis the plan describes."""
        config = self.config
        assert config is not None  # settle() takes this route only with one
        logger.info(
            "Building RESAMPLED mass axis (resampling enabled) - "
            "will NOT iterate through all spectra"
        )
        plan = self.resolve()
        refuse_wide(plan.target_bins)

        min_mz, max_mz = plan.min_mz, plan.max_mz
        axis_type, target_bins = plan.axis_type, plan.target_bins

        # Determine reference parameters for physics generators
        reference_width, ref_mz = self.reference_params(axis_type)

        # Kept for provenance: the processing step must declare what was
        # actually done, and with "auto" settings the requested config
        # says nothing about the method, axis and bin width the decision
        # tree resolved to.  See the uns assembler's processing provenance.
        provenance: Dict[str, Any] = {
            "method": self._method,
            "axis_type": axis_type,
            "target_bins": target_bins,
            "min_mz": min_mz,
            "max_mz": max_mz,
            "mass_width_da": reference_width,
            "reference_mz": ref_mz,
        }
        tof_law: Optional[Tuple[float, float]] = None
        if axis_type is AxisType.TOF:
            a, b, k = tof_plan(config, detected_tof_law=self._detected_tof_law)
            tof_law = (a, b)
            provenance.update({"tof_a": a, "tof_b": b, "bins_per_fwhm": k})

        if config.axis_type is not None:
            logger.info(f"Using manually specified axis type: {axis_type}")
        else:
            logger.info(f"Auto-detected axis type: {axis_type}")

        logger.info(
            f"Building resampled mass axis: {min_mz:.2f} - {max_mz:.2f} m/z, "
            f"{target_bins} bins"
        )

        # Build the physics-based axis
        builder = CommonAxisBuilder()

        if hasattr(axis_type, "value") and axis_type.value != "constant":
            # Use physics-based generator with reference parameters
            mass_axis = builder.build_physics_axis(
                min_mz=min_mz,
                max_mz=max_mz,
                num_bins=target_bins,
                axis_type=axis_type,
                reference_mz=ref_mz,
                reference_width=reference_width,
                tof_law=tof_law,
            )
            logger.info(
                f"Built physics-based {axis_type} mass axis with "
                f"{len(mass_axis.mz_values)} points"
            )
        else:
            # Fall back to uniform axis
            mass_axis = builder.build_uniform_axis(min_mz, max_mz, target_bins)
            logger.info(
                f"Built uniform mass axis with " f"{len(mass_axis.mz_values)} points"
            )

        axis = mass_axis.mz_values.astype(np.float64)

        # Bin sizes for the log line, in chunks. ``np.diff`` over the whole
        # axis is another float64 array of its length -- 1.6 GB on a 200M
        # bin axis, allocated for two numbers in one INFO line (#251).
        min_bin_size, max_bin_size = bin_width_range(axis)

        # A positive narrowest gap is "strictly ascending", which the
        # closed-form bin index needs and this log line already computed.
        linearisation = usable_linearisation(
            mass_axis.linearisation, axis, min_bin_size
        )

        logger.info(
            f"Resampled mass axis created: {len(axis)} bins, "
            f"range {axis[0]:.2f}-{axis[-1]:.2f} m/z, "
            f"bin sizes {min_bin_size*1000:.2f}-{max_bin_size*1000:.2f} mDa "
            f"({axis_type})"
        )
        logger.info(f"Built resampled mass axis with {len(axis)} bins")

        return SettledAxis(
            axis=axis,
            # The range the bins were laid across, kept because a physics
            # axis reports bin *centres* and so stops half a bin short of
            # it at either end. It, not the axis's own span, is what
            # decides whether a peak is in range -- see kept_mz_range()
            # (issue #239).
            axis_range=(float(min_mz), float(max_mz)),
            linearisation=linearisation,
            method=self._method,
            provenance=provenance,
        )

    def _source_axis(self, refuse_wide: Callable[[int], None]) -> SettledAxis:
        """Take the source's own m/z values as the axis (``--no-resample``)."""
        if self.reader is None:
            raise ValueError("A raw mass axis needs a reader to read it from.")
        if self.reader.has_shared_mass_axis:
            logger.info(
                "Loading RAW mass axis (no resampling) - "
                "continuous mode, reading m/z from first spectrum only"
            )
        else:
            logger.warning(
                "Building RAW mass axis (no resampling) - "
                "processed mode, iterating ALL spectra to collect unique m/z values. "
                "This is slow for large datasets!"
            )
        axis = self.reader.get_common_mass_axis()
        if len(axis) == 0:
            raise ConversionRefused(
                "Common mass axis is empty. Cannot proceed with conversion."
            )
        refuse_wide(len(axis))
        logger.info(f"Using raw mass axis with {len(axis)} unique m/z values")

        # A raw union axis was laid by no generator, so there is no
        # coordinate it is uniform in and every placement onto it is a
        # search -- which is what _map_mass_to_indices does, and no
        # strategy is built for it.
        return SettledAxis(
            axis=axis,
            axis_range=None,
            linearisation=None,
            method=None,
            provenance=None,
        )
