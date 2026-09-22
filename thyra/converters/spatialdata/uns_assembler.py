# thyra/converters/spatialdata/uns_assembler.py

"""The provenance block every table's ``uns`` carries, built from arguments.

One conversion owns one :class:`UnsAssembler`, and it is the only thing
that composes ``uns``. It exists because the converters once had two
write paths and they drifted: the in-memory route handed the table to
``anndata``'s writer, the streaming route hand-wrote the Zarr layout and
composed its own, much smaller block. There is one write path now, and
one place the block is assembled.

**It is handed what it needs, not the converter.** Everything the block
reads that the converter decides *during* a conversion -- the sibling
keys, the mobility grid, the resolved resampling plan, the region
summary, where the pixel pitch came from -- arrives per call in an
:class:`UnsContext`, because every one of those is assigned after the
converter is constructed and several change between one table and the
next. Only the four things settled at construction time (the reader, the
pitch accessor, the pixel-size detection info and the resampling config)
are constructor arguments.

**It lives with the converter, not with the metadata package**, because
only a conversion runs it: its inputs are converter types (a pixel-size
source, a resampling config, a mobility grid, a reader) and its one caller
is :class:`~.base_spatialdata_converter.BaseSpatialDataConverter`. It
landed under ``thyra/metadata/`` first because it assembles metadata;
by dependency it is a converter collaborator, and the document side of
that package (``schema/``, ``ontology/``, ``types.py``) imports nothing
from here or from any converter (design decision D25).

**Exceptions are logged as text here, never as objects**, for the reason
``base_spatialdata_converter``'s module docstring gives at length: a
retained log record drags an exception's traceback along, and through it
every CSC scratch memmap the finalize path's frames still hold (issue
#249).
"""

import json
import logging
from dataclasses import dataclass
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Tuple,
)

import numpy as np

from ...errors import MALFORMED_METADATA
from ...metadata.schema import (
    MSI_METADATA_UNS_KEY,
    ProcessingStep,
    SoftwareRef,
    build_msi_metadata,
    forget_resolved_table,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ...core.base_converter import PixelSizeSource
    from ...core.base_reader import BaseMSIReader
    from ...resampling.mobility_grid import MobilityGrid
    from ...resampling.types import ResamplingConfig

logger = logging.getLogger(__name__)


#: How far the heatmap's total may sit from the stored mean spectrum's and
#: still be called equal. Looser than the sibling tables'
#: ``_MARGINAL_TOLERANCE`` because the two are not the same sum reordered:
#: the heatmap coarsens the mass axis by an integer factor and clips
#: mobility into the edge channels, so float32 storage of the counts sets
#: the floor. Measured at 4e-8 under scan_sum and 0.15 under the vendor
#: centroid, so the gap either side of this is four orders of magnitude
#: wide.
_HEATMAP_TOLERANCE = 1e-4


def _numeric_only(value: Any) -> bool:
    """True when a list round-trips through AnnData/zarr as a numeric array."""
    if isinstance(value, list):
        return all(_numeric_only(item) for item in value)
    if isinstance(value, np.ndarray):
        return value.dtype.kind in "iufb"
    return isinstance(value, (bool, int, float, np.integer, np.floating, np.bool_))


def _json_fallback(obj: Any) -> Any:
    """Last-resort encoder for values ``json.dumps`` does not know."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    return str(obj)


def _jsonify_string_lists(obj: Any) -> Any:
    """Replace any list that is not purely numeric with its JSON encoding.

    AnnData/zarr cannot round-trip such lists: a list of dicts is
    stringified entry by entry into Python ``repr`` strings, and any list
    of strings comes back as a numpy string array -- whose ``deepcopy``
    segfaults outright on numpy 2.1-2.2 (numpy#28609). Since every table
    copy deepcopies ``uns`` (``AnnData.copy``, spatial queries, joins),
    one such array in ``uns`` kills the reader's process with no
    traceback. JSON side-steps both: it stores a single scalar string,
    and hands consumers the actual structure back through
    ``json.loads`` instead of ``repr`` output.

    Purely numeric lists (nested included) are kept: they become plain
    numeric arrays, which are safe and more useful as arrays.
    """
    if isinstance(obj, dict):
        return {key: _jsonify_string_lists(value) for key, value in obj.items()}
    if isinstance(obj, list):
        if _numeric_only(obj):
            return obj
        return json.dumps(obj, default=_json_fallback)
    return obj


def serialize_for_zarr(obj: Any) -> Any:
    """Recursively convert tuples to lists for Zarr serialization.

    Dict keys are coerced to strings: Zarr group members must be named,
    and a non-string key otherwise fails at write time, after the whole
    conversion has already been done.
    """
    if isinstance(obj, dict):
        return {str(k): serialize_for_zarr(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [serialize_for_zarr(item) for item in obj]
    elif hasattr(obj, "__dict__"):
        # Convert dataclass/object to dict
        return serialize_for_zarr(vars(obj))
    else:
        return obj


def record_heatmap_current(adata: Any) -> None:
    """Say how much of the stored mean spectrum the heatmap holds.

    The heatmap's marginal over mobility is the stored mean spectrum
    coarsened to its own m/z bins -- exactly, under the default
    ``--tdf-spectrum scan_sum``, and not at all under the vendor
    centroid, which is a peak-picked spectrum over the same scans
    while the heatmap is built from the raw points. Measured at 15%
    off on every vendor_centroid store.

    docs/output-format.md says the identity holds only under
    scan_sum, but nothing in the store said which case a given store
    was, so a consumer plotting the heatmap next to the mean spectrum
    had no way to tell the 15% from a bug (issue #253). The grid's
    ``uns["mobility_marginal"]`` has recorded its ratio since it
    existed; this is the same number for the heatmap.

    Both quantities are per-pixel means over the same pixels, so
    their totals compare directly: the ratio is one scalar and needs
    nothing materialised.
    """
    block = adata.uns.get("mobility_heatmap")
    spectrum = adata.uns.get("average_spectrum")
    if not isinstance(block, dict) or spectrum is None:
        return
    try:
        stored = float(np.asarray(spectrum, dtype=np.float64).sum())
        held = float(np.asarray(block["counts"], dtype=np.float64).sum())
    except MALFORMED_METADATA as e:  # pragma: no cover - defensive
        # A block without ``counts``, or counts that will not become a
        # float array (issue #280). Two numpy sums cannot fail any
        # other way that is not a defect.
        # str(e), not e: this is the finalize path, where a retained
        # record pins the memmaps the table is built on (issue #249).
        logger.debug("Could not compare the mobility heatmap: %s", str(e))
        return
    if stored <= 0:
        return
    ratio = held / stored
    block["current_ratio"] = ratio
    if abs(ratio - 1.0) > _HEATMAP_TOLERANCE:
        logger.warning(
            "The mass-mobility heatmap holds %.4fx the stored mean "
            "spectrum's ion current. Its marginal over mobility is that "
            "spectrum only under --tdf-spectrum scan_sum; the vendor "
            "centroid keeps only the current inside the peaks it picks, "
            "while the heatmap reads raw points. "
            'uns["mobility_heatmap"]["current_ratio"] records it.',
            ratio,
        )


@dataclass(frozen=True)
class UnsContext:
    """What the converter knows at the moment one table's ``uns`` is built.

    Every field is state the converter assigns after it is constructed,
    and several of them change between one table and the next:
    ``mobility_table_key`` and ``msms_table_key`` are decided per slice
    and set back to ``None`` when a sibling builder declines,
    ``resolved_resampling_plan`` is filled when the resampled axis is
    built, ``mobility_grid`` when the grid is resolved, ``region_info``
    when the conversion is initialised, and ``pixel_size_source`` is
    reassigned when a pitch is auto-detected from the reader's metadata.
    Passing them per call is what keeps the assembler from reaching back
    into the converter for any of it.

    Attributes:
        mobility_table_key: Element key of the mobility-resolved sibling
            this slice gets, or ``None``.
        msms_table_key: Element key of the demultiplexed MS/MS sibling
            this slice gets, or ``None``.
        mobility_grid: The resolved common mobility grid, or ``None``.
            Only ``to_schema_report()`` is read.
        resolved_resampling_plan: What the resampling decision tree
            actually resolved to, which wins over the requested config in
            the ``mass axis resampling`` step.
        region_info: The acquisition region summary, stored as JSON.
        pixel_size_source: Where the in-plane pitch came from; read as
            ``.value``.
        mobility_heatmap: The conversion's mass-mobility heatmap block, as
            a callable so it is built at the point the block is recorded
            rather than at the point the context is packed.
    """

    mobility_table_key: Optional[str]
    msms_table_key: Optional[str]
    mobility_grid: Optional["MobilityGrid"]
    resolved_resampling_plan: Optional[Dict[str, Any]]
    region_info: Optional[list]
    pixel_size_source: "PixelSizeSource"
    mobility_heatmap: Callable[[], Optional[Dict[str, Any]]]


class UnsAssembler:
    """The ``uns`` block of one conversion: composed, applied, unnamed again.

    Three things the converter asks of it:

    1. :meth:`build` -- the mapping of ``uns`` key to value that every
       table of the conversion renders, the summed one and its siblings
       alike, so a section added here reaches every store.
    2. :meth:`apply` -- the same mapping written onto an ``AnnData`` that
       is about to be stored, plus the heatmap's ion-current check, which
       can only run once the mean spectrum is on the object.
    3. :meth:`unname_declined_siblings` -- taking a sibling's key back out
       of every table that already named it, when its builder declined.

    :meth:`fragmentation` is the fourth entry point and the odd one: the
    converter's own MS/MS planning asks for the schedule before any
    ``uns`` is built, and the schedule must be read once per conversion
    with its chimera warning said once. The cache lives here because this
    is where both readers of it meet.
    """

    def __init__(
        self,
        reader: "BaseMSIReader",
        *,
        pixel_size_xy: Callable[[], Tuple[float, float]],
        pixel_size_detection_info: Optional[Dict[str, Any]],
        resampling_config: Optional["ResamplingConfig"],
    ) -> None:
        """Build the ``uns`` assembler of one conversion.

        Args:
            reader: The MSI reader. Only ``get_comprehensive_metadata``,
                ``has_ion_mobility``, ``get_mobility_axis``,
                ``has_fragmentation``, ``get_fragmentation`` and -- through
                ``getattr`` -- ``tdf_spectrum`` and ``file_type`` are read,
                so a reader here is whatever satisfies that interface, not
                necessarily a
                :class:`~thyra.core.base_reader.BaseMSIReader` subclass.
            pixel_size_xy: The conversion's in-plane pitch as
                ``(x_um, y_um)``, called when the block is built rather
                than read here: the converter settles its pitch from the
                reader's metadata after this object is built.
            pixel_size_detection_info: What pitch detection reported, or
                ``None``. Only ``source_format`` is read.
            resampling_config: The requested resampling configuration, or
                ``None`` when the mass axis is left as the source has it.
                Serialised field by field into the processing provenance.
        """
        self.reader = reader
        self._pixel_size_xy = pixel_size_xy
        self._pixel_size_detection_info = pixel_size_detection_info
        self._resampling_config = resampling_config
        # The reader's fragmentation schedule, read once per conversion.
        # ``_read`` distinguishes "not yet" from "asked, nothing to say".
        self._fragmentation_schedule: Any = None
        self._fragmentation_read = False

    #: The ``uns`` block outside the versioned schema that names each
    #: sibling, keyed by the ``ms_analysis`` section that names it inside.
    SIBLING_UNS_BLOCK = {
        "ion_mobility": "mobility_axis",
        "fragmentation": "msms_schedule",
    }

    def build(self, ctx: UnsContext) -> Dict[str, Any]:
        """The provenance block every write path must persist, identically.

        Single source of truth for what lands in the table's ``uns``:

        - ``essential_metadata`` -- dimensions, mass range, source path,
          spectrum type and the Thyra version that wrote the store.
        - ``format_specific`` -- vendor metadata (FlexImaging areas,
          teaching points, imzML file mode, ...).
        - ``acquisition_params`` / ``instrument_info`` -- when the reader
          has them.
        - ``raw_metadata`` -- the source metadata as read.
        - ``regions`` -- the acquisition region summary, as JSON.

        This exists because the converters once had two write paths and
        they drifted. The in-memory converters handed the table to
        ``anndata``'s writer, which serialises whatever is in
        ``adata.uns``; the streaming path hand-wrote the Zarr layout and
        composed its own, much smaller block -- with ``spectrum_type``
        hardcoded to ``"processed"``, which is not even a value the
        extractors produce. Routing was on a size threshold at the time,
        so a dataset large enough to reach the streaming path came out
        claiming a spectrum representation it did not have, and without
        any of the other sections, while a slightly smaller one from the
        same instrument came out complete. There is one write path now,
        through ``anndata``'s writer, and every table -- the summed one
        and its siblings -- renders this mapping, so a section added here
        reaches every store.

        Sections the reader has nothing for are omitted rather than
        written empty, so consumers can tell "not available from this
        format" from "available and empty".

        Args:
            ctx: What the converter knows about the table being built.

        Returns:
            Mapping of ``uns`` key to the value to store. Empty if the
            reader cannot produce comprehensive metadata at all.
        """
        try:
            comp_meta = self.reader.get_comprehensive_metadata()
        except Exception as e:
            # Non-fatal: a store without provenance still holds the
            # spectra. But it is not a debug-level event -- the whole
            # point of the block is that a consumer can say where the
            # data came from, so losing it has to be visible in the log.
            logger.warning("Could not read metadata for uns provenance: %s", str(e))
            return {}

        uns: Dict[str, Any] = {}
        try:
            self._collect_essential_metadata(uns, comp_meta)
            self._collect_optional_sections(uns, comp_meta)
            self._collect_region_info(uns, ctx)
        except (*MALFORMED_METADATA, RecursionError) as e:
            # Three Thyra methods reading an object a Thyra extractor built,
            # so the breadth the reader call above has is not earned here
            # (issue #280). ``RecursionError`` is in the set because
            # serialize_for_zarr walks ``vars()`` and a vendor object that
            # points back at its parent has no bottom.
            logger.warning("Could not build the full uns provenance block: %s", str(e))

        self._collect_msi_metadata_block(uns, comp_meta, ctx)
        self._collect_mobility_axis(uns, ctx)
        self._collect_mobility_heatmap(uns, ctx)
        self._collect_msms_schedule(uns, ctx)

        return uns

    def apply(self, adata: Any, ctx: UnsContext) -> None:
        """Apply :meth:`build` to an AnnData about to be written."""
        uns = self.build(ctx)
        adata.uns.update(uns)
        record_heatmap_current(adata)
        logger.debug("Added MSI metadata to AnnData .uns: %s", sorted(uns))

    def fragmentation(self) -> Any:
        """The reader's fragmentation schedule, read once and cached.

        ``None`` when the reader cannot say. Asked for through the base
        reader contract, so a format that learns to report it later needs
        no change here.
        """
        if not self._fragmentation_read:
            self._fragmentation_read = True
            if self.reader.has_fragmentation:
                try:
                    self._fragmentation_schedule = self.reader.get_fragmentation()
                except Exception as e:  # pragma: no cover - reader-defined
                    logger.warning("Could not describe the fragmentation: %s", str(e))
                    self._fragmentation_schedule = None
            self._warn_if_precursors_merge()
        return self._fragmentation_schedule

    def _fragmentation_report(self) -> Any:
        """The schedule in the shape the schema builder reads, or ``None``."""
        schedule = self.fragmentation()
        return None if schedule is None else schedule.to_extractor_report()

    def _warn_if_precursors_merge(self) -> None:
        """Say out loud when a stored spectrum sums several precursors.

        The stored spectrum of such a pixel holds fragments of every
        precursor the frame isolated, with nothing marking which came
        from which. That is not visible in the output -- it looks like an
        ordinary spectrum -- so it is said once, at WARNING, rather than
        left for a reader of the peaks to work out.
        """
        schedule = self._fragmentation_schedule
        if schedule is None or not schedule.merges_precursors:
            return
        # Each precursor once: a mass isolated at two collision energies
        # is two windows and one precursor (issue #383).
        distinct = list(dict.fromkeys(f"{w.target:g}" for w in schedule.windows))
        targets = ", ".join(distinct[:6])
        if len(distinct) > 6:
            targets += ", ..."
        logger.warning(
            "This acquisition isolates %d precursors per pixel (%s) over %d "
            "isolation windows. Thyra sums them into one spectrum per pixel, "
            "so the stored spectrum holds fragments of all of them and cannot "
            "be attributed to a single precursor. The schedule is recorded in "
            "uns['msms_schedule'].",
            schedule.n_precursors,
            targets,
            len(schedule.windows),
        )

    def _collect_msms_schedule(self, uns: Dict[str, Any], ctx: UnsContext) -> None:
        """Add ``msms_schedule`` when the source fragmented anything.

        Written on the summed MSI table so a consumer can tell fragment
        m/z from intact m/z, and see which precursors a chimeric spectrum
        merges. Kept out of the versioned ``msi_metadata`` block for the
        same reason ``mobility_axis`` is: that block is versioned, this
        one carries arrays.
        """
        schedule = self.fragmentation()
        if schedule is None or not schedule.is_msms:
            return
        block = schedule.to_uns()
        if ctx.msms_table_key is not None:
            block["resolved_table"] = ctx.msms_table_key
        uns["msms_schedule"] = _jsonify_string_lists(serialize_for_zarr(block))

    def _collect_mobility_axis(self, uns: Dict[str, Any], ctx: UnsContext) -> None:
        """Add ``mobility_axis`` when the source has a mobility dimension.

        Written on the summed MSI table so a consumer can tell "summed over
        mobility" from "never had any", and so it can find the
        mobility-resolved sibling (``resolved_table``) when one was written.
        Kept out of the ``msi_metadata`` schema block: that block is
        versioned, this one carries arrays.
        """
        try:
            if not self.reader.has_ion_mobility:
                return
            axis = self.reader.get_mobility_axis()
        except Exception as e:  # pragma: no cover - reader-defined
            logger.warning("Could not describe the mobility axis: %s", str(e))
            return
        if axis is None:
            return
        block = axis.to_uns()
        if ctx.mobility_table_key is not None:
            block["resolved_table"] = ctx.mobility_table_key
        uns["mobility_axis"] = _jsonify_string_lists(serialize_for_zarr(block))

    def _collect_mobility_heatmap(self, uns: Dict[str, Any], ctx: UnsContext) -> None:
        """Add ``mobility_heatmap`` when the source has a mobility dimension.

        The mean (m/z, mobility) frame of the whole dataset, the surface
        a consumer looks at to decide whether mobility separates anything
        before asking for a mobility-resolved table. Arrays, not lists,
        and plain-name keys; nothing in it is per pixel.
        """
        block = ctx.mobility_heatmap()
        if block is not None:
            uns["mobility_heatmap"] = block

    def _collect_msi_metadata_block(
        self, uns: Dict[str, Any], comp_meta: Any, ctx: UnsContext
    ) -> None:
        """Add the versioned ``msi_metadata`` schema block.

        Built in its own try so a schema failure cannot take the other
        provenance sections down with it.  See docs/metadata-schema.md
        for the storage contract and ``thyra validate`` for the
        consumer side.
        """
        try:
            info = self._pixel_size_detection_info or {}
            meta = build_msi_metadata(
                comp_meta,
                pixel_size_um=self._pixel_size_xy(),
                pixel_size_source=ctx.pixel_size_source.value,
                source_format=info.get("source_format"),
                processing=self._processing_provenance(ctx),
                mobility_resolved_table=ctx.mobility_table_key,
                mobility_grid=(
                    None
                    if ctx.mobility_grid is None
                    else ctx.mobility_grid.to_schema_report()
                ),
                fragmentation=self._fragmentation_report(),
                msms_resolved_table=ctx.msms_table_key,
            )
            uns[MSI_METADATA_UNS_KEY] = meta.to_uns_dict()
        except MALFORMED_METADATA as e:
            # The schema builder is Thyra's, and so is everything handed to
            # it, so its breadth was not earned (issue #280). What remains
            # is the vendor-shaped part: a comp_meta section the builder
            # reads positionally or by key and this source spells
            # differently.
            logger.warning("Could not build the msi_metadata block: %s", str(e))

    def _processing_provenance(self, ctx: UnsContext) -> List[Any]:
        """The processing steps this conversion performed, oldest first.

        Modeled on mzQC provenance.  The list describes what was done to
        the data, so nothing about how the store was written -- the sparse
        layout, the number of passes -- belongs here.
        """
        from thyra import __version__

        thyra_ref = SoftwareRef(name="thyra", version=__version__)
        conversion_parameters: Dict[str, Any] = {}
        # A TDF frame's mobility ramp is collapsed into one spectrum, and
        # the two correct ways to do that differ in TIC by up to a fifth
        # (vendor centroid vs. lossless scan sum), so the store must say
        # which one it holds.
        tdf_spectrum = getattr(self.reader, "tdf_spectrum", None)
        if (
            tdf_spectrum is not None
            and getattr(self.reader, "file_type", None) == "tdf"
        ):
            conversion_parameters["tdf_spectrum"] = str(tdf_spectrum)
        steps = [
            ProcessingStep(
                name="conversion",
                software=thyra_ref,
                parameters=conversion_parameters,
            )
        ]

        config = self._resampling_config
        if config is not None:
            parameters: Dict[str, Any] = {}
            for field_name, value in vars(config).items():
                if value is None:
                    continue
                parameters[field_name] = getattr(value, "value", value)
            # The resolved plan wins over the requested config: with
            # "auto" settings the config says nothing about the method
            # and axis the decision tree actually picked, and the step
            # must declare what was done, not what was asked for.
            for field_name, value in (ctx.resolved_resampling_plan or {}).items():
                if value is None:
                    continue
                parameters[field_name] = getattr(value, "value", value)
            steps.append(
                ProcessingStep(
                    name="mass axis resampling",
                    software=thyra_ref,
                    parameters=parameters,
                )
            )
        return steps

    def _collect_essential_metadata(self, uns: Dict[str, Any], comp_meta: Any) -> None:
        """Add ``essential_metadata`` (tuples become lists for Zarr)."""
        essential = getattr(comp_meta, "essential", None)
        if essential is None:
            return

        dims = essential.dimensions
        mrange = essential.mass_range
        from thyra import __version__

        uns["essential_metadata"] = {
            "source_path": str(essential.source_path),
            "dimensions": list(dims) if dims else None,
            "mass_range": list(mrange) if mrange else None,
            "spectrum_type": getattr(essential, "spectrum_type", None),
            "thyra_version": __version__,
        }

    def _collect_optional_sections(self, uns: Dict[str, Any], comp_meta: Any) -> None:
        """Add the vendor sections the reader actually populated.

        Lists that are not purely numeric (imzML ``cvParams``, or any
        string list a vendor extractor reports) are stored as JSON
        strings -- see :func:`_jsonify_string_lists` for why letting
        them reach the writer as lists corrupts them and crashes
        readers on numpy 2.1-2.2.
        """
        for key in ("format_specific", "acquisition_params", "instrument_info"):
            value = getattr(comp_meta, key, None)
            if value:
                uns[key] = _jsonify_string_lists(serialize_for_zarr(value))

        raw_metadata = getattr(comp_meta, "raw_metadata", None)
        if raw_metadata:
            uns["raw_metadata"] = _jsonify_string_lists(
                serialize_for_zarr(raw_metadata)
            )

    def _collect_region_info(self, uns: Dict[str, Any], ctx: UnsContext) -> None:
        """Add the acquisition region summary as JSON.

        Stored as a JSON string because AnnData/zarr cannot round-trip
        a list of dicts (they get stringified individually). JSON
        preserves the structure and can be parsed with json.loads().

        Always written for a consistent schema. Single-region datasets
        get a single-entry list with region_number=1.
        """
        region_info = ctx.region_info
        if region_info:
            uns["regions"] = json.dumps(serialize_for_zarr(region_info))

    def unname_declined_siblings(
        self,
        tables: Mapping[str, Any],
        keys: Iterable[Optional[str]],
        declined: List[str],
    ) -> None:
        """Take a declined sibling's key back out of every table that names it.

        The alternative was to decide before naming -- hoist whatever makes
        a builder decline up into ``SiblingTables._plan_mobility_table``
        and ``_plan_msms_table``, so a table that will not be built is
        never named. That is the cleaner shape and it is not the one
        taken: the two builders decline at eight separate points across
        ``mobility_table`` and ``msms_table``, several of them knowable
        only once the pass has run (a feature listing that comes back
        empty, a var count over the ceiling), and a ninth added later
        would re-open the defect silently. Reacting to ``None`` in one
        place closes all of them, including the ones nobody has written
        yet.

        Every table written for this slice is cleaned, not just the summed
        one: the converter takes the siblings' ``uns`` from :meth:`build`
        *before* either builder runs, so a surviving sibling carries the
        same stale pointer. A mobility table that declines while the MS/MS
        table is written would otherwise leave the MS/MS table naming it
        too.

        Args:
            tables: The tables written for this slice, keyed by element
                key. Looked up, never replaced.
            keys: The element keys to clean -- the summed table's and each
                sibling's. A ``None`` key is skipped, as is one naming a
                table that has no ``uns``.
            declined: The ``ms_analysis`` sections whose sibling declined.
        """
        for key in keys:
            if key is None:
                continue
            uns = getattr(tables.get(key), "uns", None)
            if uns is None:
                continue
            for section in declined:
                block = uns.get(self.SIBLING_UNS_BLOCK[section])
                if isinstance(block, dict):
                    block.pop("resolved_table", None)
                meta = uns.get(MSI_METADATA_UNS_KEY)
                if isinstance(meta, dict):
                    forget_resolved_table(meta, section)


__all__ = [
    "UnsAssembler",
    "UnsContext",
    "record_heatmap_current",
    "serialize_for_zarr",
]
