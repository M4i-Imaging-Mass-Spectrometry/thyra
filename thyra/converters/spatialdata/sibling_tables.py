# thyra/converters/spatialdata/sibling_tables.py

"""The tables written *beside* the summed MSI table, and what decides them.

A conversion writes one summed table per slice and, for a source that
carries more than a spectrum per pixel, up to two siblings of it: the
mobility-resolved ``{table}_mobility`` (see :mod:`.mobility_table`) and
the demultiplexed ``{table}_msms`` (see :mod:`.msms_table`). The builders
of both, and the pass engine that can feed them from the summed table's
own reads (:mod:`.fused_passes`), have always been their own modules.
This is the orchestration between them, which was not: seventeen
converter methods over seventeen attributes, deciding *whether* each
sibling is written, resolving the grid one of them needs, running or
fusing the raw pass they share, and taking a name back out of the store
when a builder declines.

**It exists because none of that was reachable without a conversion.**
Issue #343 is the shape of defect that produces: the summed table's
``uns`` is composed before either builder runs, so it names siblings that
have not been built yet, and a builder that declines by returning ``None``
used to leave its name behind -- a store whose summed table pointed at an
element nobody wrote, reported as a success. It was found, and could only
be tested, by driving ``convert()`` end to end with a stub reader. The
same decision is now one method call on a collaborator that needs a
reader and the store it would write and nothing else.

**Three orderings here are load-bearing and all three read like
oversights**, so each is stated where it happens:

* :meth:`SiblingTables._plan_msms_table` asks
  :func:`~thyra.converters.spatialdata.msms_table.demultiplex_refusal`
  *before* the reader's capability check, or an acquisition with a single
  precursor is told its reader cannot separate precursors at all.
* :meth:`SiblingTables.release_scratch` empties the whole tables mapping
  before unlinking anything: an AnnData over a memmap keeps the file
  mapped, and Windows will not delete a mapped file.
* :meth:`SiblingTables.attach` composes the siblings' ``uns`` once, up
  front, from the summed table's -- so every table of the slice carries
  the same pointers, and the unnaming after a decline has to reach all of
  them (issue #343).

**The sibling modules are imported inside the methods that use them**,
as they were on the converter. The cost is not the reason -- the whole
chain is 2.5 ms against this package's 2.8 s import -- the reason is that
a deferred import is the seam the ``#343`` tests patch
(``mobility_table.build_mobility_table``) to make a builder decline.
Binding those names at module import time would make those patches
silently ineffective.

**Exceptions are logged as text here, never as objects**, for the reason
:mod:`.base_spatialdata_converter`'s module docstring gives at length: a
retained log record holds the exception, the exception holds its
traceback, and the traceback's frames reach every CSC scratch memmap the
finalize path built a table over -- which is the scratch this module is
also responsible for removing (issue #249).
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from ...core.base_reader import BaseMSIReader
from ...errors import MALFORMED_METADATA
from ...resampling.mobility_grid import (
    MOBILITY_CHANNELS,
    MobilityGrid,
    build_mobility_grid,
    report_channel_width,
)
from ...resampling.types import AxisLinearisation

logger = logging.getLogger(__name__)


#: How far a marginal may sit from the column it mirrors, relative to the
#: largest value in the summed table, and still be called exact. Summing
#: the same float64 values in a different order is the only difference
#: under ``--tdf-spectrum scan_sum``.
_MARGINAL_TOLERANCE = 1e-9


def _current_ratio_block(
    table: Any, summed_key: str, summed: Any
) -> Optional[Dict[str, Any]]:
    """How much of the summed table's ion current a sibling holds, per pixel.

    A sibling's row sum over every one of its columns is that pixel's ion
    current, and so is the summed table's, so the two row sums compare
    directly. Both matrices sit on memmaps, and a row sum is one pass
    over each -- nothing the size of a matrix is held in RAM, which is
    why this is the whole comparison: a cell-by-cell deviation between
    the grid table's marginal and the summed table needs the product and
    the difference materialised, each as large as the summed table, and
    the route that writes both never holds either. ``None`` when the two
    cannot be compared (a row count mismatch, no ion current at all),
    which is not a disagreement.
    """
    if summed is None:
        return None
    totals = np.asarray(summed.X.sum(axis=1)).ravel().astype(np.float64)
    split = np.asarray(table.X.sum(axis=1)).ravel().astype(np.float64)
    if split.size != totals.size or not totals.any():
        return None
    per_pixel = split / np.where(totals == 0, np.nan, totals)
    return {
        "summed_table": summed_key,
        "current_ratio": float(split.sum() / totals.sum()),
        "current_ratio_pixel_min": float(np.nanmin(per_pixel)),
        "current_ratio_pixel_max": float(np.nanmax(per_pixel)),
    }


@dataclass(frozen=True)
class SiblingContext:
    """What the conversion decided after the converter was constructed.

    The same split the ``uns`` assembler draws, for the same reason: what
    is settled when the converter is built is a constructor argument of
    :class:`SiblingTables`, and what a conversion decides later arrives
    per call. Here that is the common mass axis and everything downstream
    of it -- none of which exists until ``_setup_mass_axis`` has run,
    which is after every collaborator is constructed -- plus the three
    converter-side delegates the orchestration has to call back into.

    The three delegates are packed fresh on every call rather than
    captured, which is what lets a test replace the converter's ``uns``
    assembler and have this see the replacement (the failure a captured
    collaborator reproduced in issue #288).

    Attributes:
        mass_axis: The common mass axis of this conversion, or ``None``
            before it is built. Every sibling is binned onto it, so a
            conversion that has not laid one writes no sibling at all.
        linearisation: ``mass_axis``'s own closed-form index mapping when
            it has one, so the sinks place points by arithmetic rather
            than by search (design decision D21). Belongs to the axis,
            not to the strategy that bins with it: an interpolated
            conversion still writes siblings onto the same axis.
        dimensions: The source's ``(n_x, n_y, n_z)`` grid, or ``None``
            before the reader's essential metadata is loaded. Only the
            fused route needs it, to size one row space.
        n_spectra: How many spectra the source holds, for the raw passes'
            progress bars. A callable because on a Bruker source it is a
            query, and the passes that want it may never run.
        build_uns: The provenance block the summed table of this slice
            carries, which the siblings inherit. A callable because it
            must be composed at the moment :meth:`SiblingTables.attach`
            runs, not when the context is packed.
        fragmentation: The source's fragmentation schedule, or ``None``.
            Read through the ``uns`` assembler, which caches it, because
            ``demultiplex_refusal`` and the ``uns`` block ask the same
            question of the same reader.
        unname_declined: Takes a declined sibling's key back out of every
            table of the slice that names it. ``uns`` surgery, so it
            lives on the assembler; the decision to call it is here.
    """

    mass_axis: Optional[NDArray[np.float64]]
    linearisation: Optional[AxisLinearisation]
    dimensions: Optional[Tuple[int, int, int]]
    n_spectra: Callable[[], int]
    build_uns: Callable[[], Dict[str, Any]]
    fragmentation: Callable[[], Any]
    unname_declined: Callable[
        [Mapping[str, Any], Iterable[Optional[str]], List[str]], None
    ]


class SiblingTables:
    """The sibling tables of one conversion: planned, scanned, built, named.

    Five things the converter asks of it, in this order:

    1. :meth:`fused_passes` -- when the reader hands its frames over as
       records, the sinks to feed from the summed table's own two passes
       instead of scanning again. Plans the siblings on the way, because
       the sinks are theirs.
    2. :meth:`plan` -- which siblings this slice gets, decided *before*
       the summed table's ``uns`` is composed so the two agree. Idempotent
       per table key, so the route that already planned in step 1 does not
       re-decide and the route that did not still plans every slice.
    3. :meth:`prepare_scans` -- the raw mobility pass, once, for
       everything that still needs it.
    4. :meth:`attach` -- build each planned sibling, add it to the slice's
       tables, and take back the name of any that declined.
    5. :meth:`release_scratch` -- drop every table's memmaps and remove
       the scratch directories, on every way out of a conversion.

    :meth:`heatmap` sits beside them: the mass-mobility navigator block is
    a sibling concern (it is built from the same raw pass) that lands in
    the summed table's own ``uns`` rather than in a table of its own.
    """

    def __init__(
        self,
        reader: BaseMSIReader,
        store_path: Callable[[], Path],
        *,
        write_mobility_table: bool = True,
        mobility_heatmap: bool = True,
        mobility_grid: bool = False,
        mobility_bins: int = MOBILITY_CHANNELS,
        mobility_min: Optional[float] = None,
        mobility_max: Optional[float] = None,
        write_msms_table: bool = True,
    ) -> None:
        """Build the sibling tables of one conversion.

        Every argument is settled when the converter is constructed;
        everything a conversion decides later arrives per call in a
        :class:`SiblingContext`.

        Args:
            reader: The MSI reader this conversion reads. Only its
                mobility and fragmentation capabilities are asked about
                here (``has_ion_mobility``, ``has_shared_mobility_axis``,
                ``get_mobility_axis``, ``has_precursor_spectra``,
                ``has_frame_scans``); the builders read the rest. A reader
                here is whatever satisfies that interface, not necessarily
                a :class:`~thyra.core.base_reader.BaseMSIReader` subclass.
            store_path: The store this conversion is writing. A table's
                scratch is made in its parent, so the memmaps and the store
                they are written into share a filesystem. Called when a
                scratch is made rather than read here, because the output
                path a converter was *given* is not always the one it
                writes: ``prepare_zarr_output_path`` shortens it on
                Windows, and a caller that stands a converter up itself
                assigns the shortened path afterwards. Reading it once at
                construction would put the scratch beside a path nothing
                was written to.
            write_mobility_table: Whether a mobility-resolved sibling is
                wanted at all.
            mobility_heatmap: Whether the mass-mobility heatmap block is
                wanted.
            mobility_grid: Whether to bin a per-pixel mobility point cloud
                onto a common grid (opt in; a source that already shares a
                feature axis needs no grid and ignores this).
            mobility_bins: Channels that grid divides the mobility range
                into.
            mobility_min: Lower edge of the grid, or ``None`` for the
                smallest value the source's axis holds.
            mobility_max: Upper edge, or ``None`` for the largest.
            write_msms_table: Whether a demultiplexed MS/MS sibling is
                wanted.
        """
        self.reader = reader
        self.store_path = store_path

        # The mobility-resolved sibling (see mobility_table.py): whether to
        # write one, and -- once a finalize step has decided for its slice --
        # the element key it gets, so the MSI table's uns can name it.
        self._write_mobility_table = bool(write_mobility_table)
        self._mobility_table_key: Optional[str] = None
        # The common mobility grid (see resampling/mobility_grid.py): the
        # second way to fill the same sibling, for a source whose pixels
        # each carry their own mobility values. Resolved once by
        # _plan_mobility_table, so the MSI table's metadata block and the
        # sibling describe the same grid.
        self._mobility_grid_enabled = bool(mobility_grid)
        self._mobility_bins = int(mobility_bins)
        self._mobility_bounds = (mobility_min, mobility_max)
        self._mobility_grid: Optional[MobilityGrid] = None
        self._mobility_grid_resolved = False
        # The grid's discovery pass for the slice being written, when the
        # converter ran it fused with the heatmap's pass (see
        # prepare_scans); consumed by attach().
        self._grid_discovery: Any = None
        # The MS/MS table's accumulator when the converter fed it from the
        # summed table's own passes (see fused_passes.py); consumed by
        # attach() like the grid's discovery.
        self._msms_accumulator: Any = None
        # Whether the sibling sinks were fed from the summed table's passes
        # already, so prepare_scans has nothing left to scan; and the table
        # key the siblings were last planned for, so the route that plans in
        # fused_passes() and the route that plans per slice can both simply
        # call plan() (see there).
        self._sibling_scans_done = False
        self._planned_key: Optional[str] = None
        # Scratch directories holding the memmapped matrices of every table
        # until it is written; released by release_scratch.
        self._table_scratch: List[Tuple[Any, Path]] = []
        # The mass-mobility heatmap (see mobility_heatmap.py): built once
        # per conversion, on first demand, and shared by every uns block
        # that asks for it. ``_built`` distinguishes "not yet" from
        # "tried, nothing to write".
        self._mobility_heatmap_enabled = bool(mobility_heatmap)
        self._mobility_heatmap_block: Optional[Dict[str, Any]] = None
        self._mobility_heatmap_built = False
        # The demultiplexed MS/MS sibling (see msms_table.py): on by default,
        # and -- once a finalize step has decided for its slice -- the element
        # key it gets, so the MSI table's uns can name it.
        self._write_msms_table = bool(write_msms_table)
        self._msms_table_key: Optional[str] = None

    # -- What the rest of the conversion reads ---------------------------

    @property
    def mobility_table_key(self) -> Optional[str]:
        """Element key of this slice's mobility sibling, or ``None``.

        ``None`` both before :meth:`plan` has decided and after a builder
        declined: a key here is a promise that the element exists, which
        is what the summed table's ``uns`` records (issue #343).
        """
        return self._mobility_table_key

    @property
    def msms_table_key(self) -> Optional[str]:
        """Element key of this slice's demultiplexed sibling, or ``None``."""
        return self._msms_table_key

    @property
    def mobility_grid(self) -> Optional[MobilityGrid]:
        """The resolved common mobility grid, or ``None``.

        Resolved once per conversion, by :meth:`plan`, so the summed
        table's metadata block and the sibling describe the same grid.
        """
        return self._mobility_grid

    # -- The heatmap -----------------------------------------------------

    def heatmap(self, ctx: SiblingContext) -> Optional[Dict[str, Any]]:
        """Build the heatmap the first time it is asked for; cache the result.

        Needs the common mass axis, so it can only run after
        ``_initialize_conversion``. A failure is logged and leaves the
        summed table untouched: the block is additive, and a store
        without it is still complete.
        """
        if self._mobility_heatmap_built:
            return self._mobility_heatmap_block
        self._mobility_heatmap_built = True
        if not self._mobility_heatmap_enabled:
            return None
        try:
            if not self.reader.has_ion_mobility:
                return None
        except Exception as e:  # pragma: no cover - reader-defined
            logger.warning("Could not inspect the mobility axis: %s", str(e))
            return None
        if ctx.mass_axis is None:
            logger.warning(
                "No mass-mobility heatmap: the common mass axis is not built yet"
            )
            return None
        from .mobility_heatmap import build_mobility_heatmap

        try:
            self._mobility_heatmap_block = build_mobility_heatmap(
                self.reader,
                ctx.mass_axis,
                n_spectra=ctx.n_spectra(),
                linearisation=ctx.linearisation,
            )
        except Exception as e:
            logger.error("Could not build the mass-mobility heatmap: %s", str(e))
            self._mobility_heatmap_block = None
        return self._mobility_heatmap_block

    # -- Planning --------------------------------------------------------

    def plan(self, ctx: SiblingContext, table_key: str) -> None:
        """Decide which siblings of ``table_key`` this conversion writes.

        Called by a finalize step before the summed table's ``uns`` is
        composed, so the two agree about what exists, and by
        :meth:`fused_passes`, whose sinks belong to the siblings it plans.

        Idempotent per key. The fused route plans the single table it
        serves and the finalize step then asks again for that same key;
        a multi-plane conversion, which the fused route never serves,
        asks once per plane and gets a decision per plane. The flag this
        replaced (``_siblings_planned``) said only *that* something had
        been planned, and every caller had to remember to test it before
        planning -- an obligation with nothing enforcing it, on a path
        where planning twice re-names a sibling a builder has already
        declined.
        """
        if self._planned_key == table_key:
            return
        self._mobility_table_key = self._plan_mobility_table(ctx, table_key)
        self._msms_table_key = self._plan_msms_table(ctx, table_key)
        self._planned_key = table_key

    def _plan_mobility_table(
        self, ctx: SiblingContext, table_key: str
    ) -> Optional[str]:
        """The key of the mobility table this slice gets, or ``None``.

        Decided before the MSI table's ``uns`` is built so the two agree,
        which for the grid route also means the grid itself is resolved
        here: the summed table's metadata block names it, and the sibling
        must be binned onto the very grid that was named.

        Two mechanisms can fill the same key -- a shared feature axis
        (nothing to decide) or a common grid (opt in) -- and a source that
        allows neither gets no table, said by name rather than silently.
        """
        if not self._write_mobility_table:
            return None
        try:
            if not self.reader.has_ion_mobility:
                return None
            shared = bool(self.reader.has_shared_mobility_axis)
        except Exception as e:  # pragma: no cover - reader-defined
            logger.warning("Could not inspect the mobility axis: %s", str(e))
            return None
        from .mobility_table import mobility_table_key

        if shared:
            return mobility_table_key(table_key)
        if not self._mobility_grid_enabled:
            logger.info(
                "No mobility-resolved table: %s carries mobility per pixel "
                "rather than as a shared feature axis. Pass --mobility-grid "
                "to bin it onto a common mobility grid.",
                type(self.reader).__name__,
            )
            return None
        if not self._mobility_grid_resolved:
            self._mobility_grid_resolved = True
            self._mobility_grid = self._resolve_mobility_grid(ctx)
        if self._mobility_grid is None:
            return None
        return mobility_table_key(table_key)

    def _resolve_mobility_grid(self, ctx: SiblingContext) -> Optional[MobilityGrid]:
        """The common mobility grid this conversion bins onto, or ``None``.

        Bounds come from the axis values unless the caller overrode them:
        the per-scan 1/K0 of a real file overhangs its declared
        acquisition range, and the mass-mobility heatmap already bins over
        the values, so anything else breaks the index-for-index mapping
        between the two. Every refusal is said by name.
        """
        from .mobility_table import (
            MAX_GRID_VAR_ENTRIES,
            grid_refusal,
            grid_var_bound,
            mobility_grid_range,
        )

        if ctx.mass_axis is None:
            logger.warning(
                "No mobility-resolved table: the common mass axis is not "
                "built yet, so the grid's m/z bins are unknown"
            )
            return None
        try:
            measured = mobility_grid_range(self.reader)
        except Exception as e:  # pragma: no cover - reader-defined
            logger.warning("Could not read the mobility axis values: %s", str(e))
            return None
        lower, upper = self._mobility_bounds
        if measured is None and (lower is None or upper is None):
            logger.warning(
                "No mobility-resolved table: the source's mobility axis "
                "carries no per-scan values to bin over (a reader opened "
                "without its vendor library cannot supply them). Give "
                "--mobility-min and --mobility-max to bin over a stated range."
            )
            return None
        span = measured or (0.0, 0.0)
        try:
            grid = build_mobility_grid(
                span[0] if lower is None else float(lower),
                span[1] if upper is None else float(upper),
                self._mobility_bins,
            )
        except ValueError as e:
            logger.warning("No mobility-resolved table: %s", str(e))
            return None
        refusal = grid_refusal(self.reader, ctx.mass_axis, grid)
        if refusal is not None:
            logger.warning("No mobility-resolved table: %s", refusal)
            return None
        unit = None
        axis = self.reader.get_mobility_axis()
        if axis is not None and axis.unit_name:
            unit = str(axis.unit_name)
        logger.info(
            "Mobility grid: %d %s channels over [%.5f, %.5f]%s",
            grid.n_channels,
            grid.law,
            grid.lower,
            grid.upper,
            "" if unit is None else f" {unit}",
        )
        report_channel_width(grid)
        bound = grid_var_bound(ctx.mass_axis, grid)
        if bound > MAX_GRID_VAR_ENTRIES:
            # A bound above the ceiling settles nothing -- real occupancy
            # runs an order of magnitude below it -- so it is said and the
            # source is read; the count decides (see var_ceiling_refusal).
            logger.info(
                "The mobility grid spans %s (m/z bin, channel) pairs, above "
                "the var ceiling of %s. Most of them will be empty; the "
                "table is refused only if the pairs that carry signal pass "
                "the ceiling too.",
                f"{bound:,}",
                f"{MAX_GRID_VAR_ENTRIES:,}",
            )
        return grid

    def _plan_msms_table(self, ctx: SiblingContext, table_key: str) -> Optional[str]:
        """The key of the demultiplexed MS/MS table this slice gets, or ``None``.

        Decided before the MSI table's ``uns`` is built so the two agree.
        Every refusal is said by name: writing no table is the right
        answer for an acquisition whose precursors cannot be told apart,
        but silently writing none is not.
        """
        if not self._write_msms_table:
            return None
        from .msms_table import demultiplex_refusal, msms_table_key

        # Order matters: the schedule-specific refusal is the one a user can
        # act on ("a single precursor", "the windows carry no scan range"),
        # while the capability check can only name the reader. Ask for the
        # specific reason first and fall back to the generic one, or an
        # acquisition whose precursors are merely uninteresting gets told
        # its reader is incapable.
        refusal = demultiplex_refusal(ctx.fragmentation())
        if refusal is not None:
            logger.info("No demultiplexed MS/MS table: %s", refusal)
            return None
        if not self.reader.has_precursor_spectra:
            logger.info(
                "No demultiplexed MS/MS table: %s cannot separate the "
                "precursors of a pixel",
                type(self.reader).__name__,
            )
            return None
        # The fragment axis is the summed table's axis whatever that axis
        # is (design decision D6): resampled, it is the grid the user chose
        # for the whole store; raw, it is the union of the very fragment
        # m/z values the split re-reads, so the mapping is exact either way.
        return msms_table_key(table_key)

    # -- The raw pass ----------------------------------------------------

    def prepare_scans(
        self, ctx: SiblingContext, obs: pd.DataFrame, z_value: Optional[int] = None
    ) -> None:
        """Run the raw mobility pass once for everything that needs it.

        Called by a finalize step right after the sibling tables are
        planned and before the summed table's ``uns`` is built, with the
        ``obs`` the siblings will mirror. When a grid table is planned its
        discovery pass -- which occupied cells there are, and how many
        rows each holds -- is fused into the heatmap's pass, so the two
        share one read *and* one mapping of every point onto the mass
        axis; the mapping is the larger cost of the two. The heatmap is
        built here once and cached for every ``uns`` block that asks; a
        route that never calls this still gets it from :meth:`heatmap` on
        first demand.

        A no-op when the sinks were already fed from the summed table's
        own passes (a reader that hands its frames over as records; see
        ``fused_passes.py``).
        """
        if self._sibling_scans_done:
            return
        self._grid_discovery = None
        try:
            if not self.reader.has_ion_mobility:
                return
        except Exception as e:  # pragma: no cover - reader-defined
            logger.warning("Could not inspect the mobility axis: %s", str(e))
            return
        if ctx.mass_axis is None:
            return
        from .mobility_heatmap import finish_mobility_heatmap, scan_mobility

        heatmap = self._pending_heatmap(ctx)
        discovery = self._pending_grid_discovery(ctx, obs, z_value)
        sinks = [sink for sink in (heatmap, discovery) if sink is not None]
        if not sinks:
            return
        try:
            scan_mobility(
                self.reader,
                ctx.mass_axis,
                *sinks,
                n_spectra=ctx.n_spectra(),
                description=(
                    "Mobility heatmap + grid"
                    if discovery is not None
                    else "Mobility heatmap"
                ),
                linearisation=ctx.linearisation,
            )
        except Exception as e:
            logger.error("Could not scan the mobility spectra: %s", str(e))
            if heatmap is not None:
                self._mobility_heatmap_built = True
                self._mobility_heatmap_block = None
            return
        if heatmap is not None:
            self._mobility_heatmap_built = True
            self._mobility_heatmap_block = finish_mobility_heatmap(heatmap)
        if discovery is not None:
            discovery.finish()
            self._grid_discovery = discovery

    def _pending_heatmap(self, ctx: SiblingContext) -> Any:
        """An empty heatmap accumulator, when one is wanted and not yet built."""
        if not self._mobility_heatmap_enabled or self._mobility_heatmap_built:
            return None
        from .mobility_heatmap import prepare_mobility_heatmap

        return prepare_mobility_heatmap(self.reader, ctx.mass_axis)

    def _pending_grid_discovery(
        self, ctx: SiblingContext, obs: pd.DataFrame, z_value: Optional[int]
    ) -> Any:
        """The grid's discovery accumulator, when a grid table is planned."""
        if self._mobility_table_key is None or self._mobility_grid is None:
            return None
        from .mobility_table import GridDiscovery, row_lookup

        try:
            return GridDiscovery(
                ctx.mass_axis,
                self._mobility_grid,
                row_lookup(obs, z_value, None),
                int(len(obs)),
            )
        except MemoryError as e:
            logger.warning("No mobility-resolved table: %s", str(e))
            return None

    def fused_passes(self, ctx: SiblingContext, table_key: str) -> Any:
        """The sibling sinks to feed from the summed table's own passes, or ``None``.

        Only for a reader that hands its frames over as records
        (:attr:`~thyra.core.base_reader.BaseMSIReader.has_frame_scans`);
        plans the siblings of ``table_key`` first, since the sinks are
        theirs. ``None`` when nothing wants the frames, in which case the
        passes read the summed spectra as they always did and
        :meth:`prepare_scans` scans on its own later.
        """
        if not self.reader.has_frame_scans:
            return None
        if ctx.mass_axis is None or ctx.dimensions is None:
            return None
        self.plan(ctx, table_key)
        from .fused_passes import SiblingPasses
        from .mobility_table import GridDiscovery
        from .msms_table import new_msms_accumulator

        n_x, n_y, n_z = ctx.dimensions
        n_grid = int(n_x * n_y * n_z)
        heatmap = self._pending_heatmap(ctx)
        discovery = None
        if self._mobility_table_key is not None and self._mobility_grid is not None:
            try:
                # Rows are handed to the sinks by the passes themselves,
                # so the lookup a standalone pass would use is not needed.
                discovery = GridDiscovery(
                    ctx.mass_axis,
                    self._mobility_grid,
                    lambda coords: None,
                    n_grid,
                )
            except MemoryError as e:
                logger.warning("No mobility-resolved table: %s", str(e))
        msms = None
        if self._msms_table_key is not None:
            msms = new_msms_accumulator(
                self.reader,
                ctx.mass_axis,
                n_grid,
                ctx.linearisation,
            )
        if heatmap is None and discovery is None and msms is None:
            return None
        self._sibling_scans_done = True
        return SiblingPasses(
            ctx.mass_axis,
            heatmap=heatmap,
            discovery=discovery,
            msms=msms,
            linearisation=ctx.linearisation,
        )

    def take_fused_results(self, passes: Any) -> None:
        """Keep what the fused passes built for the finalize step to write."""
        if passes.heatmap_wanted:
            self._mobility_heatmap_built = True
            self._mobility_heatmap_block = passes.heatmap_block
        self._grid_discovery = passes.discovery
        self._msms_accumulator = passes.msms

    # -- Scratch ---------------------------------------------------------

    def _new_scratch(self, prefix: str) -> Path:
        """A scratch directory for one sibling's memmaps, next to the output."""
        from .csc_assembly import scratch_directory

        return scratch_directory(f".thyra_{prefix}_", parent=self.store_path().parent)

    def register_scratch(self, prefix: str, assembly: Any) -> Path:
        """A scratch directory for ``assembly``, released with the others once written."""
        scratch = self._new_scratch(prefix)
        self._table_scratch.append((assembly, scratch))
        return scratch

    def release_scratch(self, tables: Optional[Dict[str, Any]] = None) -> None:
        """Drop every table's memmaps and remove their scratch directories.

        Called once the store is written. ``tables`` is the mapping that
        still holds the tables; it is emptied first, because a mapped file
        cannot be deleted on Windows and the AnnData is what keeps it
        mapped. Idempotent, so the ``finally`` of ``convert`` can call it
        too.
        """
        from .csc_assembly import remove_scratch

        if not self._table_scratch:
            return
        if tables is not None:
            tables.clear()
        pending = self._table_scratch
        self._table_scratch = []
        for assembly, path in pending:
            if assembly is not None:
                assembly.release()
            remove_scratch(path)

    # -- Building --------------------------------------------------------

    def attach(
        self,
        ctx: SiblingContext,
        tables: Dict[str, Any],
        table_key: str,
        region_key: str,
        obs: pd.DataFrame,
        z_value: Optional[int] = None,
    ) -> None:
        """Build the sibling tables of ``table_key`` and add them to ``tables``.

        No-op for a sibling :meth:`plan` did not name for this slice.

        A sibling that *was* named is not additive, and this used to say
        it was. The summed table's ``uns`` is built before the siblings
        are, and it carries their keys, so a failure swallowed here leaves
        a store whose summed table points at an element nobody wrote. A
        failure therefore propagates now (issue #280). A builder may still
        decline by returning ``None`` -- a decision, not a failure -- and
        the name is taken back out when it does (issue #343).
        """
        if ctx.mass_axis is None:
            return
        if self._mobility_table_key is None and self._msms_table_key is None:
            return
        # The siblings carry the same provenance as the summed table,
        # minus the heatmap: that block is the summed table's navigator
        # over the very data the siblings hold resolved or split.
        sibling_uns = ctx.build_uns()
        sibling_uns.pop("mobility_heatmap", None)
        summed = tables.get(table_key)
        declined: List[str] = []
        if self._mobility_table_key is not None:
            table = self._build_mobility_sibling(
                ctx, obs, table_key, region_key, dict(sibling_uns), z_value
            )
            if table is None:
                declined.append("ion_mobility")
                self._mobility_table_key = None
            else:
                if self._mobility_grid is not None:
                    self._record_mobility_marginal(table, summed, table_key)
                tables[self._mobility_table_key] = table
        if self._msms_table_key is not None:
            table = self._build_msms_sibling(
                ctx, obs, table_key, region_key, dict(sibling_uns), z_value
            )
            if table is None:
                declined.append("fragmentation")
                self._msms_table_key = None
            else:
                self._record_demultiplexed_current(table, summed, table_key)
                tables[self._msms_table_key] = table
        if declined:
            # Every table written for this slice, not just the summed one:
            # each sibling's uns was taken from the summed table's before
            # either builder ran, so each carries the same stale pointer.
            ctx.unname_declined(
                tables,
                (table_key, self._mobility_table_key, self._msms_table_key),
                declined,
            )

    def _build_mobility_sibling(
        self,
        ctx: SiblingContext,
        obs: pd.DataFrame,
        table_key: str,
        region_key: str,
        uns: Dict[str, Any],
        z_value: Optional[int],
    ) -> Optional[Any]:
        """The mobility sibling of one slice, on a scratch directory of its own."""
        from .mobility_table import build_mobility_table

        discovery, self._grid_discovery = self._grid_discovery, None
        # The fused passes allocate and register their own scratch; a
        # discovery without one is scattered by the builder on a new one.
        scratch = None if discovery is None else discovery.scratch
        if scratch is None:
            scratch = self._new_scratch("mobility")
            self._table_scratch.append(
                (None if discovery is None else discovery.assembly, scratch)
            )
        # Deliberately uncaught (issue #280). By the time this runs, the
        # summed table's uns has already named this table, twice: in
        # ``mobility_axis["resolved_table"]`` and in the versioned block at
        # ``msi_metadata.ms_analysis.ion_mobility.resolved_table``. Turning
        # a failure into a log line does not leave "a store without the
        # sibling, still complete". It leaves a store pointing at an
        # element that was never written (issue #343 has the measurement),
        # and
        # nothing downstream checks that the pointer resolves: no validator
        # rule, no test. A ConversionRefused from the builder is a sentence
        # addressed to whoever ran the conversion (csc_assembly's
        # two-passes-disagree, mobility_table's off-axis pair); anything
        # else is an invariant break whose traceback is the explanation.
        # Neither is served by being swallowed here. Nothing is on disk yet
        # at this point, so propagating costs a store that would have lied,
        # not a store that was written.
        return build_mobility_table(
            self.reader,
            obs,
            ctx.mass_axis,
            table_key,
            region_key,
            uns,
            z_value=z_value,
            grid=self._mobility_grid,
            discovery=discovery,
            scratch=scratch,
            linearisation=ctx.linearisation,
        )

    def _build_msms_sibling(
        self,
        ctx: SiblingContext,
        obs: pd.DataFrame,
        table_key: str,
        region_key: str,
        uns: Dict[str, Any],
        z_value: Optional[int],
    ) -> Optional[Any]:
        """The demultiplexed sibling of one slice, on a scratch directory of its own."""
        from .msms_table import build_msms_table

        accumulator, self._msms_accumulator = self._msms_accumulator, None
        scratch = None if accumulator is None else accumulator.scratch
        if scratch is None:
            scratch = self._new_scratch("msms")
            self._table_scratch.append((None, scratch))
        # Deliberately uncaught, for the reason given in
        # :meth:`_build_mobility_sibling`: the summed table's uns already
        # names this table, in ``msms_schedule["resolved_table"]`` and at
        # ``msi_metadata.ms_analysis.fragmentation.resolved_table``, which
        # schema 0.5.0 added for exactly that purpose (issue #280).
        return build_msms_table(
            self.reader,
            obs,
            ctx.mass_axis,
            table_key,
            region_key,
            uns,
            z_value=z_value,
            scratch=scratch,
            accumulator=accumulator,
            linearisation=ctx.linearisation,
        )

    # -- What the two tables say about each other ------------------------

    @staticmethod
    def _record_mobility_marginal(table: Any, summed: Any, summed_key: str) -> None:
        """Say how far the grid table's marginal is from the summed table.

        Summing a grid table's channels within one m/z bin must reproduce
        that bin's column of the summed table: both are the same points,
        binned the same way, differing only in whether mobility was kept.
        Under ``--tdf-spectrum scan_sum`` that holds exactly; under the
        vendor centroid it cannot, because the centroid is a peak-picked
        spectrum over the same ramp and keeps only the current inside the
        peaks it picks (87-96% on measured acquisitions) while the grid
        reads raw scans. A store whose two tables disagree
        must say by how much rather than leave a reader to find it by
        subtraction.

        The comparison is per pixel: the grid table's row sums against
        the summed table's, which is one bounded pass over each memmap.
        The per-cell deviation the in-memory converters used to record
        needed the marginal and its difference from the summed table
        materialised, each the size of the summed table, and went with
        them (design decision D11); the current ratio is the same number
        it always was.
        """
        try:
            block = _current_ratio_block(table, summed_key, summed)
        except MALFORMED_METADATA as e:  # pragma: no cover - defensive
            # Two row sums over two memmaps (issue #280): a released memmap
            # or a matrix that is not the shape it should be raises one of
            # these. A comparison is the whole job here, so anything wider
            # would hide the defect that broke it in a debug line.
            logger.debug("Could not compare the mobility marginal: %s", str(e))
            return
        if block is None:
            return
        table.uns["mobility_marginal"] = block
        exact = abs(float(block["current_ratio"]) - 1.0) <= _MARGINAL_TOLERANCE
        log = logger.info if exact else logger.warning
        log(
            "The mobility grid table holds %.4fx the summed table's ion "
            "current (per pixel %.4f to %.4f). They agree exactly only "
            "under --tdf-spectrum scan_sum; the vendor centroid is a "
            "peak-picked spectrum over the same scans and keeps less of "
            "the current.",
            block["current_ratio"],
            block["current_ratio_pixel_min"],
            block["current_ratio_pixel_max"],
        )

    @staticmethod
    def _record_demultiplexed_current(table: Any, summed: Any, summed_key: str) -> None:
        """Say how much of the summed table's ion current the split holds.

        The two tables agree exactly under ``--tdf-spectrum scan_sum``,
        which writing this table now selects; under an explicit
        ``vendor_centroid`` they do not, because the vendor peak picker
        drops the index bins it assigns to no peak while the split reads
        raw scans, so the
        demultiplexed table holds *more*. That is not a defect, but a
        store whose two tables disagree must say so rather than leave a
        reader to find it by subtraction.
        """
        try:
            block = _current_ratio_block(table, summed_key, summed)
        except MALFORMED_METADATA as e:  # pragma: no cover - defensive
            # Same two row sums, same reasoning as _record_mobility_marginal.
            logger.debug("Could not compare the demultiplexed current: %s", str(e))
            return
        if block is None:
            return
        table.uns["demultiplexed_current"] = block
        if abs(float(block["current_ratio"]) - 1.0) > _MARGINAL_TOLERANCE:
            # WARNING, not INFO. It is the same disagreement between the
            # same two tables that the mobility grid reports at WARNING,
            # and a reader who has to know about one has to know about
            # the other (issue #253).
            logger.warning(
                "The demultiplexed table holds %.4fx the summed table's ion "
                "current (per pixel %.4f to %.4f). Above 1 under an explicit "
                "--tdf-spectrum vendor_centroid, which discards counts the "
                "raw scans keep; below 1 when the schedule's windows do not "
                "cover every scan that carries current.",
                block["current_ratio"],
                block["current_ratio_pixel_min"],
                block["current_ratio_pixel_max"],
            )


__all__ = ["SiblingContext", "SiblingTables"]
