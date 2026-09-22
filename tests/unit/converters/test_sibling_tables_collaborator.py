"""The sibling tables decided, and declined, without a conversion (#353).

Every test here stands :class:`~thyra.converters.spatialdata.sibling_tables.SiblingTables`
up on a reader stub and an accessor for the store it would be writing -- no
converter, no store, no source read -- and asks it the questions a
conversion asks it.

That is the point of the module's existence. The same decisions used to be
seventeen private methods on the base converter, reachable only by driving
``convert()`` end to end, and issue #343 is what that cost: a sibling
builder that declines by returning ``None`` left its name behind in the
summed table's ``uns``, so a store was written, reported as a success, and
pointed at an element nobody wrote. :class:`TestADeclinedSiblingIsUnnamed`
is that scenario in isolation -- it is the test the issue asked for, and
``tests/unit/converters/test_internal_failures_are_not_warnings.py``
still pins the same behaviour through a whole conversion.

The ``uns`` half is not stubbed: the assembler that composes the block and
the assembler that takes a pointer back out are the real one, built on the
same stub reader, because the defect lived in the seam between the two
collaborators rather than inside either.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("spatialdata")

from tests.unit.converters.test_mobility_grid import (  # noqa: E402
    MASS_AXIS,
    PIXELS,
    GridStubReader,
)
from thyra.converters.spatialdata.sibling_tables import (  # noqa: E402
    SiblingContext,
    SiblingTables,
)
from thyra.converters.spatialdata.uns_assembler import (  # noqa: E402
    UnsAssembler,
    UnsContext,
)
from thyra.core.base_converter import PixelSizeSource  # noqa: E402
from thyra.core.msms import FragmentationSchedule, IsolationWindow  # noqa: E402
from thyra.metadata.schema import MSI_METADATA_UNS_KEY  # noqa: E402

_TABLE_KEY = "stub_z0"
_REGION_KEY = "stub_z0_pixels"


class _SharedButEmpty(GridStubReader):
    """Reports a shared feature axis it cannot list.

    The reachable route into issue #343 and the one it was measured on:
    ``_plan_mobility_table`` names the table on the shared-axis route with
    no further check, and ``build_mobility_table`` declines whenever the
    listing comes back empty.
    """

    @property
    def has_shared_mobility_axis(self) -> bool:
        return True

    def get_shared_mobility_features(self):
        return None


def _schedule(n_windows: int) -> FragmentationSchedule:
    """A separable PASEF schedule over ``n_windows`` disjoint scan ranges."""
    windows = tuple(
        IsolationWindow(
            300.0 + 100.0 * i,
            0.5,
            0.5,
            34.3,
            scan_begin=20 + 50 * i,
            scan_end=60 + 50 * i,
        )
        for i in range(n_windows)
    )
    return FragmentationSchedule(ms_level=2, windows=windows, source="stub")


def _fragmenting(reader_cls, n_windows: int = 2):
    """``reader_cls`` with a schedule of ``n_windows`` it can also split."""
    schedule = _schedule(n_windows)

    class _Reader(reader_cls):  # type: ignore[valid-type, misc]
        @property
        def has_fragmentation(self) -> bool:
            return True

        @property
        def has_precursor_spectra(self) -> bool:
            return True

        def get_fragmentation(self):
            return schedule

        def iter_precursor_spectra(self):
            for p, (x, y) in enumerate(PIXELS):
                for window in range(n_windows):
                    yield (
                        (x, y, 0),
                        window,
                        MASS_AXIS[: 2 + window].copy(),
                        np.arange(2 + window, dtype=np.float64) + p + 1.0,
                    )

    return _Reader()


def _obs() -> pd.DataFrame:
    """The ``obs`` of the summed table the siblings mirror."""
    return pd.DataFrame(
        {
            "x": [x for x, _ in PIXELS],
            "y": [y for _, y in PIXELS],
            "instance_id": [str(i) for i in range(len(PIXELS))],
        }
    )


def _context(siblings: SiblingTables, reader: Any) -> SiblingContext:
    """The context a converter would pack, with the real ``uns`` assembler.

    ``build_uns`` and ``unname_declined`` are the assembler's own, so the
    pointers a decline has to remove are the ones a conversion writes.
    """
    assembler = UnsAssembler(
        reader,
        pixel_size_xy=lambda: (10.0, 10.0),
        pixel_size_detection_info={"source_format": "stub"},
        resampling_config=None,
    )

    def build_uns() -> Dict[str, Any]:
        return assembler.build(
            UnsContext(
                mobility_table_key=siblings.mobility_table_key,
                msms_table_key=siblings.msms_table_key,
                mobility_grid=siblings.mobility_grid,
                resolved_resampling_plan=None,
                region_info=None,
                pixel_size_source=PixelSizeSource.DEFAULT,
                mobility_heatmap=lambda: None,
            )
        )

    return SiblingContext(
        mass_axis=MASS_AXIS.copy(),
        linearisation=None,
        dimensions=(4, 4, 1),
        n_spectra=lambda: len(PIXELS),
        build_uns=build_uns,
        fragmentation=assembler.fragmentation,
        unname_declined=assembler.unname_declined_siblings,
    )


def _pointers(uns: Dict[str, Any]) -> Dict[str, str]:
    """Every place in one table's ``uns`` that names a sibling element."""
    found: Dict[str, str] = {}
    for key in ("mobility_axis", "msms_schedule"):
        plain = uns.get(key)
        if isinstance(plain, dict) and "resolved_table" in plain:
            found[f"{key}.resolved_table"] = plain["resolved_table"]
    sections = uns.get(MSI_METADATA_UNS_KEY, {}).get("ms_analysis", {})
    for section in ("ion_mobility", "fragmentation"):
        block = sections.get(section)
        if isinstance(block, dict) and "resolved_table" in block:
            found[f"{section}.resolved_table"] = block["resolved_table"]
    return found


def _summed_table(uns: Dict[str, Any]) -> Any:
    """A stand-in for the summed table: its ``uns`` and a matrix to sum."""

    class _Table:
        def __init__(self) -> None:
            self.uns = uns
            self.X = np.ones((len(PIXELS), MASS_AXIS.size), dtype=np.float64)

    return _Table()


class TestPlanning:
    """Which siblings a source gets, asked of the source alone."""

    def test_a_shared_feature_axis_names_the_mobility_table(self, tmp_path: Path):
        reader = _SharedButEmpty()
        siblings = SiblingTables(reader, lambda: tmp_path / "out.zarr")

        siblings.plan(_context(siblings, reader), _TABLE_KEY)

        assert siblings.mobility_table_key == f"{_TABLE_KEY}_mobility"
        assert siblings.msms_table_key is None

    def test_a_per_pixel_cloud_needs_the_grid_and_says_so(
        self, tmp_path: Path, thyra_logs
    ):
        """``GridStubReader`` carries mobility per pixel, not as a feature axis."""
        reader = GridStubReader()
        siblings = SiblingTables(reader, lambda: tmp_path / "out.zarr")

        with thyra_logs("thyra", logging.INFO) as records:
            siblings.plan(_context(siblings, reader), _TABLE_KEY)

        assert siblings.mobility_table_key is None
        assert any("--mobility-grid" in r.getMessage() for r in records)

    def test_the_grid_route_names_it_and_resolves_one_grid(self, tmp_path: Path):
        reader = GridStubReader()
        siblings = SiblingTables(
            reader, lambda: tmp_path / "out.zarr", mobility_grid=True
        )

        siblings.plan(_context(siblings, reader), _TABLE_KEY)

        assert siblings.mobility_table_key == f"{_TABLE_KEY}_mobility"
        assert siblings.mobility_grid is not None

    def test_turning_it_off_names_nothing(self, tmp_path: Path):
        reader = _SharedButEmpty()
        siblings = SiblingTables(
            reader, lambda: tmp_path / "out.zarr", write_mobility_table=False
        )

        siblings.plan(_context(siblings, reader), _TABLE_KEY)

        assert siblings.mobility_table_key is None

    def test_a_source_without_the_dimension_is_asked_nothing_else(self, tmp_path: Path):
        class _Flat(GridStubReader):
            @property
            def has_ion_mobility(self) -> bool:
                return False

        reader = _Flat()
        siblings = SiblingTables(
            reader, lambda: tmp_path / "out.zarr", mobility_grid=True
        )

        siblings.plan(_context(siblings, reader), _TABLE_KEY)

        assert siblings.mobility_table_key is None
        assert siblings.mobility_grid is None


def _incapable(schedule: FragmentationSchedule):
    """A reader that reports ``schedule`` and cannot split a pixel by it.

    Both refusals in ``_plan_msms_table`` apply to such a source, so which
    one it is told is decided by the order they are asked in -- and only
    by that. It is the ordinary shape for a fragmentation source Thyra
    has no precursor reader for.
    """

    class _Reader(GridStubReader):
        @property
        def has_fragmentation(self) -> bool:
            return True

        @property
        def has_precursor_spectra(self) -> bool:
            return False

        def get_fragmentation(self):
            return schedule

    return _Reader()


class TestTheMsmsRefusalOrder:
    """The ordering the #276 review called load-bearing, asked directly.

    ``_plan_msms_table`` puts ``demultiplex_refusal()`` before the
    reader's ``has_precursor_spectra`` check. Both fire on a source that
    states a schedule Thyra has no precursor reader for, so the order is
    the whole of what decides which sentence that source gets: reversed,
    an acquisition that isolates one precursor -- a fact about the method,
    which its own schedule states and a user can act on -- is told instead
    that its reader cannot separate precursors at all. Until this module
    existed the order could only be checked by reading it.
    """

    def test_a_single_precursor_is_told_about_its_schedule(
        self, tmp_path: Path, thyra_logs
    ):
        reader = _incapable(_schedule(1))
        siblings = SiblingTables(reader, lambda: tmp_path / "out.zarr")

        with thyra_logs("thyra", logging.INFO) as records:
            siblings.plan(_context(siblings, reader), _TABLE_KEY)

        assert siblings.msms_table_key is None
        messages = [r.getMessage() for r in records]
        assert any("isolates a single precursor" in m for m in messages)
        assert not any("cannot separate the precursors" in m for m in messages)

    def test_a_reader_that_cannot_split_is_told_that(self, tmp_path: Path, thyra_logs):
        """The fallback still fires for a schedule with nothing wrong with it."""
        reader = _incapable(_schedule(2))
        siblings = SiblingTables(reader, lambda: tmp_path / "out.zarr")

        with thyra_logs("thyra", logging.INFO) as records:
            siblings.plan(_context(siblings, reader), _TABLE_KEY)

        assert siblings.msms_table_key is None
        assert any("cannot separate the precursors" in r.getMessage() for r in records)

    def test_a_separable_schedule_names_the_table(self, tmp_path: Path):
        reader = _fragmenting(GridStubReader, n_windows=2)
        siblings = SiblingTables(reader, lambda: tmp_path / "out.zarr")

        siblings.plan(_context(siblings, reader), _TABLE_KEY)

        assert siblings.msms_table_key == f"{_TABLE_KEY}_msms"


class TestADeclinedSiblingIsUnnamed:
    """Issue #343, reproduced with no converter and no store.

    The mobility builder is named on the shared-axis route and then
    declines because the listing comes back empty. Nothing fails; the
    conversion would report success. What must not survive is the name.
    """

    def _attach(self, tmp_path: Path, reader: Any):
        siblings = SiblingTables(reader, lambda: tmp_path / "out.zarr")
        ctx = _context(siblings, reader)
        siblings.plan(ctx, _TABLE_KEY)
        assert siblings.mobility_table_key == f"{_TABLE_KEY}_mobility"

        # The summed table's uns, composed while the sibling was still
        # named -- which is the order a conversion writes it in.
        tables: Dict[str, Any] = {_TABLE_KEY: _summed_table(ctx.build_uns())}
        assert (
            _pointers(tables[_TABLE_KEY].uns)["mobility_axis.resolved_table"]
            == f"{_TABLE_KEY}_mobility"
        )

        siblings.attach(ctx, tables, _TABLE_KEY, _REGION_KEY, _obs())
        return siblings, tables

    def test_the_key_comes_back_out(self, tmp_path: Path):
        siblings, tables = self._attach(tmp_path, _SharedButEmpty())

        assert siblings.mobility_table_key is None
        assert set(tables) == {_TABLE_KEY}

    def test_the_summed_table_stops_naming_it(self, tmp_path: Path):
        _, tables = self._attach(tmp_path, _SharedButEmpty())

        assert _pointers(tables[_TABLE_KEY].uns) == {}

    def test_the_acquisition_fact_survives_the_unnaming(self, tmp_path: Path):
        """Only the pointer is withdrawn: the source still had the dimension."""
        _, tables = self._attach(tmp_path, _SharedButEmpty())

        uns = tables[_TABLE_KEY].uns
        assert uns[MSI_METADATA_UNS_KEY]["ms_analysis"]["ion_mobility"]["present"]
        assert "mobility_axis" in uns

    def test_a_surviving_sibling_does_not_name_the_declined_one(self, tmp_path: Path):
        """The siblings' ``uns`` is taken before either builder runs.

        So the MS/MS table inherits the summed table's pointer to the
        mobility table, and cleaning only the summed table would move the
        dangling pointer rather than remove it.
        """
        reader = _fragmenting(_SharedButEmpty, n_windows=2)
        siblings, tables = self._attach(tmp_path, reader)

        msms_key = f"{_TABLE_KEY}_msms"
        assert set(tables) == {_TABLE_KEY, msms_key}
        for key in (_TABLE_KEY, msms_key):
            pointers = _pointers(tables[key].uns)
            assert "ion_mobility.resolved_table" not in pointers
            assert "mobility_axis.resolved_table" not in pointers
            # The one that was written is still named, on both tables.
            assert pointers["fragmentation.resolved_table"] == msms_key

        siblings.release_scratch(tables)


class TestPlanningIsIdempotentPerKey:
    """What replaced the ``_siblings_planned`` flag every caller had to test.

    The fused route plans the one table it serves, and the finalize step
    then asks again for that same key; a multi-plane conversion, which the
    fused route never serves, asks once per plane. Keying the decision on
    the table it is about makes both right without a caller remembering
    which case it is in -- and planning twice after a decline would put
    the name the decline just removed straight back.
    """

    def test_a_second_ask_for_the_same_table_re_names_nothing(self, tmp_path: Path):
        reader = _SharedButEmpty()
        siblings = SiblingTables(reader, lambda: tmp_path / "out.zarr")
        ctx = _context(siblings, reader)

        siblings.plan(ctx, _TABLE_KEY)
        tables: Dict[str, Any] = {_TABLE_KEY: _summed_table(ctx.build_uns())}
        siblings.attach(ctx, tables, _TABLE_KEY, _REGION_KEY, _obs())
        assert siblings.mobility_table_key is None

        siblings.plan(ctx, _TABLE_KEY)

        assert siblings.mobility_table_key is None

    def test_the_next_slice_gets_its_own_decision(self, tmp_path: Path):
        reader = _SharedButEmpty()
        siblings = SiblingTables(reader, lambda: tmp_path / "out.zarr")
        ctx = _context(siblings, reader)

        siblings.plan(ctx, _TABLE_KEY)
        siblings.plan(ctx, "stub_z1")

        assert siblings.mobility_table_key == "stub_z1_mobility"


class TestTheScratchGoes:
    """The directories a table's memmaps live in, and the order they go in."""

    def test_a_registered_directory_is_made_next_to_the_output(self, tmp_path: Path):
        output = tmp_path / "store" / "out.zarr"
        output.parent.mkdir()
        siblings = SiblingTables(GridStubReader(), lambda: output)

        scratch = siblings.register_scratch("summed", None)

        assert scratch.parent == output.parent
        assert scratch.name.startswith(".thyra_summed_")
        assert scratch.is_dir()

    def test_the_store_is_read_when_the_scratch_is_made(self, tmp_path: Path):
        """The scratch follows the store, not the path handed over first.

        ``convert_msi`` shortens the output path before it builds a
        converter, but a caller standing one up itself assigns the
        shortened path afterwards -- and the scratch a table's memmaps
        live in has to land beside the store that is actually written,
        on the same filesystem as it.

        Reading the accessor once and keeping the answer puts the scratch
        beside ``first``, which this fails on.
        """
        first = tmp_path / "given"
        second = tmp_path / "written"
        first.mkdir()
        second.mkdir()
        output = first / "out.zarr"
        siblings = SiblingTables(GridStubReader(), lambda: output)

        output = second / "out.zarr"
        scratch = siblings.register_scratch("summed", None)

        assert scratch.parent == second
        assert not any(p.name.startswith(".thyra_") for p in first.iterdir())

    def test_releasing_removes_it(self, tmp_path: Path):
        siblings = SiblingTables(GridStubReader(), lambda: tmp_path / "out.zarr")
        scratch = siblings.register_scratch("mobility", None)

        siblings.release_scratch()

        assert not scratch.exists()

    def test_the_tables_are_dropped_before_anything_is_unlinked(self, tmp_path: Path):
        """Windows will not delete a mapped file, and the table is the map."""
        siblings = SiblingTables(GridStubReader(), lambda: tmp_path / "out.zarr")
        scratch = siblings.register_scratch("summed", None)
        tables = {_TABLE_KEY: object()}

        siblings.release_scratch(tables)

        assert tables == {}
        assert not scratch.exists()

    def test_releasing_twice_is_harmless(self, tmp_path: Path):
        """``convert``'s ``finally`` runs after ``_save_output`` already did."""
        siblings = SiblingTables(GridStubReader(), lambda: tmp_path / "out.zarr")
        siblings.register_scratch("summed", None)

        siblings.release_scratch()
        siblings.release_scratch()
