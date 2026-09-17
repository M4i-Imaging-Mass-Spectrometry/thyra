"""A failure inside Thyra must not come out as a warning line (issue #280).

The converter carried 28 ``except Exception`` catches. Most of them sit at
a boundary -- a reader call, a TIFF decode, a zarr write -- and keep their
breadth, because what surfaces there was raised by a vendor SDK, an XML
parser or sqlite, and no narrower tuple can be written honestly. The ones
these tests pin are the other kind: Thyra wrapping Thyra, where a failure
is a defect and the broad catch turned it into a log line.

Two halves, and both matter:

* the narrowed catches still tolerate what they were written for -- a
  source whose metadata is not the shape the block expected costs that one
  block, as before;
* an invariant break now propagates instead, because a store quietly
  missing a section is a worse outcome than a traceback.

The third class has a consequence that reached the store. A sibling table
is named in the summed table's ``uns`` *before* it is built, so swallowing
its failure left a store pointing at an element nobody wrote -- and nothing
downstream checks that the pointer resolves.

That naming has two ways to go wrong and this module pins both. A builder
that *fails* is issue #280, above. A builder that *declines* -- returns
``None``, which it is entitled to do -- reached the same dangling store
without anything failing, and is issue #343 (:class:`TestADeclinedSibling`).
"""

import logging
from pathlib import Path
from typing import Any, Dict

import pytest

pytest.importorskip("spatialdata")

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader  # noqa: E402
from thyra.errors import ConversionRefused  # noqa: E402

_DATASET_ID = "mock"


def _converter(tmp_path: Path, reader=None, **kwargs):
    from thyra.converters.spatialdata.streaming_converter import (
        StreamingSpatialDataConverter,
    )

    if reader is None:
        reader = MockMSIReader(
            MockMSIConfig(n_x=3, n_y=3, n_mz_bins=64, peaks_per_spectrum=(5, 10))
        )
    return StreamingSpatialDataConverter(
        reader,
        tmp_path / "out.zarr",
        dataset_id=_DATASET_ID,
        pixel_size_um=10.0,
        **kwargs,
    )


def _FragmentingGridReader():
    """A source whose mobility sibling declines while its MS/MS one is built.

    Built on the grid stub -- the only fixture here with a mobility
    dimension -- reporting a shared feature axis it cannot list (so
    ``build_mobility_table`` declines) and a two-precursor schedule it
    *can* demultiplex. One slice, two named siblings, one of them written:
    the shape that tells whether unnaming reaches the sibling as well as
    the summed table.
    """
    import numpy as np

    from tests.unit.converters.test_mobility_grid import (
        MASS_AXIS,
        PIXELS,
        GridStubReader,
    )
    from thyra.core.msms import FragmentationSchedule, IsolationWindow

    schedule = FragmentationSchedule(
        ms_level=2,
        windows=(
            IsolationWindow(313.275, 0.5, 0.5, 34.3, scan_begin=150, scan_end=200),
            IsolationWindow(936.578, 0.5, 0.5, 49.6, scan_begin=20, scan_end=60),
        ),
        source="stub",
    )

    class _Reader(GridStubReader):
        @property
        def has_shared_mobility_axis(self) -> bool:
            return True

        def get_shared_mobility_features(self):
            return None

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
                for window in (0, 1):
                    yield (
                        (x, y, 0),
                        window,
                        MASS_AXIS[: 2 + window].copy(),
                        np.arange(2 + window, dtype=np.float64) + p + 1.0,
                    )

    return _Reader()


class TestTheProvenanceBlock:
    """``build_uns_metadata``: one catch for the reader, one for us."""

    def test_metadata_it_does_not_recognise_still_costs_only_that_block(
        self, tmp_path, thyra_logs
    ):
        """The tolerance the broad catch provided is the tolerance kept."""
        converter = _converter(tmp_path)
        # Patched on the assembler, which is what composes the block now;
        # an attribute on the converter would be silently ignored.
        converter.uns._collect_optional_sections = lambda *a, **k: (
            _ for _ in ()
        ).throw(KeyError("format_specific"))

        with thyra_logs("thyra", logging.WARNING) as records:
            uns = converter.build_uns_metadata()

        # The section that raised is gone; the blocks composed after the
        # catch are not, and the failure was said out loud.
        assert "format_specific" not in uns
        assert "msi_metadata" in uns
        assert any(
            "Could not build the full uns provenance block" in r.getMessage()
            for r in records
        )

    def test_an_invariant_break_is_not_downgraded_to_a_warning(self, tmp_path):
        """A RuntimeError is a defect, and a defect's traceback is the point."""
        converter = _converter(tmp_path)
        converter.uns._collect_region_info = lambda *a, **k: (_ for _ in ()).throw(
            RuntimeError("the region table is not built yet")
        )

        with pytest.raises(RuntimeError, match="not built yet"):
            converter.build_uns_metadata()

    def test_the_reader_call_above_it_keeps_its_breadth(self, tmp_path, thyra_logs):
        """The boundary is where an unnameable exception is legitimate."""

        class _RefusingReader(MockMSIReader):
            def get_comprehensive_metadata(self):
                raise RuntimeError("the vendor library said no")

        converter = _converter(
            tmp_path,
            reader=_RefusingReader(
                MockMSIConfig(n_x=3, n_y=3, n_mz_bins=64, peaks_per_spectrum=(5, 10))
            ),
        )

        with thyra_logs("thyra", logging.WARNING) as records:
            uns = converter.build_uns_metadata()

        assert uns == {}
        assert any(
            "Could not read metadata for uns provenance" in r.getMessage()
            for r in records
        )


class TestANamedSibling:
    """A sibling the summed table already named is not additive."""

    def test_the_summed_table_names_the_sibling_it_expects(self, tmp_path):
        """Why the catch below it was harmful, stated as a fact about a store.

        Read off a *successful* conversion: the key is in the summed
        table's ``uns`` whether or not the sibling was ever built.
        """
        from tests.unit.converters.test_mobility_grid import GridStubReader
        from thyra.utils.windows_paths import prepare_zarr_output_path

        converter = _converter(tmp_path, reader=GridStubReader(), mobility_grid=True)
        converter.output_path = prepare_zarr_output_path(
            converter.output_path, _DATASET_ID
        )
        assert converter.convert() is True

        import spatialdata

        from thyra.utils.windows_paths import prepare_zarr_read_path

        sdata = spatialdata.read_zarr(prepare_zarr_read_path(converter.output_path))
        summed = sdata.tables[f"{_DATASET_ID}_z0"]
        named = summed.uns["mobility_axis"]["resolved_table"]

        assert named == f"{_DATASET_ID}_z0_mobility"
        assert named in sdata.tables

    def test_a_refused_mobility_table_fails_the_conversion(self, tmp_path, thyra_logs):
        """It used to be an ERROR line and a store that named a missing table."""
        from tests.unit.converters.test_mobility_grid import GridStubReader
        from thyra.converters.spatialdata import mobility_table
        from thyra.utils.windows_paths import prepare_zarr_output_path

        def _refuse(*args: Any, **kwargs: Dict[str, Any]):
            raise ConversionRefused(
                "2 columns were counted for more entries than were scattered"
            )

        converter = _converter(tmp_path, reader=GridStubReader(), mobility_grid=True)
        converter.output_path = prepare_zarr_output_path(
            converter.output_path, _DATASET_ID
        )

        original = mobility_table.build_mobility_table
        mobility_table.build_mobility_table = _refuse
        try:
            with thyra_logs("thyra.core", logging.ERROR) as records:
                assert converter.convert() is False
        finally:
            mobility_table.build_mobility_table = original

        # Presented as the refusal it is, by convert()'s ConversionRefused
        # branch -- not as "Could not build the mobility-resolved table".
        assert any("were counted for more entries" in r.getMessage() for r in records)


def _store(converter):
    """Convert and read the store back."""
    import spatialdata

    from thyra.utils.windows_paths import (
        prepare_zarr_output_path,
        prepare_zarr_read_path,
    )

    converter.output_path = prepare_zarr_output_path(converter.output_path, _DATASET_ID)
    assert converter.convert() is True
    return spatialdata.read_zarr(prepare_zarr_read_path(converter.output_path))


def _pointers(uns: Dict[str, Any]) -> Dict[str, Any]:
    """Every place in one table's ``uns`` that can name a sibling."""
    analysis = uns.get("msi_metadata", {}).get("ms_analysis", {})
    found: Dict[str, Any] = {}
    for open_block, section in (
        ("mobility_axis", "ion_mobility"),
        ("msms_schedule", "fragmentation"),
    ):
        block = uns.get(open_block)
        if isinstance(block, dict) and block.get("resolved_table"):
            found[f"{open_block}.resolved_table"] = block["resolved_table"]
        versioned = analysis.get(section)
        if isinstance(versioned, dict):
            if versioned.get("resolved_table"):
                found[f"{section}.resolved_table"] = versioned["resolved_table"]
            if versioned.get("grid"):
                found[f"{section}.grid"] = "<present>"
    return found


class TestADeclinedSibling:
    """A sibling that says no leaves nothing naming it (issue #343).

    The other half of :class:`TestANamedSibling`. There a builder *failed*
    and the failure was swallowed; here nothing fails -- the builder
    returns ``None``, which is a decision it is entitled to make -- and the
    name it was given before it ran has to come back out.
    """

    def test_a_shared_axis_it_cannot_list_names_nothing(self, tmp_path):
        """The reachable case, and the one the issue was measured on.

        ``_plan_mobility_table`` names the table on the shared-axis route
        with no further check, and ``_build_from_shared_axis`` declines
        whenever the listing comes back empty.
        """
        from tests.unit.converters.test_mobility_grid import GridStubReader

        class SharedButEmpty(GridStubReader):
            @property
            def has_shared_mobility_axis(self) -> bool:
                return True

            def get_shared_mobility_features(self):
                return None

        sdata = _store(_converter(tmp_path, reader=SharedButEmpty()))

        assert set(sdata.tables) == {f"{_DATASET_ID}_z0"}
        assert _pointers(sdata.tables[f"{_DATASET_ID}_z0"].uns) == {}

    def test_the_grid_description_goes_with_the_name(self, tmp_path):
        """``grid`` describes the binning of the table that was not written.

        Left behind it would say a table was binned onto a grid nobody
        wrote -- the same lie as the name, one field over.
        """
        from tests.unit.converters.test_mobility_grid import GridStubReader
        from thyra.converters.spatialdata import mobility_table

        converter = _converter(tmp_path, reader=GridStubReader(), mobility_grid=True)
        original = mobility_table.build_mobility_table
        mobility_table.build_mobility_table = lambda *a, **k: None
        try:
            sdata = _store(converter)
        finally:
            mobility_table.build_mobility_table = original

        summed = sdata.tables[f"{_DATASET_ID}_z0"]
        assert set(sdata.tables) == {f"{_DATASET_ID}_z0"}
        assert _pointers(summed.uns) == {}
        # The acquisition still had the dimension; only the table's
        # description is withdrawn.
        assert summed.uns["msi_metadata"]["ms_analysis"]["ion_mobility"]["present"]
        assert "mobility_axis" in summed.uns

    def test_a_surviving_sibling_does_not_name_the_declined_one(self, tmp_path):
        """The siblings' own ``uns`` is taken before either builder runs.

        So the MS/MS table inherits the summed table's pointer to the
        mobility table, and cleaning only the summed table would move the
        dangling pointer rather than remove it.
        """
        sdata = _store(_converter(tmp_path, reader=_FragmentingGridReader()))

        assert set(sdata.tables) == {
            f"{_DATASET_ID}_z0",
            f"{_DATASET_ID}_z0_msms",
        }
        msms_key = f"{_DATASET_ID}_z0_msms"
        for key in (f"{_DATASET_ID}_z0", msms_key):
            pointers = _pointers(sdata.tables[key].uns)
            assert "ion_mobility.resolved_table" not in pointers
            assert "mobility_axis.resolved_table" not in pointers
            # The one that was written is still named, on both tables.
            assert pointers["fragmentation.resolved_table"] == msms_key

    def test_a_written_sibling_keeps_its_name(self, tmp_path):
        """The guard: unnaming must not reach a table that was written."""
        from tests.unit.converters.test_mobility_grid import GridStubReader

        sdata = _store(
            _converter(tmp_path, reader=GridStubReader(), mobility_grid=True)
        )

        summed = sdata.tables[f"{_DATASET_ID}_z0"]
        pointers = _pointers(summed.uns)
        assert pointers["mobility_axis.resolved_table"] == f"{_DATASET_ID}_z0_mobility"
        assert pointers["ion_mobility.resolved_table"] == f"{_DATASET_ID}_z0_mobility"
        assert pointers["ion_mobility.grid"] == "<present>"
        assert f"{_DATASET_ID}_z0_mobility" in sdata.tables
