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


class TestTheProvenanceBlock:
    """``build_uns_metadata``: one catch for the reader, one for us."""

    def test_metadata_it_does_not_recognise_still_costs_only_that_block(
        self, tmp_path, thyra_logs
    ):
        """The tolerance the broad catch provided is the tolerance kept."""
        converter = _converter(tmp_path)
        converter._collect_optional_sections = lambda *a, **k: (_ for _ in ()).throw(
            KeyError("format_specific")
        )

        with thyra_logs("thyra.converters", logging.WARNING) as records:
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
        converter._collect_region_info = lambda *a, **k: (_ for _ in ()).throw(
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

        with thyra_logs("thyra.converters", logging.WARNING) as records:
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
