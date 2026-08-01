"""The provenance block in ``uns`` must survive, and agree, on every path.

A converted store is supposed to be able to say where it came from and how
it was interpreted: ``uns["essential_metadata"]`` carries ``source_path``,
``dimensions``, ``mass_range``, ``spectrum_type`` and the ``thyra_version``
that wrote it, and the sections beside it carry the vendor metadata.  Ousia
reads those back out of the store; nothing else records them.

There are three write paths and they had drifted:

* the in-memory converter and the streaming-COO path hand the table to
  anndata's writer, which serialises ``adata.uns`` as-is, and
* the streaming-PCS path hand-writes the Zarr layout, and used to compose
  its own, much smaller block -- ``spectrum_type`` hardcoded to
  ``"processed"``, ``mass_range`` taken from the *resampled target axis*
  rather than the source, ``source_path`` from a different attribute, and
  ``format_specific`` / ``instrument_info`` / ``raw_metadata`` / ``regions``
  dropped entirely.

Which path a dataset takes is decided by size, in
``StreamingSpatialDataConverter._should_use_pcs`` -- so two acquisitions
off the same instrument landed on opposite sides of the split and came out
described differently.  On ``test_data/``: ``pea.imzML`` estimates 29.9 GB
(COO, complete provenance) and ``bellini.imzML`` 43.3 GB (PCS, wrong
``spectrum_type``, everything else missing).  Ousia's import wizard forces
``use_csc=True``, so *every* wizard conversion took the PCS path.

Nothing caught it because ``MockMSIReader`` reported
``spectrum_type="processed"`` too -- the fixture agreed with the hardcoded
literal, so a test could compare them and pass.  The fixture now reports a
value from the vocabulary the real extractors produce; see the comment on
it.

These tests convert the same mock dataset down all three paths and compare
what each one stored.
"""

from typing import Any, Callable, Dict

import numpy as np
import pytest

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.converters.spatialdata.base_spatialdata_converter import (
    SPATIALDATA_AVAILABLE,
)
from thyra.converters.spatialdata.spatialdata_2d_converter import SpatialData2DConverter
from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)

pytestmark = pytest.mark.skipif(
    not SPATIALDATA_AVAILABLE,
    reason="SpatialData dependencies not available",
)

_DATASET_ID = "mock"
_TABLE_NAME = f"{_DATASET_ID}_z0"
_N_X = 6
_N_Y = 6
_N_MZ_BINS = 500

# What the mock reader reports, and therefore what must come back out of
# every store. Kept as literals rather than read off the fixture: the point
# is that the stored value tracks the source, and a test that derives the
# expectation from the same place the writer does cannot show that.
_EXPECTED_SPECTRUM_TYPE = "centroid spectrum"
_EXPECTED_SOURCE_PATH = "mock_msi_data"
_EXPECTED_ESSENTIAL_KEYS = {
    "source_path",
    "dimensions",
    "mass_range",
    "spectrum_type",
    "thyra_version",
}


def _config() -> MockMSIConfig:
    return MockMSIConfig(
        n_x=_N_X, n_y=_N_Y, n_mz_bins=_N_MZ_BINS, peaks_per_spectrum=(20, 40)
    )


def _in_memory(output_path):
    """The default converter -- anything under the streaming threshold."""
    return SpatialData2DConverter(
        reader=MockMSIReader(_config()),
        output_path=output_path,
        dataset_id=_DATASET_ID,
        pixel_size_um=10.0,
    )


def _streaming_coo(output_path):
    """``streaming=True, use_csc=False`` -- writes via ``SpatialData.write()``."""
    return StreamingSpatialDataConverter(
        reader=MockMSIReader(_config()),
        output_path=output_path,
        dataset_id=_DATASET_ID,
        pixel_size_um=10.0,
        use_csc=False,
    )


def _streaming_pcs(output_path):
    """``streaming=True, use_csc=True`` -- the hand-written Zarr layout."""
    return StreamingSpatialDataConverter(
        reader=MockMSIReader(_config()),
        output_path=output_path,
        dataset_id=_DATASET_ID,
        pixel_size_um=10.0,
        use_csc=True,
    )


WRITE_PATHS: Dict[str, Callable] = {
    "in_memory": _in_memory,
    "streaming_coo": _streaming_coo,
    "streaming_pcs": _streaming_pcs,
}


def _convert(tmp_path_factory, path_name: str):
    """Convert once through one write path; hand back the table group path."""
    output_path = tmp_path_factory.mktemp(f"prov_{path_name}") / "out.zarr"
    converter = WRITE_PATHS[path_name](output_path)
    assert converter.convert() is True, f"{path_name} conversion failed"
    return output_path / "tables" / _TABLE_NAME


@pytest.fixture(scope="module")
def stores(tmp_path_factory) -> Dict[str, Any]:
    """One conversion per write path. Module-scoped: converting is the slow part."""
    return {name: _convert(tmp_path_factory, name) for name in WRITE_PATHS}


def _read_uns(table_path) -> Dict[str, Any]:
    """Read ``uns`` back eagerly, as a plain dict."""
    import anndata as ad
    import zarr

    return ad.io.read_elem(zarr.open_group(str(table_path), mode="r")["uns"])


def _plain(value: Any) -> Any:
    """Normalise for comparison across paths.

    A list written by one path and an ndarray by another are the same
    stored value; ``==`` on ndarrays is not a bool. Recurse so nested
    sections compare too.
    """
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


@pytest.mark.parametrize("path_name", list(WRITE_PATHS))
def test_essential_metadata_is_not_empty(stores, path_name):
    """The reported symptom: the block present but with nothing in it.

    Deliberately the weakest assertion here, and deliberately first --
    a store whose provenance is empty cannot say what it came from at
    all, which is worse than one that says something inaccurate.
    """
    uns = _read_uns(stores[path_name])

    assert "essential_metadata" in uns, f"{path_name} wrote no essential_metadata"
    essential = uns["essential_metadata"]
    assert essential, f"{path_name} wrote an EMPTY essential_metadata block"
    assert set(essential) == _EXPECTED_ESSENTIAL_KEYS


@pytest.mark.parametrize("path_name", list(WRITE_PATHS))
def test_stored_spectrum_type_is_the_readers(stores, path_name):
    """``spectrum_type`` must come from the data, not from a literal.

    The PCS path wrote ``"processed"`` for everything. That is not a
    value any extractor produces, so a consumer could neither trust it
    nor tell it apart from a real reading -- which is the whole failure:
    losing provenance *silently*.
    """
    essential = _read_uns(stores[path_name])["essential_metadata"]

    assert essential["spectrum_type"] == _EXPECTED_SPECTRUM_TYPE, (
        f"{path_name} stored spectrum_type={essential['spectrum_type']!r}; "
        f"the reader reports {_EXPECTED_SPECTRUM_TYPE!r}"
    )
    assert essential["spectrum_type"] != "processed"


@pytest.mark.parametrize("path_name", list(WRITE_PATHS))
def test_essential_metadata_values_track_the_source(stores, path_name):
    """The rest of the block must describe the source, not the output."""
    essential = _read_uns(stores[path_name])["essential_metadata"]
    cfg = _config()

    assert essential["source_path"] == _EXPECTED_SOURCE_PATH
    assert _plain(essential["dimensions"]) == [cfg.n_x, cfg.n_y, cfg.n_z]
    # mass_range is the SOURCE range. The PCS path used to take it from
    # the resampled target axis, which is a different quantity whenever
    # resampling clips or extends the range.
    assert _plain(essential["mass_range"]) == pytest.approx([cfg.mz_min, cfg.mz_max])
    assert essential["thyra_version"]


@pytest.mark.parametrize("path_name", list(WRITE_PATHS))
def test_vendor_sections_are_stored(stores, path_name):
    """The sections beside essential_metadata, which PCS dropped entirely.

    ``acquisition_params`` is deliberately not asserted: the mock leaves
    it empty and empty sections are omitted, so that a consumer can tell
    "this format has none" from "this one has none recorded".
    """
    uns = _read_uns(stores[path_name])

    assert uns.get("format_specific") == {"format": "mock"}
    assert uns.get("instrument_info") == {"instrument": "mock"}
    assert uns.get("raw_metadata") == {"source": "mock"}
    assert "acquisition_params" not in uns
    assert uns.get("regions"), f"{path_name} stored no region summary"


def test_every_path_stores_the_same_provenance(stores):
    """The invariant the individual assertions above are instances of.

    Compares whole blocks rather than named keys, so a section added to
    one path and forgotten on another fails here without anyone having
    to remember to extend this file.
    """
    provenance_keys = (
        "essential_metadata",
        "format_specific",
        "acquisition_params",
        "instrument_info",
        "raw_metadata",
        "regions",
    )
    blocks = {
        name: {k: _plain(v) for k, v in _read_uns(path).items() if k in provenance_keys}
        for name, path in stores.items()
    }

    reference_name = "in_memory"
    reference = blocks[reference_name]
    for name, block in blocks.items():
        if name == reference_name:
            continue
        assert (
            block == reference
        ), f"{name} stored different provenance from {reference_name}"


@pytest.mark.parametrize("path_name", list(WRITE_PATHS))
def test_read_lazy_sees_the_block(stores, path_name):
    """The access mode a consumer actually uses.

    ``read_lazy`` is how Ousia opens a converted store without pulling
    the matrix into memory, and it is where the loss was noticed. Array
    entries come back as dask, so this checks the keys are reachable --
    the values are covered eagerly above.
    """
    anndata = pytest.importorskip("anndata")

    lazy = anndata.experimental.read_lazy(str(stores[path_name]))
    essential = lazy.uns["essential_metadata"]

    assert set(essential) == _EXPECTED_ESSENTIAL_KEYS
    assert essential["spectrum_type"] == _EXPECTED_SPECTRUM_TYPE
