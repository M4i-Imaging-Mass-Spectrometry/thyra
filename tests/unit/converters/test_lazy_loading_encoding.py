"""anndata encoding metadata on the hand-written PCS Zarr layout.

The streaming PCS path (``use_csc=True``) does not go through anndata's
writer: it hand-writes the AnnData layout into Zarr itself.  Every
``encoding-type`` / ``encoding-version`` attribute in that table is therefore
written by Thyra, and nothing upstream will notice if one is forgotten -- the
store simply stops being readable by
``anndata.experimental.read_lazy()``, which is how downstream consumers open
Thyra output without loading it into memory.

Nothing tested those attributes until now: ``git grep encoding-type tests/``
came back empty.  The coverage here started life on
``feature/lazy-loading-support``, whose *implementation* ``main`` superseded
and improved on -- ``main`` distinguishes ``"string"`` for the 0-d scalars in
``uns`` from ``"string-array"`` for real arrays, which that branch's binary
``_set_encoding_attrs(is_string=...)`` helper could not express.  The branch
is gone; the test it carried lives here, rewritten against the layout ``main``
actually writes.

The tree walk is the load-bearing one: it fails when a new array is added to
the hand-written layout without its encoding attrs, which is the way this
breaks in practice.
"""

from typing import Dict, Iterator, Tuple

import pytest
import zarr

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.converters.spatialdata.streaming_converter import (
    SPATIALDATA_AVAILABLE,
    StreamingSpatialDataConverter,
)

pytestmark = pytest.mark.skipif(
    not SPATIALDATA_AVAILABLE,
    reason="SpatialData dependencies not available",
)

_DATASET_ID = "mock"
_TABLE_NAME = f"{_DATASET_ID}_z0"
_N_X = 4
_N_Y = 4
_N_MZ_BINS = 40

# The attributes anndata's spec requires, for the nodes whose *value* matters
# rather than merely being present. Scalars written into ``uns`` are
# ``"string"``; anything with a length is ``"string-array"``. Getting that
# split wrong is what stops read_lazy reconstructing the table.
_EXPECTED_ENCODINGS: Dict[str, Tuple[str, str]] = {
    "": ("anndata", "0.1.0"),
    "X": ("csc_matrix", "0.1.0"),
    "X/data": ("array", "0.2.0"),
    "X/indices": ("array", "0.2.0"),
    "X/indptr": ("array", "0.2.0"),
    "obs": ("dataframe", "0.2.0"),
    "obs/instance_id": ("string-array", "0.2.0"),
    "obs/instance_key": ("string-array", "0.2.0"),
    "obs/region": ("categorical", "0.2.0"),
    "obs/region/categories": ("string-array", "0.2.0"),
    "obs/region/codes": ("array", "0.2.0"),
    "obs/x": ("array", "0.2.0"),
    "obs/y": ("array", "0.2.0"),
    "obs/spatial_x": ("array", "0.2.0"),
    "obs/spatial_y": ("array", "0.2.0"),
    "var": ("dataframe", "0.2.0"),
    "var/_index": ("string-array", "0.2.0"),
    "var/mz": ("array", "0.2.0"),
    "uns": ("dict", "0.1.0"),
    "uns/average_spectrum": ("array", "0.2.0"),
    "uns/essential_metadata": ("dict", "0.1.0"),
    "uns/essential_metadata/dimensions": ("array", "0.2.0"),
    "uns/essential_metadata/mass_range": ("array", "0.2.0"),
    "uns/essential_metadata/source_path": ("string", "0.2.0"),
    "uns/essential_metadata/spectrum_type": ("string", "0.2.0"),
    "uns/essential_metadata/thyra_version": ("string", "0.2.0"),
    "uns/spatialdata_attrs": ("dict", "0.1.0"),
    "uns/spatialdata_attrs/region": ("string", "0.2.0"),
    "uns/spatialdata_attrs/region_key": ("string", "0.2.0"),
    "uns/spatialdata_attrs/instance_key": ("string", "0.2.0"),
}


@pytest.fixture(scope="module")
def pcs_store(tmp_path_factory):
    """Convert once via the PCS path and hand back the output path.

    Module-scoped: every test here reads the same immutable store, and the
    conversion is the slow part.
    """
    output_path = tmp_path_factory.mktemp("pcs_encoding") / "pcs.zarr"
    converter = StreamingSpatialDataConverter(
        reader=MockMSIReader(
            MockMSIConfig(
                n_x=_N_X, n_y=_N_Y, n_mz_bins=_N_MZ_BINS, peaks_per_spectrum=(5, 10)
            )
        ),
        output_path=output_path,
        dataset_id=_DATASET_ID,
        pixel_size_um=10.0,
        use_csc=True,  # force the hand-written PCS layout
    )
    assert converter.convert() is True, "PCS conversion should succeed"
    return output_path


def _table_group(output_path) -> zarr.Group:
    return zarr.open_group(str(output_path), mode="r")["tables"][_TABLE_NAME]


def _walk(node, path: str = "") -> Iterator[Tuple[str, object]]:
    """Yield ``(path, node)`` for the table group and everything beneath it."""
    yield path, node
    if isinstance(node, zarr.Group):
        for name in sorted(node):
            child_path = f"{path}/{name}" if path else name
            yield from _walk(node[name], child_path)


def test_read_lazy_opens_the_pcs_store(pcs_store):
    """The end-to-end capability the encoding attrs exist to serve.

    Deliberately coarse: read_lazy tolerates a fair amount, and dropping a
    single scalar's attrs does not stop it opening the store (measured). It
    covers the shape of the layout; the per-node tests below cover the parts
    it would let through.
    """
    anndata = pytest.importorskip("anndata")

    lazy = anndata.experimental.read_lazy(str(pcs_store / "tables" / _TABLE_NAME))

    assert lazy.shape == (_N_X * _N_Y, _N_MZ_BINS)
    assert list(lazy.var.index[:3]) == ["mz_0", "mz_1", "mz_2"]
    # region must come back as a categorical, which only works if the
    # categories/codes pair carries the right encodings.
    assert str(lazy.obs["region"].dtype) == "category"
    assert {"spatialdata_attrs", "essential_metadata"} <= set(lazy.uns)


def test_every_table_node_declares_its_encoding(pcs_store):
    """No array or group in the hand-written layout may omit the attrs.

    This is the regression guard for adding something to the PCS layout and
    forgetting its encoding metadata.
    """
    missing = [
        path or "<table root>"
        for path, node in _walk(_table_group(pcs_store))
        if "encoding-type" not in node.attrs or "encoding-version" not in node.attrs
    ]

    assert not missing, f"nodes without encoding attrs: {missing}"


@pytest.mark.parametrize("node_path", sorted(_EXPECTED_ENCODINGS))
def test_node_declares_the_expected_encoding(pcs_store, node_path):
    """Each known node carries the exact encoding anndata's spec wants."""
    table = _table_group(pcs_store)
    node = table[node_path] if node_path else table

    expected_type, expected_version = _EXPECTED_ENCODINGS[node_path]
    assert node.attrs["encoding-type"] == expected_type
    assert node.attrs["encoding-version"] == expected_version


def test_uns_scalars_are_string_not_string_array(pcs_store):
    """The distinction ``main`` draws that the superseded branch could not.

    A 0-d scalar in ``uns`` is ``"string"``; anything with a length is
    ``"string-array"``. Writing ``"string-array"`` for a scalar is exactly the
    bug a binary is-it-a-string helper produces, and it makes read_lazy
    reconstruct ``uns`` wrongly.
    """
    table = _table_group(pcs_store)

    for path in (
        "uns/spatialdata_attrs/region",
        "uns/spatialdata_attrs/region_key",
        "uns/spatialdata_attrs/instance_key",
        "uns/essential_metadata/source_path",
        "uns/essential_metadata/spectrum_type",
        "uns/essential_metadata/thyra_version",
    ):
        node = table[path]
        assert node.shape == (), f"{path} should be 0-d, got shape {node.shape}"
        assert node.attrs["encoding-type"] == "string", path

    for path in ("obs/instance_id", "obs/region/categories", "var/_index"):
        node = table[path]
        assert node.shape != (), f"{path} should have a length"
        assert node.attrs["encoding-type"] == "string-array", path
