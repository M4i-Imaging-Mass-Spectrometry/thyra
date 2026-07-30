"""Write-path coverage under pandas 3's string dtypes.

Under pandas 3.0 -- or pandas 2.x with ``future.infer_string=True``, which is
what 3.0 makes the default -- the table's obs index (``instance_id``), the
``instance_key`` column and the var index (``mz_*``) are inferred as pandas'
``str`` dtype, backed by ``ArrowStringArrayNumpySemantics``.  anndata's IO
registry has a writer for ``pandas.core.arrays.string_.StringArray`` but not
for the ``*NumpySemantics`` subclasses pandas 3 actually produces, and it
matches exact types rather than subclasses, so ``SpatialData.write()`` raises::

    No method registered for writing
    <class 'pandas.core.arrays.string_arrow.ArrowStringArrayNumpySemantics'>
    into <class 'zarr.core.group.Group'>

``BaseSpatialDataConverter._coerce_table_strings_to_object`` is the fix.  Both
write paths that reach anndata's writer fail without it: the in-memory
converter (the default for anything under the 10 GB streaming threshold, i.e.
most conversions) and the streaming COO path.  The streaming PCS path
hand-writes the AnnData layout straight to Zarr and never reaches anndata's
writer, so it passes either way -- it is covered here to keep it that way.

The rest of the suite runs with pandas' default inference, so nothing else
would notice if the coercion were dropped; every test in this module turns
``future.infer_string`` on for its duration and restores it afterwards.

``test_anndata_can_write_arrow_backed_strings`` is the removal trigger: see
``_coerce_table_strings_to_object`` for the full anndata #2221 story.
"""

from typing import Callable, Dict

import numpy as np
import pandas as pd
import pytest

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.converters.spatialdata.base_spatialdata_converter import (
    SPATIALDATA_AVAILABLE,
    BaseSpatialDataConverter,
)
from thyra.converters.spatialdata.spatialdata_2d_converter import SpatialData2DConverter
from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)

pytestmark = pytest.mark.skipif(
    not SPATIALDATA_AVAILABLE,
    reason="SpatialData dependencies not available",
)

_N_X = 6
_N_Y = 6
_N_MZ_BINS = 500
_N_PIXELS = _N_X * _N_Y


@pytest.fixture(autouse=True)
def infer_string():
    """Run every test in this module with pandas 3 string inference on.

    ``pd.option_context`` restores the previous value on exit, so the option
    cannot leak into the rest of the suite.  The option still exists on pandas
    3 (where it defaults to ``True``), so this is a no-op there rather than a
    behaviour change.
    """
    with pd.option_context("future.infer_string", True):
        yield


def _small_config() -> MockMSIConfig:
    return MockMSIConfig(
        n_x=_N_X, n_y=_N_Y, n_mz_bins=_N_MZ_BINS, peaks_per_spectrum=(20, 40)
    )


def _in_memory(output_path):
    """The default (``streaming=False``) converter."""
    return SpatialData2DConverter(
        reader=MockMSIReader(_small_config()),
        output_path=output_path,
        dataset_id="mock",
        pixel_size_um=10.0,
    )


def _streaming_pcs(output_path):
    """``streaming=True, use_csc=True`` -- the hand-written Zarr layout."""
    return StreamingSpatialDataConverter(
        reader=MockMSIReader(_small_config()),
        output_path=output_path,
        dataset_id="mock",
        pixel_size_um=10.0,
        use_csc=True,
    )


def _streaming_coo(output_path):
    """``streaming=True, use_csc=False`` -- writes via ``SpatialData.write()``."""
    return StreamingSpatialDataConverter(
        reader=MockMSIReader(_small_config()),
        output_path=output_path,
        dataset_id="mock",
        pixel_size_um=10.0,
        use_csc=False,
    )


WRITE_PATHS: Dict[str, Callable] = {
    "in_memory": _in_memory,
    "streaming_pcs": _streaming_pcs,
    "streaming_coo": _streaming_coo,
}


def _string_obs_table() -> pd.DataFrame:
    """An obs table shaped like the converters build it, with pandas 3 dtypes.

    Mirrors what ``_create_coordinates_dataframe`` plus the ``region`` /
    ``instance_key`` fixups produce: a string index, a string column, a
    string-backed categorical with more than one category, and numeric
    columns that must be left alone.
    """
    return pd.DataFrame(
        {
            "x": np.arange(3, dtype=np.int32),
            "spatial_x": np.arange(3, dtype=np.float64) * 10.0,
            "region": pd.Categorical(np.array(["mock_z0_pixels", "other", "other"])),
            "instance_key": np.array(["0", "1", "2"]),
        },
        index=pd.Index(np.array(["0", "1", "2"]), name="instance_id"),
    )


@pytest.mark.parametrize("path_name", list(WRITE_PATHS))
def test_write_path_converts_under_infer_string(tmp_path, path_name):
    """All three write paths must convert with pandas 3 string dtypes.

    ``in_memory`` and ``streaming_coo`` both fail here without the coercion.
    """
    # Own output path per case: reusing one filename makes the second
    # conversion fail with "Destination already exists", which reads like a
    # real regression.
    output_path = tmp_path / f"{path_name}.zarr"

    converter = WRITE_PATHS[path_name](output_path)
    assert converter.convert() is True, f"{path_name} conversion failed"

    import spatialdata

    sdata = spatialdata.read_zarr(str(output_path))
    assert len(sdata.tables) == 1
    table = list(sdata.tables.values())[0]
    assert table.shape == (_N_PIXELS, _N_MZ_BINS)
    assert "mock_z0_tic" in sdata.images
    assert "mock_z0_pixels" in sdata.shapes


def test_written_store_reads_back_with_pandas_string_index(tmp_path):
    """The coercion is write-time only -- the on-disk format is unchanged.

    Readers still get pandas' native string dtype for the obs and var indices,
    exactly as they would from a pandas 2 write, so the workaround does not
    leak ``object`` dtypes into the stored layout.
    """
    output_path = tmp_path / "readback.zarr"
    assert _in_memory(output_path).convert() is True

    import spatialdata

    table = list(spatialdata.read_zarr(str(output_path)).tables.values())[0]

    assert isinstance(table.obs.index.dtype, pd.StringDtype)
    assert isinstance(table.var.index.dtype, pd.StringDtype)
    assert table.obs.index.tolist() == [str(i) for i in range(_N_PIXELS)]
    assert table.var.index.tolist() == [f"mz_{i}" for i in range(_N_MZ_BINS)]
    assert isinstance(table.obs["region"].dtype, pd.CategoricalDtype)
    assert table.obs["region"].tolist() == ["mock_z0_pixels"] * _N_PIXELS


def test_coercion_changes_only_dtypes():
    """Same order, same content, only the string dtypes move to ``object``."""
    df = _string_obs_table()

    # Precondition: without these the rest of the test is vacuous.
    assert isinstance(df.index.dtype, pd.StringDtype)
    assert isinstance(df["instance_key"].dtype, pd.StringDtype)

    index_before = df.index.tolist()
    instance_key_before = df["instance_key"].tolist()
    x_before = df["x"].tolist()
    spatial_x_before = df["spatial_x"].tolist()

    BaseSpatialDataConverter._coerce_table_strings_to_object(df)

    assert df.index.dtype == object
    assert df["instance_key"].dtype == object
    assert df.index.tolist() == index_before
    assert df.index.name == "instance_id"
    assert df["instance_key"].tolist() == instance_key_before

    # Non-string columns are untouched, dtype included.
    assert df["x"].dtype == np.int32
    assert df["spatial_x"].dtype == np.float64
    assert df["x"].tolist() == x_before
    assert df["spatial_x"].tolist() == spatial_x_before

    assert list(df.columns) == ["x", "spatial_x", "region", "instance_key"]


def test_coercion_preserves_categorical_codes_and_categories():
    """A string-backed categorical keeps its codes and its category values."""
    df = _string_obs_table()

    # Precondition: the categorical must be string-backed to be worth testing.
    assert isinstance(df["region"].dtype, pd.CategoricalDtype)
    assert isinstance(df["region"].cat.categories.dtype, pd.StringDtype)

    codes_before = df["region"].cat.codes.tolist()
    categories_before = df["region"].cat.categories.tolist()
    values_before = df["region"].tolist()
    assert len(set(codes_before)) > 1, "need >1 category for codes to mean anything"

    BaseSpatialDataConverter._coerce_table_strings_to_object(df)

    assert isinstance(df["region"].dtype, pd.CategoricalDtype)
    assert df["region"].cat.categories.dtype == object
    assert df["region"].cat.codes.tolist() == codes_before
    assert df["region"].cat.categories.tolist() == categories_before
    assert df["region"].tolist() == values_before


def test_coercion_is_a_noop_on_object_dtypes():
    """On pandas < 3 the same columns are already ``object`` -- nothing moves.

    Built by coercing once, so this also covers a table that has already been
    through ``_save_output``.
    """
    df = _string_obs_table()
    BaseSpatialDataConverter._coerce_table_strings_to_object(df)

    # Precondition: this is the pandas < 3 shape.
    assert df.index.dtype == object
    assert df["instance_key"].dtype == object
    assert df["region"].cat.categories.dtype == object

    before = df.copy(deep=True)
    BaseSpatialDataConverter._coerce_table_strings_to_object(df)

    pd.testing.assert_frame_equal(df, before)


@pytest.mark.xfail(
    strict=False,
    reason=(
        "The pinned anndata cannot serialize pandas' str dtype -- see anndata "
        "#2221 ('Pandas 3.0 compatibility', milestone 0.14.0). An XPASS here "
        "means the installed anndata can write Arrow-backed strings and "
        "_coerce_table_strings_to_object can be deleted."
    ),
)
def test_anndata_can_write_arrow_backed_strings(tmp_path):
    """Removal trigger: XPASSes when the workaround becomes unnecessary.

    Deliberately non-strict, so an anndata that gains pandas 3 support shows up
    as an XPASS in the run summary rather than as a red build.  The point is
    that the coercion announces its own obsolescence instead of silently
    pessimising every write forever.
    """
    import anndata as ad
    import zarr

    df = _string_obs_table()
    assert isinstance(df.index.dtype, pd.StringDtype)  # precondition

    group = zarr.open_group(str(tmp_path / "probe.zarr"), mode="w")
    ad.io.write_elem(group, "obs", df)

    written_back = ad.io.read_elem(group["obs"])
    assert written_back.index.tolist() == df.index.tolist()
    assert written_back["instance_key"].tolist() == df["instance_key"].tolist()
