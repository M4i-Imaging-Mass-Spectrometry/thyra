"""Foundation tests for the Thyra chunk/shard policy.

Two halves, for the two things Thyra writes:

**Rasters.** Locks the current chunk tuples (so the future zarr v3 sharding
wiring, gated on spatialdata PR #1106, cannot silently change output) and
documents the cross-repo INNER_TILE_EDGE contract with Ousia's serving tile
size.

**The MSI table.** Locks the shard size the table's ``X`` arrays land in. This
half exists because of a real, shipped regression: anndata >= 0.13 writes every
zarr v3 array with ``shards="auto"``, and zarr < 3.1.4 sized the auto chunk with
``max_bytes=1024`` where it meant 1 MiB. Nothing in Thyra changed -- a
dependency bump alone was enough to turn one table into ~400,000 files of
roughly 1 KB each, and the zarr write into 98% of conversion wall-clock
(pea.imzML: 434 s of 445 s).

That failure is invisible to every other test in this suite. The store still
round-trips, ``read_lazy`` still opens it, and the values are all correct --
only the geometry is catastrophic. Hence a test that asserts geometry, on a
store big enough for the geometry to be non-degenerate: below ~1 KB of ``data``
the broken and fixed layouts agree, so a 4-pixel fixture would prove nothing.
"""

from __future__ import annotations

import contextlib
from pathlib import Path

import numpy as np
import pytest
import zarr

from thyra.converters.spatialdata import _chunking
from thyra.converters.spatialdata._chunking import (
    IMAGE_CHUNK_EDGE,
    INNER_TILE_EDGE,
    MIN_TABLE_SHARD_BYTES,
    TABLE_SHARD_TARGET_BYTES,
    image_chunks,
    table_write_config,
)

_ZARR_CONFIG_KEY = "array.target_shard_size_bytes"


class TestRasterChunkPolicy:
    def test_image_chunks_2d_reproduces_todays_literal(self):
        # No-op-refactor lock: the helper must return exactly the (1, 4096, 4096)
        # literal it replaced at base_spatialdata_converter.py's optical parse.
        assert image_chunks(2) == (1, 4096, 4096)

    def test_image_chunks_3d_prepends_a_z_chunk(self):
        assert image_chunks(3) == (1, 1, 4096, 4096)

    def test_image_chunks_rejects_other_ndim(self):
        with pytest.raises(ValueError):
            image_chunks(4)

    def test_inner_tile_edge_matches_ousia_serving_contract(self):
        # Duplicated-by-contract with Ousia DEFAULT_TILE_SIZE = 512 (Thyra cannot
        # import Ousia). If Ousia's serving tile changes, update this in lockstep.
        assert INNER_TILE_EDGE == 512
        # Outer chunk must tile exactly into inner tiles (no ragged remainder).
        assert IMAGE_CHUNK_EDGE % INNER_TILE_EDGE == 0


class TestTableWriteConfig:
    """The zarr-config seam itself."""

    def test_publishes_thyras_budget_while_open(self):
        with table_write_config():
            assert zarr.config.get(_ZARR_CONFIG_KEY) == TABLE_SHARD_TARGET_BYTES

    def test_restores_the_previous_value_on_exit(self):
        before = zarr.config.get(_ZARR_CONFIG_KEY, None)
        with table_write_config():
            pass
        assert zarr.config.get(_ZARR_CONFIG_KEY, None) == before

    def test_budget_is_an_int_so_anndata_yields_to_it(self):
        # anndata only skips its own 1 GB default when the configured value is
        # an int (anndata._io.specs.methods.zarr_v3_sharding). A float here
        # would silently hand the budget back to anndata.
        assert isinstance(TABLE_SHARD_TARGET_BYTES, int)
        assert TABLE_SHARD_TARGET_BYTES > MIN_TABLE_SHARD_BYTES


@pytest.fixture(scope="module")
def converted_store(tmp_path_factory) -> tuple[Path, int]:
    """A real converted store, plus how many times the policy context opened.

    Sized so ``X/data`` is a few hundred KB: comfortably past the ~1 KB cliff
    where the broken and fixed layouts diverge, and still a couple of seconds
    to convert.
    """
    from pyimzml.ImzMLWriter import ImzMLWriter

    from thyra.convert import convert_msi

    n_side, n_mz = 16, 200
    tmp_path = tmp_path_factory.mktemp("table_chunks")
    imzml_path = tmp_path / "grid.imzML"

    mzs = np.linspace(100.0, 1000.0, n_mz)
    with ImzMLWriter(str(imzml_path), mode="processed") as writer:
        for x in range(1, n_side + 1):
            for y in range(1, n_side + 1):
                # Every channel non-zero, so nnz is n_side**2 * n_mz and the
                # array size the assertions rely on is not a function of how
                # aggressively the converter drops zeros.
                intensities = np.full(n_mz, 1.0) + x + y
                writer.addSpectrum(mzs, intensities, (x, y, 1))

    entered = 0
    real = _chunking.table_write_config

    @contextlib.contextmanager
    def counting_config():
        nonlocal entered
        entered += 1
        with real():
            yield

    output_path = tmp_path / "out.zarr"
    with pytest.MonkeyPatch.context() as mp:
        # Patch the name the converter module bound at import time.
        mp.setattr(
            "thyra.converters.spatialdata.base_spatialdata_converter."
            "table_write_config",
            counting_config,
        )
        assert convert_msi(
            str(imzml_path),
            str(output_path),
            format_type="spatialdata",
            dataset_id="ds",
            pixel_size_um=10.0,
        ), "conversion fixture failed"

    return output_path, entered


class TestConvertedTableGeometry:
    def test_save_output_holds_the_policy_open(self, converted_store):
        # Catches the wiring being dropped from _save_output. The geometry
        # assertion below cannot: with the zarr floor in place, anndata's own
        # 1 GB fallback also yields large shards on a fixture this size.
        _, entered = converted_store
        assert entered >= 1, "_save_output wrote the table outside table_write_config()"

    @pytest.mark.parametrize("member", ["data", "indices", "indptr"])
    def test_x_arrays_are_not_kilobyte_sharded(self, converted_store, member):
        store, _ = converted_store
        arr = zarr.open_group(str(store), mode="r")[f"tables/ds_z0/X/{member}"]

        # One shard is one file; with no sharding the chunk is the file.
        unit = arr.shards[0] if arr.shards else arr.chunks[0]
        file_bytes = unit * arr.dtype.itemsize
        array_bytes = arr.shape[0] * arr.dtype.itemsize

        # An array smaller than the floor can only ever be one short file.
        assert file_bytes >= min(array_bytes, MIN_TABLE_SHARD_BYTES), (
            f"X/{member}: {file_bytes} B per file over {array_bytes} B of data "
            f"-> {int(np.ceil(arr.shape[0] / unit)):,} files. Below "
            f"{MIN_TABLE_SHARD_BYTES} B this is the zarr < 3.1.4 "
            "auto-shard defect; check the zarr floor in pyproject.toml."
        )

    def test_data_array_is_large_enough_for_that_to_mean_anything(
        self, converted_store
    ):
        # Guards the fixture, not the code. The assertion above is
        # `file_bytes >= min(array_bytes, MIN_TABLE_SHARD_BYTES)`, so it says
        # nothing at all unless X/data is bigger than the floor -- which it
        # stops being if the converter ever drops these spectra to a handful of
        # non-zeros. Today it is ~400 KB.
        store, _ = converted_store
        data = zarr.open_group(str(store), mode="r")["tables/ds_z0/X/data"]
        assert data.shape[0] * data.dtype.itemsize > MIN_TABLE_SHARD_BYTES

    def test_whole_table_stays_in_a_handful_of_files(self, converted_store):
        store, _ = converted_store
        x_dir = store / "tables" / "ds_z0" / "X"
        n_files = sum(1 for p in x_dir.rglob("*") if p.is_file())
        # Fixed this is single digits (each array fits one shard). Broken, the
        # same table is ~400 files: X/data alone shards at 200 elements.
        assert n_files < 100, f"{n_files} files under tables/ds_z0/X"
