"""Foundation tests for the Thyra chunk/shard policy.

Two halves, for the two things Thyra writes:

**Rasters.** Locks the current chunk tuples (so the future zarr v3 sharding
wiring, gated on spatialdata PR #1106, cannot silently change output) and
documents the cross-repo INNER_TILE_EDGE contract with Ousia's serving tile
size.

**The MSI table.** Locks the shard size the table's ``X`` arrays land in. This
half exists because of a real, shipped regression: anndata >= 0.13 writes every
zarr v3 array with ``shards="auto"``, and zarr < 3.1.6 sized the auto chunk with
``max_bytes=1024`` where it meant 1 MiB (zarr-python#3603, first released in
3.1.6). On a zarr that also predates ``array.target_shard_size_bytes`` (< 3.1.4)
there is no budget to override that chunk, so it sizes the shard too. Nothing in
Thyra changed -- a dependency bump alone was enough to turn one table into
~400,000 files of roughly 1 KB each, and the zarr write into 98% of conversion
wall-clock (pea.imzML: 434 s of 445 s).

The shard, i.e. the file, is only half of it. 3.1.4 and 3.1.5 gained the budget
knob while keeping the bad ``max_bytes=1024``, so there the file count looks fine
and it is the chunk *inside* each file that collapses. Both halves are asserted
below, against ``MIN_TABLE_SHARD_BYTES`` and ``MIN_TABLE_CHUNK_BYTES``.

That failure is invisible to every other test in this suite. The store still
round-trips, ``read_lazy`` still opens it, and the values are all correct --
only the geometry is catastrophic. Hence a test that asserts geometry, on a
store big enough for the geometry to be non-degenerate: below ~1 KB of ``data``
the broken and fixed layouts agree, so a 4-pixel fixture would prove nothing.

**The storage profile of a table.** Locks which parts of Zarr a table uses:
format 3, little-endian, zstd, the sharding codec, variable-length strings, and
nothing else (design decision D91, and "Promised: how the tables are stored" in
docs/for-tool-developers.md). Thyra chooses none of these itself. They are the
defaults of zarr and anndata, so a dependency bump can change them with no
change in Thyra, and a reader written against a subset of Zarr would then stop
reading new stores. The check reads the ``zarr.json`` files as JSON, the way
such a reader does, and never opens the store with zarr.
"""

from __future__ import annotations

import contextlib
import json
from pathlib import Path

import numpy as np
import pytest
import zarr
from zarr.codecs import BloscCodec, BytesCodec, GzipCodec, ZstdCodec

from thyra.converters.spatialdata import _chunking
from thyra.converters.spatialdata._chunking import (
    IMAGE_CHUNK_EDGE,
    INNER_TILE_EDGE,
    MIN_TABLE_CHUNK_BYTES,
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
            f"{MIN_TABLE_SHARD_BYTES} B this is the pre-3.1.4 auto-shard "
            "collapse (no array.target_shard_size_bytes to override the "
            "1 KB auto chunk); check the zarr floor in pyproject.toml, "
            "which is >= 3.1.6."
        )

    @pytest.mark.parametrize("member", ["data", "indices", "indptr"])
    def test_x_arrays_are_not_kilobyte_chunked(self, converted_store, member):
        # The other half of the same defect, and the half the shard assertion
        # above is blind to: with Thyra's budget held open, zarr 3.1.4/3.1.5
        # still land a large shard and only the chunk inside it collapses.
        store, _ = converted_store
        arr = zarr.open_group(str(store), mode="r")[f"tables/ds_z0/X/{member}"]

        chunk_bytes = arr.chunks[0] * arr.dtype.itemsize
        array_bytes = arr.shape[0] * arr.dtype.itemsize

        # An array smaller than the floor can only ever be one short chunk.
        assert chunk_bytes >= min(array_bytes, MIN_TABLE_CHUNK_BYTES), (
            f"X/{member}: {chunk_bytes} B per chunk over {array_bytes} B of data "
            f"-> {int(np.ceil(arr.shape[0] / arr.chunks[0])):,} chunks per array, "
            f"each separately compressed and indexed. Below "
            f"{MIN_TABLE_CHUNK_BYTES} B this is the zarr < 3.1.6 auto-chunk "
            "defect (zarr-python#3603); check the zarr floor in pyproject.toml."
        )

    def test_data_array_is_large_enough_for_that_to_mean_anything(
        self, converted_store
    ):
        # Guards the fixture, not the code. Both assertions above are
        # `<bytes> >= min(array_bytes, <floor>)`, so they say nothing at all
        # unless X/data is bigger than the floor -- which it stops being if the
        # converter ever drops these spectra to a handful of non-zeros. Today it
        # is ~400 KB.
        store, _ = converted_store
        data = zarr.open_group(str(store), mode="r")["tables/ds_z0/X/data"]
        assert data.shape[0] * data.dtype.itemsize > max(
            MIN_TABLE_SHARD_BYTES, MIN_TABLE_CHUNK_BYTES
        )

    def test_whole_table_stays_in_a_handful_of_files(self, converted_store):
        store, _ = converted_store
        x_dir = store / "tables" / "ds_z0" / "X"
        n_files = sum(1 for p in x_dir.rglob("*") if p.is_file())
        # Fixed this is single digits (each array fits one shard). Broken, the
        # same table is ~400 files: X/data alone shards at 200 elements.
        assert n_files < 100, f"{n_files} files under tables/ds_z0/X"


# The storage profile a table is promised to stay inside (D91). The names are
# the Zarr v3 ones, as they appear in ``zarr.json``.
_PROFILE_DATA_TYPES = frozenset(
    {
        "bool",
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
        "float32",
        "float64",
        "string",
    }
)
_ONE_BYTE_DATA_TYPES = frozenset({"bool", "int8", "uint8"})
_ZARR_V2_METADATA = (".zarray", ".zgroup", ".zattrs", ".zmetadata")


def _codec(entry) -> tuple[str, dict]:
    """One codec of a ``zarr.json`` as (name, configuration)."""
    if isinstance(entry, str):
        return entry, {}
    return entry["name"], entry.get("configuration") or {}


def _chain_violations(codecs: list, data_type: str) -> list[str]:
    """What is wrong with the codecs that encode one chunk."""
    names = [_codec(entry)[0] for entry in codecs]
    first = "vlen-utf8" if data_type == "string" else "bytes"
    if names != [first, "zstd"]:
        return [f"codecs {names}, expected {[first, 'zstd']}"]
    endian = _codec(codecs[0])[1].get("endian")
    if first == "bytes" and data_type not in _ONE_BYTE_DATA_TYPES:
        if endian != "little":
            return [f"byte order {endian!r}, expected 'little'"]
    return []


def _shard_violations(configuration: dict) -> list[str]:
    """What is wrong with a sharding codec's own settings."""
    problems = []
    location = configuration.get("index_location", "end")
    if location != "end":
        problems.append(f"shard index at the {location}, expected the end")
    index = [_codec(entry) for entry in configuration["index_codecs"]]
    if [name for name, _ in index] != ["bytes", "crc32c"]:
        problems.append(f"shard index codecs {[name for name, _ in index]}")
    elif index[0][1].get("endian") != "little":
        problems.append("shard index is not little-endian")
    return problems


def _array_violations(meta: dict) -> list[str]:
    """What takes one array's ``zarr.json`` outside the profile."""
    problems = []
    data_type = meta["data_type"]
    # zarr writes some types (fixed-length bytes, datetimes, records) as a
    # JSON object, not a name; none of them is in the profile.
    if not isinstance(data_type, str) or data_type not in _PROFILE_DATA_TYPES:
        return [f"data type {data_type!r}"]
    if meta.get("storage_transformers"):
        problems.append(f"storage transformers {meta['storage_transformers']}")
    if meta["chunk_grid"]["name"] != "regular":
        problems.append(f"chunk grid {meta['chunk_grid']['name']!r}")
    key_name, key_config = _codec(meta["chunk_key_encoding"])
    if key_name != "default" or key_config.get("separator", "/") != "/":
        problems.append(f"chunk keys {meta['chunk_key_encoding']}")
    codecs = meta["codecs"]
    if [_codec(entry)[0] for entry in codecs] == ["sharding_indexed"]:
        configuration = _codec(codecs[0])[1]
        problems += _shard_violations(configuration)
        codecs = configuration["codecs"]
    return problems + _chain_violations(codecs, data_type)


def _profile_violations(root: Path) -> list[str]:
    """Every way the Zarr nodes under ``root`` leave the storage profile."""
    problems = [
        f"{path.relative_to(root).as_posix()}: Zarr 2 metadata file"
        for name in _ZARR_V2_METADATA
        for path in sorted(root.rglob(name))
    ]
    for path in sorted(root.rglob("zarr.json")):
        meta = json.loads(path.read_text(encoding="utf-8"))
        where = path.parent.relative_to(root).as_posix()
        if meta.get("zarr_format") != 3:
            problems.append(f"{where}: zarr_format {meta.get('zarr_format')!r}")
        elif meta.get("node_type") == "array":
            problems += [f"{where}: {text}" for text in _array_violations(meta)]
    return problems


def _table_arrays(store: Path) -> dict[str, dict]:
    """The ``zarr.json`` of every array under ``tables/``, by path."""
    root = store / "tables"
    found = {}
    for path in sorted(root.rglob("zarr.json")):
        meta = json.loads(path.read_text(encoding="utf-8"))
        if meta.get("node_type") == "array":
            found[path.parent.relative_to(root).as_posix()] = meta
    return found


@pytest.fixture(scope="module")
def mobility_store(tmp_path_factory) -> Path:
    """The committed mobility fixture, converted: a summed and a sibling table.

    It adds what ``converted_store`` has none of: a second table, and the one
    table array with two dimensions (``uns/mobility_heatmap/counts``, float32).
    """
    from thyra.convert import convert_msi

    source = Path(__file__).parents[2] / "data" / "fixtures"
    output_path = tmp_path_factory.mktemp("table_profile") / "out.zarr"
    assert convert_msi(
        str(source / "mobility_continuous.imzML"),
        str(output_path),
        format_type="spatialdata",
        dataset_id="ds",
    ), "conversion fixture failed"
    return output_path


_LEFT_THE_PROFILE = (
    "A table left the promised storage profile. If a dependency bump did this, "
    "a reader built on a subset of Zarr stops reading new stores: pin the "
    "dependency or change the profile on purpose (D91, "
    "docs/for-tool-developers.md).\n"
)

# One array's metadata, inside the profile, to break one field at a time.
_SHARDED_FLOAT64 = {
    "zarr_format": 3,
    "node_type": "array",
    "shape": [16],
    "data_type": "float64",
    "chunk_grid": {"name": "regular", "configuration": {"chunk_shape": [8]}},
    "chunk_key_encoding": {"name": "default", "configuration": {"separator": "/"}},
    "fill_value": 0.0,
    "storage_transformers": [],
    "codecs": [
        {
            "name": "sharding_indexed",
            "configuration": {
                "chunk_shape": [4],
                "codecs": [
                    {"name": "bytes", "configuration": {"endian": "little"}},
                    {"name": "zstd", "configuration": {"level": 0}},
                ],
                "index_codecs": [
                    {"name": "bytes", "configuration": {"endian": "little"}},
                    {"name": "crc32c"},
                ],
                "index_location": "end",
            },
        }
    ],
}


def _broken(path: tuple, value) -> dict:
    """``_SHARDED_FLOAT64`` with the field at ``path`` set to ``value``."""
    meta = json.loads(json.dumps(_SHARDED_FLOAT64))
    node = meta
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    return meta


class TestConvertedTableStorageProfile:
    """A table stays inside the storage profile the docs promise (D91)."""

    def test_every_table_array_is_inside_the_profile(self, converted_store):
        store, _ = converted_store
        problems = _profile_violations(store / "tables")
        assert not problems, _LEFT_THE_PROFILE + "\n".join(problems)

    def test_the_sibling_table_is_inside_the_profile_too(self, mobility_store):
        problems = _profile_violations(mobility_store / "tables")
        assert not problems, _LEFT_THE_PROFILE + "\n".join(problems)

    def test_the_stores_exercise_every_part_of_the_profile(
        self, converted_store, mobility_store
    ):
        # Guards the fixtures, not the code: the checks above only mean
        # something for the kinds of array the stores actually hold. Which
        # arrays are sharded is not promised, so none is named as sharded.
        store, _ = converted_store
        arrays = _table_arrays(store)
        sharded = {
            path
            for path, meta in arrays.items()
            if _codec(meta["codecs"][0])[0] == "sharding_indexed"
        }
        strings = {p for p, meta in arrays.items() if meta["data_type"] == "string"}
        wide = {
            p for p, m in arrays.items() if m["data_type"] not in _ONE_BYTE_DATA_TYPES
        }

        # The arrays docs/for-tool-developers.md names for a reader.
        named = ["X/data", "X/indices", "X/indptr", "var/mz", "obs/x", "obs/y"]
        assert {f"ds_z0/{name}" for name in named} <= set(arrays)
        assert sharded & (wide - strings), "no sharded numeric array to check"
        assert set(arrays) - sharded, "no unsharded array to check"
        assert strings, "no string array to check"

        siblings = _table_arrays(mobility_store)
        assert len({path.split("/")[0] for path in siblings}) == 2, "no sibling table"
        assert any(len(meta["shape"]) == 2 for meta in siblings.values())

    @pytest.mark.parametrize(
        ("path", "value", "expected"),
        [
            (("data_type",), {"name": "numpy.datetime64"}, "data type"),
            (("data_type",), "float16", "data type 'float16'"),
            (("storage_transformers",), [{"name": "x"}], "storage transformers"),
            (("chunk_grid", "name"), "rectilinear", "chunk grid 'rectilinear'"),
            (("chunk_key_encoding", "name"), "v2", "chunk keys"),
            (("chunk_key_encoding", "configuration", "separator"), ".", "chunk keys"),
            (
                ("codecs", 0, "configuration", "index_location"),
                "start",
                "shard index at the start",
            ),
            (
                ("codecs", 0, "configuration", "index_codecs"),
                [{"name": "bytes", "configuration": {"endian": "little"}}],
                "shard index codecs ['bytes']",
            ),
            (
                ("codecs", 0, "configuration", "codecs", 0, "configuration"),
                {"endian": "big"},
                "byte order 'big'",
            ),
            (
                ("codecs", 0, "configuration", "codecs"),
                [{"name": "transpose"}, {"name": "bytes"}, {"name": "zstd"}],
                "codecs ['transpose', 'bytes', 'zstd']",
            ),
        ],
    )
    def test_the_check_reports_each_field_it_reads(self, path, value, expected):
        assert _array_violations(_SHARDED_FLOAT64) == []
        problems = _array_violations(_broken(path, value))
        assert len(problems) == 1 and expected in problems[0], problems

    @pytest.mark.parametrize(
        ("kwargs", "expected"),
        [
            ({"compressors": GzipCodec()}, "codecs ['bytes', 'gzip']"),
            ({"compressors": BloscCodec()}, "codecs ['bytes', 'blosc']"),
            ({"compressors": None}, "codecs ['bytes']"),
            (
                {"serializer": BytesCodec(endian="big")},
                "byte order 'big'",
            ),
            ({"dtype": "complex128"}, "data type 'complex128'"),
            ({"zarr_format": 2, "compressors": None}, "Zarr 2 metadata file"),
            (
                {"shards": (8,), "compressors": GzipCodec()},
                "codecs ['bytes', 'gzip']",
            ),
        ],
    )
    def test_the_check_reports_an_array_outside_the_profile(
        self, tmp_path, kwargs, expected
    ):
        # The check is the whole guarantee, so it must not be able to pass on
        # the layouts a changed upstream default would produce.
        settings = {"dtype": "float64", "compressors": ZstdCodec()}
        settings.update(kwargs)
        zarr.create_array(
            store=str(tmp_path / "a"), shape=(16,), chunks=(4,), **settings
        )
        problems = _profile_violations(tmp_path)
        # A Zarr 2 array is two metadata files, so two lines; the rest are one.
        assert problems and all(expected in line for line in problems), problems

    @pytest.mark.parametrize("shards", [None, (8,)])
    @pytest.mark.parametrize("dtype", ["float64", "int8", "str"])
    def test_the_check_accepts_what_zarr_writes_by_default(
        self, tmp_path, shards, dtype
    ):
        zarr.create_array(
            store=str(tmp_path / "a"),
            shape=(16,),
            chunks=(4,),
            shards=shards,
            dtype=dtype,
        )
        assert _profile_violations(tmp_path) == []
