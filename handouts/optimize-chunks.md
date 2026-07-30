# B. `--optimize-chunks` has never worked

> **RESOLVED** -- Option 2 (remove) was taken. `optimize_zarr_chunks` and
> `thyra/utils/data_processors.py` are gone; the flag survives as a hidden,
> deprecated no-op that warns, so scripts passing it keep running.
>
> Three corrections to the brief below, found while investigating:
>
> * `_chunking.py` covers **raster images only**, not tables -- but the real
>   finding is that *four* writers produce `tables/<id>/X` and only the
>   streaming-CSC one sets chunks deliberately; the rest inherit zarr's
>   `_guess_chunks`. Both layouts are sensible, which is why removal won.
> * The function was broken on **dense arrays too**: it reads
>   `array.compressor`, which zarr 3 removed.
> * A third defect, unlisted below: the CLI built `tables/<dataset_id>/X`, but
>   the 2D and streaming converters write `tables/<dataset_id>_z0/X`, so via the
>   CLI the key never resolved. (The 3D converter *does* use the bare id, which
>   is the one path where the documented `.shape` error is reachable.)
> * The attrs claim at line ~123 is wrong: `encoding-type`/`encoding-version` on
>   `data`/`indices`/`indptr` are **not** load-bearing and are absent on three of
>   four write paths. The load-bearing attrs are on the `X` group and the table
>   group. See `tests/unit/test_read_lazy_contract.py`.

**Branch:** `fix/optimize-chunks-sparse`
**Worktree:** `../Thyra-chunks`
**Priority:** independent, shares no files with the other handouts.

```bash
git worktree add -b fix/optimize-chunks-sparse ../Thyra-chunks main
```

---

## The defect

The documented `--optimize-chunks` CLI flag fails on every conversion.

Reproduce on `main` at v1.27.0, with any imzML and a perfectly ordinary
short output path:

```python
from thyra.convert import convert_msi
from thyra.utils.data_processors import optimize_zarr_chunks

convert_msi(imzml, out, dataset_id="ds", pixel_size_um=2.5)
optimize_zarr_chunks(str(out), "tables/ds_z0/X")
```

```
ERROR Error optimizing chunks: 'Group' object has no attribute 'shape'
-> returns False
```

### Why

`tables/<dataset_id>/X` is a Zarr **Group**, not an Array. It holds a SciPy
sparse matrix:

```
encoding-type attr: csc_matrix
members: ['data', 'indices', 'indptr']
```

`optimize_zarr_chunks` (`thyra/utils/data_processors.py:38-39`) does
`zarr_store[array_path]` and then reads `array.shape`. Its chunk-selection
logic branches on `len(shape) == 4` for a `(c, z, y, x)` layout, which says
plainly that it was written for a dense image layout predating the sparse
CSC tables the converter now writes. It is stale code from an earlier
architecture.

### A second, separate defect

`thyra/__main__.py::_handle_post_conversion` calls it and **ignores the
return value**:

```python
if success:
    if optimize_chunks and format == "spatialdata":
        optimize_zarr_chunks(str(output), f"tables/{dataset_id}/X")
    logger.info(f"Conversion completed successfully. Output stored at {output}")
```

So `--optimize-chunks` logs an error, does nothing, and the CLI still
reports success and exits 0. This is the same shape of bug as the exit-code
defect fixed in `1584f10`, one layer down: a step fails and nothing
propagates.

---

## What to do

Decide between fixing and removing. Investigate before choosing.

### Option 1: make it work on the sparse layout

`data`, `indices`, and `indptr` are 1-D, so the existing 4-D and 2-D
branches do not apply at all — the whole chunk-selection body needs
rewriting, not patching. Rechunking a CSC store means picking chunk sizes
for three 1-D arrays whose lengths are `nnz`, `nnz`, and `n_cols + 1`.

Before writing any of it, answer: **does this actually buy anything?**
Compare against what `thyra/converters/spatialdata/_chunking.py` already
does. That module was introduced as "the sharding seam" and may already set
sensible chunks at write time, making a post-hoc pass pointless.

There is a live reason it might buy something. spatialdata PR #1055 — see
[handout E](upstream-lazy-table-pr.md) — adds `read_zarr(..., lazy=True)`,
and it already works on Thyra output. Under lazy reading the chunk sizes of
`data`, `indices`, and `indptr` **directly determine read performance**,
because a lazy consumer pulls chunks rather than the whole array. So
"chunking is right for the access pattern" stops being cosmetic.

If you take Option 1, benchmark against the lazy access pattern, not just
the eager one: an ion image over an m/z window, and a single pixel's
spectrum, via `read_zarr(..., lazy=True)`. That is how Ousia will read
these stores. Note that whole-array dask reductions over sparse `X` do not
work at all (a dask/scipy-sparse limitation described in handout E), so do
not build a benchmark around `X.sum()`.

### Option 2: remove the flag and the function

If `_chunking.py` already covers it, deleting is the better answer: a
documented flag that silently does nothing is worse than no flag. Remove
`--optimize-chunks` from the CLI, drop `optimize_zarr_chunks` and its tests,
and update `docs/cli.md` (Performance section) and `docs/getting-started.md`
if it appears there.

Removing a public CLI flag is a breaking change for anyone scripting it.
Since it has never done anything, nobody can be depending on its effect —
but a script passing the flag would start failing on an unknown option.
Consider accepting it as a deprecated no-op that warns, or note the break
clearly in the commit message. Your call; state the reasoning.

### Either way

`_handle_post_conversion` must not report success when the step failed.

---

## Constraints

**Do not break the AnnData encoding attributes.** This is the important one.

`main` writes `encoding-type` / `encoding-version` attrs throughout the
store, and they are what makes `anndata.experimental.read_lazy()` work — the
lazy path the Ousia software depends on. Verified working on `main`:

```
read_lazy OK: 36 x 40, X=Array
  obs columns: ['y', 'region', 'region_number', 'x', 'spatial_y', 'spatial_x']
```

Relevant attrs on the table's `X` group and its members:

```
X_group.attrs["encoding-type"]  = "csc_matrix"   (or "csr_matrix")
X_group.attrs["encoding-version"]
data/indices/indptr each: encoding-type = "array", encoding-version = "0.2.0"
```

If you rewrite those arrays to rechunk them, **the attrs must survive**.
`zarr` will not carry them across a create-and-copy. A regression here would
silently break lazy reading, which is exactly the class of change to avoid.

Add a test asserting `read_lazy` still works after whatever you do, not just
that the arrays exist.

**If you keep the flag, it needs the extended-length path treatment.** The
CLI passes the user's unprefixed path, while `convert_msi` may have written
through a `\\?\` path — see `prepare_zarr_output_path` in
`thyra/utils/windows_paths.py`. Opening deep keys inside a long-path store
with a plain path will fail on Windows. Use the same helper.

## Verification

```bash
cd ../Thyra-chunks
PYTHONPATH=$(pwd) poetry run pytest -q
poetry run black . && poetry run isort . && poetry run flake8
poetry run mkdocs build --strict
```

Plus, whichever option you take, a test that fails on `main`:

- If fixed: `optimize_zarr_chunks` returns `True` on a real converted store,
  chunk shapes actually changed, and `read_lazy` still works.
- If removed: the CLI no longer advertises the flag, and
  `_handle_post_conversion` surfaces failures.
