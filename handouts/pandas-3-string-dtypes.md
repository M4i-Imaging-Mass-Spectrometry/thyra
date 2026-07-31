# A. pandas 3.0 string dtypes break SpatialData writes

**Branch:** `fix/pandas3-string-dtypes`
**Worktree:** `../Thyra-pandas3`
**Priority:** land this first. It is the only handout fixing something a
user can hit today, and handout C waits on it.

```bash
git worktree add -b fix/pandas3-string-dtypes ../Thyra-pandas3 main
```

---

## The defect

Under pandas 3.0 — or pandas 2.x with `future.infer_string=True`, which is
what 3.0 makes the default — writing a converted dataset fails:

```
Error saving SpatialData: No method registered for writing
<class 'pandas.core.arrays.string_arrow.ArrowStringArrayNumpySemantics'>
into <class 'zarr.core.group.Group'>
```

`convert_msi` catches it and returns `False`. Since v1.27.0 the CLI exits 1
rather than silently reporting success, so this now fails loudly — but it
still fails.

### Reproduced on main at v1.27.0

Setting `pd.set_option("future.infer_string", True)` before converting:

| Path | Result |
|---|---|
| in-memory (`streaming=False`) | **FAIL** |
| streaming PCS (`streaming=True, use_csc=True`) | PASS |
| streaming COO (`streaming=True, use_csc=False`) | **FAIL** |

The PCS path survives because it hand-writes the AnnData layout straight to
Zarr and never goes through anndata's writer. Everything that reaches
`SpatialData.write()` fails.

Note this is **not** COO-specific, contrary to what the existing branch's
commit message implies. The in-memory path — the default for any dataset
under the 10 GB streaming threshold, i.e. most conversions — fails too.

### Why it happens

With `infer_string=True`, the table's obs index (`instance_id`), the
`instance_key` column, and the var index (`mz_*`) are inferred as pandas'
`str` dtype, backed by `ArrowStringArrayNumpySemantics`. anndata's IO
registry has a writer for `pandas.core.arrays.string_.StringArray` but not
for the `*NumpySemantics` subclasses pandas 3 actually produces, and the
registry matches exact types rather than subclasses.

### Upstream has already fixed this, in a release the ceiling excludes

**Corrected 2026-07-30.** This section previously said anndata #2221 was open
and milestoned 0.14.0, and concluded that upstream would not fix this in a
usable timeframe. That was wrong. The fix below was written on that mistaken
basis, which does not make it wrong — it was and is needed on the pinned
anndata — but it does change it from a permanent fixture into removable work.

Checked directly:

- **anndata #2377** is our exact error. Closed as a duplicate of #2221.
- **anndata #2221** ("Pandas 3.0 compatibility") was **closed as completed**
  on 2025-12-11 by anndata #2133, and shipped in **anndata 0.13.0**
  (2026-07-07). Its milestone label still reads 0.14.0, which is what caused
  the original misreading — the label was never updated to match the release
  the fix actually landed in.
- **spatialdata-io #364** is the same failure hitting spatialdata's own
  Xenium writer, so this is an ecosystem gap, not a Thyra bug. Note the right
  advice for that issue is now "upgrade to anndata >= 0.13", not the coercion
  below.

Verified by direct test on 2026-07-30, writing an obs table with a string
index, a string column, and a string-backed categorical:

| pandas | pyarrow | anndata | `write_zarr` |
|---|---|---|---|
| 2.3.2 + `infer_string` | yes | 0.12.2 | `IORegistryError` on `ArrowStringArrayNumpySemantics` |
| 3.0.5 | yes | 0.12.2 | `IORegistryError` on `pandas.arrays.ArrowStringArray` |
| 3.0.5 | no | 0.12.2 | `RuntimeError` on `allow_write_nullable_strings` |
| 2.3.2 + `infer_string` | yes | **0.13.2** | **PASS** |
| 3.0.5 | yes | **0.13.2** | **PASS** |

The concrete failing class varies with pandas version and whether pyarrow is
installed, which is why a fix keyed to any one class name would be wrong.

`pyproject.toml` pins `anndata = ">=0.11.0, <0.13"`, which excludes the fix.
spatialdata 0.7.3 requires only `anndata>=0.9.1`, so raising that ceiling is
**not** blocked by the spatialdata pin, and not blocked by
[handout E](upstream-lazy-table-pr.md) either.

**anndata's documented escape hatches do not work on our pinned version.**
Verified against anndata 0.12.2 with `infer_string=True`:

| Attempt | Result |
|---|---|
| default | `IORegistryError` on `ArrowStringArrayNumpySemantics` |
| `ad.settings.allow_write_nullable_strings = True` | still `IORegistryError` |
| `pd.set_option("mode.string_storage", "python")` | still fails, now on `StringArrayNumpySemantics` |
| both of the above together | still fails, on `StringArrayNumpySemantics` |

That last row is worth stating explicitly, because **on pandas 3 proper the
combination does work.** Measured on pandas 3.0.5 + pyarrow 25.0.0 + anndata
0.12.2: `mode.string_storage="python"` together with
`allow_write_nullable_strings = True` writes cleanly, with values, index and
categorical dtype all intact. Neither setting alone is enough there either --
anndata 0.12.x has a `StringArray` writer gated behind that setting and no
`ArrowStringArray` writer at all, so you have to leave Arrow storage *and*
unlock nullable strings.

It does not help under `infer_string` on pandas 2, which is what CI runs,
because there the arrays are the `*NumpySemantics` subclasses that no writer
matches.

Either way it is the wrong tool for a library: both are global process state,
and Thyra should not flip pandas and anndata settings on its users' behalf. A
local coercion at the write chokepoint stays the right layer.

So a Thyra-side fix is the right option until the `anndata` ceiling moves to
`>= 0.13`.

### It is reachable today

`pyproject.toml` line 74 is `pandas = ">=2.0.0"` with **no upper bound**. A
fresh `pip install thyra` can resolve pandas 3.x and every non-PCS
conversion breaks.

---

## What to do

There is an existing commit that does the right thing:
`b23e636` on `origin/fix/streaming-coo-arrow-string-write`. It adds
`_coerce_table_strings_to_object()` and calls it in the base `_save_output()`
for every table's `obs` and `var` before `SpatialData.write()`.

**It cherry-picks cleanly onto current main** (verified) and **it works**:
with it applied, all three write paths pass under `infer_string=True`.

So:

1. `git cherry-pick b23e636`
2. Add the test coverage it lacks (below).
3. Decide the pandas bound (below).
4. Correct the commit message's claim that this is COO-specific.

### Why this is the right layer, not a patch

Worth stating in the commit message, because it is the crux:

- It sits at the single serialization chokepoint, so the 2D, 3D, and COO
  converters are all covered by one call.
- **It does not change the on-disk format.** Verified: after the coercion,
  `sd.read_zarr()` returns `obs.index.dtype == str` (i.e. pandas' native
  string dtype). The coercion is write-time only; readers get the modern
  dtype. Nothing about the stored SpatialData/AnnData layout differs from a
  pandas 2 write.
- It is a no-op on pandas < 3, where those columns are already `object`.

### Required: a removal trigger

This is dead code the moment the `anndata` ceiling moves to `>= 0.13`, which
already exists and is not blocked on anything. Leave something that makes that
obvious rather than letting it calcify:

- A module-level comment naming **anndata #2221** and the 0.13.0 release that
  fixed it, stating the coercion can be deleted when the `anndata` ceiling
  moves to `>=0.13`.
- Ideally a test that fails loudly when the installed anndata *can* write
  Arrow-backed strings, so the workaround announces its own obsolescence
  instead of silently pessimising. Something like: if
  `ad.io.write_elem` succeeds on a `str`-dtype frame, the coercion is no
  longer needed. Mark it `xfail(strict=False)` or skip-with-reason so it is
  informative rather than a CI blocker.

### Required: tests

The branch has none, which is the main thing wrong with it. The commit
message claims validation via
`test_streaming_pcs_roundtrip.py::test_coo_path_roundtrips` under pandas
3.0, but nothing in CI runs under `infer_string=True`, so the fix can
regress silently.

Add a test module that, per test, sets `future.infer_string=True` and
restores it afterwards (a fixture with `pd.option_context` is cleanest), and
covers:

- all three write paths (`streaming=False`, `streaming=True` with
  `use_csc=True` and `use_csc=False`) convert successfully;
- the store reads back via `sd.read_zarr()` with the expected shape;
- `_coerce_table_strings_to_object` leaves values unchanged — same order,
  same content, only dtype differs;
- a string-backed **categorical** (the `region` column) keeps its codes and
  category values;
- it is a no-op when dtypes are already `object`.

Give each conversion its own output path. My first attempt at this shared
one filename between the two streaming cases and the second run failed with
"Destination already exists", which looked like a real failure for a while.

### Decide: pandas upper bound

Two defensible options. Pick one and say why in the commit message.

- **Coercion only.** Thyra keeps working on pandas 3. Preferred if you want
  users on current pandas.
- **Coercion plus `pandas = ">=2.0.0, <4"`** or similar. The coercion
  handles the write; a bound documents what has actually been tested. Note
  the neighbouring dependency block already carries deliberate upper bounds
  with a long comment explaining why, so adding one here is consistent with
  house style.

Do **not** pin `pandas < 3` as the whole fix. It papers over a defect that
is already solved by a commit sitting in the repo, and it will strand users
whose other packages want pandas 3.

---

## Constraints

- Do not touch the PCS hand-written Zarr path in
  `thyra/converters/spatialdata/streaming_converter.py`. It already works
  and its `encoding-type` / `encoding-version` attrs are what makes
  `anndata.experimental.read_lazy()` work for Ousia. Breaking those breaks
  lazy reading.
- Do not "simplify" by building obs/var with `dtype=object` at construction
  instead. It would scatter the workaround across the 2D, 3D, and streaming
  converters instead of keeping it at one chokepoint, and it would fight
  pandas' inference in more places than it fixes.

## Verification

```bash
cd ../Thyra-pandas3
PYTHONPATH=$(pwd) poetry run pytest -q
PYTHONPATH=$(pwd) poetry run python -c "
import pandas as pd; pd.set_option('future.infer_string', True)
# then run a conversion on each of the three paths
"
poetry run black . && poetry run isort . && poetry run flake8
poetry run mkdocs build --strict
```

`PYTHONPATH=$(pwd)` matters — see the note in [README.md](README.md).

Once merged, delete the superseded branch and worktree:

```bash
git worktree remove ../Thyra-coo-fix
git push origin --delete fix/streaming-coo-arrow-string-write
```
