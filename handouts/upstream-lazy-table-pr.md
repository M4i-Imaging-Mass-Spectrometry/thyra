# E. Upstream: land spatialdata PR #1055 (lazy table loading)

**Repo:** `scverse/spatialdata`, not this one
**PR:** https://github.com/scverse/spatialdata/pull/1055
**Branch:** `Tomatokeftes:spatialdata:feature/lazy-table-loading`
**Head:** `4b1da50`, rebased onto `main` (`eb4fb3d`) on 2026-07-30
**Priority:** highest leverage for Ousia. Independent of everything else here.

This is the only handout whose work happens outside the Thyra repo. It is
included because it is on the critical path for Ousia consuming Thyra
output lazily, and because its state changes how handouts B and C should be
judged.

---

## State as of 2026-07-30 (after the rebase)

| | |
|---|---|
| Status | **open, mergeable, green, zero human reviews** |
| Opened | 2026-01-27 (six months) |
| Last touched | 2026-07-30 — rebased, fixed, review requested |
| Size | 7 files, +345 / -63, in 5 commits (squashed from 13) |
| Divergence | 5 ahead, **0 behind** `main` |
| CI | all 10 checks pass |
| Comments | the codecov bot, plus your review request |

The 10 green checks are ubuntu / macOS / Windows across Python 3.12, 3.13
and 3.14, the min-dask job, pre-commit.ci, Read the Docs, and codecov.

Nothing has been requested by a maintainer. It is not blocked on changes; it
is simply unattended. What has changed is that there is no longer any
housekeeping reason to leave it alone: it is rebased, green, documented, and
carries a regression test.

## What it does

Adds `lazy: bool = False` to `SpatialData.read()` / `read_zarr()` and
`_read_table()`, backing tables with `anndata.experimental.read_lazy()`
instead of loading them into memory. Plus a `_is_lazy_anndata()` helper,
validation that skips eager checks for lazy tables, and fixes to
`_filter_table_by_element_names`, `_filter_table_by_elements`,
`get_values`, and `_inplace_fix_subset_categorical_obs` so
`bounding_box_query` and `aggregate` work on lazy tables.

Reported benchmark, 100,000 pixels x 100,000 m/z bins, ~296M non-zeros:
15.4 MB versus 2,270.7 MB, 0.13s versus 1.57s.

---

## Found while rebasing: the PR was silently broken

Upstream #1131 rewrote the join helpers to call `reset_index()` and
`groupby()` directly on `table.obs`. For a lazily-read table, `obs` is an
xarray `Dataset2D`, not a `DataFrame`, so those calls raise
`AttributeError`. That took out `bounding_box_query`, `get_values`, and all
five `join_spatialelement_table` modes — the exact paths the PR exists to
support.

This was reproduced on the **pre-rebase** head `bfd2b5f`, so it was
pre-existing rather than rebase-induced: the PR had been quietly broken by
upstream movement for some time, and nothing in CI or the review queue
surfaced it.

The PR now carries the fix: an `_obs_as_dataframe()` helper in
`relational_query.py`, called at the five sites that touch `obs`, plus
computing `X` in `get_values`. A regression test
(`test_lazy_table_relational_queries_match_eager`) asserts lazy and eager
relational queries agree, so the next upstream refactor of these helpers
fails loudly instead of silently.

If a maintainer asks why a single-feature PR carries a bug fix, that is the
answer.

---

## Verified: Thyra output already works with it

This was tested directly, against the pre-rebase head (`bfd2b5f`) and
pointed at stores produced by every Thyra write path:

```
in_memory      X=Array  obs=Dataset2D  block_equal=True  nnz=12  var=60
streaming_pcs  X=Array  obs=Dataset2D  block_equal=True  nnz=12  var=60
streaming_coo  X=Array  obs=Dataset2D  block_equal=True  nnz=12  var=60
```

`block_equal=True` means a lazy `X[:6, :30].compute()` is identical to the
eager read of the same block, and `obs["x"]` resolves to the right values.
That result is about the on-disk format and still holds; the relational-query
breakage above was in spatialdata's own join helpers, not in anything Thyra
writes.

**So no change is needed in Thyra for the lazy path to work.** The
`encoding-type` / `encoding-version` attrs `main` already writes are
sufficient. This also confirms the Thyra branch
`feature/lazy-loading-support` is genuinely redundant, not merely stale.

### One practical caveat worth knowing

Dask reductions over the sparse `X` do not work — and this is broader than
`sum`. All four common reductions fail, each in its own way:

```python
lazy.X.sum()    # TypeError: _cs_matrix.sum() got an unexpected keyword argument 'keepdims'
lazy.X.mean()   # IndexError: Index dimension must be 1 or 2
lazy.X.std()    # IndexError: Index dimension must be 1 or 2
lazy.X.max()    # TypeError about 'ndmin'
```

Note the absence of `.compute()` in those lines: they raise while dask is
**building the graph**, not at compute time. Dask derives a reduction's
result metadata by calling the matching NumPy reduction on one
`scipy.sparse` block, and `scipy.sparse` rejects the `keepdims` / `ndmin`
arguments NumPy forwards. Confirmed on both dask 2026.7.1 / scipy 1.18.0 and
dask 2026.1.1 / scipy 1.16.0, so it is not a single-version accident.

This is a dask/scipy-sparse interop limitation, not something the PR or
Thyra causes, and not something either can fix. Slicing and `.compute()`
work fine, which is the realistic MSI access pattern (extract an ion image,
pull one pixel's spectrum). Ousia should avoid whole-array dask reductions
on lazy tables, reduce a materialized slice (`X[:1000].compute().sum()`), or
use `map_blocks` with a sparse-aware function.

The limitation is now documented in the docstrings of `read_zarr()` and
`SpatialData.read()`, with those workarounds spelled out.

---

## The real blocker for Ousia: dependency ceilings

Even once this merges, **Thyra users cannot reach it** without Thyra moving
three pins. Current spatialdata `main` declares:

```
dask>=2026.3.0
ome_zarr>=0.16.0
distributed>=2026.3.0
```

Those floors come from upstream `main`, not from this PR — the PR touches no
`pyproject.toml` at all. They are the ambient cost of tracking spatialdata,
and they would apply to any release carrying #1055.

Thyra's `pyproject.toml` pins:

```
dask         = ">=2025.12.0, <2026.2"   # installed: 2026.1.1
ome-zarr     = ">=0.14.0, <0.16"        # installed: 0.15.0
spatialdata  = ">=0.7.3, <0.8"          # installed: 0.7.3
```

spatialdata is **already at v0.8.0** (released 2026-07-02), which Thyra's
ceiling excludes. Whatever release carries #1055 will be further out still.

Those ceilings are deliberate and well-commented — there is a long note in
`pyproject.toml` explaining that this cluster releases frequently with
breaking cross-version interactions, and that the clean-venv CI job exists
to keep them honest. They are not accidental and should not be bumped
casually.

Encouraging data point: the PR branch **imported and ran correctly against
Thyra's installed dask 2026.1.1 and ome_zarr 0.15.0**, below the declared
floors. So the floors may be conservative and the real compatibility window
may be wider than the metadata claims. Worth establishing rather than
assuming.

---

## Suggested work

Three of the four items below were done on 2026-07-30 and are recorded here
so the reasoning survives rather than being deleted.

1. ~~**Rebase onto `main`.**~~ Done. 13 commits squashed to 5, replayed onto
   `eb4fb3d`, force-pushed to `4b1da50`. All 10 checks green.
2. ~~**Ping for review.**~~ Done. A comment with the benchmark and Thyra
   named as the concrete downstream consumer (MSI data, for Ousia) is
   posted. Still zero human reviews — the scverse Zulip remains the next
   escalation if it stays quiet.
3. ~~**Cover the caveat.**~~ Done. The sparse-reduction limitation is
   documented in the docstrings of `read_zarr()` and `SpatialData.read()`,
   including which reductions fail, that they fail at graph-build time, and
   what to do instead.
4. **Once it lands**, open a follow-up in Thyra to raise the
   `spatialdata` / `dask` / `ome-zarr` ceilings together, exercised by the
   existing clean-venv CI job. That is the step that actually delivers lazy
   reading to Ousia. **This is the only item still open here**, and it is
   gated on a maintainer, not on you.

## Related upstream gap you are positioned to fix

`spatialdata-io` issue #364 reports `spatialdata.write` failing on Xenium
data under pandas 3.0, with the same `IORegistryError` on
`ArrowStringArray` that handout A works around inside Thyra. anndata #2221
tracks the root cause and is milestoned 0.14.0.

Handout A fixes this downstream in Thyra because Thyra cannot wait. But the
same coercion belongs in spatialdata's own table write path, and you already
have commit rights to a PR there. If you fix it upstream, Thyra's workaround
eventually deletes itself. Worth raising on #364 at least.
