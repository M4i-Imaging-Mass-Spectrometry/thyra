# E. Upstream: land spatialdata PR #1055 (lazy table loading)

**Repo:** `scverse/spatialdata`, not this one
**PR:** https://github.com/scverse/spatialdata/pull/1055
**Branch:** `Tomatokeftes:spatialdata:feature/lazy-table-loading`
**Priority:** highest leverage for Ousia. Independent of everything else here.

This is the only handout whose work happens outside the Thyra repo. It is
included because it is on the critical path for Ousia consuming Thyra
output lazily, and because its state changes how handouts B and C should be
judged.

---

## State as of 2026-07-30

| | |
|---|---|
| Status | **open, mergeable, zero human reviews** |
| Opened | 2026-01-27 (six months) |
| Last touched | 2026-07-13 |
| Size | 8 files, +208 / -58 |
| Divergence | 13 ahead, **5 behind** `main` |
| Patch coverage | 96.55% (codecov, 1 line uncovered) |
| Only comment | the codecov bot |

Nothing has been requested by a maintainer. It is not blocked on changes; it
is simply unattended. Five commits behind is nothing — this is very
rebasable.

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

## Verified: Thyra output already works with it

This was tested directly, by checking out the PR head
(`bfd2b5f`) and pointing it at stores produced by every Thyra write path:

```
in_memory      X=Array  obs=Dataset2D  block_equal=True  nnz=12  var=60
streaming_pcs  X=Array  obs=Dataset2D  block_equal=True  nnz=12  var=60
streaming_coo  X=Array  obs=Dataset2D  block_equal=True  nnz=12  var=60
```

`block_equal=True` means a lazy `X[:6, :30].compute()` is identical to the
eager read of the same block, and `obs["x"]` resolves to the right values.

**So no change is needed in Thyra for the lazy path to work.** The
`encoding-type` / `encoding-version` attrs `main` already writes are
sufficient. This also confirms the Thyra branch
`feature/lazy-loading-support` is genuinely redundant, not merely stale.

### One practical caveat worth knowing

Dask reductions over the sparse `X` do not work:

```python
lazy.X.sum().compute()
# TypeError: _cs_matrix.sum() got an unexpected keyword argument 'keepdims'
```

Dask passes `keepdims` into `scipy.sparse`'s `sum`, which does not accept
it. This is a dask/scipy-sparse interop limitation, not something the PR or
Thyra causes, and not something either can fix. Slicing and `.compute()`
work fine, which is the realistic MSI access pattern (extract an ion image,
pull one pixel's spectrum). Ousia should avoid whole-array dask reductions
on lazy tables, or call `.to_memory()` on a slice first.

---

## The real blocker for Ousia: dependency ceilings

Even once this merges, **Thyra users cannot reach it** without Thyra moving
three pins. The PR branch declares:

```
dask>=2026.3.0
ome_zarr>=0.16.0
distributed>=2026.3.0
```

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
Thyra's installed dask 2026.1.1 and ome_zarr 0.15.0**, below its declared
floors. So the floors may be conservative and the real compatibility window
may be wider than the metadata claims. Worth establishing rather than
assuming.

---

## Suggested work

1. **Rebase onto `main`** (only 5 commits behind) so it is trivially
   mergeable and CI is green against current `main`.
2. **Ping for review.** Six months with no maintainer response usually means
   it fell off the queue rather than that anyone objects. A comment with the
   benchmark table and a concrete downstream consumer named (MSI data via
   Thyra, for Ousia) is the useful nudge. Consider the scverse Zulip.
3. **Cover the caveat.** A note in the docstring or docs that dask
   reductions over sparse `X` are unsupported would save the next person the
   confusion, and pre-empts a "lazy loading is broken" issue.
4. **Once it lands**, open a follow-up in Thyra to raise the
   `spatialdata` / `dask` / `ome-zarr` ceilings together, exercised by the
   existing clean-venv CI job. That is the step that actually delivers lazy
   reading to Ousia.

## Related upstream gap you are positioned to fix

`spatialdata-io` issue #364 reports `spatialdata.write` failing on Xenium
data under pandas 3.0, with the same `IORegistryError` on
`ArrowStringArray` that handout A works around inside Thyra. anndata #2221
tracks the root cause and is milestoned 0.14.0.

Handout A fixes this downstream in Thyra because Thyra cannot wait. But the
same coercion belongs in spatialdata's own table write path, and you already
have commit rights to a PR there. If you fix it upstream, Thyra's workaround
eventually deletes itself. Worth raising on #364 at least.
