# C. Write amplification and resampling memory

**Branch:** `perf/write-amplification`
**Worktree:** `../Thyra-perf`
**Priority:** start after handout A merges, then branch from the updated
`main`. A touches `_save_output` in the same file you will be working in.

```bash
git fetch origin && git worktree add -b perf/write-amplification ../Thyra-perf origin/main
```

---

## Framing

This is an **investigation**, not a list of fixes. Thyra's stated purpose is
converting 100+ GB datasets, and there are four measured or observed leads
suggesting the conversion path does substantially more work than it needs
to. Measure each at scale first. Report numbers before proposing changes,
and do not assume a lead is real because it is written down here — two of
the four are observations rather than measurements, and are marked as such.

There are existing benchmarks in `benchmarks/`. Check what they already
cover before writing new ones.

---

## Lead 1: metadata rewritten dozens of times per conversion (measured)

Instrumenting `zarr.storage._local._atomic_write` during **one conversion of
a 4-spectrum dataset** showed 40 distinct Zarr keys written but 22 of them
written more than once:

| key | times written |
|---|---|
| `<store>/zarr.json` | **28** |
| `<store>/tables/zarr.json` | 21 |
| `<store>/tables/ds_z0/zarr.json` | 20 |
| `<store>/tables/ds_z0/obs/zarr.json` | 15 |
| `<store>/images/ds_z0_tic/zarr.json` | 8 |
| `<store>/tables/ds_z0/X/zarr.json` | 7 |

Roughly 430 metadata renames total, for four pixels.

This is not cosmetic. It is the direct cause of the Windows `WinError 5`
race that v1.27.0 works around with a bounded retry — see the analysis in
`thyra/utils/zarr_atomic_write.py`. Reducing the rewriting would remove the
cause rather than retrying past it.

Find out:

- Is this `SpatialData.write()` re-serialising the store root once per
  element added? Compare against `write_element` used incrementally.
- **How does it scale?** The suspicion is that it scales with *element
  count*, not pixel count, which would make it a fixed overhead rather than
  a large-dataset problem. Confirm or kill that. If it is fixed overhead of
  a few hundred renames, it matters for robustness but not for throughput,
  and the finding is "keep the retry, close the issue".
- Whether the streaming converter's own `zarr.open_group` /
  `consolidate_metadata` calls add to it.

## Lead 2: default bin counts upsample heavily (measured)

The default is 5 mDa at m/z 1000 for every axis type except `linear_tof`
(17 mDa at m/z 300). On the tutorial's 250-1200 Da synthetic dataset that is
**190,000 bins from 4,000 source m/z values**, a 45x upsample.

The store stays small (~29 MB) because it is sparse, so this is not a disk
problem. It costs every consumer that materialises the table, and it feeds
lead 4 below.

See `_calculate_bins_from_width` in
`thyra/converters/spatialdata/base_spatialdata_converter.py` and the
bin-count table in `docs/resampling.md`.

**Treat this lead sceptically.** It looked stronger before the upstream
lazy-loading work was taken into account. spatialdata PR #1055 — see
[handout E](upstream-lazy-table-pr.md) — benchmarks itself at **100,000
pixels x 100,000 m/z bins** and works on Thyra output today. Wide `var`
axes are the design target of the lazy path, not a pathology it struggles
with. If the intended consumption mode is `read_zarr(..., lazy=True)`, a
190k-wide axis is largely free to the consumer and shrinking the default
buys much less than it appears to.

So the question to answer is not "are 190,000 bins too many" but:

- who actually pays for the width, once lazy reading is in play? Measure a
  realistic Ousia-style access (ion image over an m/z window, single-pixel
  spectrum) against a wide and a narrow axis.
- is 45x upsampling scientifically meaningful, or is it inventing
  resolution the instrument never had? That is a different argument from
  the performance one and may be the stronger reason to change.

Deriving the default from observed source spacing — `peak density` is
already collected in `DataCharacteristics` — remains a reasonable idea. But
justify it on the numbers above, not on the raw bin count looking large.

**This is a behaviour change with compatibility impact.** Thyra is a library
with an external consumer (Ousia). Changing default bin counts changes the
`var` axis of every future conversion. Do not change it quietly: if the
numbers support a change, propose it with before/after figures and let the
owner decide.

## Lead 3: processed-mode imzML scans every spectrum twice (observed)

With `--no-resample` on processed data, the raw-axis path iterates all
spectra to collect unique m/z values. The code already says so:

```
Building RAW mass axis (no resampling) - processed mode, iterating ALL
spectra to collect unique m/z values. This is slow for large datasets!
```

There is also a separate "Scanning mass range and counting peaks" pass
before it.

Measure both passes on something large. Then see whether they can be merged
into one, and whether the unique-m/z collection can use something cheaper
than a Python set over every value.

## Lead 4: the dense resampling path allocates per pixel (observed)

`_resample_spectrum_to_indices` returns `np.arange(len(common_mass_axis))`
for **every spectrum**, and the non-nearest-neighbor branch of
`_process_single_spectrum` builds a dense array of the full axis width per
pixel (`base_spatialdata_converter.py`, roughly lines 890-970).

`nearest_neighbor` has a sparse fast path that returns only non-zero bins.
`tic_preserving` does not. At 190k bins that is ~1.5 MB allocated per pixel
before sparsification, so a 10,000-pixel dataset with `tic_preserving` does
15 GB of allocation churn.

Check whether `tic_preserving` can produce sparse output directly. Note the
constraint that makes it non-trivial: TIC preservation is defined over the
whole spectrum, so the rescaling factor needs the full interpolated
spectrum's sum. That may be computable without materialising the dense
array.

---

## Constraints

**Do not change the on-disk format in a way that breaks lazy reading.**
`anndata.experimental.read_lazy()` currently works on Thyra output and Ousia
depends on it. Verified on `main`:

```
read_lazy OK: 36 x 40, X=Array
  obs columns: ['y', 'region', 'region_number', 'x', 'spatial_y', 'spatial_x']
```

That works because the streaming PCS path hand-writes `encoding-type` /
`encoding-version` attrs for every array, group, categorical, and 0-d
scalar. `main` is deliberately granular here: `"string"` for scalars in
`uns` versus `"string-array"` for arrays. If you touch how the store is
written — and lead 1 is squarely about that — **assert `read_lazy` still
works**, not merely that the files exist.

**Do not "fix" the WinError 5 retry by removing it** even if lead 1 reduces
the rewriting. Fewer renames lowers the probability of the race; it does not
eliminate it, since the race needs only one contended rename.

**Coordinate systems are contract-tested.** `tests/unit/converters/
test_coordinate_systems.py` pins that the TIC image and pixel polygons agree
at `"global"`. If a write-path change moves transforms, that is a real
regression, not a stale test.

## Verification

```bash
cd ../Thyra-perf
PYTHONPATH=$(pwd) poetry run pytest -q
poetry run black . && poetry run isort . && poetry run flake8
poetry run mkdocs build --strict
```

For lead 1, the instrumentation that produced the numbers above is worth
recreating: wrap `zarr.storage._local._atomic_write`, count writes per key
path, and report the distribution. Be aware that instrumenting it perturbs
timing enough to change the `WinError 5` failure rate substantially, so do
not use write counts and race frequency from the same run to draw
conclusions about both.

## Deliverable

A findings report with measurements, then separate commits per change that
the measurements justify. It is a perfectly good outcome for this handout to
conclude that two of the four leads are not worth acting on — say so with
the numbers that show it.
