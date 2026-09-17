# Closed-form bin index: the measurements behind D21

**Issue:** thyra#295. **Decision:** `docs/design-decisions.md`, D21.
**Verdict:** GO, implemented on `perf/closed-form-bin-index`.
**Measured:** 2026-09-15 to 2026-09-16, one machine, 268 conversions,
12.8 hours of wall clock, every one of them bit-identical to the current
tree. This note holds what the decision page summarises: the derivation
that settled correctness, the harness, the environment, and every number,
so nobody re-derives them.

## The question

`_nn_map_to_bins` maps each peak of a spectrum to its nearest bin on the
resampled axis by `np.searchsorted` plus a left/right comparison. Six call
sites share it: two in the converter (`_nearest_neighbor_resample`,
`_build_nn_shared_cache`) and four in the sibling tables (the
mass-mobility heatmap, the mobility grid's discovery and scatter, the
MS/MS split). The issue proposed computing the index instead -- every
generator lays its bins uniformly in some analytic coordinate -- and a
first challenge (the issue's first comment, 2026-09-14) found the proposal
sound but its numbers wrong on three counts, and left three questions
open. The validation brief that followed asked for:

- **Q1 (blocking):** is the +-1 repair provably sufficient? The physics
  generators report m/z midpoints of edges laid uniformly in `u = f(m/z)`,
  so `f(centre_i) != u0 + i du`; the empirical worst deviation, 0.478 of a
  bin, sat 4.5 percent from the half-bin cliff past which the repair
  returns a wrong bin silently.
- **Q2:** does the win hold across the corpus and the laws, cold and warm,
  with spread, and on the TIMS route in particular?
- **Q3:** two call sites or six?

## Q1: resolved by derivation

**The bound, in one line.** The centre `c_i = (m_i + m_{i+1}) / 2` lies
strictly between its edges, and `f` is strictly monotone, so `f(c_i)` lies
strictly between `u_i` and `u_{i+1}`: strictly within half a step of the
linearised centre `(u_i + u_{i+1}) / 2`. This holds for every bin of every
axis any of the six generators can lay, with no condition on the
configuration. The `constant` axis is the linspace itself and deviates by
float rounding only.

**Per law**, with `t = m_{i+1} / m_i > 1` the edge ratio of one bin, the
deviation as a fraction of the step is:

| law | `f` | deviation | limit `t -> inf` | `1/2 - deviation` |
|---|---|---|---|---|
| constant | `m` | 0 | 0 | 1/2 |
| linear_tof | `sqrt m` | `(sqrt((1 + t^2) / 2) - (1 + t) / 2) / (t - 1)` | `1/sqrt(2) - 1/2 = 0.2071` | `> 0.29` |
| reflector_tof | `ln m` | `ln(cosh h) / (2h)`, `h = ln(t) / 2` | 1/2 | `(h - ln cosh h) / (2h) > 0` |
| orbitrap | `1 / sqrt m` | `((1 + t) / 2 - sqrt(2) t / sqrt(1 + t^2)) / (t - 1)` | 1/2 | `(sqrt(2) t / sqrt(1 + t^2) - 1) / (t - 1) > 0` |
| fticr | `1 / m` | `(t - 1) / (2 (t + 1))` | 1/2 | `1 / (t + 1)` |
| tof | `(2 / sqrt B) asinh(sqrt(B m / A))` | no closed form needed | 1/2 | the general bound |

Each closed form was derived by hand and checked with sympy
(`q1_derivation.py`, appendix A). The FT-ICR entry is the one the brief
sketched, confirmed. The margin `1/2 - d` is strictly positive for every
`t > 1` and vanishes only as one bin spans an unbounded ratio: two bins
over `[1, 10^8]` gives `d = 0.5 - 2 x 10^-8`; the configuration that
produced 0.478 was 1,000 bins over `[1, 10^6]`. On the axes the corpus
resolves to, `d = 6.25 x 10^-7`; over the realistic configurations swept,
the worst per law was `2.8 x 10^-11` (constant), `3.0 x 10^-4`
(linear_tof), `2.2 x 10^-4` (reflector_tof), `8.9 x 10^-4` (orbitrap),
`4.0 x 10^-3` (fticr), `2.9 x 10^-4` (tof, three laws).

**The composed index error.** A value's true nearest bin `j` has the value
in `[c_{j-1}, c_{j+1}]`; `f` is monotone, so its position `x = (f(mz) -
u0) / du` lies between those two centres' positions, each within `d < 1/2`
of its own index; so `x` is in `(j - 3/2, j + 3/2)` and `rint(x)` is one of
`j - 1, j, j + 1`. The repair compares those three in m/z exactly as the
search's caller compares its two, moving left only on a strict
improvement, so a tie resolves to the right as before. A value in the
half-bin skirt beyond either end rounds to `-1` or `n`; the clip sends it
to the edge bin, which is where the search's clip sends it too. The
probe's own rounding does not enter: the argument bounds the interval `x`
lies in, not `x`'s distance from `j`.

**Numerically**, on the real generators: 305 configurations (8 realistic
ranges x 4 bin counts x 8 laws, plus 7 absurd configurations), 541,518,423
probes -- every axis point, every midpoint, both float neighbours of each,
the half-bin skirt, random values, sorted and unsorted -- zero
mismatches against `_nn_map_to_bins`. One defect was found and fixed in
the closed form during this sweep, not in the argument: the first draft
compared with `<=` from the right and took the left of two exactly
equidistant bins, which the `constant` law's exact midpoints expose
(2,787,257 mismatches on the constant law, zero elsewhere). Strict
comparison fixes it; the shipped code and its tests pin the tie rule.

**The guard.** The proof is about the generators; the converter should not
trust that a `MassAxis` came from one. So `_usable_linearisation` measures
the deviation on the axis actually built and requires it under
`NN_LINEARISATION_MARGIN = 0.25` bins on a strictly ascending axis
(`np.searchsorted` maps a duplicated value to its first occurrence; the
repair would not). A quarter is not the bound but half of it, so float
rounding in `x` -- a few ULP of an index up to `10^8` -- has nowhere to
matter. It refuses two bins over six decades and nothing an instrument
asks for. The check is chunked like the bin-width log line beside it, for
the 200-million-bin reason of #251.

## Q2 and Q3: the harness

`q2_harness.py` (appendix B) runs one conversion per subprocess from a
plan, appends each finished run to a JSONL file at once and skips runs
already recorded, so a crash costs one conversion. It converts through
`thyra.convert_msi` with the resampling configuration the CLI builds for
`thyra convert` with no flags (the API's `None` means no resampling), and
patches two things in the process:

- `_build_resampled_mass_axis`, to derive `(f, u0, du)` from the resolved
  plan after the axis is built and measure its deviation;
- `_nn_map_to_bins`, replaced by a dispatcher that computes the closed form
  when the axis is the resampled one and the variant says so, else calls
  the original.

Variants: `current` (the tree as it is), `closed2` (closed form at the two
converter sites only), `closed6` (all six). Modes: `time` (the variant
alone, wall clock and process CPU of the conversion), `verify` (every
closed-form call also computed by the original and asserted equal element
for element), `profile` (every call timed and attributed to its call-site
group). In every mode a running checksum of every mapping (sum and count
of the indices) is kept, and after every run each table matrix in the
store (`data`, `indices`, `indptr` of every `tables/*/X`) is hashed.
Within one (dataset, law) every run of every variant and mode must carry
the same hashes and the same checksum, or the report says MISMATCH.

The source under test was a snapshot of `main` at `25d22e7` copied out of
the worktree before the run, so the implementation could proceed in the
worktree while the measurement ran. The committed implementation was then
run without any patching on five datasets and its store hashes compared
with the recorded ones (section "The committed code" below).

### Environment

Windows 11 Enterprise 10.0.26100, Intel Core i9-14900K, 128 GB RAM, NVMe
system disk holding the datasets and the output. Python 3.13.3, numpy
2.4.6, spatialdata 0.8.0, anndata 0.13.3, zarr 3.3.0, in a venv that
borrows the user site through a `.pth` and pins spatialdata (the recipe in
the brief). Every dataset fits in the page cache (the largest is 27 GB),
so **run 1 of each dataset is the first touch of the session and runs 2
and up are warm**; the machine cannot drop the cache without elevation,
so no run after the first is cold and none is reported as such. Runs
executed one at a time, in the order: run 1 of every dataset and variant,
then the profiles, then runs 2 and 3 everywhere, then runs 4 and 5 on the
nine smaller sets, then the forced laws. **One contamination is known**:
the unit suite (3,171 tests, 12 minutes) ran alongside the driver from
18:27 to 18:39 on 2026-09-15, overlapping the first-touch runs of
`pea_imzml`, `bellini_imzml` and `tims03_biofilm_20um`; their run-1
numbers are inflated and the warm medians are unaffected.

### Datasets

Tier A is `<thyra-checkout>/test_data`, Tier B is the six-set TIMS corpus,
named as in the brief. Every one of them resolved to `reflector_tof` and
`nearest_neighbor` on the auto path -- the whole corpus is one law -- so
the other five laws were forced on one imzML and one TDF set with
`axis_type` in the resampling configuration (the `tof` law with the
timsTOF pair `A = 0.0877, B = 8.74e-4`). `bellini.imzML` was converted
with `spectrum_type="centroid"`, as the prior measurements did. `mock1273`
is the issue comment's mock at Xenium density (150 x 150 pixels, 1,273
peaks per spectrum) through the mock reader, kept as the one sanity point
the brief asked for. Tier C (the lab share) was not needed: no law went
unexercised.

### Results: every dataset, auto law (wall clock, s)

| dataset | axis | bins | mean pts/spec | current r1 (first touch) | current warm median (n) | closed2 warm median (n) | gain | closed6 warm median (n) | gain |
|---|---|---|---|---|---|---|---|---|---|
| bellini_imzml | reflector_tof | 2,373,513 | 2,220 | 54.6 | 18.8 (4) | 9.3 (4) | +50.5% (2.02x) | 9.4 (4) | +50.0% (2.00x) |
| mock1273 | reflector_tof | 599,146 | 1,272 | 33.8 | 17.7 (4) | 13.4 (4) | +24.1% (1.32x) | 13.5 (4) | +23.9% (1.31x) |
| pea_imzml | reflector_tof | 460,495 | 4,722 | 40.0 | 15.1 (4) | 9.7 (4) | +35.8% (1.56x) | 9.7 (4) | +35.6% (1.55x) |
| pea_nedc_d | reflector_tof | 541,610 | 2,321 | 56.6 | 31.7 (4) | 22.5 (4) | +29.1% (1.41x) | 22.3 (4) | +29.8% (1.42x) |
| tims01_msms_315px | reflector_tof | 599,146 | 424 | 50.4 | 27.4 (4) | 26.9 (4) | +2.0% (1.02x) | 26.4 (4) | +3.9% (1.04x) |
| tims02_longramp_1465px | reflector_tof | 21,072 | 6,246 | 31.1 | 27.0 (4) | 26.6 (4) | +1.4% (1.01x) | 26.5 (4) | +2.1% (1.02x) |
| tims03_biofilm_20um | reflector_tof | 92,861 | 6,868 | 110.4 | 80.0 (4) | 72.4 (4) | +9.5% (1.11x) | 69.4 (4) | +13.2% (1.15x) |
| tims04_ratbrain_71k | reflector_tof | 183,258 | 1,433 | 146.6 | 127.6 (4) | 118.2 (4) | +7.3% (1.08x) | 113.7 (4) | +10.9% (1.12x) |
| tims05_shortramp_26k | reflector_tof | 138,629 | 34,243 | 620.8 | 365.8 (4) | 326.8 (4) | +10.7% (1.12x) | 313.7 (4) | +14.2% (1.17x) |
| tims06_glycans_66k | reflector_tof | 240,794 | 66,488 | 2,041.0 | 2,001.8 (2) | 1,715.9 (2) | +14.3% (1.17x) | 1,632.6 (2) | +18.4% (1.23x) |
| xenium_d | reflector_tof | 313,723 | 1,273 | 679.3 | 625.5 (2) | 477.3 (2) | +23.7% (1.31x) | 477.5 (2) | +23.7% (1.31x) |
| xenium_imzml | reflector_tof | 313,713 | 1,273 | 425.4 | 414.3 (2) | 276.3 (2) | +33.3% (1.50x) | 275.9 (2) | +33.4% (1.50x) |

### Spread of the warm runs (min-max, s; n runs)

| dataset | law | current | closed2 | closed6 |
|---|---|---|---|---|
| bellini_imzml | auto | 13.6-19.2 (4) | 9.1-9.6 (4) | 9.3-9.4 (4) |
| mock1273 | auto | 17.6-18.0 (4) | 13.3-13.7 (4) | 13.3-13.6 (4) |
| pea_imzml | auto | 14.2-16.2 (4) | 9.6-9.9 (4) | 9.5-9.8 (4) |
| pea_imzml | constant | 13.0-13.4 (2) | 9.0-9.1 (2) | n/a |
| pea_imzml | fticr | 16.2-16.3 (2) | 10.2-10.3 (2) | n/a |
| pea_imzml | linear_tof | 11.4-11.6 (2) | 9.0-9.0 (2) | n/a |
| pea_imzml | orbitrap | 15.1-15.1 (2) | 10.1-10.1 (2) | n/a |
| pea_imzml | reflector_tof | 13.9-14.1 (2) | 9.7-9.7 (2) | n/a |
| pea_imzml | tof | 12.8-13.0 (2) | 9.9-10.0 (2) | n/a |
| pea_nedc_d | auto | 29.5-32.5 (4) | 22.0-22.8 (4) | 21.2-23.3 (4) |
| tims01_msms_315px | auto | 27.3-28.2 (4) | 26.6-27.0 (4) | 26.1-26.8 (4) |
| tims02_longramp_1465px | auto | 26.7-27.2 (4) | 26.5-26.8 (4) | 26.1-26.8 (4) |
| tims03_biofilm_20um | auto | 79.9-81.4 (4) | 71.4-75.6 (4) | 67.8-69.8 (4) |
| tims03_biofilm_20um | constant | 75.8-79.0 (2) | 69.6-72.2 (2) | n/a |
| tims03_biofilm_20um | fticr | 82.6-84.2 (2) | 74.5-77.1 (2) | n/a |
| tims03_biofilm_20um | linear_tof | 72.0-74.7 (2) | 70.0-70.9 (2) | n/a |
| tims03_biofilm_20um | orbitrap | 80.0-81.6 (2) | 73.3-74.0 (2) | n/a |
| tims03_biofilm_20um | reflector_tof | 79.7-80.3 (2) | 72.9-73.1 (2) | n/a |
| tims03_biofilm_20um | tof | 77.1-77.9 (2) | 75.1-76.2 (2) | n/a |
| tims04_ratbrain_71k | auto | 126.5-128.6 (4) | 117.3-120.1 (4) | 112.8-115.9 (4) |
| tims05_shortramp_26k | auto | 363.5-367.4 (4) | 326.2-328.8 (4) | 308.3-316.7 (4) |
| tims06_glycans_66k | auto | 1963.6-2040.0 (2) | 1706.9-1724.9 (2) | 1616.5-1648.7 (2) |
| xenium_d | auto | 622.5-628.4 (2) | 475.8-478.8 (2) | 477.3-477.7 (2) |
| xenium_imzml | auto | 408.3-420.4 (2) | 275.8-276.7 (2) | 275.6-276.2 (2) |

### Every law forced (warm medians, s)

| dataset | law | bins | deviation (bins) | current median (n) | closed2 median (n) | gain |
|---|---|---|---|---|---|---|
| pea_imzml | constant | 179,997 | 0.00e+00 | 13.1 (3) | 9.0 (3) | +31.2% (1.45x) |
| pea_imzml | fticr | 1,799,780 | 1.25e-06 | 16.2 (3) | 10.3 (3) | +36.0% (1.56x) |
| pea_imzml | linear_tof | 44,059 | 6.13e-06 | 11.4 (3) | 9.0 (3) | +21.6% (1.28x) |
| pea_imzml | orbitrap | 864,841 | 9.38e-07 | 15.1 (3) | 10.1 (3) | +33.3% (1.50x) |
| pea_imzml | reflector_tof | 460,495 | 6.25e-07 | 13.9 (3) | 9.7 (3) | +30.5% (1.44x) |
| pea_imzml | tof | 200,257 | 1.31e-06 | 12.8 (3) | 10.0 (3) | +22.3% (1.29x) |
| tims03_biofilm_20um | constant | 26,000 | 0.00e+00 | 76.0 (3) | 70.9 (3) | +6.6% (1.07x) |
| tims03_biofilm_20um | fticr | 337,662 | 4.37e-07 | 83.1 (3) | 75.1 (3) | +9.5% (1.11x) |
| tims03_biofilm_20um | linear_tof | 7,897 | 4.14e-06 | 72.9 (3) | 70.9 (3) | +2.7% (1.03x) |
| tims03_biofilm_20um | orbitrap | 176,679 | 5.55e-07 | 81.6 (3) | 73.7 (3) | +9.7% (1.11x) |
| tims03_biofilm_20um | reflector_tof | 92,861 | 6.25e-07 | 79.7 (3) | 73.1 (3) | +8.3% (1.09x) |
| tims03_biofilm_20um | tof | 40,348 | 1.25e-06 | 77.9 (3) | 75.1 (3) | +3.6% (1.04x) |

### Profile: mapping time by call-site group (current tree, warm, one run each)

| dataset | wall s | converter sites s (share) | sibling sites s (share) | converter calls / probes | sibling calls / probes |
|---|---|---|---|---|---|
| tims01_msms_315px | 28.6 | 0.1 (0.3%) | 0.4 (1.3%) | 1,426 / 604,212 | 22,103 / 907,027 |
| tims02_longramp_1465px | 27.5 | 0.6 (2.0%) | 0.3 (0.9%) | 2,930 / 18,300,488 | 1,465 / 9,150,857 |
| tims03_biofilm_20um | 78.8 | 10.2 (12.9%) | 5.0 (6.3%) | 37,676 / 258,747,608 | 18,838 / 129,377,751 |
| pea_nedc_d | 32.1 | 10.2 (31.7%) | 0.0 (0.0%) | 67,380 / 156,404,432 | 0 / 0 |
| tims04_ratbrain_71k | 127.4 | 14.6 (11.4%) | 6.5 (5.1%) | 142,498 / 204,219,086 | 71,249 / 102,110,835 |
| tims05_shortramp_26k | 365.1 | 55.1 (15.1%) | 26.6 (7.3%) | 52,174 / 1,786,581,580 | 26,087 / 893,296,892 |
| xenium_d | 641.9 | 189.7 (29.6%) | 0.0 (0.0%) | 1,837,710 / 2,338,675,584 | 0 / 0 |
| tims06_glycans_66k | 1,947.1 | 304.5 (15.6%) | 144.5 (7.4%) | 131,712 / 8,757,254,600 | 65,856 / 4,378,632,250 |
| xenium_imzml | 422.1 | 168.7 (40.0%) | 0.0 (0.0%) | 1,837,710 / 2,338,675,584 | 0 / 0 |
| pea_imzml | 16.0 | 6.2 (38.7%) | 0.0 (0.0%) | 25,474 / 120,299,250 | 0 / 0 |

### Reading the tables

**Identity.** 268 of 268 runs succeeded. Within every (dataset, law) group
-- 24 groups, 16 runs each on the auto law, 7 on each forced law -- every
variant and mode produced the same store hashes for every table (summed,
MS/MS where present) and the same mapping checksum. The verify runs
compared 4,388,000 closed-form calls against the search element for
element: `bellini` 32,768, `mock1273` 45,000, `pea_imzml` 25,474,
`pea_nedc_d` 67,380, `tims01` 23,529, `tims02` 4,395, `tims03` 56,514,
`tims04` 213,747, `tims05` 78,261, `tims06` 197,568, `xenium_d`
1,837,710, `xenium_imzml` 1,837,710, plus the forced laws. No call fell
back to the search on any resampled axis (`fallback_calls = 0`
throughout), which also shows the linearisation reached every call site
through the axis identity the harness keyed on.

**The sanity point.** The brief asked for the xenium-density mock to be
reproduced: 17.7 s to 13.4 s warm, +24.1 percent, against the +27.3
percent reported on 2026-09-14 (19.27 to 14.01 s, a different session and
a first-touch pair). Same direction, same size.

**Where the win is.** The gain scales with mean peaks per spectrum times
log of the bin count, which is the search it removes, divided by the
conversion's other work:

- imzML and TDF-without-mobility sets, 1,300 to 4,700 peaks per spectrum
  on 300k to 2.4M bins: **+24 to +50 percent** of wall clock
  (`xenium_imzml` +33, `xenium_d` +24, `pea` +36, `pea_nedc` +29,
  `bellini` +50). The mapping was 30 to 40 percent of these conversions
  (profile table) and is now 4 to 5 times cheaper.
- TIMS sets with the heatmap fed from the passes, 1,400 to 66,500 points
  per frame on 93k to 241k bins: **+7 to +14 percent** at the two
  converter sites, **+11 to +18 percent** at all six. The mapping is a
  smaller share here (11 to 16 percent at the converter sites, 5 to 7 at
  the sibling sites) because the frames are read from a Bruker SDK and
  scattered through a 10 GB memmap, and both cost more than the search.
- the two tiny sets, 315 and 1,465 pixels: +1 to +4 percent, within the
  spread, as expected on a 27-second conversion whose mapping is under a
  second.

**Every law.** Forced on `pea_imzml` and `tims03`, all six laws convert
bit-identically through the closed form and gain in proportion to their
bin count: `linear_tof` and `tof` lay the fewest bins here (44k and 200k
on pea; 7.9k and 40k on the biofilm) and gain the least (+22 percent and
+3 to +4), `fticr` and `orbitrap` the most (1.8M and 865k on pea, +36 and
+33). The measured deviation of every forced axis is under `10^-5` bins.

**Spread.** Warm runs of one variant on one dataset differ by under 2
percent on every set except `bellini` (13.6 to 19.2 s on `current`, four
runs; the imzML parser's first pass dominates this small file and varies)
and the two giants where only two warm runs exist. The gains above are
five to fifty times the spread.

### Q3: two sites, and what the sibling sites are worth

On the TIMS route the sibling sites map each frame's unique m/z once in
pass 1 (for the heatmap), while the summed table maps the same array twice
(once per pass): the sibling probes are exactly half the converter probes
on every TIMS set in the profile table. Their share of a whole-slide wall
clock is **5.1 to 7.4 percent**; the closed form recovers about four
fifths of that (`closed6` minus `closed2`: +3.6 points on `tims04`, +3.5
on `tims05`, +4.1 on `tims06`, +3.7 on `tims03`). The issue body's "~110 s
of a 171 s pass" is not what this slide does now: on `tims05` the whole
mapping, both groups, is 82 s of a 365 s warm conversion.

So the prior assessment holds in substance: the two converter sites carry
the imzML and mobility-free TDF win whole, and two thirds of the TIMS win.
The sibling sites are worth 3 to 4 points on a TIMS whole slide, not
nothing, and the plumbing is smaller than the "five container classes"
estimate -- `SiblingPasses` and `MsmsAccumulator` each take the axis in
their constructor and could take the linearisation beside it, and the
three mapping helpers would take one optional argument -- but it is a
separate change. The better follow-up is different: the summed spectrum
of a TDF frame *is* the frame's `unique_mz` (`spectrum()` returns
`self.unique_mz` under `scan_sum`), so the sibling sinks re-map an array
the summed table has just mapped. Sharing that mapping per frame would
remove the sibling search entirely rather than make it four times
cheaper. Neither is in this change; both are recorded on the issue.

### The committed code

The harness measured a monkeypatched snapshot of `main` (`25d22e7`). The
implementation that shipped (`perf/closed-form-bin-index`, first commit
`f5c095d`) was then run without any patching, through `thyra.convert_msi`
with the CLI's default resampling configuration, on five datasets, and its
store hashes compared with the ones recorded above for the same datasets:

| dataset | wall (s) | log says closed form | store hash identical to the recorded one |
|---|---|---|---|
| `pea_imzml` | 10.2 | yes | yes |
| `bellini_imzml` | 9.5 | yes | yes |
| `pea_nedc_d` | 21.5 | yes | yes |
| `tims03_biofilm_20um` | 71.2 | yes | yes |
| `tims01_msms_315px` | 26.9 | yes | yes |

The wall clocks sit on the `closed2` medians (9.7, 9.3, 22.5, 72.4,
26.9 s). The shipped code adds `forward()` to every generator, an
`AxisLinearisation` on the `MassAxis` it returns, the closed-form route
and its guard in `base_spatialdata_converter.py`, and 192 tests
(`tests/unit/resampling/test_axis_linearisation.py`,
`tests/unit/converters/test_nn_closed_form_bins.py`); the full unit suite
passes (3,171 tests). Micro-benchmarked, the shipped `_nn_map_to_bins`
with a linearisation is 4.2x the search at 1,273 sorted peaks on 600k
bins, 4.9x at 4,722, 3.5x at 15,000 on 1.05M bins, and 1.8x at 400.

### Things learned, for whoever comes next

- The API's `resampling_config=None` disables resampling; the CLI's
  default is a dict with `method="auto"`, `axis_type="auto"` and every
  other key `None`. A harness that passes `None` measures
  `_map_mass_to_indices`, a different function.
- `pea_nedc_d` (`20231109_PEA_NEDC.d`) has no mobility dimension: zero
  sibling calls. It is a TDF acquisition with TIMS off, which is why it
  behaves like an imzML set in every table.
- A tie between two candidate bins must be resolved with strict
  comparisons from the right; `<=` from the right takes the left bin, and
  only the `constant` law's exact midpoints show it.
- Warm-run spread on a 27-second conversion is 1 to 2 percent; a single
  first-touch pair says nothing at that scale. `pea_nedc_d` run 1 gave
  `closed6` 31.0 s against `closed2` 40.0 s on identical code paths (it
  has no sibling sites).

## Appendix A: `q1_derivation.py`

```python
"""Q1 of thyra#295: is the +-1 repair provably sufficient?

Three parts.

1. The general lemma, stated and checked symbolically per law:
   every non-uniform generator lays edges u_i = u0 + i*du in u = f(m/z)
   and takes bin centres c_i = (m_i + m_{i+1}) / 2 in m/z. Because
   m_i < c_i < m_{i+1} and f is strictly monotone, f(c_i) lies strictly
   between u_i and u_{i+1}, so the deviation from the linearised centre
   (u_i + u_{i+1}) / 2 is strictly less than |du| / 2 for EVERY bin of
   EVERY axis any of these generators can build. Per law the deviation has
   a closed form in the edge ratio; sympy verifies each is < 1/2.

2. The composed index error. A probe's true nearest bin j* has the probe
   in [c_{j*-1}, c_{j*+1}], so its linearised position x lies in
   [j*-1-delta, j*+1+delta]; with delta < 1/2, rint(x) is within one of
   j*, and the +-1 repair (which compares in m/z space with the same
   comparisons searchsorted's caller makes) returns j* exactly.

3. Numerics: the deviation measured on real generator output for a sweep
   of configurations including the absurd ones, and an exhaustive probe
   test of the closed form against the reference on every axis point,
   every midpoint, its float neighbours, the half-bin skirt and random
   values -- bit-identical or the script fails.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import sympy as sp

sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else ".")

from thyra.converters.spatialdata.base_spatialdata_converter import (  # noqa: E402
    _nn_map_to_bins,
)
from thyra.resampling.mass_axis.fticr_generator import FTICRAxisGenerator  # noqa: E402
from thyra.resampling.mass_axis.linear_generator import LinearAxisGenerator  # noqa: E402
from thyra.resampling.mass_axis.linear_tof_generator import (  # noqa: E402
    LinearTOFAxisGenerator,
)
from thyra.resampling.mass_axis.orbitrap_generator import (  # noqa: E402
    OrbitrapAxisGenerator,
)
from thyra.resampling.mass_axis.reflector_tof_generator import (  # noqa: E402
    ReflectorTOFAxisGenerator,
)
from thyra.resampling.mass_axis.tof_generator import (  # noqa: E402
    MRT_TOF_LAW,
    PHI_NANOTOF_LAW,
    TIMSTOF_TOF_LAW,
    TOFAxisGenerator,
)

# --------------------------------------------------------------------------
# Part 1: symbolic. d(t) is the deviation / |du| for one bin whose edges
# have ratio t = m_{i+1} / m_i > 1 (or the equivalent parameter). Show
# 1/2 - d(t) > 0 for all t > 1.
# --------------------------------------------------------------------------

t = sp.symbols("t", positive=True)  # t > 1 assumed in the argument below
h = sp.symbols("h", positive=True)


def sym_fticr():
    # u = 1/m; edges a = 1/m_i, b = 1/m_{i+1}; t = m_{i+1}/m_i = a/b
    m1 = sp.Integer(1)
    m2 = t
    c = (m1 + m2) / 2
    f = lambda m: 1 / m
    ubar = (f(m1) + f(m2)) / 2
    d = sp.simplify((ubar - f(c)) / (f(m1) - f(m2)))
    margin = sp.simplify(sp.Rational(1, 2) - d)
    return d, margin


def sym_orbitrap():
    m1 = sp.Integer(1)
    m2 = t
    c = (m1 + m2) / 2
    f = lambda m: 1 / sp.sqrt(m)
    ubar = (f(m1) + f(m2)) / 2
    d = sp.simplify((ubar - f(c)) / (f(m1) - f(m2)))
    margin = sp.simplify(sp.Rational(1, 2) - d)
    return d, margin


def sym_linear_tof():
    m1 = sp.Integer(1)
    m2 = t
    c = (m1 + m2) / 2
    f = lambda m: sp.sqrt(m)
    ubar = (f(m1) + f(m2)) / 2
    d = sp.simplify((f(c) - ubar) / (f(m2) - f(m1)))
    margin = sp.simplify(sp.Rational(1, 2) - d)
    return d, margin


def sym_reflector_tof():
    # u = ln m; edges a, b with b - a = 2h; centre c = (e^a + e^b)/2
    a = sp.symbols("a", real=True)
    b = a + 2 * h
    f = sp.log
    c = (sp.exp(a) + sp.exp(b)) / 2
    ubar = (a + b) / 2
    d = sp.simplify((f(c) - ubar) / (b - a))
    margin = sp.simplify(sp.Rational(1, 2) - d)
    return d, margin


def part1():
    print("== Part 1: symbolic deviation per law, 1/2 - d(t) ==")
    for name, fn in [
        ("fticr", sym_fticr),
        ("orbitrap", sym_orbitrap),
        ("linear_tof", sym_linear_tof),
        ("reflector_tof", sym_reflector_tof),
    ]:
        d, margin = fn()
        print(f"  {name:14s} d = {d}")
        print(f"  {'':14s} 1/2 - d = {margin}")
        # Evaluate the margin on a log-spaced grid of the ratio to show it is
        # positive everywhere and how it tends to zero.
        sym = t if name != "reflector_tof" else h
        grid = [1 + 10.0**k for k in range(-9, 7)] if sym is t else [10.0**k for k in range(-6, 4)]
        vals = [float(margin.subs(sym, g)) for g in grid]
        assert all(v > 0 for v in vals), (name, vals)
        lim = sp.limit(margin, sym, sp.oo)
        lim0 = sp.limit(d, sym, 1 if sym is t else 0)
        print(f"  {'':14s} margin > 0 on the grid; limit of margin at infinity = {lim}; "
              f"d at a vanishing bin = {lim0}")
    print("  tof (two-term): f(m) = (2/sqrt(B)) asinh(sqrt(B m / A)) is strictly increasing;")
    print("     the general lemma (m_i < c_i < m_{i+1}  =>  u_i < f(c_i) < u_{i+1}) applies verbatim.")
    print("  constant: the axis IS the linspace, deviation 0 up to the last ULP.")
    print()


# --------------------------------------------------------------------------
# Part 3: numerics on the real generators.
# --------------------------------------------------------------------------


def law_forward(name, tof_law=None):
    if name == "constant":
        return lambda m: np.asarray(m, dtype=np.float64)
    if name == "linear_tof":
        return lambda m: np.sqrt(np.asarray(m, dtype=np.float64))
    if name == "reflector_tof":
        return lambda m: np.log(np.asarray(m, dtype=np.float64))
    if name == "orbitrap":
        return lambda m: 1.0 / np.sqrt(np.asarray(m, dtype=np.float64))
    if name == "fticr":
        return lambda m: 1.0 / np.asarray(m, dtype=np.float64)
    if name == "tof":
        return TOFAxisGenerator(*tof_law).cumulative
    raise ValueError(name)


def generator_for(name, tof_law=None):
    return {
        "constant": LinearAxisGenerator,
        "linear_tof": LinearTOFAxisGenerator,
        "reflector_tof": ReflectorTOFAxisGenerator,
        "orbitrap": OrbitrapAxisGenerator,
        "fticr": FTICRAxisGenerator,
    }[name]() if name != "tof" else TOFAxisGenerator(*tof_law)


def linearisation(name, min_mz, max_mz, n, tof_law=None):
    """(f, f0, du) such that f(centre_i) ~= f0 + i*du."""
    f = law_forward(name, tof_law)
    if name == "constant":
        f0 = float(min_mz)
        du = (float(max_mz) - float(min_mz)) / (n - 1)
        return f, f0, du
    u_lo = float(f(min_mz))
    u_hi = float(f(max_mz))
    du = (u_hi - u_lo) / n
    f0 = u_lo + 0.5 * du
    return f, f0, du


def deviation(name, axis, f, f0, du):
    i = np.arange(axis.size, dtype=np.float64)
    return float(np.max(np.abs(f(axis) - (f0 + i * du)) / abs(du)))


def closed_form(axis, mzs, f, f0, du):
    """The issue's proposal: rint of the linearised position, +-1 repair, ties right."""
    n = axis.size
    k = np.rint((f(mzs) - f0) / du).astype(np.int64)
    np.clip(k, 0, n - 1, out=k)
    lo = np.maximum(k - 1, 0)
    hi = np.minimum(k + 1, n - 1)
    d_lo = np.abs(axis[lo] - mzs)
    d_k = np.abs(axis[k] - mzs)
    d_hi = np.abs(axis[hi] - mzs)
    # Ties resolve to the right: the reference takes the left neighbour
    # only when strictly closer, so scan from the right and move left
    # only on a strict improvement.
    best = hi
    best_d = d_hi
    take = d_k < best_d
    best = np.where(take, k, best)
    best_d = np.where(take, d_k, best_d)
    take = d_lo < best_d
    best = np.where(take, lo, best)
    return best


def probes_for(axis, min_mz, max_mz, rng, n_random=200_000):
    mids = (axis[:-1] + axis[1:]) / 2
    parts = [
        axis,
        mids,
        np.nextafter(mids, np.inf),
        np.nextafter(mids, -np.inf),
        np.nextafter(axis, np.inf),
        np.nextafter(axis, -np.inf),
        rng.uniform(min_mz, max_mz, n_random),
        # the half-bin skirt beyond the first and last centre
        rng.uniform(min_mz, axis[0], 2000),
        rng.uniform(axis[-1], max_mz, 2000),
        np.array([min_mz, max_mz, axis[0], axis[-1]]),
    ]
    return np.sort(np.concatenate(parts))


def part3():
    print("== Part 3: numerics on the real generators ==")
    rng = np.random.default_rng(295)
    # (name, min, max, n, tof_law) -- realistic first, then the absurd.
    configs = []
    realistic = [
        (50.0, 1000.0), (100.0, 2000.0), (200.0, 500.0), (300.0, 3000.0),
        (400.0, 800.0), (900.0, 3000.0), (12.0, 400.0), (540.0, 600.0),
    ]
    for lo, hi in realistic:
        for n in (2000, 50_000, 190_000, 1_050_000):
            for name in ("constant", "linear_tof", "reflector_tof", "orbitrap", "fticr"):
                configs.append((name, lo, hi, n, None))
            for law in (MRT_TOF_LAW, TIMSTOF_TOF_LAW, PHI_NANOTOF_LAW):
                configs.append(("tof", lo, hi, n, law))
    absurd = [
        (1.0, 1e6, 1000), (0.001, 1e6, 1000), (1.0, 1e6, 2), (1.0, 1e6, 3),
        (1.0, 1e8, 2), (10.0, 1e5, 10), (1.0, 100.0, 2),
    ]
    for lo, hi, n in absurd:
        for name in ("constant", "linear_tof", "reflector_tof", "orbitrap", "fticr"):
            configs.append((name, lo, hi, n, None))
        for law in (TIMSTOF_TOF_LAW, PHI_NANOTOF_LAW):
            configs.append(("tof", lo, hi, n, law))

    worst = {}
    n_probes = 0
    n_mismatch = 0
    t0 = time.perf_counter()
    for name, lo, hi, n, law in configs:
        gen = generator_for(name, law)
        axis = gen.generate_axis(lo, hi, n).mz_values.astype(np.float64)
        f, f0, du = linearisation(name, lo, hi, n, law)
        dev = deviation(name, axis, f, f0, du)
        key = name if law is None else f"tof{law}"
        prev = worst.get(key, (0.0, None))
        if dev > prev[0]:
            worst[key] = (dev, (lo, hi, n))
        assert dev < 0.5, (name, lo, hi, n, dev)
        if not bool(np.all(np.diff(axis) > 0)):
            print(f"  SKIP probes: axis not strictly ascending for {name} {lo} {hi} {n}")
            continue
        probes = probes_for(axis, lo, hi, rng, n_random=min(200_000, 40 * n))
        ref = _nn_map_to_bins(axis, probes)
        got = closed_form(axis, probes, f, f0, du)
        bad = int(np.count_nonzero(ref != got))
        n_probes += probes.size
        n_mismatch += bad
        if bad:
            print(f"  MISMATCH {name} {lo} {hi} {n} law={law}: {bad} of {probes.size}, dev={dev:.4f}")
    dt = time.perf_counter() - t0
    print(f"  {len(configs)} configurations, {n_probes:,} probes, {n_mismatch} mismatches, {dt:.1f} s")
    print("  worst deviation / |du| per law (0.5 is the cliff):")
    for key, (dev, cfg) in sorted(worst.items()):
        print(f"    {key:40s} {dev:.6f}  at min={cfg[0]:g} max={cfg[1]:g} n={cfg[2]}")
    print()
    # Realistic-only worst, which is what a build-time guard would see.
    print("  worst deviation on the realistic configurations only:")
    worst_r = {}
    for name, lo, hi, n, law in configs:
        if (lo, hi) not in realistic:
            continue
        gen = generator_for(name, law)
        axis = gen.generate_axis(lo, hi, n).mz_values.astype(np.float64)
        f, f0, du = linearisation(name, lo, hi, n, law)
        dev = deviation(name, axis, f, f0, du)
        key = name if law is None else f"tof{law}"
        worst_r[key] = max(worst_r.get(key, 0.0), dev)
    for key, dev in sorted(worst_r.items()):
        print(f"    {key:40s} {dev:.3e}")
    assert n_mismatch == 0


if __name__ == "__main__":
    part1()
    part3()
    print("Q1: every law's deviation is provably < 1/2; closed form bit-identical on every probe.")
```

## Appendix B: `q2_harness.py`

Run as `THYRA_WORKTREE=<source> THYRA_TEST_DATA=<test_data> THYRA_TIMS_DATA=<tims corpus> python q2_harness.py --plan plan.json`; the plan is a JSON list of `{dataset, law, variant, mode, run[, extra]}`.

```python
"""Q2/Q3 of thyra#295: does the closed-form bin index pay on real data?

One conversion per subprocess, driven by a plan. Every completed run is
appended to a JSONL file the moment it finishes, and a run already in the
file is skipped, so a crash costs one conversion, not the night.

Variants (what ``_nn_map_to_bins`` does when the axis is the resampled one):

    current   the tree as it is: np.searchsorted + left/right comparison
    closed2   closed form at the two converter call sites only
              (_nearest_neighbor_resample, _build_nn_shared_cache); the
              four sibling-table sites keep searchsorted
    closed6   closed form at all six call sites

Modes:

    time      the variant alone, wall clock and CPU time of convert_msi();
              a cheap running checksum of every mapping (sum of the
              indices and their count) is kept in both variants
    verify    every closed-form call is also computed by the reference and
              the two are asserted equal, element for element
    profile   every call is timed and attributed to its call-site group,
              converter or sibling, for the Q3 split

The output store's every table matrix is hashed after each run, so a
variant that changed a stored value is caught even where the per-call
checksum would not see it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKTREE = Path(os.environ.get("THYRA_WORKTREE", HERE)).resolve()
TEST_DATA = Path(os.environ["THYRA_TEST_DATA"])
TIMS = Path(os.environ["THYRA_TIMS_DATA"])
RESULTS = HERE / "q2_results.jsonl"
OUT_ROOT = HERE / "out"
LOGS = HERE / "logs"

#: name -> (input path, reader options)
DATASETS = {
    "xenium_imzml": (TEST_DATA / "20240826_xenium_0041899.imzML", {}),
    "bellini_imzml": (TEST_DATA / "bellini.imzML", {"spectrum_type": "centroid"}),
    "pea_imzml": (TEST_DATA / "pea.imzML", {}),
    "pea_nedc_d": (TEST_DATA / "20231109_PEA_NEDC.d", {}),
    "xenium_d": (TEST_DATA / "20240826_Xenium_0041899.d", {}),
    "tims01_msms_315px": (TIMS / "01_tiny_msms_315px" / "220425_MSMS_pos_brain1.d", {}),
    "mock1273": (None, {}),  # the issue comment's mock at xenium density, via the mock reader
    "tims02_longramp_1465px": (None, {}),  # filled from a glob below
    "tims03_biofilm_20um": (TIMS / "03_biofilm_maldi_vs_maldi2_20um" / "20241023_MALDI2_TIMS_20um.d", {}),
    "tims04_ratbrain_71k": (TIMS / "04_ratbrain_71k_px_9990_scans" / "20200915_Ratbrain_B2S3_R.d", {}),
    "tims05_shortramp_26k": (TIMS / "05_maldi2_shortramp_26k_px" / "250213_SlideV_Retest.d", {}),
    "tims06_glycans_66k": (None, {}),
}


def _fill_globs() -> None:
    for name, sub in (("tims02_longramp_1465px", "02_tiny_longramp_1465px"),
                      ("tims06_glycans_66k", "06_glycans_highmass_66k_px")):
        if DATASETS[name][0] is None:
            hits = sorted((TIMS / sub).glob("*.d"))
            if hits:
                DATASETS[name] = (hits[0], DATASETS[name][1])


_fill_globs()

# --------------------------------------------------------------------------
# The variant machinery (only imported inside the child process)
# --------------------------------------------------------------------------

STATE: dict = {"law": None, "delta": None, "plan": None}
COUNT: dict = {
    "conv_calls": 0, "conv_probes": 0, "sib_calls": 0, "sib_probes": 0,
    "fallback_calls": 0, "closed_calls": 0, "checksum": 0, "n_mapped": 0,
    "t_conv": 0.0, "t_sib": 0.0, "verified_calls": 0,
}
CONVERTER_SITES = {"_nearest_neighbor_resample", "_build_nn_shared_cache"}


def linearisation(plan, axis):
    """(f, f0, du, delta): f(axis[i]) ~= f0 + i*du, delta the worst deviation / |du|."""
    import numpy as np
    from thyra.resampling.mass_axis.tof_generator import TOFAxisGenerator

    name = plan["axis_type"].value
    lo, hi, n = float(plan["min_mz"]), float(plan["max_mz"]), int(plan["target_bins"])
    if name == "constant":
        f = lambda m: m  # noqa: E731
        f0 = lo
        du = (hi - lo) / (n - 1)
    else:
        f = {
            "linear_tof": np.sqrt,
            "reflector_tof": np.log,
            "orbitrap": lambda m: 1.0 / np.sqrt(m),
            "fticr": lambda m: 1.0 / m,
            "tof": (TOFAxisGenerator(plan["tof_a"], plan["tof_b"]).cumulative if name == "tof" else None),
        }[name]
        u_lo, u_hi = float(f(lo)), float(f(hi))
        du = (u_hi - u_lo) / n
        f0 = u_lo + 0.5 * du
    i = np.arange(axis.size, dtype=np.float64)
    delta = float(np.max(np.abs(f(axis) - (f0 + i * du)) / abs(du)))
    return f, f0, du, delta


def closed_form(axis, mzs, f, f0, du):
    import numpy as np

    n = axis.size
    k = np.rint((f(mzs) - f0) / du).astype(np.int64)
    np.clip(k, 0, n - 1, out=k)
    lo = np.maximum(k - 1, 0)
    hi = np.minimum(k + 1, n - 1)
    d_lo = np.abs(axis[lo] - mzs)
    d_k = np.abs(axis[k] - mzs)
    d_hi = np.abs(axis[hi] - mzs)
    best = hi
    best_d = d_hi
    take = d_k < best_d
    best = np.where(take, k, best)
    best_d = np.where(take, d_k, best_d)
    take = d_lo < best_d
    return np.where(take, lo, best)


def install(variant: str, mode: str) -> None:
    import numpy as np
    from thyra.converters.spatialdata import base_spatialdata_converter as base
    from thyra.converters.spatialdata import mobility_heatmap as mh
    from thyra.converters.spatialdata.base_spatialdata_converter import (
        BaseSpatialDataConverter,
    )

    orig_map = base._nn_map_to_bins
    orig_build = BaseSpatialDataConverter._build_resampled_mass_axis

    def build(self):
        orig_build(self)
        plan = self._resolved_resampling_plan
        axis = self._common_mass_axis
        f, f0, du, delta = linearisation(plan, axis)
        STATE["delta"] = delta
        STATE["plan"] = {
            "axis_type": plan["axis_type"].value,
            "target_bins": int(plan["target_bins"]),
            "min_mz": float(plan["min_mz"]),
            "max_mz": float(plan["max_mz"]),
            "mass_width_da": float(plan["mass_width_da"]),
            "reference_mz": float(plan["reference_mz"]),
            "method": str(getattr(plan["method"], "value", plan["method"])),
            "tof_law": [plan["tof_a"], plan["tof_b"]] if "tof_a" in plan else None,
            "strictly_ascending": bool(np.all(np.diff(axis) > 0)),
        }
        # The build-time guard: fall back to searchsorted when the axis
        # deviates from its linearisation by more than a quarter bin (the
        # proof needs < 1/2; a quarter leaves the float error nowhere to go)
        # or is not strictly ascending.
        if delta < 0.25 and STATE["plan"]["strictly_ascending"] and variant != "current":
            STATE["law"] = (axis, f, f0, du)
        else:
            STATE["law"] = None

    BaseSpatialDataConverter._build_resampled_mass_axis = build

    def make(group):
        def mapped(axis, mzs):
            law = STATE["law"]
            g = group
            if g == "auto":
                g = "conv" if sys._getframe(1).f_code.co_name in CONVERTER_SITES else "sib"
            use_closed = (
                law is not None
                and axis is law[0]
                and (variant == "closed6" or (variant == "closed2" and g == "conv"))
            )
            if mode == "profile":
                t0 = time.perf_counter()
            if use_closed:
                idx = closed_form(axis, mzs, law[1], law[2], law[3])
                COUNT["closed_calls"] += 1
                if mode == "verify":
                    ref = orig_map(axis, mzs)
                    if not np.array_equal(ref, idx):
                        bad = int(np.count_nonzero(ref != idx))
                        raise AssertionError(f"closed form differs from reference on {bad} of {idx.size}")
                    COUNT["verified_calls"] += 1
            else:
                idx = orig_map(axis, mzs)
                if law is not None and axis is not law[0]:
                    COUNT["fallback_calls"] += 1
            if mode == "profile":
                COUNT["t_" + g] += time.perf_counter() - t0
            COUNT[g + "_calls"] += 1
            COUNT[g + "_probes"] += int(mzs.size)
            COUNT["checksum"] += int(idx.sum())
            COUNT["n_mapped"] += int(idx.size)
            return idx

        return mapped

    base._nn_map_to_bins = make("auto")
    mh._nn_map_to_bins = make("sib")


# --------------------------------------------------------------------------
# Store hashing
# --------------------------------------------------------------------------


def store_hash(path: Path) -> dict:
    import numpy as np
    import zarr

    root = zarr.open_group(str(path), mode="r")
    digests = {}
    tables = root["tables"]
    for name, table in sorted(tables.groups()):
        h = hashlib.blake2b(digest_size=16)
        x = table["X"]
        keys = sorted(k for k, _ in x.arrays()) if isinstance(x, zarr.Group) else None
        if keys is None:
            arr = x
            step = 1 << 22
            for s in range(0, arr.shape[0], step):
                h.update(np.ascontiguousarray(arr[s : s + step]).tobytes())
        else:
            for k in keys:
                arr = x[k]
                h.update(k.encode())
                step = 1 << 24
                for s in range(0, arr.shape[0], step):
                    h.update(np.ascontiguousarray(arr[s : s + step]).tobytes())
        digests[name] = h.hexdigest()
    return digests


# --------------------------------------------------------------------------
# One run
# --------------------------------------------------------------------------


def run_one(dataset: str, law: str, variant: str, mode: str, run: int, extra: dict) -> dict:
    sys.path.insert(0, str(WORKTREE))
    LOGS.mkdir(exist_ok=True)
    OUT_ROOT.mkdir(exist_ok=True)
    tag = f"{dataset}__{law}__{variant}__{mode}__r{run}"
    logging.basicConfig(
        filename=str(LOGS / f"{tag}.log"), filemode="w", level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    import thyra  # noqa: F401  (registers readers/converters)

    install(variant, mode)
    src, reader_options = DATASETS[dataset]
    reader_options = dict(reader_options)
    reader_options.update(extra.get("reader_options", {}))
    out = OUT_ROOT / f"{tag}.zarr"
    if out.exists():
        shutil.rmtree(out)
    # Exactly what the CLI builds for `thyra convert` with no resampling
    # flags (``_build_resampling_config`` in thyra/__main__.py); the API's
    # None means no resampling at all.
    rc = {
        "method": "auto", "axis_type": law, "target_bins": None, "min_mz": None,
        "max_mz": None, "width_at_mz": None, "reference_mz": 1000.0,
        "gap_tolerance_da": None, "tof_a": None, "tof_b": None,
    }
    if extra.get("resampling_config"):
        rc.update(extra["resampling_config"])
    kwargs = dict(extra.get("kwargs", {}))
    if dataset == "mock1273":
        sys.path.insert(0, str(WORKTREE / "tests" / "fixtures"))
        from mock_msi_generator import MockMSIConfig, MockMSIReader
        from thyra.converters.spatialdata.streaming_converter import (
            StreamingSpatialDataConverter,
        )

        cfg = MockMSIConfig(n_x=150, n_y=150, n_mz_bins=50000, peaks_per_spectrum=(1273, 1274))
        reader = MockMSIReader(cfg)
        t_wall = time.perf_counter()
        t_cpu = time.process_time()
        converter = StreamingSpatialDataConverter(
            reader=reader, output_path=out, dataset_id="mock", pixel_size_um=cfg.pixel_size_um,
            resampling_config=rc,
        )
        ok = converter.convert()
        wall = time.perf_counter() - t_wall
        cpu = time.process_time() - t_cpu
        reader.close()
    else:
        t_wall = time.perf_counter()
        t_cpu = time.process_time()
        ok = thyra.convert_msi(
            str(src), str(out), resampling_config=rc, reader_options=reader_options or None,
            **kwargs,
        )
        wall = time.perf_counter() - t_wall
        cpu = time.process_time() - t_cpu
    rec = {
        "dataset": dataset, "law": law, "variant": variant, "mode": mode, "run": run,
        "ok": bool(ok), "wall_s": round(wall, 3), "cpu_s": round(cpu, 3),
        "plan": STATE["plan"], "delta": STATE["delta"],
        "counts": {k: (round(v, 4) if isinstance(v, float) else v) for k, v in COUNT.items()},
        "store": store_hash(out) if ok else None,
        "extra": extra,
        "finished": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "numpy": __import__("numpy").__version__,
        "thyra_file": thyra.__file__,
    }
    if ok and not extra.get("keep"):
        shutil.rmtree(out, ignore_errors=True)
    return rec


# --------------------------------------------------------------------------
# The driver
# --------------------------------------------------------------------------


def key_of(r: dict) -> tuple:
    return (r["dataset"], r["law"], r["variant"], r["mode"], r["run"], json.dumps(r.get("extra", {}), sort_keys=True))


def load_results() -> list:
    if not RESULTS.exists():
        return []
    return [json.loads(line) for line in RESULTS.read_text().splitlines() if line.strip()]


def drive(plan_items: list, python: str) -> None:
    done = {key_of(r) for r in load_results()}
    for item in plan_items:
        item.setdefault("extra", {})
        if key_of(item) in done:
            continue
        args = [python, str(Path(__file__)), "--one", json.dumps(item)]
        print(f"[{datetime.now():%H:%M:%S}] {item['dataset']} {item['law']} {item['variant']} "
              f"{item['mode']} r{item['run']}", flush=True)
        t0 = time.perf_counter()
        proc = subprocess.run(args, capture_output=True, text=True, env={**os.environ, "THYRA_WORKTREE": str(WORKTREE)})
        if proc.returncode != 0:
            print(f"   FAILED after {time.perf_counter() - t0:.0f} s\n{proc.stderr[-3000:]}", flush=True)
            rec = {**item, "ok": False, "error": proc.stderr[-3000:], "wall_s": None,
                   "finished": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        else:
            rec = json.loads(proc.stdout.strip().splitlines()[-1])
        with RESULTS.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
        done.add(key_of(rec))
        print(f"   {'ok' if rec.get('ok') else 'FAIL'} wall={rec.get('wall_s')} s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--one", help="JSON of one run (internal)")
    ap.add_argument("--plan", help="JSON file holding a list of runs")
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()
    if a.one:
        item = json.loads(a.one)
        rec = run_one(item["dataset"], item["law"], item["variant"], item["mode"], item["run"], item.get("extra", {}))
        print(json.dumps(rec))
        return
    plan = json.loads(Path(a.plan).read_text())
    drive(plan, a.python)


if __name__ == "__main__":
    main()
```

## Appendix C: `q2_report.py`

```python
"""Tables from q2_results.jsonl: identity, timing, and the Q3 split."""

from __future__ import annotations

import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

RESULTS = Path(sys.argv[1] if len(sys.argv) > 1 else "q2_results.jsonl")
rows = [json.loads(line) for line in RESULTS.read_text().splitlines() if line.strip()]
ok = [r for r in rows if r.get("ok")]
failed = [r for r in rows if not r.get("ok")]

print(f"{len(rows)} runs recorded, {len(ok)} ok, {len(failed)} failed")
for r in failed:
    print(f"  FAILED {r['dataset']} {r['law']} {r['variant']} {r['mode']} r{r['run']}: "
          f"{(r.get('error') or '')[-300:]!r}")
print()

# ---- identity: within one (dataset, law, extra), every run must carry the
# same store hashes and the same mapping checksum, across variants and modes.
print("== Identity across variants (store hashes and mapping checksum) ==")
groups = defaultdict(list)
for r in ok:
    groups[(r["dataset"], r["law"], json.dumps(r.get("extra", {}), sort_keys=True))].append(r)
for key, rs in sorted(groups.items()):
    stores = {json.dumps(r["store"], sort_keys=True) for r in rs}
    sums = {(r["counts"]["checksum"], r["counts"]["n_mapped"]) for r in rs}
    variants = sorted({r["variant"] for r in rs})
    verified = sum(r["counts"]["verified_calls"] for r in rs if r["mode"] == "verify")
    fallback = sum(r["counts"]["fallback_calls"] for r in rs)
    flag = "OK " if len(stores) == 1 and len(sums) == 1 else "MISMATCH"
    print(f"  {flag} {key[0]:24s} {key[1]:14s} runs={len(rs):2d} variants={variants} "
          f"verified_calls={verified:,} fallback_calls={fallback}")
print()

# ---- the plan each dataset resolved to
print("== Resolved axis per dataset (auto law) ==")
seen = set()
for r in ok:
    if r["law"] != "auto" or r["dataset"] in seen or not r["plan"]:
        continue
    seen.add(r["dataset"])
    p = r["plan"]
    c = r["counts"]
    print(f"  {r['dataset']:24s} {p['axis_type']:14s} {p['target_bins']:>10,} bins "
          f"[{p['min_mz']:.1f}, {p['max_mz']:.1f}] {p['method']}; delta={r['delta']:.2e}; "
          f"conv {c['conv_calls']:,} calls / {c['conv_probes']:,} probes; "
          f"sib {c['sib_calls']:,} calls / {c['sib_probes']:,} probes")
print()

# ---- timing
print("== Wall clock, time mode: run 1 (first touch) and median of runs >= 2 (warm), s ==")
print(f"  {'dataset':24s} {'law':14s} {'variant':8s} {'r1':>8s} {'n':>2s} {'median':>8s} {'min':>8s} {'max':>8s} {'vs current':>11s}")
timing = defaultdict(dict)
for key, rs in sorted(groups.items()):
    per_variant = defaultdict(list)
    for r in rs:
        if r["mode"] == "time":
            per_variant[r["variant"]].append(r)
    base_med = None
    if "current" in per_variant:
        warm = [r["wall_s"] for r in per_variant["current"] if r["run"] >= 2]
        base_med = statistics.median(warm) if warm else None
    for v in ("current", "closed2", "closed6"):
        if v not in per_variant:
            continue
        rs_v = per_variant[v]
        r1 = [r["wall_s"] for r in rs_v if r["run"] == 1]
        warm = [r["wall_s"] for r in rs_v if r["run"] >= 2]
        med = statistics.median(warm) if warm else None
        timing[key][v] = (med, warm)
        rel = ""
        if med and base_med:
            rel = f"{100 * (base_med - med) / base_med:+.1f}% ({base_med / med:.2f}x)"
        print(f"  {key[0]:24s} {key[1]:14s} {v:8s} "
              f"{(r1[0] if r1 else float('nan')):8.1f} {len(warm):2d} "
              f"{(med if med else float('nan')):8.1f} "
              f"{(min(warm) if warm else float('nan')):8.1f} {(max(warm) if warm else float('nan')):8.1f} {rel:>11s}")
print()

# ---- Q3: the split of mapping time between the converter and sibling sites
print("== Profile mode: mapping time by call-site group (current variant) ==")
print(f"  {'dataset':24s} {'wall':>8s} {'conv s':>8s} {'conv %':>7s} {'sib s':>8s} {'sib %':>7s} {'conv calls':>11s} {'sib calls':>10s} {'conv probes':>13s} {'sib probes':>13s}")
for r in ok:
    if r["mode"] != "profile":
        continue
    c = r["counts"]
    w = r["wall_s"]
    print(f"  {r['dataset']:24s} {w:8.1f} {c['t_conv']:8.1f} {100 * c['t_conv'] / w:6.1f}% "
          f"{c['t_sib']:8.1f} {100 * c['t_sib'] / w:6.1f}% {c['conv_calls']:11,} {c['sib_calls']:10,} "
          f"{c['conv_probes']:13,} {c['sib_probes']:13,}")


# ---- markdown tables for the write-up
def _med(xs):
    return statistics.median(xs) if xs else None


def _fmt(x, nd=1):
    return "n/a" if x is None else f"{x:,.{nd}f}"


if "--md" in sys.argv:
    print("\n\n## Markdown: per-dataset, auto law\n")
    print("| dataset | axis | bins | mean pts/spec | current r1 (first touch) | current warm median (n) | closed2 warm median (n) | gain | closed6 warm median (n) | gain |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for key, rs in sorted(groups.items()):
        if key[1] != "auto":
            continue
        plan = next((r["plan"] for r in rs if r["plan"]), None)
        c = next((r["counts"] for r in rs if r["mode"] == "time"), None)
        if not plan or not c:
            continue
        pts = c["conv_probes"] / max(c["conv_calls"], 1)
        tv = {v: [r["wall_s"] for r in rs if r["mode"] == "time" and r["variant"] == v and r["run"] >= 2]
              for v in ("current", "closed2", "closed6")}
        r1 = [r["wall_s"] for r in rs if r["mode"] == "time" and r["variant"] == "current" and r["run"] == 1]
        cur = _med(tv["current"])
        cells = [key[0], plan["axis_type"], f"{plan['target_bins']:,}", f"{pts:,.0f}", _fmt(r1[0]) if r1 else "n/a",
                 f"{_fmt(cur)} ({len(tv['current'])})"]
        for v in ("closed2", "closed6"):
            m = _med(tv[v])
            gain = f"{100 * (cur - m) / cur:+.1f}% ({cur / m:.2f}x)" if (m and cur) else "n/a"
            cells += [f"{_fmt(m)} ({len(tv[v])})", gain]
        print("| " + " | ".join(cells) + " |")

    print("\n## Markdown: spread of the warm runs (min-max), s\n")
    print("| dataset | law | current | closed2 | closed6 |")
    print("|---|---|---|---|---|")
    for key, rs in sorted(groups.items()):
        cells = [key[0], key[1]]
        for v in ("current", "closed2", "closed6"):
            xs = [r["wall_s"] for r in rs if r["mode"] == "time" and r["variant"] == v and r["run"] >= 2]
            cells.append(f"{min(xs):.1f}-{max(xs):.1f} ({len(xs)})" if xs else "n/a")
        print("| " + " | ".join(cells) + " |")

    print("\n## Markdown: every law forced\n")
    print("| dataset | law | bins | deviation (bins) | current median (n) | closed2 median (n) | gain |")
    print("|---|---|---|---|---|---|---|")
    for key, rs in sorted(groups.items()):
        if key[1] == "auto":
            continue
        plan = next((r["plan"] for r in rs if r["plan"]), None)
        delta = next((r["delta"] for r in rs if r["delta"] is not None), None)
        cur = [r["wall_s"] for r in rs if r["mode"] == "time" and r["variant"] == "current"]
        c2 = [r["wall_s"] for r in rs if r["mode"] == "time" and r["variant"] == "closed2"]
        mc, m2 = _med(cur), _med(c2)
        gain = f"{100 * (mc - m2) / mc:+.1f}% ({mc / m2:.2f}x)" if (mc and m2) else "n/a"
        print(f"| {key[0]} | {key[1]} | {plan['target_bins']:,} | {delta:.2e} | {_fmt(mc)} ({len(cur)}) | {_fmt(m2)} ({len(c2)}) | {gain} |" if plan else f"| {key[0]} | {key[1]} | failed | | | | |")

    print("\n## Markdown: profile split (current variant, warm)\n")
    print("| dataset | wall s | converter sites s (share) | sibling sites s (share) | converter calls / probes | sibling calls / probes |")
    print("|---|---|---|---|---|---|")
    for r in ok:
        if r["mode"] != "profile":
            continue
        c = r["counts"]; w = r["wall_s"]
        print(f"| {r['dataset']} | {w:,.1f} | {c['t_conv']:,.1f} ({100 * c['t_conv'] / w:.1f}%) | {c['t_sib']:,.1f} ({100 * c['t_sib'] / w:.1f}%) | {c['conv_calls']:,} / {c['conv_probes']:,} | {c['sib_calls']:,} / {c['sib_probes']:,} |")
```

## Appendix D: every run

| dataset | law | variant | mode | run | wall s | cpu s | finished (UTC) |
|---|---|---|---|---|---|---|---|
| mock1273 | auto | current | time | 1 | 33.8 | 33.5 | 16:20:13 |
| mock1273 | auto | closed6 | verify | 1 | 33.8 | 33.7 | 16:20:51 |
| mock1273 | auto | closed2 | time | 1 | 19.6 | 20.8 | 16:21:14 |
| mock1273 | auto | closed6 | time | 1 | 18.0 | 19.5 | 16:21:36 |
| tims01_msms_315px | auto | current | time | 1 | 50.4 | 82.4 | 16:22:38 |
| tims01_msms_315px | auto | closed6 | verify | 1 | 30.5 | 55.8 | 16:23:15 |
| tims01_msms_315px | auto | closed2 | time | 1 | 34.4 | 61.6 | 16:23:53 |
| tims01_msms_315px | auto | closed6 | time | 1 | 36.0 | 63.2 | 16:24:33 |
| tims02_longramp_1465px | auto | current | time | 1 | 31.1 | 62.1 | 16:25:09 |
| tims02_longramp_1465px | auto | closed6 | verify | 1 | 60.5 | 81.5 | 16:26:13 |
| tims02_longramp_1465px | auto | closed2 | time | 1 | 29.3 | 59.2 | 16:26:46 |
| tims02_longramp_1465px | auto | closed6 | time | 1 | 35.8 | 68.2 | 16:27:26 |
| pea_imzml | auto | current | time | 1 | 40.0 | 38.3 | 16:28:11 |
| pea_imzml | auto | closed6 | verify | 1 | 26.6 | 29.2 | 16:28:44 |
| pea_imzml | auto | closed2 | time | 1 | 18.7 | 21.5 | 16:29:18 |
| pea_imzml | auto | closed6 | time | 1 | 20.1 | 24.7 | 16:29:43 |
| bellini_imzml | auto | current | time | 1 | 54.6 | 35.0 | 16:30:45 |
| bellini_imzml | auto | closed6 | verify | 1 | 54.5 | 36.8 | 16:31:59 |
| bellini_imzml | auto | closed2 | time | 1 | 14.7 | 16.6 | 16:32:18 |
| bellini_imzml | auto | closed6 | time | 1 | 11.5 | 13.2 | 16:32:34 |
| tims03_biofilm_20um | auto | current | time | 1 | 110.4 | 355.9 | 16:34:31 |
| tims03_biofilm_20um | auto | closed6 | verify | 1 | 171.2 | 470.3 | 16:37:30 |
| tims03_biofilm_20um | auto | closed2 | time | 1 | 90.8 | 310.3 | 16:39:07 |
| tims03_biofilm_20um | auto | closed6 | time | 1 | 99.9 | 305.8 | 16:40:54 |
| pea_nedc_d | auto | current | time | 1 | 56.6 | 171.2 | 16:41:55 |
| pea_nedc_d | auto | closed6 | verify | 1 | 52.3 | 174.1 | 16:42:52 |
| pea_nedc_d | auto | closed2 | time | 1 | 40.0 | 124.2 | 16:43:40 |
| pea_nedc_d | auto | closed6 | time | 1 | 31.0 | 110.4 | 16:44:16 |
| tims04_ratbrain_71k | auto | current | time | 1 | 146.6 | 257.0 | 16:46:47 |
| tims04_ratbrain_71k | auto | closed6 | verify | 1 | 142.2 | 251.4 | 16:49:15 |
| tims04_ratbrain_71k | auto | closed2 | time | 1 | 122.6 | 222.7 | 16:51:23 |
| tims04_ratbrain_71k | auto | closed6 | time | 1 | 138.2 | 253.7 | 16:53:47 |
| tims05_shortramp_26k | auto | current | time | 1 | 620.8 | 2183.7 | 17:04:37 |
| tims05_shortramp_26k | auto | closed6 | verify | 1 | 613.1 | 2153.5 | 17:15:16 |
| tims05_shortramp_26k | auto | closed2 | time | 1 | 389.2 | 1430.4 | 17:22:06 |
| tims05_shortramp_26k | auto | closed6 | time | 1 | 329.6 | 1198.0 | 17:27:56 |
| xenium_imzml | auto | current | time | 1 | 425.4 | 445.2 | 17:35:28 |
| xenium_imzml | auto | closed6 | verify | 1 | 471.3 | 495.3 | 17:43:53 |
| xenium_imzml | auto | closed2 | time | 1 | 278.5 | 305.5 | 17:48:58 |
| xenium_imzml | auto | closed6 | time | 1 | 279.6 | 307.9 | 17:54:11 |
| xenium_d | auto | current | time | 1 | 679.3 | 2570.1 | 18:05:57 |
| xenium_d | auto | closed6 | verify | 1 | 694.4 | 2673.6 | 18:18:04 |
| xenium_d | auto | closed2 | time | 1 | 479.9 | 1829.4 | 18:26:30 |
| xenium_d | auto | closed6 | time | 1 | 470.0 | 1785.9 | 18:34:47 |
| tims06_glycans_66k | auto | current | time | 1 | 2041.0 | 7678.6 | 19:10:38 |
| tims06_glycans_66k | auto | closed6 | verify | 1 | 2032.4 | 7709.6 | 19:46:20 |
| tims06_glycans_66k | auto | closed2 | time | 1 | 1695.2 | 6383.0 | 20:16:25 |
| tims06_glycans_66k | auto | closed6 | time | 1 | 1639.4 | 6152.2 | 20:45:33 |
| tims01_msms_315px | auto | current | profile | 1 | 28.6 | 48.5 | 20:46:12 |
| tims02_longramp_1465px | auto | current | profile | 1 | 27.5 | 55.2 | 20:46:42 |
| tims03_biofilm_20um | auto | current | profile | 1 | 78.8 | 261.6 | 20:48:06 |
| pea_nedc_d | auto | current | profile | 1 | 32.1 | 117.2 | 20:48:43 |
| tims04_ratbrain_71k | auto | current | profile | 1 | 127.4 | 228.0 | 20:50:55 |
| tims05_shortramp_26k | auto | current | profile | 1 | 365.1 | 1337.5 | 20:57:20 |
| xenium_d | auto | current | profile | 1 | 641.9 | 2445.6 | 21:08:36 |
| tims06_glycans_66k | auto | current | profile | 1 | 1947.1 | 7314.1 | 21:42:44 |
| xenium_imzml | auto | current | profile | 1 | 422.1 | 444.0 | 21:50:14 |
| pea_imzml | auto | current | profile | 1 | 16.0 | 17.3 | 21:50:41 |
| mock1273 | auto | current | time | 2 | 17.7 | 19.3 | 21:51:02 |
| mock1273 | auto | closed2 | time | 2 | 13.3 | 14.3 | 21:51:19 |
| mock1273 | auto | closed6 | time | 2 | 13.3 | 14.8 | 21:51:35 |
| tims01_msms_315px | auto | current | time | 2 | 27.5 | 48.3 | 21:52:06 |
| tims01_msms_315px | auto | closed2 | time | 2 | 26.8 | 47.4 | 21:52:36 |
| tims01_msms_315px | auto | closed6 | time | 2 | 26.8 | 47.7 | 21:53:06 |
| tims02_longramp_1465px | auto | current | time | 2 | 27.2 | 54.4 | 21:53:36 |
| tims02_longramp_1465px | auto | closed2 | time | 2 | 26.6 | 54.0 | 21:54:06 |
| tims02_longramp_1465px | auto | closed6 | time | 2 | 26.4 | 51.8 | 21:54:35 |
| pea_imzml | auto | current | time | 2 | 14.2 | 16.5 | 21:54:53 |
| pea_imzml | auto | closed2 | time | 2 | 9.6 | 12.5 | 21:55:07 |
| pea_imzml | auto | closed6 | time | 2 | 9.5 | 13.2 | 21:55:20 |
| bellini_imzml | auto | current | time | 2 | 18.6 | 14.8 | 21:55:42 |
| bellini_imzml | auto | closed2 | time | 2 | 9.1 | 11.0 | 21:55:55 |
| bellini_imzml | auto | closed6 | time | 2 | 9.4 | 11.1 | 21:56:08 |
| tims03_biofilm_20um | auto | current | time | 2 | 79.9 | 264.3 | 21:57:33 |
| tims03_biofilm_20um | auto | closed2 | time | 2 | 75.6 | 249.7 | 21:58:53 |
| tims03_biofilm_20um | auto | closed6 | time | 2 | 69.8 | 228.8 | 22:00:08 |
| pea_nedc_d | auto | current | time | 2 | 32.5 | 114.6 | 22:00:45 |
| pea_nedc_d | auto | closed2 | time | 2 | 22.7 | 80.8 | 22:01:21 |
| pea_nedc_d | auto | closed6 | time | 2 | 23.3 | 82.3 | 22:01:48 |
| tims04_ratbrain_71k | auto | current | time | 2 | 128.0 | 227.8 | 22:04:01 |
| tims04_ratbrain_71k | auto | closed2 | time | 2 | 117.3 | 211.6 | 22:06:03 |
| tims04_ratbrain_71k | auto | closed6 | time | 2 | 112.8 | 204.8 | 22:08:01 |
| tims05_shortramp_26k | auto | current | time | 2 | 366.7 | 1337.7 | 22:14:28 |
| tims05_shortramp_26k | auto | closed2 | time | 2 | 327.3 | 1186.2 | 22:20:16 |
| tims05_shortramp_26k | auto | closed6 | time | 2 | 316.7 | 1142.0 | 22:25:52 |
| xenium_imzml | auto | current | time | 2 | 408.3 | 439.0 | 22:33:13 |
| xenium_imzml | auto | closed2 | time | 2 | 276.7 | 306.3 | 22:38:16 |
| xenium_imzml | auto | closed6 | time | 2 | 276.2 | 306.6 | 22:43:19 |
| xenium_d | auto | current | time | 2 | 622.5 | 2374.8 | 22:54:08 |
| xenium_d | auto | closed2 | time | 2 | 478.8 | 1818.1 | 23:02:41 |
| xenium_d | auto | closed6 | time | 2 | 477.7 | 1813.9 | 23:11:05 |
| tims06_glycans_66k | auto | current | time | 2 | 2040.0 | 7660.0 | 23:46:46 |
| tims06_glycans_66k | auto | closed2 | time | 2 | 1706.9 | 6428.2 | 00:17:06 |
| tims06_glycans_66k | auto | closed6 | time | 2 | 1648.7 | 6189.3 | 00:46:23 |
| mock1273 | auto | current | time | 3 | 17.6 | 18.4 | 00:46:52 |
| mock1273 | auto | closed2 | time | 3 | 13.7 | 15.1 | 00:47:09 |
| mock1273 | auto | closed6 | time | 3 | 13.5 | 14.6 | 00:47:25 |
| tims01_msms_315px | auto | current | time | 3 | 27.3 | 46.9 | 00:47:55 |
| tims01_msms_315px | auto | closed2 | time | 3 | 26.9 | 48.7 | 00:48:25 |
| tims01_msms_315px | auto | closed6 | time | 3 | 26.3 | 47.1 | 00:48:55 |
| tims02_longramp_1465px | auto | current | time | 3 | 27.0 | 53.4 | 00:49:25 |
| tims02_longramp_1465px | auto | closed2 | time | 3 | 26.7 | 53.5 | 00:49:55 |
| tims02_longramp_1465px | auto | closed6 | time | 3 | 26.1 | 52.3 | 00:50:24 |
| pea_imzml | auto | current | time | 3 | 15.8 | 18.0 | 00:50:43 |
| pea_imzml | auto | closed2 | time | 3 | 9.6 | 13.0 | 00:50:57 |
| pea_imzml | auto | closed6 | time | 3 | 9.7 | 14.5 | 00:51:10 |
| bellini_imzml | auto | current | time | 3 | 19.2 | 15.5 | 00:51:33 |
| bellini_imzml | auto | closed2 | time | 3 | 9.6 | 11.2 | 00:51:46 |
| bellini_imzml | auto | closed6 | time | 3 | 9.4 | 11.0 | 00:51:59 |
| tims03_biofilm_20um | auto | current | time | 3 | 80.0 | 264.2 | 00:53:24 |
| tims03_biofilm_20um | auto | closed2 | time | 3 | 71.4 | 232.6 | 00:54:41 |
| tims03_biofilm_20um | auto | closed6 | time | 3 | 67.8 | 220.1 | 00:55:54 |
| pea_nedc_d | auto | current | time | 3 | 31.7 | 115.0 | 00:56:30 |
| pea_nedc_d | auto | closed2 | time | 3 | 22.8 | 80.5 | 00:56:57 |
| pea_nedc_d | auto | closed6 | time | 3 | 21.2 | 75.2 | 00:57:23 |
| tims04_ratbrain_71k | auto | current | time | 3 | 126.5 | 227.9 | 00:59:34 |
| tims04_ratbrain_71k | auto | closed2 | time | 3 | 118.5 | 214.4 | 01:01:37 |
| tims04_ratbrain_71k | auto | closed6 | time | 3 | 114.5 | 206.5 | 01:03:37 |
| tims05_shortramp_26k | auto | current | time | 3 | 364.8 | 1333.6 | 01:10:10 |
| tims05_shortramp_26k | auto | closed2 | time | 3 | 326.2 | 1190.2 | 01:15:56 |
| tims05_shortramp_26k | auto | closed6 | time | 3 | 313.1 | 1130.9 | 01:21:29 |
| xenium_imzml | auto | current | time | 3 | 420.4 | 440.3 | 01:29:04 |
| xenium_imzml | auto | closed2 | time | 3 | 275.8 | 307.3 | 01:34:06 |
| xenium_imzml | auto | closed6 | time | 3 | 275.6 | 307.2 | 01:39:16 |
| xenium_d | auto | current | time | 3 | 628.4 | 2398.6 | 01:50:11 |
| xenium_d | auto | closed2 | time | 3 | 475.8 | 1808.7 | 01:58:40 |
| xenium_d | auto | closed6 | time | 3 | 477.3 | 1815.3 | 02:07:04 |
| tims06_glycans_66k | auto | current | time | 3 | 1963.6 | 7378.5 | 02:41:29 |
| tims06_glycans_66k | auto | closed2 | time | 3 | 1724.9 | 6502.5 | 03:12:00 |
| tims06_glycans_66k | auto | closed6 | time | 3 | 1616.5 | 6065.4 | 03:40:44 |
| mock1273 | auto | current | time | 4 | 18.0 | 19.2 | 03:41:13 |
| mock1273 | auto | closed2 | time | 4 | 13.4 | 14.4 | 03:41:30 |
| mock1273 | auto | closed6 | time | 4 | 13.4 | 14.9 | 03:41:47 |
| tims01_msms_315px | auto | current | time | 4 | 27.3 | 47.5 | 03:42:17 |
| tims01_msms_315px | auto | closed2 | time | 4 | 26.6 | 47.2 | 03:42:46 |
| tims01_msms_315px | auto | closed6 | time | 4 | 26.1 | 47.5 | 03:43:15 |
| tims02_longramp_1465px | auto | current | time | 4 | 27.0 | 54.1 | 03:43:46 |
| tims02_longramp_1465px | auto | closed2 | time | 4 | 26.5 | 53.3 | 03:44:15 |
| tims02_longramp_1465px | auto | closed6 | time | 4 | 26.8 | 54.3 | 03:44:45 |
| pea_imzml | auto | current | time | 4 | 16.2 | 18.4 | 03:45:05 |
| pea_imzml | auto | closed2 | time | 4 | 9.9 | 13.5 | 03:45:19 |
| pea_imzml | auto | closed6 | time | 4 | 9.8 | 13.2 | 03:45:33 |
| bellini_imzml | auto | current | time | 4 | 19.0 | 15.3 | 03:45:55 |
| bellini_imzml | auto | closed2 | time | 4 | 9.5 | 11.4 | 03:46:08 |
| bellini_imzml | auto | closed6 | time | 4 | 9.4 | 11.5 | 03:46:21 |
| tims03_biofilm_20um | auto | current | time | 4 | 81.4 | 272.9 | 03:47:47 |
| tims03_biofilm_20um | auto | closed2 | time | 4 | 71.9 | 235.5 | 03:49:05 |
| tims03_biofilm_20um | auto | closed6 | time | 4 | 69.7 | 226.7 | 03:50:20 |
| pea_nedc_d | auto | current | time | 4 | 31.7 | 113.8 | 03:50:56 |
| pea_nedc_d | auto | closed2 | time | 4 | 22.0 | 78.5 | 03:51:22 |
| pea_nedc_d | auto | closed6 | time | 4 | 21.8 | 77.3 | 03:51:48 |
| tims04_ratbrain_71k | auto | current | time | 4 | 127.1 | 227.1 | 03:54:00 |
| tims04_ratbrain_71k | auto | closed2 | time | 4 | 118.0 | 212.9 | 03:56:03 |
| tims04_ratbrain_71k | auto | closed6 | time | 4 | 112.8 | 205.8 | 03:58:00 |
| tims05_shortramp_26k | auto | current | time | 4 | 363.5 | 1325.7 | 04:04:24 |
| tims05_shortramp_26k | auto | closed2 | time | 4 | 326.3 | 1183.9 | 04:10:18 |
| tims05_shortramp_26k | auto | closed6 | time | 4 | 314.2 | 1135.1 | 04:15:52 |
| mock1273 | auto | current | time | 5 | 17.6 | 18.9 | 04:16:14 |
| mock1273 | auto | closed2 | time | 5 | 13.5 | 15.0 | 04:16:31 |
| mock1273 | auto | closed6 | time | 5 | 13.6 | 14.5 | 04:16:47 |
| tims01_msms_315px | auto | current | time | 5 | 28.2 | 47.4 | 04:17:18 |
| tims01_msms_315px | auto | closed2 | time | 5 | 27.0 | 46.7 | 04:17:48 |
| tims01_msms_315px | auto | closed6 | time | 5 | 26.5 | 47.8 | 04:18:18 |
| tims02_longramp_1465px | auto | current | time | 5 | 26.7 | 54.0 | 04:18:48 |
| tims02_longramp_1465px | auto | closed2 | time | 5 | 26.8 | 55.0 | 04:19:18 |
| tims02_longramp_1465px | auto | closed6 | time | 5 | 26.5 | 52.2 | 04:19:47 |
| pea_imzml | auto | current | time | 5 | 14.5 | 17.1 | 04:20:06 |
| pea_imzml | auto | closed2 | time | 5 | 9.8 | 13.2 | 04:20:19 |
| pea_imzml | auto | closed6 | time | 5 | 9.8 | 13.8 | 04:20:33 |
| bellini_imzml | auto | current | time | 5 | 13.6 | 15.4 | 04:20:50 |
| bellini_imzml | auto | closed2 | time | 5 | 9.1 | 11.2 | 04:21:03 |
| bellini_imzml | auto | closed6 | time | 5 | 9.3 | 11.0 | 04:21:15 |
| tims03_biofilm_20um | auto | current | time | 5 | 80.0 | 269.1 | 04:22:40 |
| tims03_biofilm_20um | auto | closed2 | time | 5 | 72.8 | 239.4 | 04:23:58 |
| tims03_biofilm_20um | auto | closed6 | time | 5 | 69.2 | 222.8 | 04:25:13 |
| pea_nedc_d | auto | current | time | 5 | 29.5 | 107.5 | 04:25:47 |
| pea_nedc_d | auto | closed2 | time | 5 | 22.3 | 79.3 | 04:26:13 |
| pea_nedc_d | auto | closed6 | time | 5 | 22.7 | 81.5 | 04:26:40 |
| tims04_ratbrain_71k | auto | current | time | 5 | 128.6 | 231.8 | 04:28:54 |
| tims04_ratbrain_71k | auto | closed2 | time | 5 | 120.1 | 219.6 | 04:30:59 |
| tims04_ratbrain_71k | auto | closed6 | time | 5 | 115.9 | 210.9 | 04:32:59 |
| tims05_shortramp_26k | auto | current | time | 5 | 367.4 | 1341.4 | 04:39:34 |
| tims05_shortramp_26k | auto | closed2 | time | 5 | 328.8 | 1201.8 | 04:45:31 |
| tims05_shortramp_26k | auto | closed6 | time | 5 | 308.3 | 1114.4 | 04:50:59 |
| pea_imzml | constant | closed6 | verify | 1 | 15.4 | 18.3 | 04:51:19 |
| pea_imzml | constant | current | time | 1 | 13.1 | 17.5 | 04:51:36 |
| pea_imzml | constant | closed2 | time | 1 | 8.9 | 12.4 | 04:51:48 |
| pea_imzml | constant | current | time | 2 | 13.4 | 16.7 | 04:52:06 |
| pea_imzml | constant | closed2 | time | 2 | 9.1 | 12.5 | 04:52:19 |
| pea_imzml | constant | current | time | 3 | 13.0 | 16.5 | 04:52:35 |
| pea_imzml | constant | closed2 | time | 3 | 9.0 | 13.5 | 04:52:48 |
| pea_imzml | linear_tof | closed6 | verify | 1 | 12.9 | 16.0 | 04:53:05 |
| pea_imzml | linear_tof | current | time | 1 | 11.3 | 14.8 | 04:53:20 |
| pea_imzml | linear_tof | closed2 | time | 1 | 8.6 | 11.5 | 04:53:32 |
| pea_imzml | linear_tof | current | time | 2 | 11.4 | 14.8 | 04:53:48 |
| pea_imzml | linear_tof | closed2 | time | 2 | 9.0 | 12.2 | 04:54:00 |
| pea_imzml | linear_tof | current | time | 3 | 11.6 | 14.5 | 04:54:16 |
| pea_imzml | linear_tof | closed2 | time | 3 | 9.0 | 11.4 | 04:54:28 |
| pea_imzml | reflector_tof | closed6 | verify | 1 | 16.3 | 18.5 | 04:54:49 |
| pea_imzml | reflector_tof | current | time | 1 | 13.9 | 17.1 | 04:55:06 |
| pea_imzml | reflector_tof | closed2 | time | 1 | 9.6 | 13.4 | 04:55:20 |
| pea_imzml | reflector_tof | current | time | 2 | 13.9 | 18.0 | 04:55:37 |
| pea_imzml | reflector_tof | closed2 | time | 2 | 9.7 | 13.7 | 04:55:51 |
| pea_imzml | reflector_tof | current | time | 3 | 14.1 | 17.4 | 04:56:09 |
| pea_imzml | reflector_tof | closed2 | time | 3 | 9.7 | 13.2 | 04:56:22 |
| pea_imzml | tof | closed6 | verify | 1 | 16.3 | 20.0 | 04:56:42 |
| pea_imzml | tof | current | time | 1 | 12.6 | 16.6 | 04:56:59 |
| pea_imzml | tof | closed2 | time | 1 | 10.3 | 14.1 | 04:57:13 |
| pea_imzml | tof | current | time | 2 | 12.8 | 15.4 | 04:57:30 |
| pea_imzml | tof | closed2 | time | 2 | 10.0 | 13.3 | 04:57:43 |
| pea_imzml | tof | current | time | 3 | 13.0 | 16.5 | 04:58:00 |
| pea_imzml | tof | closed2 | time | 3 | 9.9 | 13.6 | 04:58:14 |
| pea_imzml | orbitrap | closed6 | verify | 1 | 16.7 | 21.2 | 04:58:35 |
| pea_imzml | orbitrap | current | time | 1 | 15.3 | 19.5 | 04:58:54 |
| pea_imzml | orbitrap | closed2 | time | 1 | 10.4 | 14.9 | 04:59:08 |
| pea_imzml | orbitrap | current | time | 2 | 15.1 | 19.6 | 04:59:27 |
| pea_imzml | orbitrap | closed2 | time | 2 | 10.1 | 13.5 | 04:59:41 |
| pea_imzml | orbitrap | current | time | 3 | 15.1 | 18.9 | 05:00:00 |
| pea_imzml | orbitrap | closed2 | time | 3 | 10.1 | 14.8 | 05:00:14 |
| pea_imzml | fticr | closed6 | verify | 1 | 18.0 | 20.8 | 05:00:36 |
| pea_imzml | fticr | current | time | 1 | 16.0 | 21.3 | 05:00:56 |
| pea_imzml | fticr | closed2 | time | 1 | 10.5 | 13.4 | 05:01:10 |
| pea_imzml | fticr | current | time | 2 | 16.3 | 20.3 | 05:01:30 |
| pea_imzml | fticr | closed2 | time | 2 | 10.2 | 13.7 | 05:01:44 |
| pea_imzml | fticr | current | time | 3 | 16.2 | 19.1 | 05:02:04 |
| pea_imzml | fticr | closed2 | time | 3 | 10.3 | 14.1 | 05:02:18 |
| tims03_biofilm_20um | constant | closed6 | verify | 1 | 79.7 | 263.5 | 05:03:43 |
| tims03_biofilm_20um | constant | current | time | 1 | 76.0 | 250.8 | 05:05:11 |
| tims03_biofilm_20um | constant | closed2 | time | 1 | 70.9 | 233.1 | 05:06:27 |
| tims03_biofilm_20um | constant | current | time | 2 | 75.8 | 251.7 | 05:07:48 |
| tims03_biofilm_20um | constant | closed2 | time | 2 | 69.6 | 226.3 | 05:09:02 |
| tims03_biofilm_20um | constant | current | time | 3 | 79.0 | 265.7 | 05:10:26 |
| tims03_biofilm_20um | constant | closed2 | time | 3 | 72.2 | 236.8 | 05:11:43 |
| tims03_biofilm_20um | linear_tof | closed6 | verify | 1 | 76.6 | 254.7 | 05:13:04 |
| tims03_biofilm_20um | linear_tof | current | time | 1 | 72.9 | 238.3 | 05:14:22 |
| tims03_biofilm_20um | linear_tof | closed2 | time | 1 | 71.9 | 234.7 | 05:15:38 |
| tims03_biofilm_20um | linear_tof | current | time | 2 | 72.0 | 236.6 | 05:16:54 |
| tims03_biofilm_20um | linear_tof | closed2 | time | 2 | 70.9 | 228.5 | 05:18:10 |
| tims03_biofilm_20um | linear_tof | current | time | 3 | 74.7 | 245.4 | 05:19:36 |
| tims03_biofilm_20um | linear_tof | closed2 | time | 3 | 70.0 | 229.1 | 05:20:50 |
| tims03_biofilm_20um | reflector_tof | closed6 | verify | 1 | 84.5 | 283.8 | 05:22:20 |
| tims03_biofilm_20um | reflector_tof | current | time | 1 | 79.6 | 266.3 | 05:23:45 |
| tims03_biofilm_20um | reflector_tof | closed2 | time | 1 | 73.8 | 242.1 | 05:25:04 |
| tims03_biofilm_20um | reflector_tof | current | time | 2 | 79.7 | 264.8 | 05:26:29 |
| tims03_biofilm_20um | reflector_tof | closed2 | time | 2 | 72.9 | 239.9 | 05:27:47 |
| tims03_biofilm_20um | reflector_tof | current | time | 3 | 80.3 | 268.2 | 05:29:13 |
| tims03_biofilm_20um | reflector_tof | closed2 | time | 3 | 73.1 | 239.2 | 05:30:31 |
| tims03_biofilm_20um | tof | closed6 | verify | 1 | 85.3 | 287.5 | 05:32:02 |
| tims03_biofilm_20um | tof | current | time | 1 | 79.6 | 267.0 | 05:33:27 |
| tims03_biofilm_20um | tof | closed2 | time | 1 | 74.4 | 245.3 | 05:34:46 |
| tims03_biofilm_20um | tof | current | time | 2 | 77.9 | 259.2 | 05:36:09 |
| tims03_biofilm_20um | tof | closed2 | time | 2 | 75.1 | 248.4 | 05:37:29 |
| tims03_biofilm_20um | tof | current | time | 3 | 77.1 | 256.3 | 05:38:52 |
| tims03_biofilm_20um | tof | closed2 | time | 3 | 76.2 | 254.2 | 05:40:13 |
| tims03_biofilm_20um | orbitrap | closed6 | verify | 1 | 88.8 | 301.0 | 05:41:47 |
| tims03_biofilm_20um | orbitrap | current | time | 1 | 82.5 | 275.5 | 05:43:15 |
| tims03_biofilm_20um | orbitrap | closed2 | time | 1 | 73.7 | 244.0 | 05:44:34 |
| tims03_biofilm_20um | orbitrap | current | time | 2 | 81.6 | 274.5 | 05:46:02 |
| tims03_biofilm_20um | orbitrap | closed2 | time | 2 | 74.0 | 243.4 | 05:47:21 |
| tims03_biofilm_20um | orbitrap | current | time | 3 | 80.0 | 268.8 | 05:48:46 |
| tims03_biofilm_20um | orbitrap | closed2 | time | 3 | 73.3 | 242.8 | 05:50:05 |
| tims03_biofilm_20um | fticr | closed6 | verify | 1 | 88.4 | 300.1 | 05:51:38 |
| tims03_biofilm_20um | fticr | current | time | 1 | 83.1 | 279.2 | 05:53:07 |
| tims03_biofilm_20um | fticr | closed2 | time | 1 | 75.1 | 248.5 | 05:54:27 |
| tims03_biofilm_20um | fticr | current | time | 2 | 82.6 | 278.9 | 05:55:55 |
| tims03_biofilm_20um | fticr | closed2 | time | 2 | 77.1 | 254.1 | 05:57:17 |
| tims03_biofilm_20um | fticr | current | time | 3 | 84.2 | 284.5 | 05:58:47 |
| tims03_biofilm_20um | fticr | closed2 | time | 3 | 74.5 | 245.5 | 06:00:07 |
