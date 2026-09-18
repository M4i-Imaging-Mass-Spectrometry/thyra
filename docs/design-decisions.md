# Design Decisions

This page records the decisions behind Thyra's defaults where a reasonable
person could have chosen otherwise. Each entry gives the decision, the
reason, the strongest objection that was raised against it and why it did
not win, and the known limit. The goal is that someone who disagrees can see
exactly which premise to argue with, and open an issue against that premise
rather than against a number in a table.

Every measurement quoted here was taken on real acquisitions on the date
given, so a future maintainer can repeat it.

**Status vocabulary.** *Implemented* means the code on `main` does this.
*Accepted* means the decision is made and recorded but the change has not
shipped yet. *Deferred* means no change, with the condition that would reopen
it stated.

---

## D1. Which spectrum a reader takes

**Decision.** Thyra reads the instrument's sparse record when the file holds
one, otherwise the vendor's picked spectrum. It applies no peak picking of
its own.

**Status:** Implemented 2026-09-07. The Bruker TDF default changed from
`vendor_centroid` to `scan_sum`; the centroid stays available through
`--tdf-spectrum vendor_centroid`. Every other row of the table below was
already what the code does.

| Source | What the file holds | What Thyra reads | Why |
|---|---|---|---|
| Bruker TDF (TIMS on) | per-scan digitizer-index counts for every mobility scan, plus a vendor centroid on request | the scans, summed per index (`scan_sum`) | the sparse record; equals Bruker's own frame total |
| Bruker TSF (TIMS off) | a vendor line spectrum, plus a raw digitizer trace | the line spectrum | the trace is a continuum with a baseline, not a peak list |
| Bruker solariX | centroided peak lists in `peaks.sqlite`; raw transients | the peak lists | the only sparse record; transients need an FT |
| Bruker rapifleX | profile spectra | the profile, resampled | nothing else exists; the resampler conserves current |
| Waters SELECT SERIES MRT | a zero-suppressed digitiser trace, plus a vendor centroid on demand | the trace, onto a digitiser-matched axis | the picker merges 8 to 9 mDa doublets the trace resolves; see [Supported Formats](supported-formats.md#waters-masslynx) |
| every other Waters instrument | the same trace, plus a vendor centroid on demand | the centroid (profile on request) | measured on a Synapt G2-Si, the picker keeps everything the trace resolves, so the trace would cost 2 to 3 times the store for nothing |
| imzML, mzPeak, PHI | whatever was exported | as exported | the export already made the choice |

The line between TDF and TSF is the whole decision, so it was measured
rather than assumed.

**Measured on TDF** (2026-09-07; 30 frames of a 26k-pixel MALDI-2 slide and
10 frames of a 20 um biofilm set):

| Quantity | Slide | Biofilm |
|---|---|---|
| `scan_sum` total against Bruker's quasi-profile export (`tims_extract_profile_for_frame`) | identical | identical |
| `scan_sum` total against Bruker's per-frame TIC column (`Frames.SummedIntensities`) | 1.0000 | 0.9992 |
| vendor centroid total against the same TIC column | 0.8748 | 0.9644 |
| centroid intensity of a strong peak against the summed indices under it | 0.97 to 0.99 | 0.98 to 0.99 |
| points per frame, `scan_sum` over centroid | 2.9x | 1.5x |
| raw pair intensities equal to 1 | none | none |

Three things are therefore true at once. The vendor centroid reports peak
*areas* and conserves the current inside each peak it picks. What it drops
is every index bin that its picker assigned to no peak: 12.5 percent of the
slide's current and 3.6 percent of the biofilm's. And Bruker's own
definition of a frame's total ion current is the scan sum, not the centroid.
The earlier description of the loss as "single-count noise" was wrong: these
files contain no single counts, and the dropped current is sub-threshold
signal.

**Measured on TSF** (2026-09-07; two files from different acquisitions): the
line spectrum's intensity equals the maximum of the digitizer trace under
the peak exactly (ratio 1.000 on every peak checked), so TSF intensities are
peak *heights*, and `Frames.SummedIntensities` is the sum of those heights.
The trace has a baseline of 18 on every sample, peaks 15 to 150 samples
wide, and a total 3.4 to 3.7 times the line sum. Summing it would import a
baseline and change the meaning of intensity from height to area. That is a
different kind of data, so TSF stays on the vendor's line spectrum. The
Waters trace is not the same kind: it is zero-suppressed, so it carries no
baseline, and it is the digitiser's record. Which side of the rule a Waters
instrument falls on is therefore decided by whether its peak picker has been
shown to lose what the trace keeps, which is the MRT measurement in the
next row of the table.

**Why `scan_sum` is the TDF default.** Thyra's one invariant everywhere
else is that ion current is conserved: the resampler is TIC-preserving and
every sibling table is checked against an exact identity. The vendor
centroid was the single place where a closed, unversioned algorithm removed
current before the store was written. Under `scan_sum` every identity holds
by construction: the mobility heatmap's marginal is the mean spectrum, the
mobility grid's marginal is the summed table, and the MS/MS blocks add back
up to the summed table. And the number a reader calls TIC is the number the
vendor calls TIC.

**Objections considered.**

- *The dropped current is noise, and keeping it raises the noise floor.*
  Some of it is. Removing it is analysis, reproducible downstream by any
  method the user chooses and documented in their own pipeline. Removing it
  at conversion is irreversible and undocumented.
- *Users comparing with SCiLS Lab or TSF exports will see different
  numbers.* True. The mode is recorded in
  `msi_metadata.processing[0].parameters.tdf_spectrum`, and
  `--tdf-spectrum vendor_centroid` restores the vendor numbers exactly.
- *It changes the numbers of every existing TDF pipeline.* True, and the
  real cost. It ships as a minor release with a changelog entry, and an old
  store can be told from a new one by the provenance field.
- *Then TSF should sum its profile for consistency.* No. The TSF profile
  is a continuum with a baseline, and its line intensities are heights by
  Bruker's own definition. The rule is about sparse records, not about
  reading the largest array available.
- *Then Waters should stay on its centroid, as this page first said.* The
  first version of this table put every Waters instrument on the centroid
  and called its profile "the same kind of data" as the TSF trace. That was
  wrong on both counts, and the correction came from the measurements
  behind the MRT profile default rather than from this page: the Waters
  trace has no baseline, and on the MRT the picker merges near-isobars the
  trace resolves. The rule did not change; the row did.

**Known limits.** The summed table holds about three times the points per
pixel on a short-ramp slide, half again on a long-ramp one. Measured on
the whole 26,087-pixel slide (2026-09-07, default mass axis, optical image
left out, same machine):

| | `vendor_centroid` | `scan_sum` | ratio |
|---|---|---|---|
| non-zeros in the summed table | 282,223,487 | 893,289,573 | 3.2x |
| store on disk | 1.36 GB | 2.50 GB | 1.8x |
| conversion, warm | 438 s | 528 s | 1.2x |

The store grows less than the point count because the sharded chunks
compress the low counts well. That is the price of the default, and the
flag buys the old size back.

---

## D2. The MS/MS table is written by default when the schedule qualifies

**Status:** Implemented 2026-09-07. `--msms-table` is the default; the
table is written for a Bruker PASEF acquisition with a constant,
non-overlapping schedule of at least two precursors and refused, with the
reason logged, on everything else. `--no-msms-table` opts out.

**Reason.** A scheduled PASEF acquisition fragments each precursor in turn
at every pixel. The one-spectrum-per-pixel summed table of such a pixel is
therefore a mixture of unrelated fragment spectra. It is not a spectrum of
anything, and writing only that mixture is the inaccurate representation.
The split is parameter-free, since the schedule is the instrument's own
list, and lossless, since every recorded point falls in exactly one
isolation window; under D1 its blocks add back to the summed table exactly.
The refusals for a schedule that varies across pixels, overlapping windows,
or a single precursor already fall back silently to the summed table, so the
default never breaks a conversion. The prior-art survey in
[Output Format](output-format.md#fragmentation-msms) found no open
analysis-layer convention to defer to.

**Objections considered.**

- *It has run on one dataset: 713 pixels, 15 precursors, one schedule
  shape.* True when it shipped, and stated as a limit. Answered on
  2026-09-07 by four more acquisitions, below: a second schedule of 13
  precursors in the other polarity, sharing no precursor with the first,
  found by walking the lab share. The assembly engine underneath is the one
  verified on the 26k-pixel mobility grid, and the MS/MS count span is
  tiny, so scale is not where it would fail; the schedule shape was, and
  now there are two of them measured plus a third made by deleting a window
  from a copy.
- *A default-on writes an element the user did not ask for and consumers
  must recognise.* The sibling is discriminated by its `uns` block exactly
  as the mobility sibling already is, and the summed table remains the
  primary element.

### What five more acquisitions showed (2026-09-07)

Every file in `TIMS-test-data\01_tiny_msms_315px`, and the two negative-mode
images the share search below turned up, converted with the defaults and
`--no-optical` on both write routes (`--streaming true` and `false`). The
two routes agreed exactly on every table they wrote -- same `var` labels,
same `obs`, and the same CSR `indices`, `indptr` and `data` arrays -- so
only one set of numbers is given:

| acquisition | pixels | `MsMsType` | precursor table | tables written | windows | `current_ratio` | MS/MS `var` |
|---|---|---|---|---|---|---|---|
| `220425_MSMS_pos_brain1.d` | 713 | 8 | `PasefFrameMsMsInfo`, 10,695 rows | summed + `_msms` | 15 | 1.0 (per pixel 1.0 to 1.0) | 713 x 264,006, 302,460 nnz |
| `220425_MSMS_pos_brain1_2.d` | 668 | 8 | `PasefFrameMsMsInfo`, 10,020 rows | summed + `_msms` | 15 | 1.0 (per pixel 1.0 to 1.0) | 668 x 284,235, 335,840 nnz |
| `220425_MSMS_neg.d` | 487 | 8 | `PasefFrameMsMsInfo`, 6,331 rows | summed + `_msms` | 13 | 1.0 (per pixel 1.0 to 1.0) | 487 x 162,829, 169,079 nnz |
| `220425_MSMS_neg_test.d` | 519 | 8 | `PasefFrameMsMsInfo`, 6,747 rows | summed + `_msms` | 13 | 1.0 (per pixel 1.0 to 1.0) | 519 x 93,926, 96,046 nnz |
| `220425_MSMS_pos.d` | 315 | 2 | `FrameMsMsInfo`, 315 rows | summed only | 1 | -- | -- |

Conversion took 6 to 19 s per file per route. The summed axis is 599,146
bins on the positive files and 635,610 on the negative ones, which acquired
to m/z 1200 rather than 1000 -- the default grid follows the acquisition
range, so **two runs share an axis only when they share that range.**

**The refusal on `220425_MSMS_pos.d` is the right answer.** Its 315 frames
are all `MsMsType = 2`, single-precursor MS/MS, and it has no
`PasefFrameMsMsInfo` at all: `FrameMsMsInfo` carries one row per frame,
every one of them `TriggerMass` 1046.54 at `IsolationWidth` 1.5 and
`CollisionEnergy` 57.327. The whole mobility ramp of every frame fragments
that one precursor, so the summed table already *is* its fragment spectrum
and there is nothing to separate. This is the acquisition the "one
precursor" refusal was written for, met for the first time.

**Two pixel sets of one method align on the labels.** `brain1` and
`brain1_2` share 38,453 of their 264,006 and 284,235 feature labels, and
every shared label agrees on `precursor_mz`, `mz`, `mz_index` *and*
`precursor_index`; `anndata.concat` gives 1,381 x 509,788 with per-pixel
totals preserved and no column holding two precursors. `precursor_mobility`
does **not** agree between those two: they carry their own TIMS
calibration and the same window comes out up to 4.0e-3 apart in 1/K0. The
two negative-mode runs, acquired back to back, agree on it *bitwise*. So
the coordinate is reproducible when the calibration is and not otherwise,
which is why aligning on `(precursor_mz, precursor_mobility)` means
matching the pair rather than comparing it with `==`.

**A second schedule, acquired as such.** The negative-mode pair isolates 13
precursors from 519.182 to 1179.731 -- a list with **nothing in common**
with the positive-mode 15 from 313.275 to 936.578. Concatenating a negative
store with a positive one shares **zero** labels, which is the whole point:
two unrelated acquisitions must not merge. Had the labels been named after
the rank, the two would have shared **5,262** labels and *every one of them*
would have named two different precursors -- `p0_mz31749` is 519.182 in one
and 313.275 in the other -- and the different mass axes would not have
saved it, because both axes start at m/z 50 and their low bins line up. The
two negative runs, which do share a schedule, share 7,155 labels and agree
on every column of `var` including the rank.

**And a third shape, made from the data at hand.** A copy of `brain1_2`
with the lowest-m/z window (313.275) deleted from all 668 frames -- from
every frame, so `constant_across_pixels` still holds and only the shape
changes. The result is 14 precursors, each one rank lower than in `brain1`:
353.32 is `precursor_index` 1 there and 0 here. The alignment survives it.
36,912 labels are shared and all of them agree on `precursor_mz`, `mz` and
`mz_index`; `precursor_index` does not, and
`anndata.concat(join="inner", merge="unique")` drops that column while
keeping the other three -- anndata itself finding the rank disagrees. The
15,741 labels of the deleted precursor are `brain1`-only and the reduced
sample contributes exactly 0 to them; a rank-named scheme would have merged
all 27,145 labels the two stores would then have shared.

**And the current block says what was removed.** The reduced copy's summed
table is untouched, so its split now misses exactly the deleted window's
ion current: `current_ratio` 0.9717425865942857 against
1 - 0.028257413405714304 = 0.9717425865942857, equal to the last digit.

### The share census

Every `analysis.tdf` on the lab share was opened read-only and asked what it
is: **1,021** of them, plus the nine in the local corpus. The share's timsTOF
instrument directory holds 25 -- most `.d` folders there are TSF, TIMS off --
and the other twenty instrument directories hold none.
Of the 1,021, 230 are imaging acquisitions (they have `MaldiFrameInfo`) and
95 carry `PasefFrameMsMsInfo`, but only **six are both**, and those six are
the four positive-mode files above (two of them duplicate copies) plus the
two negative-mode ones. Every other PASEF file on the share is an LC run
with no MALDI geometry at all, so the converter stops there long before the
tables -- a clean conversion failure, not a crash.

Their *schedules* are still worth reading, and reading them exercised two
refusals no real file had reached before.

- **Data-dependent PASEF.** The datasets registered as `tims_dda_pro2_2026`
  are fourteen DDA runs on a timsTOF Pro 2. One
  (`tims_dda_pro2_2026_741`) has 1,840 survey and 9,912 fragment
  frames and **14,952 distinct isolation windows, each in 1 to 8 frames**;
  across the share the count reaches 90,915. The precursors are chosen per
  frame, so there is no global feature axis, and the windows overlap on the
  ramp as well. Refused on the first of those, which is the right one: the
  overlap is a consequence of the design, the varying schedule *is* the
  design.
- **diaPASEF.** 654 files interleave survey and DIA fragment frames -- 1,786
  and about 28,566 of them in one -- with their windows in
  `DiaFrameMsMsWindows`.

**Reading a diaPASEF one found a refusal that lied.**
`DiaFrameMsMsWindows` is not one of the two tables the TDF reader looks in,
so such a file yields a schedule with an empty window list, and
`constant_across_pixels` false because the frame types are mixed.
`demultiplex_refusal` asked the window count before the constancy and
answered *"the acquisition isolates a single precursor"* about a run that
isolates 32. It now asks the conditions in order of what each says about
the acquisition -- not MS/MS, then not constant, then no precursor
recorded, then one precursor, then overlapping -- so that file is refused
as non-constant, which is true, and an empty window list has a reason of
its own. The `INFO` line about a `current_ratio` off 1 was likewise
attributing every deviation to `vendor_centroid`; a schedule whose windows
do not cover every scan carrying current is the other direction, as the
reduced copy above measures.

**Known limits.** Tested on four MALDI PASEF images across **two
independently acquired schedules** -- 15 precursors in positive mode over
m/z 50-1000, 13 in negative over 50-1200 -- at 487 to 713 pixels, plus a
third shape made by deleting one window from a copy, the single-precursor
refusal on a real `MsMsType = 2` acquisition, and the varying-schedule
refusal on real DDA-PASEF. All four come from the same instrument and the
same method family, and none exceeds 713 pixels: the assembly engine
underneath is the one verified on the 26k-pixel mobility grid, but this
table has not itself been run at that scale.

---

## D3. The mobility grid stays opt in

**Status:** Implemented. `--mobility-grid` is off by default.

**Reason.** This is the opposite case to D2. Summing a TIMS pixel's
mobility scans gives an ordinary MS1 spectrum; nothing about it is
inaccurate. The grid is the same data unfolded along a second axis with
choices the summed table does not need: how many mobility channels, and a
channel width that then depends on each acquisition's ramp. Anything with
free parameters is a derived product, and derived products are written when
asked, with their parameters recorded, so that nobody later finds a 256 in a
store and wonders who chose it. The lossless summary of the mobility
dimension, the mass-mobility heatmap plus the mobility axis metadata, is
written by default, and under D1 the heatmap's marginal equals the stored
mean spectrum exactly.

The sentence that decides D2 and D3 together: the default store contains
what is accurate without any parameter beyond the mass axis, plus a
lossless summary of every extra dimension. Grids over an extra dimension
are opt in and carry their parameters.

**Objections considered.**

- *The mass axis is also a parameterised, derived choice, and it is
  default.* A mass axis is necessary to build any table at all; the grid is
  optional on top of a complete table. This is the thinnest argument on the
  page and is recorded as such.
- *Users will not discover it.* The conversion log says the source has
  mobility, and this page exists.
- *Practical costs are secondary but real.* On 400 frames the store is six
  times larger and the conversion about five times slower warm; at the
  default mass axis the whole 26k-pixel slide is refused (see D4). A
  default that refuses the flagship dataset is not a default anyone can
  rely on.

**Known limits.** None beyond D4.

---

## D4. The grid's feature ceiling is a memory guard, not a format limit

**Status:** Implemented 2026-09-07. The fixed ceiling of 20,000,000 occupied
(m/z bin, mobility channel) pairs became a projection of the `var` frame's
memory (330 bytes per feature, measured) against the machine's free memory,
warning past a quarter of it and refusing past half, with an absolute cap of
100,000,000 kept as a statement of what downstream tools can be expected to
open rather than as the operative guard.

**Reason.** Both sibling tables are built out of core, so the remaining
peak is the `var` frame, measured at roughly 330 bytes per feature including
the AnnData copies. The ceiling exists to protect that build. Raising the
constant from 20M to 25M because the flagship slide occupies 21.1M pairs at
the default axis would be fitting a constant to one dataset, and the
objection to that is immediate. The honest guard is the one the grid used
once before and that was accepted then: project bytes from the count and
compare with what the machine has. Fractions rather than sizes, because a
table that is routine on a workstation is fatal on a laptop.

**Objections considered.**

- *The per-feature constant is version and machine dependent.* It is
  measured and documented, and a guess wrong by a factor of two still
  refuses at half of free memory.

**Measured on the whole slide** (2026-09-07, 26,087 pixels at the default
138,629-bin axis, `--mobility-grid --no-optical`, warm, no optical image).
The conversion the fixed ceiling refused now runs: 21,101,151 occupied
`(m/z bin, mobility channel)` pairs, 1,627,659,585 non-zeros through 18.19 GB
of scratch, 239 of 256 channels occupied, a marginal of exactly 1.0 in every
pixel, an 8.0 GB store (2.5 GB summed, 5.6 GB grid). The guard neither warned
nor refused, which is what a projection of 6.5 GB against 90 GB free should
do.

The constant held: peak private bytes were 6.38 GB, which over 21,101,151
features is **302 bytes per feature**, 8 percent under the 330 the guard
projects with and measured the way that number was, as the whole process's
private peak over the feature count. The `var` frame's own step is narrower
-- private sat at 2.49 GB through both passes and both column sorts and rose
to 6.38 GB over three seconds while the frame was built, so 184 bytes per
feature is what the frame alone costs -- and the guard is conservative under
either reading, on a grid 1.6 times the one the constant was measured on. It
stays at 330.

**Known limits.** The `var` frame is still the peak (D8), and the grid at the
default axis is 5.6 GB of store for a slide whose summed table is 2.5 GB.

---

## D5. One raw read per frame per pass serves every table

**Status:** Implemented 2026-09-07, the same day it was deferred. It was
deferred on a 400-frame warm measurement that put the grid's extra read at
about 4 percent of wall time; the whole-slide logs from the D1 measurement
then showed the mobility heatmap's own pass at 201 s of a 528 s conversion
under `scan_sum` and 272 s of 438 s under the vendor centroid. The extra
reads were 38 to 62 percent of a default conversion, not 4, and the reopen
condition was met on the same afternoon. The mapping this left as the next
lever was taken the same week, by the frame's own unique indices; see below.

**Decision.** A Bruker TDF reader hands each frame over once, as a record
of its raw scan read, from which the summed spectrum, the mobility point
cloud and the fragment spectrum of each precursor are all derived. On the
streaming route the two passes the summed table already takes, count then
scatter, feed the heatmap, the mobility grid and the MS/MS table from that
same read. Two raw reads of the source per conversion, whatever is written;
a default conversion used to take three, a grid four, a grid with the MS/MS
split six. The engine's two passes are still inherent, so "one pass" was
never the right name; "one read per pass" is.

**What makes it safe.** Every derivation goes through the very helpers the
reader's three iterators use, and every sink is the same accumulator the
standalone passes feed, so the fused route cannot see different numbers.
That was checked rather than assumed: stores written by the committed code
and by the fused code were compared table by table on five configurations
(the 400-frame slide with and without the grid on both write routes, and
the 713-pixel PASEF set with its split on both routes) and were identical
in every array, every `var` and `obs` frame and every `uns` block, index
dtypes included. A stub source that answers both the iterators and the
records pins the same identity in the unit tests, and pins that the fused
route reads it exactly twice and never through the iterators.

**Measured** (2026-09-07, warm, optical image left out):

| conversion | before | after |
|---|---|---|
| 400 frames, default (heatmap) | 13 s | 11 s |
| 400 frames, `--mobility-grid` | 32 s | 21 s |
| 713-pixel PASEF, MS/MS split | 12 s | 9 s |
| whole 26,087-pixel slide, default | 528 s | 454 s |
| whole slide, `--mobility-grid --resample-bins 40000` | 1,400 s | 858 s |

The grid row's "before" was measured on an earlier commit, before the
engine's column sort was made six times faster, so part of that gain is
the sort's; the 400-frame row above is the clean comparison for the grid.
On the whole slide the heatmap's own pass (201 s) became part of the count
pass, which grew from about 120 s to 266 s: what was saved is the read and
the index-to-m/z conversion of every frame, about 55 s, plus the second
scatter-side read. What remains is the mapping of every raw point onto the
mass axis, 51,000 points per frame on this slide, which the heatmap and
the grid need and the summed spectrum does not. That mapping was the next
lever, not another read, and it was taken the same week.

**The mapping, by the frame's unique indices** (implemented 2026-09-07). A
TDF frame is read as digitizer indices, and the reader already converts the
unique ones: measured on this slide, 28,560 unique of 49,070 points, 58
percent. Both halves of the mapping -- nearest bin, and in range or not --
are elementwise in the m/z, so mapping the unique values and gathering by
the frame's own inverse gives the same bins, the same mask and the same drop
count for 58 percent of the binary searches. The record offers that view as
`mobility_points_indexed`, the fused passes prefer it, and a record that
does not offer it is mapped point by point as before; the iterators and the
standalone passes are untouched.

**Measured** (2026-09-07, the same slide, warm, local disk, no optical
image; each pair one session, v3.19.0 against the branch):

| whole 26,087-pixel slide | v3.19.0 | indexed mapping |
|---|---|---|
| default, wall | 420 s | 343 s |
| default, pass 1 (maps every point, for the heatmap) | 250 s | 171 s |
| default, pass 2 (no mobility sink) | 137 s | 139 s |
| `--mobility-grid` at the default axis, wall | 1,201 s | 724 s |
| grid, pass 1 (heatmap and discovery) | 401 s | 237 s |
| grid, pass 2 (grid scatter) | 627 s | 357 s |

The pass that maps points is a third faster, and where both passes map, the
conversion is 40 percent faster; the default's second pass, which has no
mobility sink to map for, is unchanged, which is the control. Taking 42
percent of the searches out of the default's mapping pass took 79 s off it,
which puts the searches at roughly 190 s of that pass's 250 and roughly 110 s
of its 171 now -- still the largest single item in it, and the gather and
the mask that remain are proportional to the points however few searches
they cost.

Identity was checked the way the fused passes were: stores written by
v3.19.0 and by the branch, compared array by array (X data, indices and
indptr with dtypes, `var`, `obs` and `uns`), on the 400-frame slice with and
without the grid and on the 713-pixel PASEF set, and then on the whole slide
in both configurations by hashing every array of every table in chunks.
Identical throughout, the heatmap's `counts` included.

**Cold cache on the lab share** (2026-09-07, the deferral D5 left open).
The same slide copied to the SMB share, default conversion, branch code, the
file pushed out of the machine's standby list before the first run so the
read is genuinely cold: 365 s cold against 342 s warm, pass 1 178 s against
170 s, pass 2 143 s against 138 s. The first read of a never-read file costs
about 5 percent of the pass; a conversion's second pass is warm either way,
since the first pulls the acquisition into the cache. The share is not the
bottleneck the deferral suspected -- on this machine a cold network read is
within seconds of a warm local one (171 s for the same pass), so the mapping
and not the read is what a conversion of this shape spends its time on,
wherever the file lives.

**Objections considered.**

- *The in-memory route still runs the standalone passes.* True at the
  time, and left so: that route held the whole matrix in RAM and was taken
  by small files, where the passes cost seconds. Moot since D11 folded it
  into the streaming route, which feeds the sinks from its own passes.
- *Under `vendor_centroid` the summed spectrum cannot be derived from the
  raw scans.* Correct; in that mode the record asks the library for the
  centroid as a second call per frame, and the raw read still serves the
  sinks. That mode is the opt-in.
- *The record contract grew a method for one reader's convenience.* The
  indexed view is reached with `getattr` and a record that does not offer
  it is mapped point by point, which the stub records in the tests pin. It
  is a view of the read every record already has, not a new obligation.
- *Mapping every unique value maps ones that are then dropped as out of
  range.* It does, and it costs a few searches more on a frame whose points
  overhang the axis -- 8,328 points of 1.65 billion on this slide. The
  alternative, masking before mapping, would need the mask expanded to the
  points first, which is the work the factoring exists to avoid.

**Known limits.** Only the Bruker TDF reader hands frames over as records;
every other source keeps its iterators and its standalone passes, which for
an imzML mobility export is one pass over an already sparse file.

---

## D6. The MS/MS fragment axis is the MS1 mass axis

**Status:** Implemented. A condition attached to it on 2026-09-07, to refuse
the table on an unresampled axis, was withdrawn the same day; see below.

**Reason.** Fragments and precursors pass through the same TOF and have the
same resolving power, so one mass axis per store is the accurate
representation. A separate fragment axis would invent a second resolution
for the same analyser. The coupling is also what makes `mz_index`
meaningful and the conservation check exact.

**Objections considered.**

- *On a raw, unresampled axis fragments snap to MS1-only m/z values, and
  any fragment outside the axis range is dropped.* This objection was
  raised, accepted, implemented as a refusal, and then found to be wrong
  by the existing PASEF test, which relies on the raw axis for an exact
  equality. The premise fails because on a scheduled MS/MS acquisition the
  summed table holds no intact ions: its spectra are the same fragment
  peaks the split re-reads, so the raw axis is the union of the fragment
  m/z values themselves and the mapping is exact. The refusal was removed
  and a test now pins the raw axis as exact on both write routes.

**Known limits.** The resampling grid chosen for the summed table sets
fragment resolution. That is a consequence of the decision, not an
accident.

---

## D7. Waters: which functions hold the image

**Status:** Implemented 2026-09-07, then **restated the same day against
real files, which contradicted its premise.** The first version is kept
below because the correction is the point.

### What the first version decided, and why it was wrong

The Waters reader classifies acquisition functions and then yields every
scan of every MS-classified function as a pixel spectrum, keyed by laser
coordinate. A two-function acquisition -- MSe low and high energy, or a
data-dependent run -- would therefore emit two spectra at the same
coordinate, one MS1 and one MS2, and the per-scan MS level MassLynx reports
was parsed and then ignored. Summing an intact-ion and a fragment spectrum
into one pixel makes a spectrum of nothing, so the fix was to convert the
functions MassLynx labelled MS level 1 and record the rest under
`excluded_functions`. No Waters MS/MS imaging file was available; the
objection *no test data* was answered with "the filter is on a field already
parsed and is testable with a synthetic scan record".

That answer was wrong, and the entry said so itself without noticing: a
synthetic record can only confirm that the code reads the field it was
written to read. It cannot say what the field **means**.

### What the real files say

7,486 Waters `.raw` directories were found on the lab share (its six Waters
instrument folders, plus one user directory);
1,396 hold more than one `_FUNC*.DAT`, and the multi-function ones were
opened. **Every multi-function MALDI imaging run is a single-function raster
that MassLynx split across functions**, because it caps a `_FUNC*.DAT` file
at about 1.6 GB and opens a new *function* when a long run reaches it. In
each of those files:

- `_extern.inf` declares exactly **one** acquisition function -- "MALDI TOF
  MS FUNCTION", or "MALDI MOBILITY TOF MS FUNCTION" on the G2-Si -- however
  many functions the file holds. `_FUNCTNS.INF` carries one 416-byte record
  per stored chunk (10,400 bytes for the 25-function file).
- The chunks **tile** the stage and the run: consecutive, non-overlapping y
  bands and retention-time ranges, and **zero** shared pixels between any
  two functions. Their positioned scans sum exactly to the file's distinct
  laser positions (7,682 on `180814_EVO_Fresh_image.raw`, which is also
  exactly its 167 x 46 grid).
- MassLynx reports MS level 1 for the first chunk, 2 for the middle ones and
  **0** for the last, and `getLockmassFunction` names that last chunk as the
  file's lockmass function (it returns -1 only when the file has one
  function). `isMsFunction` is 1 for every chunk, including the one it calls
  lockmass, and **no chunk carries a precursor m/z**.

So on these files the level filter dropped every chunk after the first, and
the lockmass classification -- which predates this decision -- dropped the
last one on top of that. Pixels converted, against the pixels the file
holds:

| File | Instrument | Functions | v3.19.0 | Now (centroid) | In the file |
|---|---|---|---|---|---|
| `20140509_ZF_No13.raw` | Synapt G1 | 25 | 1,657 | 36,633 | 37,033 |
| `20140513_ZF_No16.raw` | Synapt G1 | 11 | 2,496 | 20,438 | 21,788 |
| `20141003 MTB_04.raw` | Synapt G1 | 8 | 1,700 | 12,110 | 12,972 |
| `140903_trypsin_vs_no.raw` | Synapt G1 | 5 | 2,284 | 8,843 | 9,221 |
| `180814_EVO_Fresh_image.raw` | Synapt G1 | 3 | 3,200 | 6,408 | 7,682 |
| `20201209_BCtumor_left Analyte 2.raw` | Synapt G2-Si | 4 | 15,790 | 49,096 | 61,108 |
| `20201223_BCTumor_Eva MDA_468.raw` | Synapt G2-Si | 3 | 5,837 | 11,413 | 17,550 |
| `20201218_MDA468_slide20201208 Analyte 5.raw` | Synapt G2-Si | 2 | 7,355 | 7,355 | 8,547 |

The two-function G2-Si run is the mildest case and still lost 14 percent of
the image; the 25-function one kept 4.5 percent of it. Nothing in the store
said so: the reader logged a warning about "MS level 2" functions and the
conversion looked healthy. The middle column is what the rule below converts
by default; the gap that remains is one chunk per file, and the next section
is why.

**One real multi-function acquisition was found**, and it is what makes the
rule below decidable rather than a guess:
the dataset registered as `waters_dda_neg_16func`, a fast-DDA run on a Xevo
DESI. Its
`_extern.inf` declares **16** functions -- one "TOF FAST DDA FUNCTION" and
15 "TOF SURVEY FUNCTION"s -- against the one function the chunked files
declare. Function 0 is MS1 with no precursor; functions 1 to 15 are level 2
and report a *different* precursor per scan, 28 distinct values each. And
every one of the 16 lands on the **same** position. So the two cases
separate cleanly on the positions: a chunked raster tiles them, a real
parallel acquisition repeats them.

(That file grids to 1x1 -- one unique x, one unique y -- and the reader now
refuses it rather than converting the whole acquisition onto a single pixel
with a 0.0 um pitch (issue #213). The refusal is right for this file: it is a
single-spot DDA run, not an image. It is *not* a DESI gap, as that issue
first supposed. Surveying the whole share settled it: real DESI images do
carry stage coordinates in the same laser fields and grid normally --
401x401, 247x140, 200x63 -- while 76 of the 105 `.raw` dirs there collapse
to one pixel and every one of those is a calibration, a tuning run, a lysis
test or a single-spot acquisition. Not one is an image.)

### Decision

Decide from the **laser positions**, the same measurement the pixel grid is
already built from.

- A function landing on pixels no earlier function covers **extends the
  raster** and is converted, whatever level MassLynx reports for it and
  whether or not MassLynx calls it the lockmass function. Only ever widens a
  file that already has an MS-classified function, so a run with no MS
  function is still refused.
- Functions **competing for the same pixels** were acquired in parallel.
  Only one of them can be the pixel's spectrum, so among those the MS1 ones
  win when the group has any, and the rest are listed with their level,
  precursor m/z, scan count and the reason they stayed out under
  `excluded_functions`. That is the original decision, kept -- scoped to the
  functions it was actually about.

`format_specific.function_types` keeps MassLynx's own classification next to
`format_specific.ms_functions`, so a rescued chunk is visible in the store.

### The catch: the tail chunk cannot be centroided

The library will not centroid the function `getLockmassFunction` names.
Measured on `180814_EVO_Fresh_image.raw`, sampling scans from each chunk:

| Function | `setCentroid(1)` | `setCentroid(0)` | `isRawSpectrumContinuum` |
|---|---|---|---|
| 0 (MS) | 7,408 points, TIC 3.35e4 | 82,876 points | continuum |
| 1 (MS, "level 2") | 8,508 points, TIC 3.60e4 | 91,274 points | continuum |
| 2 (named lockmass) | **112,594 points, TIC 1.11e5** | 112,594 points | continuum |

All three are acquired as continuum, and the request is honoured for the
first two and ignored for the third: it returns the same profile trace
either way. `ScanInfo.isProfile` reports this faithfully (0, 0, 1), since it
is read after `setCentroid`. There is no per-function centroid entry point
in the library to work around it.

Converting that chunk into a store of centroids therefore lays a band of
profile rows across the top of the image. It was converted that way once, by
accident, and the TIC image shows it: rows 0 to 37 average 3.1e4 to 3.7e4 and
rows 39 to 45 -- exactly the rescued chunk -- average 8.6e4 to 9.3e4, a
sharp 2.3x step at the chunk boundary and not a feature of the sample.

So the rule takes the representation into account: while the run is read as
centroids, a chunk the library will not centroid **stays out**, and
`excluded_functions` records it with its scan count, the pixels it would
have added and the reason. Reading the run as the profile trace
(`--waters-spectrum profile`, and the default on an MRT) selects every
chunk, because then all of them come back the same way. The warning names
the cost and the flags:

    Function(s) 2 hold 1274 pixels (16.6% of the image) that no other
    function covers, but MassLynx names them the lockmass function and will
    not centroid them. They stay out rather than put profile rows in a table
    of centroids: pass --waters-spectrum profile to convert the whole image.

`--streaming true` was in that sentence until D11 because, when it was written,
`--streaming auto` did not notice how large the profile store is: its
estimate assumed 10,000 peaks per spectrum whatever the source, which is
0.57 GB for this run, while the streaming converter's own estimate once
running is **74.1 GB** (7,682 pixels x 2,590,447 bins). Left on `auto` the
conversion stayed in memory and died at 72 percent asking for a 24.5 GiB
array on a 128 GB machine.

That was a general defect, not a Waters one, and it is fixed (issue #214).
Two things were wrong with the old estimate and both had to go, which the
measurements below settled:

| what `auto` scores this run at | GB |
|---|---|
| the old fixed 10,000 peaks per spectrum, 8 bytes each | 0.57 |
| the source's own measured width (85,117 points per spectrum), 16 bytes each | 9.7 |
| the same, against the axis the run is resampled onto (2,590,447 bins) | 296.5 |

The first fix is to stop guessing the width. `total_peaks` divided by the
spectrum count is *measured* by every extractor, and it is the profile bin
count on profile data and the peak count on centroided data. The second is
to count the 16 bytes a value costs the standard converter's COO arrays
(`int32` row, `int32` column, `float64` value) rather than 8.

Those two alone give 9.7 GB, which is still under the 10 GB threshold --
so they do not fix this run, and that is why the axis is consulted as well.
Interpolating a contiguous trace onto an axis 30x finer than it (2,590,447
bins against 85,117 points) fills the bins in between, and the array the
conversion died on says by how much: 24.5 GiB of `float64` is 3.3e9 values,
about 428,000 per spectrum, five times the source width. So a *profile*
source is sized at the axis it will be written onto and a centroid source
at its own peak count, since peaks stay peaks and a finer axis does not
multiply them. The bin count comes from the converter's own planner rather
than a second copy of that arithmetic, which is what issue #87 was.

This run then scored 296.5 GB and upgraded itself, so `--streaming true`
became belt-and-braces rather than required. That over-stated what the
conversion really needed -- those 3.3e9 values are about 49 GB across the
three COO arrays -- and deliberately so: the decision was one-sided, since
under-estimating kept a conversion in memory that had to stream and it
died there, while over-estimating cost at most a streaming run that would
also have fit. The estimate lasted one day: D11 removed the in-memory
converter, and with it the gate, so every conversion streams and the flag
selects nothing.

The alternative -- forcing the whole run to the profile trace whenever a
chunk cannot be centroided -- would keep the image whole automatically, but
it silently overrides D1's measured choice for a whole vendor and turns a
241 MB store into a 74 GB one, on every chunked file, of which this share
holds 1,396. Better to convert what is consistent, say what is missing, and
leave the trade to the flags that already exist.

The `fragmentation` block follows the **precursor**, not the level: a
converted function holds fragment spectra when MassLynx reports a precursor
m/z for it. A reported level with no precursor behind it is the chunk
artefact above, and reads as MS1.

**Verified on real data** (2026-09-07):

| Check | Result |
|---|---|
| `180814_EVO_Fresh_image.raw` converted, default centroid | 6,408 pixels against v3.19.0's 3,200, no duplicate coordinates, no empty pixel, `ms_functions [0, 1]`, `excluded_functions["2"]` carrying `n_unique_pixels 1274` and its reason, `fragmentation` MS1. The TIC image is continuous across the chunk boundary at row 19 |
| the same run with `--waters-spectrum profile --streaming true --resample-bins 60000` | all three chunks, **7,682 pixels** -- the whole 167 x 46 raster, no duplicate coordinate, no empty pixel -- and the TIC image is continuous across both chunk boundaries, so the band above was the representation and nothing else. 109M non-zeros, 870 MB. The axis was coarsened only to keep the store off a full disk; the default axis is the 74.1 GB estimate above |
| `20170818_08.raw`, a real single-precursor MS/MS run | `fragmentation` MS level 2, one window at m/z 377.4, collision energy 35.0, CID; the converter declines a demultiplexed table because one precursor needs none |
| `waters_dda_neg_16func`, the real DDA run | functions 1 to 15 recorded under `excluded_functions` with their per-scan precursors, function 0 converted, `fragmentation` MS1 |
| every file in the table above | no two converted functions share a pixel, and converted plus recorded pixels equal the file's distinct laser positions |

### Objections considered

- *The level is what MassLynx says; the chunking is its bug to report, not
  ours to work around.* It is not a bug that can be worked around later:
  there is no other field that separates a chunk from a high-energy
  function, and the level is wrong in both directions at once (0 for a
  chunk that holds image data, 2 for one that holds MS1 data). The
  positions are a direct measurement of the thing the rule is about --
  whether two spectra land on the same pixel -- so they are the better
  signal even if MassLynx were fixed tomorrow.
- *A lockmass function could legitimately carry laser positions, and would
  now be converted.* A reference function is acquired **alongside** the
  image, so its scans land on pixels the MS function already covers and it
  is excluded by the same rule. The one real parallel lockmass function
  found (`050517_BILE ACIDS 01.raw`) sits at a single constant position for
  the whole run, which no raster chunk does.
- *Two functions covering complementary mass ranges at one pixel should be
  summed, not filtered.* They still are: functions competing for a pixel and
  reporting the same level are all kept and summed, exactly as before. Only
  a mixed-level group is filtered.

**Known limit.** A raster chunk is recognised by covering new pixels, so a
chunked acquisition whose stage revisits a position -- a re-scan of the same
area in a later function -- would have that function read as a parallel one
and excluded. No such file was found, and the acquisition would be
ambiguous anyway: the store holds one spectrum per pixel.

### No Waters demultiplexed table yet

The second half of the original entry -- "a demultiplexed Waters table waits
until a scheduled Waters acquisition exists to look at" -- still waits. Of
the 506 methods whose `_extern.inf` was read, **113 declare a "MALDI TOF
MSMS FUNCTION"**, and every one of them has exactly **one** function isolating
**one** precursor, 224 KB to 3.9 MB, 9 to 60 scans on five or fewer distinct
positions: spot acquisitions, not images. A Waters demultiplexer needs a
file whose MS/MS functions each isolate a *different* constant precursor
over one raster, and no such file exists here, so it was not built.

What those files do establish, which the synthetic tests could not: MassLynx
reports the precursor faithfully when there is one (377.4, 482.0, 476.16,
each matching `Set Mass` in `_extern.inf`) together with a real collision
energy (35.0, and a 6 -> 30 eV ramp on one), so `precursor_mz` is the field
to trust. `quadIsolationStart`/`End` stay 0.0 even on these, so the
isolation window has no offsets from this API and only a target.

---

## D8. The `var` frame is the memory peak, and stays so

**Status:** Accepted limit.

With both sibling tables out of core, memory is proportional to the number
of features, not to the number of non-zeros. That is the intended shape;
D4 is the guard on it.

---

## D9. One streaming route

**Status:** Implemented 2026-09-08. The COO route went with PR #220 on the
same day this was accepted; `thyra/converters/spatialdata/` has held one
route since, and the code refers to the old one in the past tense.

The streaming converter had two write routes. PCS (pre-calculated scatter)
counts entries per column in one pass and scatters straight into
memory-mapped CSC arrays in the second, so the matrix is never a scipy
object in RAM. COO counted non-zeros per row, wrote CSR components to a
temporary Zarr, read them back whole into a `scipy.sparse.csr_matrix` and
called `.tocsc()` on it. A size threshold sent anything estimated under
30 GB to COO on the assumption that PCS bought memory safety at a cost in
speed.

Measured on the mock reader, peak process RSS sampled at 20 ms, one
subprocess per route, that trade did not exist:

| nnz | PCS | COO |
|---|---|---|
| 16M | 7.4 s / 540 MB | 10.4 s / 655 MB |
| 64M | 35.6 s / 1.1 GB | 58.5 s / 1.9 GB |

Both routes read the source twice, so the whole gap was the
materialise-then-convert step, and it widened with the dataset. v3.19
made PCS the unconditional default and left COO as an opt-in escape hatch.
Nothing reachable from `convert()` or the CLI could select it, Ousia's
wizard pinned `use_csc=True`, and its own bug class (the 79.5 GiB temp
directory leak, the missing out-of-grid guard, the pandas 3 string
failure) was the only thing it still produced. It is gone: about 750 lines
of `streaming_converter.py`, the `chunk_size` and `temp_dir` constructor
arguments that served only it, and its half of nine test modules.

**What stayed, for a day.** `use_csc` stays as a keyword because Ousia
passes it; `True` and `"auto"` mean the one route, `False` raises and says
why. The in-memory 2D and 3D converters stayed too, and for reasons rather
than inertia: only they wrote depth (the PCS scatter had no z term and
refused `n_z > 1`), they made one pass where PCS resamples every spectrum
twice, and they went through spatialdata's own writer, which is what caught
the five layout drifts the hand-written PCS store had (phantom rows,
missing root attrs, missing obs column, invented provenance, missing
encoding attrs). Deleting the route that used the real writer would have
left the hand-written layout with no oracle. Each of those three reasons
was a piece of work rather than a permanent argument; D11 did the work and
removed them.

**One thing this settled by accident.** `sparse_format="csr"` was only
honoured by COO; PCS ignored it and stored CSC. Since COO was unreachable
from `convert()`, a streaming CSR request had been silently producing CSC
for a release. The streaming converter now refuses `csr` and names the
in-memory converter as the one that writes it.

**Not decided here.** Whether PCS should grow a z term and route its table
through spatialdata's writer, at which point the in-memory converters and
the parity tests both become removable, and whether `--sparse-format`
earns its place as a flag at all. Both were filed as issues and both are
decided below: D10 (#219) and D11 (#218).

---

## D10. CSC is the only layout, and `sparse_format` is gone

**Status:** Implemented (2026-09-08).

**Decision.** Every converter writes CSC. `--sparse-format` and the
`sparse_format` keyword on `convert_msi` and the converters are removed,
and passing the keyword is answered rather than ignored. A caller who
wants row-major access calls `.tocsr()` on the matrix they read back.

**How it is answered, exactly.** The converter class raises
`ConversionRefused`. `convert_msi` does not re-raise it: its
`except ConversionRefused` logs the message once at `ERROR` and returns
`False`, which is what it does for *every* refusal since issue #234, and
what the CLI turns into exit 1. This paragraph used to say the keyword
"raises", full stop, which was true of the class and not of the front
door most callers use -- the discrepancy is issue #261's first item.
Both behaviours are pinned by tests
(`tests/unit/test_api_argument_gaps.py`,
`tests/unit/converters/test_spatialdata_converter.py`). What the decision
requires is that the removal is not *silent*; whether the caller learns
it from an exception or from `False` plus a logged reason is
`convert_msi`'s convention, not this decision's.

**Why.** After D9 the flag was honoured by the in-memory converters only.
That made it a flag whose effect depended on a second flag: `csr` did
something with `--streaming false` and was refused with it on. One option
per concept, and no option whose meaning changes with another, is the rule
the flag review of 2026-09-07 set; this was the clearest violation of it
left in the CLI. The layout it selected is also not one anybody here asked
for. Every consumer of a Thyra store reads columns -- an ion image is one
m/z across all pixels, which is one contiguous CSC column -- and Ousia,
the only shipped consumer, never passed the keyword at all.

**The alternative, and why it lost.** The honest version of keeping the
flag is to teach the streaming route a row-wise scatter, mirroring the
column one, so `csr` means the same thing on both routes. That is a second
two-pass write path, with its own count-then-scatter arithmetic and its own
half of every parity test, carried for a layout with no requester. The
conversion it replaces costs one `.tocsr()` in memory on data the caller
has already read.

**Why the keyword is answered instead of being accepted as a no-op.**
Unknown keyword arguments fall through `**kwargs` into
`BaseMSIConverter.options` without a word. Dropping `sparse_format` from the signature and stopping
there would mean `sparse_format="csr"` silently producing CSC -- which is
exactly the failure D9 found and this decision is meant to end. So the
base converter names it: the message says the keyword is gone, that CSC is
what is written, and what to call for rows. `"csc"` is refused on the same
terms, because a keyword that is accepted and does nothing is the shape of
option being cleared out here.

**The CLI differs, and deliberately.** `--sparse-format` is deleted
outright rather than kept as a hidden no-op the way `--optimize-chunks`
was. The precedents point in the same direction once the difference is
named: `--optimize-chunks` never had an effect, so accepting it could not
mislead anyone, while `--tof-a`/`--tof-b`/`--bins-per-fwhm` did have
effects and were removed outright. `--sparse-format csr` had an effect, and
the docs told people to pass it together with `--streaming false`; those
command lines have to fail loudly rather than quietly write the opposite
layout. Click's unknown-option error is that failure, and it arrives
before any data is read.

**Known limit.** A script passing `--sparse-format csc`, which asked for
what it would have got anyway, also stops working and needs the argument
deleted. That is the cost of not having a value that is accepted and
ignored.

---

## D11. One converter

**Status:** Implemented (2026-09-08).

The in-memory 2D and 3D converters are folded into the streaming one.
There is one converter: two passes over the source, every table scattered
into memory-mapped CSC arrays (`csc_assembly.CscAssembly`, the engine the
sibling tables already used, with every m/z bin a column) and written
through spatialdata's writer as an AnnData over those memmaps. `handle_3d`
decides the shape of the store, one table per z plane or one for the
volume, exactly as it used to pick between the two in-memory converters.
`SpatialData2DConverter` and `SpatialData3DConverter` are gone;
`SpatialDataConverter` is the registered name of the one converter and
`StreamingSpatialDataConverter` the class, kept because Ousia imports it.

**The measurement D9 asked for.** Whether small datasets convert faster in
memory had never been measured. On the real files in `test_data/` and the
registry, warm, `--no-optical`, one subprocess per conversion, peak RSS
summed over the process tree and sampled at 50 ms (2026-09-08, before the
fold, so both routes are the shipped code):

| dataset | spectra | one pass, in memory | two passes, streaming |
|---|---|---|---|
| `pea.imzML` | 12,737 | 12.2 s / 4.0 GB | 15.7 s / 1.0 GB |
| `bellini.imzML` (36M nnz) | | 8.3 s / 1.1 GB | 10.8 s / 0.6 GB |
| a TSF slide | 33,800 | 19.2 s / 5.1 GB | 28.2 s / 1.3 GB |
| `tims_msms_pos_brain1` (TDF, PASEF) | 713 | 10.7 s / 0.42 GB | 10.0 s / 0.40 GB |

So the claim was true, and small: the one-pass route is about 1.3x faster
on small imzML and TSF files, at three to four times the peak memory, and
on a TDF the fused passes (D5) had already made the two equal. Where the
time goes on `pea.imzML`: the in-memory route spends 5.3 s in its one pass,
0.8 s converting COO to CSC and 3.2 s writing; the streaming route 4.6 s
counting, 5.9 s scattering and 2.2 s writing. The second pass costs more
than the first because the scatter writes each entry to a random position
of the memmap, where the count only increments a dense array. That is the
lever if the second pass ever matters, and it lives inside the one route.

The same files through the one converter after the fold, same harness,
best of four runs: `pea.imzML` 15.8 s / 1.15 GB, `bellini.imzML` 10.7 s /
0.67 GB, the TSF slide 29.1 s / 1.48 GB, the PASEF TDF 8.4 s / 0.45 GB.
The wall time is the old streaming route's to within a second and the
PASEF set is faster; peak memory sits 10 to 15 percent above the old
streaming route, which is the writer's 128 MiB shard buffer and the
`var` frame the hand-written layout wrote as bare arrays. Store sizes are
the old in-memory route's to the megabyte, since it is the same writer.

**Why one route anyway.** Three seconds on a fifteen-second conversion
does not pay for a second code path with its own bug class (a COO
pre-allocation sized from `total_peaks` that grew by 50 percent when the
count was wrong, a 4 GB RSS on a 700 MB file) and its own half of every
converter test. The three reasons D9 kept the in-memory converters were
each a piece of work, and they are done:

- *Depth.* The scatter indexes rows within a table unit -- `y * n_x + x`
  on a plane's table, `z * n_x * n_y + y * n_x + x` on the volume's -- and
  a multi-plane source gets one table per plane or one volume, with the
  same element keys, `obs` columns, `instance_id` convention (the plane-
  local grid index per plane, the volume grid index for a volume), TIC
  images and shapes the in-memory converters wrote. The refusal of
  `n_z > 1` is gone; the measurement that justified it now asserts the
  planes are apart (`test_pcs_phantom_pixel_rows`).
- *The writer.* Nothing composes a Zarr layout by hand any more. Encoding
  attributes, root attributes, `obs`, `var` and the provenance block all
  come from `_save_output`, which every table goes through, so the five
  drifts D9 lists cannot recur -- there is nothing left to drift from.
  The seam the issue asked about is `CscAssembly.matrix()`: scipy's
  constructor takes the memmaps without copying when the index dtype is
  the one it would have chosen, and anndata's sparse writer copies each
  array into the store shard by shard, so the write stays out of core.
- *The oracle.* The real writer is not an oracle when it is the only
  writer; it is the layout. The value-level oracles stay and are the
  guard now: `test_stored_pixel_spectrum_oracle` (computed expectation
  per pixel, never read back), `test_out_of_grid_guard`, the 3D probe
  readers and the sibling-table stubs. The cross-path parity suite went,
  as the issue said it could; what survives of it checks that both table
  shapes describe the source the same way.

**What changed on disk.** Nothing in the values. Measured rather than
claimed: the four files above were converted with the code before the
fold down both routes and with the code after it, and every store was
read back through spatialdata and compared element by element -- root
attributes minus the timestamp, every table's `X` arrays, `obs`, `var`
and `uns`, every shapes element's index, geometry and transform, every
image's pixels and transform. Against the old in-memory route all four
new stores are identical. Against the old streaming route they differ in
layout details that were the in-memory route's all along -- `indptr`
follows scipy's index-dtype rule (int32 while it fits) instead of always
int64; `obs` string columns are anndata's nullable string arrays; the
summed table's shards come from `table_write_config()` (128 MiB) instead
of hand-rolled one-million-entry chunks; a raw-axis dense spectrum no
longer stores its explicit zeros -- and in two drifts of the hand-written
layout that nobody had listed: it stamped `spatialdata_software_version`
with a literal `0.6.1` from whenever it was written, where the writer
records the installed version, and it gave the sibling tables an `obs` of
`x`, `y`, `region` and `instance_key` only, where the in-memory route
gave them the summed table's full `obs` (`spatial_x`, `spatial_y`,
`region_number` included), which every sibling now gets. Two things
changed the other way. The grid table's `mobility_marginal` block
carries the per-pixel ratio only: the per-cell deviation the in-memory
route recorded needed the marginal and its difference from the summed
table materialised, each the size of the summed table, and the route that
was kept never holds either. And a single-plane FlexImaging acquisition
converted with `handle_3d=True` now places its TIC image by the alignment
affine like every other 2D image, where the old 3D converter scaled it in
micrometres while putting its shapes in optical pixels, so the two
disagreed at `"global"`.

**What went with them.** `--streaming` is a hidden no-op and `streaming=`
selects nothing (though it is still checked: only `True`, `False` and
`"auto"` are accepted, because an argument that selects nothing still has
to say when it was misspelled -- issue #261): with one route the `auto`
estimate has nothing to decide,
so the sizing PR #216 landed the day before (`_values_per_spectrum`, the
profile-versus-centroid rule, `Thresholds.STREAMING_SIZE_GB`) is gone with
the gate it served. The refusal of a broken file still happens at the
first metadata read, outside any try, which `_create_converter` does
itself now. `sparse_format` went the same day (D10), which this makes
moot twice over: the only writer of CSR was in-memory.

**Objections considered.**

- *Small files got slower.* By the numbers above, and said here rather
  than hidden. Conversions of a few seconds are not what the memory
  budget is for; conversions of hours are, and those already streamed.
- *A second pass reads the source twice.* It always did on the route
  that was kept, and on a cold network share a pass is bounded by the
  mapping, not the read (D5's cold-cache measurement).
- *The in-memory route dropped explicit zeros and the streaming one kept
  them.* It drops them now, on the raw-axis path where they arose;
  the resampling paths never produced any.

**Known limit.** The `var` frame is still O(bins) in RAM (D8), and the
string index it carries is built by pandas at about 125 bytes per bin: a
raw, unresampled axis of millions of bins is the one thing that still
scales with the source, and it is the case to resample.

---

## D12. An anisotropic raster is carried, not refused

**Status:** Implemented (2026-09-09), issue #228.

**Decision.** When a source declares a different pitch on each axis, the
store carries both — root attrs, `coordinate_systems.global` and its
`raster_to_global_affine`, the image and shapes transforms, the pixel
footprints, `obs["spatial_x"]`/`["spatial_y"]` and the `msi_metadata`
block. Thyra does not refuse the acquisition, and `--pixel-size` keeps its
meaning of one number applied to both axes.

**The defect this replaces.** The converter carried a single float. Pixel
size detection kept the source's x pitch (`# Use X size`), the root attrs
wrote that same value as `pixel_size_y_um`, and it went on into the affine
and the footprints; only `msi_metadata.ms_analysis.pixel_size_um` read the
true pair back out of the detection info. A raster acquired at 30 x 50 um
was stored as 30 x 30 in one block and (30, 50) in another. So the store
contradicted itself and rendered squashed by y/x.

**Why not refuse.** Refusing was the alternative the issue named, and it
loses on both counts that decide it:

- *The format is already per-axis.* `pixel_size_um` is `{x, y}` and
  **required** by the metadata schema, mapped to `IMS:1000046` and
  `IMS:1000047`. Nothing needed adding and no schema version moves; the
  converter was the only thing carrying one number where the format
  carries two. Refusing would have been declining to write a store the
  format already describes.
- *The only escape would be a false statement.* `--pixel-size` overrides
  detection with one number on both axes, so the way out of a refusal
  would be to declare a 30 x 50 um raster square. That is a worse store
  than the one the refusal was protecting against, because nothing in it
  records that a number was invented.

**What it cost.** `self.pixel_size_um` had 31 call sites across
`base_spatialdata_converter.py`, `streaming_converter.py` and
`core/base_converter.py`, and each had to be read to decide which axis it
meant. Most are x, or a single number that is written twice. Eleven meant
y and now take `pixel_size_y_um`: the root attrs, `pixel_size_um_y` and
the `raster_to_global_affine`, `stage_offset_um`, the optical-to-micrometre
scale matrix, the half-pixel footprint, `obs["spatial_y"]` on all three obs
builders, and the TIC image's `Scale` on the plane and volume paths. Two
provenance echoes gained a `pixel_size_y_um` key beside the x one. Two were
judgement calls:

- *`_resolve_z_spacing`'s fallback* ("nobody said, so reuse the in-plane
  pitch") takes x. There is no single in-plane pitch on an anisotropic
  raster and the choice is arbitrary — which is the point, since nothing
  about the raster predicts the section thickness either way, and
  `z_spacing_source` already marks the number as assumed rather than
  measured.
- *`stage_offset_um`* multiplies the source's raw acquisition-index
  offsets by the pitch, so each offset takes its own axis.

**The tolerance.** The anisotropy is reported in the log only when the two
pitches differ by more than float noise (`rel_tol=1e-9`). Nothing is
refused, so the tolerance decides a log line rather than whether a file
converts. It is deliberately tight: a real anisotropy is a percent or
more, and the near-misses a detector used to produce were not noise but a
bug — until issue #217 the Waters extractor returned 29.66 x 28.93 for
every square raster, which this code would have carried faithfully into
the affine. That is fixed at the extractor, which is where it belonged.

**Known limit.** No dataset in the registry is anisotropic, so this is
tested synthetically: a reader is asked for a 30 x 50 um pitch and the
written store is read back. A DESI method with `DesiXStep != DesiYStep`
would be the first real case, and nothing here has been run against one.

---

## D13. A peak on a declared mass-range bound lands in the edge bin

**Status:** Implemented (2026-09-09), issue #239.

**Decision.** A peak is in range when it is inside the **declared**
`[min_mz, max_mz]`, not when it is inside `[axis[0], axis[-1]]`. Both
resampling methods use that range: nearest-neighbour keeps the peak and
maps it to its nearest bin, TIC-preserving measures the preserved share
against it. Peaks outside the declared range are still dropped, not
clamped.

**The defect this replaces.** Every physics generator lays
`target_bins + 1` bin *edges* across the requested range and returns the
midpoints, so the first and last axis point sit half a bin inside what was
asked for. A source declaring 50-1000 m/z built the axis
`[50.0001, 999.9975]`, and a peak sitting exactly on a declared bound was
below `axis[0]` and discarded. It costs most where a source declares its
range *as* its first and last sample: `phi_extractor` takes `mass_range`
from the first and last detector channel, so both were dropped in every
pixel — the mock fixture stored 12 counts against 14 in the source, where
`--no-resample` stored all 14.

**Why not widen the axis.** Widening so that
`axis[0] <= min_mz <= max_mz <= axis[-1]` was the other candidate the
issue named. It changes the axis itself, so `var["mz"]` and the bin count
move on **every** resampled store, and two stores of the same acquisition
written either side of the change no longer share a feature axis. Edge
bins change only which bin a boundary peak lands in — nothing else about
any store moves, and a store with no peak on a bound is byte-identical.

**Why it cannot become the clamp again.** `_nearest_neighbor_resample`
used to clip every out-of-range index into the axis and accumulate, so
narrowing the range piled the discarded part of the spectrum onto two
bins: on real `pea.imzML` resampled to 400-800 m/z, bin 0 held 654,158
counts where a real peak there is around 80, with the total conserved
exactly so no TIC check could see it. The new rule reaches at most half a
bin beyond the first and last centre, because the declared range *is* the
outer bin edges. A peak further out is still dropped.

**Where nothing changes.** A uniform axis is
`np.linspace(min_mz, max_mz, n)`, whose end points already are the
declared bounds, so `constant` stores are untouched. `--no-resample` is
untouched too: no axis is built, and the range falls back to the axis's
own span, which is what the rule was before.

**What it is worth, measured.** On `tims_msms_pos_brain1` (713 PASEF
frames, 302,107 peaks, declared 50-1000 m/z) **zero** peaks lie in the
half-bin skirt; the one peak the warning reports is genuinely above 1000
and is still dropped. So on data whose peaks are interior this is a
correctness fix that recovers nothing. The sources it pays on are the ones
whose declared range is a detector channel rather than an acquisition
setting.

---

## D14. A 0-based imzML folds its base down; a cropped one does not move

**Status:** Implemented (2026-09-09), issue #244.

**Decision.** x and y are rebased on `min(observed_minimum, 1)`. A file
whose smallest coordinate is 0 is 0-based and is rebased on 0; a file
starting at 1, or at 5, keeps the specification's base of 1. z keeps its
own rule — the smallest value present. What was subtracted is reported as
`EssentialMetadata.coordinate_offsets` and written to
`coordinate_systems.global.coordinate_offsets_px`.

**The defect this replaces.** Three sites subtracted a constant 1 (the
reader's cached coordinate array, its per-spectrum fallback, and the
extractor's grid sizing). On a file written 0-based that produced
`x = -1` and `y = -1` for the first row and column, which the converter's
`_locate` guard dropped. A 3x3 file at coordinates 0..2 previewed as
`grid (2, 2)`, warned that "5 spectra sat outside the declared 2x2x1
grid", stored 4 rows, and exited 0.

**Why not rebase on the observed minimum, as z does.** `_z_base` measures
its base off the file and its docstring argues the case well — imzML
guarantees nothing about the base and pyimzml is inconsistent about it.
The argument does not carry over, because **z has no physical origin and
x and y do.** An acquisition cropped to a region of the slide legitimately
starts at `x = 5`; rebasing on the observed minimum would slide it to
`x = 0`, changing the grid width, every `obs["spatial_x"]`, the TIC image
extent and the pixel footprint — on a file that converts correctly today.
Nothing in the file distinguishes that acquisition from a 0-based export
whose first column happens to be empty. Folding only a 0 down fixes the
reported defect and moves nothing else, which is the whole point: the only
files whose stored coordinates change are the ones that were losing a row
and a column.

**Why not the declared pixel counts.** `IMS:1000042` / `IMS:1000043` were
the third option. Nothing in Thyra reads them today except
`mzpeak_extractor`, which carries a comment about a declared extent
disagreeing with the coordinates it ships with. Trusting a declared extent
over the coordinates is a larger change with a wider blast radius than the
defect it would fix.

**What the warning says now.** Nothing, on this path: a 0-based file is no
longer off-grid, so the out-of-grid warning does not fire for it and keeps
its meaning for a reader whose coordinates genuinely disagree with the
grid it declares. The base itself is reported once, at INFO, naming the
smallest coordinate found.

**Known limit.** No 0-based imzML is in the registry, so this is tested
against a written fixture (`zero_based_imzml`), alongside a 1-based one
and a cropped 1-based one — the last being the file the rejected
alternative would have moved.

---

## D15. An override that contradicts the detector warns; it does not refuse or self-correct

**Status:** Implemented (2026-09-09), issue #246.

**Decision.** When `--resample-method` is given explicitly, the detector
is asked what it would have chosen and a WARNING is logged if the two
differ. Nothing stored changes. For a `tic_preserving` override the
warning names `--resample-gap-tolerance`, or reports the tolerance already
in force.

**The defect this replaces.** The detector has a verdict for every source
and only the `auto` path ever asked for it, so an explicit method was
applied with nothing checked and nothing said.
`--resample-method tic_preserving` on a Bruker TDF — for which the
detector chooses nearest-neighbour — interpolates across the gaps of a
sparse centroid list and fills the axis. Measured on
`tims_msms_pos_brain1`:

| | nearest_neighbor (default) | tic_preserving |
|---|---|---|
| stored non-zeros | 302,106 | **423,386,757** |
| table | 7.6 MB | 583 MB |
| peak RSS | 0.5 GB | 5.9 GB |
| wall | 26 s | 57 s |

Per-pixel TIC is identical either way, so the TIC identity cannot see it.
What breaks is the siblings, quietly: MS/MS blocks against the summed TIC
come to 0.998711 per pixel, and the heatmap marginal against the stored
mean spectrum to rel 68 with nothing recorded at all.

This is the same bug class as #168 on PHI ToF-SIMS, which was fixed *by*
adding a detector — and that is precisely why a detector is not enough on
its own. A detector only ever steers `auto`.

**Why not refuse.** A refusal would have to fire on "points per spectrum
is a small fraction of the axis", which is true of *any* centroid source
against a fine axis, including the ones where `tic_preserving` with a gap
tolerance is exactly what the user wants. The threshold would be arbitrary
and the refusal would break working commands.

**Why not an automatic gap tolerance.** It changes stored values on every
sparse `tic_preserving` conversion, including ones somebody is relying on,
and there is no principled value to choose: half the widest source gap is
per-spectrum and data-dependent, so the stored numbers would depend on a
heuristic no flag records. The remedy already has a flag, and one flag per
concept is the rule.

**That the remedy works, measured.** The same conversion with
`--resample-gap-tolerance 0.01`: **4,801,946** stored non-zeros, an 88x
reduction against the unguarded override.

**Where the size is already reported.** Batch 3's `_matrix_size_gb` logs
the counted non-zeros after pass 1 — 423,386,757 entries, before pass 2
writes anything. The axis guard (`AXIS_BYTES_PER_BIN`) does not fire and
should not: 599,146 columns is an ordinary axis. What exploded is the
matrix.

---

## D16. A preview that cannot count spectra says so rather than reporting the raster

**Status:** Implemented (2026-09-09), issue #240.

**Decision.** `MsiPreview.n_pixels` is `Optional[int]`, and is `None` when
the format cannot count spectra without decoding them.
`EssentialMetadata` gains `n_spectra_counted`, false only on that path,
where `n_spectra` and `total_peaks` are 0 meaning "not counted". PHI is
the only format that takes it, and only under `metadata_only=True`.

**The defect this replaces.** `preview_msi` promises "No spectra are
decoded" and passes `metadata_only=True` for that purpose.
`PhiReader.__init__` swallowed the kwarg through `**kwargs`, and
`PhiMetadataExtractor` called `occupied_channel_counts()`, which
aggregates every 8-byte event in the stream. Measured at 0.16 s for a
16 MB file with 2.02 M events and 0.14 s for a 14 MB one — linear, so a
multi-gigabyte SmartSoft acquisition previewed as slowly as it converted.

**Why not report the raster size.** From the header alone PHI knows the
tile geometry (`n_x * n_y`) but not which pixels carry events: it stores a
stream of ion arrivals, not a list of spectra. Every other reader reports
spectra *present*, and cheaply — imzML `len(coords)`, Bruker a SQL count —
so putting `n_x * n_y` in `n_pixels` would make that field mean "positions
the raster covers" for one format and "spectra present" for every other.
The raster extent is already reported, in `grid_dims`, so nothing is lost
by declining to answer twice.

**Why a flag rather than `Optional[int]` on `n_spectra`.** `n_spectra` is
read on the conversion path in half a dozen places where it is always a
real count; widening its type would push a `None` check into all of them
to describe a state none of them can reach. A default-true boolean beside
it says the same thing and is inert everywhere else.

**Precedent.** Bruker already does this for the other half:
`skip_total_peaks=self._metadata_only` reports `total_peaks` as 0 without
scanning the Frames table. `skip_event_aggregate` is the same idea, and
`n_spectra_counted` retro-fits the missing part — a way to tell
0-not-counted from 0-none-present.

**What the preview still answers.** Dimensions, coordinate bounds, m/z
range, pixel size and both detector verdicts, all from the header and the
block chain. `PhiToFSIMSDetector` matches on the format flag rather than
on peak density, so zeroing the counts does not cost the preview its
nearest-neighbour verdict — which would have been a regression of #168.

---

## D17. The docs deploy builds `main`, and the `.ibd` UUID warns rather than refuses

**Status:** Implemented (2026-09-09), issues #223 and #261.

Two independent calls, both of the same shape: a check that could be made
stricter, and the measurement that says how strict it should be.

### The docs deploy pins its ref

`docs.yml` ends every run in `mkdocs gh-deploy --force`, which force-pushes
the whole built site over `gh-pages`. Two runs overlapping is a lost-update
race. #154 removed the duration-dependent half of it with a `pages-deploy`
concurrency group, and deliberately left a narrower half: GitHub's
workflow-syntax documentation says runs in a group are processed by the time
each started waiting, and then adds, verbatim, that **"ordering is not
guaranteed"**. An older run admitted second builds its own older checkout and
force-pushes it over newer docs -- both runs green, nothing in either log.

**Decision.** The deploy checks out
`${{ github.event_name == 'push' && 'main' || github.ref }}`. Whichever run
executes last then builds current `main`, so the published site is right
regardless of admission order.

**Why this rather than tighter serialisation.** There is nothing tighter
available: the group is already unkeyed (one `gh-pages` branch, so a
`workflow_dispatch` and a push must contend), and `cancel-in-progress: true`
would reintroduce the race it was added to remove, because cancellation is not
synchronous and the older run's force-push can still be in flight. The
ordering guarantee simply is not offered. Not depending on it is the only fix
that does not depend on it.

**The cost, accepted.** Re-running an old docs run now publishes today's docs
rather than reproducing that run's build, and `gh-deploy`'s `Deployed <sha>`
message stops matching the run's trigger. For a workflow whose only job is
"deploy the site" that is the behaviour wanted; it would be wrong for a
workflow meant to reproduce a historical build, and this is not one.

**Why conditional rather than a hardcoded `main`.** A hardcoded ref makes a
deliberate `workflow_dispatch` against another branch silently deploy `main`
instead -- a dispatch that does the opposite of what it says. The conditional
keeps that path honest and changes nothing about the release path:
`release.yml`'s `deploy-docs` job dispatches with `--ref main` already.

**Known limit.** The race is narrowed to nothing for *content*, not for
*attribution*: if two runs overlap, the gh-pages commit message names
whichever ref the surviving run checked out, which is now always `main`'s tip
rather than the merge that triggered it.

### An `.ibd` UUID that disagrees warns

The imzML specification puts the binary file's UUID in the first 16 bytes of
the `.ibd` and the same value in the XML as `IMS:1000080`. Thyra read the XML
term for the metadata store and never compared the two (issue #261, item 5).
Comparing them is the only check that can tell an `.imzML` apart from a
*different* acquisition's `.ibd` sitting beside it under the right name: every
other check in `_validate_parser_state` reads the XML's own offsets and
lengths against the binary's size, which a wrong-but-plausible pairing
satisfies.

**Decision.** Compare them, and **warn** on a disagreement. Do not refuse.

**Why not refuse, measured on 2026-09-09.** Of the three real files in the
corpus, two match byte for byte:

| file | writer | declared `IMS:1000080` | first 16 bytes of `.ibd` |
|---|---|---|---|
| `pea` | SCiLS | `9069a51f-11a8-4f15-aabb-43624def10d0` | same |
| Xenium export | SCiLS | `6c176aa6-5cf8-4027-8661-657695f2d400` | same |
| `bellini` | IONTOF SurfaceLab 7.5 | `{FC37F303-A9C0-4CD3-A28E-1D18E523C269}` | `3ad1bacd-dcc3-4f7b-aea7-f9b375dbf731` |

`bellini`'s first spectrum starts at byte 16, so the header slot is there and
populated -- the writer simply put two different values in the two places.
That file passes every other check and converts correctly. A refusal would
therefore reject data Thyra reads right today, over a disagreement between two
copies of an identifier that is never used to locate a byte. One real file in
three is not a rare shape.

**What the warning has to say, then.** Both values and both filenames, and
that reading continues -- because the likeliest reading is "this vendor fills
the two fields independently", and the person needs to be able to dismiss it
without reading the source. The second reading, a genuinely mispaired `.ibd`,
is the one worth the noise.

**Comparison rule.** The 32 hex digits only. Vendors differ on the registry
braces (IONTOF writes them, SCiLS does not) and on case, and the hyphens are
positional rather than data, so a file re-spelling its own UUID in the other
convention must not read as a disagreement with itself.

**Silent when either side is absent.** A file that declares no UUID, or an
`.ibd` shorter than its header, has nothing to compare; the extent checks
already speak for the truncated case.

---

## D18. `main` stays unprotected until a release credential exists

**Status:** Deferred (2026-09-09), issue #224.

**Decision.** Do not turn on required status checks for `main`. Keep the
existing ruleset, which forbids deletion and non-fast-forward pushes and costs
nothing.

**Why.** Required status checks and this repository's release flow cannot both
work with the credentials available. Measured on 2026-08-01 against a
throwaway branch, and re-confirmed on 2026-09-09:

| configuration | result |
|---|---|
| no protection (baseline) | `github-actions[bot]` push **allowed** |
| ruleset, required checks, bypass = admin role | push **blocked** (`GH013`) |
| ruleset + GitHub Actions app (id 15368) as bypass actor | **refused by the API**, `HTTP 422` |
| classic protection, `enforce_admins: false` | push **blocked** (`GH006`) |

`release.yml` runs `semantic-release version`, which pushes the version commit
to `main` with `secrets.GITHUB_TOKEN`. That push cannot bypass required
checks, and the Actions app cannot be added as a repo-level bypass actor --
"must be part of the ruleset source or owner organization". Turning the rule
on would leave the repository tagged-but-not-bumped on the first releasable
merge. **A repository that cannot cut a release is worse than one without
required checks.**

**The obvious wrong assumption, named.** Moving the release trigger off `push`
(the batched-release change) did **not** fix this. It changes *when*
semantic-release pushes the version commit, not *what* that push is.

**What has changed since the measurement, and what it means.** The check list
has grown from the seven contexts the issue lists to twelve: four
`test (os, py)`, one `lint`, four `integration (os, py)`, two
`clean-venv-install (py)` and one `complexity-check`. Only `tests.yml` and
`complexity-monitoring.yml` run on a pull request at all -- `docs.yml` is
`push: [main]` plus `workflow_dispatch`, and `release.yml` has no PR trigger
-- so those five jobs are the whole list. Every addition makes the rule
stricter, not more reachable. More to the point, the risk the issue was
opened for is already
closed by a different mechanism: `release.yml` waits for the **Tests**
workflow on the exact commit it is about to release, and refuses to release
when that run is absent, incomplete, or not successful. Required checks would
add "a red PR cannot be merged"; they would not add "a red commit cannot reach
PyPI", because that is already true. This is why the issue is deferred rather
than escalated.

**What would reopen it.** Any one of: a repo-admin PAT or GitHub App token in
`release.yml`'s checkout `token:` (note Ousia went deliberately no-PAT for its
own release flow, so this is a reversal, not a detail); a deploy key as a
`DeployKey` bypass actor, pushing over SSH; or an organization-level ruleset,
where the 422's wording hints the Actions app may be an acceptable actor --
untested, because it needs `gh auth refresh -s admin:org`, which is
interactive.

---

## D19. Lint can veto a release

**Status:** Implemented (2026-09-14), issue #291.

**Decision.** The `lint` job -- `pre-commit run --all-files`, which is black,
isort, flake8, mypy, bandit, pydocstyle and the repo-local path guards -- is a
job inside `tests.yml`, the workflow named **Tests**. `release.yml` polls for a
completed, successful run of that workflow on the exact SHA it is about to
release and fails closed. So a red lint on `main` stops a publish, and that is
intended rather than incidental.

**Why it is a decision and not a placement.** `release.yml`'s own comment says
the gate is on tests *specifically*, and gives the reason: waiting on "every
check for this commit" would let an **advisory** workflow veto a release, and
`complexity-monitoring.yml` is exactly such a workflow. Putting lint inside
`tests.yml` hands lint the veto that comment was written to withhold from
complexity monitoring. The two are not the same kind of check. Complexity is a
number that a legitimate refactor can move; a lint failure is a file that does
not meet the rules this repository has written down, including the two guards
that keep colleague names and lab-share layouts out of a **public**
repository. Refusing to publish from it is the right answer.

**The alternative, and why not.** A separate `lint.yml` with the same `on:`
block would enforce lint on pull requests and leave the release gate untouched.
It was not taken: it buys a second workflow file, a second cache key and a
second place for the trigger block to drift, to preserve a distinction between
"blocks the merge" and "blocks the publish" that nobody wants for this check.

**The cost, accepted.** A failure with no bearing on the code can now block a
release: pre-commit clones each hook's repository, so a GitHub outage or a
deleted hook tag turns into a red **Tests** run, and that SHA stays unreleasable
until a new commit lands. The window is small -- releases are batched and cut
deliberately (`gh workflow run release.yml`), not on every merge -- and a
blocked release is not a lost one, since the commits stay on `main` and the
next run releases them together.

**Known limit, and what would reverse this.** If a hook-repo failure ever
actually blocks a release, the answer is to split `lint.yml` out with the same
`on:` block. It is **not** to weaken `release.yml`'s gate: that gate is the only
thing standing between a red commit and PyPI, and it fails closed on purpose.

---

## D20. The PHI target axis follows measured peak width, fitted at its narrow edge

**Status:** Implemented (2026-09-14).

**Decision.** `PhiToFSIMSDetector` asks for `AxisType.TOF` with the measured
pair `A = 0.454`, `B = 0.0284`, laid at the usual three bins per peak width.
It keeps `nearest_neighbor`, and keeps reporting `linear_tof` as its
`source_grid_law`. The two are deliberately different laws.

**Why the source law is not the target law.** `PhiMassAxis` steps at a
constant flight time, so the detector's *channels* are spaced as `sqrt(m/z)`.
That is a fact about the digitiser and is what `source_grid_law` reports. The
*peaks* do not follow it: measured over 311 isolated peaks from twelve
acquisitions, resolving power climbs from about 2,200 at m/z 10 and levels off
near 4,250 above m/z 60. A rise then a plateau is neither single-term limit --
`linear_tof` has R climbing as `sqrt(m)` without end, `reflector_tof` has it
flat from the start -- and it is exactly the knee `sqrt(A m + B m^2)`
describes. Only peak width has any business setting a bin width.

**What the old answer cost, measured 2026-09-14** on a 512x512 negative-mode
nanoTOF acquisition against the instrument software's own `.bif6` peak-image
export:

| | `--no-resample` | old `linear_tof` default | this decision |
|---|---|---|---|
| bins | 863,670 | 86,204 | 103,070 |
| store | 65.8 MB | 46.5 MB | 51.6 MB |
| TIC vs the vendor export | exact, 262,144/262,144 px | exact | exact |
| per-peak recovery | 97.9-103.6% | 88.6-113.1% | 92.8-105.1% |
| worst correlation | 0.9509 | 0.9024 | 0.9695 |
| bins in the 6 mDa window at m/z 27 | 12 | 1 | 3 |

Inheriting the generic `linear_tof` width of 17 mDa at m/z 300 put 42% of the
corpus's peaks under two bins per peak width and one single bin inside the
window that acquisition's own peak list uses to separate C15N- from 13CN-,
6.3 mDa apart. **No total-ion check could see any of it**: nearest-neighbour
conserves counts, so the total ion image is bit-exact against the vendor's
export on all three axes. Only windowed numbers move, and on the old axis they
moved in both directions at once -- 113.1% for one ion and 88.6% for another --
which is the signature of one bin per peak.

**Why the fit is a 5th-percentile quantile regression, not least squares.**
This is the part worth arguing with. Least squares runs the law through the
middle of the corpus, so half of every acquisition's peaks are narrower than it
predicts; on this corpus a least-squares pair leaves 12% of the 311 measured
peaks under two bins per width and 70% under three. The two errors are not
symmetric. A peak
**broader** than the law gets more bins than it needs, which costs store and
nothing else. A peak **narrower** than the law is the one case that loses
shape, and that is the defect being fixed. So the pair is placed at the narrow
edge: quantile regression at `q = 0.05` on `FWHM^2 = A m + B m^2`, which puts
94% of measured peaks at three bins per width or better and none below two, for
20% more bins than the least-squares pair. The old `linear_tof` default left
87% under three bins and 42% under two. (Every figure here counts each of the
311 peaks once; weighting one vote per acquisition instead gives 54%, 4% and
83% / 35%, with the same ordering.)

**Why one law for an instrument whose tune varies.** Resolving power varies by
a factor of three between acquisitions in the corpus, and a C60 primary beam
resolves about five times worse again (R 280-1,200 against 1,700-7,100 for
Bi3). The header does not distinguish these reliably -- `AcqPulseWidth` is 16.0
on both the sharpest and among the broadest files measured -- so a
mode-switching detector would be guessing. One law, placed at the sharp end,
converts every acquisition correctly; a poorly resolved one is simply
oversampled. This matches how the MRT and timsTOF pairs were derived.

**Why `get_reference_width` stays `None`.** A `tof` axis is sized in bins per
peak width, so the law and `DEFAULT_BINS_PER_FWHM` already fix the bin width at
every m/z. Declaring a width as well would be a second spelling of the same
quantity, and `reference_params` answers for a `tof` axis before it ever
consults a detector's width -- it would be dead code. Same reasoning as
`WatersMRTCentroidDetector`.

**Why resampling stays on by default.** `--no-resample` was the other
candidate, and on the *old* axis it was the better store. On this one it is
not: the resampled store is 22% smaller, has 8.4x fewer columns, opens about
1.8x faster (0.61 s against 1.08 s, medians of four alternating warm reads --
the first cold read of the native store took 3.0 s, which is not a fair
comparison and is not the number used here), and comes out better correlated
with the vendor's peak images than the unresampled one. Its remaining recovery
spread is window-edge quantisation -- the vendor's windows span 3 to 10 bins
and their edges fall mid-bin, and the windows with the fewest bins are exactly
the ones with the worst recovery -- rather than lost peak shape. `--no-resample`
stays the way to get `var["tof_us"]` and the raw channel grid.

**The known limit.** The corpus is one instrument. Every acquisition came from
the same nanoTOF, so the pair describes that instrument's range of tunes and
not the model line. A second instrument would be worth measuring before
treating these two numbers as a nanoTOF constant. The high-mass end rests on
few peaks -- 8 above m/z 200, all from polystyrene standards -- so `B` is
anchored mostly by the m/z 60-200 plateau.

---

## D21. The nearest-neighbour bin index is computed, not searched

**Status:** Implemented (2026-09-15).

**Decision.** Every axis generator records the coordinate it lays its bins
in -- `forward()` on `BaseAxisGenerator`, carried as an `AxisLinearisation`
on the `MassAxis` it returns -- and the converter keeps that beside the
common mass axis. At the two converter call sites of `_nn_map_to_bins`
(`_nearest_neighbor_resample` and `_build_nn_shared_cache`) the nearest bin
of a peak is then computed: round its position `(forward(mz) - u0) / du`,
take the nearest of that bin and its two neighbours with the search's own
comparison and tie rule. The four sibling-table sites (the heatmap, the
mobility grid's discovery and scatter, the MS/MS split) keep
`np.searchsorted`. At axis-build time the converter measures how far the
axis deviates from its own linearisation and uses the closed form only
under a quarter of a bin on a strictly ascending axis; otherwise, and for
every axis that has no law (`--no-resample`, an axis passed bare), it
searches. Issue #295.

**Why the repair is exact, for every law.** The objection raised against
the issue's proposal was right on its facts: only `constant` is a
`np.linspace` in its own coordinate. The five physics generators lay a
uniform grid of bin *edges* `u_i = u_0 + i du` in `u = forward(m/z)` and
report centres `c_i = (m_i + m_{i+1}) / 2` -- arithmetic midpoints in m/z,
not in `u` -- so `forward(c_i)` is not `u_0 + (i + 1/2) du`, and a
rounded position can be off by one. The empirical worst deviation found
then, 0.478 of a bin on an FT-ICR axis of 1,000 bins over `[1, 10^6]`, sat
4.5 percent from the half-bin cliff past which a +-1 repair returns a
wrong bin, and nothing in the conversion would notice.

The cliff is unreachable, and the reason is one line. `m_i < c_i <
m_{i+1}` strictly, and `forward` is strictly monotone, so `forward(c_i)`
lies strictly between `u_i` and `u_{i+1}`: strictly within half a step of
the linearised centre, for every bin of every axis any of these generators
can lay. Per law the deviation has a closed form in the edge ratio
`t = m_{i+1} / m_i`, verified symbolically (sympy) and numerically on
541 million probes over 305 configurations:

| law | `forward` | deviation / `du` | limit as `t -> inf` |
|---|---|---|---|
| constant | `m` | 0 (float rounding only) | 0 |
| linear_tof | `sqrt(m)` | `(sqrt((1 + t^2) / 2) - (1 + t) / 2) / (t - 1)` | `1/sqrt(2) - 1/2 = 0.207` |
| reflector_tof | `ln(m)` | `ln(cosh(h)) / (2h)`, `h = ln(t) / 2` | 1/2 |
| orbitrap | `1 / sqrt(m)` | `((1 + t) / 2 - sqrt(2) t / sqrt(1 + t^2)) / (t - 1)` | 1/2 |
| fticr | `1 / m` | `(t - 1) / (2 (t + 1))` | 1/2 |
| tof | `(2 / sqrt(B)) asinh(sqrt(B m / A))` | the general bound | 1/2 |

Every entry is strictly below 1/2 for `t > 1`, and approaches it only as
one bin spans an unbounded ratio -- two bins over `[1, 10^8]` gives
`0.5 - 2 x 10^-8`. On a real axis the deviation is second order in one
bin's relative width: `4 x 10^-3` at worst over the realistic
configurations swept, `6 x 10^-7` on the reflector-TOF axes the corpus
resolves to.

The composed index error follows. A value's true nearest bin `j` has the
value in `[c_{j-1}, c_{j+1}]`; `forward` is monotone, so its position lies
between those two centres' positions, each within `d < 1/2` of its own
index; so the position lies in `(j - 3/2, j + 3/2)` and rounds to `j - 1`,
`j` or `j + 1`. Among those three the closed form makes the reference's
own comparison in m/z, moving left only on a strict improvement, so a tie
resolves to the right as it always has. A value in the half-bin skirt
beyond either end rounds to `-1` or `n` and is clipped to the edge bin,
which is where the search's clip sends it too.

**Why a guard anyway.** The proof is about the generators; the converter
should not have to trust that a `MassAxis` came from one. So
`_usable_linearisation` measures the deviation on the axis actually built
and requires it under `NN_LINEARISATION_MARGIN = 0.25` bins, on a strictly
ascending axis (a duplicated value must still map to its first occurrence,
which is what `np.searchsorted` does and the repair does not). The margin
is not the proof's bound but a quarter of it, so float rounding in the
position has nowhere to matter; it refuses two bins over six decades and
nothing an instrument asks for. The check is chunked like the bin-width
log line beside it, for the same 200-million-bin reason (#251).

**Measured** (2026-09-15/16, this machine, every dataset in the local
corpus, medians of warm runs, first run of each pair excluded as the page
cache's first touch; every run's store hashed and identical across
variants; `handouts/closed-form-bin-index.md` has the harness, the
spread and every run):

| dataset | bins | mean pts/spectrum | current | two sites | gain | six sites | gain |
|---|---|---|---|---|---|---|---|
| `20240826_xenium_0041899.imzML` (Xenium, centroid) | 313,713 | 1,273 | 414 s | 276 s | +33% | 276 s | +33% |
| `20240826_Xenium_0041899.d` (timsTOF, TIMS off) | 313,723 | 1,273 | 626 s | 477 s | +24% | 478 s | +24% |
| `pea.imzML` (flex, centroid) | 460,495 | 4,722 | 15.1 s | 9.7 s | +36% | 9.7 s | +36% |
| `20231109_PEA_NEDC.d` (timsTOF, TIMS off) | 541,610 | 2,321 | 31.7 s | 22.5 s | +29% | 22.3 s | +30% |
| `bellini.imzML` (ToF-SIMS, centroid) | 2,373,513 | 2,220 | 18.8 s | 9.3 s | +50% | 9.4 s | +50% |
| `06_glycans_highmass_66k_px` (TIMS slide) | 240,794 | 66,488 | 2,002 s | 1,716 s | +14% | 1,633 s | +18% |
| `05_maldi2_shortramp_26k_px` (TIMS slide) | 138,629 | 34,243 | 366 s | 327 s | +11% | 314 s | +14% |
| `04_ratbrain_71k_px_9990_scans` (TIMS slide) | 183,258 | 1,433 | 128 s | 118 s | +7% | 114 s | +11% |
| `03_biofilm_maldi_vs_maldi2_20um` (TIMS) | 92,861 | 6,868 | 80 s | 72 s | +10% | 69 s | +13% |
| `02_tiny_longramp_1465px`, `01_tiny_msms_315px` | 21k, 599k | | 27 s | 27 s | +1 to +2% | 26 s | +2 to +4% |

Every set resolved to `reflector_tof`, so the other five laws were forced
on `pea.imzML` and the biofilm set: all six convert bit-identically, and
gain in proportion to their bin count (+22 percent for `linear_tof` at
44k bins on pea, +36 for `fticr` at 1.8M). The issue's own mock at Xenium
density reproduces the 2026-09-14 sanity point: 17.7 s to 13.4 s, +24
percent against the +27 reported then.

**Why not the four sibling sites.** Profiled on the current tree, warm,
the two converter sites are 11 to 16 percent of a TIMS whole-slide wall
clock and 30 to 40 percent of an imzML or mobility-free TDF conversion;
the sibling sites are 5 to 7 percent of a TIMS slide and nothing anywhere
else, because they map each frame's unique m/z once (pass 1, for the
heatmap) where the summed table maps it twice. The closed form there is
worth 3 to 4 points on a slide (the six-site column above), through one
optional argument on `SiblingPasses`, `MsmsAccumulator` and the three
mapping helpers. It is left out of this change for a better reason than
size: under `scan_sum` a TDF frame's summed spectrum *is* its `unique_mz`,
so the sinks re-map an array the summed table has just mapped, and
sharing that mapping per frame would remove the sibling search rather
than speed it up. Both are recorded on #295.

**Strongest objection.** That the win was measured on a mock ten to forty
times sparser than the corpus, and was a net loss there. It was: at 124
peaks per spectrum the closed form's fixed cost -- three gathers and two
comparisons over every probe -- exceeds the search it replaces, and the
issue's own "10 percent of conversion CPU" came from that mock. The numbers
above are the corpus, and at the corpus's peak counts the search is 85
percent of the mapping and the mapping is a quarter to a third of the
conversion. The premise objection -- "not a linspace, and 4.5 percent from
a cliff" -- is answered above by derivation rather than by more sampling.

**Known limit.** A `MassAxis` built by hand carries no linearisation and
searches. The sibling tables search. `--no-resample` maps through a
different function altogether (`_map_mass_to_indices`, an exact-match
search over the raw union axis) and is untouched.

---

## D22. The converter holds one resampling strategy

**Status:** Implemented (2026-09-17).

**Decision.** A conversion that resamples builds exactly one
`ResamplingStrategy` when it builds its mass axis, and calls it once per
spectrum per pass. The strategy is the whole per-spectrum operator:
`NearestNeighborStrategy` bins, `TICPreservingStrategy` interpolates and
rescales, both answer `resample(mzs, intensities)` with the bins that
spectrum fills and nothing else, and both carry the out-of-range counter
and its warn-once line. The converter reaches them through one factory,
`build_strategy(method, axis, axis_range, linearisation, gap_tolerance_da)`
in `thyra.resampling.strategies`, and holds the result as `_resampler`.
Six private methods and five pieces of state leave
`base_spatialdata_converter.py` with them. `--no-resample` is not a
strategy: it maps onto the reader's own axis through
`_map_mass_to_indices`, which stays on the converter. Issues #277 and #352.

**Why the package that already bore these names could not be wired in.**
It did the opposite thing with intensity. The old `NearestNeighborStrategy`
*interpolated*: each target point took the intensity of the nearest source
point, so a lone source peak was copied onto every target point closer to
it than to any other, and the summed intensity scaled with how densely the
target axis was laid. The conversion has always *binned*: each source peak
lands in exactly one target bin and the total is preserved. The two classes
shared a name, a module and nothing else -- the converter never imported
them, `docs/api.md` never rendered them, and the old `TICPreservingStrategy`
was a second copy of the interpolating operator that had already drifted
(it lacked the declared `axis_range`, so it still carried #239's edge-bin
bug after the converter's copy was fixed). Making the strategies real
therefore meant replacing the class bodies with the converter's, not
teaching the converter to call what was there. The meaning of the name
changed; with zero consumers that is a changelog line rather than a break.

**Why the operator moved before the axis planner.** #276 splits this
converter along two seams: what decides the axis, and what places a
spectrum on it. The operator is the smaller half and the one with a hard
acceptance test -- every stored byte unchanged -- so it went first, in two
steps: PR #366 moved the pure functions into `thyra/resampling/`, this one
moved the state and the methods that read it. What the split buys is the
contract above. `build_strategy`'s middle three arguments are exactly the
triple the `AxisPlanner` of #352 returns, so the planner could be lifted
out without the operator noticing, and the operator became testable
without standing up a converter: the test files that drove it stopped
posing as one through `SimpleNamespace` and `MethodType` and construct a
strategy instead.

**Why the contract is sparse.** `resample` returns unique ascending bin
indices and their non-zero values, never an array of the axis's length.
The dense form of the TIC-preserving method was 35x the cost of reading a
zero-suppressed profile source, in the words of the comment that survives
on `_process_spectrum`: a Waters MRT pixel stores about 15,000 samples
against a 1.05M-bin axis, and every other bin interpolates to exactly
zero. The dense array is available from `to_dense` for the callers that
want one, so no strategy has to build it for a caller that does not.

**What was measured.** Store identity, with `tests/tools/store_identity.py`
(#349), against `origin/main` at `8f16d89`: `pea` (imzML, nearest
neighbour, shared-axis cache), a TDF with mobility sibling tables, a
TIMS-off TDF, the same imzML at `--mass-axis-type constant` (the route
#356 made linearised), and the TDF again at `--resample-method
tic_preserving`. Both strategies, both mapping routes -- closed form and
search -- and the sibling sinks under an interpolated conversion are
covered. Five dataset-runs, and zero Zarr elements differ in any pair.
That is the acceptance criterion for the whole #276 decomposition, and it
is what makes the change of meaning above safe to assert rather than hope.

Cost on the hot path: `pea.imzML` with `--no-optical`, five warm runs per
side, interleaved so machine drift cancels -- 12.13 s median here against
12.28 s on `main`, a 1.2% difference in this change's favour. That is
noise, and the honest way to say so is to report what measuring it badly
gave: run in blocks instead of interleaved, minutes apart, the same pair
came out 1.8% the other way. It should be a wash either way. The call the
strategy replaced was already one Python-level dispatch per spectrum --
the `_nn_route` flag `__init__` resolved once for exactly that reason --
so the bound method on `_resampler` costs the same lookup and everything
inside it is the code that was already there. The one real difference is
that the kept range is computed once when the strategy is built rather
than once per spectrum.

**What the strategy is not given.** The axis's linearisation is not the
strategy's to own, even though the binning strategy is the only thing that
uses it to resample. It is a property of the axis -- one of the triple
`AxisPlanner` returns -- and the sibling sinks place peaks onto that
same axis whichever method the spectra took, including an interpolated
one. So the converter keeps it as `_axis_linearisation`, assigned and
cleared beside the axis itself, and passes a copy to `build_strategy`.
Deriving it back out of the strategy instead was tried and rejected: it
reads as tidier and quietly makes a sibling table written beside a
`tic_preserving` conversion search where it used to compute, which is the
same bin (D21) at a cost no measurement would explain. The sinks still map
peaks themselves rather than through `NearestNeighborStrategy.map_to_bins`;
#353 closes that by giving them the strategy.
