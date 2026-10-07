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
it stated. *Proposed* means the shape is written down but nothing is decided:
the entry states what it waits on, and until that arrives it binds nothing.

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

---

## D23. A metadata document does not need a pixel size (PROPOSED)

**Proposal.** Split `msi_metadata` into a **core** every mass-spectrometry
acquisition can fill and an **imaging profile** that adds what only a raster
has. `ms_analysis.pixel_size_um`, today the schema's one required field,
moves into the profile. A document that states no profile is a valid core
document; a document that states the imaging profile must carry the pitch,
exactly as every document does now.

**Status:** Proposed 2026-09-18, not decided. It waits on a partner schema
(below) so that the split lands where a layered model already in progress
expects it, rather than being guessed at and then reconciled.

### What forced it

`thyra metadata` can now build the block from a source that will never be
converted, and the first such source read was a Bruker `.d` from a run that
imaged nothing. Everything the schema asks about the mass spectrometry was
there -- instrument model, positive polarity, a 100 to 8000 mass range, one
precursor at m/z 3888 isolated over a 10 Da window at two collision
energies, the acquisition timestamp with its UTC offset, the method file
name. What was not there was a pixel size, because there is no raster, and
no number would be one: the document is invalid by construction, on its one
required field, for a reason that is a property of the acquisition rather
than a gap in it.

That is not a bug in the document and it is not fixed by making the field
optional. "Optional" would say a pitch may be missing from an *image*, which
is exactly the thing conversion refuses to do (see
[CLI Reference](cli.md#exit-status)) and exactly the field METASPACE
requires. The distinction the schema is missing is not present/absent, it is
*which kind of document this is*.

### The shape proposed

| Section | Core | Imaging profile adds |
|---|---|---|
| `sample` | whole section | -- |
| `preparation` | whole section | -- |
| `ms_analysis` | polarity, ionisation source, analyzer, instrument model, manufacturer, serial number, spectrum count, resolving power, ion mobility, fragmentation | `pixel_size_um` |
| `acquisition` | `acquisition_datetime`, `method_file` | `laser_power_percent`, `laser_frequency_hz`, `shots_per_pixel` |
| `calibration` | whole section | -- |
| `alignment` | -- | whole section: it registers a raster onto an optical image, and a run with no raster has neither |
| `processing` | whole section | -- |
| `provenance` | whole section | -- |

Three properties the split has to keep:

1. **Every document Thyra writes into a store is an imaging document.** A
   conversion refuses without a pitch, so nothing already on disk changes
   meaning and nothing already written stops validating.
2. **The METASPACE export stays a mirror of the imaging profile.** Its
   `Pixel_Size` comes straight off `ms_analysis.pixel_size_um`, and
   METASPACE is a platform for imaging data; an export of a core-only
   document is not a submission with a field missing, it is not a
   submission. The export should refuse it and say which profile it needs.
3. **The shots-per-pixel fields go with the raster, not with the laser.**
   `shots_per_pixel` is a count *per position*, which a run with no
   positions cannot state; `laser_power_percent` and `laser_frequency_hz`
   are settings an untargeted MALDI run has as much as an imaging one, but
   they are reported today only by the three formats' imaging tables, and
   promoting them into the core would put fields in it that nothing can
   fill. They stay in the profile until a source fills them without a
   raster.

### The objection, and why it does not win yet

**"Just make the field optional and be done."** One line of model change
against a table of sections. It loses on what the document then says: an
optional pitch makes "this is an image and nobody recorded the pitch"
indistinguishable from "this is not an image", and the first of those is a
defect a consumer must not silently accept. A profile makes the two
different documents. The cost is that a consumer must look at which profile
a document declares before reading `ms_analysis`, which is the cost of the
distinction being real.

### Why it is not decided here

The community this came from is drafting a layered metadata guideline for
native and ion-mobility MS -- a base schema for mass spectrometry with
technology extensions on top -- and Thyra's schema is a candidate to be
referenced from it. A split designed against that draft is one schema in two
layers; a split designed ahead of it is two schemas that will have to be
reconciled, and the reconciliation would land on stored documents. The draft
has not arrived. So this entry records the shape and the constraints, and
the change waits.

**Known limit while it waits.** `thyra metadata` on a source with no raster
writes a document, reports `ms_analysis.pixel_size_um: Field required`, and
exits 1. That is honest and it is not usable in a pipeline that gates on the
exit status. Anyone who needs such a document to pass validation today has
to supply a pitch that is not a measurement, which is the thing this schema
exists to avoid.

---

## D24. The converter is the two-pass loop and four collaborators

**Status:** Implemented (2026-09-22), issues #276 and #349 to #354. The
one question the decomposition left open, whether the two converter
modules fold into one, is **deferred** below with the condition that
reopens it.

**Decision.** The SpatialData converter is the lifecycle and the write
path (`base_spatialdata_converter.py`), the two-pass loop
(`streaming_converter.py`), and four collaborators the converter holds
rather than inherits, each constructible and testable without a
converter:

| collaborator | home | takes | issue |
|---|---|---|---|
| `OpticalImages` | `converters/spatialdata/optical_image.py` | the reader, a store-path accessor, the dataset id, a pitch accessor | #350 |
| `UnsAssembler` and `RootAttrsBuilder` | `converters/spatialdata/uns_assembler.py`, `converters/spatialdata/root_attrs.py` | the reader and what is settled at construction; everything decided during a conversion arrives per call in a context | #351 |
| `AxisPlanner` | `resampling/axis_planner.py` | a reader, or the metadata dict the detectors read, and a `ResamplingConfig` | #352 |
| `SiblingTables` | `converters/spatialdata/sibling_tables.py` | the reader, a store-path accessor, the sibling options; the axis and grid per call in a `SiblingContext` | #353 |

Beside them the per-spectrum operator became a `ResamplingStrategy`
(D22, #277). What is left on the converter is what the plan on #276 said
should be: `__init__`, `convert`, `_initialize_conversion`, the mass-axis
guards, the `var` frame, the pixel shapes and region numbers, the three
context packers, `_save_output`, and the two-pass loop with its table
units.

Two rules the collaborators follow, each learned from a defect on the
way. **A collaborator the converter can replace is read per call, never
captured**: `OpticalImages` captured the converter's pitch accessor in
its constructor and a test that swapped `converter.optical` reproduced
issue #288, so every delegate travels in a context packed fresh at the
call. **A deferred import is a patch seam, not a cost**: the sibling
builders are imported inside the methods that call them because the #343
tests patch `mobility_table.build_mobility_table`, and a module-level
import binds the name before the patch can (measured, the whole import
chain is 2.5 ms; cost was never the reason).

**What was measured.** Every step ran `tests/tools/store_identity.py`
(#349) against `main` at its base, every Zarr array and attribute
hashed, and the recount below was taken from the AST at each merge, not
from reading the modules.

| step | merged as | `base_spatialdata_converter.py` lines / methods / broad catches | new home (lines) | store identity against its base |
|---|---|---|---|---|
| filed, `97184f5` | | 4,174 / 85 / 27 | | |
| planned, `f4c11b7` | | 4,628 / 90 / 20 | | |
| 1, optical | PR #361 | 4,098 / 79 / 17 | `optical_image.py` 710 to 1,451 | 13 stores, 0 differences |
| #277 A, the pure operator | PR #366 | 3,666 / 79 / 17 | `resampling/binning.py` 307 | 5 dataset-runs, 0 |
| #277 B, the strategies | PR #374 | 3,383 / 73 / 17 | `resampling/strategies/` 507 | 5 dataset-runs, 0 |
| 2a, the `uns` block | PR #375 | 2,941 / 60 / 14 | `metadata/uns_assembler.py` 668 | 3 datasets, 0 |
| 2b, the root attrs | PR #376 | 2,648 / 53 / 13 | `metadata/root_attrs.py` 374 | 3 datasets, 3 differences, all the root attrs #67 item 1 added on purpose |
| 3, the axis planner | PR #377 | 1,760 / 40 / 9 | `resampling/axis_planner.py` 1,198 | 13 dataset-runs, 0 |
| 4, the sibling tables | PR #378 | 1,125 / 26 / 3 | `sibling_tables.py` 977 | 10 stores, 0 |

`streaming_converter.py` is 1,078 lines and 24 methods at `f4c11b7` and
1,079 and 24 after step 4: the loop was already what it should be. The
broad `except Exception` catches went from 27 in one file to 3 there and
20 across the modules above, which is the 20 the plan started from: they
moved with their methods, none was added and none was removed by the
decomposition (narrowing them was #280). The lines across the same
modules grew, 6,715 to 7,996, and the growth is the contexts, the
docstrings and the tests that build each collaborator alone. Line count
was never the point; that those tests exist is.

**The fold, deferred.** `BaseSpatialDataConverter` is an ABC with one
subclass (D11), and the plan asked whether the two modules should still
be two once the collaborators existed. The measurement that decides it is
what the subclass still reaches into: every `self.<name>` in
`streaming_converter.py` defined on the base chain and not in the
subclass. At `f4c11b7` that was 27 names over 47 sites; at `88d926a` it
is **14 names over 35 sites**, and the 13 names that went were all
collaborator extractions. Of the 35 sites, 23 are four pieces of
per-conversion state, `_dimensions` (10, unchanged since `f4c11b7`),
`_common_mass_axis` (5), `_resampler` (5) and `_region_map` (3), and most
of the other twelve are template-method calls a subclass is meant to make
(`_create_pixel_shapes`, `_create_mass_dataframe`,
`_drop_unusable_intensities`, `_coalesce_duplicate_bins`). The coupling
that remains is data, not behaviour: two peers of the same size sharing
four variables that both genuinely need.

So the two modules stay two. The reach-ins that made the split look wrong
when #276 was filed were symptoms of the missing collaborators and went
without a fold; what a fold would buy now is the removal of two abstract
hooks and one `super().__init__`, at the cost of a 2,200-line module
whose two halves have different reasons to change (the lifecycle and
the write path change when spatialdata does; the loop changes when the
scatter does). The condition that reopens this: a change that needs a
fifth piece of shared state, or a second subclass, or a reach-in count
that grows past 35 again. The response then is to lift the state into a
collaborator of its own (the four variables are the shape of one: the
axis, the grid and the regions of one conversion), not to fold.

**Objections considered.**

- *Mixins would have cut the file at a tenth of the cost.* They cut the
  file, not the coupling; nothing becomes constructible alone, and the
  test each step added, the collaborator built from arguments with no
  converter anywhere, is the thing the decomposition was for.
- *The store-writing assemblers landed in `thyra/metadata/` but only a
  conversion runs them.* True by dependency: their inputs are converter
  types and their only caller is the converter. They landed there because
  they assemble metadata, and they were moved beside the converter once
  #382 stated the boundary (D25). The boundary that matters, that nothing
  on the document side imports them, held before the move and is asserted
  by `tests/unit/test_import_boundaries.py`; the move was a rename against
  that test, not a design question.

**Known limit.** The four shared variables are assigned in the base's
`_initialize_conversion` and read in the subclass's passes, so a test of
the loop still needs a converter with a reader; the collaborators are
testable alone, the loop is not yet.

---

## D25. The metadata document depends on nothing that reads spectra

**Status:** Implemented (2026-09-22) for the boundary and the address;
the two structural changes it names are **deferred**, each with the
condition that triggers it. Issue #382, on the findings of #381, #387 and
#379.

**Decision.** The schema (`thyra.metadata.schema`, `thyra.metadata.ontology`,
`thyra.metadata.types`) imports nothing from the readers, the converters
or the resampling package, and `tests/unit/test_import_boundaries.py`
asserts it in a fresh interpreter, at import time and again while a
document is built. Anything that can describe an acquisition may produce
a document: what the builder consumes is the four dictionaries of
`ComprehensiveMetadata`, an optional pixel size with its provenance, a
source format name and an optional fragmentation report, not a
`BaseMSIReader`. The store-writing assemblers, `UnsAssembler` and
`RootAttrsBuilder`, are converter collaborators and live with the
converter under `converters/spatialdata/`. The schema stays inside the
`thyra` distribution until a second program commits to writing the
document; then, and not before, it becomes a distribution of its own,
which the boundary test makes a mechanical move.

**Measurement.** Through the package front door the schema imported in
2.7 s with 2,822 modules, 32 of them readers and 8 converters, plus
spatialdata, dask and anndata. With the front door out of the way (#381,
PR #386) the same import takes 0.14 s, 215 modules, sixteen of them
Thyra's, and one third-party package, pydantic. The document builder
reads exactly one field of `EssentialMetadata`, the source path. Nothing
outside the repository imports the schema API: the one downstream
application reads the block out of the store as a dictionary. The
committed JSON Schema had no `$id` and the LinkML file no address; both
are published at a versioned URL now (#387), which is what a guideline
cites and what a non-Python implementer consumes.

**Three layers by dependency, two by packaging.** Data (readers,
converters, resampling), metadata (extractors, `ComprehensiveMetadata`,
`document.py`) and document (`schema/`, `ontology/`) point one way, and
the test pins the direction. They are one wheel. The middle layer's
contract, a `MetadataSource` protocol of four methods that every reader
already satisfies (`get_comprehensive_metadata`, `has_fragmentation`,
`get_fragmentation`, `close`) and a `document_from_source` entry point
beside `read_metadata_document`, is **deferred until the first non-reader
source is written**: a protocol with one implementer is a class with a
docstring, and every abstraction in this repository built before its
second caller has been deleted or restated (D7, D10, the routers of D11).
The vendor extractors stay where they are for the same reason: none of
them needs the schema or a reader class, so moving 3,333 lines buys
nothing the boundary test does not.

**Which axis each change moves.** This entry is about how the *code* is
layered. How the *document* is layered, a core every acquisition can fill
and an imaging profile that adds the raster, is D23 and moves on a
`schema_version` bump, not a package change. The one place they touch: a
core-only document is what a non-imaging source produces, so D23 is what
makes the protocol above worth having. A dependency, not a merger.

**The objection, and why it does not win.** Split the distribution now,
because a guideline will not reference a subpackage of a converter. What a
guideline references is a versioned schema at a stable address, which
exists, and what a non-Python implementer consumes is the JSON Schema and
the LinkML source, neither of which is a wheel. A second wheel before a
second writer is an orphan package with one dependant and two release
trains for one maintainer: the release automation stamps one version into
one `pyproject`, and a workspace needs a changelog, a compatibility policy
across packages and someone to approve a schema change that is not a
Thyra change. The `schema_version` inside the document is already
versioned independently of the package, so the rule exists; the
governance does not, and it is not free.

**Known limit.** `ComprehensiveMetadata.essential` is required, and its
docstring promises fields a source with no spectra has no honest value
for (dimensions, bounds, counts). The first non-reader source makes it
optional on the document path, one branch in a builder that already
tolerates `comprehensive is None`. Until then `thyra metadata` on a
non-imaging source describes it correctly and fails validation on the
pixel size, which is D23's limit, not this one's.

---

## D26. The acquisition order is an `obs` column, when the source has one

**Status:** Implemented (2026-09-28). Schema 0.10.0.

**Decision.** Every table's `obs` carries `acquisition_order` when the reader
knows in which order the spectra were acquired, and leaves it out when it does
not. The column is `int64`, grows with acquisition time and differs from row to
row. It is the source's own number for a spectrum where the source numbers them
(Bruker `Frames.Id`, solariX `Spectra.Id`), else the spectrum's 0-based position
in the order the source lists them (the imzML spectrum list, mzPeak's
`spectrum_index`, Waters scans numbered on through the converted functions). A
row summed from two measurements takes the earlier one. The rows stay in grid
order. A reader reports the order through `has_acquisition_order` and
`iter_spectra_with_acquisition_order`, and the TDF frame record carries it as
`acquisition_order`, so the fused passes of D5 keep it too.

**Why.** The rows follow the raster, and nothing in a store said when each pixel
was measured. A serpentine scan comes back along every second row and a slide
of several regions is measured one region after another, so a value plotted
against the row number is plotted against the raster, not against time. Two
consumers need time. A QC view plots per-pixel values against acquisition order
to show drift during the run, as SCiLS Lab does; Ousia's QC view labels its axis
"raster order" because the store offered nothing better. And a per-spectrum
correction that falls back to "the previous successful spectrum" needs
"previous" to mean previously acquired. The Bruker reader already walked the
frames by `Frames.Id` and dropped the id at the last step.

**Why the source's own number and not a rank.** A rank from 0 is what a plot
wants, and a consumer gets it with one `argsort`. The source's number gives the
same order, finds the spectrum again in the source, and survives a `--region`
conversion unchanged: a region's frame ids are not contiguous in the file, and
renumbering them would make two stores of one slide disagree about the same
frame. The cost is that the column is not a row index and, for Bruker, not
0-based; the docs say so where the column is described.

**Why `int64` and not a nullable integer.** pandas' nullable `Int64` was the
obvious alternative, and its null would be unreachable. A row exists only
because a spectrum reached it, and every spectrum from a reader that knows the
order carries one; the converter refuses a spectrum without it rather than
leave a blank. "Unknown" is therefore a property of a whole table, and the
store says it by leaving the column out, as it says "not stated" everywhere
else. Plain `int64` is also what `x`, `y` and `region_number` are, and it is
stored as one flat array, where anndata stores a nullable integer as values
plus a mask that any reader of the Zarr without anndata has to know about. The
price: concatenating a table that has the column with one that has not gives
pandas' `float64` with NaN, unless the caller asks for `Int64`.

**Why this name.** `spectrum_index` was the other candidate. mzML and mzPeak use
it for a 0-based list position, while the Bruker value is a 1-based id with
gaps, so it would mislead exactly the readers who know the name.
`acquisition_order` says what the column is for and claims nothing about where
it starts. The spec reserves it beside the `var` names
(`MSI_OBS_ACQUISITION_ORDER_COLUMN`), and adding an optional name moves
`schema_version` to 0.10.0, a minor bump by the versioning rule. The document
itself is unchanged: the 0.10.0 JSON Schema differs from 0.9.0 only in its
version string.

**What was measured.** Each "grows with time" was checked read-only on real
acquisitions, 2026-09-28:

| Source | Value | Measured |
|---|---|---|
| Bruker timsTOF | `Frames.Id` | `Frames.Time` never decreases along it, on two TDF (713 and 11,752 frames) and two TSF (33,800 and 918,855 frames) |
| Bruker solariX | `Spectra.Id` | On 50 acquisitions of 4,270 to 30,824 spectra, the per-scan `DateTime` never decreases along it, nor the `minutes` of `ImagingInfo.xml` along its own list |
| Waters | scan number, carried on through the converted functions | Retention time never decreases along it on three runs of 2,436, 5,400 and 7,683 scans; the last is one raster split across three functions, whose time ranges follow each other |
| imzML | position in the spectrum list | Not measurable, as the file records no time. The one export available beside its source (918,855 spectra) lists them in frame order, one to one, but that run was rastered row by row, so it cannot tell "as acquired" from "raster order" |
| mzPeak | `spectrum_index` | Not measured; it is the archive's own order |
| PHI | absent | Each frame passes over the whole raster and a pixel sums every frame |
| Bruker rapifleX | absent | See below |

**What a store gains, measured.** Fifteen datasets were converted before and
after the change: the eight fixtures, two imzML exports, a TSF, a PASEF TDF,
two Waters runs and a PHI run. They differ in 49 places, all of three kinds:
each table of an ordered source gains the array and its name in the `obs`
column list, and every table's `schema_version` moves. The matrix, `var`, the
other `obs` columns, the images, the shapes and the root attributes are
byte-identical. On the TDF and the TSF every row holds the frame id
`MaldiFrameInfo` records at its position; the TSF has 33,690 rows for 33,800
frames, and the ids of the frames without a row are simply absent.

**Why rapifleX has none yet.** The reader walks the raster offset table of the
`.dat`. The `_poslog.txt` it also parses is a timestamped list of positions, and
very likely the acquisition order, but no real rapifleX acquisition was
available to check it, and a synthetic file can only confirm that code reads a
field, never what the field means (D7). One real acquisition reopens this.

**Why the earlier of two measurements.** A position measured twice is one row
holding the sum of both spectra (issue #241), so the row has no single time.
The earlier is when the pixel was first measured, and it does not depend on the
order a reader yields the two. It happens on real data: the split Waters run
above records one stage position twice.

**Known limit.** Two converted Waters functions that cover the same pixels are
summed (D7), and their scans are numbered one function after the other rather
than interleaved. Each such row takes the first function's number, which keeps
the rows in time order as long as the first function visited all of them. No
real file with two such functions converted was available.

## D27. One pixel size in an mzPeak archive is tested before it is believed

**Status:** Implemented (2026-09-29).

**Decision.** The mzPeak extractor reads a pixel size in one of two ways.

- **A declared pair.** `IMS:1000046` and `IMS:1000047` under their current
  names are two lengths and are taken as written, each converted from its
  unit. The declared extent is not held against them.
- **One number.** `IMS:1000046` without `IMS:1000047`, or under the name
  "pixel size", is tested against the declared pixel count (`IMS:1000042`,
  `IMS:1000043`) and extent (`IMS:1000044`, `IMS:1000045`). It is a length if
  value x count is the extent, and an area if sqrt(value) x count is the
  extent, in which case the side of the pixel is sqrt(value). The tolerance
  is 1% of the extent. Every axis that gives both a count and an extent has
  to agree. The pixel is square.

When the test cannot be made, or fits neither reading, the extractor returns
no pixel size and logs the reason at WARNING. The conversion then asks for
`--pixel-size`, as it does for the same file as imzML. An area that was read
is logged at WARNING too, with the value in the file and the side taken from
it.

A value with no unit is still read as micrometres, in both cases. Centimetre
stays out of the unit table. `IMS:1000047` without `IMS:1000046` is not a
pixel size.

The imzML path makes no test. It does not read `IMS:1000046` under the name
"pixel size" at all, with or without `IMS:1000047` beside it. The file gives
no pixel size, the log says why at WARNING, and the conversion asks for
`--pixel-size`.

**Why.** The term changed its meaning. Until commit `421481e` of the imzML
vocabulary (2017-09-07) `IMS:1000046` was named "pixel size" and gave the
area of a pixel, and `IMS:1000047` was "image shape". Since then they are
"pixel size (x)" and "pixel size y", two lengths. Files of the old form are
still published, and the reference mzPeak converter copies the parameter
into `scan_settings_list` unchanged: name, value, and a null unit.

Thyra 4.2.0 read a lone `IMS:1000046` as a length in micrometres and copied
it to the other axis. An archive that says pixel size 10000, 5 pixels on x
and extent 500 um came out with a pixel of 10000 um by 10000 um. The pixel
is 100 um. `convert_msi` returned `True` and nothing was logged.

**Why the numbers decide and not the name.** The name is what a writer
typed. The count and the extent are a second statement of the same
geometry, made by the same file. A writer may keep the old name and mean a
length, and a lone term under the new name may still hold an area, so the
name only says when to test. The test then says which vocabulary the file
speaks. If it is a length and the file gives `IMS:1000047` as well, the pair
is read as declared.

**Why the imzML path refuses on the name alone.** That path asks for both
terms, so a lone `IMS:1000046` never gave a pixel size there. One form got
through: the old name with a number on `IMS:1000047`, read as two lengths.
A file that says pixel size 10000 beside such a number was converted with a
pixel of 10000 um by 10000 um. Refusing the old name closes that and changes
nothing else: every real imzML at hand and every sampled public pair uses
the current names. Testing against the grid there as well would let a file
and its archive agree on the old form too. It is a larger change and was
not needed to stop the wrong value.

**Why a missing unit is still micrometres.** Three reasons.

1. The imzML path reads it so. A file and the archive made from it should
   give the same pixel size, and the converter adds no unit of its own.
2. For one number, the test is also a test of the unit. The extent carries
   micrometre in every old-form file that could be tested, so a value that
   passes is a value in micrometres.
3. For a declared pair, the public files agree where they can be tested.
   See the table below.

**Why centimetre is refused.** 31 public files declare the unit accession
`UO:0000015`, centimetre, and name it "micrometer" in the same parameter. A
factor of 10000 would be wrong for every one of them. Reading the number as
written would trust a unit the file contradicts. An archive keeps the
accession only, so Thyra cannot see the contradiction. The value is refused
and the log names the unit.

**What was measured.** Headers of public imzML files, read on 2026-09-29.
479 of them give a pixel size.

| Files | What they give | Result |
|---|---|---|
| 39 | `IMS:1000046` alone, named "pixel size", no unit, with count and extent | sqrt(value) x count is the extent exactly, in all 39. None is a length. Values 10000 and 100: pixels of 100 um and 10 um |
| 31 | `IMS:1000046` alone, named "pixel size", unit `UO:0000015`, count but no extent | Cannot be tested |
| 409 | Both terms | One file read per deposit, writer, unit and value (per contributing group for HuBMAP): 28 files that stand for 408. All 28 use the current names. One file could not be fetched |

Four of those 28 files give the pair with no unit. They stand for 28 files
in four MetaboLights deposits.

| Deposit | Value | Test |
|---|---|---|
| MTBLS2639 | 25 | A length, against an extent in micrometres |
| MTBLS12782 | 50 | A length, against an extent in micrometres |
| MTBLS12204 | 1.0 | The extent equals the count and has no unit. Says nothing about the unit |
| MTBLS2075 | 40 | No count and no extent |

The reference converter (HUPO-PSI/mzPeak at `bb0f307`) was run on an imzML
of the old form. The archive holds `"name": "pixel size", "value": 10000,
"unit": null`. Thyra 4.2.0 stored 10000 by 10000 from it. This branch stores
100 by 100.

**The objection.** Refuse every pixel size that has no unit, or that the
grid does not confirm, pairs included. That is the most cautious rule. It
did not win because it would send every file of a whole writer family to
`--pixel-size` for a number the file states twice, and it would make the
archive stricter than the imzML it was made from. The error that was found
is in one number read two ways, and that is what the test covers.

**Known limits.**

- One public file writes its extent as (count - 1) x size. A lone length in
  such a file fails the test and the conversion asks for `--pixel-size`.
- A pixel that is not square cannot be recovered from an area. Its two axes
  disagree in the test, and no pixel size is taken.
- A value within 2% of 1 passes both readings. They give the same side.
- The pixel size of MTBLS12204 is 1.0 with no unit, and is stored as 1 um on
  both paths. Nothing in the file says whether that is a size or a
  placeholder.
- The two paths differ on the old name. The mzPeak path tests the value, and
  the imzML path refuses it. An old-form file with a count and an extent
  needs `--pixel-size` as imzML and none as an archive.
- A writer that keeps the old name and means a length is refused on the
  imzML path. No such file was found.

## D28. The chunked mzPeak layout is decoded into the same spectra

**Status:** Implemented (2026-09-29).

**Decision.** The mzPeak reader reads the chunked layout. It decodes the
chunks into the stream of spectra the point layout gives, so nothing after
the reader knows which layout the archive had.

- **Four chunk encodings are decoded**: no compression (`MS:1000576`), delta
  encoding (`MS:1003089`), MS-Numpress linear prediction (`MS:1002312`) and
  grid encoding (`MS:1003826`) with the linear (`MS:1003824`) or the square
  root (`MS:1003825`) model.
- **The rest is refused by its CV term**, before a spectrum is read: any
  other chunk encoding, any other grid model, intensities stored under a
  transform, and chunks cut along another axis than m/z.
- **The intensity list says how many points a chunk holds.** Every decoder
  is held to that count, and a chunk that decodes to another is refused.
- **The first and last m/z of a chunk are its bounds.** `mz_chunk_start` and
  `mz_chunk_end` replace the two decoded ends where they lie within 1e-5 of
  them.
- **No shared mass axis is reported**, for the grid encoding either.
- **The store names the encodings** it was read from, in
  `format_specific.chunk_encodings`.

The member that is read is chosen as before: `data_arrays` when it holds
rows, else `peaks`. The layout is that member's.

**Why read the layout.** `mzpeak-convert` writes it unless told otherwise
(0.14.0 and 0.16.0 were checked). All eight public example archives use it,
under MS-Numpress; they were written by its 0.12.0. Asking the converter
for the point layout is not a way around: for centroid m/z on a fixed-point
lattice it still writes a chunked `peaks` member, beside a `data_arrays`
member in the point layout.

**Why four encodings at once.** The converter chooses the encoding from the
values, so one encoding would read one kind of input.

| Input | Encoding written by default |
|---|---|
| m/z that are ordinary floating point numbers (two real images, one profile and one centroid) | MS-Numpress linear, every chunk |
| profile m/z on a fixed-point lattice | delta encoding |
| centroid m/z on a fixed-point lattice | grid encoding, linear model |

Delta encoding alone would have left out both real images.

**Why MS-Numpress is decoded here.** The reference reader decodes it with
`pynumpress`, a compiled package. The byte format is short and published
(Teleman et al., Mol Cell Proteomics 2014, 13, 1537). Decoding it in numpy
costs speed, and it spares every install a compiled dependency for one
encoding of an experimental format.

A residual in the buffer is a half byte that says how many half bytes
follow, so the place of one is known only from the one before. The decoder
gives every half byte the place the next residual would have, and finds the
places that are reached by doubling the stride. The values are whole
numbers until the last division, so two running sums give them exactly.

**Why the intensity list decides the count.** The list is plain in every
encoding, so its length needs no decoding. Delta encoding needs the count.
A chunk normally leaves its first value to `mz_chunk_start`. A chunk that
begins with a null keeps the null in the list, and then the list is one
longer. The reference decoder tells the two apart by looking at the second
value, and reads a chunk that opens with a null pair one value too long.
The count gives the right answer in both cases.

**Why the bounds replace the decoded ends.** A lossy encoding returns the
first peak of a spectrum a little off. The archive declares the m/z range
of each spectrum from the exact values, and the resampled axis is built on
that range. A first peak decoded 1e-7 below the axis was left out of the
store. On a grid-encoded archive of 20 pixels, one pixel lost a peak of
intensity 512 that way before this rule.

The bounds are the source's own numbers. In 1,127,540 chunks of five
archives the first and last m/z of the point layout were equal to the two
bounds, to the bit. In all 45,525 spectra the declared lowest and highest
m/z were equal to the first and last bound.

**Why no shared axis.** A shared axis means every spectrum holds the same
m/z values. The grid encoding gives each chunk a model, and a spectrum
still lists the grid indices it holds. The converter fits one model per
spectrum:

| Input, 20 spectra | Grid models in the archive |
|---|---|
| processed centroid imzML | 20 |
| continuous centroid imzML | 1 |

One model for the whole archive came about only where the spectra held the
same m/z already. To know that they do, every index has to be read. That is
the pass a shared axis exists to save, so nothing would be gained. The
default conversion resamples and takes its axis from the declared range,
with no pass over the data in either case.

**What was measured.** `mzpeak-convert` 0.14.0 (commit `0ed311e`) converted
each input twice, in the chunked and in the point layout, on 2026-09-29.
Grid-encoded inputs were given `--no-mz-lattice` for the point layout.
Both archives were read and compared point by point.

| Input | Encoding | Points | Largest m/z difference |
|---|---|---|---|
| centroid image, 12,737 pixels | `MS:1002312` | 60,149,625 | 2.3e-7, or 0.0002 ppm |
| profile image, 16,384 pixels | `MS:1002312` | 36,371,295 | 6.6e-7, or 0.005 ppm |
| the same, with `--no-numpress` | `MS:1003089` | 36,371,295 | 7.1e-15 |
| synthetic profile, lattice m/z | `MS:1003089` | 5,048 | 0 |
| synthetic centroid, lattice m/z | `MS:1003826` | 5,594 | 9.4e-8 |
| synthetic continuous centroid | `MS:1003826` | 8,000 | 9.4e-8 |

Pixels, order, point counts and intensities were the same in all six.

The two real images were measured again on 2026-09-30, through the reader
as it now decodes, with the ends of each chunk taken from its bounds. The
largest relative difference is given beside the largest in Da. An earlier
version of this page gave 8.0e-7 and 0.11 ppm for the profile image.

None of these archives holds padding. The reference writer's own sample
does, in both layouts (HUPO-PSI/mzPeak at `bb0f307`, delta encoding). Its
chunked member decodes to 217,710 points, 39,968 of them padding. The
padding sits at the same points as in the point layout, and every m/z is
equal. The sample is not an image, so only its signal member was read.

The stores, converted with the command line's defaults:

| Input | Encoding | Store against the point layout's |
|---|---|---|
| profile image | delta | every array the same |
| synthetic centroid | grid | every array the same |
| centroid image | MS-Numpress | 287 of 60,149,625 entries in a neighbouring bin |
| profile image | MS-Numpress | one entry more among 24,267,357 |

In the two MS-Numpress stores the total ion current image was the same in
every pixel. The average spectrum differed in 10 of 460,495 channels and
in 34 of 708,847.

Three archives in the point layout convert to the same stores as before
this change.

Constructed archives cover the rest. A chunked archive and its point twin
convert to identical stores in each of the four encodings, with padding,
and for the centroid member.

Conversion time, chunked against point: 116 s against 21 s for the centroid
image, 74 s against 15 s for the profile image, 21 s against 16 s with delta
encoding. MS-Numpress decodes at about 2 million points per second.

**The objection.** Refuse the lossy encodings: a store should hold the
numbers of the instrument. It did not win because the loss is made when the
archive is written. Refusing the archive does not bring the numbers back, and
it leaves the default output of the converter unread. The largest error
measured is 0.005 ppm, and the store says which encoding it came from.

**Known limits.**

- The archives in the tables were made here. The eight published example
  archives were read later and compared with their source imzML files; the
  results are on issue #422.
- Without resampling, a lossy archive gives a mass axis far longer than its
  point twin's. The centroid image has 331,701 distinct m/z in the point
  layout and 48,116,750 under MS-Numpress. The profile image has 1,122,721,
  against 1,147,161 under delta encoding and 18,731,686 under MS-Numpress.
  The reader warns of it. Resampling, the default, is not affected.
- The square root model and the encoding without compression were read from
  constructed archives only. The converter wrote neither for these inputs.
- A writer may give the edges of the interval it cut at as the bounds of a
  chunk. Such bounds lie far from the decoded ends and are not used, so a
  lossy archive of that kind can still lose a peak at the edge of the axis.
  No such archive was found.
- Vendor grid models are refused. The converter writes them for timsTOF
  input.
- Archives of the older prototype, mzML2mzPeak, are not read.

---

## D29. An mzPeak archive converts to the store of its source

**Status:** Implemented (2026-09-30), issue #422.

**Decision.** Where an archive states the same fact as the file it was made
from, Thyra reads it the way it reads that file. Five facts:

- **Where the image sits.** Positions are rebased on the base the archive
  declares in `imaging.coordinate_base` (1 when it declares none), or on the
  smallest position when that is lower. That is D14's rule for imzML, so a
  cropped image keeps its place. A writer's shift, `imaging.position_offset`,
  is added back into `coordinate_offsets_px`.
- **The instrument.** The instrument configurations are resolved by the
  function the imzML extractor uses. The analyzer or model picks the mass
  axis, and `msi_metadata` names the model.
- **Regions.** The regions under `bruker_maldi.regions` are read when their
  boxes place every pixel exactly once and each region holds as many pixels
  as it lists frames. Otherwise the log says why and the store has one
  region.
- **MS levels.** When MS1 spectra sit on the pixels, spectra of level 2 and
  up are left out, and the log counts them. Level 0, which mzpeak-convert
  writes when the source states no level, is kept. A spectrum left out adds
  nothing to the mass axis or the mass range.
- **Embedded images** are carried into the store, unaligned.

**Why.** The eight public example archives were converted beside their
source imzML files, and one archive from a Bruker TSF run beside the `.d`.
The spectra agreed; the stores did not.

| Example | Before | After, and from its source |
|---|---|---|
| glioma, positions x 675 to 735 | 61 x 83 grid, x from 0 | 735 x 259 grid, x from 674 |
| chilli | 85 x 50 grid | 92 x 70 grid |
| Bruker TSF run | offset (1, 1) | offset (669, 109) |
| mouse bladder, LTQ Orbitrap | instrument unknown, constant axis | Orbitrap axis |
| DESI, Exactive | instrument unknown, reflector TOF axis | Orbitrap axis, 504,040 bins |

After the change, all six public imaging datasets land on the grid of their
imzML, pixel for pixel. The Bruker archive gives the grid and offset of the
`.d`: 169 x 200 at (669, 109).

**The objection.** An archive is its own format. Rebasing on the smallest
position, as before, gives the smallest grid, and a reader of mzPeak owes
imzML nothing. It did not win because the archive states the same base as
the imzML, and the imaging profile's own converter writes it. A store that
changes shape with the container turns one acquisition into two, and every
pixel coordinate a user wrote down against one no longer fits the other.

**Why regions are held to their boxes.** No scan names its region; the
archive lists each region's box of raster indices and its frame count. A box
says where a region is, not that every pixel in it belongs to it. Two
irregular regions can share a box. The count check is what makes a box an
answer.

**Why MS2 is left out.** A pixel's spectra are summed into one row. An MS2
spectrum summed into its pixel's MS1 spectrum is a spectrum the instrument
never measured. A store made from mzPeak has no table for MS2 spectra, so
leaving them out loses nothing that the store could hold.

**Why the images are not aligned.** Each image comes with an affine from
image pixels to MS pixels. In every public example it is marked
`assumed_full_extent`: the image stretched over the whole acquisition. That
is a guess about the frame, not a registration. The affine stays in the
store's raw metadata.

**Known limits.**

- No archive with two regions exists. The region rule is tested on
  constructed archives.
- The Waters and timsTOF TDF lanes of mzpeak-convert 0.16.0 were not tried.
- The chilli archive declares profile spectra and holds them in the centroid
  member. It is read as centroid, since the member decides (see
  [Supported Formats](supported-formats.md#mzpeak-experimental)); its imzML is
  read as profile.

---

## D30. A z the source does not state is recorded as 0

**Status:** Implemented (2026-10-01), issue #425.

**Decision.** `coordinate_offsets_px` gives the source's own coordinate of
index 0 on each axis. On z that is the smallest z the source states, and 0
when it states none. Every vendor reader that records offsets already
recorded 0 on z, its format having none. Two routes change:

- **imzML.** pyimzml gives a spectrum without `IMS:1000052` the z of 1. That
  1 is still subtracted, so the spectra sit on plane 0, but the store
  records 0. A z the file states is recorded as before (D14). When every z
  is 1, the first spectrum is read again to tell a stated 1 from pyimzml's.
- **mzPeak.** The reader reads position z, which mzpeak-convert writes when
  its source states z. One plane records its z. Spectra on several planes
  are refused, and so are spectra that state z beside spectra that do not:
  Thyra reads one plane of an mzPeak archive. Before, the column went unread
  and every plane became one.

**Why.** The public examples gave the same 2D data two z offsets: 1 from the
imzML and 0 from the archive. The 1 was pyimzml's, not the file's. The store
defines the field as the offsets that normalisation erases, and here it
erased nothing the file holds.

| Example | z stated by | Before: imzML / archive | After: imzML / archive |
|---|---|---|---|
| DESI, both 3 x 3 examples | neither | 1 / 0 | 0 / 0 |
| glioma, archive from 0.16.0 | both | 1 / 0 | 1 / 1 |
| glioma, archive from 0.12.0 | the imzML only | 1 / 0 | 1 / 0 |
| Bruker TSF run, archive from 0.16.0 | neither (`.d` / archive) | 0 / 0 | 0 / 0 |

**The objection.** Leave imzML alone, and record 1 for an archive without z:
the base mzPeak declares for its positions. Only mzPeak stores would change.
It did not win: the archive made from a Bruker run would then record 1
where its `.d` records 0, and that 1 would again be a value nobody stated.

**Known limits.**

- The glioma archive from 0.12.0 still records 0 against its imzML's 1. That
  converter dropped the z column; 0.16.0 writes it.
- An imzML that cannot be read a second time keeps the subtracted 1.
- No mzPeak archive with more than one plane was found, so the refusal is
  tested on constructed archives.

## D31. FlexImaging spots are placed on their lattice, not on the drawn Area

**Status:** Implemented (2026-10-05), issue #428.

**Decision.** A Bruker MALDI run is placed on its FlexImaging photo by one
affine. The teaching points map the photo to the stage. The spots sit on a
lattice one raster step apart, through the `.mis` reference point. The
raster node on that point is the one that puts every measured spot inside
its own Area, and of those the one that leaves the fewest unmeasured nodes
inside the Areas. The same affine places the TIC image, the pixel polygons,
the crop window and `raster_to_global_affine`.

A TIC cell `i` spans `[i, i + 1)`, as SpatialData draws it, so its centre
`i + 0.5` lands on the spot. A pixel polygon is that cell exactly: a
parallelogram when the photo is rotated against the stage.

**Why.** Thyra stretched each region over the bounding box of its Area
outline. The outline is drawn by hand and only selects nodes, so pixels came
out 0.87 to 1.14 raster steps wide. The TIC image was also half a cell off
its own polygons.

| Run (MassIVE MSV000088438) | Before: polygons / TIC, max | After: both, max |
|---|---|---|
| TSF, 1000 um, 276 spots | 0.57 / 1.06 steps | 0.0010 steps |
| TDF, 1000 um, 240 spots | 0.52 / 0.51 steps | 0.0014 steps |

The reference is flexImaging's own spot list, mapped through the teaching
points. On two local runs at 20 um, without a spot list, the commanded
positions in the poslog fit the lattice within 0.09 um after one
translation. Both runs need the same translation, so no node is a step off.

**What changes in a store.** Every store aligned to a FlexImaging photo: the
polygons, the TIC transform, `raster_to_global_affine` and the crop window.
With `apply_optical_alignment=False` the photo still maps a spot onto
`obs["spatial_x"]`, now from the lattice.

**The objection.** Keep the Area stretch for a run with one region. The
outline usually hugs the spots, it needs no teaching points, and the error is
a fraction of a step. It did not win: on a one-region local run at 20 um the
polygons were still up to 0.46 steps off and the TIC up to 0.99. A fraction
of a step is the scale MSI images are read at.

**Known limits.**

- Fewer than three teaching points, no raster step, or measured Areas with
  different steps fall back to the bounding-box stretch, with a warning.
- No `<ReferencePoint>`: the first teaching point is used. Every `.mis` seen
  puts the reference point there.
- In micrometre stores the TIC image and the polygons were still half a
  pixel apart. That was not FlexImaging-specific and was fixed by
  [D32](#d32-a-tic-cell-is-centred-on-its-pixels-position).

## D32. A TIC cell is centred on its pixel's position

**Status:** Implemented (2026-10-05), issue #431.

**Decision.** In a micrometre store, pixel `x` sits at
`obs["spatial_x"] = x * pixel_size_um`. Its polygon is centred there, and
now so is its TIC cell. SpatialData draws cell `x` over `[x, x + 1)`, so the
TIC transform moves each cell back half a step, then scales. A volume does
the same in z. `raster_to_global_affine` carries the same half-step
translation.

Both kinds of store now follow one rule: a TIC cell's centre is where its
pixel is. On a FlexImaging photo that is the spot (D31). In micrometres it
is the `obs` position.

**Why.** The TIC transform was a plain `Scale(pixel_size_um)`. Cell `x` was
centred on `(x + 0.5) * pixel_size_um`, half a pixel right of and below its
own polygon. The guard test allowed a whole pixel, and its comment called the
shift padding.

Each TIC cell's centre against its row, in raster steps, the worst row:

| Store | Pixels | Before: polygon / `obs` | After |
|---|---|---|---|
| DESI mzPeak example, 100 um | 17,952 | 0.5 / 0.5 | 0 / 0 |
| Waters kidney, 50 um | 5,400 | 0.5 / 0.5 | 0 / 0 |
| PHI nanoTOF, 512 x 512 | 262,142 | 0.5 / 0.5 | 0 / 0 |
| MSV000088438 TSF, alignment declined | 276 | 0.5 / 0.5 | 0 / 0 |

In the last store the photo was half a step off the TIC as well, and now
sits on it within 1e-14 steps. Between the two versions only the TIC
transform and `raster_to_global_affine` differ: every array, table and
polygon is byte-identical. The same run aligned to its photo is unchanged.

**What changes in a store.** In every micrometre store, the TIC transform
becomes a translation, then a scale. `raster_to_global_affine` gains minus
half a pixel in its third column. The polygons, `obs`, `stage_offset_um` and
every intensity stay as they were. Stores on a FlexImaging photo do not
change. With `apply_optical_alignment=False`, the photo now sits on the TIC
as well as on the polygons.

**The objection.** Move `obs` and the polygons instead, to
`(x + 0.5) * pixel_size_um`. The TIC transform would stay a pure scale, and
a tool that looks for a scale-only raster as its reference would still find
it. It did not win: every `obs` position would change, and so would the
meaning of `stage_offset_um`. A source position is where a spectrum was
taken, which is the pixel's centre, and `x * pixel_size_um` already says
that.

**Known limits.**

- The raster now starts at minus half a pixel, not at 0. A consumer that
  applied the TIC matrix to the integer index `x` now gets the cell's corner.
  Apply it to `x + 0.5` for the centre.

## D33. An mzPeak archive of a Bruker run converts to the store of its `.d`

**Status:** Implemented (2026-10-06), issue #429. Extends D29, and replaces
its "embedded images are carried into the store, unaligned" for an image
fitted to the teaching points.

**Decision.** `mzpeak-convert` 0.17 states five more facts of a Bruker run.
Thyra reads each the way it reads the `.d`.

- **m/z on the timsTOF grid.** The default TDF lane stores each point as a
  TOF bin and the frame's calibration, under the grid model `MS:9999002`.
  Thyra solves the bin through that calibration with the equation of the
  reference implementation (`mzdata`, `MzCalibrationModel2`).
- **TDF intensities.** The archive holds each point's raw count. Bruker's
  library returns `floor(count x 100 / accumulation_ms + 0.5)` for every
  point of a scan, and the `.d` reader sums the scans that fall on one TOF
  bin. Thyra does both, with the accumulation time the archive gives as each
  scan's ion injection time. The store says so in
  `format_specific.intensity_scale`.
- **The mass range** is `MzAcqRangeLower` to `MzAcqRangeUpper` from
  `vendor_metadata`, the range the `.d` takes, not the observed range.
- **Regions.** Each scan names its region in the parameter
  `acquisition region`. That name is read; the boxes are the fallback for an
  older archive.
- **The alignment image.** An image whose affine is marked
  `registration_quality: teach_points` maps image pixels to pixel positions,
  centre to centre. Its inverse, after the shift from a position to a TIC
  cell, is the lattice of D31. The TIC, the polygons and
  `raster_to_global_affine` are placed by it, as for the `.d`. The
  teaching points and the image name the archive carries from the `.mis`
  fill the metadata document's `alignment` section, as the `.d`'s do. The
  pixels are placed on the image only when the store carries it.

**Why.** Two public timsTOF fleX runs of MassIVE MSV000088438 (one TSF, one
TDF) were converted from the `.d` and from the archives `mzpeak-convert`
0.17.2 wrote of them.

| Archive | 4.3.4, against the `.d` | Now, against the `.d` |
|---|---|---|
| TDF, default | refused | every value the same: 0 of 1,352,283 entries differ |
| TDF, `--no-ims-compact` | intensities x 1.9995 (accumulation 199.953 ms) | 58 entries one bin over; every pixel's total the same |
| TSF | axis of 406,413 bins (m/z 254-1939) | the `.d`'s axis, 599,146 bins (100-2000) |
| Both, regions | read from boxes | read from the scans; the same on every pixel |
| Both, pixels on the photo | not aligned | within 0.004 photo pixels of the `.d` route; one pixel is 123 to 125 |
| Both, metadata `alignment` section | absent | the same as the `.d`'s: image, method, three teaching points |

The decoded timsTOF m/z were checked against the ends of all 8,735 chunks of
the TDF archive, which the converter evaluated itself: within 7e-16,
relative.

The rounding was measured on a copy of the TDF run with its accumulation
time set to 200 and to 40 ms. Then half of all points fall exactly on .5.
Over 439,522 points, half up matched Bruker's library on every one, and
round-half-to-even on 76%. The archive keeps the time in single precision.
On this run that changed none of the 1,784,948 rounded values.

**The objection.** Keep the raw counts: they are what the instrument stored,
and an archive is a record. It did not win. D29 makes a store independent
of the route, and the two stores of one run differed by accumulation_ms /
100, 0.1 on a 10 ms run. The archive gives the time, so the `.d`'s values
can be had exactly.

**Known limits.**

- No archive declares its intensity convention. Thyra takes a TDF source
  (`MS:1002817`) as raw counts. A writer that stored the library's values
  would be scaled twice. A TDF archive without every accumulation time keeps
  its raw counts, with a warning.
- A TSF archive's m/z are 3.6 to 5.0 ppm below the `.d`'s on its six
  strongest peaks. The archive stores m/z as numbers, from its converter's
  own TSF calibration, and carries no vendor calibration to correct them. So
  its entries sit in neighbouring bins; every pixel's total is the same.
- Only the image marked `teach_points` is aligned. A FlexImaging run's other
  images are not in the archive.

## D34. The intensity threshold tests a peak summed over its mobility scans, in every table

**Status:** Implemented (2026-10-07).

**Decision.** `--intensity-threshold` keeps a peak when its intensity, summed
over every mobility scan of its pixel, reaches the threshold. Every table
applies that one test:

- **timsTOF `.d`.** A digitizer index is kept when its sum over the whole
  ramp passes. The summed table, the `--mobility-grid` table, the
  mass-mobility heatmap and the MS/MS table keep or drop all of that index's
  points together, under either `--tdf-spectrum`.
- **imzML and mzPeak.** A spectrum that lists one point per mobility scan
  repeats a peak's m/z. Points whose m/z agree to within 1e-12, relative, are
  one peak, tested by their sum. A spectrum that lists each m/z once is
  tested point by point, as before.
- **The store says so.** A store converted with a threshold records it as
  `intensity_threshold` in the `conversion` step of
  `msi_metadata.processing`, from any source.

**Why.** Before, each table tested a different quantity. The summed table
tested the ramp sum, the grid and the heatmap each raw `(index, scan)`
point, and the MS/MS table each window sum. Intensities are never negative,
so the other tables could only lose current the summed table kept. D1 and
D3 state that they add up exactly under `scan_sum`, and with a threshold
they did not. The store did not record the threshold either.

| Store, `--mobility-grid` | Threshold | Grid | Heatmap | MS/MS | Before |
|---|---|---|---|---|---|
| `synthetic_tims` with a PASEF schedule | 50 | 1.0000 | 1.0000 | 1.0000 | 0.8427 / 0.8427 / 0.9862 |
| `synthetic_tims` | 300 | 1.0000 | 1.0000 | - | no grid table, no heatmap |
| `tims_msms_pos_brain1` | 100 | 1.0000 | 1.0000 | 1.0000 | 0.9869 / 0.9869 / 0.9997 |

The ratios are each table's ion current over the summed table's. On
`tims_msms_pos_brain1` the summed table is byte-identical to before; the
grid gains 710 entries and the MS/MS table 55. At 300 every raw point of
`synthetic_tims` is below the threshold and every summed peak above it, so
the grid used to come out empty at exit 0.

`mobility_continuous.imzML` at threshold 2.5 now holds 11 at m/z 300 and 22
at m/z 600.25 in pixel 0, as without a threshold; before it held 10 and 20.
`synthetic_tims.d` and its three `mzpeak-convert` 0.17.2 archives (default,
`--no-ims-compact`, and that with `--layout point`) give the same 3471 per
pixel at thresholds 100 and 300. Those archives already did on main,
because D33 sums a frame's points per TOF bin before the threshold.

Stores made without a threshold did not change: the store-identity harness
found 0 differences on the eight default datasets and on
`tims_msms_pos_brain1`. At threshold 1e-9 on five sources without ion
mobility, the only difference is the new `processing` parameter.

**The objection.** Test each single scan point everywhere, the summed table
too: one simple rule, and weak single points leave the grid. It did not
win. A peak far above the threshold once summed could vanish, thresholded
`.d` stores would lose current they keep today, and removing noise from the
grid is analysis, as D1 argues for the centroid.

**What changes in a store.** Stores converted with `--intensity-threshold`
from ion-mobility data: the grid, heatmap and MS/MS tables regain current,
and their values only rise. The summed table of a `.d` does not change. An
imzML or mzPeak spectrum that lists a peak once per scan keeps the points
of a peak whose sum passes, so its summed table rises to match. Every store
converted with a threshold gains the `intensity_threshold` parameter.
Stores made without a threshold do not change.

**Known limits.**

- Under `vendor_centroid` the summed table is Bruker's picked spectrum,
  thresholded per picked peak, while the other tables test raw index sums.
  These tables did not add up in that mode before, and still do not.
- Points of one peak are grouped by m/z within 1e-12, relative. An archive
  whose encoding spreads one bin's m/z further apart is tested per stored
  m/z. The MS-Numpress archives of `mzpeak-convert` 0.16.0 were seen to do
  so; none was measured for this entry.

## D35. Polarity is the value every statement in the source agrees on

**Status:** Implemented (2026-10-07).

**Decision.** `ms_analysis.polarity` is recorded when every place the
source states a scan polarity gives the same one. A source that states both
values anywhere records neither, and says so in a warning. The places read:

- **imzML:** the `fileContent` terms, every referenceable parameter group,
  and every spectrum's own terms (`MS:1000130` positive, `MS:1000129`
  negative). The parser collects the spectra's terms in the one pass it
  makes over them anyway.
- **mzPeak:** the terms in `file_description.contents`, and every
  spectrum's `scan_polarity` (1 positive, -1 negative; any other value
  states nothing).
- **timsTOF `.d`:** `Frames.Polarity`, unchanged.

`acquisition_params.polarity` holds the agreed value. It holds `mixed`
when the statements disagree, which no polarity term matches, so the file's
`fileContent` term cannot fill the field instead. A source that states
nothing gets no key.

**Why.** Writers state polarity in different places. pyimzml's writer
puts it in a parameter group; `pea` states it on each of its 12,737
spectra and nowhere else. Thyra read only `fileContent`, so both converted
with no polarity at exit 0, and `export-metaspace` then left a required
field empty. A file whose `fileContent` said positive while every spectrum
said negative was stored as positive.

An mzPeak archive was never read for polarity, so an archive converted
without the polarity of the file it was made from. D29 says the two give
the same store.

| Source | Before | Now |
|---|---|---|
| `pea.imzML` | none | positive |
| A Xenium-run imzML, negative on all 918,855 spectra | none | negative |
| `bellini.imzML` (`fileContent`) | positive | positive |
| `mzpeak-convert` 0.17.2 archives of two timsTOF fleX runs | none | positive, as their `.d` |
| `mzpeak-convert` 0.17.2 archives of `Example_Continuous`, `Example_Processed`, `desi_colad_centroid` | none | negative |
| `fileContent` positive, every spectrum negative | positive | none, with a warning |

Collecting the terms adds 15 MB of traced memory to the parse of the
918,855-spectrum file. Its time is within the noise: 37 to 65 s for either
version, over two runs each.

**The objection.** Take the majority, or the first spectrum, as pyimzml
does. It did not win. A file that alternates polarity has no single true
value, and the timsTOF rule already leaves it unset.

**What changes in a store.** Stores of imzML files that state polarity, and
of mzPeak archives that state one, gain `polarity` in `acquisition_params`
(and so the root attribute `acquisition_parameters`). Where only a
parameter group or the spectra stated it, `ms_analysis` gains `polarity` and
`polarity_term`. A file whose statements disagree loses the value it had.
Stores of sources that state no polarity do not change.

**Known limits.**

- A metadata-only read (`thyra metadata`, `read_metadata_document`, the
  preview) reads the imzML head and the first spectrum, not every spectrum.
  It can report a polarity that a conversion of a file whose later spectra
  switch would not record.
- Only the scan polarity terms are read. A writer that states polarity in
  a user parameter is not.

## D36. The imzML reader widens 32-bit arrays to float64 where it yields them

**Status:** Implemented (2026-10-07).

**Decision.** Every m/z, intensity and ion-mobility array the imzML reader
yields is float64, whatever precision the file stores. The widening happens
once, at the reader boundary: `_read_spectrum_arrays` (both the
`getspectrum` routes and the continuous fast read), the m/z-only reads in
`thyra/utils/pyimzml_direct.py`, and the points where the continuous shared
m/z block is cached. The block is cast once when cached, never per spectrum,
so the converter's identity check on the shared axis stays O(1). A file that
already stores 64-bit floats is not copied: `astype(float64, copy=False)`
returns the decoded array unchanged.

**Why.** Many imzML files store m/z and intensity as 32-bit floats. The
reader handed those arrays on in their storage dtype, and every downstream
sum and comparison ran in 32-bit. Four results came out wrong at that
precision, all from the same cause:

| Effect | Before | After |
|---|---|---|
| `--no-resample` TIC image vs its own row sums, a 32-bit processed file, 4 pixels | all 4 differ, max 1.8e-7 relative; `[2^24, 1, 1, 1]` stores TIC 16777216 beside a row sum of 16777219 | equal on every pixel |
| `tic_preserving` row total vs the pixel's measured total, a 32-bit continuous file | 2.9e-8 off, matching the 32-bit sum to 8e-15 | 9.3e-15, matching the 64-bit file |
| `--intensity-threshold 0.7` on a 32-bit processed file | keeps 0.699999988 (three pixels' worth), which the lossless mzPeak archive of the same values drops | keeps nothing below 0.7; the two routes agree |
| `--resample-min-mz 400.00001` on a 32-bit centroid file | the peak at 400.0 is kept on seven of eight pixels, and a valid file is refused with "the two passes over the source disagree" | the peak drops on every pixel; the file converts |

The mzPeak, solariX, timsTOF, Waters and PHI readers already yield float64,
which is why only the imzML route showed these. Their stores were exact;
imzML stores of the same data were not (D29 says the two give the same
store).

Cost, measured on a 918,855-spectrum Xenium-run export whose m/z array is
32-bit, default options, one subprocess per conversion, alternating sides,
two runs each, this machine:

| Side | Run 1 | Run 2 | Peak working set |
|---|---|---|---|
| Before | 366.5 s | 368.4 s | 16,704 / 16,828 MB |
| After | 363.7 s | 361.1 s | 16,826 / 16,728 MB |

The time difference is inside this machine's noise (identical code has
measured 37-65 s apart on a parse alone); the extra memory is about 10 MB,
0.07% of the peak, because the store's `X` was already float64 and the cast
replaces one in-flight array with a wider one, not a second copy.

**The objection.** Widen at each consumer instead: the two TIC sums in the
streaming converter, the threshold comparison, the resamplers' masks. It did
not win. The list of consumers is open-ended — every future sum or
comparison would have to remember the rule — and two of the four effects
(the two-pass refusal and the per-pixel disagreement) come from two code
paths disagreeing, which per-consumer fixes must keep in lockstep forever.
One cast at the boundary makes the reader's signature honest: its
annotations always said `NDArray[np.float64]`.

**What changes in a store.** Stores converted from an imzML file that
stores 32-bit floats:

- The TIC image values of a `--no-resample` store change by about 1e-7
  relative, to equal each pixel's row sum.
- `tic_preserving` rows are rescaled to the measured pixel total; stored
  values move by about 1e-7 relative.
- With `--intensity-threshold T`, values that are `float32(T)` but below
  `T` drop out; with a very small `T`, zero-intensity points drop.
- With a typed `--resample-min-mz`/`--resample-max-mz` bound, peaks stored
  at `float32(bound)` are now dropped when they lie outside the typed
  range. Ranges taken from the data are unaffected.

The stored m/z axis of a 32-bit file is float64 now; the values are
unchanged (widening is exact). Stores of 64-bit imzML files and of every
non-imzML source are unchanged.

**Known limits.**

- The arrays are wider, not more accurate: a value the file stored as
  32-bit stays the 32-bit-rounded number. Only arithmetic Thyra does on it
  improves.
- A threshold below about 7e-46 was 0 once cast to float32 and kept every
  zero-intensity point; it now drops them, as the mzPeak route always did.
- `read_metadata_document` and the preview read the imzML head and early
  spectra; the mass range they report was always computed from the stored
  values, which do not change, so metadata documents are unaffected.
