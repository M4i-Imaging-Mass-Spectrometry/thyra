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

**Measured again** (2026-10-12; centroid total against `Frames.SummedIntensities`):
0.813 (0.790 to 0.826 per frame) over 8 frames of a 1,465-pixel MALDI-2 slide,
0.98 to 0.99 on two MS/MS brain images, and 0.895 on survey frames and 0.971
on DIA frames of an LC diaPASEF run. Across the imaging acquisitions measured
the centroid keeps 79 to 99 percent of the current, not the 87 to 96 percent
first written down.

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
`DiaFrameMsMsWindows` is not one of the two tables the TDF reader looks in
(D37 reads it now),
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
`BaseMSIConverter.options` without a word (through `convert_msi` they are
refused since D76). Dropping `sparse_format` from the signature and stopping
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

**Status:** Implemented (2026-09-09), issue #244. Amended by D65: the
preview reads the declared pixel counts, and the base is logged once.

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

**Status:** Superseded by D49 (2026-10-14): releases go through a pull
request, so no release step pushes to `main` and required checks need no
bypass. Previously deferred (2026-09-09), issue #224. Closed as not planned
(2026-10-13). The organization-level ruleset route is ruled out: it needs a
paid plan, and the organization is on the Free plan. Probed with an
`admin:org` token, the org rulesets endpoint returns `HTTP 403`, "Upgrade to
GitHub Team to enable this feature". The routes left both need a new
credential. What reopens it, replacing the conditions under **What would
reopen it** below: several people merging regularly. Then use a GitHub App
owned by the organization as the ruleset's bypass actor, with `release.yml`
pushing with a token minted from it.

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

**Later.** D54 adds pyright to the same job, next to mypy.

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
else. `x`, `y` and `region_number` are plain integers too, and a plain integer
is stored as one flat array, where anndata stores a nullable integer as values
plus a mask that any reader of the Zarr without anndata has to know about. The
price: concatenating a table that has the column with one that has not gives
pandas' `float64` with NaN, unless the caller asks for `Int64`.

*Corrected 2026-10-16 (D57):* this entry first said `x`, `y` and
`region_number` are `int64`. In a 2D store all three are `int32`; in a volume
`x` and `y` are `int64` and `region_number` is `int32`. The reasoning holds:
each is a plain integer column.

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
- A value from 0.9606 to 1.0410, that is (1 - 1%)^4 to (1 - 1%)^-4, can pass
  both readings. The length is taken then, and the two sides differ by up to
  about 2%. (Corrected 2026-10-20, #523: this line said "within 2% of 1" and
  "the same side".)
- The pixel size of MTBLS12204 is 1.0 with no unit, and is stored as 1 um on
  both paths. Nothing in the file says whether that is a size or a
  placeholder.
- The two paths differ on the old name. The mzPeak path tests the value, and
  the imzML path refuses it. An old-form file with a count and an extent
  needs `--pixel-size` as imzML and none as an archive.
- A writer that keeps the old name and means a length is refused on the
  imzML path. No such file was found.

## D28. The chunked mzPeak layout is decoded into the same spectra

**Status:** Implemented (2026-09-29). Amended by D70.

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
| centroid image | MS-Numpress | 6 of 60,149,625 entries in a neighbouring bin, each one bin away |
| profile image | MS-Numpress | one entry more among 24,267,357 |

In the two MS-Numpress stores the total ion current image was the same in
every pixel. The average spectrum differed in 10 of 460,495 channels and
in 34 of 708,847.

(Corrected 2026-10-20, #523: the centroid row said 287 entries. 287 counts
the positions at which the stored `indices` arrays differ; 6 entries moved,
each by one bin.)

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
  The reader warns that it can happen: the archive does not say how many m/z
  its source held, so the growth itself is not measured. Resampling, the
  default, is not affected.
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

**Status:** Implemented (2026-09-30), issue #422. Amended by D69
(2026-10-19): an embedded image with no alignment sits only in a
coordinate system of its own. Amended by D70. Amended by D81 (2026-10-22):
the instrument comes from the configuration the converted spectra use, and
the spectrum type follows the converted spectra on both routes, which lifts
the chilli known limit below.

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
The spectra agreed; the stores did not. (Corrected 2026-10-20, #523: for the
Bruker TSF archive, positions, peak counts and intensities agree, but its m/z
sit several ppm off the `.d`'s in either calibration state. D33's known limits
give the measurement.)

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

**Status:** Implemented (2026-10-01), issue #425. Amended by D63: an
imzML's z is judged on every spectrum, not read again from the first.

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
  Corrected by D39: that converter states z in each scan's parameters list,
  which the reader now reads.
- An imzML that cannot be read a second time keeps the subtracted 1.
- No mzPeak archive with more than one plane was found, so the refusal is
  tested on constructed archives.

## D31. FlexImaging spots are placed on their lattice, not on the drawn Area

**Status:** Implemented (2026-10-05), issue #428. Amended by
[D59](#d59-a-fleximaging-run-that-stopped-early-is-placed-by-its-planned-grid)
for a run that does not fill its Area, and by
[D75](#d75-the-photo-crop-follows-the-regions-converted-and-one-affine-places-every-region)
for the fallback, which now places every region by one affine.

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
fitted to the teaching points. Amended by D61 (2026-10-16) and D69
(2026-10-19).

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

**Status:** Implemented (2026-10-07). Amended by D55 (2026-10-16) and D82 (2026-10-22).

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

**Status:** Implemented (2026-10-07). Amended by D55 (2026-10-16).

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

## D36. The imzML reader yields float64, whatever precision the file stores

**Status:** Implemented (2026-10-07).

**Decision.** Every m/z and intensity array the imzML reader hands on is
float64: each spectrum on the processed and the continuous path, the common
m/z axis, the mobility point cloud and the shared mobility features. A file
that stores 32-bit floats (or integers) is widened once, where the reader
reads it. Widening is exact; the stored numbers are the file's numbers. The
shared m/z block of a continuous file is widened once, when it is cached, so
every spectrum still yields that one object.

**Why.** pyimzml decodes an array in the precision the file declares, and
the reader passed it on unchanged. Everything after it computed in that
type: the TIC sum, the intensity threshold, the typed m/z bounds and the
total `tic_preserving` keeps. The other readers, and the mzPeak archive of
the same file (D29), yield float64. So a 32-bit imzML converted at 32-bit
precision, and four results came out wrong at exit 0:

| Small 32-bit imzML | Before | Now |
|---|---|---|
| `--no-resample`, intensities `[2^24, 1, 1, 1]` | TIC 16777216, row sum 16777219 | both 16777219 |
| `tic_preserving`, 3 pixels of 30,000 profile points | rows miss the measured total by 2.9e-8 | 9.3e-15, as the 64-bit twin |
| `--intensity-threshold 0.7`, value `float32(0.7)` | kept (0.699999988 is below 0.7) | dropped, as the mzPeak archive drops it |
| `nearest_neighbor`, typed minimum 400.00001, peak at m/z 400.0 | kept in 7 of 8 pixels | dropped from all 8 |
| the same, with the peak in the first pixel | refused: "the two passes over the source disagree" | converts |
| `tic_preserving`, four pixels whose only peak is at 400.0 | 5 rows | 1 row, as the 64-bit twin |

`pea` and a 918,855-spectrum Xenium-run export store 32-bit intensities,
but as whole counts whose pixel totals stay below 2^24. A 32-bit sum holds
those exactly, so their TIC images and `tic_preserving` totals do not
change. A threshold set one part in a billion above one of `pea`'s stored
intensity levels does:

| `pea`, `--no-resample`, that threshold | Before | Now |
|---|---|---|
| Stored values | 57,564,204 | 55,140,652 |
| Of them, below the threshold | 2,423,552 | 0 |

No file with 32-bit m/z, or with fractional 32-bit intensities, was at
hand. Those cases are measured on the constructed files above only.

**The objection.** Leave the reader as it is and compute in float64 where
it matters: the TIC sums, the threshold test, the two resamplers' masks.
It did not win. That is at least five places, each easy to miss: changing
the TIC sum alone makes every such conversion refuse, because the
second-pass check then disagrees with it. One widening where the arrays
enter covers them all, and any step added later.

**What changes in a store.** Stores of imzML files that store 32-bit
intensities or m/z:

- With `--no-resample`, the TIC image now equals each row's sum. A pixel of
  fractional intensities moves by up to about one part in ten million
  (1.8e-7 on the constructed file). Whole counts whose pixel total stays
  below 2^24 do not move.
- With `tic_preserving`, each row sums to the pixel's measured total.
- With `--intensity-threshold` at a value 32-bit floats cannot hold, a
  stored value just below it is dropped.
- With a typed `--resample-min-mz` or `--resample-max-mz`, a peak within
  half a 32-bit step outside the bound is treated as in the 64-bit twin.
  `nearest_neighbor` drops it from every pixel, and a file that was refused
  converts. `tic_preserving` drops a pixel whose only peak is outside; in a
  pixel that also has peaks inside, the peak's share is measured by area,
  as for a 64-bit file.
- With `--no-resample`, a file with 32-bit m/z stores `var["mz"]` as
  float64 instead of float32, with the same values, as every other reader.

A 32-bit file converted on the default route, without a threshold or a
typed bound, gives the same store as before. Stores of 64-bit imzML files and
of every other reader do not change.

**Known limits.**

- Each spectrum's arrays take twice the memory in flight. The store's `X`
  was already float64. On the 918,855-spectrum Xenium-run export (default
  route; old and new code alternated, two runs each) the time was 399 and
  363 s before, 388 and 365 s after. The peak working set was 16.35 GiB
  before and 16.44 GiB after; other runs of the old code gave 16.32 to
  16.44 GiB. That file stores 64-bit m/z, so only the intensity widening
  runs there.
- With `--no-resample`, a processed file with 32-bit m/z holds its native
  axis at twice the bytes. While the axis is widened both copies are held,
  so in the case where every m/z is distinct the axis build peaks near
  three times the m/z payload instead of two (61 to 91 MiB on a 38 MiB
  payload).
- The m/z-only read in `thyra/utils/pyimzml_direct.py`, and the axis
  builder it feeds, still return the file's precision. Each current caller
  widens the result or takes only its smallest and largest value, which are
  exact. Nothing stored differs, but a new caller must widen too.

## D37. Type-9 fragment frames take their isolation windows from the DIA tables

**Status:** Implemented (2026-10-12). Supersedes the statement in D2 that
`DiaFrameMsMsWindows` is not read.

**Decision.** The TDF reader reads the isolation windows of `MsMsType = 9`
frames. Each such frame names a window group in `DiaFrameMsMsInfo`, and
`DiaFrameMsMsWindows` lists the group's windows with their scan ranges,
read as `[ScanNumBegin, ScanNumEnd)` like a PASEF window. The schedule is
asked of `PasefFrameMsMsInfo`, then of the DIA tables, then of
`FrameMsMsInfo`; the first with rows is the one the frames use. Three
things stay as they were, on purpose:

- A file holding survey frames (`MsMsType = 0`) among its fragment frames
  still records its schedule as varying.
- A schedule that varies, or whose windows overlap, is never split by scan
  number. The converter already refused it; now `has_precursor_spectra` is
  false and `iter_precursor_spectra` refuses it too, so a direct caller
  cannot get one either.
- A raster position holding frames of more than one `MsMsType` is still
  summed into one pixel. The reader now says so at `WARNING`, naming the
  types.

**Why.** Bruker stores a scheduled precursor list in two layouts. PASEF
frames (type 8) list their windows per frame in `PasefFrameMsMsInfo`; the
two positive-mode MS/MS images Thyra was built on are of this kind, written
by a pre-release MALDI-prm build of the acquisition software. DIA frames
(type 9) use the DIA tables. diaPASEF writes them, and MALDI iprm-PASEF is
reported to; no such MALDI run has been measured. One LC diaPASEF run,
measured read-only: 2,755 survey and
22,039 type-9 frames, eight window groups of three windows each, the
groups taking turns frame by frame. Before this, every such file came out
as MS/MS with no windows at all.

Adjacent windows of one group share a boundary scan (`[76, 521)` and
`[521, 710)` on that run), which only a half-open range allows. The split
of a one-group type-9 schedule, through the vendor library on the
synthetic fixture, reproduces the in-window current of every frame exactly.

**The objection.** Wait for a real MALDI iprm-PASEF run before reading
anything. It did not win: the table layout is confirmed on a real file, and
the only split it enables -- one window group, no survey frames -- is
exact by the same argument as PASEF. A second objection, to stop treating
survey frames as making the schedule vary, waits for that run. Without
keeping a spot's survey and fragment frames apart first, the split would
cut survey frames into precursor blocks by scan range.

**What changes in a store.** Only stores of files with type-9 frames:

- `uns["msms_schedule"]` and `msi_metadata.ms_analysis.fragmentation` gain
  the windows (targets, offsets, collision energies, scan ranges), the
  collision-induced dissociation term, and `merges_precursors: true` when
  the windows isolate more than one precursor.
- One window group and no survey frames: `{table}_msms` is now written by
  default. It was refused before ("records no precursor").
- Survey frames among the fragment frames: still refused, for the same
  reason as before (not constant). The schedule is recorded.
- Window groups that take turns, without survey frames: refused as not
  constant, where it was refused as recording no precursor.
  `constant_across_pixels` in both blocks goes from true to false.

PASEF, single-precursor, MS1 and TSF stores do not change. The
merged-precursor `WARNING` names no per-pixel count when the schedule
varies.

**Known limits.**

- No MALDI iprm-PASEF run has been measured. Whether it writes type 9,
  keeps survey frames, or uses more than one window group is still open;
  the split is tested on a synthetic file only.
- A spot's survey and fragment frames are still summed into one pixel row.
- The half-open reading of a DIA scan range is inferred from shared
  boundaries, not compared with the vendor's own per-window output.
- `PrmFrameMsMsInfo` (prm-PASEF, `MsMsType` 10) is not read. If a MALDI
  iprm-PASEF run turns out to use it, that run records MS/MS with no
  windows and is refused, as every type-9 run was before this.
- A file with rows in both `PasefFrameMsMsInfo` and the DIA tables records
  the PASEF windows only. It is still refused, as not constant. No such
  file is known.

## D38. The mass-mobility heatmap is a mean over positions, not spectra

**Status:** Implemented (2026-10-13). Amended by D58 (the divisor is the
summed table's rows).

**Decision.** `uns["mobility_heatmap"]["counts"]` is the heatmap's total
divided by the number of distinct pixel positions it was fed, the rows
`uns["average_spectrum"]` divides by. Two spectra the source places on one
position add their points to the total and count as one pixel, as the
summed table sums them into one row. This holds on every route that builds
the heatmap: the fused frame pass, the separate mobility pass and the
fallback built on first demand. The `WARNING` raised when the heatmap and
the mean spectrum disagree names the vendor centroid as the cause only when
`--tdf-spectrum vendor_centroid` is in use.

**Why.** The heatmap divided by spectra. Its own notes promise that, under
`scan_sum`, its marginal over mobility is the stored mean spectrum; with
two spectra on one position that broke. On the synthetic iprm-PASEF
fixture (seven frames on six spots, one spot holding a survey and a
fragment frame) the heatmap held 6/7 of the mean spectrum's current under
`scan_sum`. The conversion then warned that the vendor centroid keeps only
part of the current, a centroid the conversion was not using.

No real mobility source with two spectra on one position is known. A
split Waters run records one stage position twice (D26), but Waters has no
mobility pass, so no heatmap. The nine local timsTOF imaging sets, MS1 and MS/MS (types 2 and 8), have one frame per
raster spot. Over each set's own mobility pass the accumulator saw as many
positions as spectra, so the divisor, and every stored count, is unchanged
on all of them:

| Set | Spectra | Positions |
|---|---|---|
| MS/MS, type 2 | 315 | 315 |
| MS/MS brain, type 8 | 713 | 713 |
| MS/MS brain, type 8, second | 668 | 668 |
| MALDI-2 slide, long ramp | 1,465 | 1,465 |
| biofilm, MALDI-2, 20 um | 18,838 | 18,838 |
| biofilm, MALDI, 20 um | 19,562 | 19,562 |
| MALDI-2 slide, short ramp | 26,087 | 26,087 |
| glycans, high mass | 65,856 | 65,856 |
| rat brain, 9,990 scans | 71,249 | 71,249 |

**The objection.** Divide by the row count the summed table ends with,
handed over by the converter, so the two denominators are one number. It
did not win. The separate mobility pass covers every plane of a
multi-plane store, while each plane's table has its own rows, and the
fallback has no table at all. Counting distinct positions inside the
heatmap gives the same number wherever a table exists and needs nothing
from its caller.

**What changes in a store.** Only stores whose source puts two spectra on
one position, under any spectrum mode: `counts` rises by the ratio of
spectra to positions, and `current_ratio` with it (from 0.857 to 1.0 on the
fixture). The misplaced `WARNING` goes. Stores of sources with one spectrum
per position do not change.

**Known limits.**

- A frame with no points never reaches the heatmap. A position whose
  points all fall outside the mass axis, or one off the grid, still counts
  as a pixel while the table gives it no row, as before. Neither has been
  seen on a real file. (Lifted by D58: a position with no row is no pixel.)
- On a store with more than one plane, the heatmap covers every plane
  while each table's mean spectrum covers its own; `current_ratio` there
  compares the two and is not 1. The new `WARNING` says so.
- When survey and fragment frames on one spot are kept apart, as D37
  expects to happen, the heatmap must follow whatever the summed row
  keeps. Today both go into the row and into the heatmap.

## D39. mzPeak position z is read from a scan's parameters list too

**Status:** Implemented (2026-10-13).

**Decision.** The mzPeak reader reads position z (`IMS:1000052`) from the
scans' z column when there is one, and from each scan's `parameters` list
otherwise. A scan's z is the first parameter in its list that carries that
accession or the name "position z"; it must state z as an integer, or the
scan states none. Everything D30 says about a z the
archive states then holds for both places: one plane records its z, and
spectra on several planes, or z stated for some spectra only, are refused.

**Why.** mzpeak-convert 0.12.0 and 0.14.0 do not drop z. They write no z
column, but put `IMS:1000052` into every scan's `parameters` list. Read only
from a column, an archive from those versions with two planes converted
with exit 0 into one plane, each pixel the sum of both. The same archive
with a z column was refused. D30's first known limit rested on the same
misreading: the glioma archive from 0.12.0 states z, in the list.

The scans join their spectra on `source_index`, so each spectrum takes the
z of its own scan, whatever order the scans member lists them in.

**What changes in a store.** Only archives from mzpeak-convert 0.12.0 and
0.14.0 whose scans state z in the list:

- One plane: `coordinate_offsets_px` on z goes from 0 to the stated z. The
  public archives `gbm_Test_P15_r2` and `ltpmsi-chilli` state z = 1 on
  every scan, so they record 1, as their imzML files do.
- Several planes, or z for some spectra only: refused, where the planes
  were summed into one.

Archives from 0.16.0 and 0.17.2 write a z column and no z parameter; their
stores do not change. Region numbers read from the same list (D33) are
read by the same code and do not change either; both were compared on
0.17.2 archives from Bruker runs.

**Known limits.**

- The public archives from 0.12.0 were not converted again for this entry.
  Their shape is taken from the stress test of 4.3.3, which read them.
- When an archive has both a z column and z in the list, the column wins
  and the list is not compared with it. No such archive is known.

## D40. A position measured twice is one summed row in the sibling tables too

**Status:** Implemented (2026-10-13).

**Decision.** When a source gives one raster position more than one
spectrum, the mobility-resolved table and the demultiplexed MS/MS table sum
the repeated entries into one, as the summed table already does (issue
#241). This holds on both routes that fill them: the passes fused into the
summed table's own (D5), and the standalone passes that read the source
again.

**Why.** Only the fused route was told about the repeat. The standalone
builders (grid, shared feature axis, MS/MS) stored both measurements'
entries side by side at the same (row, column). Such a matrix densifies to
the same values, so the dense rows and the current ratios agreed. But it
held more entries than cells, it was not canonical, and the two routes
wrote different stores from one source. On a four-pixel stand-in with
pixel 0 measured again last, the standalone route stored 16 entries for
12 cells in both siblings; the fused route stored 12.

The converter counts the repeats in its first pass and now tells every
sibling builder, per table. The two mobility builders also notice a row
reached twice themselves, since one spectrum is one row there. That
catches a repeat the summed table never saw: the shared feature axis keeps
features outside the mass range, so a spectrum that resamples to nothing
still fills the mobility table. A builder that knows of a repeat sums it
with the same chunked sort and merge the summed table uses. A source with
no repeat pays nothing.

**What changes in a store.** Only sources that reach the standalone passes
and measure some position twice. Today those are imzML files with a
mobility array: their `_mobility` table's `X` holds one entry per (row,
column), so `X.nnz` drops by the repeated entries. The dense values do not
change. A Bruker `.tdf` read through the vendor library always takes the
fused route and does not change. No shipped reader feeds the MS/MS table
from the standalone passes, so no `_msms` table changes.

**Known limits.**

- No real source with a repeated position was converted for this entry.
  The synthetic iprm-PASEF copy has a spot holding two frames, but the two
  share no grid cell and it writes no MS/MS table, so it stores the same
  either way.
- The MS/MS builder relies on the converter's count alone, because a row
  legitimately arrives once per precursor there. A repeat the summed table
  does not count would still be stored twice. No shipped reader reaches
  that builder outside the fused route.

## D41. A plane's table holds the mass-mobility heatmap of its own plane

**Status:** Implemented (2026-10-13). Removes the second known limit of D38
(the heatmap of a multi-plane store covering every plane).

**Decision.** A multi-slice source converted as 2D writes one table per z
plane. Each of them now holds, in `uns["mobility_heatmap"]`, the heatmap of
its own plane's pixels: the same pixels its `uns["average_spectrum"]`
averages over. The one raw mobility pass the conversion already made builds
every plane's heatmap, routing each pixel by its z. Under D38 each plane's
heatmap then divides by the distinct positions on that plane, which are the
rows its table holds. A volume (`--handle-3d`) and a single-plane source
keep one heatmap over every pixel, which is again the population of their
one table.

**Why.** Since #243 a plane's `average_spectrum` is its own plane's mean.
The heatmap beside it was still built once over every plane, and every
plane's table was handed the same dictionary. `current_ratio` was then
written into that shared dictionary once per table, comparing the
dataset-wide heatmap with that plane's mean, so the store kept the last
plane's ratio in every table. A two-plane source whose second plane is three
times brighter, converted from the raw sum (`scan_sum`), showed it:

| Table | Stored ratio, before | Ratio of its own table, before | Now |
|---|---|---|---|
| `z0` | 0.6667 | 2.0 | 1.0 |
| `z1` | 0.6667 | 0.6667 | 1.0 |

The heatmap's total was 47 in both tables, against plane means of 23.5 and
70.5. docs/output-format.md promised the marginal equals the
`average_spectrum` beside it under `scan_sum`, and it held on no plane.

**The objection.** Keep one dataset-wide heatmap, give each table its own
copy, and compute the ratio against a dataset-wide mean. It did not win.
That mean is stored nowhere since #243, so the ratio would compare against a
number no reader of the store can see, and the marginal would still differ
from the spectrum it is drawn beside. A viewer showing one plane's ion image
next to a heatmap of every plane misleads the same way the old shared mean
spectrum did. `average_spectrum_per_region` stays dataset-wide for a
different reason: a region has no z (see docs/output-format.md), whereas
every pixel the heatmap reads has one.

**Cost.** One accumulator per plane while the pass runs: about 8 MB per
plane at the default 4,000 m/z bins by 256 mobility channels. The point
buffer (up to 4 million points, 64 MB) is folded in whenever the pass moves
to another plane, so only one plane holds one at a time, as before. No extra
read; the stub conversion in the tests still makes one mobility pass for two
planes. A plane whose pixels yield no mobility spectra gets no heatmap, and
the conversion says so at `WARNING`.

**What changes in a store.** Stores of multi-slice sources with an ion
mobility dimension converted as 2D: each plane's `mobility_heatmap.counts`
becomes that plane's mean frame, and `current_ratio` becomes that plane's
ratio (1.0 under `scan_sum`). Each plane's table stores its own block where
before all held the same values. Volumes, single-plane sources (every Bruker
TDF) and sources without mobility do not change.

## D42. A store records a SHA-256 of its mass axis

**Status:** Implemented (2026-10-13).

**Decision.** Every converted store records `msi_metadata.provenance.mass_axis`,
an object with two members: `sha256` and `n_bins`. `sha256` is the SHA-256 of
the main MSI table's `var["mz"]`, taken as little-endian float64 bytes whatever
dtype the array is held or written in, as lower-case hex. `n_bins` is the
number of values. It is written whether or not the conversion resampled. There
is one value per store: the sibling tables (`_mobility`, `_msms`) carry the
main table's value, because their `mz_index` points into that axis (D6), not
because their own `var["mz"]` hashes to it.

**Why.** Showing that two stores share an axis needed both arrays loaded and
compared. The recorded settings do not answer it: a tof law is filled from
detection, numpy and the platform can change a float, and a `--no-resample`
conversion has no resampling step to read a setting from. A later plan of one
axis for several sources (issue #449) has to check a conversion against that
axis, and a hash is the cheapest thing to check against.

**Where it lives and why.** In `provenance`, which is versioned and validated,
and which every store has. Not on the "mass axis resampling" processing step:
that step is absent under `--no-resample`, and the axis exists all the same.
Not in the root attributes or in `uns["essential_metadata"]`: both are outside
the versioned, validated schema, so a validator could not check the value.
The field is optional, so a standalone document, which has no axis, stays valid
without it. It is not bound to a CV term: MS:1003151 (SHA-256) names a file
checksum, and this is a hash of an array. The gap is listed in
`CANDIDATE_CV_CONCEPTS`.

**What changes in a store.** Every new store gains `provenance.mass_axis`, and
the schema version is 0.11.0. No other value changes. Stores written earlier are
unchanged and still valid. `thyra validate` recomputes the hash from `var["mz"]`
of each plain MSI table that has the field and reports an error at
`provenance.mass_axis` when it differs; a store without the field gets no such
check.

**Known limits.**

- The hash is of the bits. Two axes that differ by one ulp are two axes. This
  is intended: the use is to show that two stores are bin for bin the same.
- When the `msi_metadata` block cannot be built, which is logged, the store
  carries no fingerprint.
- Sibling tables carry the main axis's value, not the hash of their own
  `var["mz"]` (D6).
- `thyra validate` checks every plain table that records the field. A table
  added to a store later, with a different `var["mz"]` and a copy of the main
  table's block, fails that check. Its block must drop or replace the field.
- The field is not a CV term.

## D43. A plan fixes one mass axis for several sources

**Status:** Implemented (2026-10-13). Amended by D56 (2026-10-16).

**Decision.** `plan_mass_axis(sources)` reads each source's header, as
`preview_msi` does, and returns a `MassAxisPlan`. The plan holds the axis
settings only: the axis type, the range, the bin width at one reference m/z
(or, for a `tof` axis, the width law and the bins per peak width), the bin
count, and the D42 fingerprint of the axis those settings lay. The resampling
method is not part of the axis (`AxisPlan`), so each conversion keeps its own
detector's method unless the caller forces one in the plan. The plan decides
each setting the caller leaves open:

- The range is the intersection of the sources' ranges.
- The width is the widest width any source's detector declares (else the
  axis type's default), each first stated at the reference m/z under the axis law (width grows as
  m^0 on `constant`, m^0.5 on `linear_tof`, m^1 on `reflector_tof`, m^1.5 on
  `orbitrap`, m^2 on `fticr`).
- The plan is refused when the ranges have no m/z in common, when the
  detectors want different axis types, and when `tof` sources declare
  different width laws. A caller's value is used as given.

A conversion given a plan builds the axis from the plan's settings and
refuses, before reading any spectrum, when the fingerprint of the axis it
built differs from the plan's.

**Why.** Comparing two stores bin for bin needs one `var["mz"]` for both.
Five settings rebuild the axis, and the D42 fingerprint proves they did. The
width is the widest because, on an axis finer than a source's own peaks, that
source's centroids scatter over neighbouring bins from pixel to pixel and one
ion becomes several features. A coarser axis costs the finer source
resolution instead, and the plan names it in a warning. The range is the
intersection so that every bin is measured in every store.

**What changes in a store.** Nothing. A plan only picks settings that already
exist. The new check refuses a conversion; it writes nothing.

**Known limits.**

- An imzML range is the span of the data observed, not the acquisition window,
  so an intersection bound set by an imzML source may cut signal near the edge.
  The plan warns. An mzPeak range is observed data too when the archive states
  no scan window, but that cannot be told from the header cheaply, so it is
  not flagged.
- The source with the finer width loses resolution on the shared axis, and
  under `tic_preserving` bins coarser than a source's own can merge narrow
  peaks (#442). D50 fixes the second part.
- There is no plan for "the axis of an existing store" yet.

## D44. Only a reader that can select a region accepts the region option

**Status:** Implemented (2026-10-13). Superseded in part by D57
(2026-10-16): one region now gives a region map, so a single-region store
keeps the source's region number and name. Superseded in part by D83
(2026-10-22): the Waters reader selects a region too, and the gate is a
reader attribute (`selects_regions`), not a timsTOF check.

**Decision.** `region` selects one region, and only the timsTOF reader can do
that. Every other reader refuses it, through the API as well as the CLI. The
refusal lives in `BaseMSIReader`, so a new reader refuses it until it learns
to honour it. `thyra.list_regions(path)` lists the regions of any input from
its header: each region's number, `.mis` Area name, spectrum count and box.
Two readers record more than before:

- Rapiflex builds a region map from the position log. A spot named
  `R{n}X{x}Y{y}` is region `n`, at its raster cell minus the `.dat` header
  origin. The name is the n-th `.mis` Area, the pairing D31 already uses. As
  on timsTOF, one region gives no map.
- solariX adds each region's `bounds` to its region summary, in the store's
  pixel grid, as timsTOF and mzPeak give them.

**Why.** solariX, Rapiflex, imzML, PHI and Waters dropped `region` in silence
through the API, and converted every region of an input the caller asked to
cut. On mzPeak the CLI logged it as ignored and the reader then refused it.
One rule gives one message. A caller who wants one store per region needs the
region list first, and had to open the vendor database to get it.
`list_regions` always returns at least one entry, so a loop over it works on
any input; `MsiPreview.regions` stays `None` for one region, which answers a
different question.

**What changes in a store.**

- A multi-region Rapiflex store: `obs["region_number"]` holds the
  position-log region instead of 1 everywhere, `uns["regions"]` has one entry
  per region with `n_spectra`, `bounds` and `name`, and
  `uns["average_spectrum_per_region"]` appears. A single-region Rapiflex store
  is unchanged.
- A solariX store: each entry of `uns["regions"]` and of
  `format_specific.regions` gains `bounds`.
- A call that passed `region` to another format now fails instead of writing a
  store of the whole input.

**Known limits.**

- The Rapiflex keying was checked on one real single-region acquisition, where
  every position-log spot landed on its offset-table cell. No multi-region
  Rapiflex acquisition was available. When a spectrum has no position-log spot,
  or a cell is in two regions, the reader warns and writes no map.
- solariX regions have no names. The positional Area pairing is not verified
  on solariX.
- `n_spectra` counts spectra, not pixels: a position measured twice (D40)
  counts twice.
- Region selection for solariX, Rapiflex and mzPeak, which know their regions,
  is not part of this decision.

## D45. The metadata block names the sample and the spectrum representation

**Status:** Implemented (2026-10-13).

**Decision.** Schema 0.12.0 adds three optional facts to `msi_metadata`:

- `sample.sample_id` (MS:1000001 sample number) and `sample.sample_name`
  (MS:1000002 sample name). Only the person converting fills them, through
  `convert_msi(sample_id=..., sample_name=...)` or `--sample-id`. Thyra never
  reads them from the source.
- `ms_analysis.spectrum_representation`, `"profile"` or `"centroid"`, with its
  term (MS:1000128 / MS:1000127, under MS:1000525 spectrum representation).
  It is the representation the reader delivered, the same fact
  `uns["essential_metadata"]["spectrum_type"]` records. That key stays.

The Bruker timsTOF extractor now states the representation too. It used to
state none. TSF is `centroid` (the vendor line spectrum, D1). TDF read with
`vendor_centroid` is `centroid`. TDF read with `scan_sum`, the default, is
`profile`: the summed digitizer indices are unpicked and equal Bruker's own
quasi-profile export (D1), the same reasoning that labels PHI `profile`.

**Why.** Code that combines stores of one study had two problems. Whether two
stores can be combined depended on the spectrum representation, which sat
outside the versioned block. For example, profile and centroid stores should
not share a `tic_preserving` axis. And the stores could be told apart only by
a file path, which is not a sample. The `Provenance` docstring said the root
attributes recorded the representation; they never did. That text is
corrected.

**Where it lives and why.** The representation sits in `ms_analysis`, beside
polarity, because it describes the acquired spectra. Resampling does not
change it: Thyra never picks peaks (D1), and the resampling step is recorded
in `processing`. The sample fields sit in `sample`, the section METASPACE and
SDRF both use for sample identity.

**What is left unset.** An imzML file that declares neither MS:1000127 nor
MS:1000128 gets a guess: processed mode is taken as centroid, with a warning.
The guess still picks the resampling method, and it is still written to
`uns["essential_metadata"]["spectrum_type"]`. The block leaves the field
unset, by the builder's rule that unset beats guessed. The sample fields are
never filled from a vendor `SampleName`, an imzML sample list, or
`dataset_id`. Any of these can hold a placeholder or a person's identifier,
and the block is copied into every metadata document Thyra writes.

**Objection that did not win.** "Take the sample from the file." Bruker
`GlobalMetadata.SampleName` and imzML `sampleList` exist, but both hold free
text typed at the instrument or written by the exporting software, which can
be a placeholder or a person's identifier. Promoting that text into a field
named "sample id" would present it as a curated identifier.

**What changes in a store.** Every new store has schema version 0.12.0 and
`ms_analysis.spectrum_representation` wherever the reader states one. Bruker
timsTOF stores also gain `uns["essential_metadata"]["spectrum_type"]`, which
was `None`. Their resampling defaults do not change when the instrument name
is read: the timsTOF detector matches on it before any spectrum type is
consulted. A `.d` whose name cannot be read, or lacks "timsTOF", used to fall
to the default detector (constant axis). Read as centroid (TSF, or TDF with
`vendor_centroid`), it now reaches the centroid detector, whose reflector-TOF
axis is the one the timsTOF detector would have chosen. When the Bruker
library cannot centroid, the reader falls back to `scan_sum` with a warning;
the processing step and the representation now say `scan_sum` and `profile`,
where the step used to repeat the request. A store
converted with a sample id or name carries them in `sample`. Stores written
earlier are unchanged and still valid.

**Known limits.**

- One conversion names one sample, even when the slide holds several
  regions. Each region has its own name in `uns["regions"]`.
- The sample id and name are stored as typed. Thyra cannot tell a patient
  identifier from a code, so the documentation says not to use one.
- `thyra metadata` and the preview read a timsTOF `.d` with the default
  `scan_sum`, and an imzML without `--spectrum-type`. Their document can
  therefore differ from what a conversion with other options stores.

## D46. A processed imzML with mobility gets the heatmap, and the grid when asked

**Status:** Implemented (2026-10-13). Amended by D82 (2026-10-22).

**Decision.** The range a per-pixel mobility source is binned over is the
smallest and largest finite mobility value it holds. A reader states it
through one method, `get_mobility_value_range()`. A source with a shared
axis (Bruker TDF, a continuous imzML export) reads it off the axis values,
as before. An imzML whose pixels carry their own mobility values (processed
mode, as TIMSCONVERT writes it, or a continuous file whose arrays differ)
reads it in one pass over its mobility arrays, and only those arrays, the
first time anything asks. The result is cached, so the heatmap and the grid
get the same two numbers. From this:

- Such a file gets `uns["mobility_heatmap"]` by default, like every other
  mobility source (D3).
- `--mobility-grid`, `--mobility-bins`, `--mobility-min` and
  `--mobility-max` are accepted on imzML. On a per-pixel file they write
  `{table}_mobility` with the grid's contract: `uns["mobility_grid"]`, 256
  channels on the heatmap's own edges, and `uns["mobility_marginal"]`.
- On a continuous file with a shared axis the table is read off the shared
  features, as before. `--mobility-grid` is ignored there, and a `WARNING`
  now says so.
- A point whose mobility value is not a finite number has no channel. The
  imzML reader leaves it out of the mobility point cloud and says how many
  at `WARNING`. The summed table keeps it.

**Why.** Before this, a processed imzML with a mobility array kept no
mobility at all. The log told the user to pass `--mobility-grid`, and the
CLI then called that flag ignored on imzML. Processed mode is how
TIMSCONVERT writes timsTOF imaging data, and a recent review of ion mobility
imaging names unread mobility in such files as a bottleneck of the field
(doi:10.1021/acsmeasuresciau.6c00282). Nothing in the converter was specific
to TDF except where the range came from.

The range comes from the values because only they are sure to contain every
point. TDF uses its values for the same reason: its per-scan 1/K0 overhangs
the declared acquisition range. PSI-MS has terms for a lowest and highest
observed mobility (`MS:1003437`, `MS:1003438`), but Thyra reads neither, and
a stated range narrower than the data would pile points into an edge
channel.

**Objections considered.**

- *A whole extra read for a range.* It reads the mobility arrays only: a
  third of a 64-bit `.ibd`. The heatmap itself is already one more full pass
  on any source other than TDF. `--no-mobility-heatmap` without a grid reads
  nothing extra.
- *A range stated in the file would be free.* It would be, when present and
  right. Neither can be checked without reading the values.

**What changes in a store.**

- A processed imzML (or a continuous one with differing arrays) with a
  mobility array: the summed table gains `uns["mobility_heatmap"]`.
  `uns["mobility_axis"]` is unchanged and records no range: the range is in
  the heatmap's `mobility_edges` and the grid's `edges`.
- The same file with `--mobility-grid`: a `{table}_mobility` table, its
  `uns["mobility_grid"]` and `uns["mobility_marginal"]`,
  `uns["mobility_axis"]["resolved_table"]` on the summed table, and
  `msi_metadata.ms_analysis.ion_mobility.grid`. When the array states 1/K0
  in `MS:1002814`, the table also carries `var["ccs"]` and `uns["ccs"]`, as
  every 1/K0 mobility table does (D47). That lifts the last known limit D47
  lists.
- A continuous imzML whose shared mobility axis holds a non-finite value:
  the heatmap was refused with an error before (the channel index was
  undefined) and is now written without those points.
- Bruker TDF and continuous imzML with a finite shared axis: nothing.

**Known limits.**

- The grid's marginal reproduces the summed table only when the summed
  table puts the same points on their nearest bin: no resampling, or
  `nearest_neighbor`. Under `tic_preserving` it does not. On a resampled
  axis a point exactly at the data's m/z edge can also lie just outside the
  first or last bin centre; the heatmap and the grid drop it while the
  summed table folds it into the edge bin (#347). `current_ratio` records
  both. D52 lifts the second part: the heatmap and the grid now keep the
  points the summed table keeps.
- A multi-plane imzML converted as 2D reads the whole file once per plane
  for the grid's discovery pass.
- With `--no-resample` the raw union axis times 256 channels can pass the
  grid's size ceiling (D4).

---

## D47. A 1/K0 mobility table carries each feature's CCS, for charge 1

**Status:** Implemented (2026-10-13). Reverses the rule, stated in the code
since the CCS function shipped in v3.13.0, that the converter never writes a
cross section.

**Decision.** A mobility table whose axis is 1/K0 in V s/cm2 carries
`var["ccs"]`: each feature's collision cross section in square angstroms,
float64. It is computed from the feature's own `mz` and `mobility` by the
Mason-Schamp equation, for a charge of 1, in nitrogen (28.0134 Da) at
305 K. The same table carries `uns["ccs"]`, which records the quantity and
unit (`MS:1002954`, `UO:0000324`), the method, the assumed charge, the gas,
its mass and the temperature, and that these conditions are a convention.
Keys are plain names and accessions are values, as in the other `uns`
blocks. One function adds both, where the two mechanisms that fill the
table meet, so a continuous imzML export and a TDF grid, per plane or as a
volume, all get them. The schema reserves the column name (0.13.0); the
document is unchanged.

The axis must say it is 1/K0 (`MS:1002815`) and that its unit is
`MS:1002814`. A drift-time axis gets neither the column nor the block.
Neither does an axis whose quantity or unit is not stated.

**Why.** A cross section is what IM-MSI work compares across instruments
and against libraries. 1/K0 depends on the gas, the temperature and the
calibration; a CCS is meant to be the ion's own. Every consumer of a TIMS
table would otherwise repeat the same conversion with the same assumptions,
and some would get the constants wrong: 28.013 Da instead of 28.0134
already moves the result by 7e-6, and 304.85 K instead of 305 K by 2.5e-4.

**Where the conditions come from.** From the convention, not from the run.
Bruker's `tims_oneoverk0_to_ccs_for_mz` takes no gas or temperature; it uses
nitrogen at 305 K, and the formula reproduces it with these two values. The
TDF holds no drift-gas or temperature field that Thyra reads: GlobalMetadata
gives the 1/K0 range, `TimsCalibration` gives unnamed coefficients, and the
mobility keys of `CalibrationInfo` are zero on the one imaging TDF read.
`conditions_source: "convention"` says so, so a reader does not take the
values for the run's conditions.

**Why the formula, never the vendor library.** The same formula fills the
column on every route and every platform, so a store does not depend on
whether the machine's library exports the conversion. The library also
converts one value per call, which is millions of calls on a grid table.
The two agree to 1e-6 relative from m/z 100 up and to 2.2e-6 at m/z 50,
measured on the synthetic TDF and on a probe from m/z 50 to 3000. The gap
grows as the ion gets lighter and no single gas mass removes it.

**The charge.** MALDI ions are almost always singly charged, and no source
states a feature's charge. Recording 1 as an assumption, in the block
beside the values, is not the same as storing a value the source does not
state (the objection D30 answers): the column is derived, and the block
says from what. A consumer who knows a feature is doubly charged recomputes
it with `mason_schamp_ccs`. The stored value cannot be rescaled, because
the charge also changes the ion's mass in the reduced mass.

**Alternatives.**

- *No column; leave it to consumers.* The previous rule. It keeps the store
  free of an assumption, but every consumer repeats the conversion and
  chooses its own constants, so two tools give two answers for one store.
  The function was already public and tested; the gap was that nobody had
  to call it.
- *A per-feature charge column.* Honest only if something states the
  charge. No reader does: isotope spacing would have to be read off each
  feature, which is annotation, not conversion. A column full of 1s would
  claim per-feature knowledge the converter does not have; one scalar in
  `uns` claims only what is true.
- *One column per charge (`ccs_z1`, `ccs_z2`).* Doubles the column for a
  case MALDI rarely has, and still leaves the consumer to pick.
- *Compute through the vendor library on TDF and the formula elsewhere.*
  Two routes, two sets of values for the same pairs, and a store that
  depends on the machine. Rejected for the reasons above.
- *A drift-time CCS.* Needs the cell length, field, pressure and
  temperature, or single-field calibrants. No source Thyra reads states
  them, and Waters drift-time data is not read at all.

**Cost.** 8 bytes per feature in `var`, and in the projection D4 makes
before building a grid table (`VAR_BYTES_PER_FEATURE` is 338 now). The
values are computed in chunks of about one million features, so the
formula's temporaries stay a few MB. On a synthetic 8-million-feature grid
the peak of building and wrapping `var` rose by at most 8 bytes per feature.
D8 stands: `var` is still the memory peak.

**What changes in a store.** Every mobility table on a 1/K0 axis gains
`var["ccs"]` and `uns["ccs"]`, and the schema version in every table's
`msi_metadata` is 0.13.0. Every other array is unchanged. Stores written
earlier are unchanged and still valid. `thyra validate` checks that a `ccs`
column, on whatever table carries it, is numeric, finite, positive and as
long as `mz`.

**Known limits.**

- The value is for charge 1. A multiply charged feature's stored CCS is
  wrong, and the store cannot tell which features those are.
- On a grid table the pair is the mass bin and the channel centre, so the
  value is that cell's, not one peak's. The 1/K0 it rests on can be up to
  half a channel from a peak's own.
- A 1/K0 array that states no unit gets no column. If real TIMSCONVERT or
  TIMSImaging exports turn out to omit the unit, this gate is the place to
  relax.
- Processed imzML with a mobility array keeps no mobility table, so no CCS.

---

## D48. A store in micrometres places its raster where the source did

**Status:** Implemented (2026-10-14). Amended by D69 (2026-10-19): an
optical image with no alignment is not in `"global"` at all.

**Decision.** A store whose `"global"` is in micrometres gets one more
coordinate system, `"acquisition"`. Every TIC image, pixel polygon element
and optical image placed in micrometres carries a transform to it: the
element's `"global"` transform followed by a `Translation` to the source's
origin. The root attribute `coordinate_systems` gains an `"acquisition"`
entry beside `"global"`, with `origin_um`, `origin_source` and
`raster_to_acquisition_affine`. `"global"` and the per-dataset system do not
change.

The origin is where grid position 0 sits in the source's own frame:

- Bruker timsTOF, rapifleX, solariX, imzML, mzPeak: the index the source's
  own numbering gives grid position 0, times the pitch of each axis (D12).
  `origin_source` is `"raster_index"`. That index is what the reader
  records as `coordinate_offsets_px`: the declared imaging area's or the
  selected region's minimum (timsTOF), the smallest raster index (solariX),
  the base the file numbers from, 1 or 0 (imzML, D14, so a cropped file
  keeps its own coordinates), the `.dat` header's first raster position
  (rapifleX, whose indices are already 0-based), or that base plus the
  writer's `position_offset` (mzPeak, D29).
- PHI: 0, which it states.
- Waters: the fitted lattice origin, the smallest stage reading on each axis,
  in micrometres. `origin_source` is `"stage_position"`. It is new: Waters
  stored no origin before. It reaches the converter as
  `EssentialMetadata.stage_origin_um`, because it is not an integer index and
  cannot go into `coordinate_offsets`.

On a volume the TIC volume also moves in z, by the z offset times the slice
pitch. Pixel polygons stay flat, as they are in `"global"`.

**Why.** Readers re-zero the raster. Before this, the subtracted offset
survived only as two root attributes, `coordinate_offsets_px` and
`stage_offset_um`, and no element transform carried it. Two stores cut from
one slide, such as two `--region` runs of one `.d`, landed on top of each
other in every SpatialData reader. Placing them took Thyra-specific code.
Now any reader that follows transforms places them.

**The name.** `"acquisition"`, not `"stage"`: on every format but Waters the
origin is a raster index times the pitch, a position in the source's own
grid, not a measured stage position. `"raster"` was rejected because it reads
as pixel indices.

**When it is not written.**

- A store placed in a FlexImaging scan's pixels (D31). There `"global"` is
  already the scan's frame, and the alignment affine absorbs the offset.
- A source that reports no origin, such as a Waters grid without a fitted
  lattice. A zero translation would claim an origin nobody stated.
- An optical image with no alignment. It sits at `Identity` in `"global"`,
  which says nothing about where it lies relative to the raster, so moving
  it by the raster's origin would place it falsely.

**What `stage_offset_um` means now.** It was defined as
`coordinate_offsets_px` times the pixel size (D12, D32). It is now "where
grid position 0 sits in the source's frame", the same numbers as
`acquisition.origin_um`. For every format that has index offsets the value
is unchanged. Waters stores gain it.

**Objection that did not win.** "Bump `convention_version`; the attribute
has a new top-level key." The key is additive, like `z_spacing_um` and
`raster_to_global_affine` before it. A consumer that reads only `"global"`
is unaffected, and a bump would make every existing consumer warn on every
store.

**Reserved name.** A `dataset_id` of `"acquisition"` is refused, at the
front door and by the converter. The per-dataset system is named by the id,
so it would silently become this one. `"global"` stays allowed: both systems
would carry the same transform.

**What changes in a store.** Every micrometre store whose source states an
origin gains the `"acquisition"` transform on its TIC image, pixel polygons
and placed optical images, and the `coordinate_systems.acquisition` entry.
Waters stores also gain `stage_offset_um`, and their raw metadata gains
`imaging_grid.stage_origin_um`. FlexImaging-aligned stores do not change.
The metadata schema moves to 0.14.0, which reserves the coordinate system
name, as 0.10.0 reserved an `obs` column (D26). The document's fields do not
change. Tables carry no transform, so `obs["spatial_x"]` and
`obs["spatial_y"]` stay in `"global"`.

**Known limits.**

- imzML and mzPeak number from 1 (D14), so most of their stores have an
  origin one pitch from 0. Two files exported from parts of one acquisition
  line up only if the exporter kept the acquisition's own numbering.
- A user-supplied `--pixel-size` scales an index offset like any other
  pitch. On Waters it does not move the stage origin, which is measured.
- A volume's z origin uses the slice pitch. When that pitch is
  `assumed_isotropic`, the z origin is as unreliable as the depth.

---

## D49. A release is a pull request, and `main` requires passing checks

**Status:** Implemented (2026-10-14), #474. Supersedes D18. The first release
pull request, #478, opened on the merge and ran its checks; the required-checks
ruleset and auto-merge were turned on after it.

**Decision.** release-please replaces python-semantic-release. After every
merge it keeps one release pull request open, `chore(main): release X.Y.Z`,
which bumps `pyproject.toml`, `thyra/__init__.py`, the `thyra` entry in
`uv.lock` and `CHANGELOG.md`. Merging that pull request is the release:
release-please tags the merge commit and creates the GitHub Release, then
`release.yml` builds the tag and publishes it to PyPI once **Tests** has
passed on that commit. A Monday cron merges the release pull request when
**Tests** passed on it, so releases still come weekly without anyone acting.
`main` gets a ruleset requiring the CI checks, not strict, with no bypass
list, and the repository allows auto-merge.

**Why.** D18 measured that a `GITHUB_TOKEN` push of the version commit cannot
pass a required-checks rule, and that the GitHub Actions app cannot be a
bypass actor. Both routes it left open needed a new credential. This removes
the push instead: the version commit reaches `main` through a pull request,
like any other change, so the rule needs no bypass and no credential exists to
store or rotate. It is the setup Ousia already runs.

**The alternative, and why not.** A GitHub App holding Contents write, on the
ruleset's bypass list, with `release.yml` pushing as the App. It keeps
semantic-release, but adds a private key to store, and checked against
python-semantic-release 10.6.1 it is not the one-line change it looks like:
`remote_settings.token` puts `GITHUB_TOKEN` into the push URL, so the checkout
token alone does not change who pushes. Making the App the pusher everywhere
would let its GitHub Release start a second `release: published` run.

**What `GITHUB_TOKEN` costs, and how each is handled** (all in
`release-please.yml`):

- Runs on a pull request it pushed are parked for approval and report no
  checks. The workflow approves them (`.github/scripts/approve-parked-runs.sh`).
- A tag it creates starts no workflow. The workflow dispatches `release.yml`
  at the tag.
- A merge it makes, on the cron path, starts no `push` run. The workflow runs
  release-please again in the same job and dispatches **Tests** at the tag,
  which the publish gate waits for.

**What changes for releases.**

- The Monday cron no longer releases by itself; it merges the release pull
  request if **Tests** passed on it. Checks still running defer the release
  to a dispatch (`gh workflow run release-please.yml -f merge_release_pr=true`)
  or the next Monday.
- The changelog lists `feat`, `fix` and `perf` commits only. Under
  semantic-release it also listed `docs`, `refactor`, `test`, `build` and
  `style`. release-please releases whenever a visible section is non-empty, so
  those types are hidden to keep them from triggering a release.
- New changelog sections use release-please's format. Older ones stay.
- `release.yml` runs only on a dispatch at a `vX.Y.Z` tag and refuses any
  other ref. Re-publishing a tag is `gh workflow run release.yml --ref vX.Y.Z`.
  The `publish_only` input and the `release: published` trigger are gone.
- The PR title check reads its allowed types from `release-please-config.json`.
- `python-semantic-release` left the `dev` group, and with it ten packages
  that only it needed.

**The ruleset.** Required: `lint`, `test` and `integration` on Ubuntu and
Windows for Python 3.12 to 3.14, `clean-venv-install` for 3.12 to 3.14, and
`Conventional commit title`. Not required: `Build the site (strict)`, because
`docs-check.yml` is path-filtered and a code-only pull request would wait for
it forever; and `complexity-check`, which is advisory.

---

## D50. A bin wider than the source sampling is rebinned, not sampled

**Status:** Implemented (2026-10-13), issue #442.

**Decision.** `tic_preserving` checks each spectrum for a sample that
interpolation at the bin centres would skip: a bin holding a sample with
intensity that is neither of the two samples around its centre. That is
three samples in one bin, or two on the same side of its centre; a sample
exactly on the centre is read alone, so any other in its bin counts. A skipped
zero loses nothing and does not count, so a repeated explicit zero leaves a
spectrum alone; two samples at one m/z count when either has intensity,
because interpolation reads only one of them.

- A spectrum with no such bin is resampled exactly as before.
- A spectrum with one is rebinned wherever a bin holds two or more samples.
  Each sample stands for its cell, the stretch of m/z from the midpoint with
  the sample before it to the midpoint with the sample after it. Every bin
  takes the share of each cell it overlaps. Where the spacing jumps by more
  than a factor of two, at the edge of a gap, the cell side facing the gap is
  cut back to the width of the other side.
- The spectrum's other bins are still interpolated at their centres, then
  scaled by their width over the spacing of the two samples around the
  centre, so both kinds of bin hold ion current.

The rescale to the preserved TIC then runs as before. The gap tolerance
applies only to interpolated bins: a rebinned bin holds samples and is
measured.

**The defect this replaces.** The operator evaluated the linear interpolant
at each bin centre and rescaled the result to the TIC. Once bins are wider
than the source spacing that is sampling, not binning. A peak narrower than a
bin that falls between two centres is mostly lost, and the rescale hands its
ion current to whatever the centres land on. The TIC still balances, so
nothing downstream can see it. A spectrum sampled every 0.01 Da with two
peaks of 0.05 Da FWHM, onto 0.73 Da bins, against every sample summed into
its nearest bin:

| | total | within 1 Da of 167.05 | within 1 Da of 333.10 |
|---|---|---|---|
| samples summed into bins | 56,886.4 | 20,788 | 6,536 |
| before | 56,886.4 | 421 | 415 |
| after | 56,886.4 | 20,788 | 6,536 |

On `imzml_processed_profile_a` (2,965 pixels, sampled every 0.011 Da) onto a
constant axis of 2,600 bins (0.731 Da), the mean spectrum against every
source spectrum summed into the same bins:

| | Pearson r | rel. L1 | sum ratio |
|---|---|---|---|
| before | 0.132 | 1.47 | 1.000 |
| after | 1.000 | 0.057 | 1.000 |
| after, without the gap cut | 1.000 | 0.076 | 1.000 |

The file is processed and has gaps without an explicit zero at their edges,
which is what the gap cut is for.

**Why this rule.** Synthetic profiles were built whose samples are the
integral of a known peak shape over each sample's cell, so the ion current
each bin should hold is known exactly. Error is rel. L1 against that. The
sqrt rows are a zero-suppressed trace on a `sqrt(m/z)` grid onto `linear_tof`
bins, like the Waters profile route; the constant rows are a constant grid
onto constant bins, like rapifleX.

| bins per sample | before | after | rebin every bin holding two |
|---|---|---|---|
| 0.2 (constant) | 0.2768 | 0.2768 | 0.2768 |
| 1.08 (sqrt) | 0.0244 | 0.0244 | 0.0245 |
| 1.32 (sqrt) | 0.0188 | 0.0188 | 0.0202 |
| 1.40 (constant) | 0.0920 | 0.0920 | 0.0943 |
| 1.90 (sqrt) | 0.0210 | 0.0254 | 0.0254 |
| 2.01 (constant) | 0.0778 | 0.0925 | 0.0925 |
| 2.32 (sqrt) | 0.0464 | 0.0253 | 0.0253 |
| 2.49 (constant) | 0.1773 | 0.0707 | 0.0707 |
| 3.31 (sqrt) | 0.1226 | 0.0198 | 0.0198 |
| 4.14 (sqrt) | 0.1939 | 0.0163 | 0.0163 |
| 8.02 (constant) | 0.7637 | 0.0128 | 0.0128 |
| 29.3 (constant) | 1.7784 | 0.0029 | 0.0029 |

Sampling at the centres is as good or better while no sample is skipped,
which on an even grid holds below two bins per sample. Past that it falls
apart. Rebinning every bin that holds two samples, whatever the spectrum, was
measured too (last column): it changes spectra that lose nothing, and costs
up to 0.0023 there. Rebinning only the bins that skip a sample was worse
between two and three bins per sample (0.1315 against 0.0707 at 2.49): the
interpolated bins beside them then hold the wrong share. Once a spectrum is
rebinned at all, every bin holding two samples is.

**Why cells and not trapezoids.** Spreading the trapezoid between each pair
of neighbouring samples over the bins it overlaps was measured too. In a
separate run, rebinning every bin that holds a sample, it was worse at
every width above 1 bin per sample (0.0229 against 0.0161 at 4.14).
It also smears a sample over its neighbours' bins even when the axis is the
source grid; cells reproduce the source exactly there.

**What changes in a store.**

- Conversions whose bins are wide enough to skip samples. For example
  `--resample-method tic_preserving` with a coarse `--resample-bins` or
  `--resample-width-at-mz`, the method forced onto the default width for a
  source sampled more finely than half of it, or a `plan_mass_axis` plan
  wider than a source's own width (D43). Their values change, and narrow peaks come back.
- Not the Waters profile route, the default for a SELECT SERIES MRT and what
  `--waters-spectrum profile` takes elsewhere. Its bins are 1.14 times the
  sample spacing: some hold two samples, one either side of the centre, and
  none skips one. Measured on `waters_desi_profile_a` with
  `--waters-spectrum profile` (2,436 pixels, 210,636 bins, about 50,000
  samples per spectrum): 0 of 87,507,157 stored values changed. The MRT
  reference run was not on hand; it shares the ratio and the `sqrt(m/z)`
  grid, and a synthetic zero-suppressed MRT trace at 1.14 bins per sample is
  unchanged bit for bit (`test_tic_preserving_wide_bins.py`).
- Not the default rapifleX and Bruker MALDI-TOF route. The source is sampled
  every 24.9 mDa onto a 5 mDa axis, so no bin holds two samples.
- Not `nearest_neighbor`, nor any sibling table: they never use this method.

**Cost.** Deciding that a spectrum skips nothing takes one pass over its
sample spacings, against a bound on the bin widths laid once per axis and
read once per 64 samples. On a real `waters_desi_profile_a` spectrum of
46,490 samples that is 0.18 ms of a 2.75 ms resample, and 0.08 ms of 0.80 ms
on a 15,000-sample MRT-like trace onto a 1.05M-bin axis. A spectrum that is
rebinned costs a few milliseconds more.

**Objection: the rule depends on the data (D15).** D15 refused an automatic
gap tolerance because the stored numbers would depend on a heuristic no flag
records. This rule has no setting to record. Whether a centre skips a sample
is a fact of the axis and the source, the same in every run. The method name
and the Thyra version in the processing step say which operator ran.

**Known limits.**

- Between about 1.9 and 2 bins per sample, rebinning is slightly less
  accurate than sampling (0.0254 against 0.0210; 0.0925 against 0.0778). On
  an uneven grid a few bins there already skip a sample, and the spectrum is
  rebinned.
- A spectrum that is not rebinned keeps the old operator, with its known
  error when the source and target spacing laws differ
  (`_gate_tic_preserving`).
- The gap cut uses a fixed factor of two. A source whose spacing changes by
  more than that between neighbouring samples without a gap would have its
  cells trimmed.

## D51. A mass axis plan moves an imzML header range out by the precision its text states

**Status:** Implemented (2026-10-14).

**Decision.** When `plan_mass_axis` takes a range from an imzML source, it
moves each header bound out by one unit of the last digit its text states:
`100.123457` by 1e-6, `1000` by 1, `1.2E3` by 100. When two texts give the
same extreme, the coarser step is used. A low bound that would reach 0 or
below is halved instead. This supersedes D43's range rule for imzML sources
only; the intersection, the forced bounds and every other format are as D43
describes. A bound the caller forces is used exactly as given. The plan's
`origins` for a derived bound say that an imzML range is its header text
moved out.

**Why.** An imzML header states the lowest and highest observed m/z
(MS:1000528 and MS:1000527) as rounded text. The rounded low bound can sit
above the source's first sample, and the rounded high bound below its last.
A conversion keeps the m/z inside its declared bounds, inclusive, so a plan
built from the header as it stands dropped those samples in every spectrum
and reported nothing. Measured: on `imzml_continuous_profile_a` a one-source
plan lost one sample per spectrum (620 of 3.4e7 counts); with the step, the
planned and unplanned stores both hold 34,446,401.7. A synthetic file lost
15%. The step is read from the text, not the parsed number: `1000` and
`1000.000000` parse to the same float, but only the first can hide a sample
at 1000.4. A fixed relative margin (1 ppm was tried first) covers six-decimal
headers and misses coarser ones. The same file converted without a plan reads
the exact extremes of the decoded data and loses nothing.

**What changes in a store.** A planned conversion of imzML sources lays an
axis up to one stated digit wider at a header-derived end and keeps its edge
samples; the extra bins are zero in every store. Stores converted without a
plan do not change. Nothing released changes: plans are new since 4.4.0.

**Known limits.**

- The step bounds rounding and truncation of the stated text. A header whose
  value is wrong by more than that still cuts samples; the plan cannot tell
  without decoding spectra.
- Only imzML sources are moved. Every other format's plan range comes from
  the same values its conversion uses, so a plan cuts nothing the conversion
  keeps.

---

## D52. The mobility and MS/MS tables keep the m/z range the summed table keeps

**Status:** Implemented (2026-10-14). Amended by D82 (2026-10-22).

**Decision.** The mass-mobility heatmap, the mobility grid table and the
demultiplexed MS/MS table keep a point when its m/z is inside the range the
summed table keeps (D13): the declared range on a resampled axis, which
reaches up to half a bin past the first and last bin centre. A point there
goes to the edge bin, as it does in the summed table. A point beyond the
range is still dropped. This extends D13 to every table built from the raw
points, and lifts the known limit D46 records for its grid and heatmap.

**Why.** These tables map their points onto the common mass axis with the
summed table's nearest-bin rule, but they tested membership against the
strict span `[axis[0], axis[-1]]`, the rule D13 retired for the summed
table. On a resampled axis whose end centres sit inside the data, the
lowest and highest m/z samples of a file were dropped from them while the
summed table kept them. Measured on the `mobility_processed.imzML` fixture
with the command line's default resampling and `--mobility-grid`: the grid
table held 810 of the summed table's 1000 counts (`current_ratio` 0.81,
0.6 in the worst pixel), and the heatmap the same 0.81. Both are now 1.0,
with every pixel at 1.0. The same file under `--no-resample` was already
1.0.

**What changes in a store.** On a resampled axis that is not `constant`,
points in the half-bin margin at either end now appear in the edge bins of
the heatmap `counts`, the grid table's `X` and `var`, and the MS/MS
table's `X` and `var`. `mobility_heatmap.current_ratio`,
`mobility_marginal.current_ratio` (with its pixel minimum and maximum) and
the MS/MS table's `demultiplexed_current` ratios move towards 1. Under `--no-resample` and
on a `constant` axis the kept range is the axis span, so nothing changes.
A continuous imzML's shared-axis mobility table drops nothing by m/z and
does not change. Released stores: the per-pixel imzML grid is new since
4.4.0 (D46); Bruker TDF heatmaps, grid tables and MS/MS tables change only
where a point falls in that margin.

**Known limits.**

- The heatmap's stored `mz_edges` are unchanged. On an axis whose bins
  widen with m/z the declared upper bound sits a fraction of a bin above
  the last stored edge, so a point there is counted in the last row while
  lying just outside the edges that row states.
- Under `tic_preserving` the grid's marginal still differs cell by cell
  from the summed table, as D46 says; only the points each keeps now agree.

---

## D53. Waters drift time is read as a mobility axis, and a short list of bins is its own channel set

**Status:** Implemented (2026-10-14).

**Decision.** Two parts.

1. A Waters function acquired with ion mobility on is converted. Its spectrum
   is the one MassLynx reports summed over the drift bins, and the bins
   become the reader's mobility dimension, as a TIMS frame's scans are:
   `uns["mobility_axis"]` (drift time, `MS:1002476`, in ms, `UO:0000028`;
   `values` the drift time of each bin; `source` `"waters_raw"`;
   `bin_width` the step), `uns["mobility_heatmap"]`, the `{table}_mobility`
   sibling with `--mobility-grid`, and `msi_metadata.ms_analysis.ion_mobility`.
   Each bin is read raw (`getRawMobilityScanDataPoints`), whatever
   representation the summed table was read in.
2. When a mobility axis declares itself a list of equal-width bins
   (`MobilityAxis.bin_width`) with fewer bins than the 256-channel anchor,
   the heatmap and the default grid use one channel per bin, edges half a
   bin outside the first and last value. Every other axis keeps 256
   channels over its values. `--mobility-bins` now defaults to "not given",
   so `--mobility-bins 256` is a choice that keeps 256.

**Why part 1.** Before it, the reader kept MS functions only, and a mobility
function is classified `IMS_MS`. So every Waters image recorded with
mobility on was refused with "No MS functions found": not the drift
dimension dropped, the whole image. Measured on an in-house Synapt G2-Si
tissue image. The bundled `MLReader` already exports the calls, which mzmine
binds the same way: `getNumberOfMobilityScans`, `getMobility`,
`getRawMobilityScanDataPoints`.

Measured on the public MassIVE MSV000092638 Synapt XS image (8,475 scans,
MS1) and on in-house Synapt G2-Si runs:

- 200 bins on every run, equal steps from exactly 0 ms (0.045 to 0.216 ms
  per bin across runs; 0.06625 ms on the XS image). The step matches the
  run's "ADC Pusher Frequency" in `_extern.inf`.
- `driftScanCount` on every scan equals the function's bin count.
- Summed over the bins, the raw points reproduce the profile spectrum's ion
  current exactly (ratio 1.0 on every sampled scan). The per-bin centroided
  call is about ten times slower and does not sum to the centroided
  spectrum (up to 15.6 % off), so it is not bound.

**The guard.** Over 91 Synapt runs with mobility on, read one per child
process: the 56 that declare `Acquisition Algorithm: ADC Mode` read cleanly;
the 35 that do not, all older Synapt G1 runs, crashed the library with an
access violation (32) or returned bins that do not sum to the spectrum (3).
Even the summed-spectrum read crashed on most of them. An access violation
ends the Python process with no message, so on such a run the reader leaves
the mobility functions out, decided from the side file before it reads any
spectrum, and refuses the run when nothing else is left. That is no regression: before drift time was
read, mobility functions were never converted, so those runs store what they
stored before, and a run with nothing else was refused then too, as "No MS
functions found". `getMobility` past the last
bin raises a C++ exception inside the library, so it is never asked past
`getNumberOfMobilityScans`.

**Why part 2.** 256 channels over 200 values leave 56 channels holding no
value at all, in a fixed comb (`1111011101111...`): a striped heatmap and 56
empty mobility columns per m/z bin in a grid table. The anchor exists so a
heatmap box indexes grid channels. One channel per bin keeps that, because
the heatmap and the grid take their channels from the same function
(`default_channel_plan`), and it puts each bin in exactly one channel. The
channels are still equal-width, so the grid's `law` stays `"linear"`. No
schema change.

**Why a declaration, not a test on the values.** "Equal steps, fewer than
256" also describes stubs in the TIMS tests and a synthetic TDF fixture.
A real TIMS ramp has 400 to 1,000 scans on a non-linear 1/K0 law. Only the
reader knows whether its values are bins it defines or samples of a
continuum, so only the Waters reader sets `bin_width`.

**What changes in a store.**

- Waters runs with mobility on, which were refused, now convert. They gain
  `uns["mobility_axis"]`, `uns["mobility_heatmap"]` (200 mobility channels)
  and `ion_mobility` in the versioned block.
- `uns["mobility_axis"]` gains the key `bin_width` where a reader sets one
  (Waters only).
- Waters runs without mobility, and every Bruker and imzML store, are
  unchanged: their axes declare no bins, so the plan is the 256 anchor it
  always was.
- The heatmap is on by default, so a Waters mobility run costs one extra
  read of every bin: 200 library calls per scan, about 40 ms per scan on the
  XS image. `--no-mobility-heatmap` skips it. A grid table costs two more
  passes; Waters has no single-read frame record as TDF has (D5).

**Objections that did not win.**

- "Centroid each bin to match the centroided table." It is ten times
  slower, it crashed on G1 runs, and centroids of bins still do not sum to
  the centroid of the sum. The store records the gap instead
  (`current_ratio`, with a warning naming `--waters-spectrum profile`), as
  it does for TDF `vendor_centroid` (D1).
- "Keep 256 everywhere; consumers may assume it." The docs said `k = 256`.
  A consumer that hard-codes 256 instead of reading `mobility_edges` breaks
  on Waters stores. That is the consumer's fix (the format comes first),
  and the docs now say to read `k`.
- "Probe each file in a child process instead of reading `_extern.inf`."
  It would catch crashes the field does not predict, at a process start per
  conversion. The field separated all 91 runs; a probe can be added if a
  counter-example appears.

**Known limits.**

- `getMobility` takes no function. The reader refuses mobility when two
  mobility functions report different bin counts. A run whose functions
  share a count but not a pusher period would be described by one axis;
  none has been seen.
- Functions with and without drift bins on the same pixels (a parallel
  acquisition): the ones without are converted, as before drift time was
  read, and the mobility functions are listed under `excluded_functions`.
  Functions with and without drift bins on different pixels (chunks of
  one raster) are all converted, and no mobility is recorded, since the
  bins would not cover the image. No measured run has either layout.
- Drift time is not converted to CCS. That needs the run's calibration
  (`getDriftTime_CCS` in `MassLynxRaw`), which is not bound.
- No public DESI or cyclic IMS image with its `.raw` folder was found, so
  neither is tested.

---

## D54. pyright checks types in the lint lane, next to mypy

**Status:** Implemented (2026-10-15), issue #434.

**Decision.** The `lint` job runs pyright 1.1.414 in standard mode as a
pre-commit hook, beside mypy. Its settings are `[tool.pyright]` in
`pyproject.toml`: `thyra/` minus the generated ontology tables, Python 3.12,
platform Windows, and the project's `.venv` for third-party types. A pyright
error fails `lint`, so by D19 it can veto a release. To give pyright those
types, the job now installs the project's runtime dependencies (not the
project itself) as well as the `dev` group.

**Why add it.** As configured here, mypy sees less than it appears to. The
hook runs it with `--ignore-missing-imports` and `--no-strict-optional`, and
`[tool.mypy]` makes pyimzml, spatialdata, anndata, zarr, dask and pandas
`Any`. On `08f0ea4`, with the mypy hook green, pyright in standard mode
reported 276 errors. A triage found no reachable bug among them, but some
were real defects:

- `@register_reader` returned an unannotated class, so every reader class
  read as the base class and its own attributes went unchecked (#483).
- Untyped abstract copies in the Bruker base reader hid twelve override
  mismatches (#484).
- zarr nodes were used as groups unchecked; a malformed store crashed with
  an unrelated `AttributeError` or `IndexError` (#485).
- A closed or metadata-only timsTOF reader failed on `None` with an
  `AttributeError` that named neither cause (#489).

By the release of 4.5.0 (`7845456`) the count was 279. #483 to #490 took it
to 0 on `cb060e8`. The gate keeps it there.

**Why mypy stays.** The two overlap but do not nest. mypy's settings here
report what pyright's standard mode does not: unused `# type: ignore`
comments (`warn_unused_ignores`), unreachable code (`warn_unreachable`),
returning `Any` from a typed function (`warn_return_any`) and comparisons
that can never be equal (`strict_equality`). Dropping either checker drops
checks the other does not make. Running both costs about ten seconds of
pyright per lint run.

**What would retire mypy.** pyright's own switches for those checks
(`reportUnnecessaryTypeIgnoreComment` and its reachability and `Any`
reporting, or strict mode) turned on and at 0 errors, then mypy's hook
removed in a new D-entry. Not before: the lint job would quietly lose
checks.

**The platform pin.** `pythonPlatform = "Windows"`, for the reason
`[tool.mypy] platform = "win32"` exists. Measured on the same tree, `Linux`
and `Darwin` each report four errors, all of them `winreg` attributes and
`ctypes.windll`, which exist only on Windows, behind guards that are correct
on every platform. Without the pin the verdict depends on the runner.

**How the hook is wired.** Two settings are load-bearing, and a test guards
both:

- `pass_filenames: false`. Given file names, pyright checks exactly those
  and ignores `include`, so `--all-files` would hand it tests and scripts.
  `types_or` adds `toml`, so an edit to `pyproject.toml` alone runs it too.
- No args. Settings live in `pyproject.toml`, so the hook and a bare
  `pyright` run agree, as flake8's do with `.flake8`.

**Alternatives.**

- "Replace mypy now." It drops the four checks above.
- "Run pyright as its own CI step, not a hook." That is a second command to
  keep in step with what a contributor runs locally. D19 chose one command,
  `pre-commit run --all-files`, for exactly that reason.
- "Basic mode." Fewer checks for no gain: standard is pyright's default,
  and it is the mode the 279 errors were counted and fixed in.

**Costs, accepted.**

- The lint job installs spatialdata, dask, zarr and the rest of the runtime
  tree. It is a larger install than `--only-group dev` was.
- A local run needs `uv sync` first. Without `.venv` the hook fails on
  every third-party import. A bare `pyright` instead uses whatever `python`
  is on PATH, so its verdict then depends on the machine.
- The hook env downloads the pyright package from npm when it is built.
  Like D19's GitHub dependency, an outage there can turn `lint` red.
- The `rev` is not in `uv.lock`, so nothing moves it; bump it by hand,
  re-measure, and fix what a newer pyright reports in the same PR.

---

## D55. Non-finite values decide nothing, and a mislabelled polarity term is read on spectra

**Status:** Implemented (2026-10-16). Amends D34 and D35.

**Decision.** Four rules, each about a value that is not what it claims to
be.

- **The intensity threshold sums measurements only.** A point whose m/z is
  not finite, or whose intensity is not finite or is negative, is in no
  peak's sum. Its own mask entry stays true: the converter drops, counts or
  refuses it exactly as it does with no threshold set. Every reader and
  both mobility passes call the one helper, `summed_peak_mask`. The
  timsTOF mobility and window tables test summed digitizer counts
  instead, which cannot be NaN or negative. A spectrum of only finite,
  non-negative values takes the helper's first branch, which is the old
  body unchanged, so its store is byte-identical.
- **A non-finite m/z is dropped before resampling, by either method.** The
  converter drops it with the unusable intensities and counts it in the
  same warning, on a resampled axis. `nearest_neighbor` already dropped it
  as out of range, but `tic_preserving` interpolated over it: an infinity
  gave the spectrum no bins at all, and a NaN made every bin NaN. On the
  raw axis it is `var["mz"]` itself, and the axis check refuses it.
- **The imzML source range is taken over finite m/z.** A spectrum's NaN or
  infinite m/z adds nothing to it, and a resolved range that is still not
  finite (an API bound of `inf`, or another reader's range) is refused in
  one line, as an empty range is.
- **`MS:1000128` named "positive scan" states positive polarity on a
  spectrum too.** pyimzml already corrects that pair to `MS:1000130` in
  `fileContent` and the parameter groups. Thyra applies the same correction
  to the per-spectrum terms, and to the first spectrum on the metadata-only
  route. The name decides it, so a true "profile spectrum" term is untouched.
  The spectrum-type fallback no longer reads the mislabelled term as profile.

**Why.** The converter drops a NaN, infinite or negative intensity before
resampling (`_drop_unusable_intensities`), but the threshold summed raw
values first. So one unusable copy of a peak decided it. A NaN copy removed
a peak of 100 that the store kept without a threshold; an infinite copy kept
a peak whose only usable share was 10, under a threshold of 50. The grouping
also put an infinite m/z into the largest finite peak, because the
split test `step > 1e-12 * m/z` is false for infinity.

The source range took the minimum and maximum of every value. One infinite
m/z made it infinite, and the default conversion ended in `OverflowError:
cannot convert float infinity to integer`. A NaN made its spectrum's extrema
NaN, which the running minimum then skipped, finite peaks included, so the
range could come out narrower than the data.

The polarity term was matched by accession alone on a spectrum. A file that
stated the mislabelled pair on every spectrum stored no polarity, while the
same pair in a group gave "positive". The same statement should give the
same value wherever the file puts it.

**What changes in a store.**

- Stores converted with `--intensity-threshold` from sources whose repeated
  m/z carry a NaN, infinite or negative intensity, or that hold an infinite
  m/z: they keep or drop different peaks. Stores of sources with only
  finite, non-negative values do not change.
- `--no-resample` with `--intensity-threshold` on a source with a NaN m/z
  is now refused, as it is without a threshold. Before, the threshold
  dropped the NaN point whenever that point's group did not reach the
  threshold (a NaN sum never does), and the conversion then succeeded.
- Processed imzML with a NaN m/z: `essential_metadata.mass_range`, and so
  the default axis, now covers every finite peak. Such a file converted on
  main too, on a range that could miss the spectra holding a NaN. Files
  with an infinite m/z, and continuous files with a non-finite shared m/z,
  now convert by default instead of failing. With `--no-resample` they are
  still refused, since those values would become `var["mz"]`.
- Resampled with `tic_preserving` (the default for some profile sources),
  a spectrum with a non-finite m/z keeps its finite peaks, where an
  infinity stored it empty and a NaN stored it as NaN.
- imzML files that state "positive scan" under `MS:1000128` on each
  spectrum gain `polarity` and `polarity_term` (as `MS:1000130`), and the
  root `acquisition_parameters` polarity. When `fileContent` states no
  spectrum type, such a file is no longer stored as profile.

No schema change: every field already exists.

**Known limits.**

- An mzPeak archive of a Bruker run sums raw counts per TDF-scale bin
  (`_on_tdf_scale`) before the threshold helper runs, so a NaN count still
  reaches that bin's sum there, and an infinite m/z still joins the largest
  finite bin. "Every reader calls the one helper" holds for the threshold
  test, not for that earlier grouping. PR #514 fixes it.
- When the threshold removes every measurement and only unusable points
  reach the converter, the empty-store refusal blames the unusable values
  and does not name the threshold.

**Costs, accepted.** The term correction runs one more XML lookup on each
spectrum that states no correct polarity term itself. Measured on a 60,000
spectrum file with the polarity in a group: parse 3.20 s to 3.58 s (median
of five). The hook is pyimzml's private per-spectrum step; a test converts a
file that states the pair on each spectrum, so a pyimzml release that
renames the step fails there.

**Alternatives.**

- "Drop unusable points in the threshold's mask." Then `--no-resample` with
  a threshold would silently drop an infinite m/z that it refuses without
  one. Leaving them to the converter keeps one rule for them.
- "Refuse a source with any non-finite m/z." A resampled axis has no use
  for those values, and dropping them as out of range is what happens to
  any peak outside the axis.
- "Add `MS:1000128` to the polarity accessions." Every "profile spectrum"
  term would then read as positive.

---

## D56. The resampling step records the axis it describes, and refusals name the way forward

**Status:** Implemented (2026-10-16). Amends D43.

**Decision.** Four rules, so that what a conversion says matches what it
did.

- **`bins_per_fwhm` is honoured or refused.** It is refused beside
  `width_at_mz` or `target_bins`, which each fix the bin count of a `tof`
  axis too, and on an axis that resolves to any other type. These are the
  rules `plan_mass_axis` already applied. The option exists only in the
  Python API.
- **The step records the width the axis has.** When a width set the bin
  count, `mass_width_da` is that width, as before. When `target_bins` set
  it, or the floor of 100 bins raised it, `mass_width_da` is the bin width
  at `reference_mz` that the count gives, and on a `tof` axis so is
  `bins_per_fwhm`. A width law (`tof_a`, `tof_b`, `bins_per_fwhm`) is
  recorded on a `tof` axis only. A plan keeps its width as its setting, and
  its origins and warnings state the width it lays.
- **An empty range from the source names the way through.** When both
  bounds came from the source and its range is a single m/z, the refusal
  names `--resample-min-mz` / `--resample-max-mz` and, for a source whose
  peaks all sit at that m/z, `--no-resample`. A plan of one such source
  says its range is a single m/z.
- **A Waters range widened past the acquisition range is reported where it
  sets the axis.** The extractor logs the widening at debug. The axis
  planner warns, and only for a bound the source's range set: typed bounds
  and a plan give no such warning.

**Why.** Each case was a statement that did not hold.

- `bins_per_fwhm` beside a width was dropped without a word. Beside
  `target_bins` it was recorded as applied, though the count fixed the
  axis. On a `constant` axis it changed nothing and was still recorded.
  `plan_mass_axis` refused the same requests, so a plan and a conversion
  disagreed on what is allowed.
- The step took `mass_width_da` from the request or the default, whatever
  fixed the axis. A width of 1.0 Da over 500-510 m/z recorded 1 for an axis
  of 100 points 0.101 Da apart. `target_bins` 100 on a `constant` axis over
  100-1000 recorded the 5 mDa default for points 9.09 Da apart.
- A source whose every peak sat at m/z 500 was refused as "[500, 500] is
  empty" with no remedy. `--no-resample` and typed bounds both convert it.
- The Waters extractor warned that "the resampled axis is widened" while it
  read the metadata, before typed bounds or a plan applied. With
  `--resample-min-mz 700 --resample-max-mz 900` the axis was 700-900 and
  the next warning said those peaks were dropped.

**Why the recorded width is the count's in two cases only.** Every
generator lays its bins from the range and the count alone, so the count
fixes the width. When a width set the count on a physics axis, the width the
axis has differs from it only by the truncation to a whole number of bins.
Recording the request there keeps every other store unchanged. Where the
count came from elsewhere, the request says nothing about the axis.

A `constant` axis is the exception. A width `w` sets `n = (max - min) / w`,
but the generator lays `n` points, so the step is `(max - min) / (n - 1)`,
wider by `n / (n - 1)`. At the 5 mDa default over 900 Da that is about 6 parts in
a million; at the 100-bin floor it is 1%. The step records the request
there as before, so the same `constant` axis records `w` when a width laid
it and the laid step when a count or the floor did. Recording the laid step
always would change every `constant` store's step for a difference no
reader acts on. Whether the floor applied is read from the count the width
lays, never by comparing widths, since the laid step can come out wider
than the width asked for.

**Why refuse `bins_per_fwhm` off a `tof` axis, when `--tof-law` only warns.**
`bins_per_fwhm` has no CLI flag, and `plan_mass_axis` already refuses it.
`--tof-law` is a CLI flag whose warning existing commands rely on; that
warning stays, and its values are no longer recorded.

**What changes in a store.** Only the `mass axis resampling` processing
step in `msi_metadata`, for conversions that used `target_bins`, a width
the floor raised, or a width law on a non-`tof` axis. The data, `var` and
every other store are unchanged. No schema change: the step's parameters
are free-form scalars.

**What changes for callers.** API calls that set `bins_per_fwhm` beside
`width_at_mz` or `target_bins`, or on a non-`tof` axis, are refused where
they were accepted.

**Alternatives.**

- "Leave the fields out when the count fixed the axis." The step merges
  the request under the resolved values, so a left-out field would fall
  back to the request, the value this fixes.
- "Always record the width the axis has." Every store's step would change
  by the truncation, for no reader's benefit.
- "Warn on `bins_per_fwhm` off a `tof` axis, as for `--tof-law`." A plan
  would refuse what a conversion accepts.

---

## D57. A store keeps its source's regions, and a region store describes only its own

**Status:** Implemented (2026-10-16), issue #496. Superseded in part by
D83 (2026-10-22): a Waters run with a `MultiRegions` table in its method
file keeps those regions; "Region 1 stays only for ... Waters" now holds
only for a Waters run that declares none.

**Decision.** A store describes the regions it holds, with the numbers and
names its source states. Under `region=`, everything the timsTOF reader
reports is about the selected region. Five rules follow.

- A reader that knows its regions gives a region map and a region summary
  for one region too: timsTOF (the DB RegionNumber), rapifleX (the poslog's
  `R{n}`), mzPeak (the region each scan names, or its one listing) and
  solariX, which already did. A single-region store then keeps the source's
  number, name and box instead of the placeholder region 1. This supersedes
  D44's "As on timsTOF, one region gives no map" and its "A single-region
  Rapiflex store is unchanged". Region 1 stays only for formats with no
  region concept (imzML, PHI, Waters). On timsTOF the name of a lone region
  is the `.mis` Area the alignment pairs it with, which on a run cut from a
  multi-region slide is the one its spot names give, while the number stays
  the database's, the one `region=` accepts. `region=` accepts that name too,
  and refuses another Area's name on such a run.

- `get_region_info`, and with it `uns["regions"]`, lists only the selected
  region. Its count and box are that region's, inside the stored raster.
- The precursor schedule is read from the selected region's frames: the frame
  types in `Frames`, and the windows in `PasefFrameMsMsInfo`,
  `DiaFrameMsMsInfo` and `FrameMsMsInfo`. The mixed-level check already
  filtered this way. Without `region=` the queries are unchanged.
- A position measured in two regions is one row, as any repeated position is
  (D40). On timsTOF and solariX that row takes the region of its first frame
  that holds peaks, the frame whose number is its `acquisition_order` (D26),
  and on mzPeak the region of its first scan. An empty frame gives no row, so
  it names no region. The reader warns once, naming the regions and the first
  shared position.
- `region` takes any integer type except `bool`, so a value read back from
  `obs["region_number"]`, which is a numpy integer, selects its region.

`uns["average_spectrum_per_region"]` is written when the region map holds two
regions or more: one region's mean is the table's own `average_spectrum`. A
region with no row has no entry in it, on every reader. D26's note on column
types is corrected in place, with a dated note: in a 2D store `x`, `y` and
`region_number` are `int32`, not `int64`.

**Why.** A single-region run was stored as region 1 with no name, whatever
the source said: a timsTOF run with one region numbered 3, or an mzPeak
archive naming region 3 "tumour", came out as region 1, while the optical
crop record of the same store named region 0 and Area "01". The placeholder
also broke the region loop: `list_regions` answered 1 for a single-region
`.d` whose region is 0, and `region=1` was then refused. Runs of one slide
merged on `region_number` could collide, because the placeholder 1 equals
another run's real region 1.

A `region=` store listed every region of the run, with whole-run
counts and boxes outside its raster, some negative. A region on one constant
schedule lost its MS/MS table when another region was survey only, because
the survey frames of the other region made the file "mixed". A survey region
was described as MS level 2. A shared position took the region of its last
frame and the order of its first, so one row named two different frames. A
region number read back from a store was refused with a message about
positional arguments nobody passed. A region with no row got a mean of zeros,
which is not a mean.

**What changes in a store.**

- A single-region timsTOF store: `obs["region_number"]` changes from 1 to the
  DB RegionNumber, which was 0 on all five single-region runs measured
  (2026-10-16), and the `uns["regions"]` entry
  gains `bounds` and, where the `.mis` has it, `name`. The same holds for a
  single-region rapifleX store (the poslog's number, usually 0) and a
  single-region mzPeak archive of a Bruker run.
- A single-region solariX store and a `region=` timsTOF store lose
  `uns["average_spectrum_per_region"]`, which held a copy of
  `average_spectrum`. A multi-region store keeps it, even when only one of
  its regions holds rows.
- A `region=` timsTOF store: `uns["regions"]` holds one entry. The
  `fragmentation` block and whether an `_msms` table is written follow the
  region's own frames.
- A multi-region timsTOF or solariX store with a shared position: that row's
  `region_number` is the region of its first frame (it was the last), and the
  per-region means move with it.
- A store with a region that has no row: that region's
  `average_spectrum_per_region` entry is gone.

**Alternatives.**

- "Write no region map when two regions share a position", as the rapifleX
  reader does (D44). Rejected for timsTOF and solariX: their region number is
  a column of the vendor database, so one shared position does not cast doubt
  on the others. The rapifleX map is rebuilt from spot names, and a conflict
  there can mean the keying itself is wrong.
- "Keep each region's spectrum in its own row." The store has one row per
  position (D40); two rows at one `(x, y)` would break every consumer that
  indexes by position. `region=` gives each region its own measurement.
- "An all-NaN mean for an empty region." A consumer that plots every entry
  then draws nothing under a legend label. No entry says the same thing.

**Known limits.**

- A consumer that picked a single-region store's rows with
  `region_number == 1` now finds none on timsTOF, rapifleX and mzPeak.
- On a timsTOF run cut from a multi-region slide (RegionNumber 0, spot names
  `R02...`), the optical crop record names region 2 while
  `obs["region_number"]` is 0. The name agrees; the alignment's number is
  not part of this decision.
- A rapifleX run whose poslog misses a spectrum's spot, and an older mzPeak
  archive whose one listed box does not hold every pixel, still get the
  placeholder 1, with a warning, new for a single-region run: before, one
  region was not checked at all.
- In the full store, a shared position's row sums the spectra of both
  regions, so the earlier region's mean is not the one its `region=` store
  gives. The warning says so.
- The schedule filter needs `MaldiFrameInfo`. A `region=` request on a file
  without it is already refused (issue #252).
- A frame that lists peaks but loses all of them to `--intensity-threshold`
  still counts as holding peaks, so on a shared position it can name the
  region while a later frame gives the row its `acquisition_order`.

---

## D58. The mass-mobility heatmap is a mean over the summed table's rows

**Status:** Implemented (2026-10-16), issue #495. Narrows D38's divisor and
lifts its first known limit.

**Decision.** A spectrum the summed table gives no row adds neither points
nor a pixel to `uns["mobility_heatmap"]`. On the fused route (a TDF), a
frame feeds the heatmap only when it gets a row. On the separate raw pass
(imzML, Waters), the converter tells the heatmap which positions have a row
once the passes have decided them. The heatmap's population is then the one
`uns["average_spectrum"]` divides by, as D38 meant it to be.

Three smaller rules come with it:

- The shared-axis mobility table (a continuous imzML with a mobility array)
  takes the grid route's memory guard (D4) before it builds anything, on
  the distinct `(m/z, mobility)` pairs the list holds, which are its
  `var`. Past
  half of free memory the table is left out with a `WARNING` naming
  `--no-mobility-table`, once per table; past a quarter it warns and is
  built.
- A mobility spectrum with no row counts as skipped only when it is on the
  table's own plane and carries signal, on the separate pass and the fused
  one alike. Another plane's spectrum has its row in that plane's table,
  and an all-zero one loses nothing.
- A 2D conversion of a multi-plane source plans a table for each plane the
  pre-scan meets, not for every z in the span. A single plane and a volume
  are planned up front, as before.

**Why.** On a continuous imzML with one all-zero spectrum, `--no-resample`
or not, the summed table had three rows and the heatmap four pixels:
`current_ratio` 0.75 and a warning, though no current was lost. D38 kept
the heatmap's own count of distinct positions and rejected taking the
table's rows, because one raw pass served every plane and the fallback has
no table. Both still hold, so the rows are not handed over as a number: the
heatmap asks, per position, whether the table has a row there. The
fallback, with no table, still counts every position.

The other three were seen in the same run. A feature list the grid guard
would refuse went through, and could end in a raw allocation error.
Each plane of a two-plane conversion warned that the other
plane's spectra "had no row ... and were skipped". A source stating z 1
and z 120 planned 120 tables: 118 warnings, and 14.4 s against 4.5 s for
the same six spectra on z 1 and 2.

**What changes in a store.**

- A store with a mobility heatmap and a spectrum that gets no row (all
  zero, every point outside the kept mass range, off the grid, or a TDF
  frame whose vendor centroid is empty while its raw scans are not):
  `counts` changes from (all points) / (positions fed) to (points of the
  spectra with a row) / (rows), and `current_ratio` with it. On the issue's
  file, whose rowless spectrum is all zero, that is a factor of 4/3, and
  `current_ratio` goes from 0.75 to 1.0.
- Every store without a rowless spectrum is byte-identical: measured on the
  eight default fixtures and two TDF runs, before and after.
- A store whose shared-axis feature list is too wide for the machine has no
  `_mobility` table, and `uns` does not name one.
- A TDF store under `scan_sum`, with the default mass range and no
  intensity threshold, does not change: every frame with raw points then has
  a summed spectrum, so every one already had a row. A narrowed range or a
  threshold that empties a frame changes it as the first bullet says.
- Planning tables per plane met changes no store: the planes the pre-scan
  never met had no rows, so they wrote no table before either.

**Alternatives.**

- "Divide by the table's row count." The separate pass builds every plane's
  heatmap at once (D41), so it would need every plane's count, and an
  off-grid spectrum's points would stay in the numerator.
- "Skip all-zero spectra in the imzML reader." The reader belongs to every
  route, and a spectrum off the grid or outside the mass range would still
  be counted. The converter is where rows are decided.
- "One skip warning per conversion." Counting only the table's own plane
  already counts each lost spectrum once; each plane's warning then names
  its own.

**Known limits.**

- On the separate pass, a repeated position whose second spectrum is empty
  in the summed table still adds that spectrum's raw points: the pass
  knows positions, not which spectrum of a position gave the row. The fused
  route leaves those points out.
- `VAR_BYTES_PER_FEATURE` was measured on a grid `var`. The shared-axis
  `var` is built differently and its cost per feature is not measured.
- A plane the source states but whose every spectrum is empty gets no table
  and no warning, since the pre-scan never meets it. A spectrum that turns up
  on such a plane in pass 2 is refused as a disagreement between the passes,
  as any other is.

---

## D59. A FlexImaging run that stopped early is placed by its planned grid

**Status:** Implemented (2026-10-16), issue #497. Amends
[D31](#d31-fleximaging-spots-are-placed-on-their-lattice-not-on-the-drawn-area).
D69 extends its overview-photo rule to every optical image with no
alignment.

**Decision.** D31 asks which raster node sits on the `.mis` reference
point: the one that puts every measured spot inside its Area, then the one
that leaves the fewest Area nodes unmeasured. Two things change.

- The search tries every node that keeps each region's spots inside its
  Area's box. It used to try the 3 nodes either side of a guess made from
  the spots' centre and the Area's centre.
- When the source states the grid it planned, that grid chooses among the
  nodes that keep every spot inside: the node that centres it on the
  Areas. timsTOF states it as the `ImagingArea` bounds in GlobalMetadata,
  rapifleX in the `.dat` header. A grid is used only when it holds every
  spot and its size matches the Areas' extent within 1.5 steps. With no
  usable planned grid, D31's rule chooses.
  When several nodes still fit equally, the fit keeps the one nearest the
  old guess and warns that the placement is not determined.

Three related fixes from the same issue:

- **The overview photo.** A second optical image is placed by a `Scale`
  onto the alignment image only when it is a resampled copy of it: the
  same ratio on both axes, within the rounding of the copy's sides. Any
  other image keeps one
  coordinate system of its own, named after the element, with `Identity`,
  and is not in `"global"`, the per-dataset system or `"acquisition"`. A
  warning says so.
- **The Area's edge tolerance** is at least 0.75 photo pixel, besides 0.1
  raster step. A `.mis` stores Area corners in whole pixels, and rounding
  moves an edge by up to 0.71 pixel. The pixel floor stops at 0.4 step, so
  on a raster finer than 1.9 pixels per step a whole-step shift still
  leaves a spot outside.
- **A `<ReferencePoint>` or Area `<Raster>` that is not a finite number**
  (`inf`, `nan`) is unreadable: the lattice uses the first teaching point,
  or the Area's own step is dropped.

**Why.** D31 assumed the spots fill their Area. A run that stops early
fills its first rows, and every whole-step shift along the empty part
leaves the same nodes unmeasured.

| Run | Acquired | Before | After |
|---|---|---|---|
| `rapiflex_dat_profile_a` | 6 of 118 planned rows | 112 nodes tie; kept one 59 rows low | planned grid: one node, its centre 0.4 step from the Area's |
| timsTOF TSF, same slide as `timstof_tsf_centroid_a` | 273 of 2,448 planned spots | kept a node 24 rows off | the node the complete runs on that slide use |
| timsTOF TSF, same slide | 12,867 of 29,939 planned spots | kept a node 48 rows off | the same node |
| timsTOF, 30 um, one freehand Area | 71,249 spots, rows 179-373 of 179-657 | the only fitting node lay 14 x 38 nodes from the guess; kept one with 1,550 spots outside | that node, every spot inside, found with or without the planned grid |

The three timsTOF runs also show what the `ImagingArea` bounds are: they
span every planned row while the frames cover only the first ones, and
51 x 48 nodes is exactly the 2,448 spots FlexImaging planned.

The planned-grid rule assumes that on the right node the grid's centre
lies within half a step of the centre of the Areas' box. A whole-step
shift then moves it more than half a step away. On `rapiflex_dat_profile_a`
the two centres were 0.4 step apart.

On six complete runs the affine is unchanged to the last digit:
`tims_msms_pos_brain1` and two runs beside it, `timstof_tdf_centroid_a`,
`timstof_tsf_centroid_a` and one more complete TSF run. The two public
MassIVE runs still land on flexImaging's spot list.

The overview photo of `tims_msms_pos_brain1` is the alignment scan turned
by about 0.47 degrees. Its two scales were 4.012 and 4.100, and the `Scale`
left it up to 297 photo pixels off at the corners. On
`rapiflex_dat_profile_a` the scales were 2.013 and 2.096, a turn of about
1 degree. The quarter-size `_0001` scan is 4.000 on both axes and is still
scaled.

On `timstof_tsf_centroid_a` the old tolerance flagged 28 of 1,599 spots as
outside the Area, the furthest 0.38 pixel out. Now none are, and the
affine is unchanged. So a complete run that used to warn about spots
outside its Area, from corner rounding alone, no longer does.

**What changes in a store.**

- A FlexImaging run that did not fill its Area: the TIC image's transform,
  the pixel polygons, the crop window and `raster_to_global_affine`. With
  `apply_optical_alignment=False`, the photo's place in micrometres.
- A second optical image that is not a resampled copy of the alignment
  image, such as the FlexImaging overview photo: its transformations. An
  mzPeak archive's extra images follow the same rule.
- A non-finite reference point is no longer copied into the raw `.mis`
  metadata.
- A complete run whose Area holds its spots one way only: nothing.

**The objection.** Register the overview photo onto the alignment scan
from its pixels, so it stays placed. It did not win. It means decoding both
photos before the store is declared, against the bounded memory the
optical path keeps. A feature fit is random, so two conversions would
differ. It would need scikit-image as a declared dependency. And the
overview is a coarser copy of the alignment scan, which is placed.

A second objection: with no planned grid, anchor the run at the Area's
first node, where FlexImaging starts. It did not win either. On
`rapiflex_dat_profile_a` the first planned column was never acquired, so
that anchor is a step off.

**Known limits.**

- The planned grid is taken to cover every Area in the `.mis`; a grid of
  another size is not used. Every partly acquired run measured had one
  Area.
- The size check allows 1.5 steps. A grid planned with a margin column on
  both sides of the Area fails it, and the fit falls back to D31's rule and
  its tie warning: safe, but the grid goes unused.
- Where the Areas' centre sits exactly half a step from both of two nodes,
  the planned grid ties too, and the fit warns.
- A square overview photo has the same ratio on both axes at any turn, and
  is still scaled.
- The search tries at most 250,000 nodes per pair of axis directions. Past
  that, or when no node keeps the spots inside the Areas' box, it tries the
  3 nodes either side of a seed: the node that centres the planned grid on
  the Areas when there is one. Without one the seed is the old guess, as
  before, which on a run that stopped early in a large Area can miss the
  right node; the fit then warns that spots fall outside, or that the
  placement is not determined. When no node keeps every spot inside, a
  first pass over the edge spots may miss the node with the fewest spots
  outside; the fit warns either way.
- Below 1.9 photo pixels per step, whole-pixel corners can still trip the
  warning about spots outside the Area.
- A planned grid cannot tell a mirrored stage from the usual one. With one,
  the usual axis directions are preferred before the grid is consulted.

---

## D60. A damaged isolation window is left out, and an isomer's label is its rank

**Status:** Implemented (2026-10-16), issue #498.

**Decision.** Three rules for the MS/MS schedule.

- The timsTOF reader screens every isolation window it reads, from
  `PasefFrameMsMsInfo`, the DIA tables or `FrameMsMsInfo`. The target must
  be a finite m/z above 0; a NULL target (which is also how SQLite stores a
  NaN) is no target. The width, when the row has one, must be finite and
  not negative. The collision energy must be finite, of either sign, since
  the schema sets it no bound and signed energies exist. A PASEF or DIA
  scan range must ascend from 0. A row that breaks one of these is left out
  of the schedule, with one warning that counts them. The table the frames
  use is still the first one with rows (D37): when none of its rows is
  usable, the schedule records the MS level with no windows, and the next
  table is not asked. Isomer windows are listed by m/z, then by first scan
  descending, which is 1/K0 ascending on the TIMS ramp: the order the
  `_msms` table ranks them in.
- A window left out that still occupies a usable scan range stops the
  precursor split: no `_msms` table is written, with a warning. Its scans
  still hold the ions it isolated, and splitting by the windows that are
  left would put those ions in no precursor. This qualifies D2: the split
  is lossless only when every window of the schedule is kept, so a
  schedule the screen has cut is refused like any other that does not
  qualify. A left-out window with no usable range occupies no known scans,
  and the split goes ahead.
- The metadata builder screens each window it is given the same way for
  the target, the offsets and the energy, and skips a bad one instead of
  failing the whole `msi_metadata` block. A bad scan range there drops the
  range from the window, as it did before, not the window. The Waters
  reader no longer reports an infinite precursor or energy.
- In an `_msms` table, two windows at one m/z (an isomer pair) are told
  apart by their rank in mobility among the windows of that m/z, over the
  whole schedule. Every column of the second isomer ends in `_1`, the third
  in `_2`, and so on. The first isomer has no suffix.

**Why.** The schedule was stored as read. On a damaged `analysis.tdf`, an
infinite target became a `precursor_mz` that `thyra validate` refuses, a
negative target was stored as a precursor, an infinite energy wrote the
token `Infinity` into the document, which strict JSON parsers reject, and
one negative width made the builder raise. The assembler caught that and
logged it, and the store had no `msi_metadata` at all, at exit 0.

The suffix used to come from a collision: only a mass bin that both
isomers held gave the second one `_1`. In a bin only the second isomer
held, its column had no suffix and took the first isomer's name. So
`p500_mz2209` named the isomer at 1/K0 1.092 in one store and the one at
1.821 in another store of the same schedule. A join on the label merged
two different isomers, and `anndata.concat(join="inner")` drops
`precursor_mobility`, so nothing showed it. The docstring already said the
suffix followed mobility order; now it does.

**What changes in a store.**

- A store whose source has a bad window: the window is gone from
  `ms_analysis.fragmentation.windows`, and the `msi_metadata` block is there
  where it was missing. Where the window still occupies its scans, the store
  has no `_msms` table; where it does not, the table holds the other
  windows.
- An `_msms` table with an isomer pair: the second isomer's columns in bins
  it alone holds gain `_1`. Columns in bins both hold, and every column of a
  schedule without isomers, keep their names. With three or more windows at
  one m/z, a bin held by the first and third gives the third `_2`, where it
  was `_1`.
- A schedule with an isomer pair: `uns["msms_schedule"]` and
  `fragmentation.windows` list the pair by 1/K0 ascending (first scan
  descending), where SQLite's grouping order (width, energy, then first scan
  ascending) listed them before. Window 0 of the schedule is now the isomer
  without a suffix.
- No `msi_metadata` field changes shape, so the schema version stays.

**Alternatives.**

- "Refuse the source." One damaged window out of a schedule is not a reason
  to lose the conversion, and the warning names the count.
- "Clamp a negative width to 0." That invents a value the source does not
  state.
- "Suffix by the isomer's scan range." A rank needs no other value in the
  label and is the order the docstring promised.

**Known limits.**

- The screen tests the numbers, not whether they make sense together: a
  window far outside the acquired mass range is kept.
- A left-out window whose scan range is itself unusable gives no scans to
  reason about, so the split goes ahead; if it did isolate ions, they are in
  no precursor column and `demultiplexed_current` falls below 1 with its
  usual warning.
- A value bad in some frames only splits one window into groups, and the
  schedule then reads as varying, so no `_msms` table is written.
- The stem kept six significant digits (`p{mz:g}`), so any two targets
  that rounded alike, at any m/z (798.5401 and 798.5404 were both
  `p798.54`), shared a stem and were ranked as isomers. Since D67 the stem
  writes the m/z exactly, and only targets with the same m/z are ranked.
- Leaving an isomer out (a bad window, or a `region=` store whose region
  lacks it, D57) re-ranks the other one.

---

## D61. An mzPeak TDF frame sums only measurements into a TOF bin

**Status:** Implemented (2026-10-16). Amends D33, and applies D55's rule to
the mzPeak reader's own grouping.

**Decision.** When the mzPeak reader puts a TDF archive's raw counts on the
`.d`'s scale (D33), it sums the points of one frame that fall on one TOF bin.
A point whose m/z is not finite, or whose count is not finite or is
negative, now joins no bin. It is passed on by itself, scaled but not
rounded, and the converter drops, counts or refuses it as it does any other
such point. A NaN m/z was already dropped earlier, as null-pair padding.

**Why.** The grouping split two sorted points when their step exceeded
`1e-12 x m/z`. Relative to an infinite m/z no step is large enough, so an
`inf` point joined the largest finite bin. Its count was added there and
kept at that finite m/z, with or without an intensity threshold. Two `-inf`
points were summed into one. A NaN count made its whole bin NaN, and the
converter then dropped the bin's real signal with it. A negative count
shrank the bin it joined. D55 fixed the same split test in
`summed_peak_mask`; this grouping is its own copy of the rule.

Measured on a two-pixel fixture at 100 ms, so counts keep their values. Its
first pixel holds 2, 3 and 5 at m/z 100, 101 and 102 and one extra point.
Values stored for that pixel, resampled to 99.5-102.5, no threshold. A
converter-level test holds each "Now" row, and the warning that says the
extra point was dropped:

| Extra point | Before | Now |
|---|---|---|
| 7 at m/z `inf` | 2, 3, 12 | 2, 3, 5 |
| NaN on m/z 100 | -, 3, 5 | 2, 3, 5 |
| -4 on m/z 100 (bin of 9) | 5 at m/z 100 | 9 at m/z 100 |

**What changes in a store.** Only stores of TDF-sourced mzPeak archives with
a non-finite m/z, or a non-finite or negative count, on a frame. Their bins
hold the measured counts alone. With `--intensity-threshold`, the points
passed on decide nothing either: the threshold, applied after this step,
leaves them out of every sum (D55). Without resampling, a non-finite m/z is
refused as before, since it would become `var["mz"]`. Archives whose m/z
and counts are all finite and non-negative give the same stores, value for
value. A count between -0.5 and 0 used to round to 0; it is now passed on
unrounded and dropped with a warning. Real TDF counts are unsigned, so no
real archive holds one. No schema change.

**Alternatives.**

- "Drop such points in the reader." Then the converter could not count or
  refuse them, and a non-resampled conversion would quietly lose a point it
  refuses on every other route. D55 rejected the same option.
- "Reuse `summed_peak_mask`." It returns a keep mask, not the summed bins,
  so it cannot replace this step. The two share the tolerance
  (`SAME_PEAK_RELATIVE`) and now the rule for what counts.

## D62. A NaN or an infinity from the vendor metadata is left out of the metadata block

**Status:** Implemented (2026-10-17). Applies D55's rule, that a non-finite
value decides nothing, to the `msi_metadata` block. D60 did the same for
isolation windows. No schema version change.

**Decision.** The block holds finite numbers only. JSON has no NaN or
infinity, so the published JSON Schema never allowed them; the models now
agree with it.

- **The builder leaves such a value out.** The mobility range is kept with
  both edges finite or not at all; its scan count likewise. A processing
  parameter that is not finite is left out of its step. An isolation window
  with one is left out whole, as D60 decided.
- **A pixel size is refused, since the block requires one.** `convert_msi`
  refuses a detected pitch that is not a finite number above 0, as it
  refuses a missing one, and asks for `--pixel-size`. The converter refuses
  such a `pixel_size_um` or detected y pitch. A converter built without a
  pitch keeps its default when the source's is unusable, and says so; the
  store records the pitch as `default`. `thyra metadata` leaves the pitch
  out, as D23 does for a source with no raster.
- **The Python API refuses a non-finite setting, as the CLI already does.**
  That is `intensity_threshold` and every number in `resampling_config`.
  No summed peak reaches a NaN or `+inf` threshold, so either one emptied
  the store while the step recorded the threshold. `-inf` kept everything.
- **The models refuse NaN and infinity** (`allow_inf_nan=False` on every
  schema model), and `to_uns_dict` packs its JSON strings with
  `allow_nan=False`.
- **`thyra validate` reports a stored one as an error,** once per value, at
  its location. The message says that converting the source again leaves the
  field out. The exit status is 1.

The order is the point. The uns assembler catches the models' error, so
refusing at the models alone wrote the store with no `msi_metadata` block
and one warning, while the conversion reported success. A converter test
pins that the block is still written.

**Why.** #494 refused these values where a JSON document is read or
written. A conversion still stored them, because the models took them. A
`gt=0` bound refused NaN but not `+inf`; an unbounded field took all three.
`thyra metadata` and `thyra export-metaspace` then refused to write the
store's own metadata as JSON.

**What changes in a store.** Only stores whose vendor metadata or API
settings held a NaN or an infinity in one of these fields. The field is now
absent instead. A pitch like that now stops the conversion with a message,
where it used to reach the store. The problem was found by a stress test;
no real dataset is known to hold one. Every other store is unchanged.

**Why no schema version change.** The published JSON Schema is the same
file: pydantic writes nothing for `allow_inf_nan`, and a JSON document
cannot hold NaN anyway. Every JSON document that validated still does. Only
a Zarr store, whose attributes can hold a NaN float, could have one.

**Why `thyra validate` errs rather than warns.** Such a block does not
conform to the schema it names, and a strict JSON reader of it fails. A
warning would pass it through a CI gate that exists to catch that. The
remedy is the user's own: convert again.

**Alternatives.**

- "Refuse the conversion when a vendor value is not finite." One optional
  field would then block a conversion whose spectra are fine. D60 took the
  same view for isolation windows. The pitch is the exception, because the
  block cannot be written without one.
- "Store `null`." The block's convention is that a field the source does not
  report is absent, not null.
- "Coerce to a finite stand-in, such as 0." That invents a measurement.

## D63. An imzML states z on every spectrum or on none

**Status:** Implemented (2026-10-19), issue #521. Amends the imzML half of
D30.

**Decision.** Every spectrum's z is judged the way the parser reads it: a
cvParam `IMS:1000052` with a value, on the scan itself. A term with no value
counts as no z.

- A file that states z on some spectra and not on others is refused, with
  the mzPeak route's sentence: "A plane is stated for every spectrum or for
  none."
- A file whose scans take z from a parameter group is refused by name. The
  parser would give those spectra its stand-in 1 and sum the planes.
- Every spectrum states z: D14 applies as before. None does: D30 applies,
  and the z offset is 0.

The judgement costs no second pass. A hook on the parser's own pass over the
spectra records one byte per spectrum, as the polarity hook of D55 does.
D30's sentence about reading the first spectrum again is replaced.

**Why.** Reading only the first spectrum went wrong four ways, all at exit
0. The recorded offset depended on spectrum order. A term with no value
counted as stated. The stand-in 1 formed a plane of its own: 11 rows on one
plane, 1 on another. Two planes stated through groups were summed into one
table. Meanwhile the archive of the same file was refused, or got another
offset.

**The same hook names what the parser cannot read.** A position held in a
parameter group or in a second `<scan>`, a position that is not a whole
number, and an array read through a group the parser did not pick each
ended in a raw `TypeError`, `ValueError` or `AttributeError` traceback. Each
is now one `ConversionRefused` naming the spectrum and the cause. A
document that is not well-formed XML, or is empty, gets the same, naming the
file. The position cases are refused rather than read, as the z case is.

**The objection.** Convert such files with the unstated spectra on plane 0.
It lost: that mixes planes from stated and unstated z, and the two routes
would still differ. A refusal costs only hand-made files. Every public
example and test fixture states z on all spectra or on none, and none
through a group.

**What changes in a store.** Files that mix stated and unstated z, or hold
z in a group, stop converting; no store is written. A file whose z terms
all lack a value records z offset 0 instead of 1. Every other store is
unchanged.

**Known limit.** The hook wraps a private pyimzml method, as D55's does.
The z tests convert files through it, so a pyimzml release that renames the
method fails there.

## D64. imzML arrays: precision on the first spectrum, byte ranges on every spectrum

**Status:** Implemented (2026-10-19), issue #521.

**Decision.** Thyra refuses an imzML whose binary arrays overlap in the
`.ibd`, on any spectrum. Two arrays of the same kind may point at exactly
the same bytes, as a continuous file's shared m/z axis does. Thyra also refuses a file whose
first spectrum has an array whose own precision term disagrees with its
param group, whether or not the array states `IMS:1000104`. An ion mobility
array declared as 16-bit float (`MS:1000520`) is read as 16-bit; a precision
Thyra cannot decode is refused by name.

**Why.** The parser reads each array at its group's precision and never
looks at a term on the array itself. It also reads as many values as a
spectrum states, wherever they end. A spectrum that stated four values and
stored three read into the next array, and an invented m/z of 134,217,760
joined the mass range at exit 0. A 32-bit array on spectrum 0 was stored as
514 instead of 15. A 16-bit mobility array was read as 64-bit, from the
neighbouring bytes.

**Cost.** The overlap check needs only the offsets and lengths the parser
already holds: 0.46 s on a 918,855-spectrum file whose parse takes about
115 s. Checking every array's precision would need a second pass over the
document (73.6 s on that file), and a per-array precision conflict has only
been made by hand.

**What changes in a store.** No store of a well-formed file changes. A file
whose arrays overlap, or whose first spectrum has an array that contradicts
its group, converted with invented or wrong values and is now refused. A
16-bit mobility array now stores its own values.

**Known limits.** A later array that states its own precision, and a later
spectrum that states fewer values than it stores (a gap, not an overlap),
still convert wrong at exit 0. A resampled m/z axis whose source range
starts at 0, on a logarithmic or inverse law, still ends in a raw error in
the axis planner; that belongs to the resampling fix, not this one.

## D65. The imzML head is trusted only when its declarations fit together

**Status:** Implemented (2026-10-19), issue #521. Amends D14.

**Decision.** The preview and the metadata document describe an imzML from
its head only when it declares at least one spectrum, both raster sides are
at least 1, the count is no larger than the raster, and the first spectrum
sits at x >= 1 and y >= 1. Any other file takes the full read, as a
conversion does. A count equal to the raster is a full raster and is
trusted.

Three smaller reads of the head follow the same rule, that the preview
says what a conversion would read:

- **The pixel size is read where its value is.** A size stated through a
  parameter group that the scan settings reference takes its name and unit
  from that group too. The preview reported such an area, or a nanometre
  value, as micrometres.
- **The last array ends at its offset plus its count times its width.** A
  wrong encoded length on a complete file made the preview call it
  truncated, while the conversion read it.
- **The observed m/z terms are matched as XML.** Attributes in any order,
  either quote, space around `=`, and nothing inside a comment.

A pixel size in a unit Thyra cannot convert is not read: the log names the
unit and a conversion asks for `--pixel-size`, as on the mzPeak route. It
used to be refused, so the file could not be converted even with the flag.

**Why.** pyimzml's writer declares the largest coordinate, not the count,
and the head read stops at the first spectrum. So a 0 that first appeared
later went unseen: a 4 x 3 file listed from its last pixel previewed as 3 x
2. A count larger than its raster, or a side below 1, describes no
acquisition, and checking both costs nothing. Reading every position
instead would give up the header-only speed (0.4 ms against 64 s on a 2.0
GiB file).

**What changes in a store.** None from the head rule: no conversion reads
the head. A file with a pixel size in an unknown unit now converts when
`--pixel-size` is given.

**D14, amended.** D14 says only the mzPeak extractor reads `IMS:1000042` /
`IMS:1000043`; the preview reads them too, under this rule. The rebasing
line D14 describes is now logged once per file, naming each axis's base.

**Known limit.** A wrong count that still fits its raster (too low, or a
file cut short after the count was written), and a 0-based file whose
spectra fit inside the declaration with its first spectrum away from 0,
are still described from the head. Both touch only the preview and the
metadata document, never a store.

## D66. A timsTOF store records the calibration Bruker's library applied

**Status:** Implemented (2026-10-19), issue #518. Qualifies D46 for TDF: the
axis `values` are the first frame's, and a frame with a TIMS calibration of
its own is converted with it. No schema version change.

**Decision.** Five rules, each so that the store says what the library did.

- **Each frame takes its own 1/K0.** Bruker's library turns a frame's scan
  numbers into 1/K0 through the `TimsCalibration` row it applies to that
  frame. When the frames take more than one row between them, each frame's
  points are converted with the library's call for that frame, as its m/z
  already are. The mobility table, the mass-mobility heatmap and the
  mobility grid all take these values, and the range they bin over is
  widened to every frame's own first and last scan. When every frame takes
  the same row, the axis is used as before. When which rows apply cannot be
  told (the library applied a state Thyra did not find, or a mapping cannot
  be read), every frame is converted for itself, which is exact whatever
  the rows are.
- **The recorded row is the applied one.** `uns["mobility_axis"]["calibration"]`
  holds the row the library applied to the first frame, the frame the axis
  `values` come from: the applied `calibration.sqlite` state's row, or the
  analysis database's row for that frame. New keys say which: `frame_id`,
  `stored_in` (`calibration.sqlite` or `analysis.tdf`),
  `calibration_state_id` when a state applied, and `n_frame_calibrations`,
  the number of rows the file's frames use. It used to be the analysis
  database's first row whatever was applied. When the applied state is
  unknown, no row is recorded.
- **A state is named as applied only when the library applies it.** The
  library applies the newest state whose `FrameMzCalibrationMapping` (and,
  on a TDF, `TimsCalibrationMapping`) names exactly the acquisition's
  frames. A newer state that maps one frame too few or too many is skipped
  and the one before it applied; with none left, the analysis database's
  calibration is. Thyra now reads the mappings the same way, and checks its
  answer against the library's `has_recalibrated_state` once the file is
  open. When the library applies none, no state is named. When it applies
  one that Thyra did not find, the m/z calibration step leaves
  `calibration` out, as it does for an unreadable file. When a mapping
  cannot be read, the newest state is taken, as before.
- **A missing `CalibrationInfo` costs only the software version,** and a
  `calibration.sqlite` with no state names the analysis database in the
  m/z calibration step.
- **TSF or TDF is decided in one place,** and the metadata takes the
  reader's answer. `analysis.tsf` wins when both databases hold data; a
  0-byte one gives way to the other.

**Why.** Measured with the bundled library on copies of `synthetic_tims.d`.
A state mapping one frame too few, or a frame the acquisition lacks, is
skipped: the m/z were the state before it, or the analysis database's,
while the store named the newest. On tims_msms_pos_brain1 the only state
maps the 713 frames to 300 distinct TIMS rows; a frame's 1/K0 sat up to
9e-4 (about three scans) from the first frame's. Every default conversion
of a real TDF with a `calibration.sqlite` therefore binned later frames on
the first frame's ramp. With `--no-recalibrated` every frame took one row,
and the two stores still recorded the same coefficients. The pea and
Xenium TSF folders hold a 0-byte `analysis.tdf`, so their stores said
`bruker_tdf` with ion mobility present, while the spectra were read from
`analysis.tsf`.

**What changes in a store.**

- A TDF converted with the recalibrated state on, whose state maps frames
  to more than one TIMS row (every real one checked): the mobility table,
  the heatmap, the grid and their bin edges move by up to a few scans'
  width. The summed table does not change.
- Every TDF: `uns["mobility_axis"]["calibration"]` gains the four keys, and
  its coefficients are the applied state's when a state applied.
- A source whose newest state the library skips: `format_specific.calibration`
  and the m/z calibration step name the state applied, or none.
- A TSF folder with a 0-byte `analysis.tdf`: `format_specific` says
  `bruker_tsf`, and `ion_mobility` is `{present: false}`.

The MS/MS table's `precursor_mobility` stays the first frame's 1/K0 at the
window's middle scan: a window has one value, and the first frame's stands
for the others, which differ by at most a few scans' width.

**Alternatives.**

- "Convert every frame through its own call, whatever the mappings." Exact,
  but it costs a call per frame for files whose frames all share one row,
  and changes nothing there.
- "Record only the coefficients, and leave the axis as it is." The F36
  half alone; the frames' points would still be binned on the wrong ramp.
- "Name the newest state, and say the library may have skipped it." A
  store should say what was applied; the library answers it.

## D67. An MS/MS label carries its precursor m/z exactly

**Status:** Implemented (2026-10-19), issue #518. Amends D60's known
limits. Keeps D60's rank rule for isomers.

**Decision.** The stem of an MS/MS column label is `p` followed by the
shortest decimal that reads back as exactly the stored precursor m/z
(Python's `repr`, a trailing `.0` dropped). The merged-precursor warning
spells m/z values the same way. Two precursors with different m/z always
get different stems. Only precursors with the same m/z (isomers kept apart
by mobility) share one, and they take the `_1` suffix in mobility order
(D60).

**Why.** D2 named labels after the m/z so that `anndata.concat` cannot put
two precursors in one column. At six significant digits (`p{mz:g}`),
precursors 1 to 10 ppm apart shared a stem and were ranked as isomers, so
one label could name different precursors in two stores: 798.5401 and
798.5404 were both `p798.54`, and the second took `_1` only in a store that
also held the first.

**What changes in a store.** Only labels that were rounded change:
`p1179.73_mz…` becomes `p1179.731_mz…`. Exact ones such as `p936.578` and
`p500` keep their spelling, and so do all the labels of
tims_msms_pos_brain1. Joining an old and a new store of one run on labels
drops the changed columns; align on `(precursor_mz, precursor_mobility)`,
as the docs already say.

**Alternatives.**

- "Keep six digits and document the limit." No store change, but the
  docs' promise that a label names one precursor stays untrue.
- "Fixed decimals, for example four." Uniform, but nearly every label
  changes.

## D68. A timsTOF source that Bruker's library must write into has to be writable

**Status:** Implemented (2026-10-19), issue #518. Stored data does not
change.

**Decision.** Before the library opens a `.d`, Thyra reads the header bytes
of each database it will open: the analysis database always, and
`calibration.sqlite` when the recalibrated state is on. If one is in WAL
mode and the `.d` folder cannot be written (tested by creating and removing
an empty folder; a warning names it if it cannot be removed), the
conversion stops with one line. The line names the
database and says to copy the `.d` somewhere writable, or, for
`calibration.sqlite` alone, to pass `--no-recalibrated`. The docs say that
the library's `-shm` and `-wal` files appear in the `.d` during a run and
can be left by a hard stop. `thyra metadata` opens no library and is not
checked.

**Why.** The library opens a WAL database read-write in place and puts
`-shm` and `-wal` files beside it while it runs; every real
`calibration.sqlite` checked is in WAL mode, and so is the pea
`analysis.tsf`. In a folder it cannot write, a run ended in a traceback, or
exited 0 with the recalibrated state silently not applied while the store
named it. Checked on a copy of `synthetic_tims.d` with a deny-create
permission on Windows: refused in one line, exit 1, no output and nothing
new in the `.d`; `--no-recalibrated` converted.

**Alternatives.**

- "Document only." A read-only source would still fail with a traceback,
  or convert with a state silently not applied.
- "Open a temporary copy of the sqlite files." The source would stay
  untouched, but it costs a copy per run, possibly of the multi-GB binary
  on Windows. It waits until a measurement on a real read-only share shows
  the refusal is not enough.

## D69. An optical image is stored as its file shows it, and placed only where the store aligns it

**Status:** Implemented (2026-10-19), issue #524. Generalises D59's rule for
the overview photo, and amends D29, D33 and D48 where they place an
unaligned image.

**Decision.** Four rules for every optical image a conversion stores.

- **Placement.** An image is in `"global"` and the dataset's coordinate
  system only when the store holds an alignment for it: the image the
  `.mis` (or an archive's teach points) names, once the TIC-to-image affine
  is built, and a resampled copy of it, scaled onto it. Every other image
  gets one coordinate system named after its element, with `Identity` and
  nothing else. That covers a photo beside a `.d` with no `.mis`, a `.mis`
  image whose alignment failed, every archive image with no teach-points
  affine, and the copies of an alignment image that could not be aligned.
  `alignment_element` still names the image the `.mis` designates.
- **Orientation.** An image whose file carries an orientation tag (EXIF and
  TIFF tag 274) other than 1 is stored upright: turned, and mirrored where
  the tag says so. A quarter turn swaps height and width. The alignment
  image and its resampled copies are the exception: they stay as encoded,
  and the conversion warns. `optical_images.elements.<name>.orientation`
  records `tag` and `applied` on each element whose file carries a tag
  from 2 to 8.
- **Decoding.** A file that starts like a TIFF is read by tifffile,
  whatever its suffix or declared media type. A TIFF's colour model is
  applied: a palette becomes RGB, inverted grey is turned back, CMYK and
  uncompressed YCbCr become RGB. A colour model Thyra cannot turn into grey
  or RGB is skipped with a warning. A TIFF whose channels span several pages
  (an OME-TIFF) is stored with all of them; a z stack or time series keeps
  its first plane, with a warning. A PNG with 16 bits a colour sample is
  decoded at full depth by imagecodecs, with its own channel count.
- **Names.** The part of an element name after `<dataset_id>_optical_`
  keeps only lower-case ASCII letters, digits, `_` and `.`, has no dot at
  either end, and is at most 30 characters. A second file whose name
  cleans to a taken one gets `_1`, `_2` and so on, inside the 30. A second
  file with the same stem (one picture saved as `.tif` and `.jpg`) still
  replaces the first, with a warning, unless the first is the alignment
  image: then it gets a numbered name
  ([D75](#d75-the-photo-crop-follows-the-regions-converted-and-one-affine-places-every-region)).

Three related fixes from the same issue. An archive image that cannot be
copied out (a damaged member, a name the disk refuses) is skipped with a
warning naming the member. A PNG whose header declares more pixel data
than its file can hold, or any whole-page decode larger than the free
memory, is skipped before the store declares it. Skip warnings name the
archive member, not a temporary copy.

**Why.** Identity to `"global"` drew an unaligned photo at 1 micrometre per
photo pixel over the MSI data: the public desi_colad_centroid archive's
1,920 x 1,200 px photo covered 1.92 x 1.2 mm of a 13.2 x 13.6 mm image. That
placement is a guess about the frame, the kind D29 declined for an archive's
own affine. The decoder ignored orientation on the TIFF route, and Pillow
applied it to a TIFF on some compressions only, so a camera photo was stored
sideways and the same bytes gave two stores under two suffixes. Palette
TIFFs were stored as indices, 16-bit colour PNGs as their high bytes, and
an OME-TIFF as its first channel, all at exit 0. Names could collide (one
image replaced another), end in a dot (Windows drops it, so the root attrs
named a missing element), or hold a character spatialdata refuses, which
failed the whole conversion.

**Why the alignment image is not turned.** Nobody knows whether a `.mis`
states its teaching points against the encoded or the upright pixels, and
no tagged FlexImaging file exists to check. Turning the pixels alone could
break the alignment; keeping them as encoded cannot.

**What changes in a store.** Only stores with optical images. An unaligned
image leaves `"global"`, the dataset system and `"acquisition"`; its pixels
do not move. A tagged non-alignment image is stored upright. A palette,
inverted-grey, CMYK or YCbCr TIFF, a multi-page OME-TIFF and a 16-bit colour
PNG change pixels, dtype or channel count. An element whose name held a
character outside the new set, ended in a dot, or clashed with another
changes name, in `optical_images` and `reference_element` too. Aligned grey
and RGB single-page images, the TIC, the pixel polygons and the
`coordinate_systems` attr do not change, and `convention_version` stays 1. No image in the test data, in
tims_msms_pos_brain1 or in the public mzPeak archives carries an orientation
tag.

**Alternatives.**

- "Keep the Identity to `global` and document it." No store change, but
  viewers keep drawing a misplaced overlay.
- "Stretch an unaligned image over the MSI extent." A rough overlay that
  looks right, and the frame guess D29 declined.
- "Turn the alignment image too, and its teaching points with it." Right if
  the points refer to the encoded pixels, wrong if they already refer to the
  upright ones.
- "Store every image as encoded, and warn." Alignment-safe, but sideways
  photos stay sideways.
- "Refuse a name clash." One stray file name would block the conversion.

## D70. An mzPeak archive leaves MSn out per pixel, and no chunk loses a point

**Status:** Implemented (2026-10-19), issue #522. Amends D28 and D29. Amended by
D73: the imzML route applies the same MSn rule, which ends its first known
limit.

**Decision.** Five rules for the mzPeak reader.

- **MSn spectra are left out per pixel.** A spectrum of level 2 and up is
  left out where a spectrum of a lower level sits on the same pixel. Level 0
  and a null level count as MS1 there. An MSn spectrum on a pixel of its own
  is kept. The log counts what was left out and what was kept. D29 decided
  for the whole archive: one MS1 spectrum anywhere left every MSn spectrum
  out, and only level 1 counted.
- **A region's frame count may count those MSn spectra.** A listed region
  matches when its frames equal the spectra in its box after MSn is left
  out, or before.
- **A spectrum that is not converted shapes neither the axis nor the range.**
  The raw axis is built from the spectra kept, always. Before, data rows
  of a spectrum the metadata does not list reached it. The observed mass
  range counts a spectrum only when it has points in the member that is
  read; a count column with no value at all counts every spectrum. A kept
  spectrum whose lowest or highest m/z is null has it taken from its points.
- **An MS-Numpress chunk's ends are matched within its own scale.** An end is
  replaced by its bound when it lies within one over the chunk's scale, and
  never less than 1e-5. Other encodings keep 1e-5.
- **Padding is a point with no m/z and no intensity, or a zero one.** A NaN
  m/z beside an intensity is passed on, and the converter drops or refuses
  it by D55's rule, as for imzML. A null m/z beside an intensity of 0 is
  still padding: MS-Numpress stores a true m/z of 0.0 as its null.

Three kinds of malformed chunk row are refused by name. While reading: an
MS-Numpress row whose buffer holds another number of values than it has
intensities, at any count, or ends inside its header; and a delta row whose
`mz_chunk_start` is null where its first m/z belongs. At open: a transform
column such as `intensity_numpress_slof_bytes` that holds a value beside
`intensity`. A plain row's null start is only a first point with no m/z,
and is read as any such point is.

**Why.** The 4.3.3 stress test (#427) found each of these converting to a
wrong store at exit 0.

| Finding | Before | After |
|---|---|---|
| MS1 on 3 pixels, MS2 on 4 others, two regions | 3 of 7 spectra, one region | 7 spectra, regions 0 and 1 |
| MS1 at level 0 beside MS2 on every pixel | MS2 summed in: rows of 1003 | rows of 3 |
| chunked MS-Numpress, last m/z 51,234.57 | 30 of 730 counts per pixel | 730, as its point twin |
| data rows for a spectrum not listed | empty `var` columns, `--no-resample` | the kept spectra's m/z only |
| peaks-only spectrum in a profile archive | 89,000 bins over 60 to 505 | 1,000 bins over 500 to 505 |
| null `lowest_observed_mz`, first m/z 0.0 | low peaks dropped by the default axis | row sums as the imzML |
| NaN m/z with intensity 50 | dropped as padding, at exit 0 | refused on `--no-resample`, as imzML |

An MS-Numpress value is stored as the whole number nearest to it times the
scale, so it decodes within half of one over the scale. `mzpeak-convert`
picks the scale that fits the chunk's largest m/z into 31 bits. Above m/z
42,950 that half step can pass 1e-5, and the end kept its decoded value. It
then fell outside the declared range, and the default axis dropped it.

The per-pixel rule keeps D29's reason: an MS2 spectrum summed into its
pixel's MS1 spectrum was never measured. A pixel with only MSn spectra has
no MS1 spectrum to spoil. Leaving it out lost the pixel, and then its
region's frame count no longer matched.

**What changes in a store.** Only archives of the shapes in the table. The
MSn rule changes nothing for archives with only MS1 spectra, or with MS1
and MSn on every pixel at level 1. Chunked MS-Numpress archives change only where a chunk's
last or first m/z decodes more than 1e-5 off, above m/z 42,950 or so. A
NaN m/z with an intensity is now refused with `--no-resample` and dropped
with a warning on a resampled axis.

**Alternatives.**

- "Keep the archive-wide rule, and count level 0 as MS1." That still drops
  every MSn-only pixel, which is the defect.
- "Allow every chunk 1e-4." A bound that far off may not be the value of
  that point. One step of the chunk's own scale is the error the encoding
  can make, and no more.
- "Count only the regions' spectra before the filter." No Bruker archive
  with MS1 and MSn on one pixel was at hand to say which count
  `mzpeak-convert` writes, so either is accepted.

**Known limits.**

- The imzML reader reads no MS level yet. An imzML with an MS2 spectrum on
  an MS1 pixel still sums it in, while its archive leaves it out. Reading
  `MS:1000511` there needs the parser setup and a group-stated level, and
  is left to its own change.
- A null observed bound costs one pass over the signal member when the
  essential metadata is read. Only archives with such a bound pay it.

## D71. TIC-preserving resampling keeps the current measured inside the window

**Status:** Implemented (2026-10-19), issue #516. Extends D13 and D22,
supersedes the "share is measured by area" clause of D36 and the tied-m/z
clause of D50.

**Decision.** Seven rules, and two messages.

- **A cropped spectrum keeps the sum of its samples inside the range.** When
  the declared range cuts a spectrum, the preserved total is the sum of its
  samples inside the range, bounds included: the samples `nearest_neighbor`
  keeps. The area share is used only when the range holds no sample at all,
  a window narrower than the spacing of two samples.
- **Samples outside the range are counted and reported once,** as D22 says
  both strategies do. A window that holds no peak is refused as "every peak
  fell outside the target mass range", not as an all-zero source.
- **Samples at one m/z are added up first,** sorted by value so the sum is
  the same float in any order. This is what `nearest_neighbor` does, and
  what a timsTOF `.d` folder does per TOF index.
- **When interpolation leaves no bin,** each sample inside the range goes
  to its nearest bin, as `nearest_neighbor` places it. That happens when
  every sample falls between two axis points, or when the gap tolerance
  vetoes every bin. When the range holds no sample, nothing is placed.
- **A spectrum whose largest value is above 2^900 is resampled scaled down**
  by a power of two, then scaled back. Every other spectrum is resampled
  bit for bit as before. Inside a scaled spectrum, a value more than about
  2^1020 below its largest loses precision or reads as zero.
- **The table's `average_spectrum` divides before it sums in a column whose
  sum overflows.** Every other column is summed and then divided, as before.
- **A range starting at or below m/z 0 is refused on a `reflector_tof`,
  `orbitrap` or `fticr` axis, and on a `tof` axis whose law has `A = 0`,**
  in one line that names where the bound came from. Their bin counts divide
  by the lower end or take its logarithm, and a centroid source with a peak
  at m/z 0 ended the default conversion in a `ZeroDivisionError`.
  `linear_tof` axes, and a `tof` law with `A > 0`, lay a range from 0, and
  still do; below 0 they take a square root, so a negative lower end is
  refused there too. A `constant` axis lays any range.
- **Messages.** A store emptied by `--intensity-threshold` is refused as
  "every peak was below --intensity-threshold", not "the reader yielded no
  spectra at all". The "dropping N empty pixel rows" line counts positions
  the reader gave no spectrum for apart from measured spectra left empty;
  only the first can lie outside the acquisition polygon.

**Why.** The area rule bridged gaps left by `--intensity-threshold`, or by a
profile file that stores only its peaks, with straight lines. The stored
total then followed the spacing of the samples, not the ions. On public
files the in-window current came out 0.26 to 4.9 times what was measured
inside the window, and 126 times on a constructed one, at exit 0;
`--resample-gap-tolerance` changed nothing. A sum is in the same unit as the
total it is a share of, matches `nearest_neighbor` and `docs/resampling.md`,
and the fallback keeps the window between two samples the area rule was
chosen for.

Tied m/z made a store depend on the order of the points. Interpolation read
one of the tied points and rebinning another, so two scan orders of the same
points differed in half the bins, by up to 41 over the full range, and a
narrow window stored 706.7, 660.1 or 631.9 for one input. A pixel whose
samples all fell between two axis points, under a gap tolerance, was dropped
and logged as outside the acquisition polygon. Values near the float64
maximum overflowed the rescale and dropped three of four rows, and made the
average spectrum infinite on every route although every stored value was
finite.

**What changes in a store.**

- `X`, the TIC image and the average spectrum of every `tic_preserving`
  store whose spectra reach past the declared range. That includes the
  default Waters MRT and rapifleX profile routes with a narrowed window.
  Each pixel now holds the current measured inside the range, the total
  `nearest_neighbor` stores.
- `tic_preserving` stores of sources that list one m/z more than once, such
  as an imzML with one point per mobility scan. A spectrum whose only
  repeat is a tie no longer counts as a bin to rebin (D50).
- `tic_preserving` pixels that were dropped because interpolation left
  nothing now get a row.
- Stores with values near the float64 maximum keep every row and get a
  finite average spectrum.
- A source with a peak at m/z 0 on a `reflector_tof`, `orbitrap` or `fticr`
  axis, or a `tof` axis with `A = 0`, is refused where it crashed.

A `tic_preserving` store whose spectra lie inside the range and hold no
repeated m/z does not move, and neither does any `nearest_neighbor` or
`--no-resample` store with finite column sums.

**Known limits.** A window that holds no sample still takes the area share,
which bridges a gap unless a gap tolerance is set. The per-region average
spectrum is still summed before it is divided, so it can overflow on such
values. On the plan route, the refusal of an m/z 0 bound names the
`--resample-min-mz` flags rather than the plan.

**Alternatives.**

- "Keep the area rule but stop it bridging gaps", integrating only across
  spacings within the gap tolerance. The result would depend on a setting
  users rarely set, and the default Waters route would stay wrong.
- "Keep the rule and warn when a cropped run meets gapped spectra." No value
  would change, and stores would stay wrong by up to 126 times at exit 0.
- "Divide every spectrum by its largest value before summing." That moves
  the last bits of every `tic_preserving` store; a power of two moves none.
- "Start a source-derived axis at the smallest positive m/z." The planner
  sees only the source's declared range, so this needs every reader to
  report the smallest positive m/z; the refusal names the way through.

## D72. A Cyclic IMS run is read on its own drift axis, and a declared full scan is MS1

**Status:** Implemented (2026-10-19). Amends D7 and D53. The known limit
on regions that do not share one lattice is lifted by D83 (2026-10-22).

**Decision.** Three rules for Waters runs, each read from the run's own side
files:

- **Guard.** A run whose header says `$$ Instrument: Cyclic-IMS` and whose
  `_extern.inf` names `Acquisition Device: WatersADC` has its drift bins
  read, as a run declaring `Acquisition Algorithm: ADC Mode` does (D53).
  Both fields are required. A Synapt G2-Si and an MRT also name WatersADC.
- **Drift axis.** Where a function's section of `_extern.inf` declares
  `IMS Cycle Start Delay` and `Pushes Per IMS Increment`, bin `i` is stored
  at `(delay + i * pushes) * period` ms. The period is the step of the
  library's `getMobility` axis. Both fields count pusher periods; the
  delay's "(ms)" label is wrong. Functions with different timings, a
  timing that is not two whole numbers of pushes (one field alone, a
  fraction, fewer than one push per bin), or uneven library steps record
  no mobility, with a warning, rather than an axis nobody measured. A
  function without either field keeps the library axis.
- **MS level.** A function whose section declares `MSMSMode` as `MS` or
  `SIM` holds no fragment spectra, whatever precursor its scans report. Any
  other mode, or none, keeps the D7 rule: a reported precursor means MS/MS.

The store also records what the run declares: `format_specific.
function_settings` (`msms_mode`, a SIM function's `sim_set_mass_mz`,
`ccs_calibration`), `format_specific.cyclic_sequence` (sequence name, IMS
cycle time, each step's vendor function code, array mode and time, from the
method file `<run>.xml`), `start_delay_pushes` and `pushes_per_bin` in
`format_specific.ion_mobility`, and `instrument_model` "SELECT SERIES
Cyclic IMS".

**Why the guard.** Every Cyclic run lacked `Acquisition Algorithm`, so all
of them were refused. Measured on 25 Cyclic MALDI acquisitions (cell
images at 5, 10 and 20 um, lipid standards, two HDSIM runs; every
acquisition in the set but MassLynx's own exports and two empty runs), each
read in its own process with the guard lifted: on every checked scan (all
of them, or 3,000 spread evenly over the larger runs) the raw bins summed to
the profile spectrum with the same m/z points and a largest difference of
0.0. Nothing crashed, and no file in the source folders changed. Over the
695 in-house Synapt directories with drift files, the 35 G1 runs carry
neither field (their header names only a serial); the 659 G2-Si runs
declare ADC Mode and no Cyclic label.

**Why the axis.** `getMobility` gives 0 to 13.133 ms in 0.06600 ms steps on
every Cyclic run, one push per bin. MassLynx's own drift-time export of the
same runs (`*_1_dt.raw`, one scan per bin, retention time = drift ms) gives
309.002 to 348.403 ms in 0.19799 ms steps. The function declares delay 4682
and 3 pushes per bin, and `(4682 + 3 i) * 0.0659979` reproduces the export
on all 200 bins to 1.8e-5 ms (float32 precision), on a full-scan and a SIM
run. Synapt runs declare neither field in a function section (only `ADC
Pushes Per IMS Increment 1` under the instrument configuration), so their
axis is unchanged.

**Why the MS level.** Every Cyclic full scan reports `precursorMz` 500.0,
which is `Quad.MSFixedMass`; the HDSIM runs report 524.4, the quadrupole set
mass. Both report MS level 1 and declare `MSMSMode` `MS` and `SIM`. With
the guard lifted, D7's rule stored the full scan as MS2, CID, target 500.
No Synapt or Xevo `_extern.inf` among the in-house runs declares `MSMSMode`
at all, so none of their stores changes.

**No CCS.** Thyra derives CCS from 1/K0 only (`ccs_applies` refuses drift
time). Every Cyclic run declares `CCS Calibration: No`, recorded as such.
After several passes an ion outside the targeted mobility window can lap
the ones inside it, so a multipass arrival time is a drift time only within
that window.

**What changes in a store.**

- Cyclic IMS runs, which were refused, now convert with
  `uns["mobility_axis"]` on the declared timing (309.0 to 348.4 ms, bin
  width 0.198 ms on the runs measured) and `uns["mobility_heatmap"]`. Under
  `--waters-spectrum profile` the heatmap's `current_ratio` is 1.0000 and
  the grid marginal 1.0.
- A Waters function declaring `MSMSMode` `MS` or `SIM` stores
  `ms_analysis.fragmentation` at MS level 1 with no window and no
  dissociation, and no `uns["msms_schedule"]`. An excluded such function
  lists no `precursor_mz`.
- A run declaring none of these fields stores what it stored before. That
  is every in-house Synapt and Xevo run read; no MRT run was re-read.

**D7 and D53, amended.** D7's "`precursor_mz` is the field to trust" holds
where the function declares no `MSMSMode`. D53's "equal steps from exactly
0 ms" and "no cyclic IMS image was tested" describe Synapt runs only.

**Known limits.**

- A multipass separation's pass count is not a field of the method, so it
  is not recorded; the sequence steps are kept as the vendor numbers them.
- A Cyclic run split across functions at the file-size cap (D7) has not
  been seen: the largest run measured keeps a 2.64 GB `_func001.dat` as
  one function. On a Synapt the tail of a split raster has a function
  section of its own, as the public XS image shows. A Cyclic chunk without
  one would record no mobility, since its timing differs from the first
  chunk's, and would fall back to D7's precursor rule, so a full scan
  reporting 500.0 would be stored as MS/MS. Nothing measured has this
  shape.
- A Cyclic run whose regions do not share one stage lattice is refused by
  the raster fit, as any Waters run is. One lipid-standard run with two
  line scans is.

## D73. An imzML leaves MSn out per pixel, as an mzPeak archive does

**Status:** Implemented (2026-10-19), issue #522. Completes D70 for imzML.

**Decision.** The imzML route reads each spectrum's MS level and applies
D70's rule. A spectrum of level 2 and up is left out where a spectrum of
level 1, level 0 or no stated level sits on the same pixel. The pixel is the
spectrum's x, y and z, so a spectrum on another plane is on another pixel.
An MSn spectrum on a pixel of its own is kept. The log counts what was left
out and what was kept, as the mzPeak reader's does. Both routes call one
function for the rule.

- **The level is read on the spectrum, then in its groups.** A spectrum's
  own `MS:1000511` comes first. Otherwise the groups it references are read,
  the last one listed that states a level, as pyimzml resolves them. A value
  that is not a whole number, or is below 0, counts as none. A level above
  255 is kept as 255, and the log names it so; it is MSn all the same.
- **A spectrum left out is not read.** It adds nothing to `X`, the raw mass
  axis, the mass range, `n_spectra` or the peak count. Nor does it add
  anything to the mobility point clouds or the mobility range, and it is
  not asked for a mobility array or for one as long as its m/z array. A
  continuous file takes its shared axis from the first spectrum that is
  kept.
- **A spectrum keeps its place in the list.** `obs["acquisition_order"]`
  still counts every `<spectrum>`, so the spectra after one left out do not
  move (D26).

**Why.** The 4.3.3 stress test (#427, F26) converted an imzML with an MS2
spectrum on an MS1 pixel, and its own `mzpeak-convert` archive.

| Input | Before | After |
|---|---|---|
| the imzML, level stated on each spectrum | TIC 40 and 1280 | 40 and 80 |
| the imzML, level stated in a group | TIC 40 and 1280 | 40 and 80 |
| its mzPeak archive (D70) | TIC 40 and 80 | 40 and 80 |

D70 gives the reason for the rule: an MS2 spectrum summed into its pixel's
MS1 spectrum was never measured. Reading the level through pyimzml's
per-spectrum values does not work. They never follow a group reference, and
the pyimzml writer and `mzpeak-convert` both state the level in a group.
The full per-spectrum metadata would follow it, at about 10 GiB on a
918,855-spectrum file. So the level is read in the step that already records
each spectrum's z (D63), at one byte per spectrum.

**What changes in a store.** Only files with an MSn spectrum on a pixel that
holds MS1 too. A file with no level, with level 0 or 1 only, or with MSn on
pixels of their own converts as before. Before and after were converted
from the same tree otherwise: the committed imzML and timsTOF fixtures, two
example imzML sets, two profile imzML exports, a 918,855-spectrum imzML and
four mzPeak archives. All 17 stores are identical. Reading the levels adds
about 3 s to the 46 s parse of the 918,855-spectrum file.

**Alternatives.**

- "Refuse an imzML that holds MS1 and MSn on one pixel." Its archive
  converts, so the two routes would still disagree.
- "Ask pyimzml for the level with `include_spectra_metadata`." It misses a
  level stated in a group, which is where most writers state it.

**Known limits.**

- The metadata-only route and a mass axis plan read the imzML header, not
  the spectra, when the header is enough. Their spectrum count and observed
  range then still include MSn spectra that a conversion leaves out. A
  plan's axis can then reach past the converted spectra; no point is lost.
- `format_specific.spectrum_count` stays the number of spectra the file
  holds, as on the metadata-only route.

## D74. Source keys in `uns` are stored under path-safe names, and an overflowing grid confirms no pixel size

**Status:** Implemented (2026-10-20), issue #523. Amends D27.

**Decision.** Two rules, both about values a source states that Thyra
cannot store or test as written.

- **Keys.** Every dict key that becomes a node of `uns` is stored under a
  name Zarr can hold on every platform. Each `/`, `\`, `:`, `*`, `?`, `"`,
  `<`, `>`, `|` or control character, and each dot or space at the end of a
  key, becomes `_`. A key that is safe already keeps its name; one that
  would meet it is given a number, `MS_1000529_2`. Keys inside a list that is
  stored as one JSON string are kept as written: any key is legal there.
  The rule holds for every format, since every reader's metadata reaches
  `uns` through one function.
- **Pixel size.** A declared pixel count or extent is used to test a lone
  `IMS:1000046` only when it is a finite number above zero, and so is the
  product of the side and the count. A value that is not a finite number,
  a whole number too large for a float or a JSON `Infinity` for one, is
  ignored with a warning, as a value that is not a number already was.

**Why.** The 4.3.3 stress test (#427) found both.

| Input | Before | After |
|---|---|---|
| mzPeak `run` block with a key `MS:1000529`, on Windows | save fails after the conversion, exit 1, `.failed` left | stored as `MS_1000529` |
| index metadata key `vendor:note` | the same | `vendor_note` |
| `imaging.notes` key `a/b` | save fails: Zarr refuses a `/` in a key | `a_b` |
| old-name `IMS:1000046` 10000, count 1e307, extent 1e306 mm | 10000 um pixel, logged as confirmed | no pixel size: `--pixel-size` asked for |
| the same, count 1e305, extent 1e307 mm | 10000 um pixel | no pixel size |

A key with a colon is legal on Linux, so the same archive converted there,
and its store kept the colon. The extent of 1e306 mm is infinite in
micrometres, and `isclose(inf, inf)` is true whatever the pixel size, so
the test confirmed a reading it could not have made.

**What changes in a store.** Key names, in `uns` only, for a source whose
metadata keys hold one of those characters. The root attributes keep the
keys as the source wrote them; there they are JSON keys, legal everywhere.
The metadata of 14 local datasets was read and walked as the store writes
it: five imzML, four timsTOF (TSF and TDF, one with MS/MS), three Waters
(one Cyclic IMS), a PHI ToF-SIMS run and a rapifleX run. No key was
renamed. No solariX acquisition or published mzPeak archive was at hand.
The pixel size changes only for an mzPeak archive whose declared grid
overflows, or that states a pixel size as a JSON `Infinity` or `NaN`: such
a value was taken as an infinite pixel and is now ignored, so the
conversion asks for `--pixel-size`. No writer is known to produce either.

**Alternatives.**

- "Store the free-form blocks as JSON strings, as lists already are." Every
  block would change its type in every store, to rescue a few keys.
- "Refuse a source whose keys are not path-safe." The keys are provenance;
  losing a conversion over one is the wrong trade.
- "Drop the unsafe characters instead of replacing them." `a:b` and `ab`
  would then meet, and the change would be harder to read back.

**Known limits.**

- Two keys that differ only in case are two nodes, which a file system that
  folds case, as Windows does, puts in one folder. That was so before this
  rule too. No source with such keys was found.
- The renamed key is logged at DEBUG only.

## D75. The photo crop follows the regions converted, and one affine places every region

**Status:** Implemented (2026-10-20), issue #520. Amends
[D31](#d31-fleximaging-spots-are-placed-on-their-lattice-not-on-the-drawn-area)
for the fallback without teaching points, and D69's rule for names.

**Decision.** Six rules for a FlexImaging `.mis` and the optical alignment
it gives.

- **The crop follows the regions converted.** Whether the alignment photo
  is cropped, and where, is decided by the regions the conversion converts,
  counted from the data or `--region`, never by how many Areas the `.mis`
  maps. One converted region is cropped around its cells as the
  TIC-to-image affine places them, plus the 10% margin, whether or not it
  has an Area. Its `crop` then records `region_name` null. A run converted
  whole with several regions keeps the whole photo.
- **A region with no Area is placed, and named.** Its pixels, polygons
  included, are placed by the affine every region shares, and the
  conversion warns once for each converted region without an Area.
- **One affine, with or without teaching points.** Without the D31
  lattice, the raster is stretched over the box of the mapped Areas, with
  TIC cell 0 at the raster's first position. The TIC image, every polygon
  and the crop all use that one matrix. Before, each region's polygons used
  its own Area while the TIC used the box, and a lone region ignored where
  the raster started.
- **An Area that cannot be read keeps its place.** It stays in the list
  with its name and no outline, with a warning, so the Areas after it still
  pair with their own regions. That covers a point that is not two finite
  numbers within 1e9, fewer than two points, no width or height, and a
  polygon whose points lie on one line.
  A `--region` name that two Areas share is refused, naming the
  RegionNumbers that carry it.
- **Teaching points must span the photo.** Three points that repeat or lie
  on one line fix no affine. The lattice is not fitted from them, and
  `msi_metadata.alignment.method` is not claimed. A `.mis` number is read
  as a finite decimal; a Raster or TeachPoint that is not one is left out
  with a warning naming it. When the `.mis` gives no usable Raster, a
  warning says the pitch is BeamScanSize.
- **The alignment image is the exact file.** It is the file the reader
  resolved from the `.mis` it adopted, not any file with the same stem. A
  later file never replaces it: a same-stem file gets a numbered name. When
  the named file is not there by that exact name, the one image of that
  name, else of that stem, stands in, as before (`scan.tiff` for
  `scan.tif`); two that match equally are not guessed between. A Rapiflex
  run reads a `.mis` only from its own folder, as its reader does.

`ImagingArea*IndexPos` values are read as exact integers. An empty or NULL
value counts as absent, and the frames' own bounds stand in, key by key, in
the metadata and in the alignment's first raster position alike. A missing
or empty `MzAcqRange` is a one-line refusal.

**Why.** The 4.3.3 stress test (#427) found sixteen ways a `.mis` or its
`GlobalMetadata` reached a store unchecked, all on the synthetic timsTOF and
rapifleX fixtures. Among them: one unreadable Area moved every later Area
onto the wrong region, so `--region gamma` converted another region's
spectra. With one Area for three regions, the photo was cropped to that one
Area, leaving two regions' pixels off the stored photo. Without teaching
points, a TIC cell and its polygon sat up to 0.63 steps apart, or three
steps when the raster started before its Area. Three diagonal teaching
points spread polygons over 1e17 photo pixels. A same-stem `.png` at twice
the size became the alignment image. Bounds above 2^53 grew a 3x2 raster to
5x2.

**What changes in a store.** Only stores that met one of those cases.

- A multi-region run converted whole whose `.mis` maps one Area keeps the
  whole photo and loses `crop`. A region picked with `--region` that has no
  Area is cropped, with `region_name` null.
- Without teaching points: polygons, the TIC transform,
  `raster_to_global_affine` and the crop window, for runs with several
  mapped regions or whose raster starts before its Area. Regions without an
  Area get polygons.
- Which spectra `--region <name>` converts after an unreadable Area. The
  Area's placeholder, `{"name", "unusable"}` with no outline, is now in
  the table's `uns["raw_metadata"]["mis_metadata"]` (timsTOF) and in a
  rapifleX store's `format_specific_metadata.areas`, where the Area was
  left out before.
- Which file the alignment element holds, when a same-stem twin or a
  parent-folder `.mis` was there.
- `dimensions` and `coordinate_offsets_px` for bounds above 2^53.
- `msi_metadata.alignment.method` for teaching points on one line.
- Stores that were tracebacks before.

A single-region run, a run with teaching points, and a run whose Areas all
read do not change, apart from the cases above.

**Alternatives.**

- "Crop a whole multi-region run to the box of all its regions." Every
  multi-region store would change, and only the two cases above needed to.
- "Keep cropping to the mapped Areas, and warn." The regions without an
  Area would still have no photo under their pixels.
- "Refuse a `.mis` with an unreadable Area." The other Areas still place
  their regions, and a placeholder keeps the pairing right.

**Known limits.**

- The solariX reader looks for its `.mis` beside the `.d` only, while
  optical discovery also looks inside it.

## D76. `convert_msi` refuses every keyword it does not use, and checks its values as the CLI does

**Status:** Implemented (2026-10-20), issue #519. Extends D10, which answered
only `sparse_format`. Leaves D44's refusal of `region` as it is.

**Decision.** `convert_msi` takes three kinds of argument: its named
arguments; the converter keywords its docstring lists
(`write_mobility_table`, `mobility_heatmap`, `mobility_grid`,
`mobility_bins`, `mobility_min`, `mobility_max`, `msms_table`, and `use_csc`
for older callers); and, in `reader_options`, only the documented reader
keys (`intensity_threshold`, `use_recalibrated_state`, `tdf_spectrum`,
`region`, `use_centroid`, `spectrum_type`, `max_mass_axis_length`,
`use_appended_calibration`). Anything else is refused before the source is
opened, in one ERROR line naming each keyword. A reader key given as a
plain keyword is pointed to `reader_options`.

- **A reader key the format does not use is named, then left out.** One
  WARNING per key, worded as the CLI words its flags (#260). Both routes
  read one table, `READER_OPTION_FORMATS` in `thyra/reader_options.py`.
  The CLI now hands `convert_msi` only the options that reach the detected
  format, so a plain imzML run logs no "is ignored" line although the CLI
  always sets `use_recalibrated_state`. `region` stays refused by a reader
  that cannot select one (D44): dropping it would convert every region.
- **`conversion_options` holds only values Thyra computed.** A caller's
  leftover keyword used to be merged in after them, so it could overwrite
  one, record an option nothing applied, or fail the save after both
  passes when JSON could not hold it. `handle_3d` is stored as a plain
  boolean.
- **The CLI's value checks reach the API.** `normalize_resampling_config`
  refuses, for a dict and a `ResamplingConfig` alike: a width, reference
  m/z, m/z bound, gap tolerance or `bins_per_fwhm` that is not a positive
  number; `min_mz` not below `max_mz`; a negative TOF law term or a law of
  two zeros; a bin count that is not a whole number of at least 2. A
  dataclass with a string method or axis type is read as a dict's is.
  `convert_msi` checks the config before it opens the source.
- **Numbers and flags.** A pixel size or z spacing is any finite positive
  real number, numpy scalars included; a bool is refused. `handle_3d` is a
  bool, numpy's included. The pixel size is a float from the start, so
  `obs["spatial_x"]` is float64 for `pixel_size_um=20`.
- **Paths and ids.** `convert_msi` refuses an output path under a regular
  file, or under a folder it cannot write, before it opens the source, with
  the CLI's message: the check is one function both call. A dataset id may
  not start with a dot, may not be a Windows device name (`nul`, `con.2`),
  and may hold at most 216 bytes in UTF-8 (216 plain ASCII characters), so
  its longest element name (an optical image, 39 more) fits the 255 bytes a
  folder name holds.
- **Dead resampling options are named.** Each resampling flag given beside
  `--no-resample` gets one WARNING. A gap tolerance under
  `nearest_neighbor` is named as ignored instead of logged as zeroing bins.

**Why.** The 4.3.3 stress test (#427) found ten ways the API and the CLI
disagreed (F38, F58, F60, F119, F130, F139, F150, F151, F152, F161). A
caller following the rendered API reference passed `spectrum_type` as a
keyword: the option did nothing, and the store recorded it as applied. A
negative width gave a 100-bin axis; a gap tolerance of 0 switched the check
off while the log said it was in force. An id starting with a dot wrote a
store SpatialData's reader opens as empty, since it skips names starting
with a dot. Q6 of the stress report weighed refusing, warning, or only
dropping the record; refusing was chosen because a misspelt option then
cannot give a store the caller did not ask for.

**What changes in a store.** A call that passes only documented arguments
writes the same store, with one exception: the "mass axis resampling" step
no longer records a gap tolerance when the method is `nearest_neighbor`,
where it acted on nothing. Ousia's wizard always sends a tolerance, so its
nearest-neighbour stores lose that key and its log gains the warning. A
call that passed an unknown keyword or reader key, or a value the CLI
refuses, now returns `False` and writes nothing. `obs["spatial_x"]` is
float64 where an int pixel size made it an integer column.
`conversion_options` no longer carries a caller's leftover keywords.

**Alternatives.**

- "Warn once per unused keyword and keep it out of `conversion_options`."
  Breaks nobody, but a misspelt option still gives a store the caller did
  not ask for.
- "Accept in silence and stop recording unused keywords." The store stops
  contradicting itself, but an option that does nothing passes without a
  word.

**Known limits.**

- A converter class built directly still takes leftover keywords in
  `options`; they are no longer written to the store.
- A Windows device name is refused on every platform, so a store written
  anywhere can be copied to Windows. An id ending in a dot is not refused;
  Windows drops the dot from a folder name.

## D77. A conversion that will not fit is refused before it allocates, on memory and on disk

**Status:** Implemented (2026-10-20), issue #528. Extends the mass-axis
guard (#251) and D4's reading of free memory; amends D2.

**Decision.** Five changes, each one line in the log where a traceback
used to be.

- **The count ceiling is part of the axis guard.** One table counts at most
  268,435,456 keys (a 1 GiB `uint32` count array). The mass-axis guard now
  checks that first, so an axis of more bins is refused before it is built,
  however much memory is free. The count array still raises past the
  ceiling as a backstop, now as a refusal that is also a `MemoryError`, so
  the callers that leave a sibling table out still catch it.
- **The MS/MS table is left out past the ceiling.** Precursors x mass bins
  past 268,435,456 keys is asked when the table is planned. The table is
  left out with a `WARNING` naming `--resample-bins` and `--no-msms-table`,
  and the summed table is written, as for the mobility grid. D2 said the
  default never breaks a conversion; this was the case where it did.
- **Tables and regions count.** The projection is `bins x (200 + 36 per
  table past the first + 20 per region when there are two or more)` bytes.
  The regions are known before their arrays are made, and are checked
  then. The 2D route's tables are known only after the pre-scan, which
  finds the planes, and are checked then, before the matrix is allocated.
  Every check compares with the free memory read the first time: before
  a resampled axis is built, or once a raw axis is handed over. By the
  later checks the conversion holds part of what the projection counts.
- **Free memory is what this process can take.** The smaller of the
  machine's free memory and a limit on the process: a Windows job object
  (per process or per job), a Linux cgroup (v1 or v2, every level up to the
  root) or `RLIMIT_AS`. A cgroup's reclaimable page cache (`inactive_file`)
  is not counted as used. The mobility tables' guards (D4) read the same
  function.
- **Free disk is checked after the pre-scan.** The counted bytes of every
  matrix about to be allocated (12 bytes per entry, 16 past 2^31 entries),
  summed tables and the siblings fed from the same passes, are compared
  with the free space of the output's drive, where the scratch folders go.
  More than the free space is refused, naming both sizes. More than half
  of it warns, because the store is written to the same drive while the
  files exist. A sibling table counted on its own pass meets the refusal
  when it allocates.

**Why.** The 4.3.3 stress test (#427) found each case ending in a traceback.

| Finding | Before | After |
|---|---|---|
| F59: 270M bins, 150 GB free | the axis built, then `MemoryError` | refused before the axis |
| F154: 240 precursors x 1.2M bins | conversion stopped, no store | MS/MS table left out, store written |
| F148: six regions at 5M bins | 19.5 B/bin per region not counted | counted at 20 |
| F116: ten planes, 2D route | 511 B/bin against 200 | counted at 200 + 36 x 9 |
| F194: 2 GB job limit, 18M bins | admitted, then `MemoryError` | refused, 2.0 GB free |
| F149: matrix twice the drive | `[Errno 28]` after the pre-scan | refused after the pre-scan |

The per-table figure is the slope of the stress test's measurements: one
plane 200 B/bin, ten planes 511. The per-region figure is two `float64`
arrays per region (the total and the mean, 16 B/bin) plus what building
them costs, measured at 19.5. The thresholds for disk follow the memory
guard's shape: the matrix files are an exact count, not a projection, so
the refusal is at all of the free space; the store is compressed, often
to a fraction of the matrix, so its share is a warning.

**What changes in a store.** Nothing. A conversion that fitted before
converts the same way. A conversion that crashed now ends in a refusal or,
for the MS/MS table, in a store without that table. A source with many
regions or many planes, near the memory limit, may now be refused where it
used to be attempted.

**Alternatives.**

- "Say in the docs that 200 bytes per bin is per table." It leaves the
  guard admitting what it knows will not fit.
- "Refuse the conversion when the MS/MS table is too wide." The grid in the
  same case is left out, and the summed table is the store's main content.
- "Read the cgroup's `memory.current` as used." A source read once fills
  the page cache, and the guard would refuse a conversion that fits.

**Known limits.**

- On Linux the matrix files are sparse, so free space falls only as they
  are written. The check after the pre-scan counts every file the passes
  allocate; a sibling table counted later is compared with what is left
  then, which does not subtract files allocated but not yet written.
- The fused passes' count arrays for the grid and the MS/MS table (each up
  to 1 GiB) are allocated before the passes and are not in the projection.
- The 2D route's plane tables each hold a count array during the pre-scan,
  4 bytes per bin, before the table count is known.

## D78. A folder that holds several runs is refused, not picked

**Status:** Implemented (2026-10-22), issue #525.

**Decision.** A run is a `.d` folder Thyra reads, timsTOF or solariX; other
`.d` folders beside it do not count. A folder that holds one run is read as
that run by every entry point: `convert_msi`, `preview_msi`, `thyra
metadata` and the CLI, whose `--interactive-calibration` display now shows
for it. A folder that holds several runs is refused by `convert_msi`,
`preview_msi`, `thyra metadata` and the timsTOF reader with one line that
names every run and asks for one run's path. The CLI still asks which run to
convert, and ends with the same line and exit 1 when it can read no answer,
where click printed a bare `Aborted!`. Format detection is unchanged: it
still names the first run's format, and only the entry points refuse.

**Why.** The library converted the first run in name order with one INFO
line, so a real folder of three runs gave a store of a run nobody chose, at
exit 0. A library cannot ask, and the `.mis` pick already refuses to guess
(#303).

**What changes in a store.** Nothing. A call on a folder of several runs
now writes nothing. A folder holding one solariX run now converts, where
the solariX reader refused it; its store equals the store of the `.d`
itself. Stores converted from a `.d`, or from a folder holding one timsTOF
run, are unchanged.

**Alternatives.**

- "Keep the first run in sorted order, with a warning." A script still gets
  an unchosen run at exit 0.
- "Say in the docs that the library takes the first run." The pick stays
  silent.

## D79. A store appears under its output name only when it is complete

**Status:** Implemented (2026-10-22), issue #525. Extends the `.failed`
move of #245 and #293.

**Decision.** The converter writes its store under a temporary name in the
output's folder, `<name>.incomplete-<8 hex>`. Every step that touches the
store targets that name: the write, the optical pixel streaming and the
removal of a dropped image, and the metadata consolidation. The store is
renamed to the output name as the last step, and marked complete only then.
A failure Python sees moves the temporary store to `<name>.failed`, as
before; the converter does this itself, so a caller that builds a converter
directly gets it too. A run that is killed leaves only the temporary folder
and its `.thyra_*` scratch, and the same command runs again at once. The
rename is retried for two seconds, as Windows now and then refuses one for
a moment. If it still fails, the run fails, its error names the temporary
folder and says the store in it is complete, and nothing is moved to
`.failed`. The Windows path budget (`prepare_zarr_output_path`) is measured
on the temporary name, which is 20 characters longer. `thyra validate` on a
store with no metadata block names an unfinished write as a possible cause,
and names a `.incomplete-<hex>` folder as the temporary folder of a
conversion that did not finish.

**Why.** An out-of-memory kill, or a native abort in pyarrow or pandas,
gives Python no chance to run the move to `.failed`. The partial store
stayed under the output name, the retry was refused because the output
existed, and `thyra validate` blamed an old Thyra version.

**What changes in a store.** Nothing in its contents. While a conversion
runs, its store sits under the temporary name beside the output, and it
appears under the output name only once it is complete.

**Alternatives.**

- "Keep writing in place and narrow the docs to failures the process
  survives." The retry after a kill stays blocked.
- "Mark the store unfinished and let the next run move a marked store
  aside." The path check gains a rule that moves an existing path.

**Known limits.** A Ctrl-C that lands between the consolidation and the
rename is a failed run: the complete store goes to `.failed`, because the
output name does not hold it yet.

## D80. Every output file follows the store path's rule

**Status:** Implemented (2026-10-22), issue #525.

**Decision.** `--log-file`, `thyra metadata -o`, `thyra export-metaspace -o`
(and its default beside the store) and `thyra-check-ontology --output` go
through one helper, `thyra/utils/output_files.py`, that checks the path
before any work starts. A missing folder is created when the file is
written, as it is for the store. A path that is an existing folder, or that
lies under a regular file, is refused with one line and exit 2. A
`--log-file` at the output store's path, or inside it, is refused too,
since writing it would create the output.

**Why.** Each command opened its file without a check and ended in a
traceback, after the metadata read or the ontology check had already run;
a `--log-file` in a new folder stopped the conversion before it started.
Creating rather than refusing keeps one rule for every output; the cost is
that a typo makes a stray folder, as it already does for the store.

**What changes in a store.** Nothing.

---

## D81. The spectrum type and the instrument follow the spectra that are converted

**Status:** Implemented (2026-10-22), issue #526. Amends D29 (the
instrument, and its chilli known limit); relates to D45.

**Decision.** What a store records about its spectra is read from the
spectra it holds, not from the file's summary or from list order.

- **The spectrum type** is the one every converted spectrum states. On
  imzML a spectrum states it through its own terms and the parameter
  groups it references; a group no spectrum references counts for
  nothing. On mzPeak it is the `spectrum_representation` column of the
  spectra kept: on a pixel, and not an MSn spectrum left out beside MS1
  (D70, D73). An archive read from its peaks member stays centroid, and one
  that fills both members stays profile, as the member decides. When the
  converted spectra state both types, no type is recorded and one warning
  gives the counts; the axis is then chosen as for an unknown type. When
  the spectra state one type and the file-level declaration (`fileContent`,
  `file_description.contents`) names only the other, the spectra win, with
  a warning. The declaration counts only when no spectrum states a type,
  and only if it names exactly one; naming both leaves the type unset. The
  processed-mode guess runs only when nothing is stated anywhere.
  `--spectrum-type` still overrides, and its contradiction warning compares
  with the type the spectra give. One function applies the rule on both
  routes.
- **The instrument** (model, serial number, analyzer, axis family) comes
  from the configuration the converted spectra use: on imzML the one each
  spectrum's first scan names, else the run's default; on mzPeak the scans'
  `instrument_configuration_id`, else `run.default_instrument_id`. A name
  the list does not hold counts as none, so the default decides for it
  (the pyimzml writer names one that does not exist). Each configuration
  is resolved alone. A fact the used configurations state differently is
  left unset, with a warning; a configuration that states a TOF counts
  against another's Orbitrap family.
- **A serial number** written as a JSON number is kept as text.
- **The solvent** (`IMS:1001211`) on the used configuration or one of its
  components is written to `preparation.solvent`. `export-metaspace` writes
  an unstated solvent as `""` with a warning, not `"none"`: a DESI source
  sprays one.
- **The z spacing** is written to `conversion_options` only on a volume of
  more than one plane, as in `coordinate_systems`. Elsewhere nothing
  applied it, and the log already said it is not recorded.

**Why.** The 4.3.3 stress test (#427) found five rules answering one
question: the imzML file description (centroid first), the imzML document
order (any group, used or not), the mzPeak contents order, the mzPeak
column (any one centroid label wins), and both members filled (profile).
An imzML and its own archive could be binned on different axes at exit 0,
and so could a file with an unused group or a left-out MS2 spectrum. The
type picks the axis, so this changed `var` and `X`.

| Finding | Before | After |
|---|---|---|
| F176: unused centroid group, 30 profile spectra | centroid, 219,722 bins (group first) | profile, 120,000 bins |
| F178: `fileContent` lists both, spectra profile | imzML centroid, archive first-listed | profile on both routes, either order |
| F179: 8 profile labels and 1 centroid | centroid, 1,988 bins | unset with a warning, 999 bins |
| F179: centroid MS2 left out | centroid | profile, 999 bins |
| F175: two configurations, spectra use rapifleX | model of the first listed | rapifleX in either order |
| F216: serial as a JSON number | dropped | `'1877407.00228'` |
| F72: DESI file states methanol | `Solvent: "none"` | `Solvent: "methanol"` |
| F62: one plane, `--handle-3d --z-spacing 5` | `z_spacing_um: 5` in `conversion_options` | not written |

The report measured the public chilli imzML: `fileContent` declares
profile while all 4,166 spectra state centroid. It was stored as profile on
396,849 bins and its archive as centroid on 978,495. Following the spectra
gives both the archive's axis, which lifts D29's chilli known limit. Two
options lost. Leaving the type unset whenever both are declared would strip
the type from files whose spectra all agree. Recording profile would bin an
all-centroid file on the profile axis.

**What changes.** The recorded type and axis of files whose spectra
disagree with their summary or with each other; the instrument of files
whose spectra do not use the first-listed configuration; `preparation` on
sources that state a solvent; and the two z keys leave `conversion_options`
on every store that is not a multi-plane volume. Every store with one
configuration and spectra that agree with `fileContent` keeps its type,
axis and instrument.

**Known limits.**

- The metadata-only route (`thyra metadata`, preview, `plan_mass_axis`)
  reads the head and the first spectrum, so it follows the first
  spectrum's type and configuration. It differs from a conversion only for
  a file whose later spectra say otherwise, as for polarity (D35).
- The chilli imzML and the other public example files were not converted
  for this entry; the before and after above come from the stress report
  and constructed files.
- Laser parameters (`imzmldict`) are still taken from the first
  configuration that states them, and the polarity rule still reads every
  parameter group, used or not.
- No real file with two instrument configurations has been seen. The rule
  is tested on constructed imzML files and archives.
- A pyimzml parser a caller builds without Thyra's reader records nothing
  per spectrum. `ImzMLMetadataExtractor` over it falls back to `fileContent`
  and every configuration, as it does for polarity and MS levels, so it can
  report another type or instrument than a conversion of the same file.

## D82. A mobility table holds only mobility measured inside its edges, on a shared axis proven on every spectrum

**Status:** Implemented (2026-10-22). Amends D34, D46 and D52.

**Decision.** Five rules about what a mobility-resolved table may hold.

- **A stated mobility range drops what lies outside it.** With
  `--mobility-min` or `--mobility-max`, a point beyond the stated range is
  left out of the grid table and counted in one `WARNING`. A range that
  holds no point writes no grid table, and says so. A range read off the
  values (the default) still clips its overhang into the edge channels, as
  the heatmap does. This is D13's rule for a stated m/z range, applied to
  mobility.
- **The shared-axis table keeps the summed table's m/z range.** A feature
  whose m/z lies outside the range the summed table keeps (D52) is left out
  of `var` and of every row, and the count is said at `WARNING`.
- **A shared (m/z, mobility) feature list is proven on every spectrum.** An
  imzML gets one only when every converted spectrum points at one m/z
  block, measured from the offsets, and every mobility array equals the
  first, NaN for NaN. A copy at the first one's offset is not re-read, nor
  is any block twice. Otherwise the file is per-pixel mobility and is
  converted as the same bytes declared processed are: no shared-axis table,
  no axis `values`, and the heatmap over the values (D46);
  `--mobility-grid` is the way to a mobility table. The same
  measured block test decides the extractor's mass range and the reader's
  common mass axis, so a file that declares continuous mode but writes a
  block per spectrum is read as processed, with one `WARNING`. The file's
  `file_mode` metadata still records what it declares.
- **A non-finite mobility value withholds the mobility.** A NaN or infinite
  value anywhere in an imzML's converted mobility arrays writes the summed
  table and no mobility table, heatmap, grid or `uns["mobility_axis"]`. One
  `WARNING` names the file and the count. `msi_metadata` still records that
  the file declares a mobility array (`ion_mobility.present`), which it does.
- **A TDF mzPeak archive's TOF bin is grouped within 1e-7, relative.** Its
  copies of one bin are summed into one point before the intensity
  threshold, as the `.d` reader sums them.

When a grid is asked for and not written, the reason is said once: no
mobility dimension (a Bruker TSF file, a non-finite value), or
`--no-mobility-table`, naming any sizing option given with it.
`--tdf-spectrum` and the grid options are named as ignored on a TSF file,
on the command line and through `reader_options`.

**Why.** The 4.3.3 stress test (#427) found each of these at exit 0.

| Finding | Before | After |
|---|---|---|
| F13: `synthetic_tims`, `--mobility-min 1.0 --mobility-max 1.5` | 206 points clipped into edge channels, `current_ratio` 1.0 | dropped and counted, `current_ratio` 0.6417 |
| F13: the same at 3.0 to 4.0 | all current in one channel | no grid table, with the reason |
| F129: `mobility_continuous`, m/z 400 to 500 | m/z 300 and 600.25 in the edge columns | left out; the `mz_index` marginal equals the summed table |
| F121: continuous declared, a block per spectrum | m/z 50 and 2000 dropped | row sums 6, 15, 24 on an axis of 50 to 2000 |
| F185: 144 spectra, spectrum 71 differs | stored at spectrum 0's mobility | per-pixel mobility |
| F186: the same, spectrum 9's own m/z | filed under spectrum 0's m/z, or refused with a threshold | per-pixel mobility, either way |
| F168: an infinity in the shared axis | a store `thyra validate` rejects | summed table only; validates |
| F205: equal copies holding a NaN | "differ between spectra" | the NaN named and counted |

A clipped point lies at a mobility it was never measured at, and the
marginal over channels hid it: `current_ratio` stayed 1.0. The drop makes
the ratio say how much the stated range left out.

The 64-spectrum sample was chosen for cost, but nothing recorded the choice,
and the property it tests is about every spectrum. The proof reads each
mobility block once. Measured on a continuous export of 2,000 spectra by
50,000 features with a copy per spectrum (0.80 GB of mobility, 1.60 GB
`.ibd`, warm page cache): 0.57 s for the proof against 0.01 s for the
sample. One pass over the intensities of the same file took 0.22 s, and a
conversion reads them at least twice. The file was synthetic: no real TIMS
export with a mobility array was at hand.

MS-Numpress stores each scan's row on its own scale, so one TOF bin's copies
decode up to the row's largest m/z over 2^31 apart: 7e-8 relative at m/z 20
in a row reaching 3000. TOF bins lie some 1e-6 apart or more (6.06e-6 the
closest pair on `synthetic_tims`). At 1e-12 each scan's copy was tested
alone, and a peak whose sum passes the threshold was lost.

**What changes in a store.**

- Grid tables converted with `--mobility-min` or `--mobility-max` from a
  source with points beyond them: fewer entries in `X`, no features in the
  edge channels from those points, and `mobility_marginal.current_ratio`
  below 1.
- Shared-axis tables under a resampled range narrower than the features:
  fewer features in `var` and fewer columns in `X`.
- An imzML that declares continuous mode but writes a block per spectrum:
  every peak kept, on a processed axis; `essential_metadata.mass_range`
  covers every spectrum; with a mobility array, no shared-axis table.
- A continuous imzML whose mobility copies differ only outside the old
  sample: per-pixel mobility (no shared-axis table).
- An imzML with a NaN or infinite mobility: no mobility table, heatmap or
  `uns["mobility_axis"]`. Before, it got a heatmap and grid without the
  points that had one (D46), or a store `thyra validate` rejected.
- Thresholded TDF mzPeak archives under MS-Numpress keep the peaks their
  `.d` keeps. Other archives already summed their bins and do not change.

Stores converted without stated mobility bounds, from finite imzML files
whose spectra share one block, and from `.d` folders, do not change.

**Objections considered.**

- *Clip on a stated range too, and warn.* The table would still hold
  current at a mobility it was not measured at, and a consumer reads the
  table, not the log.
- *Refuse a non-finite mobility.* The summed table is right and is most of
  the store. The m/z rule (D55) drops a non-finite m/z and keeps the rest;
  a mobility value has no table of its own to drop it from, so the mobility
  products go instead.
- *Keep the sample and warn that it is one.* A wrong axis at exit 0 is the
  failure, and the cost of the proof is one read of the mobility arrays.

**Known limits.**

- The `.d` reader's mobility comes from scan numbers and the Waters
  reader's from drift bins; neither can be NaN, and neither is checked.
- A processed-declared file whose spectra all point at one m/z block is
  still read as processed: the measured test is applied only to a
  continuous declaration.
- The metadata-only route (`preview_msi`) still takes a continuous
  declaration's range from the first spectra when their recorded extrema
  agree (issue #371).
- No MS-Numpress archive from `mzpeak-convert` 0.16.0 was converted for
  this entry; the tolerance rests on the encoding's own error bound.

## D83. A Waters run keeps the regions its method declares, on one grid as on timsTOF

**Status:** Implemented (2026-10-22). Supersedes in part D44 and D57, and
D72's known limit on regions. Amends the one-lattice rule of the Waters
raster fit.

**Decision.** A Waters run whose method file (`<run>.xml` inside the
`.raw`) holds a `MultiRegions` table has regions, and the store keeps them
the way a timsTOF store does.

- **Numbers.** Regions are numbered 1, 2, ... in table order. The table
  states no number. Starting at 1 keeps region 1 on every Waters store
  without a table, and on a run whose scans all lie in the first box.
- **Membership.** A scan belongs to the declared box (start and length,
  in mm) that holds it, widened by half the declared step on each side.
  The widening is needed: a line scan's box is 0.01 mm tall and its step
  0.02 mm. If any positioned scan lies in no box or in two, the regions are
  not used, a warning says why, and the run is read as one raster, as
  before; `region=` is then refused. A scan of a later function chunk that
  lies in no box is off the raster; one in a declared box counts as on it,
  so a region acquired only in the tail chunk of a split raster is kept.
- **Placement.** Each region holding scans gets its own lattice, so the
  off-raster refusal is tested within one region. All regions sit on one
  grid at their stage positions. The pitch is the regions' shared fitted
  pitch, which must agree within 1 %; on an axis no region can measure,
  it is the declared step. The grid starts at the smallest region origin
  on each axis. Each region's origin is rounded to the nearest grid line,
  and `format_specific.regions` records its stage origin, grid offset and
  `placement_shift_um`: how far its stage origin lies from the grid line
  it was placed on.
- **Summary.** `obs["region_number"]`, `uns["regions"]` (number, scan count,
  box), `uns["average_spectrum_per_region"]` and `thyra.list_regions` follow
  D57. A declared region with no scan is not listed.
- **Selection.** `region=N` (`--region N`) converts one declared region
  alone, and the store describes only that region (D57). A run without a
  table accepts `region=1`, the whole run, as `list_regions` reports it. A
  run whose declared regions could not be used refuses every `region=`.
  Readers now state whether they select regions (`selects_regions`); the
  base reader refuses `region` for every reader that does not.
- **Line scans.** An axis with one stage position takes its pixel size
  from the method's `MaldiXStep` or `MaldiYStep`, where declared, instead
  of asking for `--pixel-size`. `format_specific.pixel_size_source` says,
  per axis, whether the size came from the stage positions or the method
  step. `--pixel-size` still overrides. A single-spot run still has no
  pixel size and is still refused as no image.

`thyra/core/region_map.py` gains `summarize_regions`, which builds the
region entries from `(x, y, region)` per scan. The rapifleX reader uses it
too; its entries are unchanged.

**Why.** A lipid-standard run with two 25-pixel line scans, one on each
standard's spot, was refused: the two lines sit 10.16 mm apart in y, and
the second starts 170 um (8.46 pitches) to the right of the first, so no
one lattice holds both. The method file states the two boxes and the
steps, so placing the regions needs no guess beyond the rounding of one
origin, and that is recorded. The same file states the y step of every
line scan, which Thyra asked the user to type.

timsTOF's region selection was kept as is; only its gate moved to the
reader. The SQL paths of timsTOF and solariX stay, because they keep the
metadata-only preview cheap. The helper serves the readers that build
regions from a list of scans.

**What changes in a store.**

- The two-line run now converts: 33 x 509 pixels at 20.1 x 20.0 um, 50 of
  them with spectra. Region 1 is row 508, columns 0-24; region 2 is row 0,
  columns 8-32, and its stage origin lies +9.2 um (46 % of a pitch) in x
  from column 8.
  `--region 2` gives a 25 x 1 store of region 2 alone.
- A Waters line scan with a declared step now converts without
  `--pixel-size`; its y pixel size is the method's step (20.0 um on the
  lipid standards, whose x pitch is 20.1 um).
- Every Waters store gains `format_specific.pixel_size_source`, and a run
  with a method file gains `declared_step_x_um` and `declared_step_y_um`
  in `raw_metadata.imaging_grid`.
- A run whose table declares two regions with scans in only one (an HDSIM
  run) keeps `region_number` 1; its `uns["regions"]` entry gains a box, and
  `format_specific.regions` holds its placement.
- Measured against main on a 25 x 25 run without a table and on the
  HDSIM run: every array is identical, and only the keys above differ.

**Alternatives.**

- "Refuse a run whose regions do not share a lattice, and only accept
  `--region`." Two conversions for one run, and no store holding both
  standards.
- "Pack the regions next to each other." The grid stays small, but one
  affine no longer places every region (D75), and the stage geometry is
  lost.
- "Number the regions from 0, as Bruker does." Bruker numbers are the
  source's; Waters states none, and 0 would change region 1 on every
  Waters store with one region.

**Known limits.**

- Regions far apart give a large, mostly empty grid: two 1.5 mm regions
  7 mm and 43 mm apart would span about 430 x 2230 pixels. Only filled
  pixels hold spectra.
- A rounded region can sit up to half a pitch from its true place relative
  to the others. Within a region nothing moves.
- Only Cyclic IMS runs carry the method file. No Synapt or Xevo run has
  been seen with one.

## D86. A position, index or vendor header value is checked before it is used, and `region_number` is int64

**Status:** Implemented (2026-10-23), issue #527. Supersedes D57's note on
the type of `region_number` (D26 is not rewritten). Extends D77 with a
projection of the pixel grid.

**Decision.** Every value a reader casts to place a spectrum, size a grid
or build a mass axis is checked first. A value that cannot be used is
refused in one line, or skipped and counted where the reader already
skipped bad records. Nothing is cut, wrapped or invented.

- **Positions, indices and region numbers are whole numbers within
  int64.** One helper, `thyra/core/positions.py`, reads them. A float
  column, as mzPeak allows for positions, takes a value within 1e-6 of a
  whole number as that number: `0.3 / 0.1` is 3, not 2. A sqlite column
  declared INTEGER (timsTOF, solariX) takes no fraction at all: a REAL of
  `100.9999999` there is damage, not writer noise. Text, NULL, NaN,
  infinity and a `uint64` of 2**63 or more are never a position.
- **Refuse or skip, by reader.** mzPeak and solariX refuse the file and
  name the column, the count and the first bad value. timsTOF already
  skipped a frame whose position was text or NULL, counted in the "frames
  with no usable coordinates" total; a frame whose position or region
  number is not whole now joins it. Its region list, region map, region
  boxes, alignment positions and metadata bounds leave the same frames
  out.
- **`obs["region_number"]` is int64** on every store. A source may number
  its regions past 2**31, which the int32 column ended in an
  `OverflowError`.
- **A pitch must keep the far edge finite.** On each axis,
  `(|origin| + pixels + 0.5) x pitch` must be a finite number, checked
  once the grid and every pitch are settled. A finite `--pixel-size` of
  1e308 um wrote `inf` into `obs` and the pixel boxes at exit 0.
- **The pixel grid is projected before it is allocated.** The arrays held
  per grid position cost 11 bytes, plus 8 when the reader knows the
  acquisition order. Past 50 % of free memory the conversion is refused,
  past 25 % it warns: D77's fractions, with the free memory read once for
  both checks. The message names the grid against the spectrum count. The
  imzML reader reads positions as int64, rebases them, and narrows to
  int32 only when they fit, so a position of 2**31 reaches this check.
- **Waters function files.** Each `_FUNCnnn.DAT` needs its `.IDX` and
  `.STS`, and a non-empty `.IDX`, checked before the library opens the run.
  The scan count the library returns must equal the `.IDX` size divided
  by its record size: 22 bytes on the Synapt G1 runs measured, 30 on Xevo
  and Cyclic. An `.IDX` size neither record size divides is warned about,
  not refused: it may come from a MassLynx version not measured. The grid
  re-fit warning says which functions were added or left out.
- **rapifleX.** `Raster` must be exactly two numbers above zero, else it is
  treated as missing, as the `.mis` parser treats it. `Mass Start` and
  `Mass End` are read as finite numbers, so `5e2` converts; `nan` or
  `100,5` is refused naming the file. The `.dat` must hold its whole offset
  table and every listed spectrum; a short one is refused naming the spot.
- **timsTOF pixel size.** A `BeamScanSize` that is not two numbers above
  zero is treated as missing, with a warning. The Raster/BeamScanSize
  mismatch warning says which is larger; "oversampled" only when the beam
  scan is.
- **PHI.** A header number that is not finite falls back to its default; a
  required one is then refused. `HeaderSize` must be positive. An appended
  `Mass/Time` must be finite and above zero. A calibration whose m/z axis
  does not rise strictly with flight time is refused, whatever its source.
  Records that fail the validity test, events off the raster and events
  off the mass axis are counted apart. Events off the raster get a
  `WARNING`, and each tile's events are tested against that tile's own
  raster.

**Why.** Each case was a wrong store at exit 0 or a traceback, which
`docs/cli.md` calls an exception nobody planned for. A cut position put a
spectrum on its neighbour's pixel; a wrapped one dropped spectra or stored
a negative `acquisition_order`; a cut-short Waters or rapifleX file lost
pixels without a word. The real MassLynx library was probed on damaged
copies of a 9-scan Cyclic export: a `.DAT` cut to 40 % returns 3 scans and
still reads as imaging; an empty `.DAT`, an empty or missing `.IDX` and a
missing `.STS` all read as "not an imaging file". The scan-count check
catches the first; the file checks and the count, asked before the
"not imaging" refusal, name the others.

The timsTOF frame lookup reads positions at about ten SQL sites. One rule
in SQL (`col = CAST(col AS INTEGER)`) keeps every one of them on the same
frames as the frame loop, rather than a check at one site and `int()` at
the rest.

**What changes in a store.**

- `obs["region_number"]` is int64 on every store (it was int32 in 2D).
- A source with a bad position, index, header value or a cut-short file
  is refused where it used to convert or crash. A timsTOF frame with a
  position that is not whole is left out and counted.
- A multi-tile PHI store whose events lie past their tile's raster loses
  those events; they used to land on the next tile's pixels.
- Measured against main on 19 sources: the 10 default fixtures, two
  timsTOF runs, a 512-pixel PHI run, two Waters Cyclic runs, a Synapt G1
  run (22-byte index records), a Xevo DESI run (30-byte), a rapifleX run
  and a continuous imzML. Every array and attribute is identical except
  the `region_number` dtype, on every table that has the column (mobility
  and MS/MS tables included).

**Alternatives.**

- "Round every float position, also from sqlite." A REAL with a fraction
  in an INTEGER column is not writer noise; rounding it would place a
  damaged frame on a guessed pixel.
- "Refuse the timsTOF file for one bad frame." The reader already skipped
  frames without coordinates and counted them; one rule for text and for
  fractions keeps that.
- "Keep `region_number` int32 and refuse a region number past 2**31." The
  number is the source's (D57); refusing it would refuse a valid file for
  a storage width.

**Known limits.**

- A Waters `.IDX` cut exactly on a record boundary matches its own scan
  count, so a run cut that way is not caught; one cut elsewhere is only
  warned about. The `.DAT` layout is not read.
- `thyra validate` checks `var`, not `obs`, so a store written by an older
  Thyra with an infinite `spatial_x` still passes it.
