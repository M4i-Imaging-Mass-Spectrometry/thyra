# Migrating to 4.0.0

This page covers upgrading from **3.23.3 to 4.0.0**. For a pin older than
3.23.3, read [If you are coming from further back](#if-you-are-coming-from-further-back-than-3233)
at the bottom first: it is the part most likely to catch you out.

4.0.0 is a major release because seven changes break a published API. **None of
them changes the tables a conversion writes.** Convert the same file with
3.23.3 and with 4.0.0 and the matrices, the `var` and the `obs` come out the
same; two keys leave the metadata, one failure mode that used to be a log line
is now a failure, and one that used to be a numpy error now names the file.
Everything else in the release is a fix, a docs change or a performance change.

If you only call `thyra convert` from the command line, almost nothing here
applies: no flag was removed or renamed in this release. Item 7 is the one you
may notice, and only on a file that never converted — what it prints is now a
sentence rather than a traceback.

## What changed, at a glance

| # | Change | You are affected if you |
|---|--------|-------------------------|
| 1 | A sibling table that cannot be built now fails the conversion | convert sources that write a `_mobility` or `_msms` table |
| 2 | One convention for declining a capability | subclass `BaseMSIReader`, or probe readers with `getattr` |
| 3 | `batch_size` is gone | pass `batch_size` to any reader |
| 4 | `EssentialMetadata` loses two fields, and two keys leave the store | implement the reader ABC, or read `estimated_memory_gb` from a store |
| 5 | `convert()`'s stages take a typed state | subclass `BaseMSIConverter` or call its stage methods |
| 6 | A preview's mass range can be unknown | read `MsiPreview.mz_range`, or build an `EssentialMetadata` yourself |
| 7 | An imzML that declares no spectra is refused by name | catch `IndexError` around an imzML read |

Only item 4 changes what is written to disk.

---

## 1. A sibling table that cannot be built now fails the conversion

**Affects:** anyone converting a source that writes a `_mobility` or `_msms`
sibling table, and anyone treating `convert()`'s return value as the only
success signal.

The two sibling-table builders were wrapped in `except Exception`, so a
failure inside one became an ERROR line and `convert()` still returned True.
That was not "a store without the sibling, still complete", because the summed
table names its siblings in `uns` *before* they are built — at
`mobility_axis["resolved_table"]` and
`msi_metadata.ms_analysis.ion_mobility.resolved_table`, and at
`msms_schedule["resolved_table"]` and `...fragmentation.resolved_table`. The
store was written pointing at an element nobody wrote, and nothing downstream
checks that the pointer resolves — not the schema, not `thyra validate`.

A failure in a sibling builder now propagates and the conversion fails.
Nothing is on disk at that point, so the cost is a conversion you have to
re-run rather than a store that lies to you.

**What to do:** if you were relying on a conversion completing despite a
sibling failing, you will now see the traceback that was being swallowed.
Treat it as a bug report — that is what it always was.

**Not the same as a builder declining.** A builder is still allowed to say
"no table here" by returning `None` (an empty feature listing, a var count
over the ceiling). That is a decision, not a failure, the conversion still
succeeds, and as of this release the name is taken back out of every table
that carried it, so a declined sibling leaves nothing pointing at it.

## 2. One convention for declining a capability

**Affects:** anyone who has written a `BaseMSIReader` subclass, or who probes
readers with `getattr`/`hasattr`.

`BaseMSIReader` answered "unsupported" three different ways. There is now one
rule for every optional capability:

- a `has_*` predicate gates the capability;
- when it is False, a **getter** describing the capability returns `None`;
- when it is False, an **iterator** producing data raises `NotImplementedError`.

Two consequences for callers:

- `has_fragmentation` is new, and derives from `get_fragmentation()` by
  default — a reader that implements the getter and forgets the predicate
  keeps the capability rather than silently losing it.
- `has_fragmentation` and `has_precursor_spectra` are **two capabilities, not
  synonyms**. Describing the fragmentation and being able to separate a
  pixel's precursors are different things: Waters reports a schedule it cannot
  demultiplex, and an MS1 TDF reports a schedule with no windows. Both are
  `has_fragmentation` True and `has_precursor_spectra` False.

**What to do:** custom readers should follow the three-part rule. Callers can
drop defensive probing — `getattr(reader, "has_ion_mobility", False)` and
`callable(getattr(reader, "iter_precursor_spectra", None))` were guarding
against attributes the ABC already guarantees, and direct attribute access is
correct now.

## 3. `batch_size` is gone

**Affects:** anyone passing `batch_size` to a reader.

It was declared by sixteen `iter_*` methods across nine modules and referenced
by two, neither of which used it to change anything: mzPeak logged that it was
ignoring it, and imzML chose between two code paths that walked the same flat
range and called the same per-spectrum function. Verified on a real fixture
before the branch was deleted — `single(batch_size=1)`, `batch(1000)` and the
default all produced the same result.

**What to do:** delete the argument. There is no replacement because it never
selected anything.

- On the `iter_*` methods this is a plain removal. None of them takes
  `**kwargs`, so `reader.iter_spectra(batch_size=10)` already raises a
  `TypeError` naming the argument — Python's own error is loud and specific
  and no wrapper improves on it.
- `ImzMLReader.__init__` is refused explicitly instead, because it forwards
  `**kwargs` and would otherwise accept `batch_size=50` silently. There is
  also a type guard on the positional form: `ImzMLReader(path, 50)` would have
  slid into `cache_coordinates=50` — truthy, and wrong.

## 4. `EssentialMetadata` loses two fields, and two keys leave the store

**Affects:** readers and extractors implementing the ABC, and any consumer
reading `estimated_memory_gb` out of a store.

- `peak_counts_per_pixel` and `get_peak_counts_per_pixel` built the CSR
  `indptr` in a single pass. CSR was retired in 3.22 (design decision D10) and
  nothing has read either since — the method was declared on the reader ABC,
  implemented by six readers and computed by four extractors, all of it
  feeding nothing.
- `estimated_memory_gb` fed the size threshold that chose between the
  in-memory and streaming converters. D11 deleted that router. The number went
  on being computed, logged and written into every store since — twice, under
  `essential_metadata` and again under `processing_stats` — while meaning
  nothing, and not even meaning the same thing per extractor: Bruker and imzML
  estimated a dense footprint, the other four `total_peaks * 16` bytes.

**This one changes what is stored.** `essential_metadata.estimated_memory_gb`
and `processing_stats.estimated_memory_gb` are no longer written. A consumer
reading either must tolerate its absence. Nothing was reading it for a reason
— it was not a meaningful number — but code that indexes it directly will
raise.

**What to do:** two real uses were kept, each moved to where it belongs:

- PHI derives `n_spectra` and `total_peaks` from per-pixel counts, because a
  PHI file is a stream of ion events. That is now
  `PhiReader.occupied_channel_counts()` — PHI's own API rather than a contract
  every reader has to answer.
- Waters counts occupied pixels rather than scans. Only the zero/non-zero test
  survived, so the int32 accumulator became a boolean occupancy bitmap.

## 5. `convert()`'s stages take a typed state

**Affects:** anyone subclassing `BaseMSIConverter` or calling its stage methods
directly.

The four stages of `BaseMSIConverter.convert()` handed each other an untyped
dictionary. Nothing declared what a stage required or produced, the twelve keys
in use could only be learned by grepping for bracket literals, and no checker
could verify any of it in a package shipping the `Typing :: Typed` classifier.

`ConversionState` is that dictionary's keys, typed, threaded through the same
four calls. The parameter is renamed `data_structures` -> `state`.

**There is no behaviour change here** — it is a signature and an attribute
access change. Two notes:

- Optional fields are now `None` rather than absent. The `region_*`
  accumulators exist only for a dataset with a region map, and the dict said so
  by not having the key, so some call sites used `state["region_row_count"]`
  and others `state.get("region_total_intensity")`. They are `None` now, once,
  for everyone.
- `mode` (`"3d_volume"` / `"2d_slices"`) is dropped. It restated
  `self.handle_3d`, which every stage already has — read that instead.

## 6. A preview's mass range can be unknown

**Affects:** anyone reading `MsiPreview.mz_range`, and anyone building an
`EssentialMetadata` by hand.

`preview_msi` no longer parses an imzML's spectrum list or decodes its `.ibd`;
it answers from the head of the document instead (issue #360, with the
measurements in `docs/supported-formats.md`). The one thing a head cannot
always state is the mass range: an imzML written without the observed-m/z
cvParams (`MS:1000528` / `MS:1000527`) records its extrema nowhere but in the
binary, and IONTOF SurfaceLab writes none.

So `MsiPreview.mz_range` is `Optional[Tuple[float, float]]`, and `None` means
"not read" rather than "empty". It is deliberately not filled in from the first
spectrum: that was measured 28 Da narrow on a real file — a range narrow enough
to be believed and wrong enough to mislead.

`EssentialMetadata` gains `mass_range_known`, which is the distinction
`n_spectra_counted` already draws for a count (issue #240). When it is False,
`mass_range` is `(0.0, 0.0)` meaning "not read", which is not a range any
acquisition has.

**This one does not change what is stored**, despite adding a field to
`EssentialMetadata`. `uns["essential_metadata"]` is written from an explicit
five-key dict rather than from the dataclass, so a new field does not reach a
store.

**What to do:** handle `None` where you format `mz_range`, as you already have
to for `n_pixels`. A custom extractor needs no change — `mass_range_known`
defaults to True, which is what every path that reads a range for real should
report.

## 7. An imzML that declares no spectra is refused by name

**Affects:** anyone catching `IndexError` around an imzML read.

An imzML whose `<spectrumList>` holds no `<spectrum>` elements used to fail
with `IndexError: too many indices for array: array is 1-dimensional, but 2
were indexed`, logged as "Failed to initialize ImzML parser". That came out of
pyimzml's constructor, which ends by taking the z extent off the coordinates it
collected, and an empty coordinate list makes that a 1-D array. It named
neither the file nor the reason, and a preview reported it as the whole account
of the file.

Such a file is now refused before the parser is built, with a `ConversionRefused`
naming the file, the term and what is missing.

**What to do:** nothing, unless you were catching `IndexError` specifically.
`ConversionRefused` subclasses `ValueError`, so `except ConversionRefused`,
`except ValueError` and `except Exception` all keep working and gain the
message. `IndexError` was never an outcome Thyra documented or intended.

---

## If you are coming from further back than 3.23.3

The `!` markers on this release are not the whole migration story for an older
pin. At least one hard removal shipped earlier as a **minor**: 3.22.0 dropped
`--sparse-format` and the `sparse_format` keyword (design decision D10), and
passing the keyword now raises. CSC is the only layout written.

If you are upgrading across several minors, read `docs/design-decisions.md`
rather than the breaking-change markers alone.
