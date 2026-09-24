# Metadata Schema

Every store Thyra writes carries a versioned, ontology-mapped metadata
block: `table.uns["msi_metadata"]`. Its base fields mirror the
[METASPACE](https://metaspace2020.org) submission form, so the metadata a
converted dataset carries is, by construction, what a METASPACE submission
needs -- filling it costs nothing extra.

The schema exists because MSI has had no structured metadata convention
the way other spatial omics modalities do. Vendor files spell the same
fact many ways, imzML metadata is free-form cvParams, and every pipeline
re-invents its own dictionary. This block fixes the names, maps them to
ontologies, and ships a validator.

!!! info "Where this sits in the format stack"
    Raw/archival exchange is the domain of
    [mzPeak](https://github.com/HUPO-PSI/mzPeak) (the HUPO-PSI working
    draft succeeding mzML); SpatialData is the analysis-ready layer and
    Thyra's output. Thyra is the bridge between the two, and this
    schema is the analysis-layer metadata contract: PSI CV-anchored
    like the raw layer, so nothing is lost crossing the bridge. Raw
    ragged spectra stay in the raw layer; they are deliberately not
    part of Thyra's output.

```python
import spatialdata as sd

sdata = sd.read_zarr("output.zarr")
block = sdata.tables["msi_dataset_z0"].uns["msi_metadata"]

print(block["schema_version"])                       # "0.8.0"
print(block["ms_analysis"]["pixel_size_um"])         # {"x": 20.0, "y": 20.0}
print(block["ms_analysis"]["ionisation_source"])     # "MALDI"
print(block["ms_analysis"]["ionisation_source_term"])
# {"accession": "MS:1000075",
#  "name": "matrix-assisted laser desorption ionization"}
```

---

## The document

Six sections and a processing list. `ms_analysis` and `provenance` are
written by the converter for every store; `acquisition` and `calibration`
are each written when the reader reports at least one of their facts and are
absent otherwise; `sample` and `preparation` describe things no raw file
records (what the tissue was, how it was prepared) and are supplied by you --
see [Completing the metadata](#completing-the-metadata).

| Section | Field | Type | Ontology |
|---------|-------|------|----------|
| (root) | `schema_version` | `MAJOR.MINOR.PATCH` string, required | -- |
| `sample` | `organism`, `organism_term` | text + term | NCBITaxon |
| | `organism_part`, `organism_part_term` | text + term | UBERON |
| | `condition` | text | -- |
| | `sample_growth_conditions` | text | -- |
| `preparation` | `sample_stabilisation` | text | -- |
| | `tissue_modification` | text | -- |
| | `matrix`, `matrix_term` | text + term | CHEBI |
| | `matrix_application` | text | -- |
| | `solvent` | text | -- |
| `ms_analysis` | `polarity`, `polarity_term` | `"positive"`/`"negative"` + term | PSI-MS |
| | `ionisation_source`, `ionisation_source_term` | text + term | PSI-MS |
| | `analyzer`, `analyzer_term` | text + term | PSI-MS |
| | `instrument_model` | text | -- |
| | `manufacturer` | text, the source's own spelling | PSI-MS (`MS:1001269`) |
| | `serial_number` | text, identifies one physical machine | PSI-MS (`MS:1000529`) |
| | `detector_resolving_power` | `{value, at_mz}` | -- |
| | `pixel_size_um` | `{x, y}`, **required** | -- |
| | `n_spectra` | integer, the spectra the source holds as its reader counted them | -- |
| | `ion_mobility` | `{present, separation, separation_term, unit_term, range_lower, range_upper, num_scans, resolved_table, grid}` | PSI-MS (`MS:1002815` / `MS:1002476`, unit `MS:1002814`) |
| | `fragmentation` | `{present, ms_level, constant_across_pixels, merges_precursors, dissociation_term, windows, resolved_table}` | PSI-MS (`MS:1000511`; windows `MS:1000827` / `828` / `829`, `MS:1000045`, `MS:1000133`) |
| `acquisition` | `acquisition_datetime` | ISO 8601 `YYYY-MM-DDThh:mm:ss[.fff]`, with a UTC offset only when the source recorded one | -- |
| | `laser_power_percent` | number, percent of the laser's range | -- |
| | `laser_frequency_hz` | number, hertz | IMS (`IMS:1006000`) |
| | `shots_per_pixel` | integer, laser shots summed into one pixel's spectrum | IMS (`IMS:1006001`) |
| | `method_file` | text, the acquisition method's file name (never a path) | -- |
| `calibration` | `calibration_datetime` | ISO 8601 at the source's precision (`YYYY-MM-DDThh:mm[:ss[.fff]]`), with a UTC offset only when the source recorded one | -- |
| | `recalibrated` | boolean, a calibration made after the acquisition replaced the one it ran under | -- |
| | `original_calibration_datetime` | as `calibration_datetime`, for the calibration a recalibration replaced | -- |
| | `software`, `software_version` | text, what made the calibration, in the source's words | -- |
| | `n_reference_peaks` | integer, the reference peaks the calibration was fitted to | -- |
| | `mz_standard_deviation_ppm` | number, their residual in ppm (see below); never zero | -- |
| | `lock_mass_corrected` | boolean, the m/z values were corrected against a lock mass | -- |
| `processing` | list of `{name, action_term, software {name, version, uri}, parameters}` | ordered steps, oldest first | PSI-MS data processing action (`MS:1001485` for `m/z calibration`) |
| `provenance` | `thyra_version` | text, required | -- |
| | `source_format` | `"imzml"`, `"bruker"`, ... | -- |
| | `source_path` | text, the source's path in a store block and its name alone in a document | -- |
| | `pixel_size_source` | `"automatic"` / `"manual"` / `"default"` | -- |

An ontology term is always the pair
`{"accession": "MS:1000075", "name": "matrix-assisted laser desorption ionization"}` --
a [CURIE](https://www.w3.org/TR/curie/) plus the term's label, so the block
is readable without resolving anything.

`pixel_size_um` is the one acquisition field that is required: conversion
refuses to run without a pixel size, so a document without it describes no
store Thyra ever wrote. It is a pair because a raster need not be square, and
the same pair now reaches the rest of the store -- the root attrs, the
coordinate system and its affine, the element transforms and
`obs["spatial_x"]`/`["spatial_y"]`. This block used to be the only place the y
pitch survived.

!!! note "Unknown fields are rejected"
    Validation refuses keys the schema does not define, so a typo fails
    loudly instead of becoming an unread field. Anything genuinely
    vendor-specific already has a home in the sections beside this block
    (`format_specific`, `acquisition_params`, `raw_metadata` -- see
    [Output Format](output-format.md#provenance)).

### What is auto-populated

Auto-population is deliberately honest: a field the source does not report
is left unset, never guessed. The only inferences are facts that follow
from the format itself.

| Source | polarity | ionisation source | analyzer | instrument model |
|--------|----------|-------------------|----------|------------------|
| imzML | -- | -- | from the `<analyzer>` component cvParam | from the instrumentConfiguration (model term or `MS:1000031` value) |
| Bruker `.d` | from `Frames.Polarity`, when every frame agrees | MALDI, when the laser tables are present | TOF (timsTOF-family formats) | from the DB |
| PHI ToF-SIMS | from the header | SIMS | TOF | platform name |
| Waters `.raw` | -- | -- | -- | from `_HEADER.TXT` |

Bruker `.d` also fills `ion_mobility`: `present: true` for a TDF acquisition
(TIMS engaged), with the acquired 1/K0 range and the ramp length in scans,
and `present: false` for TSF. Other formats leave the field unset, which
means "not reported", not "no mobility". The MSI table is always summed over
the ramp; how it was summed is the `tdf_spectrum` parameter of the
`conversion` processing step (see [Supported Formats](supported-formats.md#bruker-timstof)).
`resolved_table` names the mobility-resolved sibling table when one was
written beside the summed table, and `grid` (`{law, lower, upper,
n_channels}`) describes the common mobility grid such a table was binned onto
when it was built from per-pixel mobility values (`--mobility-grid`); both
are unset otherwise. `grid` is the only thing in the store that says which of
the two mechanisms filled that table -- reading it off a shared feature axis
bins nothing -- and it is a description, not a structural difference: the
table is the same shape either way.
The arrays that describe the axis itself and the mass-mobility heatmap live
outside this block, in `uns["mobility_axis"]` and `uns["mobility_heatmap"]`
(see [Output Format](output-format.md#ion-mobility)): this block is versioned
and carries no arrays.

Bruker `.d` also fills `fragmentation` from `Frames.MsMsType` and whichever
precursor table the acquisition uses -- `PasefFrameMsMsInfo` for PASEF frames,
`FrameMsMsInfo` for single-precursor ones -- in a TSF acquisition as well as
a TDF one. What stays TDF-only is the demultiplexed sibling table, which
needs the mobility scan ranges only PASEF records. `present: false` records a survey
acquisition; the field is left unset when the database cannot be asked at all,
which means "not reported", not "MS1". `windows` is stored as a JSON string
(a list of objects does not round-trip through AnnData/zarr) and decoded by
`read_msi_metadata_blocks` and `thyra validate`. The field names and terms are
mzPeak's, so an archive and a store describe a precursor the same way; see
[Output Format](output-format.md#fragmentation-msms) for the array block
beside it and for what `merges_precursors` means for the stored spectrum.
Waters `.raw` fills it too, from the precursor m/z MassLynx reports per scan:
a file whose converted functions carry one reports it here, and a file that
also holds MS1 functions converts those alone, reports MS1, and lists the
MS/MS functions it left out under `excluded_functions` in the Waters-specific
block. The MS level MassLynx reports is not used on its own, because on a
raster it split across functions the level is a chunk artefact (see
[Supported Formats](supported-formats.md#which-functions-hold-the-image)).

The `acquisition` section is filled the same way, from the vendor values
the reader already keeps in `uns["acquisition_params"]`, and only where
the value's meaning and unit were checked against an acquisition:

| Source | `acquisition_datetime` | `laser_power_percent` | `laser_frequency_hz` | `shots_per_pixel` | `method_file` |
|--------|------------------------|-----------------------|----------------------|-------------------|---------------|
| Bruker tsf/tdf | `GlobalMetadata.AcquisitionDateTime`, with offset | `MaldiFrameInfo.LaserPower` | `MaldiFrameInfo.LaserRepRate` | `MaldiFrameInfo.NumLaserShots` | `GlobalMetadata.MethodName`, name only |
| Bruker solariX | `Properties.AcquisitionDateTime`, with offset | `Spectra.LaserPower` | `Spectra.LaserRepRate` | `Spectra.NumSummations` (it equals the method's `NumLaserShots`) | the `*.m` directory's name |
| Bruker rapiflex | -- (the info file's `Start Time` format is unverified) | -- (unit unverified) | -- | `Number of Shots` | `Method` |
| PHI ToF-SIMS | `AcqFileDate`, no offset | -- (ion gun) | -- | -- | -- |
| Waters `.raw` | MassLynx's acquisition date and time, no offset | -- | -- | -- | `$$ MS Method` in `_header.txt`, name only |
| imzML | -- | -- | -- | -- | -- |

A Bruker per-frame value that varied across the acquisition is reported in
`acquisition_params` as a `[min, max]` pair and leaves the schema field
unset: the section states one value per acquisition or none. A timestamp
in a format the builder does not parse stays unset too, with the raw
string untouched beside it. No time zone is ever assumed: PHI and Waters
record none, so their timestamps carry none.
`resolved_table` names the demultiplexed sibling table when one was written,
mirroring `ion_mobility.resolved_table`, so both kinds of sibling are
discoverable from this block alone.

A precursor's m/z here and in `var["precursor_mz"]` is the **isolation window
target** (`MS:1000827`) -- the m/z the quadrupole was set to -- and not
`MS:1000744`, a selected ion whose m/z was measured. mzPeak keeps the two in
separate files for the same reason; do not read either as a monoisotopic mass.

`detector_resolving_power` is filled when the source states both the value
and the m/z it is quoted at (an extractor reports them as `resolving_power`
and `resolving_power_at_mz`); one without the other is not written, because
a resolving power is not comparable without its reference m/z. No shipped
reader states the pair today, so for the formats Thyra converts it is
supplied by you like the fields below.

`n_spectra` is the reader's own count of the spectra the source holds -- the
ones present, not the positions the raster covers, and one region's when one
was read. It is left unset where the reader did not count (a PHI preview
decodes no events, so `thyra metadata` has no count for PHI), and the table can
have fewer rows than it says: a spectrum with nothing in it is not stored.

The `calibration` section describes how the source's m/z values were
calibrated. "The calibration" is the one they rest on: the most recent
recalibration when there is one, otherwise the calibration the acquisition
ran under. Every field comes from what the vendor states, and only where its
meaning was checked against real acquisitions:

| Source | `calibration_datetime` | `recalibrated` | `software`, `software_version` | `n_reference_peaks`, `mz_standard_deviation_ppm` | `lock_mass_corrected` |
|--------|------------------------|----------------|--------------------------------|--------------------------------------------------|-----------------------|
| Bruker tsf/tdf | `CalibrationInfo.CalibrationDateTime` in the analysis database, with offset; for recalibrated data the latest `calibration.sqlite` state's time, with the database's as `original_calibration_datetime` | from `calibration.sqlite`: its first state is written at acquisition, so only a second is a recalibration; unset without the file | `CalibrationSoftware`, `CalibrationSoftwareVersion`; the latest state's source and version when recalibrated | the reference masses with a corrected mass, and `MzStandardDeviationPPM`; unset for a recalibration, since no recalibrated acquisition was at hand to check what its state records | -- |
| PHI ToF-SIMS | the appended recalibration's `BlockAppendedDate`, no offset; unset for the header's calibration, which has no date | an appended block with new coefficients | -- | the calibrants of the calibration in force, from their measured and theoretical m/z | -- |
| Waters `.raw` | `$$ Cal Date` and `$$ Cal Time` in `_header.txt`, to the minute, no offset | -- | -- | -- (`Cal StdDev` is 0 on every file read) | MassLynx `isLockmassCorrected` |
| imzML, solariX, rapiflex, mzPeak | -- | -- | -- | -- | -- |

`mz_standard_deviation_ppm` is the square root of the reference peaks'
summed squared m/z errors after the calibration, in ppm, over one less than
the number of peaks. That is Bruker's `MzStandardDeviationPPM`, recomputed
from the file's own arrays to six decimals on four acquisitions, and PHI's
calibrants give the same statistic: refitting their theoretical m/z against
the flight times the stated coefficients imply reproduces the coefficients
exactly, so each calibrant's measured m/z is its position under the fit and
the difference its residual. A value is never written where it has nothing
to say. A standard deviation of zero is left unset -- a timsTOF calibrated on
two peaks fits both exactly, and the online lock-mass state a MALDI timsTOF
writes records `0.000000` against a single reference mass and no measured
mass at all -- and so is one over fewer than two peaks. The count stands on
its own.

Three things are deliberately not in the section. On a timsTOF, the first
state of `calibration.sqlite` is an online lock-mass calibration, and whether
Bruker's library applies it depends on the file: a TSF applies it when opened
with `--use-recalibrated` (the default), a TDF gives the same m/z with it or
without it, measured on one acquisition of each. Which calibration a
conversion applied is therefore a processing step (below), and
`lock_mass_corrected` is left unset for timsTOF rather than stated for data it
may not describe. The mobility calibration's keys are not read: on the one
imaging TDF with TIMS engaged its measured voltages are all zero and its
standard deviation is 3578 %. And no person reaches the section:
`CalibrationUser` and `MobilityCalibrationUser` are not read, and neither is
the reference list's name, which is text the lab chose.

Everything else -- organism, tissue, condition, matrix -- cannot come from
a raw file and stays empty until you provide it.

### Processing history

`processing` is the dataset's processing provenance, modeled on
[mzQC](https://github.com/HUPO-PSI/mzQC): an ordered list of steps, each
naming the software that performed it and the parameters it ran with, and
the PSI-MS data processing action it is (`action_term`) where one exists.
The converter records its own steps -- `conversion` always, `m/z calibration`
when the reader turned flight times or digitiser indices into m/z with a
calibration it chose, and `mass axis resampling` with the resolved resampling
parameters when resampling was enabled. Downstream tools (normalisation, peak
picking, annotation) append theirs when they modify the store.

```python
[
  {"name": "conversion",
   "software": {"name": "thyra", "version": "4.1.0"}},
  {"name": "m/z calibration",
   "action_term": {"accession": "MS:1001485", "name": "m/z calibration"},
   "software": {"name": "thyra", "version": "4.1.0"},
   "parameters": {"calibration": "appended"}},
  {"name": "mass axis resampling",
   "software": {"name": "thyra", "version": "4.1.0"},
   "parameters": {"method": "nearest_neighbor", "target_bins": 50000,
                  "reference_mz": 1000.0}}
]
```

The `m/z calibration` step is how a store says which calibration it holds,
as distinct from what the source states in `calibration`: a conversion can be
told to apply one the source does not consider current. Its parameters are
in the reader's own terms:

| Source | Parameter | Meaning |
|--------|-----------|---------|
| PHI ToF-SIMS | `calibration`: `"appended"` / `"header"` | the coefficients Thyra computed m/z from: the recalibration appended to the file whenever there is one, unless the reader's `use_appended_calibration` was turned off |
| Bruker tsf/tdf | `use_recalibrated_state`: `true` / `false` | the option Bruker's library was opened with (`--use-recalibrated` / `--no-recalibrated`). The library applies the calibration itself, and what it does with a stored state differs by file type (see the `calibration` section above), so the option is recorded as the option |

Other readers take the m/z values the source stores, apply nothing and
record no step.

On disk the list is stored as a JSON string (AnnData/zarr cannot
round-trip a list of objects -- the same reason `uns["regions"]` is
JSON); `read_msi_metadata_blocks` and `validate_document` decode it
transparently.

---

## PSI CV alignment

Beyond the per-document `*_term` values, every schema **field** that
instantiates a PSI CV concept is bound to that concept's accession in
the JSON Schema itself (a `cv` annotation on the field definition), so
the claim "Thyra's output is annotated with the same PSI CV terms as
the raw file it came from" is machine-checkable from the committed
artifact alone:

| Field | CV concept |
|-------|-----------|
| `ms_analysis.pixel_size_um.x` | `IMS:1000046` pixel size (x) |
| `ms_analysis.pixel_size_um.y` | `IMS:1000047` pixel size y |
| `ms_analysis.polarity` | `MS:1000465` scan polarity |
| `ms_analysis.ionisation_source` | `MS:1000008` ionization type |
| `ms_analysis.analyzer` | `MS:1000443` mass analyzer type |
| `ms_analysis.instrument_model` | `MS:1000031` instrument model |
| `ms_analysis.manufacturer` | `MS:1001269` instrument vendor |
| `ms_analysis.serial_number` | `MS:1000529` instrument serial number |
| `ms_analysis.detector_resolving_power` | `MS:1000800` mass resolving power |
| `acquisition.laser_frequency_hz` | `IMS:1006000` repetition rate |
| `acquisition.shots_per_pixel` | `IMS:1006001` laser shots per spectrum |

On the input side, every imzML file-description cvParam is preserved in
`uns["raw_metadata"]["cvParams"]` **with its accession** (and unit
accession where the source set one) -- the name alone cannot be resolved
back to the CV concept. The list is stored as a JSON string (see
[Output Format](output-format.md#provenance) for why); `json.loads`
hands back the list of terms. Polarity declared there (`MS:1000130` /
`MS:1000129`) auto-populates the schema field.

### Candidate CV terms

Several imaging concepts this schema needs have no CV term yet. They
are tracked in `CANDIDATE_CV_CONCEPTS` as the vocabulary to raise in
the PSI/mzPeak imaging discussions, so this schema and the future
standard converge:

- pixel size semantics (raster pitch vs laser spot vs binned size)
- pixel size provenance (measured vs user-supplied vs default)
- coordinate origin and axis handedness
- stage offset of the raster origin
- ROI / acquisition region identity
- missing / empty pixel semantics
- continuous-vs-processed source provenance after conversion
- mass axis resampling provenance (method, axis law, target bins)
- acquisition start timestamp (mzML has only the run's `startTimeStamp`
  attribute; `MS:1000747` is the completion time)
- laser power as a percentage of the instrument's range (`MS:1000846`
  pulse energy is in joules)
- acquisition method identity (`MS:1002128` names a method file format,
  not the method)
- number of spectra of any MS level (`MS:4000059` and `MS:4000060` count MS1
  and MS2 spectra separately)
- when an m/z calibration was made, whether a recalibration replaced it, and
  what software made it (`MS:1001485` m/z calibration is a processing action
  with no attributes; `MS:1003200` software version is scoped to spectral
  libraries)
- an m/z calibration's fit: its residual in ppm and how many reference peaks
  it rests on (`MS:1000014` accuracy is an analyzer attribute, `MS:4000072`
  the error of one identified ion)
- lock-mass correction of the m/z values

`MS:1001485` itself is used, as the `action_term` of the `m/z calibration`
processing step -- the one place the schema names a calibration as
something done rather than something stated.

### LinkML rendering

The schema is also rendered as LinkML at
`thyra/metadata/schema/msi_metadata.linkml.yaml` -- classes, slots,
required flags, `slot_uri` CV bindings and the polarity enum with
`meaning:` accessions. The pydantic models remain the source of truth
for what Thyra writes and validates; the YAML is the discussion artifact
for LinkML-native settings (the PSI/mzPeak imaging work, the planned
spec repository), and a unit test keeps the two from drifting. A later
migration to LinkML-as-source changes no field names and nothing on
disk.

---

## var column conventions

The MSI table's `.var` column names are fixed by the spec so every
consumer can rely on one spelling:

| Column | Written by | Meaning |
|--------|-----------|---------|
| `mz` | every converter, **required** | The common mass axis. Numeric, finite, strictly increasing -- except on a sibling table, where the pair is what is unique and sorted: `(mz, mobility)` on a mobility-resolved one, `(precursor_mz, precursor_mobility, mz)` on a demultiplexed MS/MS one, where `precursor_index` identifies the block. |
| `mobility` | the converter, on mobility-resolved tables only | The feature's ion mobility (1/K0 or drift time; see `uns["mobility_axis"]`). Its presence is what marks the table as mobility-resolved. |
| `precursor_mz` | the converter, on demultiplexed MS/MS tables only | The isolated m/z the feature's fragments came from (see `uns["msms_schedule"]`). Its presence is what marks the table as demultiplexed; such a table never carries `mobility` as well. |
| `mz_index` | the converter, on sibling tables only | Column of the feature's m/z on the summed MSI table's axis |
| `mobility_index` | the converter, on mobility-resolved tables only | The feature's position on the table's mobility axis: the rank of its mobility among the distinct values a shared-axis source lists, or the channel it fell in on a common grid (`ion_mobility.grid`) |
| `precursor_mobility` | the converter, on demultiplexed MS/MS tables only | The 1/K0 the precursor was isolated at (the middle of its mobility window). What tells two precursors sharing an m/z apart -- an isomer pair -- so they are never merged |
| `precursor_index` | the converter, on demultiplexed MS/MS tables only | The precursor's position in **this store's** precursor axis, and the identity of its column block. Means nothing outside the store: align two stores on `(precursor_mz, precursor_mobility)` |
| `formula` | annotation tools | Molecular formula of the annotation |
| `adduct` | annotation tools | Adduct, e.g. `+H`, `-H`, `+Na` |
| `annotation_source` | annotation tools | Tool/database that produced the annotation |
| `fdr` | annotation tools | False discovery rate of the annotation |

Readers may add extra per-channel columns (a non-m/z native axis such as
flight time stays alongside `mz`), and annotation tools may add columns
beyond these -- but the reserved names above must never be reused with a
different meaning. `thyra validate` checks the `mz` contract on every
table of a store.

---

## Storage contract

- The block lives at `table.uns["msi_metadata"]`, in every table the
  store has. This location is stable, like `uns["essential_metadata"]`
  and the `coordinate_systems` root attribute; consumers read it from
  here and nowhere else.
- It is written identically by every converter write path (the uns
  parity tests assert this), and contains no timestamps, so converting
  the same input twice produces the same block.
- Sections with nothing in them are omitted rather than written empty,
  following the store-wide convention.

## Versioning

`schema_version` follows semantic-version rules:

- Adding an optional field bumps the **minor** version.
- Renaming, removing, retyping, or making a field required bumps the
  **major** version.
- A validator implementing major version N rejects documents with a
  different major version, accepts older minors, and warns on newer
  minors.

Versions so far: 0.1.0 (initial), 0.2.0 (`ms_analysis.ion_mobility`
added), 0.3.0 (`ion_mobility.resolved_table` and `ion_mobility.grid` added),
0.4.0 (`ms_analysis.fragmentation` added), 0.5.0
(`fragmentation.resolved_table` added), 0.6.0 (the `acquisition` section
added), 0.7.0 (`ms_analysis.manufacturer` and `ms_analysis.serial_number`
added), 0.8.0 (the `calibration` section, `ms_analysis.n_spectra` and
the processing steps' `action_term` added).

The JSON Schema rendering is published at a fixed, versioned address,
which is also its `$id`:

```
https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/schema/0.8.0/msi_metadata.schema.json
```

Every version gets its own folder under that path and a published folder
is never edited; there is deliberately no `latest`. The LinkML source of
the same version is served beside it as `msi_metadata.linkml.yaml`. A
program that writes the document without Thyra validates against that
address; see [Writing the Metadata Document](writing-the-metadata-document.md).

The same file is committed at
`thyra/metadata/schema/msi_metadata_schema_v0_8.json` and ships in the
wheel, so a Python consumer can validate documents offline without
importing Thyra:

```python
from importlib import resources
import json

schema = json.loads(
    resources.files("thyra.metadata.schema")
    .joinpath("msi_metadata_schema_v0_8.json")
    .read_text()
)
```

Unit tests keep the package artifact in sync with the models and the
published copy in sync with the package artifact; regenerate both with
`python -m thyra.metadata.schema.generate` after a model change. A
version bump creates a new folder under `docs/schema/`; the previous
version's folder stays as it is. `docs/schema/SHA256SUMS` records the
SHA-256 of every published file, and a unit test fails when any of them
changes, so regenerating a version that is already published fails the
suite instead of rewriting it. The pull request that adds a version
folder appends that folder's lines to `SHA256SUMS`; the failing test
prints them.

---

## `thyra metadata`

```
thyra metadata INPUT [--merge USER.json] [-o OUT.json]
```

Builds the same document from a **raw** source and writes it out, without
converting anything. `INPUT` is a file or folder in any format Thyra reads;
no spectra are decoded and no vendor SDK is loaded, so it is a header read
whatever the dataset's size. The default output is stdout.

The document is the one a conversion would have stored, with three
differences. Two follow from nothing having been converted: `processing`
is absent, and no sibling table is named. The third follows from where
the document goes: `provenance.source_path` is the source's **name**,
where a store block carries the path it was converted from. A store sits
on the machine that wrote it and the path is provenance a reader can act
on; a document is written to be handed to somebody, and an absolute path
describes a filesystem they do not have while carrying a user directory
to a machine that has no use for it. Everything else -- the sections, the
ontology terms, the JSON shape -- is identical, so a document written
this way and the same dataset's block read back out of a store compare
directly.

```bash
thyra metadata raw_data.d
thyra metadata raw_data.d --merge sample.json -o meta.json
```

An acquisition with no raster has no pixel size, which this schema
version requires, so its document is reported as invalid and the command
exits 1 -- while still writing the document, which is the only kind such
an acquisition has. See
[Design Decisions](design-decisions.md#d23-a-metadata-document-does-not-need-a-pixel-size-proposed).

---

## `thyra validate`

```
thyra validate PATH [--merge USER.json] [--json]
```

`PATH` is a converted `.zarr` store or a standalone metadata `.json`
document. Only the metadata block (and, for stores, the `var` axis) is
read -- validating a 100 GB store is instant. For stores it also checks
the [var column contract](#var-column-conventions) and that every table
carries a metadata block. Exit status is `0` when every document
conforms (warnings allowed) and `1` otherwise, so it can gate CI.

On Windows, a store whose files sit past the 260-character path limit is
read through an extended-length path automatically; without that, Zarr
would return empty terms for the keys it cannot open and the document
would fail validation for the wrong reason (see
[long paths](troubleshooting.md#windows-long-paths)).

Errors mean the document does not conform: structural violations, unknown
PSI-MS/IMS/UO accessions, version incompatibility. Warnings mean it
conforms but says something suspicious: a term label that does not match
its accession, or a term from an unexpected ontology.

```bash
thyra validate output.zarr
# msi_dataset_z0: OK

thyra validate output.zarr --json   # machine-readable report on stdout
```

## `thyra export-metaspace`

```
thyra export-metaspace PATH [--merge USER.json] [--table NAME] [-o OUT.json]
```

Writes the METASPACE submission metadata JSON (default:
`<input>.metaspace.json` next to the input; `-o -` for stdout). Required
fields the store cannot know are emitted empty and reported as warnings on
stderr -- the output is a truthful starting point, never a fabricated
record. The one inference: matrix-free sources (DESI, SIMS) truthfully get
`MALDI_Matrix: "none"`.

```bash
thyra export-metaspace output.zarr
# warning: Sample_Information.Organism is required ... and is not set
# Wrote METASPACE metadata for msi_dataset_z0 to output.metaspace.json
```

## Completing the metadata

The fields only you can know go in a small JSON overlay, merged at
validation or export time with `--merge`:

```json
{
  "sample": {
    "organism": "Mus musculus",
    "organism_term": {"accession": "NCBITaxon:10090", "name": "Mus musculus"},
    "organism_part": "liver",
    "organism_part_term": {"accession": "UBERON:0002107", "name": "liver"},
    "condition": "wildtype"
  },
  "preparation": {
    "sample_stabilisation": "fresh frozen",
    "matrix": "2,5-dihydroxybenzoic acid (DHB)",
    "matrix_term": {"accession": "CHEBI:17189", "name": "2,5-dihydroxybenzoic acid"},
    "matrix_application": "ImagePrep"
  },
  "ms_analysis": {
    "detector_resolving_power": {"value": 130000, "at_mz": 400}
  }
}
```

```bash
thyra validate output.zarr --merge sample.json
thyra export-metaspace output.zarr --merge sample.json
```

The overlay merges key-by-key over the stored block, so it only needs the
fields you are adding. Keep it next to the store; it applies unchanged to
every dataset from the same study.

## Python API

```python
from thyra.metadata.schema import (
    MSIMetadata,             # the pydantic model
    read_msi_metadata_blocks,  # {table_name: block} from a store
    validate_document,       # (model | None, issues)
    to_metaspace,            # (submission_dict, warnings)
)

blocks = read_msi_metadata_blocks("output.zarr")
meta, issues = validate_document(blocks["msi_dataset_z0"])
submission, warnings = to_metaspace(meta)
```
