# CLI Reference

!!! tip "In plain words"
    [Change how Thyra converts](settings.md) explains the common options one
    task at a time. This page lists every option.

```
thyra [OPTIONS] INPUT OUTPUT
```

**INPUT** -- Path to input MSI file or directory (`.imzML`, `.d`, `.raw`, a
Rapiflex acquisition folder, or `.mzpeak`). A folder that holds one `.d` run
is read as that run. A folder that holds several is asked about: type the
run's number, or pipe it in. With no answer to read, the run ends with one
line naming the runs, and exit `1`.

**OUTPUT** -- Path for output `.zarr` directory

`thyra` also has three metadata subcommands, `thyra metadata`,
`thyra validate` and `thyra export-metaspace` -- see
[Metadata subcommands](#metadata-subcommands).

!!! tip "Grouped help"
    `thyra --help` lists every option under a category heading, in this order:
    Conversion, Logging, Ion mobility grid (advanced), Resampling (advanced),
    imzML-specific, Bruker-specific, Waters-specific, Other, and a General
    section holding `--version` and `--help`. The sections on this page follow
    the same order. The one page section with no matching help group is
    [Performance](#performance), which now records only options that were
    removed or that select nothing -- there is nothing left there for `--help`
    to list.

---

## Exit status

| Code | Meaning |
|------|---------|
| `0` | Conversion completed and the output store was written |
| `1` | Conversion failed |
| `2` | Invalid command-line arguments (click usage error) |

A refusal Thyra planned for -- a `.d` directory with no analysis files, a
Waters raster whose stage never moved, a mass range with no extent, a source
in which no pixel carries a spectrum -- prints its message once at `ERROR` and
nothing else. The traceback behind it is kept for `-v DEBUG`, and a traceback
at `ERROR` now means an exception nobody planned for, which is worth reporting
as a bug.

A conversion that would store **nothing** is one of those failures. If no
position carries a spectrum -- an all-zero source, every spectrum empty,
every peak outside a narrowed `--resample-min-mz`/`--resample-max-mz` range,
every peak below `--intensity-threshold`,
every coordinate off the declared grid -- the conversion exits `1` and
leaves no store, rather than exiting `0` with a zarr holding no table. The
message names which of those happened. A multi-slice source with *some*
empty planes is not affected: those planes are dropped with a warning and
the rest of the store is written.

The store is written under a temporary name beside the output,
`<output>.zarr.incomplete-<8 hex>`, and renamed to the output name as the
last step. So the output name holds a finished store or nothing, however the
run ends. A failed conversion renames the partial store to
`<output>.zarr.failed`, so the output path stays free for a retry. A run that
is killed outright (an out-of-memory kill, a crash) leaves only the
`.incomplete-` folder and its scratch; delete them, or run the same command
again at once. If the final rename fails, the run exits `1` and its error
names the temporary folder, whose store is complete. An
interrupted conversion is a failed one: `Ctrl-C` exits `1` and takes the
same path, so the partial store is moved aside and the scratch directories
the conversion built (`.thyra_summed_*` and its siblings, next to the
output) are removed rather than left behind. A `Ctrl-C` after the store is
saved comes too late to matter: the store is complete, stays in place, and
`thyra` exits `0`. This makes `thyra` safe to
chain in a shell script or CI job:

```bash
thyra input.imzML output.zarr && python analyse.py output.zarr
```

---

## Conversion

| Option | Default | Description |
|--------|---------|-------------|
| `--pixel-size FLOAT` | auto-detect | Pixel size in micrometers, applied to **both** axes. Auto-detection keeps the source's own x and y pitch separately, so an anisotropic raster stays anisotropic; passing this declares it square |
| `--region TEXT` | all | Convert one region of a timsTOF `.d`, by `.mis` Area Name or by DB RegionNumber. Checked against the dataset's own region list whether it has one region or several; a value that matches no Area Name is read as a RegionNumber, and says so. A name two Areas share is refused: pass the RegionNumber. On a Waters `.raw` it takes the number of a region the method file declares, counted from 1. Every other format refuses it; [`list_regions`](api.md#listing-regions) lists the regions of any input |
| `--resample / --no-resample` | enabled | Mass axis resampling. With `--no-resample`, each resampling option given beside it is logged as ignored |
| `--include-optical / --no-optical` | enabled | Include optical images in output |
| `--mobility-table / --no-mobility-table` | enabled | Also write the mobility-resolved sibling table when the source shares one set of (m/z, ion mobility) features across pixels (see [Output Format](output-format.md#ion-mobility)) |
| `--mobility-heatmap / --no-mobility-heatmap` | enabled | When the source has an ion mobility dimension, store the mean mass-mobility frame on the summed table as `uns["mobility_heatmap"]`; on a Bruker TDF it is fed from the summed table's own passes, on other sources (imzML, Waters drift time) it is one extra pass, and a processed imzML first reads its mobility values for the range (see [Output Format](output-format.md#ion-mobility)) |
| `--mobility-grid / --no-mobility-grid` | **disabled** | When the source carries ion mobility per pixel rather than as a shared feature list (Bruker TDF, processed imzML, Waters drift time), bin the point cloud onto a common mobility grid and write the same mobility-resolved sibling table, a much larger table built out of core so any acquisition fits. On a Bruker TDF it is fed from the summed table's own two passes (no extra read of the source) and its marginal reproduces the summed table exactly under the default `--tdf-spectrum scan_sum`; on a processed imzML or a Waters run it takes extra passes, and its marginal is exact when points stay on their nearest mass bin (on Waters, under `--waters-spectrum profile`). Ignored, with a warning, on a continuous imzML, which has a shared feature list, and on a Bruker TSF file, which has no ion mobility (see [Output Format](output-format.md#the-same-table-from-a-common-mobility-grid)) |
| `--msms-table / --no-msms-table` | enabled | When the source isolates several precursors per pixel in disjoint mobility slices (Bruker PASEF), also write them split apart as a demultiplexed sibling table; fed from the summed table's own two passes (no extra read of the source); the split adds back up to the summed table exactly under the default `--tdf-spectrum scan_sum` (see [Output Format](output-format.md#demultiplexed-msms-table)) |

!!! warning "`--format` is deprecated"
    `spatialdata` was the only value `--format` ever accepted, and it is the
    only output format Thyra writes, so no invocation could ever change what
    the flag selects. It is still accepted, so existing scripts keep running,
    but it is no longer listed in `--help`. Any other value still fails
    before any data is read.

### Examples

```bash
# Basic conversion -- format, pixel size, and resampling all auto-detected
thyra input.imzML output.zarr

# Specify pixel size manually (when metadata is unavailable)
thyra input.imzML output.zarr --pixel-size 25

# Convert one region, named the way flexImaging names it
thyra data.d output.zarr --region 03

# Or by the database's own region number
thyra data.d output.zarr --region 0

# Skip optical images
thyra data.d output.zarr --no-optical
```

!!! note "Naming a region"
    `--region` takes a string, and the value is matched against the `.mis`
    Area Names first. Only if nothing matches is it parsed as a DB
    RegionNumber, which starts at 0. Passing something that is neither is an
    error that lists the area names it did find. A name that two Areas share
    is refused too, with the RegionNumbers that carry it: pass one of those.

    The two numbering schemes do not have to agree, which is exactly why area
    `'03'` need not be RegionNumber 3. For a multi-region dataset Thyra logs
    the whole mapping at `INFO` during startup -- a header line
    `Region mapping (DB RegionNumber -> .mis Area Name):` followed by one
    `RegionNumber <n> -> Area '<name>' (<n> frames)` line per region -- so a
    plain run already tells you which is which, and how many frames each
    holds, before you pick one.

---

## Logging

| Option | Default | Description |
|--------|---------|-------------|
| `-v, --log-level LEVEL` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL` |
| `--log-file PATH` | none | Write logs to file. A missing folder is created; a path that is a folder, or lies under a file, is refused with exit `2` before any work |

### Examples

```bash
# Verbose output -- shows pixel size detection, resampling config, timing
thyra input.imzML output.zarr -v DEBUG

# Save logs to file for later review
thyra input.imzML output.zarr --log-file conversion.log
```

!!! tip "Debugging conversions"
    When something looks wrong in the output, re-run with `-v DEBUG --log-file
    debug.log`. The log will contain pixel size detection details, resampling
    parameters, region info, and timing for each step.

---

## Ion Mobility Grid (Advanced)

These size the grid `--mobility-grid` bins onto. They do nothing without it,
and passing one without it is logged at `WARNING`. A continuous imzML
with a shared feature list needs no grid: there `--mobility-grid` and these
are ignored, with a `WARNING`. So are they on a Bruker TSF file, which has no
ion mobility. With `--no-mobility-table` the grid is cancelled, with a
`WARNING`. A bound you state drops the points outside it.

| Option | Default | Description |
|--------|---------|-------------|
| `--mobility-bins INTEGER` | the heatmap's channels | Mobility channels the grid divides the range into. At least 1. Unset, the grid takes the heatmap's own channels: 256, or one per drift bin on a Waters run (200 bins measured) |
| `--mobility-min FLOAT` | axis minimum | Lower edge of the grid, in the axis unit (1/K0 for TIMS, ms for Waters drift time). May be zero or negative -- it is a position on an axis, not a quantity -- but must be below `--mobility-max`. Points below it are dropped and counted in a `WARNING` |
| `--mobility-max FLOAT` | axis maximum | Upper edge of the grid. Points above it are dropped and counted in a `WARNING` |

### Examples

```bash
# The mobility-resolved table for a TIMS acquisition, at the default
# 256 channels over the range the acquisition actually used
thyra tims_data.d output.zarr --mobility-grid

# A narrower range at finer resolution, giving up the heatmap alignment.
# The coarser mass axis keeps the 512-channel count array under its cap
thyra tims_data.d output.zarr --mobility-grid --mobility-bins 512 --resample-bins 100000 --mobility-min 1.05 --mobility-max 1.25
```

!!! note "The grid table is built out of core"
    Like the summed table, the grid table is built in two passes: the first counts the occupied `(m/z bin, channel)` cells and
    the second scatters every pixel straight into memmapped CSC arrays in a
    scratch directory next to the output (`.thyra_mobility_*`, removed once
    the table is written). On a Bruker TDF those are the summed table's own
    two passes: each frame is read once per pass and the read serves the
    summed spectrum, the heatmap, the grid and the MS/MS split alike (see
    [Design Decisions](design-decisions.md#d5-one-raw-read-per-frame-per-pass-serves-every-table)). Memory is the count array over the grid's span --
    4 bytes x bins x channels, capped at 1 GB -- plus one frame, whatever the
    number of pixels; disk is 12 bytes per stored non-zero while the table is
    being built. Measured on a whole 26,087-pixel timsTOF acquisition on a
    40,000-bin axis (7.3M features, 1.59 billion non-zeros, 18 GB of scratch):
    the table's own phases peaked at 3.2 GB of private memory above the
    interpreter's baseline, in the `var` frame rather than in the matrix.

!!! note "A mobility-resolved table wants a coarser mass axis than the default"
    The table's features are `(m/z bin, mobility channel)` pairs, so the mass
    axis' width is multiplied by the channel count (256 on a timsTOF). A default resampled timsTOF axis
    grows with the mass range: 599,146 bins on a usual 50 to 1000 m/z run, so
    256 channels fit but 512 do not. 200 frames of one measured acquisition filled 3.9
    million of the resulting pairs, and a whole 26,000-pixel acquisition
    occupies 21 million. The table is refused when the `var` frame those
    pairs make is projected to take more than half of the free memory (see
    below) (about 338 bytes per pair, measured), and warned about past a
    quarter. Reach for `--resample-bins` when the whole image is wanted
    resolved on a smaller machine.

!!! warning "The defaults are an alignment, and the grid changes the TIC"
    The default channels are exactly what `uns["mobility_heatmap"]` bins
    over -- 256 over the mobility axis' own value range, or one per Waters
    drift bin -- so a box drawn on the heatmap selects grid channels by
    integer index. `--mobility-bins` (even `--mobility-bins 256`),
    `--mobility-min` and `--mobility-max` each give that up.

    The grid is built from raw scans, and its marginal over channels
    reproduces the summed table only when that table was built the same
    way, which the default `--tdf-spectrum scan_sum` is. An explicit
    `vendor_centroid` is kept, said at `WARNING`, and the mismatch is
    recorded in the grid's `uns` block.

---

## Resampling (Advanced)

These options control how spectra are mapped onto a common mass axis. In most
cases the defaults work well -- Thyra auto-detects the instrument type and
chooses an appropriate method, axis law and bin width. **Start with no flags**;
the table after this one says when each option is actually needed.

| Option | Default | Description |
|--------|---------|-------------|
| `--resample-method METHOD` | `auto` | `auto`, `nearest_neighbor`, or `tic_preserving`. A value that contradicts the detector's own choice is applied, and warned about |
| `--mass-axis-type TYPE` | `auto` | `auto`, `constant`, `linear_tof`, `reflector_tof`, `tof`, `orbitrap`, `fticr` |
| `--tof-law A B` | auto | Coefficients of the `tof` width law `sqrt(A m + B m^2)` mDa (`A` in mDa<sup>2</sup>/Da, `B` dimensionless). Only for an instrument Thyra has no pair for; an MRT centroid run, a PHI run or a timsTOF supplies its own |
| `--resample-bins INTEGER` | auto | Number of bins (mutually exclusive with `--resample-width-at-mz`) |
| `--resample-min-mz FLOAT` | auto | Minimum m/z value. Inclusive: a peak sitting exactly on it lands in the first bin |
| `--resample-max-mz FLOAT` | auto | Maximum m/z value. Inclusive, as above; peaks outside the range are dropped, never folded into the edge bins |
| `--resample-width-at-mz FLOAT` | auto | Mass width in Da at reference m/z for physics-based binning |
| `--resample-reference-mz FLOAT` | `1000.0` | Reference m/z for width specification |
| `--mass-axis-plan FILE` | none | A plan file from `plan_mass_axis`, so several data sets share one mass axis. It fixes the axis, so it cannot be combined with `--no-resample` or the axis options above |
| `--resample-gap-tolerance FLOAT` | none | `tic_preserving` only: discard target bins farther than this many Da from any measured m/z, instead of interpolating across the gap. Logged as ignored when the method is `nearest_neighbor` |

### Which of these you actually need

Every instrument Thyra recognises gets a measured default: the method, the
axis law and the bin width all come from the detector that matched the file
(see [Resampling](resampling.md#which-detector-wins)). The flags exist for the
cases the defaults cannot know about.

| You want to | Use | Leave alone |
|---|---|---|
| Convert a file from a recognised instrument | nothing | everything here |
| A finer or coarser axis than the default | `--resample-width-at-mz` (with `--resample-reference-mz`) | `--mass-axis-type`: the law stays the instrument's, only the width moves |
| A fixed number of bins instead of a width | `--resample-bins` | |
| A narrower mass window | `--resample-min-mz` / `--resample-max-mz` | |
| An axis law for an instrument Thyra could not identify | `--mass-axis-type` | Note that a manual axis type also resets the width to the axis type's own default (5 mDa at m/z 1000; 17 mDa at 300 for `linear_tof`), since a width tuned for one law is not a default for another |
| The measured `tof` law on a timsTOF, or on a new TOF you have fitted | `--mass-axis-type tof`, plus `--tof-law A B` only when Thyra has no pair for the instrument | |
| Force interpolation on data Thyra would bin | `--resample-method tic_preserving`, and `--resample-gap-tolerance` if the m/z arrays are sparse | The conversion warns when the method contradicts the detector, and names the tolerance flag. Take it seriously: on a sparse source, interpolation without a tolerance filled the axis and stored 1,400x the non-zeros |
| One axis shared by several data sets | `--mass-axis-plan FILE`, written by `plan_mass_axis` | every axis option above: the plan fixes them |
| The raw m/z values, no common axis | `--no-resample` | |

The bin width of a `tof` axis is set the same way as any other: the width at
the reference m/z. Thyra derives the bins per peak width from it (3 by
default), so there is no separate flag for that quantity.

!!! info "Choosing a resampling method"
    - **`nearest_neighbor`** -- Each target bin takes the nearest original m/z
      value. Correct for **centroid** data, where peaks are discrete masses.
    - **`tic_preserving`** -- Linear interpolation, rescaled so the total ion
      current is unchanged. Bins wide enough to skip source samples are
      rebinned instead, so a coarse axis keeps narrow peaks. Correct for
      **profile** data on a target axis whose bin widths scale the same way
      the source points are spaced.
    - **`auto`** -- Picks `tic_preserving` only where the source grid's law
      is known and the target axis follows it: Bruker rapifleX and other
      Bruker MALDI-TOF data, but not timsTOF (uniform in m/z, onto
      `constant`), and the Waters profile trace (uniform in flight time, onto
      `linear_tof`). Everything else gets
      `nearest_neighbor` -- including profile data from an instrument Thyra
      cannot identify, because the interpolating method is only exact when
      the two axis laws match. See
      [Resampling](resampling.md#which-detector-wins) for the full decision
      table.

!!! warning "A very fine axis is refused before it is built"
    Everything a conversion holds per m/z bin -- the axis, the `var` frame
    with its index, the running intensity totals, the pre-scan's count array
    and anndata's copies during the write -- comes to about 200 bytes for
    one table, measured. Each further table (one per plane on the 2D route)
    adds about 36 bytes per bin, and each region about 20 when a source has
    two or more. The bin count is checked before the axis is materialised,
    and again once the regions and the planes are known. It is compared with
    the memory this process can still take: the machine's free memory, or a
    job, container or cgroup limit when that is lower. Past half of it the
    conversion is refused with the projection and the free memory printed;
    past a quarter it is attempted with a `WARNING`. With 8 GB free, a
    one-table conversion is refused from around 21 million bins. No table
    can count more than 268,435,456 bins, so a wider axis is refused however
    much memory is free. `--resample-bins` and a small
    `--resample-width-at-mz` are the two ways to ask for a fine axis; a
    wide-range source at a fine width reaches the limit without either.

!!! note "Free disk is checked after the pre-scan"
    The matrix is written to temporary files next to the output while the
    store is built: 12 bytes per stored value, 16 for a matrix of more than
    two billion values. After the pre-scan their size is known exactly. A
    conversion whose files need more than the drive's free space is refused,
    with both sizes printed. Past half of it a `WARNING` says the store,
    which is written to the same drive, may not fit.

!!! info "Choosing a mass axis type"
    The axis type determines how bin widths scale with m/z:

    - **`constant`** -- Uniform bin width (Da). Used for Bruker MALDI-TOF data, whose spectra arrive evenly spaced in m/z, and for profile data Thyra cannot identify.
    - **`linear_tof`** -- Width scales as sqrt(m/z). Matches a linear (axial) TOF, whose flight time grows as sqrt(m/z).
    - **`reflector_tof`** -- Width scales linearly with m/z (constant relative resolution). Matches reflector TOF.
    - **`tof`** -- Width follows a measured peak-width law `sqrt(A m + B m^2)`, of which the two above are the limits. See [the two-term TOF law](resampling.md#the-two-term-tof-law).
    - **`orbitrap`** -- Width scales as m/z^(3/2). Matches Orbitrap resolution.
    - **`fticr`** -- Width scales as m/z^2. Matches FTICR resolution.
    - **`auto`** -- Detected from instrument metadata.

### Examples

```bash
# An Orbitrap imzML whose analyser Thyra could not identify: name the axis
# and leave the method on auto, which bins rather than interpolates
thyra input.imzML output.zarr --mass-axis-type orbitrap

# Fixed number of bins
thyra input.imzML output.zarr --resample-bins 50000

# Restrict mass range
thyra input.imzML output.zarr --resample-min-mz 100 --resample-max-mz 1000

# Specify bin width at a reference m/z (physics-based)
thyra input.imzML output.zarr \
    --resample-width-at-mz 0.01 \
    --resample-reference-mz 500
```

---

## Performance

No options: every conversion streams, and every table is CSC.

!!! info "Every conversion streams"
    There is no streaming mode to switch on: every conversion makes two passes
    over the source -- one to count, one to scatter straight into memory-mapped
    CSC arrays -- and writes each table from those arrays, so the matrix is
    never held in RAM whatever the dataset's size. On a small file the second
    pass costs about a third of the run;
    [Design Decisions](design-decisions.md#d11-one-converter) explains why
    there is one converter.

!!! warning "Options kept or removed for old scripts"
    - `--streaming auto|true|false` is still accepted so existing scripts
      keep running, but it selects nothing.
    - `--optimize-chunks` is still accepted but does nothing, logs a
      warning, and will be dropped in a future release. Chunk sizes are
      chosen at write time.
    - `--sparse-format` was removed in v3.22: passing it fails as an unknown
      option. Every table is CSC, the layout an ion image reads down, because
      one m/z across all pixels is one contiguous column. For the row-wise
      layout, convert after reading:

        ```python
        X = sdata.tables["msi_dataset_z0"].X.tocsr()
        ```

---

## imzML-Specific

!!! note "A vendor option on another vendor's file is logged as ignored"
    Each of the format-specific groups below -- imzML, Bruker, Waters -- and
    the ion mobility grid options are accepted on any input and do nothing on
    a source they do not apply to. Passing one on another format is now
    logged at `WARNING`, naming the option and both formats, so a script that
    has been carrying a dead flag says so on the next run rather than
    producing a store byte-identical to the plain one in silence. The Python
    API names a `reader_options` key the same way.

| Option | Default | Description |
|--------|---------|-------------|
| `--spectrum-type TYPE` | `auto` | `auto`, `profile`, or `centroid` -- declare the spectrum representation instead of detecting it. imzML only |

By default Thyra records the representation every converted spectrum states
(`MS:1000127` centroid / `MS:1000128` profile). A spectrum states it itself or
in a parameter group it references. A group no spectrum references counts for
nothing. When the spectra state both, no type is recorded and a warning gives
the counts. The `fileContent` declaration counts only when no spectrum states
a type, and only if it names one. Thyra guesses only when nothing is stated.
`--spectrum-type` overrides all of that.

SCiLS Lab exposes the same control as `--rep_type PROFILE|CENTROID` (2026b User
Guide, p.81).

!!! warning "This changes stored values"
    The representation is not just a label: it feeds instrument detection,
    which picks the mass-axis type, which decides the bin spacing. Overriding
    it changes the axis a dataset is resampled onto, and therefore the numbers
    written to the store. It changes nothing for anyone who does not pass it.

    Use it when a file declares the wrong thing -- that does happen. When the
    override contradicts the type the spectra state Thyra logs a warning
    naming both values; it is in the normal output, so check it and the result.

### Examples

```bash
# The file says centroid but it is really profile data
thyra input.imzML output.zarr --spectrum-type profile

# See what detection concludes first: a plain run logs the detected
# spectrum type and the resampling it chose
thyra input.imzML output.zarr
```

---

## Bruker-Specific

Grouped here because this is where they are reached for. All but
`--intensity-threshold` only apply to Bruker timsTOF `.d` directories (`--tdf-spectrum` to a TDF one), and
`--interactive-calibration` to solariX ones as well; on any other input they
are logged as ignored. `--intensity-threshold` is honoured by every reader,
and is here because continuous-mode Bruker data is what usually needs it.

| Option | Default | Description |
|--------|---------|-------------|
| `--use-recalibrated / --no-recalibrated` | enabled | Apply the newest calibration state in `calibration.sqlite` that maps the acquisition's frames, TSF and TDF alike; off, the analysis database's own calibration. timsTOF only |
| `--interactive-calibration` | off | Display available calibration states, and the one the conversion applies. timsTOF and solariX |
| `--intensity-threshold FLOAT` | none | Minimum intensity filter. **Every format** |
| `--tdf-spectrum {scan_sum,vendor_centroid}` | `scan_sum` | How a TDF (TIMS) frame's mobility scans collapse into one spectrum per pixel. Bruker TDF only; ignored, and named, on a TSF file |

### Examples

```bash
# Use raw (non-recalibrated) m/z values
thyra data.d output.zarr --no-recalibrated

# List the calibration states the file carries. Display only: it cannot
# choose one, and the conversion then runs as it would without the flag
thyra data.d output.zarr --interactive-calibration

# Filter low-intensity signals (useful for continuous-mode Bruker data)
thyra data.d output.zarr --intensity-threshold 100

# TIMS data: reproduce Bruker's centroided spectrum instead of the raw scan sum
thyra tims_data.d output.zarr --tdf-spectrum vendor_centroid
```

!!! note "TDF: the mobility ramp is summed, not sliced"
    Every scan of a TIMS frame is read. `scan_sum` (default) sums every scan
    per digitizer index and keeps all of the ion current, which is also
    Bruker's own per-frame total; `vendor_centroid` is Bruker's frame-level
    peak picker over the full ramp (the same peaks as TSF line spectra and
    SCiLS Lab's import), at 79 to 99 percent of the current and a third to
    two thirds of the points. See [Supported Formats](supported-formats.md#bruker-timstof)
    and [Design Decisions](design-decisions.md#d1-which-spectrum-a-reader-takes).

!!! warning "Intensity threshold"
    The `--intensity-threshold` option drops all peaks below the given value
    **before** writing to zarr. This reduces file size but is irreversible.
    Use with care -- inspect the data with `-v DEBUG` first to choose an
    appropriate threshold.

    With ion mobility, a peak's intensity is its sum over all mobility scans
    of the pixel. Every table keeps or drops the peak by that one sum. The
    store records the threshold as `intensity_threshold`.

    Only real measurements count toward that sum. A NaN, an infinite or a
    negative intensity adds nothing, and neither does a point whose m/z is
    not a finite number.

---

## Waters-Specific

These options only apply when converting Waters `.raw` directories; on any
other input they are logged as ignored.

| Option | Default | Description |
|--------|---------|-------------|
| `--waters-spectrum {centroid,profile}` | `profile` on a SELECT SERIES MRT, `centroid` on every other Waters instrument | What MassLynx hands back for one pixel: the vendor peak picker, or the sampled trace behind it |
| `--waters-lock-mass MZ\|off` | keep what the run carries | Correct the m/z against the run's lock-mass reference scans, given the reference ion's m/z; `off` stores them uncorrected. The run itself is never written to. See [Lock-mass correction](supported-formats.md#lock-mass-correction) |

### Examples

```bash
# An MRT run: the profile trace on a linear_tof axis at 1.3 mDa, no flags needed
thyra mrt_run.raw output.zarr

# The same run through the vendor peak picker instead (the tof law, 1.75 mDa)
thyra mrt_run.raw output.zarr --waters-spectrum centroid

# A Synapt run's profile trace; the bin width follows its own digitiser
thyra synapt_run.raw output.zarr --waters-spectrum profile

# A MALDI run corrected against its red phosphorus lock-mass spot
thyra synapt_run.raw output.zarr --waters-lock-mass 526.5545
```

!!! note "When the profile is worth its size"
    MassLynx centroiding reports a single peak wherever the sampled trace has
    two maxima a few mDa apart, and puts it partway between them. On a
    SELECT SERIES MRT brain section, four separate 8-9 mDa doublets between
    m/z 760 and 830 -- including the <sup>13</sup>C<sub>2</sub> isotopologue of
    PC 34:1 [M+K]<sup>+</sup> against PC 34:0 [M+K]<sup>+</sup> -- came back as
    one centroid in 96-100% of the pixels that resolved them, 4-8 ppm from
    either true mass. The profile keeps them apart, for about 3 times the
    store, and is not slower: MassLynx centroids on demand, so a profile read
    skips that work. On a Synapt G2-Si the centroider was found to merge
    nothing the profile resolves, which is why only the MRT defaults to it;
    see [Supported Formats](supported-formats.md#waters-masslynx).

!!! note "The other reason to pass `profile`"
    A long imaging run outgrows one `_FUNC*.DAT` file, and MassLynx continues
    it in a new *function*. It will not centroid the last of them, so a
    centroid conversion leaves that chunk out and logs how many pixels
    (16.6% of one Synapt G1 run) and this flag. Reading the trace converts
    every chunk, because then they all come back the same way -- budget the
    disk for it, since the profile store for that run is estimated at 74 GB
    against 241 MB. See
    [Which functions hold the image](supported-formats.md#which-functions-hold-the-image).

---

## Other

| Option | Default | Description |
|--------|---------|-------------|
| `--dataset-id TEXT` | `msi_dataset` | Dataset identifier used in element keys. Letters, digits, underscores, dots and hyphens only. Not a leading `.` or `__`, not a Windows device name such as `nul`, and at most 216 plain characters: it names every element in the store and every folder that holds one. Checked before any pass over the source |
| `--sample-id TEXT` | none | Names the sample in the store's metadata (`msi_metadata.sample.sample_id`). Never read from the file, and not part of any element name. Stored as typed, so never a patient's name. Checked before any pass over the source |
| `--handle-3d` | off | Process as 3D volume instead of 2D slices |
| `--z-spacing FLOAT` | in-plane pixel size (x) | Distance between consecutive slices, in um. Only used with `--handle-3d` |

### Examples

```bash
# Custom dataset ID (affects table and image key names)
thyra input.imzML output.zarr --dataset-id hippocampus
# -> table key: hippocampus_z0, TIC key: hippocampus_z0_tic

# Combine z-slices into a single 3D table
thyra volume.imzML output.zarr --handle-3d

# 20 um sections acquired on a 50 um raster
thyra volume.imzML output.zarr --handle-3d --pixel-size 50 --z-spacing 20
```

### Set `--z-spacing` whenever you know it

`--pixel-size` is the **in-plane** raster pitch. It says nothing about how far
apart consecutive slices are: that distance is set by the microtome that cut the
sections, not by the stage that rastered them, and the two match only by
coincidence.

With no `--z-spacing`, Thyra reuses the in-plane pitch (the x one, when the
raster is anisotropic and there is no single in-plane pitch to reuse), warns,
and records `z_spacing_source: "assumed_isotropic"` in the store so the guess
stays distinguishable from a measurement. The volume's voxel values are
unaffected either way — what is wrong is its depth, so a viewer reading it in micrometres
renders the stack squashed or stretched along z.

There is no way to detect this from the data. imzML has no term for slice
spacing, and 3D MSI acquisitions are frequently non-consecutive sections, so
the number has to come from whoever cut them.

!!! warning "It does nothing without `--handle-3d`"
    Without 3D handling each slice is written as its own 2D image and there is
    no z axis to space out. Passing `--z-spacing` on its own is logged as
    ignored rather than silently accepted.

!!! warning "`--handle-3d` on a single-slice acquisition"
    One plane has no slice-to-slice distance, so no `z_spacing_um` is applied
    or recorded and a `--z-spacing` given alongside is logged as ignored. The
    store is still written as a volume of one plane, which changes the element
    keys -- the table is `{id}`, not `{id}_z0` -- but not the values.

---

## Metadata subcommands

`validate` and `export-metaspace` take either a converted `.zarr` store
or a standalone metadata `.json` document. Only the metadata block is
read, never the intensity data, so both are instant on stores of any
size. `metadata` starts from the other end: it takes a *raw* source and
writes the document a conversion would have stored. See
[Metadata Schema](metadata-schema.md) for the schema itself.

### `thyra metadata`

```
thyra metadata INPUT [OPTIONS]
```

Writes the `msi_metadata` document for a raw source without converting
it. **INPUT** is a file or folder in any format Thyra reads. For every
format but Waters, no spectra are decoded, no vendor SDK is loaded and
nothing is written beside the input, so this is a header read on a dataset
of any size. A Waters run is described through the MassLynx library, which
reads every scan: it takes a while, and needs Windows or Linux.

| Option | Default | Description |
|--------|---------|-------------|
| `--merge PATH` | none | JSON file overlaid onto the metadata before validation |
| `-o, --output PATH` | stdout | Output file; `-` is stdout too. A missing folder is created; a path that is a folder, or lies under a file, is refused with exit `2` before the source is read |

Exit status is `0` when the document conforms (warnings allowed), `1`
when it does not, `2` for usage errors. The document is written either
way, and the validation issues go to stderr -- so `thyra metadata src
-o meta.json` leaves you a document to look at whatever the exit status
was, and `thyra metadata src > meta.json` keeps stdout clean.

The document names the source rather than locating it: `source_path`
carries `raw_data.d`, not the directory it was read from. A document is
written to be sent to somebody, and nothing in it should have to be
edited out first.

```bash
# Look at what a folder says before spending hours converting it
thyra metadata raw_data.d

# Keep it, with the fields the file cannot know filled in
thyra metadata raw_data.d --merge sample.json -o meta.json
```

!!! info "A source with no image still has metadata"
    The one field the schema requires is the pixel size, because every
    document Thyra had written until now described an image. A Bruker
    `.d` from a run that imaged nothing -- an electrospray acquisition
    on a timsTOF-family instrument -- states its instrument, polarity,
    precursor schedule and acquisition time exactly as an imaging run
    does, and states no pixel size, because there is no
    raster. `thyra metadata` writes that document, reports the one
    validation error, and exits `1`. Converting such a file is refused
    with a message that says so. The split that would make the document
    valid is proposed in
    [Design Decisions](design-decisions.md#d23-a-metadata-document-does-not-need-a-pixel-size-proposed).

### `thyra validate`

```
thyra validate PATH [OPTIONS]
```

Validates the `uns["msi_metadata"]` block (or a JSON document) against
the schema: structure, schema version, and ontology terms. For a store,
it also checks each table's `var` and that its m/z axis matches the
recorded `provenance.mass_axis`.

| Option | Default | Description |
|--------|---------|-------------|
| `--merge PATH` | none | JSON file overlaid onto the metadata before validation |
| `--json` | off | Machine-readable report on stdout instead of text |

Exit status follows the conversion command's convention: `0` when every
document conforms (warnings allowed), `1` otherwise, `2` for usage
errors -- so `thyra validate` can gate a CI pipeline.

A store with no metadata block was either written by a version without
the schema or left by a conversion that stopped before it finished; the
message names both. A folder named `<name>.incomplete-<8 hex>` is the
temporary folder of a conversion that did not finish, and the message says
so.

```bash
thyra validate output.zarr
thyra validate output.zarr --merge sample.json
thyra validate metadata.json --json
```

### `thyra export-metaspace`

```
thyra export-metaspace PATH [OPTIONS]
```

Writes the [METASPACE](https://metaspace2020.org) submission metadata
JSON for a dataset. Required fields the store cannot know (organism,
condition, ...) are emitted empty and reported as warnings on stderr;
supply them with `--merge` or fill them on the submission form.

| Option | Default | Description |
|--------|---------|-------------|
| `--merge PATH` | none | JSON file overlaid onto the metadata before export |
| `--table NAME` | the only one | Table to export when the store has several |
| `-o, --output PATH` | `<input>.metaspace.json` | Output file; `-` for stdout. A missing folder is created; a path that is a folder, or lies under a file, is refused with exit `2` |

```bash
thyra export-metaspace output.zarr --merge sample.json
thyra export-metaspace output.zarr -o -          # print to stdout
```

---

## The other two commands

Installing Thyra puts three commands on your path, not one.

### `thyra-example-data`

Writes the small synthetic imzML dataset the [tutorial](tutorial.md) uses. No
download, no vendor SDK.

```
thyra-example-data [OUTPUT] [OPTIONS]
```

**OUTPUT** -- Path for the `.imzML` (default:
`example_data/synthetic_brain.imzML`). A non-`.imzML` suffix is replaced.

| Option | Default | Description |
|--------|---------|-------------|
| `--pixels-x INTEGER` | `48` | Grid width in pixels, at least 1 |
| `--pixels-y INTEGER` | `36` | Grid height in pixels, at least 1 |
| `--pixel-size FLOAT` | `25.0` | Pixel size in micrometers written to the metadata, above 0 |
| `--mz-bins INTEGER` | `4000` | Number of points on the m/z axis, at least 1 |
| `--seed INTEGER` | `0` | Random seed, at least 0 |
| `--verbose` | off | Verbose output |

```bash
# The tutorial's dataset: 48 x 36 pixels, 4,000 m/z points, seed 0
thyra-example-data example_data/synthetic_brain.imzML

# A larger phantom on a coarser axis
thyra-example-data big.imzML --pixels-x 200 --pixels-y 150 --mz-bins 2000
```

!!! note "`--seed` is what makes the tutorial reproducible"
    The seed fixes the *data*: the same `--seed` gives the same m/z axis and
    the same intensities on every machine, which is why the tutorial can
    print exact numbers and expect yours to match. Change the seed and the
    numbers change with it.

    It does not make the files byte-identical, and is not meant to. Each
    imzML gets a fresh UUID and the SHA-1 of its own `.ibd`, as a real
    acquisition would.

### `thyra-check-ontology`

Validates the CV (controlled vocabulary) terms in an imzML file against
Thyra's ontology cache, and reports any it does not recognise.

```
thyra-check-ontology INPUT [OPTIONS]
```

**INPUT** -- An `.imzML` file, or a directory to check every file in,
subfolders included. A file that cannot be read is named and skipped, and
the command then exits 1. A path that does not exist exits 2.

| Option | Default | Description |
|--------|---------|-------------|
| `--output PATH` | none | Save results as JSON. For one file, this replaces the printed summary. A missing folder is created; a path that is a folder, or lies under a file, is refused with exit `2` before the check runs |
| `--verbose` | off | Log progress at INFO level to stderr |

For a directory, the command prints the file count, the unknown terms and
any file it could not read. Each file's own results are only in the
`--output` file.

```bash
# Check one file
thyra-check-ontology sample.imzML

# Check a whole directory and keep the report
thyra-check-ontology data/ --output ontology_report.json
```
