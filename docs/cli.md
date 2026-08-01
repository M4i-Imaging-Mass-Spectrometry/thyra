# CLI Reference

```
thyra [OPTIONS] INPUT OUTPUT
```

**INPUT** -- Path to input MSI file or directory (`.imzML`, `.d`, `.raw`)

**OUTPUT** -- Path for output `.zarr` directory

!!! tip "Grouped help"
    Run `thyra --help` to see all options organised by category (Conversion,
    Logging, Resampling, Performance, Bruker-Specific, Other), and
    `thyra --version` to print the installed version.

---

## Exit status

| Code | Meaning |
|------|---------|
| `0` | Conversion completed and the output store was written |
| `1` | Conversion failed |
| `2` | Invalid command-line arguments (click usage error) |

A failed conversion renames any partially written store to
`<output>.zarr.failed`, so the output path stays free for a retry and an
incomplete store is never left where a finished one is expected. This makes
`thyra` safe to chain in a shell script or CI job:

```bash
thyra input.imzML output.zarr && python analyse.py output.zarr
```

---

## Conversion

| Option | Default | Description |
|--------|---------|-------------|
| `--format [spatialdata]` | `spatialdata` | Output format |
| `--pixel-size FLOAT` | auto-detect | Pixel size in micrometers |
| `--region INTEGER` | all | Convert a specific region number |
| `--resample / --no-resample` | enabled | Mass axis resampling |
| `--include-optical / --no-optical` | enabled | Include optical images in output |

### Examples

```bash
# Basic conversion -- format, pixel size, and resampling all auto-detected
thyra input.imzML output.zarr

# Specify pixel size manually (when metadata is unavailable)
thyra input.imzML output.zarr --pixel-size 25

# Convert only region 0 from a multi-region dataset
thyra data.d output.zarr --region 0

# Skip optical images
thyra data.d output.zarr --no-optical
```

!!! note "Region numbers"
    Region numbers start at 0. Use `-v DEBUG` to see which regions were detected
    and how many spectra each contains.

---

## Logging

| Option | Default | Description |
|--------|---------|-------------|
| `-v, --log-level LEVEL` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL` |
| `--log-file PATH` | none | Write logs to file |

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

## Resampling (Advanced)

These options control how spectra are mapped onto a common mass axis. In most
cases the defaults work well -- Thyra auto-detects the instrument type and
chooses an appropriate method and bin count.

| Option | Default | Description |
|--------|---------|-------------|
| `--resample-method METHOD` | `auto` | `auto`, `nearest_neighbor`, or `tic_preserving` |
| `--mass-axis-type TYPE` | `auto` | `auto`, `constant`, `linear_tof`, `reflector_tof`, `orbitrap`, `fticr` |
| `--resample-bins INTEGER` | auto | Number of bins (mutually exclusive with `--resample-width-at-mz`) |
| `--resample-min-mz FLOAT` | auto | Minimum m/z value |
| `--resample-max-mz FLOAT` | auto | Maximum m/z value |
| `--resample-width-at-mz FLOAT` | auto | Mass width in Da at reference m/z for physics-based binning |
| `--resample-reference-mz FLOAT` | `1000.0` | Reference m/z for width specification |
| `--resample-gap-tolerance FLOAT` | none | `tic_preserving` only: discard target bins farther than this many Da from any measured m/z, instead of interpolating across the gap |

!!! info "Choosing a resampling method"
    - **`nearest_neighbor`** -- Each target bin takes the nearest original m/z
      value. Correct for **centroid** data, where peaks are discrete masses.
    - **`tic_preserving`** -- Linear interpolation, rescaled so the total ion
      current is unchanged. Correct for **profile** data on a target axis
      whose bin widths scale the same way the source points are spaced --
      pair it only with `constant` unless you know otherwise.
    - **`auto`** -- Picks `tic_preserving` only for Bruker flexImaging /
      Rapiflex data, whose source grid is uniform in m/z, and
      `nearest_neighbor` for everything else -- including profile data from
      an instrument Thyra cannot identify, because the interpolating method
      is only exact when the two axis laws match. See
      [Resampling](resampling.md#which-detector-wins) for the full decision
      table.

!!! info "Choosing a mass axis type"
    The axis type determines how bin widths scale with m/z:

    - **`constant`** -- Uniform bin width (Da). Suitable for MALDI-TOF in linear mode.
    - **`linear_tof`** -- Width scales as sqrt(m/z). Matches TOF resolution.
    - **`reflector_tof`** -- Width scales linearly with m/z (constant relative resolution). Matches reflector TOF.
    - **`orbitrap`** -- Width scales as m/z^(3/2). Matches Orbitrap resolution.
    - **`fticr`** -- Width scales as m/z^2. Matches FTICR resolution.
    - **`auto`** -- Detected from instrument metadata.

### Examples

```bash
# Physics-based resampling for Orbitrap data
thyra input.imzML output.zarr \
    --resample-method tic_preserving \
    --mass-axis-type orbitrap

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

| Option | Default | Description |
|--------|---------|-------------|
| `--streaming [auto\|true\|false]` | `auto` | Streaming mode for large datasets |
| `--sparse-format [csc\|csr]` | `csc` | Sparse matrix storage format |

!!! info "Streaming mode"
    - **`auto`** (default) -- Thyra estimates dataset size and enables streaming
      for datasets over ~10 GB.
    - **`true`** -- Force streaming. Useful if auto-detection underestimates.
    - **`false`** -- Force standard (in-memory) conversion.

    Streaming processes spectra in chunks and writes incrementally to disk. The
    output is identical to standard mode.

### Examples

```bash
# Force streaming for a large dataset
thyra large.d output.zarr --streaming true

# Use CSR format (faster row access, slower column access)
thyra input.imzML output.zarr --sparse-format csr
```

!!! tip "CSC vs CSR"
    **CSC** (default) is optimised for extracting ion images (one m/z across all
    pixels). **CSR** is optimised for extracting spectra (one pixel across all
    m/z values). Choose based on your downstream access pattern.

!!! warning "`--optimize-chunks` is deprecated"
    The flag is still accepted, so existing scripts keep running, but it does
    nothing and now logs a warning. It will be dropped in a future release.

    It never did anything: the post-hoc pass it invoked was written for a dense
    4-D image layout and could not read the sparse table the converter actually
    writes, so it failed on every conversion and the CLI still exited 0. Chunk
    sizes are chosen at write time instead.

---

## Bruker-Specific

These options only apply when converting Bruker `.d` directories.

| Option | Default | Description |
|--------|---------|-------------|
| `--use-recalibrated / --no-recalibrated` | enabled | Use recalibrated m/z state |
| `--interactive-calibration` | off | Display available calibration states |
| `--intensity-threshold FLOAT` | none | Minimum intensity filter |

### Examples

```bash
# Use raw (non-recalibrated) m/z values
thyra data.d output.zarr --no-recalibrated

# Interactively choose calibration state
thyra data.d output.zarr --interactive-calibration

# Filter low-intensity signals (useful for continuous-mode Bruker data)
thyra data.d output.zarr --intensity-threshold 100
```

!!! warning "Intensity threshold"
    The `--intensity-threshold` option drops all peaks below the given value
    **before** writing to zarr. This reduces file size but is irreversible.
    Use with care -- inspect the data with `-v DEBUG` first to choose an
    appropriate threshold.

---

## Other

| Option | Default | Description |
|--------|---------|-------------|
| `--dataset-id TEXT` | `msi_dataset` | Dataset identifier used in element keys |
| `--handle-3d` | off | Process as 3D volume instead of 2D slices |

### Examples

```bash
# Custom dataset ID (affects table and image key names)
thyra input.imzML output.zarr --dataset-id hippocampus
# -> table key: hippocampus_z0, TIC key: hippocampus_z0_tic

# Combine z-slices into a single 3D table
thyra volume.imzML output.zarr --handle-3d
```
