# CLI Reference

```
thyra [OPTIONS] INPUT OUTPUT
```

**INPUT** -- Path to input MSI file or directory (`.imzML`, `.d`, `.raw`)

**OUTPUT** -- Path for output `.zarr` directory

!!! tip "Grouped help"
    Run `thyra --help` to see all options organised by category (Conversion,
    Logging, Resampling, Performance, Bruker-Specific, Other).

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

!!! info "Choosing a resampling method"
    - **`nearest_neighbor`** -- Fast, simple assignment to nearest bin. Good for
      data that is already close to uniformly spaced.
    - **`tic_preserving`** -- Distributes intensity proportionally across bins.
      Better for high-resolution data (Orbitrap, FTICR) where bin widths vary.
    - **`auto`** -- Picks `tic_preserving` for high-resolution instruments,
      `nearest_neighbor` otherwise.

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
| `--optimize-chunks` | off | Optimise Zarr chunks after conversion |
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

# Optimise chunk layout for downstream column-access patterns
thyra input.imzML output.zarr --optimize-chunks

# Use CSR format (faster row access, slower column access)
thyra input.imzML output.zarr --sparse-format csr
```

!!! tip "CSC vs CSR"
    **CSC** (default) is optimised for extracting ion images (one m/z across all
    pixels). **CSR** is optimised for extracting spectra (one pixel across all
    m/z values). Choose based on your downstream access pattern.

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
