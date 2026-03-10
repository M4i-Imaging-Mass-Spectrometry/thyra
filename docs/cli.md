# CLI Reference

```
thyra [OPTIONS] INPUT OUTPUT
```

**INPUT**: Path to input MSI file or directory (`.imzML`, `.d`, `.raw`)

**OUTPUT**: Path for output `.zarr` file

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
# Basic conversion with defaults
thyra input.imzML output.zarr

# Specify pixel size manually
thyra input.imzML output.zarr --pixel-size 25

# Convert only region 2 from a multi-region dataset
thyra data.d output.zarr --region 2

# Skip optical images
thyra data.d output.zarr --no-optical
```

## Logging

| Option | Default | Description |
|--------|---------|-------------|
| `-v, --log-level LEVEL` | `INFO` | Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL |
| `--log-file PATH` | none | Write logs to file |

### Examples

```bash
# Verbose output
thyra input.imzML output.zarr -v DEBUG

# Save logs to file
thyra input.imzML output.zarr --log-file conversion.log
```

## Resampling (Advanced)

These options control how spectra are mapped onto a common mass axis. In most cases the defaults (auto-detection) work well.

| Option | Default | Description |
|--------|---------|-------------|
| `--resample-method METHOD` | auto | `auto`, `nearest_neighbor`, or `tic_preserving` |
| `--mass-axis-type TYPE` | auto | `auto`, `constant`, `linear_tof`, `reflector_tof`, `orbitrap`, `fticr` |
| `--resample-bins INTEGER` | auto | Number of bins (mutually exclusive with `--resample-width-at-mz`) |
| `--resample-min-mz FLOAT` | auto | Minimum m/z value |
| `--resample-max-mz FLOAT` | auto | Maximum m/z value |
| `--resample-width-at-mz FLOAT` | auto | Mass width in Da at reference m/z for physics-based binning |
| `--resample-reference-mz FLOAT` | 1000.0 | Reference m/z for width specification |

### Examples

```bash
# Physics-based resampling for Orbitrap data
thyra input.imzML output.zarr \
    --resample-method tic_preserving \
    --mass-axis-type orbitrap

# Fixed number of bins
thyra input.imzML output.zarr --resample-bins 50000

# Specify mass range
thyra input.imzML output.zarr --resample-min-mz 100 --resample-max-mz 1000
```

## Performance

| Option | Default | Description |
|--------|---------|-------------|
| `--streaming [auto\|true\|false]` | `auto` | Streaming mode for large datasets |
| `--optimize-chunks` | off | Optimize Zarr chunks after conversion |
| `--sparse-format [csc\|csr]` | `csc` | Sparse matrix storage format |

### Examples

```bash
# Stream a large dataset to keep memory usage low
thyra large.d output.zarr --streaming true

# Optimize chunk layout for downstream access patterns
thyra input.imzML output.zarr --optimize-chunks
```

## Bruker-Specific

| Option | Default | Description |
|--------|---------|-------------|
| `--use-recalibrated / --no-recalibrated` | enabled | Use recalibrated m/z state |
| `--interactive-calibration` | off | Display available calibration states |
| `--intensity-threshold FLOAT` | none | Minimum intensity filter (useful for continuous mode) |

### Examples

```bash
# Use raw (non-recalibrated) m/z values
thyra data.d output.zarr --no-recalibrated

# Interactively choose calibration state
thyra data.d output.zarr --interactive-calibration

# Filter low-intensity signals (continuous mode Bruker data)
thyra data.d output.zarr --intensity-threshold 100
```

## Other

| Option | Default | Description |
|--------|---------|-------------|
| `--dataset-id TEXT` | `msi_dataset` | Dataset identifier used in element keys |
| `--handle-3d` | off | Process as 3D volume instead of 2D slices |
