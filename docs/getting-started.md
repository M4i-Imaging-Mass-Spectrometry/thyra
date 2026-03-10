# Getting Started

## Installation

=== "pip"

    ```bash
    pip install thyra
    ```

=== "conda"

    ```bash
    conda install -c conda-forge thyra
    ```

=== "From source"

    ```bash
    git clone https://github.com/M4i-Imaging-Mass-Spectrometry/thyra.git
    cd thyra
    poetry install
    ```

## Basic Conversion

### Command Line

```bash
# ImzML file
thyra input.imzML output.zarr

# Bruker .d folder
thyra data.d output.zarr

# With custom pixel size (when auto-detection is not available)
thyra input.imzML output.zarr --pixel-size 50
```

### Python API

```python
from thyra import convert_msi

# Simple conversion (auto-detects pixel size from metadata)
success = convert_msi("data/sample.imzML", "output/sample.zarr")

# With explicit parameters
success = convert_msi(
    "data/experiment.d",
    "output/experiment.zarr",
    dataset_id="hippocampus",
    pixel_size_um=10.0,
)
```

## Multi-Dataset Folders

When you point Thyra at a folder containing multiple `.d` datasets (e.g., a slide with multiple tissue sections), it prompts you to choose:

```
$ thyra "slide_folder/" output.zarr

Found 3 datasets in slide_folder:
  [1] C1501.d
  [2] E2501.d
  [3] E2506.d

Select dataset to convert: 3
  -> E2506.d
```

Each `.d` dataset is matched to its own `.mis` file for correct optical alignment, so the TIC overlay lands on the right tissue section in the shared optical image.

## Resampling

Mass axis resampling is **enabled by default**. This maps all spectra onto a common m/z axis, which is required for most downstream analysis.

```bash
# Resampling is on by default -- just run:
thyra input.imzML output.zarr

# Disable resampling (keeps original m/z values per spectrum)
thyra input.imzML output.zarr --no-resample

# Specify resampling parameters
thyra input.imzML output.zarr \
    --resample-method tic_preserving \
    --mass-axis-type orbitrap
```

See the [CLI Reference](cli.md) for all resampling options.

## Optical Images

For Bruker data, Thyra automatically includes optical (microscopy) images alongside the MSI data. These are aligned using teaching point calibration from the `.mis` file.

```bash
# Optical images are included by default
thyra data.d output.zarr

# Disable optical images
thyra data.d output.zarr --no-optical
```

## 3D Data

By default, Thyra treats each z-slice as a separate 2D dataset. Use `--handle-3d` to combine slices into a single 3D volume:

```bash
# Default: separate 2D slices (msi_dataset_z0, msi_dataset_z1, ...)
thyra volume.imzML output.zarr

# Combined 3D volume (single msi_dataset table)
thyra volume.imzML output.zarr --handle-3d
```

## Streaming Mode

For datasets that exceed available memory (100+ GB), use streaming mode:

```bash
thyra large_dataset.d output.zarr --streaming true
```

Streaming mode processes data in chunks and writes incrementally to disk, keeping memory usage constant regardless of dataset size.

## What Next?

- [CLI Reference](cli.md) -- all command-line options
- [Output Format](output-format.md) -- understanding the zarr output structure
- [API Reference](api.md) -- Python API documentation
