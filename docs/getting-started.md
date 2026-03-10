# Getting Started

## Installation

=== "pip (recommended)"

    ```bash
    pip install thyra
    ```

=== "From source"

    ```bash
    git clone https://github.com/M4i-Imaging-Mass-Spectrometry/thyra.git
    cd thyra
    poetry install
    ```

!!! note "Requirements"
    Python 3.11 or 3.12. Windows, macOS, and Linux are all supported.
    Bruker readers require the vendor SDK DLLs (bundled for Windows).

---

## Your First Conversion

### Command Line

```bash
# ImzML file
thyra input.imzML output.zarr

# Bruker .d folder
thyra data.d output.zarr
```

That's it -- Thyra auto-detects the format, reads the pixel size from metadata,
resamples onto a common mass axis, and writes a SpatialData/Zarr directory.

### Python API

```python
from thyra import convert_msi

success = convert_msi("data/sample.imzML", "output/sample.zarr")
```

You can pass explicit parameters when needed:

```python
success = convert_msi(
    "data/experiment.d",
    "output/experiment.zarr",
    dataset_id="hippocampus",
    pixel_size_um=10.0,
)
```

!!! tip "Checking the result"
    After conversion, load the output and inspect:
    ```python
    import spatialdata as sd

    sdata = sd.read_zarr("output/sample.zarr")
    print(list(sdata.tables.keys()))
    print(list(sdata.images.keys()))
    ```
    See [Output Format](output-format.md) for the full structure.

---

## Multi-Dataset Folders

When you point Thyra at a folder containing multiple `.d` datasets (e.g., a slide
with several tissue sections), it prompts you to choose:

```
$ thyra "slide_folder/" output.zarr

Found 3 datasets in slide_folder:
  [1] C1501.d
  [2] E2501.d
  [3] E2506.d

Select dataset to convert: 3
  -> E2506.d
```

Each `.d` dataset is matched to its own `.mis` file for correct optical alignment,
so the TIC overlay lands on the right tissue section in the shared optical image.

!!! note "Converting a specific region"
    If a single dataset contains multiple acquisition regions (e.g., tissue +
    matrix), you can select one with `--region`:
    ```bash
    thyra data.d output.zarr --region 0
    ```

---

## Resampling

Mass axis resampling is **enabled by default**. This maps all spectra onto a
common m/z axis, which is required for most downstream analysis tools.

```bash
# Default -- resampling is on, method and binning auto-detected
thyra input.imzML output.zarr

# Disable resampling (keeps original m/z values per spectrum)
thyra input.imzML output.zarr --no-resample
```

For advanced control you can specify the method and instrument type:

```bash
# Physics-based resampling for Orbitrap data
thyra input.imzML output.zarr \
    --resample-method tic_preserving \
    --mass-axis-type orbitrap
```

!!! info "When to disable resampling"
    Disable resampling (`--no-resample`) if you need the raw, unmodified spectra
    -- for example, when doing your own peak picking or centroiding downstream.
    Note that without resampling, each pixel may have a different m/z axis.

See the [CLI Reference](cli.md#resampling-advanced) for all resampling options.

---

## Optical Images

For Bruker data, Thyra automatically includes optical (microscopy) images and
aligns them to the MSI data using teaching point calibration from the `.mis` file.

```bash
# Optical images included by default
thyra data.d output.zarr

# Skip optical images
thyra data.d output.zarr --no-optical
```

The TIC image is stored with an affine transform that maps it into the optical
image's coordinate space, so overlays work out of the box.

---

## 3D Data

By default, Thyra treats each z-slice as a separate 2D dataset with `_z{i}` suffixes:

```bash
# Default: separate 2D slices (msi_dataset_z0, msi_dataset_z1, ...)
thyra volume.imzML output.zarr

# Combined 3D volume (single table with x, y, z coordinates)
thyra volume.imzML output.zarr --handle-3d
```

---

## Streaming Mode

Thyra automatically switches to streaming mode for large datasets (estimated >10 GB).
You can also force it:

```bash
# Force streaming on
thyra large_dataset.d output.zarr --streaming true

# Force streaming off
thyra large_dataset.d output.zarr --streaming false
```

!!! info "How streaming works"
    Streaming mode processes spectra in chunks and writes incrementally to Zarr
    on disk, keeping memory usage roughly constant regardless of dataset size.
    The output is identical -- only the processing strategy changes.

---

## Troubleshooting

### "WinError 5: Access is denied"

This means the output `.zarr` directory is locked by another process (e.g., a
Python session, napari, or a Jupyter notebook that loaded it).

**Fix:** Close any program that has the zarr open, or write to a different output
path.

### "No module named 'timsdata'" or Bruker SDK errors

The Bruker SDK DLLs are bundled for Windows. On Linux/macOS, Bruker data requires
the vendor's `libtimsdata.so` / `libtimsdata.dylib` to be installed separately.

### Pixel size not detected

If Thyra cannot find pixel size metadata, it will ask for a manual value. You can
also pass it explicitly:

```bash
thyra input.imzML output.zarr --pixel-size 50
```

### Memory errors on large datasets

Force streaming mode:

```bash
thyra large.d output.zarr --streaming true
```

Or reduce the number of resampling bins:

```bash
thyra large.d output.zarr --resample-bins 20000
```

### Verbose logging for debugging

```bash
thyra input.imzML output.zarr -v DEBUG --log-file conversion.log
```

---

## What Next?

- **[CLI Reference](cli.md)** -- all command-line options
- **[Output Format](output-format.md)** -- understanding the zarr output structure
- **[API Reference](api.md)** -- Python API documentation
