# Thyra

[![PyPI](https://img.shields.io/pypi/v/thyra?logo=pypi&logoColor=white)](https://pypi.org/project/thyra/)
[![Tests](https://img.shields.io/github/actions/workflow/status/M4i-Imaging-Mass-Spectrometry/thyra/tests.yml?branch=main&logo=github)](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/actions/workflows/tests.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Thyra** (from Greek *thyra*, meaning "door" or "portal") converts Mass Spectrometry Imaging (MSI) data into the standardized [SpatialData](https://spatialdata.scverse.org/)/Zarr format -- your portal to spatial omics analysis workflows.

---

## Why Thyra?

Mass spectrometry imaging produces rich spatial-molecular data, but every vendor stores it differently. Downstream tools -- napari, squidpy, scanpy -- expect a common format. Thyra bridges that gap:

```
 .imzML  ──┐                        ┌── napari visualisation
 .d      ──┼──  thyra  ──> .zarr  ──┼── squidpy / scanpy analysis
 .raw    ──┘   (SpatialData)        └── custom Python workflows
```

The output is a single SpatialData/Zarr directory containing intensity matrices, TIC images, optical images, pixel geometries, and full metadata -- ready for any tool in the scverse ecosystem.

---

## Features

| | Feature | Description |
|---|---------|-------------|
| **Formats** | Multiple inputs | ImzML, Bruker (.d timsTOF + Rapiflex), Waters (.raw) |
| **Output** | SpatialData/Zarr | Cloud-ready, chunked, standardised |
| **Scale** | Memory efficient | Streaming mode for 100+ GB datasets |
| **Optics** | Optical alignment | Automatic MSI-to-microscopy registration (Bruker) |
| **Regions** | Multi-region | Handles slides with multiple tissue sections |
| **Resampling** | Physics-aware | Instrument-specific mass axis resampling (on by default) |
| **3D** | Volume support | Process as 3D volume or separate 2D slices |
| **Platform** | Cross-platform | Windows, macOS, Linux |

---

## Quick Start

### Install

```bash
pip install thyra
```

### Convert

=== "CLI"

    ```bash
    thyra input.imzML output.zarr
    ```

=== "Python"

    ```python
    from thyra import convert_msi

    success = convert_msi("input.imzML", "output.zarr")
    ```

### Explore the output

```python
import spatialdata as sd

sdata = sd.read_zarr("output.zarr")

# Intensity matrix (pixels x m/z bins)
table = sdata.tables["msi_dataset_z0"]
print(f"Shape: {table.shape}")
print(f"m/z range: {table.var['mz'].min():.1f} -- {table.var['mz'].max():.1f}")

# TIC image
import numpy as np
tic = np.asarray(sdata.images["msi_dataset_z0_tic"])[0]
```

!!! tip "What is in the output?"
    See [Output Format](output-format.md) for the full structure: tables, TIC images, optical images, pixel shapes, regions, and metadata.

---

## Supported Formats

### Input

| Format | Extension | Instruments |
|--------|-----------|-------------|
| ImzML  | `.imzML`  | Any vendor exporting to open standard |
| Bruker | `.d`      | timsTOF fleX, Rapiflex MALDI-TOF |
| Waters | `.raw`    | MassLynx imaging (DESI, MALDI) |

### Output

| Format | Description |
|--------|-------------|
| **SpatialData/Zarr** | The [scverse](https://scverse.org/) standard for spatial omics -- cloud-ready, chunked, with coordinate transforms |

---

## Next Steps

- **[Getting Started](getting-started.md)** -- installation, first conversion, common workflows
- **[CLI Reference](cli.md)** -- every command-line option explained
- **[Output Format](output-format.md)** -- what the .zarr contains and how to use it
- **[API Reference](api.md)** -- Python API documentation
