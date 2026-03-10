# Thyra

**Thyra** (from Greek thyra, meaning "door" or "portal") -- a modern Python library for converting Mass Spectrometry Imaging (MSI) data into the standardized SpatialData/Zarr format, serving as your portal to spatial omics analysis workflows.

## Features

- **Multiple Input Formats** -- ImzML, Bruker (.d directories), Waters (.raw directories)
- **SpatialData Output** -- Modern, cloud-ready format with Zarr backend
- **Memory Efficient** -- Handles large datasets (100+ GB) through streaming processing
- **Optical Alignment** -- Automatic MSI-to-optical image registration for Bruker data
- **Multi-Region Support** -- Handles slides with multiple tissue sections
- **Resampling** -- Physics-aware mass axis resampling (enabled by default)
- **3D Support** -- Process volume data or treat as 2D slices
- **Cross-Platform** -- Windows, macOS, and Linux

## Quick Install

```bash
pip install thyra
```

## Supported Formats

### Input

| Format | Extension | Description |
|--------|-----------|-------------|
| ImzML  | `.imzML`  | Open standard for MS imaging |
| Bruker | `.d`      | Bruker timsTOF and Rapiflex |
| Waters | `.raw`    | Waters MassLynx imaging |

### Output

| Format | Description |
|--------|-------------|
| SpatialData/Zarr | Modern spatial omics standard -- cloud-ready, efficient, standardized |

## Quick Example

=== "CLI"

    ```bash
    thyra input.imzML output.zarr
    ```

=== "Python"

    ```python
    from thyra import convert_msi

    convert_msi("input.imzML", "output.zarr")
    ```

See the [Getting Started](getting-started.md) guide for more details.
