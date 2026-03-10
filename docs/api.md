# API Reference

Thyra's Python API centres on a single function: `convert_msi`. For most use
cases, that is all you need. The reader and converter base classes are
documented below for advanced users who want to extend Thyra with new formats.

---

## Converting Data

The primary entry point. Detects the input format, reads metadata, and writes
a SpatialData/Zarr directory.

### Basic usage

```python
from thyra import convert_msi

# Minimal -- auto-detects format, pixel size, resampling, and streaming
success = convert_msi("input.imzML", "output.zarr")

# With explicit parameters
success = convert_msi(
    "data/experiment.d",
    "output/experiment.zarr",
    dataset_id="hippocampus",
    pixel_size_um=10.0,
)
```

### With resampling configuration

```python
success = convert_msi(
    "input.imzML",
    "output.zarr",
    resampling_config={
        "method": "tic_preserving",
        "axis_type": "orbitrap",
        "n_bins": 50000,
    },
)
```

### Multi-region dataset (select one region)

```python
success = convert_msi(
    "data/slide.d",
    "output/tissue_only.zarr",
    region=0,  # convert only region 0
)
```

### Force streaming for large datasets

```python
success = convert_msi(
    "data/large_dataset.d",
    "output/large.zarr",
    streaming=True,
)
```

### Full signature

::: thyra.convert.convert_msi

---

## Reader Base Class

All format readers (ImzML, Bruker, Waters) inherit from this base class. If
you are writing a custom reader for a new format, subclass `BaseMSIReader`
and implement the abstract methods below.

::: thyra.core.base_reader.BaseMSIReader
    options:
      members:
        - get_essential_metadata
        - get_comprehensive_metadata
        - get_common_mass_axis
        - get_optical_image_paths
        - iter_spectra
        - get_region_map
        - get_region_info
        - has_shared_mass_axis

---

## Converter Base Class

All output converters inherit from this base class. Currently only
SpatialData output is supported, but the architecture allows adding new output
formats by subclassing `BaseMSIConverter`.

::: thyra.core.base_converter.BaseMSIConverter
    options:
      members:
        - convert
        - pixel_size_um
        - pixel_size_source
        - dataset_id
        - handle_3d
