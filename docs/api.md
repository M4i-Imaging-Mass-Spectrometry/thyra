# API Reference

Thyra's Python API centres on a single function: `convert_msi`. For most use
cases, that is all you need. The remaining sections document configuration
types, metadata objects, and base classes for advanced users who want to inspect
results or extend Thyra with new formats.

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
        "target_bins": 50000,
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

## Resampling Configuration

When you pass `resampling_config` to `convert_msi`, the dictionary keys map
to the fields of `ResamplingConfig`. You can pass a plain dict (as shown in
the examples above) or construct the dataclass directly:

```python
from thyra.resampling.types import ResamplingConfig, ResamplingMethod, AxisType

config = ResamplingConfig(
    method=ResamplingMethod.TIC_PRESERVING,
    axis_type=AxisType.ORBITRAP,
    target_bins=50000,
)

success = convert_msi("input.imzML", "output.zarr", resampling_config=config)
```

::: thyra.resampling.types.ResamplingConfig
    options:
      show_root_heading: true
      heading_level: 3

::: thyra.resampling.types.ResamplingMethod
    options:
      show_root_heading: true
      heading_level: 3
      members: true

::: thyra.resampling.types.AxisType
    options:
      show_root_heading: true
      heading_level: 3
      members: true

---

## Metadata Types

Readers expose metadata through two dataclasses. `EssentialMetadata` contains
everything needed for conversion decisions (grid size, mass range, memory
estimate). `ComprehensiveMetadata` wraps essential metadata and adds
vendor-specific details for provenance and QC.

```python
from thyra.readers.imzml import ImzMLReader

with ImzMLReader("sample.imzML") as reader:
    meta = reader.get_essential_metadata()
    print(f"Grid: {meta.dimensions}")
    print(f"m/z range: {meta.mass_range}")
    print(f"Spectra: {meta.n_spectra}")
    print(f"Est. memory: {meta.estimated_memory_gb:.1f} GB")
```

::: thyra.metadata.types.EssentialMetadata
    options:
      show_root_heading: true
      heading_level: 3

::: thyra.metadata.types.ComprehensiveMetadata
    options:
      show_root_heading: true
      heading_level: 3

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
        - close

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

---

## Format Detection and Plugin Registry

Thyra uses a registry to map file extensions and directory structures to the
correct reader and converter classes. The public functions below let you detect
formats programmatically or register your own reader/converter.

### Detecting a format

```python
from pathlib import Path
from thyra.core.registry import detect_format

fmt = detect_format(Path("experiment.imzML"))  # "imzml"
fmt = detect_format(Path("data.d"))            # "bruker" or "rapiflex"
fmt = detect_format(Path("data.raw"))          # "waters"
```

### Registering a custom reader

```python
from thyra.core.registry import register_reader
from thyra.core.base_reader import BaseMSIReader

@register_reader("my_format")
class MyFormatReader(BaseMSIReader):
    ...
```

::: thyra.core.registry.detect_format
    options:
      show_root_heading: true
      heading_level: 3

::: thyra.core.registry.register_reader
    options:
      show_root_heading: true
      heading_level: 3

::: thyra.core.registry.register_converter
    options:
      show_root_heading: true
      heading_level: 3
