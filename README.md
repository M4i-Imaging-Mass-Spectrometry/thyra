# Thyra

[![Tests](https://img.shields.io/github/actions/workflow/status/M4i-Imaging-Mass-Spectrometry/thyra/tests.yml?branch=main&logo=github)](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/thyra?logo=pypi&logoColor=white)](https://pypi.org/project/thyra/)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra)

**Thyra** (from Greek thyra, meaning "door" or "portal") -- a modern Python library for converting Mass Spectrometry Imaging (MSI) data into the standardized **SpatialData/Zarr format**, serving as your portal to spatial omics analysis workflows.

## Features

- **Multiple Input Formats**: ImzML, Bruker (.d directories), Waters (.raw directories)
- **SpatialData Output**: Modern, cloud-ready format with Zarr backend
- **Memory Efficient**: Handles large datasets (100+ GB) through streaming processing
- **Optical Alignment**: Automatic MSI-to-optical image registration for Bruker data
- **Multi-Region Support**: Handles slides with multiple tissue sections
- **Resampling**: Physics-aware mass axis resampling (enabled by default)
- **3D Support**: Process volume data or treat as 2D slices
- **Cross-Platform**: Windows, macOS, and Linux

## Installation

```bash
pip install thyra
```

## Quick Start

### Command Line

```bash
# Basic conversion (resampling enabled by default)
thyra input.imzML output.zarr

# Bruker data with verbose logging
thyra data.d output.zarr -v DEBUG

# Disable resampling
thyra input.imzML output.zarr --no-resample
```

### Python API

```python
from thyra import convert_msi

success = convert_msi("data/sample.imzML", "output/sample.zarr")
```

### Working with the Output

```python
import spatialdata as sd

sdata = sd.read_zarr("output/sample.zarr")
msi_table = sdata.tables["msi_dataset_z0"]

print(f"Shape: {msi_table.shape}")  # (pixels, m/z bins)
print(f"m/z range: {msi_table.var['mz'].min():.1f} -- {msi_table.var['mz'].max():.1f}")
```

## Documentation

Full documentation: **[M4i-Imaging-Mass-Spectrometry.github.io/thyra](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra)**

- [Getting Started](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/getting-started/) -- installation, first conversion, common workflows
- [CLI Reference](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/cli/) -- all command-line options
- [Output Format](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/output-format/) -- understanding the zarr structure
- [API Reference](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/api/) -- Python API documentation

## Supported Formats

| Input | Extension | Status |
|-------|-----------|--------|
| ImzML | `.imzML` | Full support |
| Bruker | `.d` | Full support (timsTOF + Rapiflex) |
| Waters | `.raw` | Full support |

Output: **SpatialData/Zarr** -- cloud-ready, efficient, standardized

## Development

```bash
git clone https://github.com/M4i-Imaging-Mass-Spectrometry/thyra.git
cd thyra
poetry install
poetry run pre-commit install
poetry run pytest
```

## Contributing

See [CONTRIBUTING.md](docs/contributing.md) for guidelines.

## License

MIT -- see [LICENSE](LICENSE).

## Citation

```bibtex
@software{thyra2024,
  title = {Thyra: Modern Mass Spectrometry Imaging Data Conversion},
  author = {Visvikis, Theodoros},
  year = {2024},
  url = {https://github.com/M4i-Imaging-Mass-Spectrometry/thyra}
}
```

## Acknowledgments

- Built with [SpatialData](https://spatialdata.scverse.org/) ecosystem
- Powered by [Zarr](https://zarr.readthedocs.io/) for efficient storage
- Uses [pyimzML](https://github.com/alexandrovteam/pyimzML) for ImzML parsing
