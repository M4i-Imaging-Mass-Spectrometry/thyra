<p align="center">
  <img src="docs/assets/thyra-logotype.svg" alt="Thyra" width="420">
</p>

[![Tests](https://img.shields.io/github/actions/workflow/status/M4i-Imaging-Mass-Spectrometry/thyra/tests.yml?branch=main&logo=github)](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/thyra?logo=pypi&logoColor=white)](https://pypi.org/project/thyra/)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/M4i-Imaging-Mass-Spectrometry/thyra/blob/main/notebooks/Thyra_Validation_Workflow.ipynb)

**Thyra converts imaging mass spectrometry data into one open format.** It
reads the files your instrument writes (Bruker, Waters, PHI, or imzML from any
vendor) and saves them in [SpatialData](https://spatialdata.scverse.org/),
which analysis and viewing tools can read. You do not need to know how to
program.

**[Documentation](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra)** | [Install](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/install/) | [Your first conversion](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/getting-started/) | [Tutorial](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/tutorial/) | [Technical reference](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/technical-reference/)

## Quick start

```bash
uv tool install --python 3.13 thyra                     # or: pip install thyra (Python 3.12 or 3.13)
thyra-example-data example_data/synthetic_brain.imzML   # makes a small example dataset
thyra example_data/synthetic_brain.imzML example_data/synthetic_brain.zarr
```

New to the terminal or to uv? The
[Install](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/install/) page
walks you through every step. To try Thyra with nothing installed, open the
[notebook in Google Colab](https://colab.research.google.com/github/M4i-Imaging-Mass-Spectrometry/thyra/blob/main/notebooks/Thyra_Validation_Workflow.ipynb).
The [Tutorial](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/tutorial/)
also converts the published example dataset
([10.5281/zenodo.18326569](https://doi.org/10.5281/zenodo.18326569)).

## What you can convert

| Your instrument | Point Thyra at |
|---|---|
| Bruker timsTOF or solariX | the folder ending in `.d` |
| Bruker rapifleX | the acquisition folder |
| Waters (MassLynx) | the folder ending in `.raw` |
| PHI nanoTOF (ToF-SIMS) | the file ending in `.raw` |
| Any instrument that exports imzML | the `.imzML` file, with its `.ibd` file next to it |
| mzPeak (experimental) | the `.mzpeak` file |

Thyra reads timsTOF and Waters data on Windows and Linux; everything else
also works on macOS. Shimadzu support is in development: until then, export
imzML from IMAGEREVEAL MS.

PHI mosaic, MS/MS and depth-profiling acquisitions have so far been tested
only on synthetic files. If you have real data in those modes, please get in
touch.

## What you get

One folder ending in `.zarr`, holding:

- the spectrum of every pixel, on one shared m/z axis
- an image of the total ion current (TIC)
- the optical image of the slide, when the acquisition has one
- ion mobility and MS/MS data from Bruker timsTOF, when the acquisition has them
- a validated description of the acquisition, ready for METASPACE

Thyra handles datasets larger than your computer's memory, slides with
several regions, and 3D data.

## From Python

```python
from thyra import convert_msi
import spatialdata as sd

convert_msi("data/sample.imzML", "output/sample.zarr")

sdata = sd.read_zarr("output/sample.zarr")
table = sdata.tables["msi_dataset_z0"]  # one row per pixel, one column per m/z bin
```

The acquisition description has its own commands:

```bash
thyra metadata raw_data.d                                # describe a dataset without converting it
thyra validate output.zarr                               # check a result against the schema
thyra export-metaspace output.zarr --merge sample.json   # METASPACE submission file
```

Every option, format and output detail is in the
[Technical reference](https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/technical-reference/).

## Development

```bash
git clone https://github.com/M4i-Imaging-Mass-Spectrometry/thyra.git
cd thyra
uv sync
uv run pre-commit install
uv run pytest
```

CI runs the same hooks on every pull request. Installing them locally shows
you a failure before you push.

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

### Visual identity

Logomark and logotype designed by **Nepsis Scriptorium**.

[![Instagram @nepsis.scriptorium](https://img.shields.io/badge/Instagram-%40nepsis.scriptorium-E4405F?logo=instagram&logoColor=white)](https://www.instagram.com/nepsis.scriptorium/)
[![Email nepsisscriptorium@gmail.com](https://img.shields.io/badge/Email-nepsisscriptorium%40gmail.com-EA4335?logo=gmail&logoColor=white)](mailto:nepsisscriptorium@gmail.com)
