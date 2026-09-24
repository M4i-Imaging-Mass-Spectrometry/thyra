# Technical reference

These pages hold every detail: every option, every format, and the exact
layout of the result. You do not need them to convert a dataset. For that,
start with [Install](install.md) and [Your first conversion](getting-started.md).

| Page | What it covers |
|---|---|
| [Supported formats](supported-formats.md) | How each format is recognised, what it provides, and what is special about each vendor |
| [Command-line options](cli.md) | Every option of the `thyra` command, with examples |
| [Resampling](resampling.md) | How the shared m/z axis is chosen, and how to control it |
| [Output format](output-format.md) | Everything inside the `.zarr` folder, and how to read it in Python |
| [Coordinate systems](coordinate-systems.md) | How positions in the result relate to the slide and the optical image |
| [Metadata schema](metadata-schema.md) | The structured description of each acquisition, and METASPACE export |
| [Writing the metadata document](writing-the-metadata-document.md) | For other programs that want to write the same description |
| [Python API](api.md) | Every function and class you can use from Python |
| [Design decisions](design-decisions.md) | Why each default is what it is, with the measurements behind it |
| Format notes | Background on reading [imzML](imzml-parser-notes.md), [PHI ToF-SIMS](phi-tofsims-notes.md) and [Bruker solariX](solarix-notes.md) files |
