# Supported Formats

Thyra reads five MSI formats and writes all of them into the same
SpatialData/Zarr layout. The input format is detected from the path -- there is
no format flag on the CLI, and `format_type` in the Python API selects the
*output* format, not the input.

| Format | Path shape | Detected by | Vendor SDK |
|---|---|---|---|
| **imzML** | `.imzML` file + `.ibd` | extension, `.ibd` must exist | none |
| **Bruker timsTOF** | `.d` directory | `analysis.tsf` or `analysis.tdf` | bundled DLL |
| **Bruker Rapiflex** | directory | `*.dat` + `*_poslog.txt` | none |
| **Waters MassLynx** | `.raw` **directory** | `_FUNC*.DAT` files inside | bundled DLL |
| **PHI SmartSoft-TOF** | `.raw` **file** | `SOFH` magic in first 4 bytes | none |

```bash
thyra sample.imzML       out.zarr   # imzML
thyra sample.d           out.zarr   # Bruker timsTOF
thyra rapiflex_folder/   out.zarr   # Bruker Rapiflex
thyra waters_run.raw/    out.zarr   # Waters (a directory)
thyra tofsims_run.raw    out.zarr   # PHI (a file)
```

---

## The `.raw` collision

Two vendors claim `.raw`, and they are told apart by **shape, not extension**:

- **Waters** `.raw` is a *directory* containing `_FUNC001.DAT`, `_FUNC002.DAT`, …
- **PHI** `.raw` is a *single file* whose header begins with the ASCII magic
  `SOFH`

Detection checks the directory case first, then the file magic. A `.raw` path
that is neither raises an error naming both possibilities, rather than a
confusing Waters-specific complaint:

```
Unrecognised .raw file: <path>. Expected either a Waters directory containing
_FUNC*.DAT files, or a PHI SmartSoft-TOF file beginning with the SOFH magic.
```

---

## What each format provides

Not every source carries every kind of metadata. This is what Thyra can
actually populate from each.

| | imzML | timsTOF | Rapiflex | Waters | PHI |
|---|---|---|---|---|---|
| Pixel size from metadata | yes | yes | yes | yes | yes |
| Optical image | -- | yes | yes | -- | -- |
| Optical alignment | -- | yes (`.mis`) | -- | -- | -- |
| Multi-region | -- | yes | -- | -- | mosaic tiles |
| 3D / multi-slice | yes | yes | -- | -- | -- |
| Native non-m/z axis kept | -- | -- | -- | -- | flight time |

Anything a format does not supply is simply absent from the output rather than
guessed at. Pixel size is the one exception worth knowing about: when a source
cannot report it, the CLI falls back to a default and records that it did so in
`uns` (see [Output Format](output-format.md)).

---

## imzML

The open interchange format, read through
[pyimzML](https://github.com/alexandrovteam/pyimzML). Both storage modes work:

- **continuous** -- every spectrum shares one m/z array, so the common mass
  axis is read once
- **processed** -- each spectrum carries its own m/z values, so building the
  common axis requires a pass over every spectrum

An `.imzML` without its `.ibd` beside it is rejected up front, because the XML
holds only offsets and the binary holds the data.

See [imzML Parser Notes](imzml-parser-notes.md) for the hazards Thyra works
around in the underlying library.

## Bruker timsTOF

`.d` directories containing `analysis.tsf` (TOF) or `analysis.tdf` (with ion
mobility). Reading goes through Bruker's `timsdata` library, which is bundled
for Windows and Linux. This is the richest source Thyra handles: it carries
optical microscopy images, FlexImaging `.mis` teaching points for MSI-to-optical
registration, and per-pixel region annotations for multi-region slides.

## Bruker Rapiflex

A folder of `*.dat` files with `*_poslog.txt` and `*_info.txt` alongside. The
position log supplies the pixel grid. No SDK required.

## Waters MassLynx

A `.raw` **directory** of `_FUNC*.DAT` files, read through the MassLynxRaw and
MLReader native libraries (bundled for Windows and Linux). The pixel grid is
reconstructed from the laser X/Y position recorded on each scan, and only MS
functions are converted -- lockmass, MRM and ion-mobility functions are
classified and skipped.

## PHI SmartSoft-TOF (ToF-SIMS)

A single `.raw` **file** from PHI (Physical Electronics) nanoTOF instruments.
Unlike Waters and Bruker this needs **no vendor SDK** -- the format is parsed
directly, so it behaves identically on every platform.

It differs from the other formats in one important way: it records **individual
ion arrivals**, not per-pixel spectra. Thyra aggregates those events into sparse
spectra on the detector's time-channel grid. Because flight time is what the
instrument actually measures, Thyra also stores it as `var["tof_us"]` so the
mass calibration stays reversible.

```bash
thyra tofsims_run.raw out.zarr
```

```python
from thyra.readers.phi import PhiReader

with PhiReader("tofsims_run.raw") as reader:
    print(reader.dimensions)          # (512, 512, 1)
    print(reader.pixel_size_um)       # 1.0
    print(reader.calibration_source)  # 'appended' or 'header'
```

See [PHI ToF-SIMS Notes](phi-tofsims-notes.md) for the file layout, the
calibration behaviour, and why the time axis is binned the way it is.

---

## Adding another format

A reader subclasses `BaseMSIReader` and implements four methods --
`_create_metadata_extractor`, `get_common_mass_axis`, `iter_spectra` and
`close` -- then registers itself with `@register_reader("name")`. The converter
is format-agnostic, so 2D/3D handling, resampling, chunking and the whole
SpatialData output path come for free.

Optional overrides that are worth implementing when the format allows it:

| Method | Buys you |
|---|---|
| `has_shared_mass_axis` | skips a full pass when the axis is fixed |
| `get_peak_counts_per_pixel` | single-pass CSR build instead of two passes |
| `get_mass_axis_annotations` | keeps a native non-m/z axis in `var` |
| `get_region_map` / `get_region_info` | per-pixel region annotation |
| `get_optical_image_paths` | optical images carried into the output |

See [Contributing](contributing.md).
