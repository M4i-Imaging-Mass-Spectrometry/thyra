# Your first conversion

Converting a dataset takes one command: `thyra`, then your data, then a name
for the result.

```bash
thyra your_data result.zarr
```

This page shows each part. You need Thyra installed first: see
[Install](install.md).

## Step 1: Find your data

Point Thyra at the file or folder your instrument wrote:

| Your instrument | Point Thyra at |
|---|---|
| Bruker timsTOF or solariX | the folder ending in `.d` |
| Bruker rapifleX | the folder that holds the `.dat` and `_poslog.txt` files |
| Waters (MassLynx) | the folder ending in `.raw` |
| PHI nanoTOF (ToF-SIMS) | the file ending in `.raw` |
| Any instrument that exports imzML | the `.imzML` file, with its `.ibd` file next to it |
| mzPeak (experimental) | the `.mzpeak` file |

No data at hand? This command makes a small example dataset to practise on:

```bash
thyra-example-data example_data/synthetic_brain.imzML
```

## Step 2: Convert

Type `thyra`, then the path to your data, then a name for the result ending
in `.zarr`:

```bash
thyra "D:\MSI\brain section 1.d" brain_section_1.zarr
```

Three tips:

- **Put quotes around a path that contains spaces**, as above.
- **Drag and drop.** Drag your data from the file browser onto the terminal
  window, and its full path appears. This saves typing it.
- **The result goes to the folder your terminal is in**, which is shown
  before the cursor. To save it somewhere else, give a full path, such as
  `"D:\MSI\brain_section_1.zarr"`.

Thyra prints its progress as it works. When it is done, the last line says:

```
Conversion completed successfully. Output stored at brain_section_1.zarr
```

Small datasets take seconds. Large ones can take several minutes.

## Step 3: Check the result

A new folder, `brain_section_1.zarr`, is your converted dataset. Keep the
folder whole: copy or move it like any other folder, but do not take files
out of it.

To look inside it, see [Look at the result](explore-the-output.ipynb). Any
tool that reads SpatialData can open it, including napari with the
[napari-spatialdata](https://spatialdata.scverse.org/projects/napari/en/latest/)
plugin.

## Common situations

Most datasets need nothing more than the command above. These are the
exceptions you are most likely to meet.

### Your folder holds several datasets

Point Thyra at a folder that holds several `.d` datasets, and it lists them
and asks which one to convert:

```
Found 3 datasets in slide_folder:
  [1] C1501.d
  [2] E2501.d
  [3] E2506.d

Select dataset to convert: 3
```

Type the number and press Enter.

### Your acquisition has several regions

Thyra converts all regions by default. To convert only one, add `--region`
and the region's name as flexImaging shows it:

```bash
thyra data.d region_03.zarr --region 03
```

Thyra lists the regions it finds when it starts, so a first run shows you
the names.

### Thyra cannot find the pixel size

Most files record their pixel size. If yours does not, Thyra stops before
it writes anything and tells you so. Run it again with `--pixel-size` and
the size in micrometres:

```bash
thyra input.imzML output.zarr --pixel-size 50
```

### You do not want the optical image

Bruker acquisitions often include a photo of the slide, and Thyra adds it to
the result. To leave it out, add `--no-optical`.

### You want the original m/z values

Thyra puts every spectrum on one shared m/z axis, because most analysis tools
need that. To keep each spectrum's own m/z values instead, add
`--no-resample`.

??? advanced "Advanced: How Thyra chooses the shared m/z axis"
    This step is called resampling. Thyra looks at the instrument type and
    the data to choose the spacing of the axis and the method that moves each
    peak onto it. You can set the method, the spacing and the mass range
    yourself. [Resampling](resampling.md) explains each choice and when to
    change it.

### You have a 3D dataset

Thyra saves each slice as its own 2D dataset. To keep them together as one
3D volume, add `--handle-3d`.

## Convert from Python

??? advanced "Advanced: Use Thyra inside a Python script or notebook"
    The same conversion is one function call:

    ```python
    from thyra import convert_msi

    convert_msi("brain_section_1.d", "brain_section_1.zarr")
    ```

    It takes the same options as the command, such as `pixel_size_um=50`.
    One difference: from Python, resampling is off unless you turn it on.
    The [Python API](api.md) lists every argument.

## If something goes wrong

[Troubleshooting](troubleshooting.md) lists the common errors and how to fix
each one.

## Go deeper

- [Command-line options](cli.md): every option, with examples.
- [Resampling](resampling.md): how the shared m/z axis is chosen.
- [Output format](output-format.md): everything inside the `.zarr` folder.
