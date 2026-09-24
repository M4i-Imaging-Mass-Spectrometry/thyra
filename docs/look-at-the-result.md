# Look at the result

Your result is a folder ending in `.zarr`. This page shows what is inside it
and how to open it.

## What is inside

| Part | What it holds |
|---|---|
| **The spectra** | One row per measured pixel and one column per m/z bin |
| **Average spectrum** | The mean spectrum of all measured pixels |
| **TIC image** | The total ion current of every pixel, as an image |
| **Pixel outlines** | The position and size of every measured pixel |
| **Optical image** | The photo of the slide, when the acquisition has one (Bruker) |
| **Description** | The instrument, the settings and the pixel size: see [Describe and share your data](describe-your-data.md) |

There are no ion images in the result: you make them from the spectra, as
shown below. Most pixels have no signal in most m/z bins, so Thyra stores only
the values that are there, compressed. That keeps the result small.

??? advanced "Advanced: the names of the parts"
    Every part is named after the dataset, `msi_dataset` unless you choose
    another name with `--dataset-id`, and after the slice, `_z0` for the
    first one:

    - `msi_dataset_z0`: the table of spectra
    - `msi_dataset_z0_tic`: the TIC image
    - `msi_dataset_z0_pixels`: the pixel outlines
    - `msi_dataset_optical_...`: the optical images, named after their image
      files

    Some data adds more. Bruker timsTOF data with ion mobility carries a
    mass-mobility heatmap in the table, and PASEF MS/MS data adds a table
    ending in `_msms`. [Output format](output-format.md) lists every part,
    including those for 3D data.

## Keep the folder whole

The `.zarr` folder is one dataset. Copy it, move it or share it like any
other folder, but do not take files out of it or rename files inside it. It
opens on Windows, macOS and Linux alike.

## Open it without writing code: napari

[napari](https://napari.org) is a free image viewer. With the
napari-spatialdata plugin, it opens Thyra results.

1. Install napari with the plugin. With uv (see [Install](install.md)):

    ```bash
    uv tool install --python 3.13 napari --with "napari-spatialdata[all]"
    ```

2. Start napari:

    ```bash
    napari
    ```

3. Open your `.zarr` folder with **File > Open Folder...**. If napari asks
   which plugin should open it, choose **napari-spatialdata**.
4. A **SpatialData** panel appears on the left. Click **global** in its list
   of coordinate systems.
5. The panel now lists the parts of the result. Double-click
   `msi_dataset_z0_tic` to show the TIC image.

napari shows the TIC image, the optical image and the pixel outlines well.
For an image of one m/z value, use Python, as below.

## Open it in Python

A few lines of Python give you the TIC image, or an image of any m/z value.
The [Tutorial](tutorial.md) runs them step by step on example data.

You need Thyra in the Python you work in, and matplotlib to draw the images
(see [Already use Python?](install.md#already-use-python)):

```bash
pip install thyra matplotlib
```

Then:

```python
import matplotlib.pyplot as plt
import numpy as np
import spatialdata as sd

sdata = sd.read_zarr("brain_section_1.zarr")

# The TIC image
tic = np.asarray(sdata.images["msi_dataset_z0_tic"])[0]
plt.imshow(tic)
plt.show()

# An image of m/z 760.6: add up the intensities within 0.25 of it
table = sdata.tables["msi_dataset_z0"]
mz = table.var["mz"].values
lo, hi = np.searchsorted(mz, [760.35, 760.85])
values = np.asarray(table.X[:, lo:hi].sum(axis=1)).ravel()
x = table.obs["x"].values.astype(int)
y = table.obs["y"].values.astype(int)
ion_image = np.zeros((y.max() + 1, x.max() + 1))
ion_image[y, x] = values
plt.imshow(ion_image)
plt.show()
```

??? advanced "Advanced: why a window, and not one bin"
    On the shared m/z axis, most bins are empty. The single bin nearest your
    m/z value is usually one of them, so an image of it would be blank. Add
    up the bins in a small window around the value instead, as above.

To view the result in napari from Python, install `napari-spatialdata[all]`
in the same environment and run:

```python
from napari_spatialdata import Interactive

Interactive(sdata)
```

## Go deeper

- [Explore the result in Python](explore-the-output.ipynb): optical overlays,
  single spectra, slices and metadata, as a notebook.
- [Output format](output-format.md): every part of the result, in full
  detail.
