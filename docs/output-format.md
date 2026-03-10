# Output Format

Thyra converts MSI data into [SpatialData](https://spatialdata.scverse.org/)
objects stored as Zarr directories. This page describes what the output
contains and how to work with it.

---

## Loading a Dataset

```python
import spatialdata as sd

sdata = sd.read_zarr("output.zarr")

print("Tables:", list(sdata.tables.keys()))
print("Images:", list(sdata.images.keys()))
print("Shapes:", list(sdata.shapes.keys()))
```

---

## Structure Overview

A converted dataset contains the following elements:

| Element | Key Pattern | Description |
|---------|------------|-------------|
| **Table** | `{dataset_id}_z{z}` | AnnData with intensity matrix (pixels x m/z), coordinates in `.obs`, m/z axis in `.var` |
| **TIC Image** | `{dataset_id}_z{z}_tic` | 2D total ion current image, shape `(1, y, x)` |
| **Pixel Shapes** | `{dataset_id}_z{z}_pixels` | GeoDataFrame with pixel box geometries |
| **Optical Images** | `{dataset_id}_optical_{name}` | Microscopy images (when available) |

!!! note "3D mode"
    When converted with `--handle-3d`, the `_z{z}` suffix is dropped and all
    slices are merged into a single table with `x`, `y`, `z` coordinates in `.obs`.

!!! tip "Default dataset ID"
    The default `dataset_id` is `msi_dataset`, so typical keys look like
    `msi_dataset_z0`, `msi_dataset_z0_tic`, etc. Change it with `--dataset-id`.

---

## TIC Images

The TIC (Total Ion Current) image shows the sum of all intensities per pixel.

```python
import numpy as np
import matplotlib.pyplot as plt

tic_key = [k for k in sdata.images if k.endswith("_tic")][0]
tic_array = np.asarray(sdata.images[tic_key])[0]  # drop channel dim -> (y, x)

plt.imshow(tic_array, cmap="viridis")
plt.colorbar(label="TIC Intensity")
plt.title("Total Ion Current")
plt.show()
```

---

## Optical Images

When converted with `--include-optical` (the default for Bruker data),
microscopy images are stored alongside the MSI data.

```python
optical_keys = [k for k in sdata.images if "optical" in k]
print("Optical images:", optical_keys)

if optical_keys:
    opt_image = np.asarray(sdata.images[optical_keys[0]])
    # Shape is (channels, y, x) -- transpose for display
    opt_rgb = np.moveaxis(opt_image[:3], 0, -1)
    plt.imshow(opt_rgb)
    plt.title("Optical Image")
    plt.show()
```

### TIC-to-Optical Overlay

The TIC image carries an affine transform that maps it into the optical image's
coordinate space. This means overlays work automatically in tools like napari.

To inspect the transform:

```python
from spatialdata.transformations import get_transformation

tic_element = sdata.images[tic_key]
transforms = get_transformation(tic_element, get_all=True)
affine = list(transforms.values())[0]

matrix = np.array(affine.to_affine_matrix(
    input_axes=("x", "y"), output_axes=("x", "y")
))
print(f"Scale: {matrix[0,0]:.2f}x, {matrix[1,1]:.2f}x")
print(f"Offset: ({matrix[0,2]:.0f}, {matrix[1,2]:.0f})")
```

!!! info "How alignment works"
    The optical image has an Identity transform and defines the reference
    coordinate system. The TIC image has an Affine transform (scale + offset)
    that positions it in the optical coordinate space. This comes from the
    teaching point calibration in the `.mis` file (Bruker data).

---

## Mass Spectrum Data

### Average Spectrum

Each table stores a pre-computed average spectrum in `uns`:

```python
table_key = list(sdata.tables.keys())[0]
msi_table = sdata.tables[table_key]

mz_values = msi_table.var["mz"].values
avg_spectrum = msi_table.uns["average_spectrum"]

plt.plot(mz_values, avg_spectrum, linewidth=0.5)
plt.xlabel("m/z")
plt.ylabel("Average Intensity")
plt.title("Average Mass Spectrum")
plt.show()
```

### Intensity Matrix

The intensity matrix is stored as a sparse matrix. Each row is one pixel, each
column is one m/z bin:

```python
X = msi_table.X  # sparse (pixels x m/z)
print(f"Shape: {X.shape}")
print(f"Non-zero: {X.nnz:,} ({X.nnz / (X.shape[0] * X.shape[1]) * 100:.2f}%)")
```

!!! tip "Sparse format"
    The default storage is CSC (Compressed Sparse Column), which is fast for
    extracting ion images (column = one m/z across all pixels). If you need fast
    per-pixel access, convert with `--sparse-format csr`.

### Ion Images

To visualise the spatial distribution of a specific m/z value:

```python
target_mz = 760.5
mz_idx = np.abs(mz_values - target_mz).argmin()

# Extract column from sparse matrix
ion_values = np.asarray(X[:, mz_idx].toarray()).flatten()

# Reconstruct image from pixel coordinates
x_coords = msi_table.obs["x"].values.astype(int)
y_coords = msi_table.obs["y"].values.astype(int)

ion_image = np.zeros((y_coords.max() + 1, x_coords.max() + 1))
ion_image[y_coords, x_coords] = ion_values

plt.imshow(ion_image, cmap="hot")
plt.colorbar(label="Intensity")
plt.title(f"m/z {mz_values[mz_idx]:.4f}")
plt.show()
```

### Individual Spectra

```python
pixel_idx = 0
spectrum = X[pixel_idx].toarray().flatten()

plt.plot(mz_values, spectrum, linewidth=0.5)
plt.xlabel("m/z")
plt.ylabel("Intensity")
plt.title(f"Pixel {pixel_idx}")
plt.show()
```

---

## Pixel Coordinates

Coordinates are stored in the table's `.obs` DataFrame:

```python
print("Columns:", list(msi_table.obs.columns))
```

| Column | Type | Description |
|--------|------|-------------|
| `x`, `y` | int | Raster grid coordinates (pixel indices) |
| `spatial_x`, `spatial_y` | float | Physical coordinates in micrometers |
| `region` | categorical | SpatialData region key |
| `region_number` | int | Acquisition region number |
| `instance_key` | str | Unique pixel identifier |

---

## Regions

Datasets acquired from multi-region slides (e.g., multiple tissue sections on
one slide) store region information in two places.

### Per-Pixel Region Number

The `region_number` column in `.obs` indicates which acquisition region each
pixel belongs to:

```python
print(msi_table.obs["region_number"].value_counts())
```

### Region Summary

Region metadata -- including human-readable names from the instrument's Area
definitions -- is stored as JSON in `uns`:

```python
import json

regions = json.loads(msi_table.uns["regions"])
for r in regions:
    print(f"Region {r['region_number']}: {r.get('name', 'unnamed')} "
          f"({r['n_spectra']:,} spectra)")
```

Example output:

```
Region 0: E2506 (104,321 spectra)
Region 1: Matrix (1,053 spectra)
```

!!! tip "Filtering by region"
    To work with only one region:
    ```python
    mask = msi_table.obs["region_number"] == 0
    tissue_table = msi_table[mask]
    ```

---

## 3D Data / Z-Slices

By default, Thyra stores each z-slice as a separate table and TIC image:

```python
slice_tables = sorted(k for k in sdata.tables if "_z" in k)
print(f"{len(slice_tables)} z-slices: {slice_tables}")

# Access a single slice
z0_table = sdata.tables[slice_tables[0]]
print(f"Slice 0: {z0_table.shape}")
```

When converted with `--handle-3d`, all slices are combined into a single table
with `x`, `y`, `z` coordinates in `.obs`.

---

## Dataset Metadata

### Global metadata

Stored in `sdata.attrs`:

```python
if "pixel_size_x_um" in sdata.attrs:
    print(f"Pixel size: {sdata.attrs['pixel_size_x_um']} um")

if "msi_dataset_info" in sdata.attrs:
    info = sdata.attrs["msi_dataset_info"]
    print(f"Dimensions: {info.get('dimensions_xyz')}")
    print(f"Non-empty pixels: {info.get('non_empty_pixels'):,}")
```

### Table-level metadata

Instrument info, acquisition parameters, and resampling config are in
`msi_table.uns`:

```python
print("uns keys:", list(msi_table.uns.keys()))
```

---

## Recipes

### Side-by-side TIC and ion image

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.imshow(tic_array, cmap="viridis")
ax1.set_title("TIC")

target_mz = 760.5
mz_idx = np.abs(mz_values - target_mz).argmin()
ion_values = np.asarray(X[:, mz_idx].toarray()).flatten()

x_coords = msi_table.obs["x"].values.astype(int)
y_coords = msi_table.obs["y"].values.astype(int)
ion_image = np.zeros((y_coords.max() + 1, x_coords.max() + 1))
ion_image[y_coords, x_coords] = ion_values

ax2.imshow(ion_image, cmap="hot")
ax2.set_title(f"m/z {mz_values[mz_idx]:.2f}")

plt.tight_layout()
plt.show()
```

### Export ion image to TIFF

```python
from PIL import Image

# Normalise to 0-255
ion_norm = (ion_image / ion_image.max() * 255).astype(np.uint8)
Image.fromarray(ion_norm).save("ion_image.tiff")
```

### Top N most intense m/z values

```python
avg = msi_table.uns["average_spectrum"]
top_n = 10
top_indices = np.argsort(avg)[-top_n:][::-1]

for idx in top_indices:
    print(f"  m/z {mz_values[idx]:.4f}  avg intensity {avg[idx]:.1f}")
```

### Summary statistics

```python
print(f"Dataset: {table_key}")
print(f"  Pixels: {msi_table.n_obs:,}")
print(f"  m/z bins: {msi_table.n_vars:,}")
print(f"  m/z range: {mz_values.min():.2f} -- {mz_values.max():.2f}")
print(f"  Sparsity: {(1 - X.nnz / (X.shape[0] * X.shape[1])) * 100:.1f}%")
if "pixel_size_x_um" in sdata.attrs:
    print(f"  Pixel size: {sdata.attrs['pixel_size_x_um']} um")
```
