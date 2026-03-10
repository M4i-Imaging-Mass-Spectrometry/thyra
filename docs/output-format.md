# Output Format

Thyra converts MSI data into [SpatialData](https://spatialdata.scverse.org/) objects stored as Zarr directories. This page describes the structure and how to work with it.

## Loading a Dataset

```python
import spatialdata as sd

sdata = sd.read_zarr("output.zarr")

print("Tables:", list(sdata.tables.keys()))
print("Images:", list(sdata.images.keys()))
print("Shapes:", list(sdata.shapes.keys()))
```

## Structure Overview

A converted dataset contains the following elements:

| Element | Key Pattern | Description |
|---------|------------|-------------|
| **Table** | `{dataset_id}_z{z}` | AnnData with intensity matrix (pixels x m/z), coordinates in `.obs`, m/z axis in `.var` |
| **TIC Image** | `{dataset_id}_z{z}_tic` | 2D total ion current image, shape `(1, y, x)` |
| **Pixel Shapes** | `{dataset_id}_z{z}_pixels` | GeoDataFrame with pixel box geometries |
| **Optical Images** | `{dataset_id}_optical_{name}` | Microscopy images (when available) |

For 3D mode (`--handle-3d`), the `_z{z}` suffix is dropped.

## TIC Images

The TIC (Total Ion Current) image shows the sum of all intensities per pixel:

```python
import numpy as np

tic_key = [k for k in sdata.images if k.endswith("_tic")][0]
tic_array = np.asarray(sdata.images[tic_key])[0]  # drop channel dim -> (y, x)

print(f"Shape: {tic_array.shape}")
print(f"Range: {tic_array.min():.0f} -- {tic_array.max():.0f}")
```

```python
import matplotlib.pyplot as plt

plt.imshow(tic_array, cmap="viridis")
plt.colorbar(label="TIC Intensity")
plt.show()
```

## Optical Images

When converted with `--include-optical` (the default for Bruker data), microscopy images are stored alongside the MSI data. The primary optical image has an Identity transform -- it defines the coordinate system that the TIC affine maps into.

```python
optical_keys = [k for k in sdata.images if "optical" in k]
print("Optical images:", optical_keys)

if optical_keys:
    opt_image = np.asarray(sdata.images[optical_keys[0]])
    # Shape is (channels, y, x) -- transpose for display
    opt_rgb = np.moveaxis(opt_image[:3], 0, -1)
    plt.imshow(opt_rgb)
    plt.show()
```

### TIC-to-Optical Overlay

The TIC image has an affine transform that maps it into the optical image's coordinate space:

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

## Mass Spectrum Data

### Average Spectrum

```python
table_key = list(sdata.tables.keys())[0]
msi_table = sdata.tables[table_key]

mz_values = msi_table.var["mz"].values
avg_spectrum = msi_table.uns["average_spectrum"]

plt.plot(mz_values, avg_spectrum, linewidth=0.5)
plt.xlabel("m/z")
plt.ylabel("Average Intensity")
plt.show()
```

### Raw Intensity Matrix

The intensity matrix is stored as a sparse matrix. Each row is one pixel, each column is one m/z bin:

```python
import scipy.sparse as sparse

X = msi_table.X  # sparse (pixels x m/z)
print(f"Shape: {X.shape}")
print(f"Non-zero: {X.nnz:,} ({X.nnz / (X.shape[0] * X.shape[1]) * 100:.2f}%)")
```

### Ion Images

To visualize the spatial distribution of a specific m/z value:

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

## Pixel Coordinates

Coordinates are stored in the table's `.obs` DataFrame:

```python
print("Columns:", list(msi_table.obs.columns))

# Available columns:
# x, y          -- raster grid coordinates (integers)
# spatial_x/y   -- physical coordinates in micrometers
# region        -- SpatialData region key (categorical)
# region_number -- acquisition region number (integer)
# instance_key  -- pixel identifier (string)
```

## Regions

Datasets acquired from multi-region slides store region information in two places:

### Per-Pixel Region Number

The `region_number` column in `obs` indicates which acquisition region each pixel belongs to:

```python
print(msi_table.obs["region_number"].value_counts())
```

### Region Summary

Region metadata (including names from the instrument's Area definitions) is stored as JSON in `uns`:

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

## 3D Data / Z-Slices

By default, Thyra stores each z-slice as a separate table and TIC image with `_z{i}` suffixes:

```python
slice_tables = sorted(k for k in sdata.tables if "_z" in k)
print(f"{len(slice_tables)} z-slices: {slice_tables}")

# Access a single slice
z0_table = sdata.tables[slice_tables[0]]
print(f"Slice 0: {z0_table.shape}")
```

When converted with `--handle-3d`, all slices are combined into a single table with `x`, `y`, `z` coordinates in obs.

## Dataset Metadata

Global metadata is stored in `sdata.attrs`:

```python
if "pixel_size_x_um" in sdata.attrs:
    print(f"Pixel size: {sdata.attrs['pixel_size_x_um']} um")

if "msi_dataset_info" in sdata.attrs:
    info = sdata.attrs["msi_dataset_info"]
    print(f"Dimensions: {info.get('dimensions_xyz')}")
    print(f"Non-empty pixels: {info.get('non_empty_pixels'):,}")
```

Table-level metadata (instrument info, acquisition parameters) is in `msi_table.uns`:

```python
print("uns keys:", list(msi_table.uns.keys()))
```
