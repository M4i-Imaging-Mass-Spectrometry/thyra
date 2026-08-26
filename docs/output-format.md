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
    slices are merged into a single table with `x`, `y`, `z` coordinates in
    `.obs`. The TIC image becomes a single **volume** of shape `(c, z, y, x)`
    rather than one 2D image per slice — see
    [3D Data / Z-Slices](#3d-data-z-slices).

    The pixel shapes stay two-dimensional: one flat square per pixel per slice,
    carrying no depth. Read a pixel's depth from the table's `spatial_z`
    instead. See [Pixel footprints in a volume](#pixel-footprints-in-a-volume)
    for why, and [Coordinate systems](coordinate-systems.md) for how the three
    elements line up.

!!! tip "Default dataset ID"
    The default `dataset_id` is `msi_dataset`, so typical keys look like
    `msi_dataset_z0`, `msi_dataset_z0_tic`, etc. Change it with `--dataset-id`.

!!! info "Coordinate systems"
    Every element above carries a transform to a single ``"global"``
    coordinate system, and Thyra writes a self-describing
    ``coordinate_systems`` metadata attr at the zarr top level so
    consumers know what ``"global"`` is in (micrometers or pixels).
    See [Coordinate Systems](coordinate-systems.md) for the contract
    and how to read it.

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

### Per-Region Average Spectrum

For multi-region datasets, Thyra also stores a mean spectrum per acquisition
region in `uns["average_spectrum_per_region"]`. Each key is the region number
(as a string), and the value is a 1-D array matching the m/z axis.

```python
if "average_spectrum_per_region" in msi_table.uns:
    per_region = msi_table.uns["average_spectrum_per_region"]
    for region_id, spectrum in per_region.items():
        plt.plot(mz_values, spectrum, label=f"Region {region_id}", linewidth=0.5)
    plt.xlabel("m/z")
    plt.ylabel("Average Intensity")
    plt.legend()
    plt.title("Average Spectrum per Region")
    plt.show()
```

!!! note
    This key is only present when the dataset contains multiple acquisition
    regions. Single-region datasets only have the global `average_spectrum`.

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

### var columns

`var` always carries `mz` (numeric, finite, strictly increasing).
Readers whose native axis is not m/z keep that axis alongside it, and
annotation tools add `formula`, `adduct`, `annotation_source` and
`fdr` -- names fixed by the
[metadata schema](metadata-schema.md#var-column-conventions) so they
mean the same thing in every store. `thyra validate` checks the `mz`
contract on every table.

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

The DataFrame index is `instance_id` (a string pixel identifier).

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

### 3D mode (`--handle-3d`)

When converted with `--handle-3d`, all slices are combined into a single table
with `x`, `y`, `z` coordinates in `.obs`, and the per-slice TIC images are
replaced by one **volume**:

| Element | Key | Shape / dims |
|---------|-----|--------------|
| **Table** | `{dataset_id}` | pixels x m/z, with `x`, `y`, `z` and `spatial_x`, `spatial_y`, `spatial_z` in `.obs` |
| **TIC volume** | `{dataset_id}_tic` | `(c, z, y, x)` — one channel, then the three spatial axes |
| **Pixel Shapes** | `{dataset_id}_pixels` | GeoDataFrame of 2D pixel boxes — one per pixel per slice, no z |

```python
volume = sdata.images[f"{dataset_id}_tic"]
print(volume.dims)                      # ('c', 'z', 'y', 'x')

arr = np.asarray(volume.data)[0]        # drop channel -> (z, y, x)
plt.imshow(arr[0], cmap="viridis")      # first slice
```

Note the axis order is `(c, z, y, x)`, not `(c, x, y, z)`: index a slice with
`arr[z]`, not `arr[..., z]`.

#### Voxel depth

The volume carries a `Scale` to `"global"` built from **two** distinct numbers —
the in-plane pixel pitch for `x` and `y`, and the slice spacing for `z`:

```python
from spatialdata.transformations import get_transformation

axes = ("c", "z", "y", "x")
matrix = get_transformation(volume, to_coordinate_system="global").to_affine_matrix(
    input_axes=axes, output_axes=axes
)
print(matrix[1, 1])   # um per slice step
print(matrix[3, 3])   # um per pixel in x
```

The slice spacing comes from `--z-spacing`. When nothing supplies one, Thyra
falls back to the in-plane pitch and records that it did:

```python
cs = sdata.attrs["coordinate_systems"]["global"]
cs["z_spacing_um"]      # the number used
cs["z_spacing_source"]  # "manual" | "automatic" | "assumed_isotropic"
```

A `z_spacing_source` of `"assumed_isotropic"` means **nobody supplied a spacing
and the in-plane pitch was reused** — treat the depth as unknown rather than as
measured. Section thickness is set by the microtome, not by the raster, so the
two agree only by coincidence. See
[`--z-spacing`](cli.md#set-z-spacing-whenever-you-know-it).

!!! note "These keys only appear on volumes"
    `z_spacing_um` and `z_spacing_source` are written only when the store
    actually holds a multi-slice volume. Their absence is how a 2D store says it
    has no z axis, which is why `convention_version` stays at `1` — the keys are
    additive and a consumer that never reads 3D sees the schema it already
    knows.

#### Pixel footprints in a volume

The pixel polygons are **two-dimensional**, on every route including a volume.
A slice's depth is not on the geometry; read it from the table's `spatial_z`,
or from the volume's own `Scale`:

```python
shapes = sdata.shapes[f"{dataset_id}_pixels"]
print(shapes.geometry.iloc[0].has_z)    # False, on 2D and 3D alike

# Depth per pixel, in micrometres, from the table:
obs = sdata.tables[dataset_id].obs
print(obs.loc["0", "spatial_z"])        # z_index * z_spacing_um
```

Every footprint therefore resolves to the correct place in x and y, and carries
no claim about z. Overlaying the shapes on the volume is exact in-plane and
needs `spatial_z` to pick the slice.

!!! note "Why the footprints are flat"
    They briefly were not. v3.2.0 made them `POLYGON Z` at the depth of their
    slice, which is geometrically the more honest representation, and it broke
    `spatialdata`'s spatial queries: a bounding box enclosing an entire test
    volume returned 26 of 30 footprints and 26 of 30 table rows, with no
    exception and no warning. A z-restricted query returned the same rows
    whether z was inside or far outside the data.

    `spatialdata` asks for 2D here: `ShapesModel.validate` warns that a
    3-dimensional geometry column "could led to unexpected behaviors" and names
    `force_2d()` as the remedy. 3D shapes are not on its roadmap
    ([#109](https://github.com/scverse/spatialdata/issues/109) has been idle
    since June 2023), and the live 2.5D discussion
    ([#961](https://github.com/scverse/spatialdata/issues/961)) covers points,
    images and labels only. Serial-section MSI is 2.5D in that sense.

    So this is a deliberate trade: a documented gap in z, in exchange for
    queries that return every pixel. Expressing the depth as a `Translation` on
    flat geometry does not close the gap either — the transform silently drops
    z. Both negative results are pinned by
    `tests/unit/converters/test_3d_pixel_shapes_z.py`, which will fail if
    upstream changes. See [Coordinate systems](coordinate-systems.md).

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

### Provenance

`uns["essential_metadata"]` records what the store was made from, and how
the source was interpreted. It is the only place a converted dataset says
where it came from, so it is written the same way by every converter path:

```python
provenance = msi_table.uns["essential_metadata"]

print(provenance["source_path"])     # the file this was converted from
print(provenance["dimensions"])      # source grid, [x, y, z]
print(provenance["mass_range"])      # SOURCE m/z range, not the target axis
print(provenance["spectrum_type"])   # "centroid spectrum" / "profile spectrum"
print(provenance["thyra_version"])   # the Thyra that wrote the store
```

`mass_range` describes the source, not the resampled axis -- for the axis
the data actually sits on, read `msi_table.var["mz"]`.

Beside it, when the source format provides them:

| Key | Contents |
|-----|----------|
| `format_specific` | Vendor metadata (imzML file mode and UUID, FlexImaging areas, teaching points) |
| `acquisition_params` | Polarity, scan range, laser settings |
| `instrument_info` | Instrument model, serial, software version |
| `raw_metadata` | Source metadata as read, for round-trip fidelity |
| `regions` | Acquisition region summary, as a JSON string (see [Regions](#regions)) |

A section the source format has nothing for is omitted rather than written
empty, so `"instrument_info" not in uns` means "this format does not carry
it" rather than "it was carried and lost".

### Structured metadata: `uns["msi_metadata"]`

The sections above preserve what the source said, in the source's own
vocabulary. `uns["msi_metadata"]` is the normalised, versioned view of the
same facts: fixed field names, PSI-MS/NCBITaxon/UBERON/CHEBI ontology
terms, and a schema a validator can hold it to.

```python
block = msi_table.uns["msi_metadata"]

print(block["schema_version"])                    # "0.1.0"
print(block["ms_analysis"]["pixel_size_um"])      # {"x": 20.0, "y": 20.0}
print(block["provenance"]["source_format"])       # "imzml"
```

It is written by every converter path identically, validated by
`thyra validate`, and exported to a METASPACE submission by
`thyra export-metaspace`. See [Metadata Schema](metadata-schema.md) for
the full contract.

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
