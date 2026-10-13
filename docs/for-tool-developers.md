# For tool developers

This page is for authors of MSI analysis tools who want to read a Thyra
store. It states what a store promises and what it does not. For a tour of
everything a store holds, see [Output format](output-format.md).

A Thyra store is a [SpatialData](https://spatialdata.scverse.org/) store
written as a Zarr folder. Its spectra are AnnData tables: one row per pixel,
one column per m/z bin, the intensities in a sparse matrix. Its images, such
as the TIC image, are OME-Zarr images.

## Read a store

With `spatialdata`:

```python
import numpy as np
import spatialdata as sd

sdata = sd.read_zarr("sample.zarr")
table = sdata.tables["msi_dataset_z0"]       # one row per pixel
X = table.X                                  # scipy CSC matrix, pixels x m/z bins
mz = table.var["mz"].to_numpy()              # the mass axis, one value per column
x = table.obs["x"].to_numpy()                # raster column of each row
y = table.obs["y"].to_numpy()                # raster row of each row
tic = np.asarray(sdata.images["msi_dataset_z0_tic"])[0]   # (y, x)
print(table.uns["msi_metadata"]["schema_version"])
```

Without `spatialdata`, with `anndata` and `zarr` only:

```python
import anndata as ad
import numpy as np
import zarr

root = zarr.open_group("sample.zarr", mode="r")
print(list(root["tables"]))                  # the table names
table = ad.io.read_zarr("sample.zarr/tables/msi_dataset_z0")
tic_group = root["images/msi_dataset_z0_tic"]
level = tic_group.attrs["ome"]["multiscales"][0]["datasets"][0]["path"]
tic = np.asarray(tic_group[level])[0]        # (y, x)
```

Read the image's level path from its OME metadata, as above. Do not assume
its name.

The default dataset id is `msi_dataset`. In a 2D conversion the first table
is `msi_dataset_z0` and its TIC image is `msi_dataset_z0_tic`. Users can
change the id, so list the tables rather than hard-coding a name. The
snippets on this page use the default names to stay short.

A 3D conversion (`--handle-3d`) differs: the table is `msi_dataset`, `obs`
gains `z`, and the TIC image is a volume of shape `(c, z, y, x)`. The
snippets on this page assume a 2D store. See
[Output format](output-format.md#3d-data-z-slices).

## The rules a tool can rely on

### Which table is which

A store can hold three kinds of table. Tell them apart by their `var`
columns, never by the element name:

```python
def table_kind(table):
    if "mobility" in table.var.columns:
        return "mobility"      # one column per (m/z, mobility) pair
    if "precursor_mz" in table.var.columns:
        return "msms"          # one column block per precursor
    return "msi"               # one column per m/z bin
```

No table carries both columns.

| | MSI table | Mobility table | MS/MS table |
|---|---|---|---|
| Marked by | neither column | `mobility` in `var` | `precursor_mz` in `var` |
| `var["mz"]` | strictly increasing, unique | non-decreasing, repeats | restarts at every precursor |
| Sorted by | `mz` | `(mz, mobility)` | `(precursor_mz, precursor_mobility, mz)` |
| Unique key of a column | `mz` | `(mz, mobility)` | `(precursor_index, mz)` |
| Other `var` columns | none today | `mz_index`, `mobility_index` | `precursor_index`, `mz_index`, `precursor_mobility` |

Only the MSI table has a unique, strictly increasing `var["mz"]`. Code that
binary-searches the mass axis must only be given the MSI table.

`precursor_mobility` is present only when the source states the precursors'
mobility. Without it, the MS/MS table is sorted by `(precursor_mz, mz)`.
A reader may add other `var` columns to the MSI table, such as a native
flight-time axis, so do not assume `mz` is the only one.

### How the tables connect

- A sibling table has the same rows as its MSI table, in the same order,
  with the same `obs`.
- `mz_index` is the column of the MSI table that a sibling feature belongs
  to. Use it to relate the two mass axes. On a mobility table read from
  imzML, a feature outside the MSI table's mass range is given the first or
  last column.
- The MSI table names its siblings in its metadata block:
  `uns["msi_metadata"]["ms_analysis"]["ion_mobility"]["resolved_table"]` and
  `uns["msi_metadata"]["ms_analysis"]["fragmentation"]["resolved_table"]`.
  Each key is absent when there is no such table.
- `precursor_index` identifies a precursor within one store only. Two
  precursors can share an m/z (an isomer pair), so select by
  `precursor_index`. To match precursors across stores, compare
  `precursor_mz` and `precursor_mobility` within a tolerance.

### Rows and values

- `obs["x"]` and `obs["y"]` are integer raster positions. They need not
  start at 0. Rows follow the raster, left to right and top to bottom.
  `obs["acquisition_order"]`, when present, gives the time order.
- `obs["spatial_x"]` and `obs["spatial_y"]` are in micrometres. When the
  source states no pixel size, 1 µm is assumed, and
  `uns["msi_metadata"]["provenance"]["pixel_size_source"]` is `"default"`.
- A position the instrument never measured has no row. A position measured
  twice has one row holding the sum of both spectra.
- `X` is stored as CSC on every route, with sorted row indices in each
  column. Every stored value is finite and non-negative. Zeros are not
  stored.
- The `var` and `obs` column names in the
  [metadata schema](metadata-schema.md#var-column-conventions) are reserved.
  They are never reused with a different meaning.

### What `thyra validate` checks

`thyra validate sample.zarr` exits with status 1 when a table breaks its
rule, so it can gate a pipeline. It reads only `var` and the metadata
blocks, never the matrix.

| Table | Checked |
|---|---|
| every table | `var["mz"]` present, numeric and finite; a valid `uns["msi_metadata"]` block |
| MSI table | `mz` strictly increasing; it matches the mass axis fingerprint, when one is recorded |
| Mobility table | `mz` non-decreasing; `(mz, mobility)` pairs unique and sorted |
| MS/MS table | `precursor_mz` non-decreasing; each `precursor_index` is one contiguous block; `mz` strictly increasing inside each block |

## What is stable, and how changes are versioned

### Versioned: the metadata block

Thyra writes `uns["msi_metadata"]` in every table, and it stays at that
location. If Thyra cannot build the block, it logs a warning and writes
none; `thyra validate` then fails the store. The block's `schema_version`
follows semantic versioning:

- Adding an optional field raises the minor version.
- Renaming, removing or retyping a field, or making one required, raises the
  major version.
- A validator for major version N rejects other majors, accepts older
  minors and warns on newer ones.

Each schema version from 0.6.0 on is published as JSON Schema at a fixed
address, for example
`https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/schema/0.12.0/msi_metadata.schema.json`.
A published version is never edited, and there is no `latest`. The same
file ships inside the `thyra` package. See
[Metadata schema](metadata-schema.md#versioning).

### Versioned: the package

The `thyra` package version follows
[semantic versioning](contributing.md#versioning-policy), set from the
commit messages. A new feature makes a minor release, and a fix or a
speed-up makes a patch release. Only a change marked as breaking makes a
major release.

The package version does not version the store's contents. A minor or patch
release can change stored values, for example when a default changes or a
wrong value is fixed. Read the [changelog](changelog.md) before you upgrade,
not only for a major release. [Design decisions](design-decisions.md) gives
the reason behind each default.

To see what made a store:

- `uns["msi_metadata"]["provenance"]["thyra_version"]` is the Thyra version
  that wrote it. Note which versions your tool is tested against.
- `json.loads(uns["msi_metadata"]["processing"])` lists the processing steps
  and their parameters. The list is stored as a JSON string.
- `uns["msi_metadata"]["ms_analysis"]["spectrum_representation"]` is
  `"profile"` or `"centroid"`: what the reader delivered. It is absent when
  the source does not say.
- `uns["msi_metadata"]["provenance"]["mass_axis"]` is a SHA-256 of the MSI
  table's `var["mz"]`, with its length. Two stores with equal values share
  the same mass axis, bin for bin. Sibling tables carry their MSI table's
  value. See
  [design decision D42](design-decisions.md#d42-a-store-records-a-sha-256-of-its-mass-axis).

### Not promised

- **Element names.** The dataset id is a user setting. Find tables by their
  columns and by `resolved_table`.
- **The mass axis.** The default axis follows the acquisition's mass range,
  so two stores usually differ. Compare the fingerprint, or `var["mz"]`
  itself, before reading across stores.
- **The sibling feature axes.** In a timsTOF mobility grid table and in an
  MS/MS table, only features that carry signal become columns. Two
  conversions with the same settings can give different `var`.
- **Other `uns` blocks**, such as `mobility_axis`, `mobility_heatmap` and
  `msms_schedule`. They are described in [Output format](output-format.md)
  but carry no version of their own.
- **The source's own vocabulary.** `format_specific`, `acquisition_params`,
  `instrument_info` and `raw_metadata` keep what the source said, so their
  keys differ by format.
- **Zarr storage details**, such as chunk sizes and compression.

## Recipes

The recipes below run as written against stores converted from the test
fixtures in the Thyra repository.

### An ion image

Sum the columns inside an m/z window. On a resampled axis the bins are
narrow, so a single bin is usually empty.

```python
import numpy as np
import spatialdata as sd


def ion_image(table, target_mz, half_width):
    mz = table.var["mz"].to_numpy()          # MSI table only: strictly increasing
    lo, hi = np.searchsorted(mz, [target_mz - half_width, target_mz + half_width])
    values = np.asarray(table.X[:, lo:hi].sum(axis=1)).ravel()
    x = table.obs["x"].to_numpy()
    y = table.obs["y"].to_numpy()
    image = np.full((y.max() + 1, x.max() + 1), np.nan)   # NaN = not measured
    image[y, x] = values
    return image


sdata = sd.read_zarr("sample.zarr")
image = ion_image(sdata.tables["msi_dataset_z0"], 600.25, 0.05)
```

### A mobility window

Select the m/z window on the MSI table, then keep the mobility table's
features inside it whose mobility falls in a range.

```python
import numpy as np
import spatialdata as sd

sdata = sd.read_zarr("sample.zarr")
msi = sdata.tables["msi_dataset_z0"]
name = msi.uns["msi_metadata"]["ms_analysis"]["ion_mobility"]["resolved_table"]
mob = sdata.tables[name]

mz = msi.var["mz"].to_numpy()
lo, hi = np.searchsorted(mz, [600.25 - 0.05, 600.25 + 0.05])
mz_index = mob.var["mz_index"].to_numpy()
mobility = mob.var["mobility"].to_numpy()        # 1/K0 for TIMS

selected = (mz_index >= lo) & (mz_index < hi) & (mobility >= 1.3) & (mobility < 1.4)
values = np.asarray(mob.X[:, selected].sum(axis=1)).ravel()   # one value per pixel
```

On a timsTOF grid table with the default settings, the same window over all
mobility values gives the MSI table's ion image.
`uns["mobility_marginal"]["current_ratio"]` on the grid table records how
closely the two tables agree.

### A precursor block

Look the precursor up once, then select its block by `precursor_index`.
Never compare `precursor_mz` with `==`.

```python
import numpy as np
import spatialdata as sd

sdata = sd.read_zarr("sample.zarr")
msi = sdata.tables["msi_dataset_z0"]
name = msi.uns["msi_metadata"]["ms_analysis"]["fragmentation"]["resolved_table"]
msms = sdata.tables[name]

precursor_mz = msms.var.groupby("precursor_index")["precursor_mz"].first()
index = (precursor_mz - 500.0).abs().idxmin()  # nearest precursor to m/z 500
# Two precursors can share an m/z: then compare precursor_mobility as well.

block = (msms.var["precursor_index"] == index).to_numpy()
fragment_mz = msms.var["mz"].to_numpy()[block]
precursor_image = np.asarray(msms.X[:, block].sum(axis=1)).ravel()  # per pixel
fragment_images = msms.X[:, block]           # one column per fragment
```

`uns["demultiplexed_current"]["current_ratio"]` on the MS/MS table says how
much of the MSI table's ion current the blocks hold together.

### Do two stores share a mass axis?

```python
import spatialdata as sd


def mass_axis(path):
    table = sd.read_zarr(path).tables["msi_dataset_z0"]
    return table.uns["msi_metadata"]["provenance"].get("mass_axis")


a, b = mass_axis("first.zarr"), mass_axis("second.zarr")
same_axis = a is not None and a == b
```

A store without the field predates the fingerprint. Compare `var["mz"]`
directly for those.
