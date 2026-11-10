# Analyse the result with scanpy and squidpy

[scanpy](https://scanpy.readthedocs.io) and
[squidpy](https://squidpy.readthedocs.io) are free Python tools for tables
like the one in a Thyra result. scanpy groups pixels with similar spectra.
squidpy measures patterns across the tissue. This page runs both on the
example data.

The example data is synthetic and small. Real datasets have many more
pixels.

## Install the packages

Install three packages beside Thyra, in the Python you work in (see
[Already use Python?](install.md#already-use-python)):

```bash
pip install thyra "scanpy[leiden]" squidpy spatialdata-plot
```

spatialdata-plot draws a result on the tissue. `[leiden]` adds the
clustering method this page uses.

## Read the result

Make the example data and convert it:

```bash
thyra-example-data example_data/synthetic_brain.imzML
thyra example_data/synthetic_brain.imzML example_data/synthetic_brain.zarr
```

Then open the result in Python:

```python
import scanpy as sc
import spatialdata as sd
import spatialdata_plot  # adds sdata.pl, used below
import squidpy as sq

sdata = sd.read_zarr("example_data/synthetic_brain.zarr")
table = sdata.tables["msi_dataset_z0"]
```

`table` holds the spectra: one row per pixel and one column per m/z bin.
scanpy and squidpy were written for genes, so they call a column a "gene".
Here, a gene is an m/z bin.

## Normalise, reduce and cluster

```python
sc.pp.normalize_total(table)
sc.pp.log1p(table)
sc.pp.highly_variable_genes(table, n_top_genes=500)
sc.pp.pca(table)
sc.pp.neighbors(table)
sc.tl.leiden(table, resolution=0.2, flavor="igraph")
sc.tl.umap(table)
sc.tl.rank_genes_groups(table, groupby="leiden")
```

Line by line:

1. `normalize_total` scales every spectrum to the same TIC (total ion
   current), the sum of its intensities.
2. `log1p` takes the logarithm, so a few strong peaks count for less.
3. `highly_variable_genes` marks the 500 m/z bins that differ most between
   pixels.
4. `pca` reduces those bins to 50 numbers per pixel.
5. `neighbors` finds the pixels with the most similar spectra, for every
   pixel.
6. `leiden` puts similar pixels into groups, called clusters. A higher
   `resolution` gives more clusters.
7. `umap` places the pixels on a flat map, similar pixels close together.
   Draw it with `sc.pl.umap(table, color="leiden")`.
8. `rank_genes_groups` ranks the m/z bins that set each cluster apart.

## Plot the clusters on the tissue

```python
sdata.pl.render_shapes("msi_dataset_z0_pixels", color="leiden").pl.show(
    coordinate_systems="global"
)
```

Every pixel has the colour of its cluster. The example gives three
clusters: the background, the outer region and the inner region.

## Spatial neighbours and Moran's I

The spatial neighbours of a pixel are the pixels beside it on the tissue.
squidpy finds them from the pixel outlines in the result:

```python
sq.gr.spatial_neighbors_grid(
    sdata,
    elements_to_coordinate_systems={"msi_dataset_z0_pixels": "global"},
    table_key="msi_dataset_z0",
    n_neighs=4,
)
```

squidpy may warn that an index is not sorted. The result is not affected.

Moran's I is a number for one m/z bin. Near 1, neighbouring pixels have
similar intensities, so the bin shows a clear structure. Near 0, the bin
shows no pattern.

```python
bins = table.var_names[table.var["highly_variable"]]
sq.gr.spatial_autocorr(table, mode="moran", genes=bins)
print(table.uns["moranI"].head())
```

The bins with the highest Moran's I come first.

## One ion image

An ion image shows the intensity of one m/z bin in every pixel. Draw the
bin with the highest Moran's I:

```python
best = table.uns["moranI"].index[0]
print("m/z", table.var.loc[best, "mz"])
sdata.pl.render_shapes("msi_dataset_z0_pixels", color=best).pl.show(
    coordinate_systems="global"
)
```

The values are the normalised intensities from the first step. In the
example, this bin is brightest around the section.

## Where to go next

- [scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials/index.html)
  and [squidpy tutorials](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/index.html):
  more analyses on the same table.
- [Explore the result in Python](explore-the-output.ipynb): optical
  overlays, single spectra and metadata.
- [For tool developers](for-tool-developers.md): what a tool can rely on
  in a result.
