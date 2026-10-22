# Ion mobility and MS/MS

Some instruments measure more than m/z. A Bruker timsTOF or a Waters Synapt
or Cyclic IMS can also separate ions by their ion mobility. A timsTOF can break selected
ions into fragments (MS/MS). Thyra keeps both, next to the usual spectra.

## Which data has them

| Your data | What Thyra keeps |
|---|---|
| Bruker timsTOF `.d` folder, with TIMS switched on | ion mobility |
| Bruker timsTOF PASEF or iprm-PASEF imaging run | ion mobility and MS/MS |
| imzML in continuous mode with a mobility array, as TIMSImaging writes | ion mobility |
| imzML in processed mode with a mobility array, as TIMSCONVERT writes | ion mobility |
| Waters Synapt `.raw` folder, with ion mobility on | ion mobility (drift time) |
| Waters Cyclic IMS `.raw` folder | ion mobility (drift time) |

You do not need to switch anything on. Thyra finds the mobility and the
MS/MS data by itself.

## What you get

| Part | What it holds | When |
|---|---|---|
| **The spectra**, `msi_dataset_z0` | every pixel's spectrum, added up over mobility | always |
| **Mass-mobility heatmap**, inside the spectra table | the whole image in one picture: m/z against mobility | when the data has mobility |
| **Mobility table**, `msi_dataset_z0_mobility` | one column for each pair of m/z and mobility | continuous imzML: always. timsTOF, processed imzML and Waters: add `--mobility-grid` |
| **MS/MS table**, `msi_dataset_z0_msms` | the fragments, grouped by the ion they came from | when every pixel fragments the same list of ions |

The spectra table looks the same as for any other data. Tools that know
nothing about mobility can still read it.

## Look at the heatmap first

The heatmap shows whether mobility separates anything in your sample. Two
bright spots at the same m/z but at different mobility are two different
ions. The mobility table can show each of them as its own image.

```python
import matplotlib.pyplot as plt
import numpy as np
import spatialdata as sd

sdata = sd.read_zarr("my_slide.zarr")
heat = sdata.tables["msi_dataset_z0"].uns["mobility_heatmap"]
plt.pcolormesh(heat["mobility_edges"], heat["mz_edges"], np.log1p(heat["counts"]))
plt.xlabel("1/K0 (V s/cm2)")
plt.ylabel("m/z")
plt.show()
```

## Get the mobility table from timsTOF or processed imzML

For timsTOF data and processed imzML, ask for the mobility table when you
convert:

```bash
thyra my_slide.d my_slide.zarr --mobility-grid
```

The result is larger and slower to make. Continuous imzML
files get the table without asking, if all spectra share one mobility array.

??? advanced "Advanced: how this table is built"
    Every pixel measures mobility at slightly different values. So Thyra puts
    them on shared mobility channels, the same channels as the heatmap: 256
    on a timsTOF, one per drift bin on a Waters run.
    A box drawn on the heatmap then picks the same channels in the table.
    `--mobility-bins`, `--mobility-min` and `--mobility-max` change the
    channels. See
    [Output format](output-format.md#the-same-table-from-a-common-mobility-grid).

## Two ions at one m/z

This makes two images from one m/z value: one for each mobility range.

```python
mob = sdata.tables["msi_dataset_z0_mobility"]
mz = mob.var["mz"].to_numpy()
k0 = mob.var["mobility"].to_numpy()
near = np.abs(mz - 600.25) < 0.01
low = np.asarray(mob.X[:, near & (k0 < 1.3)].sum(axis=1)).ravel()
high = np.asarray(mob.X[:, near & (k0 >= 1.3)].sum(axis=1)).ravel()
```

`low` and `high` hold one value per pixel. Turn them into images with the
pixel positions, as in [Look at the result](look-at-the-result.md#open-it-in-python).

## Collision cross section

The mobility table also gives each column's collision cross section (CCS),
in square angstroms:

```python
print(mob.var[["mz", "mobility", "ccs"]].head())
```

Thyra assumes every ion has charge 1, which is almost always true for MALDI.
`mob.uns["ccs"]` lists this and the other assumptions. See
[Output format](output-format.md#collision-cross-section) for how to redo it
for another charge.

## One ion's fragments

In an MS/MS run, the instrument picks a list of ions, called precursors. It
breaks each one into fragments at every pixel. The MS/MS table keeps each
precursor's fragments together, so every fragment has its own image.

```python
msms = sdata.tables["msi_dataset_z0_msms"]
print(msms.var.groupby("precursor_index")["precursor_mz"].first())
block = (msms.var["precursor_index"] == 0).to_numpy()
image = np.asarray(msms.X[:, block].sum(axis=1)).ravel()
```

The first line lists the precursors. The last line adds up all fragments of
precursor 0, one value per pixel.

??? advanced "Advanced: two precursors at the same m/z"
    A method can pick the same m/z twice, at two mobility values. That is how
    isomers are targeted. Thyra keeps them apart, so always select a
    precursor by `precursor_index`, not by its m/z. `precursor_mobility`
    tells the two apart. See
    [Output format](output-format.md#demultiplexed-msms-table).

??? advanced "Advanced: reading these tables in your own tool"
    Recognise the tables by their columns, not by their names. A table with
    a `mobility` column is a mobility table. A table with a `precursor_mz`
    column is an MS/MS table. In both, `var["mz"]` repeats, so code that
    expects every m/z once must not be pointed at them.
    [For tool developers](for-tool-developers.md) lists the rules.

## Not covered yet

- Older Waters Synapt G1 runs with ion mobility on are refused. The Waters
  library cannot read them safely.
- CCS is computed from 1/K0 only, never from a drift time.
- On a Cyclic IMS run with several passes, ions outside the targeted mass
  range can lap. Their drift time then says little about their mobility.
- diaPASEF and data-dependent PASEF runs give the spectra table only. Their
  list of precursors changes from pixel to pixel.
- MALDI iprm-PASEF runs from the released Bruker software have been tested
  on a synthetic file only. If you have one, please
  [get in touch](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/issues).

## Go deeper

[Output format](output-format.md#ion-mobility) describes every part of these
tables in full detail.
