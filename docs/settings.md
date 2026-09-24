# Change how Thyra converts

Thyra chooses its settings from your data, and most datasets need no options
at all. This page explains the options you are most likely to want, one task
at a time.

Add an option to the end of the command, after the result's name:

```bash
thyra my_slide.d my_slide.zarr --no-optical
```

## Set the pixel size

Thyra reads the pixel size from your files. If they do not record it, the
conversion stops before it writes anything. Give the size in micrometres:

```bash
thyra input.imzML output.zarr --pixel-size 20
```

The same size is used for both directions.

??? advanced "Advanced: pixels that are not square"
    A pixel size read from the file can differ between x and y, and Thyra
    keeps both. `--pixel-size` always gives square pixels.

## Convert one region

A Bruker timsTOF slide with several regions is converted as one result, with
every pixel labelled with its region. To convert only one region, add
`--region` and the region's name as flexImaging shows it:

```bash
thyra slide.d region_03.zarr --region 03
```

Thyra lists the regions it finds when it starts, so a first run shows you the
names. If the name matches none of them, Thyra stops and lists the names it
did find.

??? advanced "Advanced: region numbers, and the cropped photo"
    If the value matches no region name, Thyra reads it as the region's
    number in the Bruker database, which starts at 0. The names and the
    numbers often differ: region `03` need not be number 3. When you convert
    one region, Thyra also crops the slide photo to that region, plus a 10%
    margin.

    For other formats, `--region` is ignored with a warning.

## Leave out the optical image

For Bruker data, Thyra adds every image it finds in the data folder and the
folder above it. The image named in the flexImaging `.mis` file is the one
the MSI data is lined up with. To leave all images out, add `--no-optical`.

## Choose the m/z axis

### What Thyra does by default

Each spectrum is recorded at slightly different m/z values. To compare
pixels, Thyra puts every spectrum on one shared set of m/z positions, called
bins. Each column of the result then means the same m/z in every pixel. Most
analysis tools need this. The step is called **resampling**.

Thyra chooses the width of the bins, and the way peaks are moved onto them,
from your instrument. For most instruments the bins are 5 mDa wide at
m/z 1000.

??? advanced "Advanced: what Thyra chooses for each instrument"
    | Your data | Method | Axis type | Bin width at m/z 1000 |
    |---|---|---|---|
    | Bruker timsTOF | nearest_neighbor | reflector_tof | 5 mDa |
    | Bruker rapifleX | tic_preserving | constant | 5 mDa |
    | Bruker solariX, or imzML from an FT-ICR | nearest_neighbor | fticr | 5 mDa |
    | imzML from an Orbitrap | nearest_neighbor | orbitrap | 5 mDa |
    | PHI ToF-SIMS | nearest_neighbor | tof | about 57 mDa, 3 bins per peak |
    | Waters SELECT SERIES MRT (profile, the default) | tic_preserving | linear_tof | 1.3 mDa |
    | Waters SELECT SERIES MRT, centroid | nearest_neighbor | tof | 1.75 mDa |
    | Other Waters (centroid, the default) | nearest_neighbor | reflector_tof | 2 mDa |
    | Anything else, centroid | nearest_neighbor | reflector_tof | 5 mDa |
    | Anything else, profile | nearest_neighbor | constant | 5 mDa |

    The axis type sets how the bin width changes across the mass range.
    [Resampling](resampling.md) explains each type and each method.

### Keep the original m/z values

Add `--no-resample`. Each spectrum then keeps its own m/z values. Use it when
you want to do your own peak picking or calibration.

### Make the bins finer or coarser

Set the bin width at m/z 1000, in Da. A smaller width gives finer bins:

```bash
thyra input.imzML output.zarr --resample-width-at-mz 0.002
```

Finer bins make a larger result and need more memory. If an axis would need
more memory than your computer has free, Thyra refuses it before the
conversion starts.

??? advanced "Advanced: set the number of bins instead"
    `--resample-bins 100000` sets the number of bins directly. Use one of the
    two options, not both. `--resample-reference-mz` moves the m/z at which
    `--resample-width-at-mz` applies (1000 by default).

### Keep only part of the mass range

```bash
thyra input.imzML output.zarr --resample-min-mz 300 --resample-max-mz 1200
```

Both limits are included. Peaks outside the range are left out, and Thyra
warns you once.

??? advanced "Advanced: choose the method and the axis type yourself"
    `--resample-method` chooses how peaks are moved onto the bins:

    - `nearest_neighbor` puts each peak in the nearest bin. It is safe on
      every axis type.
    - `tic_preserving` spreads the signal between bins and keeps the total
      ion current. It is exact only on an axis spaced like your instrument's
      own sampling.

    `--mass-axis-type` chooses how the bin width changes across the mass
    range. If you change it and leave the method on automatic, Thyra
    switches to `nearest_neighbor` where it has to. If you ask for
    `tic_preserving` yourself, Thyra does what you ask. Read
    [Resampling](resampling.md) before combining the two options.

## Keep 3D data together

A dataset with several slices becomes one result part per slice, named
`_z0`, `_z1` and so on. To keep the slices together as one 3D volume, add
`--handle-3d`. Give the distance between slices, in micrometres, with
`--z-spacing`:

```bash
thyra volume.imzML volume.zarr --handle-3d --z-spacing 20
```

Without `--z-spacing`, Thyra assumes the slices are one pixel apart, and warns
you.

## Name the dataset

Every part of the result is named after the dataset, `msi_dataset` unless you
choose another name:

```bash
thyra slide.d slide.zarr --dataset-id mouse_brain
```

Use letters, digits, `_`, `.` and `-`.

## Keep a record of the conversion

To save everything Thyra reports into a file, add `--log-file`:

```bash
thyra slide.d slide.zarr --log-file conversion.log
```

Add `-v DEBUG` for more detail. The file is added to, not replaced, on every
run.

??? advanced "Advanced: exit status, for scripts"
    `thyra` exits with status 0 when the result was written, 1 when the
    conversion failed or was refused, and 2 when the command itself is
    wrong, for example an option that does not exist or a result that
    already exists.

## Options for one kind of instrument

??? advanced "Advanced: options for imzML, Bruker and Waters data"
    - `--spectrum-type profile` or `centroid` (imzML): use it when the file
      declares the wrong type of spectrum. It changes which m/z axis Thyra
      chooses.
    - `--no-recalibrated` (Bruker timsTOF): use the calibration recorded
      during the acquisition, instead of the newest recalibration.
    - `--tdf-spectrum vendor_centroid` (Bruker timsTOF with ion mobility):
      use Bruker's own peak-picked spectrum instead of the sum of all scans.
    - `--mobility-grid` (Bruker timsTOF with ion mobility): add a table that
      keeps the ion mobility of every pixel.
    - `--no-msms-table` (Bruker timsTOF PASEF MS/MS): leave out the table
      that splits the fragments by precursor.
    - `--waters-spectrum profile` or `centroid` (Waters): profile keeps close
      masses apart, at about three times the size. It is the default on a
      SELECT SERIES MRT; other Waters instruments default to centroid.
    - `--intensity-threshold` (any data): leave out peaks below this
      intensity. They cannot be recovered from the result.

## Convert from Python

??? advanced "Advanced: the same settings in Python"
    `convert_msi` takes the same settings as the command:

    ```python
    from thyra import convert_msi

    convert_msi(
        "slide.d",
        "slide.zarr",
        pixel_size_um=20,
        resampling_config={},
    )
    ```

    From Python, resampling is off unless you pass `resampling_config`.
    An empty `{}` gives the same automatic choice as the command. The
    [Python API](api.md) page lists every argument. `convert_msi` returns
    `True` when the result was written and `False` when it was not.

## Go deeper

- [Command-line options](cli.md): every option, with examples.
- [Resampling](resampling.md): the shared m/z axis in full detail.
