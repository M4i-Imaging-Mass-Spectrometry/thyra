# Resampling

Mass axis resampling maps every spectrum in a dataset onto one shared m/z axis.
It is **enabled by default**, and by default every parameter is chosen
automatically from the acquisition metadata.

This page explains what that automatic choice is, how it is made, and how to
override it. For the list of flags alone, see the
[CLI Reference](cli.md#resampling-advanced).

---

## Why resampling exists

An imzML file comes in one of two storage modes:

| Mode | m/z arrays | Consequence |
|---|---|---|
| **continuous** | one shared array for the whole dataset | pixels are already comparable bin-for-bin |
| **processed** | a separate array per spectrum | pixel 1's column 500 and pixel 2's column 500 are *different masses* |

For processed data -- which includes most centroided vendor exports -- you
cannot build a pixel x m/z matrix at all without first agreeing on a common
axis. Thyra reports this as `needs_resampling`, and it is the reason resampling
defaults to on rather than off.

Even for continuous data, resampling is often still what you want: it lets you
impose a physically sensible bin spacing, restrict the mass range, or match the
axis convention of another tool.

To keep the original per-spectrum axes untouched, disable it:

```bash
thyra input.imzML output.zarr --no-resample
```

Note that with `--no-resample` on processed data, each pixel keeps its own m/z
values, and most downstream tools will not be able to treat the table as a
matrix of comparable features.

---

## How the automatic choice is made

Three things are decided: the **method** (how intensity moves onto the new
axis), the **axis type** (how bin widths scale with m/z), and the **bin count**.

```
  acquisition metadata
          │
          ▼
  DataCharacteristics        storage mode, spectrum type, instrument
          │                  name/type/manufacturer, peak density
          ▼
  InstrumentDetectorChain    first matching detector wins
          │
          ├──────────────► resampling method   (nearest_neighbor | tic_preserving)
          ├──────────────► axis type           (constant | *_tof | orbitrap | fticr)
          └──────────────► bin count           from a target width at a reference m/z
```

### What Thyra looks at

`DataCharacteristics` collects, from the metadata:

- whether there is a shared mass axis (continuous) or not (processed)
- spectrum type -- `centroid spectrum` (`MS:1000127`) or `profile spectrum` (`MS:1000128`)
- instrument name, instrument type, and manufacturer
- average peaks per spectrum; above **5000** the data is treated as high-density profile
- format flags for Rapiflex and timsTOF

The timsTOF flag is a case-insensitive substring match on the instrument name,
because Bruker names the family many ways (`timsTOF fleX MALDI-2`,
`timsTOF Pro 2`, `timsTOF SCP`, and so on).

### Which detector wins

Detectors are tried in a fixed priority order and the first match wins:
timsTOF, Rapiflex, FT-ICR, Orbitrap, generic centroid, then a catch-all
default. This table is the actual observed behaviour of that chain:

| Metadata | Detector | Method | Axis type |
|---|---|---|---|
| timsTOF, centroid | timsTOF | `nearest_neighbor` | `reflector_tof` |
| timsTOF, profile high-density | timsTOF | `nearest_neighbor` | `reflector_tof` |
| Rapiflex, profile | Rapiflex MALDI-TOF | `tic_preserving` | `constant` |
| Bruker MALDI-TOF | Rapiflex MALDI-TOF | `tic_preserving` | `constant` |
| unknown vendor, profile high-density | Rapiflex MALDI-TOF | `tic_preserving` | `constant` |
| unknown, centroid | ImzML Centroid | `nearest_neighbor` | `reflector_tof` |
| no usable metadata | Unknown (default) | `nearest_neighbor` | `constant` |

!!! warning "The Orbitrap and FT-ICR detectors cannot currently fire"
    Both match on an `instrument_type` of `"Orbitrap"` or `"FT-ICR"`, but no
    reader populates that key: the only place it is set is the Rapiflex reader,
    which hardcodes `"MALDI-TOF"`. Neither the imzML nor the Bruker `.d`
    metadata extractor emits it at all.

    So Orbitrap and FT-ICR data does **not** get the `orbitrap` or `fticr`
    axis automatically -- it falls through to the generic centroid or default
    row above and gets `reflector_tof` or `constant`. Set
    `--mass-axis-type orbitrap` or `--mass-axis-type fticr` explicitly if you
    want the matching physics.

!!! note "`tic_preserving` is selected for profile MALDI-TOF, not for high resolution"
    It is easy to assume the high-resolution analysers get the more elaborate
    method. They do not. `tic_preserving` is chosen for **profile** data, where
    a peak is spread over many points and rebinning would otherwise change the
    total ion count. Orbitrap and FT-ICR data is normally centroided, so it gets
    `nearest_neighbor`, which is the correct choice for discrete peaks. If you
    have *profile* Orbitrap or FT-ICR data, set
    `--resample-method tic_preserving` yourself.

!!! danger "Do not combine `tic_preserving` with a non-uniform axis type"
    `tic_preserving` interpolates onto the target axis and then applies a
    single scaling factor to the whole spectrum. A single factor cannot
    account for bin widths that vary across the mass range, so pairing it with
    `linear_tof`, `reflector_tof`, `orbitrap` or `fticr` suppresses high-m/z
    ions relative to low-m/z ones. Measured across 300-1100 m/z, two ions of
    equal abundance come back with their ratio distorted by 1.9x on
    `linear_tof`, 3.7x on `reflector_tof`, 7.0x on `orbitrap` and 13.4x on
    `fticr`.

    Auto-selection never produces these pairings -- `tic_preserving` is only
    ever chosen alongside `constant`, which is uniform and therefore exact.
    You have to ask for the combination with two explicit flags.

    If you want a non-uniform axis, use `nearest_neighbor`, which moves each
    peak into a single bin and is unaffected by bin width.

The chosen detector, method, and axis type are all logged:

```
INFO - Detected instrument type: timsTOF
INFO - Selected resampling method: NEAREST_NEIGHBOR
INFO - Selected axis type: REFLECTOR_TOF
```

---

## Methods

### `nearest_neighbor`

Each target bin takes the intensity of the nearest original m/z value. Peaks
stay sharp and stay put, and no intensity is invented between them. This is the
right choice for **centroid** data, where a peak is a single discrete mass and
smearing it across neighbouring bins would be wrong.

Consequence worth knowing: on a target axis finer than the source spacing, most
bins receive nothing. A dataset resampled from 4,000 source points onto 190,000
bins has exactly 4,000 populated bins per spectrum, and the rest are zero. This
is why you should extract ion images by **summing over an m/z window** rather
than picking the single nearest bin -- see
[the tutorial](tutorial.md#step-7-ion-images).

### `tic_preserving`

Linear interpolation onto the target axis, followed by rescaling so the
spectrum's total ion current matches the original. This is the right choice for
**profile** data: a peak spans many points, interpolation reconstructs its shape
on the new grid, and the rescaling step stops rebinning from quietly changing
quantitation.

Use it whenever the total ion count per pixel has to stay comparable before and
after conversion.

The total that is preserved is the share of the spectrum inside the axis range.
If you crop with `--resample-min-mz` or `--resample-max-mz`, intensity outside
the window is dropped rather than redistributed over the bins you kept, so the
per-pixel total falls by whatever you cropped out -- a cropped window should not
claim ions from the parts that were cut away. That share is measured by area, so
a window narrower than the spacing between source points still keeps a
proportionate amount instead of collapsing to zero.

With the default axis, which spans the dataset's own mass range, nothing is
dropped and the total is preserved exactly.

---

## Axis types

The axis type sets how bin width grows with m/z. Each corresponds to the
physics of a mass analyser, so that bins track the instrument's real resolving
power instead of over-sampling the low end and under-sampling the high end.

| Axis type | Bin width scales as | Rationale |
|---|---|---|
| `constant` | constant Da | Equidistant bins; the SCiLS Lab convention for MALDI-TOF profile data |
| `linear_tof` | `sqrt(m/z)` | Linear TOF: flight time `t ∝ sqrt(m/z)`, so equal time bins give `sqrt(m/z)` mass bins |
| `reflector_tof` | `m/z` | Constant *relative* resolution `R = m/Δm`; bins are uniform in `ln(m/z)` |
| `orbitrap` | `m/z^1.5` | Orbitrap frequency `f ∝ 1/sqrt(m/z)`, so equal frequency bins give `m/z^1.5` mass bins |
| `fticr` | `m/z^2` | Cyclotron frequency `f ∝ 1/(m/z)`, so equal frequency bins give `m/z^2` mass bins |

`reflector_tof` is the most broadly useful of these: constant relative
resolution means constant relative mass accuracy across the whole range, which
is what most MS workflows assume.

---

## Bin count

You can set the bin count directly, or specify a target bin **width at a
reference m/z** and let Thyra derive the count from the axis physics. The two
are mutually exclusive.

When you specify neither, these defaults apply:

| Axis type | Default width | At reference m/z |
|---|---|---|
| `linear_tof` | 17 mDa | 300 |
| everything else | 5 mDa | 1000 |

The 17 mDa / m/z 300 pairing for `linear_tof` reproduces the axis SCiLS Lab
produces for FlexImaging data. The count is then derived per axis type:

| Axis type | Bin count |
|---|---|
| `reflector_tof` | `ln(max/min) * (ref_mz / width)` |
| `linear_tof` | `(2/k) * (sqrt(max) - sqrt(min))`, `k = width / sqrt(ref_mz)` |
| `orbitrap` | `2 * (1/sqrt(min) - 1/sqrt(max)) * (ref_mz^1.5 / width)` |
| `fticr` | `(1/min - 1/max) * (ref_mz^2 / width)` |
| `constant` | `(max - min) / width` |

Each is the integral of `1 / width(m)` across the mass range, so the bin width
the axis actually realizes at your reference m/z is the width you asked for.
The result is floored at 100 bins.

A worked example, from the [tutorial](tutorial.md)'s synthetic dataset --
250-1200 Da, `constant` axis, default 5 mDa width:

```
(1200 - 250) / 0.005 = 190000 bins
```

which is what the log reports:

```
INFO - Calculating bins for 5.0 mDa width at m/z 1000.0
INFO - Calculated 190000 bins for constant axis type
INFO - Building resampled mass axis: 250.00 - 1200.00 m/z, 190000 bins
```

190,000 bins from 4,000 source points sounds extravagant, but the matrix is
stored sparsely, so the store stays small -- roughly 29 MB here. Pass
`--resample-bins` or `--no-resample` if you would rather not upsample.

### The 10 million bin cap

`--no-resample` on a processed-mode imzML does not build a physics axis at all.
It builds the **raw** axis: every distinct m/z value in the whole dataset, one
column each. When the peak lists genuinely share m/z values that is small. When
they do not -- which is the normal case for centroided peak picking, where two
pixels almost never report a peak at the same double -- it grows to roughly one
column per peak in the entire dataset, and there is no natural stopping point.

Thyra gives up once that axis passes **10 million** unique m/z values, and says
what to do instead:

```
ValueError: Common mass axis exceeded 10,000,000 unique m/z values after
12,400 of 918,855 spectra (10,004,112 so far). The peak lists in this dataset
do not share m/z values, so a raw axis grows to roughly one column per peak,
which is not usable downstream. Convert with resampling instead (it is the
default; --no-resample disables it), or raise max_mass_axis_length.
```

10 million is SCiLS Lab's own limit on the same quantity: "Data sets in SCiLS
Lab are limited to a maximum of 10 million bins on the common mass axis" (2026b
User Guide, p.76).

The cap is Python-API only, through `reader_options`, and applies to the imzML
reader:

```python
from thyra import convert_msi

convert_msi(
    "data.imzML",
    "out.zarr",
    reader_options={"max_mass_axis_length": 50_000_000},  # or None for no cap
)
```

It has no effect when resampling is on, which is the default: a resampled axis
has the bin count you asked for.

---

## Overriding the automatic choice

Every part is independently overridable; anything you leave alone stays
automatic.

=== "CLI"

    ```bash
    # Method and axis type. nearest_neighbor is the pairing to use with a
    # non-uniform axis -- see the warning under "Selection" above.
    thyra input.imzML output.zarr \
        --resample-method nearest_neighbor \
        --mass-axis-type orbitrap

    # Fixed bin count
    thyra input.imzML output.zarr --resample-bins 50000

    # Target width at a reference m/z instead of a bin count
    thyra input.imzML output.zarr \
        --resample-width-at-mz 0.01 \
        --resample-reference-mz 500

    # Restrict the mass range
    thyra input.imzML output.zarr \
        --resample-min-mz 400 --resample-max-mz 1000

    # Off entirely
    thyra input.imzML output.zarr --no-resample
    ```

=== "Python"

    ```python
    from thyra import convert_msi

    convert_msi(
        "input.imzML",
        "output.zarr",
        resampling_config={
            "method": "tic_preserving",
            "axis_type": "orbitrap",
            "width_at_mz": 0.01,
            "reference_mz": 500.0,
        },
    )
    ```

!!! warning "The Python API does not resample by default"
    `convert_msi()` called without `resampling_config` keeps the original mass
    axis, whereas the `thyra` command-line tool resamples by default. Pass a
    `resampling_config` explicitly if you want the CLI's behaviour from Python.

---

## Checking what was used

The decision is recorded in the output store, so a converted dataset is
self-describing:

```python
import spatialdata as sd

sdata = sd.read_zarr("output.zarr")
table = sdata.tables["msi_dataset_z0"]

mz = table.var["mz"].values
print(f"{len(mz)} bins, {mz.min():.2f} -- {mz.max():.2f} m/z")

# Bin width across the range tells you the axis type in practice:
# flat -> constant, growing linearly -> reflector_tof, and so on.
import numpy as np
d = np.diff(mz)
print(f"bin width: {d.min()*1000:.3f} -- {d.max()*1000:.3f} mDa")
```

For the tutorial dataset this prints a flat 5.000 mDa across the range,
confirming a `constant` axis at the default width.

---

## Practical guidance

- **Leave it alone** unless you have a reason. The detected settings are right for the common Bruker and imzML cases.
- **Disable it** (`--no-resample`) when you want the untouched vendor axis -- your own peak picking, centroiding, or calibration downstream.
- **Set `tic_preserving`** if quantitation matters and your data is profile, especially profile data from a high-resolution analyser, which the automatic choice will not pick it for.
- **Set `--resample-bins`** any time you need to match another tool's axis exactly.
- **Narrow the mass range** (`--resample-min-mz` / `--resample-max-mz`) before increasing bin counts; it is usually the cheaper way to get the resolution you need where you need it.

---

## See also

- **[CLI Reference](cli.md#resampling-advanced)** -- the flags
- **[Tutorial](tutorial.md)** -- resampling in the context of a full conversion
- **[Output Format](output-format.md)** -- where the m/z axis lives in the store
- **[API Reference](api.md)** -- `ResamplingConfig`, `ResamplingDecisionTree`
