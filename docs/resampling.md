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
| Orbitrap, centroid | Orbitrap | `nearest_neighbor` | `orbitrap` |
| FT-ICR, centroid | FT-ICR | `nearest_neighbor` | `fticr` |
| unknown vendor, profile high-density | Rapiflex MALDI-TOF | `tic_preserving` | `constant` |
| unknown, centroid | ImzML Centroid | `nearest_neighbor` | `reflector_tof` |
| no usable metadata | Unknown (default) | `nearest_neighbor` | `constant` |

!!! note "`tic_preserving` is selected for profile MALDI-TOF, not for high resolution"
    It is easy to assume the high-resolution analysers get the more elaborate
    method. They do not. `tic_preserving` is chosen for **profile** data, where
    a peak is spread over many points and rebinning would otherwise change the
    total ion count. Orbitrap and FT-ICR data is normally centroided, so it gets
    `nearest_neighbor`, which is the correct choice for discrete peaks. If you
    have *profile* Orbitrap or FT-ICR data, set
    `--resample-method tic_preserving` yourself.

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

---

## Overriding the automatic choice

Every part is independently overridable; anything you leave alone stays
automatic.

=== "CLI"

    ```bash
    # Method and axis type
    thyra input.imzML output.zarr \
        --resample-method tic_preserving \
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
