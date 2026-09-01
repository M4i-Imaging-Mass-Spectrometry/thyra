# solariX Notes

Thyra reads Bruker solariX / MRMS (FT-ICR) imaging `.d` directories natively
by consuming `peaks.sqlite`, the processed peak store ftmsControl writes into
every imaging acquisition. There is no public specification for this file;
the layout below was established by inspecting real acquisitions spanning
2018-2024 (ftmsControl 2.2, flexImaging 5.0.80-5.0.89) and decode-verifying
the binary blobs against their own metadata. The same file is read by at
least one independent open-source tool for the solariX XR, so the store is
not a single-instrument quirk.

What this route deliberately is **not**: a transient processor. The `.d`
also holds the raw time-domain transients (`ser`, tens of GB per region) and
Bruker's proprietary `.mcf` containers. Thyra never opens either -- no FT,
no apodization, no peak picking. The peaks were picked by the acquisition
software at acquisition time, and that processing is frozen in the data.
Anyone needing different peak picking must export via the vendor software.

---

## Directory layout

```
sample.d/
  peaks.sqlite            centroided per-pixel peak lists  <- what Thyra reads
  ImagingInfo.xml         per-scan index (spot name, TIC, base peak)
  ser                     raw transients, fixed-length per scan (ignored)
  <guid>_N.mcf/.mcf_idx   Bruker container stores (ignored)
  parameterChanges.sqlite per-scan parameter history (ignored)
  Storage.mcf_idx         container index (ignored)
  *.m/                    acquisition method directory (ignored)
sample.mis                flexImaging sequence NEXT TO the .d, same stem
sample.dat                flexImaging cache (ignored)
```

Detection requires `peaks.sqlite` **and** `ImagingInfo.xml` inside a `.d`
directory. Two traps are handled explicitly:

- **Pre-scan directories.** flexImaging writes sibling `.d` directories such
  as `sample_0_C17_000001.d` holding `fid` + `analysis.baf` + `calib.bin` --
  single-spectrum pre-scans, not imaging runs. They contain none of the
  detection files and are rejected as unrecognised Bruker data.
- **Peaks-less acquisitions.** A `.d` with `ser` and `ImagingInfo.xml` but no
  `peaks.sqlite` is still recognisably solariX-family; the error names the
  imzML export fallback (DataAnalysis, SCiLS Lab, or flexImaging) instead of
  claiming the directory is not Bruker data.

## peaks.sqlite

Three tables matter. `Properties` is a key/value table carrying the schema
identity (`SchemaType = Imaging`, `SchemaVersionMajor = 1`,
`SchemaVersionMinor = 2` in every file inspected), the acquisition software
name/version, the instrument vendor, family, source type and serial number,
the operator, the acquisition timestamp and the acquired m/z range
(`MzAcqRangeLower`/`MzAcqRangeUpper`).

The reader checks `SchemaType` and `SchemaVersionMajor` at open and refuses
anything it has not been verified against, quoting the values it found. A
newer ftmsControl that bumps the major version fails loudly rather than
decoding blobs on an unknown layout.

`AcquisitionKeys` holds one row per acquisition mode with `Polarity`,
`ScanMode`, `AcquisitionMode` and `MsLevel`. Only the polarity enum has been
verified (against a paired positive/negative acquisition of the same
sample): **0 = positive, 1 = negative**. The other enums are undocumented,
so Thyra stores the raw integers without inventing labels. Observed values:
`ScanMode 0 / MsLevel 1` on full scans, `ScanMode 2 / MsLevel 2` on a CASI
(quadrupole isolation) run.

`Spectra` holds one row per pixel: stage-raster coordinates
(`XIndexPos`/`YIndexPos`), `NumPeaks`, per-scan acquisition values
(`NumSummations`, `LaserPower`, `LaserRepRate`), `RegionNumber`, and the
peak blobs.

### Blob encoding

| Column | Encoding | Length |
|---|---|---|
| `PeakMzValues` | little-endian float64, ascending | `NumPeaks * 8` |
| `PeakIntensityValues` | little-endian float32 | `NumPeaks * 4` |
| `PeakFwhmValues` | little-endian float32 | `NumPeaks * 4` |
| `PeakSnrValues` | little-endian float32 | `NumPeaks * 4` |

The m/z values are already calibrated -- the frequency-to-m/z law (ML1/ML2/
ML3 in the acquisition method) has been applied by the time the peaks are
stored, so the reader does no calibration math. Every blob length is checked
against `NumPeaks * itemsize` before decoding; a mismatch is refused as
corruption rather than silently truncated. `PeakFlags`, `BH` and `BC` are
undecoded (flags were all-zero in every file inspected) and nothing is
derived from them.

The database is opened strictly read-only through a URI (`mode=ro`), so
sqlite never creates journal or WAL side files next to it. That matters
because solariX data commonly lives on read-only network shares.

## Coordinates and sparsity

`XIndexPos`/`YIndexPos` are **absolute** stage-raster indices -- a real
region started at X 1239, not 0. The reader subtracts the per-dataset
minimum to produce 0-based coordinates and reports the subtracted origin in
`coordinate_offsets` and the absolute extents in the format metadata.

Rasters are sparse: flexImaging acquires the polygon drawn on the optical
image, so the bounding box contains cells that were never measured (81.5 %
fill on the reference region). Unmeasured cells are simply absent from the
output, exactly like sparse pixels in any other format. A measured scan can
also record zero peaks (seen at the tissue-free edges of one acquisition);
such pixels are skipped like unmeasured ones.

One fidelity check worth recording: on the reference region the sum of a
pixel's stored peak intensities reproduces the per-scan `tic` value in
`ImagingInfo.xml` -- an independent record written at acquisition time --
to a ratio of 0.9999. TIC images computed from the peak store are therefore
faithful to what the instrument reported.

## Pixel size

The `.d` itself does not record the raster pitch anywhere -- the acquisition
method file was swept for raster/pixel tokens and carries none. The only
source is the `<Raster>x,y</Raster>` element (micrometres) of the
flexImaging `.mis` file sitting next to the `.d` with the same stem.

Resolution order: a stem-matching `.mis` in the parent directory; failing
that, a *lone* `.mis` in the parent is accepted as unambiguous. Multiple
non-matching candidates are refused -- on a multi-region slide each region
has its own `.mis`, and guessing one would silently assign a wrong pixel
size. When no source resolves, the pixel size is reported as unknown and
conversion asks for `--pixel-size` instead of defaulting.

## Resampling

The stored spectra are centroided per-pixel peak lists, so every spectrum
has its own m/z values. The metadata extractor stamps
`instrument_type = "FT-ICR"` from the instrument identity in `Properties`
(`InstrumentFamily 513` additionally maps to the model name `solariX`;
unknown family numbers are surfaced raw). The resampling decision tree
therefore lands on `nearest_neighbor` onto the `fticr` axis with no flags --
the same pair a vendor imzML export of the same data reaches, minus the
export step. See [Resampling](resampling.md) for why TIC-preserving is
refused for centroid data.

## What the native route adds over an imzML export

Everything an export carries, plus metadata no exporter writes: per-pixel
laser power, repetition rate and summation counts, the acquisition software
identity and version, the raw acquisition-mode enums, per-scan region
numbers, the absolute stage-raster extents, and the acquisition m/z range
as acquired. And no manual vendor-software step per dataset.
