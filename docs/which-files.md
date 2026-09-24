# Which files can I convert?

Thyra reads the files your instrument writes. Point it at the file or folder
in the table, and it recognises the format by itself.

| Your instrument | Point Thyra at | Windows | Linux | macOS |
|---|---|---|---|---|
| Any instrument that exports imzML | the `.imzML` file | yes | yes | yes |
| Bruker timsTOF | the folder ending in `.d` | yes | yes | no |
| Bruker solariX | the folder ending in `.d` | yes | yes | yes |
| Bruker rapifleX | the acquisition folder | yes | yes | yes |
| Waters (MassLynx) | the folder ending in `.raw` | yes | yes | no |
| PHI nanoTOF (ToF-SIMS) | the file ending in `.raw` | yes | yes | yes |
| mzPeak (experimental) | the `.mzpeak` file | yes | yes | yes |

**Your instrument is not listed?** Most instrument software can export imzML,
and Thyra converts that. For Shimadzu data, export imzML from IMAGEREVEAL MS;
native support for Shimadzu files is in development.

## Check a dataset before you convert it

To see what Thyra finds in a dataset without converting it, run:

```bash
thyra metadata my_slide.d
```

It prints what Thyra can tell about the acquisition: the instrument, the
polarity, the pixel size and more. If no pixel size appears, the conversion
will need `--pixel-size`.

??? advanced "Advanced: What `thyra metadata` prints"
    The output is a JSON document in the same form Thyra stores inside every
    result. Add `-o description.json` to save it to a file instead. For Waters
    data, Thyra reads every spectrum to describe the run, so it can take a
    while. [Describe and share your data](describe-your-data.md) explains the
    document.

## imzML

- Keep the `.ibd` file next to the `.imzML` file, with the same name. Thyra
  needs both.
- Thyra reads the pixel size from the file when the file records it. If it
  does not, give it with `--pixel-size`.
- A file with several z-slices becomes one dataset per slice. Add
  `--handle-3d` to keep them together as one 3D volume.
- If Thyra refuses the file, see
  [Troubleshooting](troubleshooting.md#an-imzml-file-is-refused-before-the-conversion-starts).

## Bruker timsTOF

- Point Thyra at the `.d` folder. You can also point it at the folder above:
  if that holds several `.d` folders, Thyra asks which one to convert.
- Keep the `.mis` file from flexImaging next to the `.d` folder, ideally with
  the same name. It gives the pixel size and lines the data up with the
  optical image.
- Thyra adds every image (TIFF, JPEG, PNG or BMP) it finds in the `.d`
  folder, in the folder you pointed at, and in the folder above the `.d`.
  Unrelated pictures in those folders end up in the result too. To leave all
  images out, add `--no-optical`.
- A slide with several regions is converted as one dataset, with every pixel
  labelled with its region. To convert only one region, see
  [Change how Thyra converts](settings.md#convert-one-region).
- Ion mobility (TIMS) and PASEF MS/MS data are kept automatically.
- Close Bruker DataAnalysis before converting. While it has the dataset open,
  Thyra cannot read it.

??? advanced "Advanced: ion mobility and MS/MS options"
    For TIMS data, Thyra sums each spectrum over the mobility dimension and
    also stores a mass-mobility heatmap. It can add a table that keeps
    mobility per pixel (`--mobility-grid`). PASEF MS/MS acquisitions get a
    separate table with the fragments split by precursor. The
    [command-line options](cli.md) page describes these options.

## Bruker solariX

- Point Thyra at the `.d` folder itself. It must hold the files
  `peaks.sqlite` and `ImagingInfo.xml`.
- The pixel size comes from a `.mis` file next to the `.d` folder. Without
  one, give it with `--pixel-size`.
- A `.d` folder that holds only raw transients, with no `peaks.sqlite`, cannot
  be read. Export imzML from DataAnalysis, SCiLS Lab or flexImaging, and
  convert that instead.
- Optical images are added, but they are not lined up with the MSI data.

## Bruker rapifleX

- Point Thyra at the acquisition folder, not at the `.d` folder inside it.
  The folder must hold the `.dat`, `_poslog.txt` and `_info.txt` files.
- The pixel size comes from the `_info.txt` file.
- A `.mis` file in the folder lines the data up with the optical image.

## Waters (MassLynx)

- Point Thyra at the `.raw` folder, and copy the folder whole: if any
  `_FUNC` file is missing, Thyra stops.
- The pixel size comes from the recorded laser positions.
- Ion mobility functions are skipped. When a run records low and high
  collision energy over the same pixels (MSe), Thyra keeps the low-energy
  data.

??? advanced "Advanced: profile and centroid data from Waters"
    Waters runs can be read as profile spectra or as centroids. For SELECT
    SERIES MRT data, Thyra reads the profile by default, which gives a larger
    result. The [command-line options](cli.md) page shows how to choose.

## PHI nanoTOF (ToF-SIMS)

- Point Thyra at the `.raw` file. Unlike Waters data, it is a single file.
- The pixel size comes from the file.
- A mosaic of tiles is joined into one image, and the frames of a depth
  profile are added into one 2D image.
- Mosaic, MS/MS and depth-profile acquisitions have so far been tested only
  on synthetic files. If you have real data of these kinds, please
  [get in touch](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/issues).

## mzPeak (experimental)

- Point Thyra at the `.mzpeak` file. Only 2D images are read.
- Thyra reads the pixel size only when the file records it.

## Go deeper

[Supported formats](supported-formats.md) explains how each format is
recognised and read, in full detail.
