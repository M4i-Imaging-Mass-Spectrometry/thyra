# Describe and share your data

Every result carries a description of the acquisition: the instrument, the
settings, the pixel size, and what Thyra did to the data. This page shows how
to read that description, check it, and turn it into the description
METASPACE asks for.

## What Thyra records for you

Thyra fills in whatever your files state:

- the pixel size, and whether it came from the file or from `--pixel-size`
- the instrument: polarity, ion source, analyser, model, maker and serial
  number
- the acquisition settings, such as the start time, laser power, repetition
  rate, shots per pixel and method file
- the calibration, for Bruker timsTOF, PHI and Waters data
- what Thyra did: its version, the number of spectra, and the processing
  steps, such as building the shared m/z axis

How much of this a file states depends on the instrument. From imzML files,
for example, Thyra reads no acquisition settings.

Thyra cannot know about your sample: the organism, the organ or tissue, the
condition, and the matrix and how it was applied. You add those when you
prepare a submission (see below).

!!! note "Before you share a result"
    Thyra leaves out people's names and e-mail addresses, and shortens
    folder paths recorded by the instrument to the file name. It keeps two
    things you may want to check: the location of the source file on your
    computer, and any free text typed into the instrument software, such as
    sample or method names.

## See the description of a dataset

To see the description without converting anything, point
`thyra metadata` at your data:

```bash
thyra metadata my_slide.d
```

It prints the description on the screen. To save it to a file, add
`-o my_slide.json`.

## Check a result

```bash
thyra validate brain_section_1.zarr
```

For a result in good order, it prints one line per table:

```
msi_dataset_z0: OK
```

A warning is printed under that line and does not make the check fail. An
error means the description does not follow the rules, and the line says
`FAILED`.

??? advanced "Advanced: what `thyra validate` checks"
    It checks that the description follows the published rules for its
    version, that every value is of the right kind, and that the codes from
    the mass spectrometry vocabularies match their names. For a result, it
    also checks that every table carries a description and that its m/z
    values are numbers in increasing order. It never reads the intensities.
    It exits with status 0 when everything conforms, 1 when something does
    not, and 2 when the command itself is wrong. Add `--json` for output a
    program can read.

## Prepare a METASPACE submission

[METASPACE](https://metaspace2020.org) annotates the molecules in imaging mass
spectrometry data. Its upload form asks for a description of the sample and
the measurement. Thyra fills in what your files state, and you add the rest.

**Step 1.** Write the details only you know into a file, for example
`sample.json`:

```json
{
  "sample": {
    "organism": "Mus musculus",
    "organism_part": "brain",
    "condition": "wildtype"
  },
  "preparation": {
    "matrix": "2,5-dihydroxybenzoic acid (DHB)",
    "matrix_application": "automated sprayer"
  },
  "ms_analysis": {
    "detector_resolving_power": {"value": 30000, "at_mz": 400}
  }
}
```

**Step 2.** Combine it with the description in your result:

```bash
thyra export-metaspace brain_section_1.zarr --merge sample.json
```

This writes `brain_section_1.metaspace.json` next to the result. If a field
METASPACE requires is still empty, Thyra names it, for example
`Sample_Information.Organism`. Add it to `sample.json` and run the command
again.

**Step 3.** Upload to METASPACE. METASPACE needs the imzML files themselves,
not the `.zarr` result. Upload your `.imzML` and `.ibd` files at
[metaspace2020.org](https://metaspace2020.org), and fill in the form with the
values from `brain_section_1.metaspace.json`. If your data is not imzML yet,
export imzML from your instrument software first.

Your `sample.json` is only combined for the export. The result itself does
not change, so keep `sample.json` next to your data.

??? advanced "Advanced: more fields, several tables, and submitting from Python"
    - **Fields your files do not state.** If Thyra says polarity, ion source
      or analyser is missing, add them under `ms_analysis`, for example
      `"polarity": "positive"`, `"ionisation_source": "MALDI"` and
      `"analyzer": "TOF"`. Write the polarity in lower case. A misspelt key
      is refused with an error that names it.
    - **Several tables.** A result with several slices, or with an extra
      table for ion mobility or MS/MS, needs `--table` to say which table to
      export, for example `--table msi_dataset_z0`.
    - **Submitting from Python.** The `metaspace2020` Python client can
      submit a dataset with this file as its metadata. It needs a METASPACE
      account and an API key.

??? advanced "Advanced: what \"versioned\" and \"ontology-mapped\" mean"
    Each description names the edition of the rules it follows, and a
    published edition is never changed afterwards. Key values carry a
    standard code next to the plain word: MALDI is `MS:1000075`, and mouse
    is `NCBITaxon:10090`. Software can then read them without guessing how
    they were spelt.

## Go deeper

- [Metadata schema](metadata-schema.md): every field, and how each format
  fills it.
- [Writing the metadata document](writing-the-metadata-document.md): for
  programs that want to write the same description.
