# Metadata Schema

Every store Thyra writes carries a versioned, ontology-mapped metadata
block: `table.uns["msi_metadata"]`. Its base fields mirror the
[METASPACE](https://metaspace2020.org) submission form, so the metadata a
converted dataset carries is, by construction, what a METASPACE submission
needs -- filling it costs nothing extra.

The schema exists because MSI has had no structured metadata convention
the way other spatial omics modalities do. Vendor files spell the same
fact many ways, imzML metadata is free-form cvParams, and every pipeline
re-invents its own dictionary. This block fixes the names, maps them to
ontologies, and ships a validator.

```python
import spatialdata as sd

sdata = sd.read_zarr("output.zarr")
block = sdata.tables["msi_dataset_z0"].uns["msi_metadata"]

print(block["schema_version"])                       # "0.1.0"
print(block["ms_analysis"]["pixel_size_um"])         # {"x": 20.0, "y": 20.0}
print(block["ms_analysis"]["ionisation_source"])     # "MALDI"
print(block["ms_analysis"]["ionisation_source_term"])
# {"accession": "MS:1000075",
#  "name": "matrix-assisted laser desorption ionization"}
```

---

## The document

Four sections. `ms_analysis` and `provenance` are written by the converter
for every store; `sample` and `preparation` describe things no raw file
records (what the tissue was, how it was prepared) and are supplied by
you -- see [Completing the metadata](#completing-the-metadata).

| Section | Field | Type | Ontology |
|---------|-------|------|----------|
| (root) | `schema_version` | `MAJOR.MINOR.PATCH` string, required | -- |
| `sample` | `organism`, `organism_term` | text + term | NCBITaxon |
| | `organism_part`, `organism_part_term` | text + term | UBERON |
| | `condition` | text | -- |
| | `sample_growth_conditions` | text | -- |
| `preparation` | `sample_stabilisation` | text | -- |
| | `tissue_modification` | text | -- |
| | `matrix`, `matrix_term` | text + term | CHEBI |
| | `matrix_application` | text | -- |
| | `solvent` | text | -- |
| `ms_analysis` | `polarity`, `polarity_term` | `"positive"`/`"negative"` + term | PSI-MS |
| | `ionisation_source`, `ionisation_source_term` | text + term | PSI-MS |
| | `analyzer`, `analyzer_term` | text + term | PSI-MS |
| | `instrument_model` | text | -- |
| | `detector_resolving_power` | `{value, at_mz}` | -- |
| | `pixel_size_um` | `{x, y}`, **required** | -- |
| `provenance` | `thyra_version` | text, required | -- |
| | `source_format` | `"imzml"`, `"bruker"`, ... | -- |
| | `source_path` | text | -- |
| | `pixel_size_source` | `"automatic"` / `"manual"` / `"default"` | -- |

An ontology term is always the pair
`{"accession": "MS:1000075", "name": "matrix-assisted laser desorption ionization"}` --
a [CURIE](https://www.w3.org/TR/curie/) plus the term's label, so the block
is readable without resolving anything.

`pixel_size_um` is the one acquisition field that is required: conversion
refuses to run without a pixel size, so a document without it describes no
store Thyra ever wrote.

!!! note "Unknown fields are rejected"
    Validation refuses keys the schema does not define, so a typo fails
    loudly instead of becoming an unread field. Anything genuinely
    vendor-specific already has a home in the sections beside this block
    (`format_specific`, `acquisition_params`, `raw_metadata` -- see
    [Output Format](output-format.md#provenance)).

### What is auto-populated

Auto-population is deliberately honest: a field the source does not report
is left unset, never guessed. The only inferences are facts that follow
from the format itself.

| Source | polarity | ionisation source | analyzer | instrument model |
|--------|----------|-------------------|----------|------------------|
| imzML | -- | -- | -- | from `instrument model` |
| Bruker `.d` | -- | MALDI, when the laser tables are present | TOF (timsTOF-family formats) | from the DB |
| PHI ToF-SIMS | from the header | SIMS | TOF | platform name |
| Waters `.raw` | -- | -- | -- | from `_HEADER.TXT` |

Everything else -- organism, tissue, condition, matrix, resolving power --
cannot come from a raw file and stays empty until you provide it.

---

## Storage contract

- The block lives at `table.uns["msi_metadata"]`, in every table the
  store has. This location is stable, like `uns["essential_metadata"]`
  and the `coordinate_systems` root attribute; consumers read it from
  here and nowhere else.
- It is written identically by every converter write path (the uns
  parity tests assert this), and contains no timestamps, so converting
  the same input twice produces the same block.
- Sections with nothing in them are omitted rather than written empty,
  following the store-wide convention.

## Versioning

`schema_version` follows semantic-version rules:

- Adding an optional field bumps the **minor** version.
- Renaming, removing, retyping, or making a field required bumps the
  **major** version.
- A validator implementing major version N rejects documents with a
  different major version, accepts older minors, and warns on newer
  minors.

The JSON Schema rendering is committed at
`thyra/metadata/schema/msi_metadata_schema_v0_1.json` and ships in the
wheel, so non-Python consumers can validate documents without importing
Thyra:

```python
from importlib import resources
import json

schema = json.loads(
    resources.files("thyra.metadata.schema")
    .joinpath("msi_metadata_schema_v0_1.json")
    .read_text()
)
```

A unit test keeps the artifact in sync with the models; regenerate it with
`python -m thyra.metadata.schema.generate` after a model change.

---

## `thyra validate`

```
thyra validate PATH [--merge USER.json] [--json]
```

`PATH` is a converted `.zarr` store or a standalone metadata `.json`
document. Only the metadata block is read -- validating a 100 GB store is
instant. Exit status is `0` when every document conforms (warnings
allowed) and `1` otherwise, so it can gate CI.

Errors mean the document does not conform: structural violations, unknown
PSI-MS/IMS/UO accessions, version incompatibility. Warnings mean it
conforms but says something suspicious: a term label that does not match
its accession, or a term from an unexpected ontology.

```bash
thyra validate output.zarr
# msi_dataset_z0: OK

thyra validate output.zarr --json   # machine-readable report on stdout
```

## `thyra export-metaspace`

```
thyra export-metaspace PATH [--merge USER.json] [--table NAME] [-o OUT.json]
```

Writes the METASPACE submission metadata JSON (default:
`<input>.metaspace.json` next to the input; `-o -` for stdout). Required
fields the store cannot know are emitted empty and reported as warnings on
stderr -- the output is a truthful starting point, never a fabricated
record. The one inference: matrix-free sources (DESI, SIMS) truthfully get
`MALDI_Matrix: "none"`.

```bash
thyra export-metaspace output.zarr
# warning: Sample_Information.Organism is required ... and is not set
# Wrote METASPACE metadata for msi_dataset_z0 to output.metaspace.json
```

## Completing the metadata

The fields only you can know go in a small JSON overlay, merged at
validation or export time with `--merge`:

```json
{
  "sample": {
    "organism": "Mus musculus",
    "organism_term": {"accession": "NCBITaxon:10090", "name": "Mus musculus"},
    "organism_part": "liver",
    "organism_part_term": {"accession": "UBERON:0002107", "name": "liver"},
    "condition": "wildtype"
  },
  "preparation": {
    "sample_stabilisation": "fresh frozen",
    "matrix": "2,5-dihydroxybenzoic acid (DHB)",
    "matrix_term": {"accession": "CHEBI:17189", "name": "2,5-dihydroxybenzoic acid"},
    "matrix_application": "ImagePrep"
  },
  "ms_analysis": {
    "detector_resolving_power": {"value": 130000, "at_mz": 400}
  }
}
```

```bash
thyra validate output.zarr --merge sample.json
thyra export-metaspace output.zarr --merge sample.json
```

The overlay merges key-by-key over the stored block, so it only needs the
fields you are adding. Keep it next to the store; it applies unchanged to
every dataset from the same study.

## Python API

```python
from thyra.metadata.schema import (
    MSIMetadata,             # the pydantic model
    read_msi_metadata_blocks,  # {table_name: block} from a store
    validate_document,       # (model | None, issues)
    to_metaspace,            # (submission_dict, warnings)
)

blocks = read_msi_metadata_blocks("output.zarr")
meta, issues = validate_document(blocks["msi_dataset_z0"])
submission, warnings = to_metaspace(meta)
```
