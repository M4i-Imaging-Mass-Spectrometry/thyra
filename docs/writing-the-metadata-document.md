# Writing the Metadata Document from Another Program

The `msi_metadata` document is not tied to Thyra. It is a JSON object
with a published schema, and any program that can describe a mass
spectrometry acquisition can write one and check it, in any language,
without installing Thyra. This page is for that program's author. The
field-by-field reference is the [Metadata Schema](metadata-schema.md)
page; this one covers the address, the contract and the versioning rule.

## The address

Every published schema version is served at a fixed URL:

```
https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/schema/<version>/msi_metadata.schema.json
https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/schema/<version>/msi_metadata.linkml.yaml
```

The current version is `0.6.0`, so the JSON Schema is at
<https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/schema/0.6.0/msi_metadata.schema.json>
and its `$id` is that same address. The LinkML source of the same version
sits beside it.

Three rules make the address safe to cite:

- **A published version is never edited.** The bytes served under
  `0.6.0` today are the bytes served under `0.6.0` in five years. A
  test in the repository fails if a published folder stops naming itself.
- **A new schema version is a new folder.** Nothing is served under a
  moving name such as `latest`. A document names its version in
  `schema_version`, and a validator fetches exactly that version.
- **The wheel ships the same file.** `importlib.resources` on
  `thyra.metadata.schema` yields byte-identical content, so a Python
  consumer offline and a validator online check against the same schema.

## What a document must contain

The schema closes every object (`additionalProperties: false`), so a
key it does not define is an error, not an extension. Where a source
has nothing to say, the field is left out; an absent field means "not
reported", never "none".

Required at the root:

| Field | What it is |
|---|---|
| `schema_version` | The version this document conforms to, `MAJOR.MINOR.PATCH`. Defaults to the current version in Thyra's models; another writer states it explicitly. |
| `ms_analysis` | How the data was acquired. Its one required member today is `pixel_size_um` (`{"x": ..., "y": ...}` in micrometres); everything else is optional. |
| `provenance` | Who wrote the document. Its one required member is `thyra_version`, which for another program is the writing software's own version string; the name is a historical accident of the field and its meaning is "the version of whatever wrote this". |

Everything else, the `sample`, `preparation`, `acquisition` and
`processing` sections and the optional members of `ms_analysis`, is
filled where known and omitted where not.

!!! note "Documents with no pixel size"
    An acquisition with no raster has no pixel size and no honest number
    for one, so under `0.6.0` its document does not validate. The split
    into a core every acquisition can fill and an imaging profile that
    adds the pitch is design decision D23 in
    [Design Decisions](design-decisions.md) and will arrive as a new
    schema version. Until then, a writer for non-imaging data should
    produce the document anyway and expect exactly that one validation
    error; nothing else about the document changes when the split lands.

## Controlled-vocabulary terms

Fields that name an instrument fact in free text have a sibling
`*_term` field carrying the controlled-vocabulary term for it, as a
CURIE accession plus its label:

```json
"analyzer": "orbitrap",
"analyzer_term": {"accession": "MS:1000484", "name": "orbitrap"}
```

Which vocabulary each field binds to is stated in the schema itself:
every bound field carries a `cv` object in the published JSON Schema,
for example `acquisition.laser_frequency_hz` binds to `IMS:1006000`.
The vocabularies in use are PSI-MS (`MS:`), the imaging MS vocabulary
(`IMS:`), the units ontology (`UO:`), NCBI Taxonomy, UBERON and ChEBI.
`thyra validate` resolves `MS`, `IMS` and `UO` accessions against the
tables it bundles: an unknown accession is an error, a label that does
not match the bundled one is a warning, and an accession from a
vocabulary other than the one designated for that field is a warning.
Accessions from the other vocabularies are checked for CURIE shape
only. The bindings are listed under
[PSI CV alignment](metadata-schema.md#psi-cv-alignment). A concept the
vocabulary does not yet have is left unset and recorded in the
[candidate terms](metadata-schema.md#candidate-cv-terms) list, and the
right fix is a term request to the vocabulary's maintainers, not a
private accession.

## Validating a document

With any JSON Schema validator, against the address:

```bash
# Python, no Thyra installed
pip install jsonschema requests
python -c "
import json, sys, jsonschema, requests
schema = requests.get('https://M4i-Imaging-Mass-Spectrometry.github.io/thyra/schema/0.6.0/msi_metadata.schema.json').json()
jsonschema.Draft202012Validator(schema).validate(json.load(open(sys.argv[1])))
print('ok')
" document.json
```

With Thyra, which additionally resolves the `MS`, `IMS` and `UO`
accessions and checks the version compatibility rule:

```bash
thyra validate document.json
```

`thyra validate` exits non-zero on an error and prints warnings for a
label that no longer matches the bundled vocabulary or a document minor
version newer than the implementation's; see
[CLI Reference](cli.md).

## The versioning rule

`schema_version` follows semantic versioning, stated in terms of the
schema and not of any code:

- **Minor**: an optional field, section or enum value is added. A
  document written against an older minor stays valid.
- **Major**: a field or enum value is removed, renamed or retyped, or an
  optional field becomes required.
- **Patch**: descriptions and documentation only; no document changes
  validity.

A validator implementing major version N rejects a document of a
different major version, accepts an older minor and warns on a newer
one. The history of versions so far is in
[Versioning](metadata-schema.md#versioning).

## Writing one from Thyra

For a format Thyra reads, `thyra metadata INPUT` writes the document
from the raw file without converting it, and
`thyra.read_metadata_document(path)` is the library form. A converted
store carries the same document in `table.uns["msi_metadata"]`, and
`read_msi_metadata_blocks` returns it in the same shape, so a document
written from a raw file and the block read out of that file's store
compare directly.
