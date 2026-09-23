# thyra/metadata/schema/generate.py
"""Regenerate the committed schema artifacts.

Two copies of the JSON Schema rendering of the pydantic models exist,
and this script writes both:

* next to the models, under the name in ``SCHEMA_JSON_FILENAME``, so the
  wheel ships it and a Python consumer can load it with
  ``importlib.resources``;
* under ``docs/schema/<version>/`` in the repository, together with a
  copy of the LinkML source, so the documentation site serves both at
  the addresses ``MSI_METADATA_SCHEMA_ID`` and ``MSI_METADATA_LINKML_ID``
  name.  That folder is append-only: a published version is never edited
  again, and a new ``schema_version`` is a new folder.

``docs/schema/SHA256SUMS`` records the hash of every published file, and
a unit test fails when one stops matching, so regenerating a version
that is already published fails the suite instead of rewriting what the
site serves under that number.  This script never writes the list:
regenerating a recorded hash would hide exactly the edit it exists to
catch.  The pull request that adds a version folder appends that
folder's lines by hand, and the failing test prints them.

Unit tests assert the committed files match the models and each other;
when one fails, rerun::

    python -m thyra.metadata.schema.generate

and commit the result together with the model change (and the version
bump in ``MSI_METADATA_SCHEMA_VERSION`` that the change warrants, with
the new folder's lines in ``docs/schema/SHA256SUMS``).
"""

import json
import shutil
from pathlib import Path

from .models import MSI_METADATA_SCHEMA_VERSION, SCHEMA_JSON_FILENAME, MSIMetadata

LINKML_FILENAME = "msi_metadata.linkml.yaml"
PUBLISHED_JSON_FILENAME = "msi_metadata.schema.json"
PUBLISHED_LINKML_FILENAME = "msi_metadata.linkml.yaml"


def render_json_schema() -> str:
    """The JSON Schema text for the current models, as committed."""
    schema = MSIMetadata.model_json_schema()
    return json.dumps(schema, indent=2, sort_keys=True) + "\n"


def published_dir(repo_root: Path, version: str = MSI_METADATA_SCHEMA_VERSION) -> Path:
    """The docs folder that serves one schema version."""
    return repo_root / "docs" / "schema" / version


def main() -> None:
    """Write the package artifact, then the published copies when in a checkout."""
    here = Path(__file__).resolve().parent
    text = render_json_schema()
    target = here.with_name(here.name) / SCHEMA_JSON_FILENAME
    target.write_text(text, encoding="utf-8", newline="\n")
    print(f"Wrote {target}")

    repo_root = here.parents[2]
    docs = repo_root / "docs"
    if not docs.is_dir():
        print("No docs/ beside the package; published copies not written.")
        return
    out = published_dir(repo_root)
    out.mkdir(parents=True, exist_ok=True)
    (out / PUBLISHED_JSON_FILENAME).write_text(text, encoding="utf-8", newline="\n")
    shutil.copyfile(here / LINKML_FILENAME, out / PUBLISHED_LINKML_FILENAME)
    print(f"Wrote {out / PUBLISHED_JSON_FILENAME}")
    print(f"Wrote {out / PUBLISHED_LINKML_FILENAME}")


if __name__ == "__main__":
    main()
