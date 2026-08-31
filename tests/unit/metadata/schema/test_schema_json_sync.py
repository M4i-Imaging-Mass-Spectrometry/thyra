"""The committed JSON Schema artifact must match the pydantic models.

The artifact exists so non-Python consumers can validate documents
without importing Thyra; it is only trustworthy if it cannot drift from
the models.  When this test fails, regenerate with::

    python -m thyra.metadata.schema.generate

and commit the result together with the model change (plus the
``MSI_METADATA_SCHEMA_VERSION`` bump the change warrants: additive
optional fields bump minor, anything else bumps major).
"""

import json
from pathlib import Path

from thyra.metadata.schema.models import SCHEMA_JSON_FILENAME, MSIMetadata

_SCHEMA_DIR = Path(__file__).resolve().parents[4] / "thyra" / "metadata" / "schema"


def test_committed_json_schema_matches_the_models():
    committed = json.loads(
        (_SCHEMA_DIR / SCHEMA_JSON_FILENAME).read_text(encoding="utf-8")
    )
    assert committed == MSIMetadata.model_json_schema()


def test_schema_artifact_ships_in_the_package():
    # importlib.resources is how consumers are told to read it, so the
    # file must be importable package data, not just a repo file.
    from importlib import resources

    data = (
        resources.files("thyra.metadata.schema")
        .joinpath(SCHEMA_JSON_FILENAME)
        .read_text(encoding="utf-8")
    )
    assert json.loads(data)["title"] == "MSIMetadata"
