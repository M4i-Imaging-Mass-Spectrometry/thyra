"""Every published schema version is served at the address it names.

``docs/schema/<version>/`` is copied verbatim onto the documentation
site (issue #385).  These tests pin three things: the artifact names its
own address, the published copy of the current version is the same
bytes as the copy the wheel ships, and every version folder is
self-consistent -- its ``$id`` names that folder -- so an edit to an
already-published version fails here rather than quietly changing what
the site serves under an old number.
"""

import json
import re
from pathlib import Path

import pytest

from thyra.metadata.schema import generate, models

_SCHEMA_DIR = Path(models.__file__).resolve().parent
_REPO_ROOT = Path(__file__).resolve().parents[4]
_PUBLISHED_ROOT = _REPO_ROOT / "docs" / "schema"
_SEMVER = re.compile(r"^\d+\.\d+\.\d+$")


def _needs_checkout():
    if not _PUBLISHED_ROOT.is_dir():
        pytest.skip("published schema folders exist only in a source checkout")


def test_emitted_schema_names_its_published_address():
    schema = models.MSIMetadata.model_json_schema()
    assert schema["$id"] == models.MSI_METADATA_SCHEMA_ID
    assert schema["$schema"] == models._JSON_SCHEMA_DIALECT
    assert models.MSI_METADATA_SCHEMA_VERSION in schema["$id"]
    assert schema["$id"].startswith(models.MSI_METADATA_SCHEMA_URL_BASE + "/")


def test_linkml_source_names_its_published_address():
    yaml = pytest.importorskip("yaml")
    linkml = yaml.safe_load(
        (_SCHEMA_DIR / generate.LINKML_FILENAME).read_text(encoding="utf-8")
    )
    assert linkml["id"] == models.MSI_METADATA_LINKML_ID
    assert linkml["version"] == models.MSI_METADATA_SCHEMA_VERSION


def test_published_copy_of_the_current_version_matches_the_package():
    _needs_checkout()
    out = generate.published_dir(_REPO_ROOT)
    assert out.is_dir(), f"run python -m thyra.metadata.schema.generate ({out} missing)"
    package_json = (_SCHEMA_DIR / models.SCHEMA_JSON_FILENAME).read_bytes()
    package_linkml = (_SCHEMA_DIR / generate.LINKML_FILENAME).read_bytes()
    assert (out / generate.PUBLISHED_JSON_FILENAME).read_bytes() == package_json
    assert (out / generate.PUBLISHED_LINKML_FILENAME).read_bytes() == package_linkml


def test_every_published_version_folder_is_self_consistent():
    _needs_checkout()
    folders = sorted(p for p in _PUBLISHED_ROOT.iterdir() if p.is_dir())
    assert folders, "no published schema versions"
    for folder in folders:
        assert _SEMVER.match(folder.name), folder
        schema = json.loads(
            (folder / generate.PUBLISHED_JSON_FILENAME).read_text(encoding="utf-8")
        )
        expected = (
            f"{models.MSI_METADATA_SCHEMA_URL_BASE}/{folder.name}"
            f"/{generate.PUBLISHED_JSON_FILENAME}"
        )
        assert schema["$id"] == expected, folder
        assert (folder / generate.PUBLISHED_LINKML_FILENAME).is_file(), folder
