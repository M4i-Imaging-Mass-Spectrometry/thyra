"""Every published schema version is served at the address it names.

``docs/schema/<version>/`` is copied verbatim onto the documentation
site (issue #385).  These tests pin four things: the artifact names its
own address, the published copy of the current version is the same
bytes as the copy the wheel ships, every version folder names itself in
its ``$id``, and every published file still has the bytes recorded for
it in ``docs/schema/SHA256SUMS`` (issue #395).

The last is the one that keeps a published version unchanged.  The
``$id`` check passes any edit that leaves the ``$id`` alone -- a
reworded description, a widened constraint, a field added or dropped --
which is how a corrected description in the already-published 0.6.0
once passed the full suite.  A folder is published the moment it
reaches ``main``, because the documentation deploy builds ``main``, not
a release; so the pull request that creates a version folder appends its
lines to ``SHA256SUMS``, and from then on a changed byte under that
folder fails here.  Nothing inside a repository can make such an edit
impossible, but it can no longer be silent: getting past this test means
rewriting a recorded hash, and that line shows in the diff.
"""

import hashlib
import json
import re
from pathlib import Path

import pytest

from thyra.metadata.schema import generate, models

_SCHEMA_DIR = Path(models.__file__).resolve().parent
_REPO_ROOT = Path(__file__).resolve().parents[4]
_PUBLISHED_ROOT = _REPO_ROOT / "docs" / "schema"
_SEMVER = re.compile(r"^\d+\.\d+\.\d+$")

# ``sha256sum`` format -- the digest, two spaces, the path relative to
# ``docs/schema`` -- so ``sha256sum -c`` checks a copy kept in the same
# layout, and the site serves the list beside the folders it describes.
_MANIFEST = "SHA256SUMS"
_MANIFEST_LINE = re.compile(r"^([0-9a-f]{64})  (\d+\.\d+\.\d+)/([^/\s]+)$")


def _needs_checkout():
    if not _PUBLISHED_ROOT.is_dir():
        pytest.skip("published schema folders exist only in a source checkout")


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _version_key(version):
    return tuple(int(part) for part in version.split("."))


def _read_manifest(root):
    """The recorded digest of each published path, and what is wrong with the list."""
    recorded, order, problems = {}, [], []
    lines = (root / _MANIFEST).read_text(encoding="utf-8").splitlines()
    for number, line in enumerate(lines, start=1):
        match = _MANIFEST_LINE.match(line)
        if match is None:
            problems.append(
                f"{_MANIFEST} line {number} is not '<sha256>  <version>/<file>': "
                f"{line!r}"
            )
            continue
        digest, version, name = match.groups()
        relative = f"{version}/{name}"
        if relative in recorded:
            problems.append(f"{relative} is recorded twice in {_MANIFEST}")
        recorded[relative] = digest
        order.append((_version_key(version), name))
    if order != sorted(order):
        problems.append(
            f"{_MANIFEST} is out of order: it is sorted by version, then file "
            "name, and a new version's lines go at the end"
        )
    return recorded, problems


def _recorded_file_problems(root, recorded):
    """Published files that are gone or no longer hash to their recorded digest."""
    problems = []
    for relative, digest in recorded.items():
        path = root / relative
        version = relative.split("/")[0]
        if not path.is_file():
            problems.append(
                f"{relative} was published and no longer exists: a published "
                "file is never removed"
            )
        elif _sha256(path) != digest:
            problems.append(
                f"{relative} no longer has the bytes recorded when {version} was "
                "published. A published version is never edited: restore it "
                f"(git checkout -- docs/schema/{version}) and make the change in "
                "a new schema version -- a description-only change is a PATCH "
                "bump, see docs/metadata-schema.md#versioning"
            )
    return problems


def _unrecorded_file_problems(root, recorded):
    """Entries under a version folder that the list does not record."""
    published = {relative.split("/")[0] for relative in recorded}
    problems = []
    for folder in sorted(p for p in root.iterdir() if p.is_dir()):
        for path in sorted(folder.iterdir()):
            relative = f"{folder.name}/{path.name}"
            if relative in recorded:
                continue
            if not path.is_file():
                problems.append(
                    f"{relative} is not a file; a version folder holds files"
                )
            elif folder.name in published:
                problems.append(
                    f"{relative} was added to {folder.name}, which is already "
                    "published: new content is a new schema version"
                )
            else:
                problems.append(
                    f"{relative} is not recorded in docs/schema/{_MANIFEST}. Once "
                    f"{folder.name} is final, append:\n    {_sha256(path)}  {relative}"
                )
    return problems


def _published_bytes_problems(root):
    """Every way the version folders under ``root`` differ from what was published."""
    if not (root / _MANIFEST).is_file():
        return [f"{root / _MANIFEST} is missing"]
    recorded, problems = _read_manifest(root)
    problems += _recorded_file_problems(root, recorded)
    problems += _unrecorded_file_problems(root, recorded)
    return problems


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


def test_every_published_version_folder_names_itself():
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


def test_every_published_file_still_has_its_recorded_bytes():
    _needs_checkout()
    problems = _published_bytes_problems(_PUBLISHED_ROOT)
    assert not problems, "\n".join(problems)


def _publish(root, version, files):
    """Write one version folder and record it, as merging it does."""
    folder = root / version
    folder.mkdir()
    for name, content in files.items():
        (folder / name).write_bytes(content)
    with (root / _MANIFEST).open("a", encoding="utf-8", newline="\n") as manifest:
        for name in sorted(files):
            manifest.write(f"{_sha256(folder / name)}  {version}/{name}\n")


class TestTheBytesGuard:
    """The guard reports each way a published version can change.

    A guard that never fails proves nothing, so each case publishes a
    small tree, changes it in the one way the test names, and checks
    that the report names it.
    """

    @pytest.fixture
    def root(self, tmp_path):
        root = tmp_path / "schema"
        root.mkdir()
        _publish(
            root,
            "0.1.0",
            {"a.schema.json": b'{"description": "old"}\n', "a.yaml": b"id: a\n"},
        )
        _publish(root, "0.2.0", {"a.schema.json": b'{"description": "new"}\n'})
        return root

    def test_an_untouched_tree_passes(self, root):
        assert _published_bytes_problems(root) == []

    def test_a_reworded_description_is_caught(self, root):
        path = root / "0.1.0" / "a.schema.json"
        path.write_bytes(path.read_bytes().replace(b"old", b"corrected"))
        (problem,) = _published_bytes_problems(root)
        assert problem.startswith("0.1.0/a.schema.json no longer has the bytes")
        assert "git checkout -- docs/schema/0.1.0" in problem

    def test_a_removed_file_is_caught(self, root):
        (root / "0.1.0" / "a.yaml").unlink()
        (problem,) = _published_bytes_problems(root)
        assert problem.startswith("0.1.0/a.yaml was published and no longer exists")

    def test_a_file_added_to_a_published_version_is_caught(self, root):
        (root / "0.1.0" / "b.json").write_bytes(b"{}\n")
        (problem,) = _published_bytes_problems(root)
        assert problem.startswith("0.1.0/b.json was added to 0.1.0")

    def test_a_new_version_is_told_the_lines_to_append(self, root):
        folder = root / "0.10.0"
        folder.mkdir()
        (folder / "a.schema.json").write_bytes(b"{}\n")
        (problem,) = _published_bytes_problems(root)
        line = f"{_sha256(folder / 'a.schema.json')}  0.10.0/a.schema.json"
        assert problem.endswith("append:\n    " + line)

    def test_a_new_version_is_ordered_by_number_not_by_spelling(self, root):
        # 0.10.0 sorts before 0.2.0 as text; the list must not care.
        _publish(root, "0.10.0", {"a.schema.json": b"{}\n"})
        assert _published_bytes_problems(root) == []

    def test_lines_out_of_order_are_caught(self, root):
        manifest = root / _MANIFEST
        lines = manifest.read_text(encoding="utf-8").splitlines()
        manifest.write_text(
            "\n".join(reversed(lines)) + "\n", encoding="utf-8", newline="\n"
        )
        (problem,) = _published_bytes_problems(root)
        assert problem.startswith(f"{_MANIFEST} is out of order")

    def test_a_path_recorded_twice_is_caught(self, root):
        manifest = root / _MANIFEST
        last = manifest.read_text(encoding="utf-8").splitlines()[-1]
        with manifest.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(last + "\n")
        (problem,) = _published_bytes_problems(root)
        assert problem == f"0.2.0/a.schema.json is recorded twice in {_MANIFEST}"

    def test_a_malformed_line_is_caught(self, root):
        with (root / _MANIFEST).open("a", encoding="utf-8", newline="\n") as handle:
            handle.write("0.1.0/a.yaml  d41d8cd98f00b204e9800998ecf8427e\n")
        (problem,) = _published_bytes_problems(root)
        assert problem.startswith(f"{_MANIFEST} line 4 is not")

    def test_a_missing_list_is_caught(self, root):
        (root / _MANIFEST).unlink()
        (problem,) = _published_bytes_problems(root)
        assert problem.endswith(f"{_MANIFEST} is missing")
