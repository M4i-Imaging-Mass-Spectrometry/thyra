"""A store-identity harness: does a change move any converted store at all?

Step 0 of the #276 decomposition plan (thyra#349). Every later step in that
plan -- and thyra#348 -- is meant to be behaviour-preserving, and
"preserving" has to be measured against a real store, not asserted from
reading the diff. This script is the one instrument every step runs instead
of improvising a comparison by hand.

It has two subcommands.

``run`` converts a set of named datasets, one conversion per subprocess, and
writes a JSON hash manifest of the resulting Zarr store next to a JSONL
checkpoint. A checkpoint entry marked ``ok`` is skipped on a re-run unless
``--force`` is given, so a crash partway through costs one conversion, not
the whole set. The skip only applies when the earlier run used the same
settings: the same input (and the size and time of its files), ``thyra``
options, source checkout and the ``thyra`` files on disk (uncommitted
edits included), interpreter and its installed packages, and mask keys.
A name run before under other settings is converted again. ``--source``
must hold the ``thyra`` package the child then imports; ``run`` checks.

The child ``thyra`` runs with its standard input closed, so a conversion
that would ask a question fails at once instead of waiting for an answer
nobody sees. A failed record keeps the ERROR lines from both output
streams: Thyra logs to stdout, so stderr alone rarely says why.

``compare`` diffs the manifests of two labelled runs and prints one line per
Zarr element or other file that differs, then exits 1 if anything did.
These count as differences too, so a label never answers for a store it
did not make: a label whose manifests were made under different settings,
a manifest whose latest checkpoint record is not ``ok``, and a name that
failed or was skipped with a store in neither label. ``compare`` exits 2
when neither label holds a manifest, or when both name one folder:
comparing nothing is not the same as finding no difference.

Comparing a branch against main
--------------------------------

Run the same dataset set through both checkouts under different labels, then
compare the labels::

    python tests/tools/store_identity.py run --label main \\
        --source /path/to/main-worktree --datasets pea bellini
    python tests/tools/store_identity.py run --label mybranch \\
        --source /path/to/branch-worktree --datasets pea bellini
    python tests/tools/store_identity.py compare main mybranch

``--source`` defaults to the repository this script lives in and
``--interpreter`` to the running Python, so the common case (same checkout,
comparing two option sets, e.g. before/after a resampling flag) needs
neither::

    python tests/tools/store_identity.py run --label before --datasets pea
    python tests/tools/store_identity.py run --label after --datasets pea \\
        -- --resample-bins 500
    python tests/tools/store_identity.py compare before after

How dataset names resolve
--------------------------

A ``NAME`` on the command line is never a path. It resolves in this order,
and a name that resolves nowhere is skipped with a message -- never a
failure:

1. The machine-local registry at ``~/.msi-datasets.toml`` (sections
   flattened, ``name = "path"``). A name that *is* registered but whose
   path does not exist here (a lab share that is offline, say) is a
   different skip from "not found anywhere" -- it says so. A registered
   folder that is not itself a ``.d`` but holds ``*.d`` runs resolves to
   its only ``.d``; with several, the name is skipped with a message, as
   ``thyra`` would otherwise ask which one to convert.
2. Each directory listed in the ``THYRA_DATASET_ROOTS`` environment
   variable (``os.pathsep``-separated), matched by basename against
   ``<root>/<name>``, ``<root>/<name>.imzML`` and ``<root>/<name>.d``. A
   directory match that is not itself a ``.d`` and holds exactly one
   ``*.d`` resolves to that ``.d``. A name matching *both*
   ``<name>.imzML`` and ``<name>.d`` in the same root (a case-insensitive
   filesystem can hold both for stems differing only in letter case) is
   also skipped rather than guessed at -- give the name with its suffix to
   disambiguate.
3. ``tests/data/fixtures`` in this repository, matched the same way,
   searched last.

The default dataset set
------------------------

With no ``--datasets``, ``run`` converts every fixture under
``tests/data/fixtures`` that actually converts with ``--no-optical`` --
checked by hand once, see ``DEFAULT_DATASETS`` below. That default is a
smoke set, not a regression set: fixtures are small and synthetic, built to
exercise one edge case each, not to stand in for a real acquisition on any
of Thyra's real routes.

For a real comparison, pass dataset names explicitly by *role* -- an imzML
export for the shared-axis nearest-neighbour path, another for
TIC-preserving, a Bruker TDF with TIMS off, a TIMS slide with mobility
engaged, an MS/MS PASEF schedule (the registry's ``tims_msms_neg_13win`` is
one), and a FlexImaging acquisition carrying a ``.mis`` alignment image.
None of those names are fixed here; look them up in
``~/.msi-datasets.toml`` or a ``THYRA_DATASET_ROOTS`` directory on the
machine at hand and pass them with ``--datasets``.

What is hashed
---------------

Every Zarr array (shape, dtype, data and its own metadata), every group's
attrs, and every other file in the store -- the pixel boxes in
``shapes/*/shapes.parquet`` are not Zarr arrays, so they are hashed as
files, byte for byte. Attrs are compared as sorted JSON text, so ``1``,
``1.0`` and ``true`` are three different values, as they are on disk.

What is masked before hashing
------------------------------

Four keys are dropped because they vary without the store's data changing:
``conversion_timestamp`` changes on every run, and ``thyra_version``,
``msi_converter_version`` and ``produced_by`` carry the version number, so
two checkouts that differ only in ``__version__`` compare equal. The same
version also sits in the ``software`` entries of the ``processing`` list;
an entry naming ``thyra`` has its ``version`` dropped, and any other
software keeps it. Masking works two ways, because ``uns`` is not attrs:

- **Attrs dicts** (a group's own ``.attrs``, an array's ``.metadata``
  ``attributes`` block) are masked by key, recursively, at any depth.
- **anndata writes ``uns`` as zarr groups and arrays, not attrs** -- a
  scalar such as ``thyra_version`` is a 0-d array literally named
  ``thyra_version``, not an attrs entry -- so ``build_manifest`` also drops
  any array or group with a mask key anywhere in its path, descendants
  included, wherever it sits. And a string *array element* that parses as JSON
  (``uns`` payloads such as the ``cvParams`` / ``regions`` blocks are
  stored as JSON strings -- see ``thyra/metadata/uns_compat.py``) is
  parsed and masked the same way (see ``_hash_block``). A payload with
  nothing to mask is hashed exactly as written. One that held a masked key
  is re-serialised with its key order kept, so spacing and the spelling
  of a number (``1.00`` against ``1.0``) inside it are not compared.

``--mask KEY`` adds more keys to drop from both paths at once. Reach for it
when the same dataset is read from two different locations. The input path
is recorded under six keys, and all six must be masked, in both runs::

    python tests/tools/store_identity.py run --label here \\
        --mask source_path --mask ibd_file --mask data_path \\
        --mask database_path --mask binary_file --mask mis_file \\
        --datasets pea

``source_path`` is the absolute input path recorded on every table;
``ibd_file`` (imzML), ``data_path`` (Bruker and Waters), ``database_path``
and ``binary_file`` (Bruker) and ``mis_file`` (solariX) are the locating
keys ``thyra/metadata/personal_data.py`` keeps in the vendor metadata.

Convert the same datasets under two labels with the same settings and
compare them to find out whether the mask list is complete for a given
store shape -- a non-zero ``compare`` there means it is not. A blanket
ISO-8601-in-string regex was tried and dropped: it would mask a genuinely
changed ``acquisition_datetime`` (real source metadata, not a run artefact)
exactly as it masks a volatile one, and measuring against both the fixture
round and a real round found nothing left over that needed it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import tomllib
import uuid
from pathlib import Path
from typing import Any, Optional

import numpy as np
import zarr

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent

DEFAULT_MASK_KEYS = {
    "conversion_timestamp",
    "thyra_version",
    "msi_converter_version",
    "produced_by",
}

# The software entry of a ``processing`` step that names Thyra itself: its
# ``version`` is the running version, like the keys above, and is dropped.
_THYRA_SOFTWARE_NAME = "thyra"

# Every tests/data/fixtures file that actually converts with --no-optical,
# checked by hand: each of tests/data/fixtures/*.imzML and synthetic_tims.d
# was run through `thyra` once. Four are left out because they are fixtures
# *for* a refusal and exit non-zero on purpose, not by defect:
# two_precision_terms (two precision terms declared for the m/z array), and
# pixel_size_area_old_name, pixel_size_area_old_name_no_extent and
# pixel_size_unit_contradiction (no pixel size Thyra accepts without
# --pixel-size).
DEFAULT_DATASETS = (
    "iontof_sparse",
    "mobility_continuous",
    "mobility_processed",
    "pixel_size_unit_absent",
    "pixel_size_unit_declared",
    "solarix_fticr",
    "synthetic_tims",
    "timstof_flex_export",
    "two_scansettings",
    "unit_nanometre",
)

# Zarr's own metadata documents. Everything else in a store that is not an
# array's chunk data is a file Zarr knows nothing about and is hashed whole.
_ZARR_METADATA_NAMES = {"zarr.json", ".zattrs", ".zgroup", ".zarray", ".zmetadata"}

# What a ``failed`` record keeps of the child's output.
_ERROR_MARKERS = (" - ERROR - ", " - CRITICAL - ")
_MAX_ERROR_LINES = 20

_HASH_CHUNK_BYTES = 64 * 1024 * 1024  # 64 MiB, chunked so a big table array
# is hashed without holding it whole in memory.

# StringDType ('T', what anndata/zarr write today), fixed-width unicode
# ('U'), bytes ('S') and object ('O') arrays all hold their text out of
# line from the fixed-size element .tobytes() would hash -- a StringDType
# element is either a 16-byte inline buffer or an arena pointer, neither of
# which is the string. Those kinds are hashed element-wise instead.
_STRING_DTYPE_KINDS = {"T", "U", "S", "O"}


# ---------------------------------------------------------------------------
# Dataset name resolution
# ---------------------------------------------------------------------------


class DatasetPathUnreachable(Exception):
    """NAME is in the registry but the path it names does not exist on
    this machine (e.g. a lab share that is offline)."""


def _d_runs(folder: Path) -> list:
    """The ``.d`` run folders directly inside ``folder``, counted the way
    ``thyra`` counts them: folders only, the suffix in any letter case."""
    return sorted(
        child
        for child in folder.iterdir()
        if child.is_dir() and child.suffix.lower() == ".d"
    )


def _registry_lookup(name: str) -> Optional[Path]:
    registry = Path.home() / ".msi-datasets.toml"
    if not registry.is_file():
        return None
    try:
        with registry.open("rb") as handle:
            data = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError):
        return None
    flat: dict[str, Any] = {}
    for section in data.values():
        if isinstance(section, dict):
            flat.update(section)
    value = flat.get(name)
    if not isinstance(value, str):
        return None
    path = Path(value)
    if not path.exists():
        raise DatasetPathUnreachable(
            f"{name}: registered but its path is not reachable on this machine"
        )
    if path.is_dir() and path.suffix.lower() != ".d":
        # A folder of runs, the way _match_in_root reads one: passed as it
        # is, thyra would ask on stdin which .d to convert.
        hits = _d_runs(path)
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            listing = ", ".join(hit.name for hit in hits)
            raise AmbiguousDatasetName(
                f"{name}: the registered folder holds {len(hits)} .d runs "
                f"({listing}); register the one to convert"
            )
    return path


class AmbiguousDatasetName(Exception):
    """A NAME matches more than one file in a root (case-insensitive
    filesystems make ``name.imzML`` and ``name.d`` both exist for a stem
    that differs only in letter case, as the Xenium pair does)."""


def _match_in_root(root: Path, name: str) -> Optional[Path]:
    # An exact match first, so a name may carry its own suffix and sidestep
    # the ambiguity below entirely.
    exact = root / name
    if exact.exists():
        if exact.is_dir() and exact.suffix.lower() != ".d":
            hits = _d_runs(exact)
            if len(hits) == 1:
                return hits[0]
            # Zero or more-than-one .d child: not a resolvable exact match;
            # fall through to the suffix-guessing candidates below.
        else:
            return exact

    candidates = [c for c in (root / f"{name}.imzML", root / f"{name}.d") if c.exists()]
    if len(candidates) > 1:
        listing = ", ".join(c.name for c in candidates)
        raise AmbiguousDatasetName(
            f"{name}: ambiguous in {root} ({listing}); give the name with its suffix"
        )
    return candidates[0] if candidates else None


def resolve_dataset(name: str) -> Optional[Path]:
    """Resolve a dataset NAME to a path, or None. See the module docstring.

    Raises ``AmbiguousDatasetName`` rather than guessing when a name
    matches more than one file in the same root, and
    ``DatasetPathUnreachable`` when the name is registered but its path
    does not exist here -- neither is the same as "not found", and the
    caller reports them differently.
    """
    found = _registry_lookup(name)
    if found is not None:
        return found
    env_roots = os.environ.get("THYRA_DATASET_ROOTS", "")
    roots = [Path(entry) for entry in env_roots.split(os.pathsep) if entry]
    roots.append(REPO_ROOT / "tests" / "data" / "fixtures")
    for root in roots:
        found = _match_in_root(root, name)
        if found is not None:
            return found
    return None


# ---------------------------------------------------------------------------
# Masking
# ---------------------------------------------------------------------------


def _maybe_json(text: str) -> Any:
    stripped = text.strip()
    if not stripped or stripped[0] not in "{[":
        return None
    try:
        return json.loads(stripped)
    except (json.JSONDecodeError, ValueError):
        return None


def mask_value(value: Any, mask_keys: set) -> Any:
    """Recursively drop ``mask_keys``, masking JSON-string payloads the same
    way (see the module docstring for what was tried and measured out)."""
    if isinstance(value, dict):
        names_thyra = value.get("name") == _THYRA_SOFTWARE_NAME
        return {
            key: mask_value(item, mask_keys)
            for key, item in value.items()
            if key not in mask_keys and not (names_thyra and key == "version")
        }
    if isinstance(value, list):
        return [mask_value(item, mask_keys) for item in value]
    if isinstance(value, str):
        parsed = _maybe_json(value)
        if parsed is None:
            return value
        masked = mask_value(parsed, mask_keys)
        if json.dumps(masked) == json.dumps(parsed):
            # Nothing was masked: keep the text exactly, so a change in key
            # order or spacing inside a stored payload still shows.
            return value
        # Key order is kept; only the spacing is redone.
        return json.dumps(masked, default=str)
    return value


# ---------------------------------------------------------------------------
# Store hashing
# ---------------------------------------------------------------------------


def _hash_block(digest: "hashlib._Hash", block: np.ndarray, mask_keys: set) -> None:
    if block.dtype.kind in _STRING_DTYPE_KINDS:
        # Element-wise: str(x), UTF-8 encoded, length-prefixed so that
        # ("ab", "c") and ("a", "bc") hash differently. Each element also
        # goes through mask_value so a JSON payload stored as a string
        # (uns's cvParams/regions blocks -- see thyra/metadata/uns_compat.py)
        # is masked the same way an attrs string is.
        for element in block.reshape(-1):
            encoded = mask_value(str(element), mask_keys).encode("utf-8")
            digest.update(len(encoded).to_bytes(4, "big"))
            digest.update(encoded)
    else:
        digest.update(np.ascontiguousarray(block).tobytes())


def _hash_array(array: "zarr.Array", mask_keys: set) -> str:
    digest = hashlib.sha256()
    shape = array.shape
    if not shape:
        # A 0-d array's [...] unboxes to a bare Python scalar rather than a
        # 0-d ndarray (StringDType does this reliably), so this is handled
        # directly against the array's own declared dtype kind rather than
        # routed through _hash_block, which needs an ndarray to inspect.
        value = array[...]
        if array.dtype.kind in _STRING_DTYPE_KINDS:
            encoded = mask_value(str(value), mask_keys).encode("utf-8")
            digest.update(len(encoded).to_bytes(4, "big"))
            digest.update(encoded)
        else:
            digest.update(np.asarray(value, dtype=array.dtype).tobytes())
        return digest.hexdigest()
    row_bytes = array.dtype.itemsize
    for dim in shape[1:]:
        row_bytes *= dim
    rows_per_chunk = max(1, _HASH_CHUNK_BYTES // row_bytes) if row_bytes else shape[0]
    total_rows = shape[0]
    for start in range(0, total_rows, rows_per_chunk):
        block = array[start : start + rows_per_chunk]
        _hash_block(digest, block, mask_keys)
    return digest.hexdigest()


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def _hash_other_files(store_path: Path, array_paths: set, mask_keys: set) -> dict:
    """Hash every file in the store that is neither Zarr metadata nor an
    array's chunk data -- the pixel boxes' ``shapes.parquet``, for one.

    The Zarr walk never reaches these, so without this a store whose
    pixel boxes moved compared equal to one whose boxes did not.
    """
    array_parts = {tuple(path.split("/")) for path in array_paths}
    files: dict[str, Any] = {}
    for path in sorted(store_path.rglob("*")):
        if not path.is_file() or path.name in _ZARR_METADATA_NAMES:
            continue
        parts = path.relative_to(store_path).parts
        if any(part in mask_keys for part in parts):
            continue
        if any(parts[:end] in array_parts for end in range(1, len(parts))):
            continue
        files["/".join(parts)] = {
            "size": path.stat().st_size,
            "sha256": _hash_file(path),
        }
    return files


def build_manifest(store_path: Path, mask_keys: set) -> dict:
    """Walk the whole store: every array's shape/dtype/hash/metadata,
    every group's (masked) attrs, including the root, and every other file.

    anndata writes ``uns`` as zarr groups and arrays, not attrs -- a scalar
    like ``thyra_version`` is a 0-d array named ``thyra_version``, not an
    attrs key -- so masking there works by path component: any array or
    group with a component of its path in ``mask_keys`` is left out of the
    manifest entirely, at whatever depth it sits, descendants included.

    An array's own ``.metadata`` (chunk grid, shard/codec configuration,
    fill value, and its own ``attributes`` -- anndata's ``encoding-type`` /
    ``encoding-version`` live there, not on the parent group) is masked and
    kept too, so a chunking or encoding change that leaves the bytes
    unchanged still shows up as a difference.

    Files the Zarr walk does not reach are hashed byte for byte under
    ``files`` (see ``_hash_other_files``).

    The walk reads each node's own metadata, not the consolidated copy in
    the root ``zarr.json``, and that copy is kept on its own under
    ``consolidated``: a store whose consolidated block went stale, or
    went missing, differs from one whose block did not.
    """
    root = zarr.open_group(str(store_path), mode="r", use_consolidated=False)
    groups: dict[str, Any] = {"": mask_value(dict(root.attrs), mask_keys)}
    arrays: dict[str, Any] = {}
    all_array_paths = set()
    for path, node in root.members(max_depth=None):
        if not isinstance(node, zarr.Group):
            all_array_paths.add(path)
        # The walk is flat, so a masked group's descendants arrive as their
        # own entries: test every component, not just the last one.
        if any(part in mask_keys for part in path.split("/")):
            continue
        if isinstance(node, zarr.Group):
            groups[path] = mask_value(dict(node.attrs), mask_keys)
        else:
            arrays[path] = {
                "shape": list(node.shape),
                "dtype": str(node.dtype),
                "sha256": _hash_array(node, mask_keys),
                "metadata": mask_value(node.metadata.to_dict(), mask_keys),
            }
    files = _hash_other_files(Path(store_path), all_array_paths, mask_keys)
    return {
        "groups": groups,
        "arrays": arrays,
        "files": files,
        "consolidated": _consolidated_block(Path(store_path), mask_keys),
    }


def _consolidated_block(store_path: Path, mask_keys: set) -> Any:
    """The root ``zarr.json``'s consolidated metadata, masked like the
    walk: an entry under a masked path component is left out, and every
    other entry has its mask keys dropped. ``None`` when there is none."""
    root_document = store_path / "zarr.json"
    if not root_document.is_file():
        return None
    block = json.loads(root_document.read_text(encoding="utf-8")).get(
        "consolidated_metadata"
    )
    if not isinstance(block, dict):
        return block
    entries = block.get("metadata", {})
    kept = {
        path: mask_value(entry, mask_keys)
        for path, entry in entries.items()
        if not any(part in mask_keys for part in path.split("/"))
    }
    return mask_value({**block, "metadata": kept}, mask_keys)


# ---------------------------------------------------------------------------
# ``run``
# ---------------------------------------------------------------------------


def _git_rev(source: Path) -> Optional[str]:
    try:
        proc = subprocess.run(
            ["git", "-C", str(source), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return proc.stdout.strip()


def _code_hash(source: Path) -> str:
    """A hash of the ``thyra`` package the child imports from ``source``.

    The commit alone does not say what runs: the child imports the files
    on disk, uncommitted and untracked edits included, so two runs at
    one commit with different edits must not share a cached store.
    """
    digest = hashlib.sha256()
    package = source / "thyra"
    for path in sorted(package.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        name = path.relative_to(package).as_posix().encode("utf-8")
        digest.update(len(name).to_bytes(4, "big") + name)
        digest.update(_hash_file(path).encode("ascii"))
    return digest.hexdigest()


_LIST_DISTRIBUTIONS = (
    "import importlib.metadata as m, json, sys; "
    "print(sys.version); "
    "print(json.dumps(sorted(f\"{d.metadata.get('Name')}=={d.version}\" "
    "for d in m.distributions())))"
)

# Variables that change what the child runs: an inherited PYTHONPATH can
# put another copy of a dependency first, and the SDK paths override the
# vendor libraries Thyra ships.
_ENVIRONMENT_VARIABLES = ("PYTHONPATH", "BRUKER_SDK_PATH", "WATERS_SDK_PATH")


def _environment_hash(interpreter: str) -> str:
    """A hash of the interpreter's Python version and installed packages,
    and of ``_ENVIRONMENT_VARIABLES``, so a venv rebuilt under the same
    path, or a shell with another ``PYTHONPATH``, does not reuse a store.

    When the interpreter cannot list its packages, the value is unique to
    this run, so it matches no other: an unknown environment never lets a
    cached store, or a mixed label, pass.
    """
    try:
        proc = subprocess.run(
            [interpreter, "-c", _LIST_DISTRIBUTIONS],
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        print(f"[warn] cannot list the packages of {interpreter}: {exc}")
        return f"unknown-{uuid.uuid4().hex}"
    variables = {name: os.environ.get(name) for name in _ENVIRONMENT_VARIABLES}
    text = proc.stdout + json.dumps(variables, sort_keys=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _fingerprint_files(files: list, base: Path, digest: Any) -> None:
    for path in sorted(files):
        stat = path.stat()
        line = (
            f"{path.relative_to(base).as_posix()}\0{stat.st_size}\0{stat.st_mtime_ns}\n"
        )
        digest.update(line.encode("utf-8"))


def _input_fingerprint(resolved: Path, thyra_options: list, source: Path) -> str:
    """Size and modification time of the files a conversion may read.

    That is every file of a folder input (a ``.d``, a ``.raw``), every
    file beside the input (an imzML's ``.ibd``; the ``.mis`` and optical
    image a Bruker reader looks for next to a ``.d``), and any file an
    option names, except the ``--log-file`` the child itself writes. The
    path alone keyed the cache, so a dataset re-exported to the same place
    reused the old store. Content hashes would read gigabytes per run;
    size and time catch a rewrite.
    """
    digest = hashlib.sha256()
    if resolved.is_dir():
        files = [path for path in resolved.rglob("*") if path.is_file()]
        _fingerprint_files(files, resolved, digest)
    beside = [path for path in resolved.parent.iterdir() if path.is_file()]
    _fingerprint_files(beside, resolved.parent, digest)
    previous = None
    for option in thyra_options:
        flag, _, value = option.partition("=")
        if previous == "--log-file" or flag == "--log-file":
            previous = option
            continue
        previous = option
        # ``--flag=value`` names its file after the "="; the child runs in
        # ``source``, so a relative path is read there.
        text = value if option.startswith("--") and value else option
        named = Path(text) if Path(text).is_absolute() else source / text
        if named.is_file():
            digest.update(option.encode("utf-8") + b"\0")
            _fingerprint_files([named], named.parent, digest)
    return digest.hexdigest()


def _child_env(source: Path) -> dict:
    """The child's environment: ``source`` first on ``PYTHONPATH``."""
    env = dict(os.environ)
    env["PYTHONPATH"] = str(source) + (
        os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else ""
    )
    return env


class SourceNotUsed(Exception):
    """``--source`` is not the ``thyra`` the child would import."""


def _check_source(source: Path, interpreter: str) -> None:
    """Make sure the child imports ``thyra`` from ``source``.

    A folder without ``thyra/`` (the package folder itself, say, or a
    checkout's parent) let the child import whatever ``thyra`` the
    interpreter had installed -- an editable install of another checkout
    -- while the label recorded ``source``.
    """
    if not (source / "thyra" / "__init__.py").is_file():
        raise SourceNotUsed(f"{source} holds no thyra package (thyra/__init__.py)")
    try:
        proc = subprocess.run(
            [interpreter, "-c", "import thyra; print(thyra.__file__)"],
            cwd=str(source),
            env=_child_env(source),
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SourceNotUsed(f"{interpreter} cannot import thyra from {source}: {exc}")
    imported = Path(proc.stdout.strip().splitlines()[-1]).resolve()
    if not imported.is_relative_to((source / "thyra").resolve()):
        raise SourceNotUsed(
            f"the child imports thyra from {imported}, not from {source}"
        )


def _load_checkpoint(checkpoint: Path) -> list:
    if not checkpoint.is_file():
        return []
    records = []
    for line in checkpoint.read_text().splitlines():
        line = line.strip()
        if line:
            records.append(json.loads(line))
    return records


def _append_checkpoint(checkpoint: Path, record: dict) -> None:
    with checkpoint.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")


def _latest_by_name(records: list) -> dict:
    """The last checkpoint record for each name wins -- a name reconverted
    after an earlier failure (or under --force) is judged on that attempt,
    not on whatever ran first."""
    latest: dict[str, dict] = {}
    for record in records:
        latest[record["name"]] = record
    return latest


# The record fields a cached ok run must share with the current call to be
# reused. A record written before these were recorded lacks some of them,
# and so never matches: it is converted again. The commit (``source_rev``)
# is recorded but not compared: ``code`` covers what the child runs, and a
# docs-only commit would otherwise redo every name.
_SETTINGS_KEYS = (
    "input",
    "input_state",
    "options",
    "source",
    "code",
    "interpreter",
    "environment",
    "mask_keys",
)


def _settings_differ(record: dict, settings: dict, keys=_SETTINGS_KEYS) -> list:
    """The names of the settings ``record`` was made under that differ
    from ``settings``; empty when the cached run can stand for this one."""
    return [key for key in keys if record.get(key) != settings.get(key)]


def _error_lines(*outputs: str) -> list:
    """The child's ERROR and CRITICAL log lines, from every stream.

    Thyra's log handler writes to stdout, so stderr alone (all a failed
    record kept before) holds the traceback at most, not the reason.
    """
    lines = [
        line.rstrip()
        for output in outputs
        for line in (output or "").splitlines()
        if any(marker in line for marker in _ERROR_MARKERS)
    ]
    return lines[-_MAX_ERROR_LINES:]


def _resolve_or_skip(
    name: str, manifest_path: Path, settings: dict, checkpoint: Path, results: list
) -> Optional[Path]:
    """Resolve ``name``; when it does not resolve, record the skip and
    return ``None``."""
    try:
        resolved = resolve_dataset(name)
    except (AmbiguousDatasetName, DatasetPathUnreachable) as exc:
        resolved, reason = None, str(exc)
    else:
        reason = (
            f"{name}: not found in the registry, THYRA_DATASET_ROOTS, "
            "or tests/data/fixtures"
        )
    if resolved is not None:
        return resolved

    print(f"[skip] {reason}")
    # This attempt made no store, so an earlier manifest must not answer
    # for it in `compare`; the record says it was skipped.
    if manifest_path.exists():
        manifest_path.unlink()
        print(f"[skip] {name}: removed its manifest from an earlier run")
    known = [key for key in _SETTINGS_KEYS if not key.startswith("input")]
    record = {"name": name, **{key: settings.get(key) for key in known}}
    record.update(status="skipped", reason=reason)
    _append_checkpoint(checkpoint, record)
    results.append(record)
    return None


def _mark_failed(record: dict, proc: subprocess.CompletedProcess) -> None:
    """Fill in a failed conversion's record and print why it failed."""
    record["status"] = "failed"
    record["returncode"] = proc.returncode
    record["error_lines"] = _error_lines(proc.stdout, proc.stderr)
    record["stderr_tail"] = proc.stderr[-4000:]
    # The first ERROR line is the reason; the last is the CLI's closing
    # "Conversion failed."
    why = record["error_lines"][0] if record["error_lines"] else "no ERROR line"
    print(
        f"[failed] {record['name']} after {record['wall_seconds']}s "
        f"(exit {proc.returncode}): {why}"
    )


def run_datasets(
    label: str,
    datasets: list,
    source: Path,
    interpreter: str,
    scratch: Path,
    thyra_options: list,
    mask_keys: set,
    force: bool,
) -> list:
    """Convert ``datasets`` under ``label``, returning one result record per
    name in ``datasets`` (status ``ok`` / ``failed`` / ``hash_failed`` /
    ``skipped``). Prints a one-line summary before returning.

    The checkpoint file is never deleted, ``--force`` or not: it holds one
    line per attempt, and only the latest attempt per name decides whether
    that name is done. A prior ok is reused only when it was made under the
    same settings (``_SETTINGS_KEYS``: the input and the size and time of
    its files, options, source and the ``thyra`` files on disk,
    interpreter and its installed packages, mask keys); otherwise the
    name is converted again.
    ``--force`` reconverts the requested names regardless of a prior ok, and
    appends a fresh record rather than discarding the checkpoint's record of
    every other name under this label.
    """
    label_dir = scratch / label
    label_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = label_dir / "checkpoint.jsonl"

    latest = _latest_by_name(_load_checkpoint(checkpoint))
    done = {name: rec for name, rec in latest.items() if rec.get("status") == "ok"}

    source = source.resolve()
    _check_source(source, interpreter)
    settings = {
        "options": list(thyra_options),
        "source": str(source),
        "source_rev": _git_rev(source),
        "code": _code_hash(source),
        "interpreter": interpreter,
        "environment": _environment_hash(interpreter),
        "mask_keys": sorted(mask_keys),
    }
    results = []

    for name in datasets:
        manifest_path = label_dir / f"{name}.hashes.json"
        resolved = _resolve_or_skip(name, manifest_path, settings, checkpoint, results)
        if resolved is None:
            continue

        settings["input"] = str(resolved)
        settings["input_state"] = _input_fingerprint(resolved, thyra_options, source)
        if not force and name in done:
            differing = _settings_differ(done[name], settings)
            if not differing and manifest_path.exists():
                print(f"[skip-cached] {name} (label {label!r} already has an ok run)")
                results.append(done[name])
                continue
            if differing:
                print(
                    f"[rerun] {name}: label {label!r} was made with other "
                    f"{', '.join(differing)}; converting again"
                )

        # Delete any manifest from a previous attempt before this one starts,
        # so a failed conversion can never leave a stale ok manifest behind
        # for `compare` to silently read.
        if manifest_path.exists():
            manifest_path.unlink()

        out_store = label_dir / f"{name}.zarr"
        if out_store.exists():
            shutil.rmtree(out_store)

        cmd = [
            interpreter,
            "-m",
            "thyra",
            str(resolved),
            str(out_store),
            *thyra_options,
        ]
        env = _child_env(source)

        print(f"[run] {name} <- {resolved}")
        started = time.perf_counter()
        # stdin closed: a conversion that would ask a question (a folder of
        # several .d runs, say) fails at once instead of waiting forever.
        proc = subprocess.run(
            cmd,
            cwd=str(source),
            env=env,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
        )
        wall_seconds = round(time.perf_counter() - started, 3)

        record = {"name": name, **settings, "wall_seconds": wall_seconds}

        # A branch switch or an edit in the source during a long run would
        # otherwise record later names under the code the run started with.
        code_now = _code_hash(source)
        if code_now != settings["code"]:
            settings["code"] = code_now
            settings["source_rev"] = _git_rev(source)
            if proc.returncode == 0:
                record["status"] = "failed"
                record["error_lines"] = [
                    "the thyra code in the source changed during the run"
                ]
                print(f"[failed] {name}: {record['error_lines'][0]}; run it again")
                _append_checkpoint(checkpoint, record)
                results.append(record)
                continue

        if proc.returncode != 0:
            _mark_failed(record, proc)
            _append_checkpoint(checkpoint, record)
            results.append(record)
            continue

        try:
            manifest = build_manifest(out_store, mask_keys)
        except (
            Exception
        ) as exc:  # noqa: BLE001 - report any hashing failure, then move on
            record["status"] = "hash_failed"
            record["error"] = str(exc)
            print(f"[hash_failed] {name}: {exc}")
            _append_checkpoint(checkpoint, record)
            results.append(record)
            continue

        manifest_path.write_text(json.dumps(manifest, sort_keys=True, indent=2))
        record["status"] = "ok"
        print(f"[ok] {name} {wall_seconds}s")
        _append_checkpoint(checkpoint, record)
        results.append(record)

    ok = sum(1 for rec in results if rec.get("status") == "ok")
    failed = sum(1 for rec in results if rec.get("status") in ("failed", "hash_failed"))
    skipped = sum(1 for rec in results if rec.get("status") == "skipped")
    print(f"run {label!r}: {ok} ok, {failed} failed, {skipped} skipped (not found)")

    return results


# ---------------------------------------------------------------------------
# ``compare``
# ---------------------------------------------------------------------------


def _label_dataset_names(scratch: Path, label: str) -> set:
    label_dir = scratch / label
    if not label_dir.is_dir():
        return set()
    return {
        path.name[: -len(".hashes.json")] for path in label_dir.glob("*.hashes.json")
    }


def _load_manifest(scratch: Path, label: str, name: str) -> dict:
    path = scratch / label / f"{name}.hashes.json"
    return json.loads(path.read_text())


def _json_text(value: Any) -> str:
    """Sorted JSON text: ``1``, ``1.0`` and ``true`` stay three values.

    Compared as parsed Python values they are one -- ``1 == 1.0 == True``
    -- so a change of JSON type in an attr passed as no change.
    """
    return json.dumps(value, sort_keys=True, default=str)


def _diff_files(name: str, manifest_a: dict, manifest_b: dict) -> int:
    """Diff the non-Zarr files and the consolidated metadata block."""
    sections = ("files", "consolidated")
    if any(
        key not in manifest for key in sections for manifest in (manifest_a, manifest_b)
    ):
        # A manifest written before these were recorded cannot vouch for
        # them; an absent section is not an empty one.
        print(f"{name} <files>: a manifest predates file hashing; run it again")
        return 1
    differences = 0
    if _json_text(manifest_a["consolidated"]) != _json_text(manifest_b["consolidated"]):
        print(f"{name} <consolidated metadata>")
        differences += 1
    files_a, files_b = manifest_a["files"], manifest_b["files"]
    for path in sorted(set(files_a) | set(files_b)):
        if path not in files_a:
            print(f"{name} {path}: only in second run")
        elif path not in files_b:
            print(f"{name} {path}: only in first run")
        elif files_a[path] != files_b[path]:
            print(f"{name} {path} hash")
        else:
            continue
        differences += 1
    return differences


def _diff_manifests(name: str, manifest_a: dict, manifest_b: dict) -> int:
    differences = 0

    group_paths = sorted(set(manifest_a["groups"]) | set(manifest_b["groups"]))
    for path in group_paths:
        element = path or "<root>"
        if path not in manifest_a["groups"]:
            print(f"{name} {element} attrs: only in second run")
            differences += 1
        elif path not in manifest_b["groups"]:
            print(f"{name} {element} attrs: only in first run")
            differences += 1
        elif _json_text(manifest_a["groups"][path]) != _json_text(
            manifest_b["groups"][path]
        ):
            print(f"{name} {element} attrs")
            differences += 1

    array_paths = sorted(set(manifest_a["arrays"]) | set(manifest_b["arrays"]))
    for path in array_paths:
        if path not in manifest_a["arrays"]:
            print(f"{name} {path}: only in second run")
            differences += 1
            continue
        if path not in manifest_b["arrays"]:
            print(f"{name} {path}: only in first run")
            differences += 1
            continue
        entry_a, entry_b = manifest_a["arrays"][path], manifest_b["arrays"][path]
        if entry_a["shape"] != entry_b["shape"]:
            print(f"{name} {path} shape")
            differences += 1
        elif entry_a["dtype"] != entry_b["dtype"]:
            print(f"{name} {path} dtype")
            differences += 1
        elif entry_a["sha256"] != entry_b["sha256"]:
            print(f"{name} {path} hash")
            differences += 1
        elif _json_text(entry_a.get("metadata")) != _json_text(entry_b.get("metadata")):
            print(f"{name} {path} metadata")
            differences += 1

    return differences + _diff_files(name, manifest_a, manifest_b)


def _mixed_settings(label: str, names: set, records: dict) -> int:
    """1 when the label's manifests were made under different settings.

    ``run`` only revisits the names it is given, so re-running a label
    for one dataset after a code change leaves the others' manifests
    from before it. Compared as one label, the old stores would answer
    for the new code. The input and its fingerprint are left out: they differ by name.
    """
    keys = [key for key in _SETTINGS_KEYS if not key.startswith("input")]
    groups: dict[str, list] = {}
    for name in sorted(names):
        record = records.get(name, {})
        made_under = _json_text({key: record.get(key) for key in keys})
        groups.setdefault(made_under, []).append(name)
    if len(groups) < 2:
        return 0
    listing = "; ".join(", ".join(group) for group in groups.values())
    print(
        f"label {label!r}: its manifests were made under {len(groups)} different "
        f"settings ({listing}); run the label again with --force"
    )
    return 1


def _vouched_names(scratch: Path, label: str, records: dict) -> tuple:
    """The label's manifest names whose latest checkpoint record is ``ok``,
    and how many manifests had no such record.

    A manifest is only as good as the record that says which run made it:
    one whose record is missing (a deleted checkpoint) or says the latest
    attempt failed or was skipped is reported, not compared.
    """
    vouched, unvouched = set(), 0
    for name in sorted(_label_dataset_names(scratch, label)):
        status = records.get(name, {}).get("status")
        if status == "ok":
            vouched.add(name)
        else:
            print(f"{name} <unvouched>: the manifest in {label!r} has no ok record")
            unvouched += 1
    return vouched, unvouched


def compare_labels(label_a: str, label_b: str, scratch: Path) -> int:
    """Diff two labels: 0 when every manifest matches, 1 on any difference,
    2 when neither label holds a manifest or both name one folder, so
    there was nothing to compare.

    Besides the stores themselves, four things count as a difference: a
    name whose latest attempt failed or was skipped and that has a store
    in neither label, a label whose manifests were made under different
    settings, a manifest whose latest record is not ``ok``, and a
    manifest from before file hashing.
    """
    # By folder, not by name: "main" and "Main" (Windows) or "main/" are one
    # folder, and a label compared with itself matches by construction.
    folder_a, folder_b = (scratch / label_a).resolve(), (scratch / label_b).resolve()
    same = folder_a == folder_b or (
        folder_a.exists() and folder_b.exists() and folder_a.samefile(folder_b)
    )
    if same:
        print(f"compare {label_a!r} vs {label_b!r}: one label; give two different ones")
        return 2

    records_a = _latest_by_name(
        _load_checkpoint(scratch / label_a / "checkpoint.jsonl")
    )
    records_b = _latest_by_name(
        _load_checkpoint(scratch / label_b / "checkpoint.jsonl")
    )
    differences = 0
    names_a, unvouched_a = _vouched_names(scratch, label_a, records_a)
    names_b, unvouched_b = _vouched_names(scratch, label_b, records_b)
    differences += unvouched_a + unvouched_b
    all_names = sorted(names_a | names_b)

    for label, names in ((label_a, names_a), (label_b, names_b)):
        if not (scratch / label).is_dir():
            print(f"label {label!r}: no such label under {scratch}")
        elif not names:
            print(f"label {label!r}: no manifest under {scratch / label}")

    # A name that failed or was skipped, with a store in neither label, was
    # not compared, and must not drop out of the count unseen.
    not_compared = ("failed", "hash_failed", "skipped")
    for label, records in ((label_a, records_a), (label_b, records_b)):
        for name, record in sorted(records.items()):
            status = record.get("status")
            if name not in all_names and status in not_compared:
                print(f"{name} <{status}>: in {label!r}, and no store in either label")
                differences += 1

    if not all_names:
        print(
            f"compare {label_a!r} vs {label_b!r}: nothing to compare, "
            "so nothing is shown to be identical"
        )
        return 2

    differences += _mixed_settings(label_a, names_a, records_a)
    differences += _mixed_settings(label_b, names_b, records_b)
    for name in all_names:
        if name not in names_a:
            print(f"{name} <missing>: present only in {label_b!r}")
            differences += 1
            continue
        if name not in names_b:
            print(f"{name} <missing>: present only in {label_a!r}")
            differences += 1
            continue
        manifest_a = _load_manifest(scratch, label_a, name)
        manifest_b = _load_manifest(scratch, label_b, name)
        differences += _diff_manifests(name, manifest_a, manifest_b)

    print(
        f"compare {label_a!r} vs {label_b!r}: {len(all_names)} dataset(s) compared, "
        f"{differences} difference(s)"
    )
    return 1 if differences else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_run_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="store_identity.py run",
        description="Convert named datasets once per subprocess and hash the resulting stores.",
    )
    parser.add_argument(
        "--label",
        required=True,
        help="Name for this run, e.g. 'main' or a branch name.",
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=REPO_ROOT,
        help="Checkout to convert with (default: this repo).",
    )
    parser.add_argument(
        "--interpreter", default=sys.executable, help="Python to run `-m thyra` with."
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=list(DEFAULT_DATASETS),
        help="Dataset names to convert.",
    )
    parser.add_argument(
        "--scratch",
        type=Path,
        default=Path(tempfile.gettempdir()) / "thyra-store-identity",
        help="Root scratch directory for stores and checkpoints.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Reconvert even if the checkpoint already has an ok run under the "
            "same settings (a run under other settings is always redone)."
        ),
    )
    parser.add_argument(
        "--mask",
        action="append",
        default=[],
        metavar="KEY",
        help=(
            "Extra attrs key(s) to drop before hashing. For one dataset read "
            "from two folders, mask source_path, ibd_file, data_path, "
            "database_path, binary_file and mis_file."
        ),
    )
    parser.add_argument(
        "thyra_options", nargs="*", help="Passed through to `thyra` after a literal --."
    )
    return parser


def _build_compare_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="store_identity.py compare",
        description="Diff two labelled runs' hash manifests.",
    )
    parser.add_argument("label_a")
    parser.add_argument("label_b")
    parser.add_argument(
        "--scratch",
        type=Path,
        default=Path(tempfile.gettempdir()) / "thyra-store-identity",
        help="Root scratch directory the labels were written under.",
    )
    return parser


def main(argv: Optional[list] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in ("run", "compare"):
        print("usage: store_identity.py {run,compare} ...", file=sys.stderr)
        return 2

    command, rest = argv[0], argv[1:]
    if command == "run":
        args = _build_run_parser().parse_args(rest)
        mask_keys = DEFAULT_MASK_KEYS | set(args.mask)
        try:
            results = run_datasets(
                label=args.label,
                datasets=args.datasets,
                source=args.source,
                interpreter=args.interpreter,
                scratch=args.scratch,
                thyra_options=args.thyra_options,
                mask_keys=mask_keys,
                force=args.force,
            )
        except SourceNotUsed as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        failed = any(rec.get("status") in ("failed", "hash_failed") for rec in results)
        return 1 if failed else 0

    args = _build_compare_parser().parse_args(rest)
    return compare_labels(args.label_a, args.label_b, args.scratch)


if __name__ == "__main__":
    raise SystemExit(main())
