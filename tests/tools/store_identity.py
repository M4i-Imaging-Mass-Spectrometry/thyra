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
the whole set.

``compare`` diffs the manifests of two labelled runs and prints one line per
Zarr element that differs, then exits 1 if anything did.

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
   different skip from "not found anywhere" -- it says so.
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

What is masked before hashing
------------------------------

Two keys are dropped because they vary run to run without the store's data
changing: ``conversion_timestamp`` and ``thyra_version``. Masking works two
ways, because ``uns`` is not attrs:

- **Attrs dicts** (a group's own ``.attrs``, an array's ``.metadata``
  ``attributes`` block) are masked by key, recursively, at any depth.
- **anndata writes ``uns`` as zarr groups and arrays, not attrs** -- a
  scalar such as ``thyra_version`` is a 0-d array literally named
  ``thyra_version``, not an attrs entry -- so ``build_manifest`` also drops
  any array or group with a mask key anywhere in its path, descendants
  included, wherever it sits. And a string *array element* that parses as JSON
  (``uns`` payloads such as the ``cvParams`` / ``regions`` blocks are
  stored as JSON strings -- see ``thyra/metadata/uns_compat.py``) is
  parsed, masked the same way, and re-serialised before it is hashed
  (see ``_hash_block``).

``--mask KEY`` adds more keys to drop from both paths at once. Reach for it
when the same dataset is read from two different locations: ``--mask
source_path`` drops ``uns/essential_metadata/source_path``, the absolute
input path recorded on every table, which otherwise differs between the two
runs for a reason that has nothing to do with the data.

Run the fixture round in the acceptance checks below (the same label
converted twice) to find out whether the mask list is complete for a given
store shape -- a non-zero ``compare`` there means it is not. A blanket
ISO-8601-in-string regex was tried and dropped: it would mask a genuinely
changed ``acquisition_datetime`` (real source metadata, not a run artefact)
exactly as it masks a volatile one, and measuring against both the fixture
round and a real round found nothing left over that needed it -- the two
named keys are the whole list, confirmed rather than assumed.
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
from pathlib import Path
from typing import Any, Optional

import numpy as np
import zarr

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent

DEFAULT_MASK_KEYS = {"conversion_timestamp", "thyra_version"}

# Every tests/data/fixtures file that actually converts with --no-optical,
# checked by hand: each of tests/data/fixtures/*.imzML and synthetic_tims.d
# was run through `thyra` once. two_precision_terms.imzML is deliberately
# excluded -- it is a fixture *for* a refusal (two precision terms declared
# for the m/z array) and exits non-zero on purpose, not by defect.
DEFAULT_DATASETS = (
    "iontof_sparse",
    "mobility_continuous",
    "mobility_processed",
    "solarix_fticr",
    "synthetic_tims",
    "timstof_flex_export",
    "two_scansettings",
    "unit_nanometre",
)

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
    if path.exists():
        return path
    raise DatasetPathUnreachable(
        f"{name}: registered but its path is not reachable on this machine"
    )


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
            hits = sorted(exact.glob("*.d"))
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
        return {
            key: mask_value(item, mask_keys)
            for key, item in value.items()
            if key not in mask_keys
        }
    if isinstance(value, list):
        return [mask_value(item, mask_keys) for item in value]
    if isinstance(value, str):
        parsed = _maybe_json(value)
        if parsed is not None:
            return json.dumps(
                mask_value(parsed, mask_keys), sort_keys=True, default=str
            )
        return value
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


def build_manifest(store_path: Path, mask_keys: set) -> dict:
    """Walk the whole Zarr store: every array's shape/dtype/hash/metadata,
    and every group's (masked) attrs, including the root.

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
    """
    root = zarr.open_group(str(store_path), mode="r")
    groups: dict[str, Any] = {"": mask_value(dict(root.attrs), mask_keys)}
    arrays: dict[str, Any] = {}
    for path, node in root.members(max_depth=None):
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
    return {"groups": groups, "arrays": arrays}


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
    that name is done. ``--force`` reconverts the requested names regardless
    of a prior ok, and appends a fresh record rather than discarding the
    checkpoint's record of every other name under this label.
    """
    label_dir = scratch / label
    label_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = label_dir / "checkpoint.jsonl"

    latest = _latest_by_name(_load_checkpoint(checkpoint))
    done = {name: rec for name, rec in latest.items() if rec.get("status") == "ok"}

    source = source.resolve()
    source_rev = _git_rev(source)
    results = []

    for name in datasets:
        if not force and name in done:
            print(f"[skip-cached] {name} (label {label!r} already has an ok run)")
            results.append(done[name])
            continue

        try:
            resolved = resolve_dataset(name)
        except (AmbiguousDatasetName, DatasetPathUnreachable) as exc:
            print(f"[skip] {exc}")
            results.append({"name": name, "status": "skipped"})
            continue

        if resolved is None:
            print(
                f"[skip] {name}: not found in the registry, THYRA_DATASET_ROOTS, "
                "or tests/data/fixtures"
            )
            results.append({"name": name, "status": "skipped"})
            continue

        # Delete any manifest from a previous attempt before this one starts,
        # so a failed conversion can never leave a stale ok manifest behind
        # for `compare` to silently read.
        manifest_path = label_dir / f"{name}.hashes.json"
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
        env = dict(os.environ)
        env["PYTHONPATH"] = str(source) + (
            os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else ""
        )

        print(f"[run] {name} <- {resolved}")
        started = time.perf_counter()
        proc = subprocess.run(
            cmd, cwd=str(source), env=env, capture_output=True, text=True
        )
        wall_seconds = round(time.perf_counter() - started, 3)

        record = {
            "name": name,
            "options": thyra_options,
            "source": str(source),
            "source_rev": source_rev,
            "wall_seconds": wall_seconds,
        }

        if proc.returncode != 0:
            record["status"] = "failed"
            record["returncode"] = proc.returncode
            record["stderr_tail"] = proc.stderr[-4000:]
            print(f"[failed] {name} after {wall_seconds}s (exit {proc.returncode})")
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
        elif manifest_a["groups"][path] != manifest_b["groups"][path]:
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
        elif json.dumps(
            entry_a.get("metadata"), sort_keys=True, default=str
        ) != json.dumps(entry_b.get("metadata"), sort_keys=True, default=str):
            print(f"{name} {path} metadata")
            differences += 1

    return differences


def compare_labels(label_a: str, label_b: str, scratch: Path) -> int:
    names_a = _label_dataset_names(scratch, label_a)
    names_b = _label_dataset_names(scratch, label_b)
    all_names = sorted(names_a | names_b)

    differences = 0
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
        help="Reconvert even if the checkpoint already has an ok run.",
    )
    parser.add_argument(
        "--mask",
        action="append",
        default=[],
        metavar="KEY",
        help="Extra attrs key(s) to drop before hashing.",
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
        failed = any(rec.get("status") in ("failed", "hash_failed") for rec in results)
        return 1 if failed else 0

    args = _build_compare_parser().parse_args(rest)
    return compare_labels(args.label_a, args.label_b, args.scratch)


if __name__ == "__main__":
    raise SystemExit(main())
