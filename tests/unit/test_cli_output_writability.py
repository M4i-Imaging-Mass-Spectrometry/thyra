"""The output-writability refusal, on a directory that really is denied.

``os.access(path, os.W_OK)`` consults the read-only *attribute* on
Windows, not the ACL, so it answers ``True`` for a directory the user
cannot write. Measured on a directory denied ``(WD,AD)`` via ``icacls``:

    os.access(W_OK)   -> True
    st_mode           -> 0o40777
    mkdir             -> PermissionError, Access is denied
    open(w)           -> PermissionError errno 13

So the refusal could never fire on the platform this project is
developed on. The failure surfaced instead from the CSC scratch
allocation deep inside the conversion, as "Error during conversion" plus
a traceback, after the metadata scan (issue #312).

A monkeypatched probe would pass against the old code too and would
prove nothing, so this uses a real ACL. It is skipped off Windows, where
``os.access`` is not broken and there is nothing to guard.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import click
import pytest

from thyra.__main__ import _validate_output_path

pytestmark = pytest.mark.skipif(
    os.name != "nt", reason="os.access(W_OK) is only a no-op on Windows"
)


@pytest.fixture
def denied_dir(tmp_path: Path):
    """A subdirectory the current user may not create anything in.

    Never ``tmp_path`` itself -- pytest's own teardown would trip on it.
    """
    target = tmp_path / "denied"
    target.mkdir()
    user = os.environ.get("USERNAME", "")
    result = subprocess.run(
        ["icacls", str(target), "/deny", f"{user}:(WD,AD)"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.skip(f"could not deny write on {target}: {result.stderr.strip()}")
    try:
        yield target
    finally:
        subprocess.run(
            ["icacls", str(target), "/remove:d", user],
            capture_output=True,
            text=True,
        )


def test_the_platform_assumption_still_holds(denied_dir):
    """If this ever fails, os.access has been fixed and the probe can go."""
    assert os.access(denied_dir, os.W_OK) is True


def test_an_unwritable_destination_is_refused(denied_dir):
    with pytest.raises(click.BadParameter, match="not writable"):
        _validate_output_path(denied_dir / "out.zarr")


def test_the_message_names_the_directory(denied_dir):
    with pytest.raises(click.BadParameter) as excinfo:
        _validate_output_path(denied_dir / "out.zarr")
    assert str(denied_dir) in str(excinfo.value)


def test_a_writable_destination_is_accepted(tmp_path):
    _validate_output_path(tmp_path / "out.zarr")


def test_the_probe_leaves_nothing_behind(tmp_path):
    before = sorted(p.name for p in tmp_path.iterdir())
    _validate_output_path(tmp_path / "out.zarr")
    assert sorted(p.name for p in tmp_path.iterdir()) == before


def test_a_nested_destination_probes_the_nearest_existing_ancestor(denied_dir):
    """The store's parents do not exist yet; the refusal is about where
    they would have to be created."""
    with pytest.raises(click.BadParameter, match="not writable"):
        _validate_output_path(denied_dir / "a" / "b" / "out.zarr")
