"""Bounded retry for Zarr's atomic metadata writes on Windows.

Zarr writes every metadata document by streaming it to a
``zarr.<32-hex>.partial`` sibling, closing that file, and then renaming it
over the real name (``zarr/storage/_local.py::_atomic_write``). On POSIX
the rename is atomic and silently replaces whatever is there. On Windows
there is no atomic replace: ``MoveFileEx`` has to delete the destination
first, and if any handle is open on it at that instant the call fails with
``[WinError 5] Access is denied`` instead of waiting.

Thyra writes several metadata documents repeatedly during one conversion
-- the store root ``zarr.json`` is rewritten once per element added, and a
single conversion issues a few hundred metadata renames -- so the
destination is genuinely contended. Measured on Windows 11 with
zarr 3.1.3: across 150 conversions and roughly 64,500 metadata renames, 40
failed with WinError 5, and every one of them succeeded on an immediate
retry with no delay. Left unhandled each of those aborts the whole
conversion, which is why it surfaced as intermittent test failures and as
``Error saving SpatialData: [WinError 5] Access is denied: ...partial ->
...zarr.json`` in real runs.

The file handle is already closed before the rename, so there is nothing
for Thyra to close earlier; the contention is on the destination. A
bounded retry is the fix.

Measured again 2026-09-15, because the paragraph above and the error text
below both blamed a handle, and a handle is not what this is. With **pure
stdlib, one thread, and zarr not imported**, a replace onto a destination
that already exists fails 720 times in 20,000 (3.6%). At the instant of
failure both paths still open (``rb`` and ``r+b`` on the destination,
``rb`` on the ``.partial``), and ``psutil.Process().open_files()`` shows
this process holding neither. Four arms then separate the conditions, 4,000
replaces each:

===================================================  ==============
arm                                                  failures
===================================================  ==============
A same destination name, destination exists          221 (5.53%)
B fresh destination name, destination exists         0
C same destination name, destination absent          0
D same name, exists, 20 ms before the replace        3 of 500 (0.6%)
===================================================  ==============

So the destination merely existing is not enough (B is clean), and the
write pattern alone is not enough (C is clean). It takes **replacing the
same name over and over**, and it eases when the freshly written source is
given a moment (D). That is the signature of NTFS *file system
tunnelling*, the cache that holds a name for a short while after it is
deleted or renamed away -- and it is the explanation Microsoft's own thread
on this reaches, having first ruled antivirus out by reproducing with it
disabled and found the failure gone on Windows Server 2025:
https://learn.microsoft.com/en-us/answers/questions/5559596/

That matters because the antivirus story is the one everybody downstream
tells -- zarr's own v2 fix is literally titled "Make DirectoryStore
__setitem__ resilient against antivirus file locking" (#698) -- and it
predicts arm B would fail too. It does not. The retry is still exactly the
right fix; only the explanation was wrong.

The same 20,000 replaces with the retry below: 0 unrecovered, 498 needing
one, and never more than four attempts.

Only a key written TWICE is exposed, which is why first writes never fail.
CI is not affected and never has been -- 0 ``test (windows-latest)``
failures in 60 runs. This is a developer-machine and end-user-machine fix.

Deferring to a fixed Zarr
-------------------------

Zarr has since taken this retry upstream, in zarr-developers/zarr-python#4358
(merged 2026-09-15, shipped in **v3.4.0** the same day). On a Zarr that has
it, ``install_windows_atomic_write_retry`` patches nothing and Zarr's own
retry runs. This project does not pin such a Zarr yet -- the ceiling is still
``<3.2`` -- so on the pinned Zarr the patch below is what runs.

That check is not politeness, it is necessary. This module does not *wrap*
``_atomic_write``, it **substitutes** a copy of it -- it has to, because the
rename to retry happens inside that function, after the caller has already
written into the file it yielded. Substituting means it would go on overriding
whatever else upstream does to that function, and something else is already
queued: zarr-developers/zarr-python#4173 moves the ``.partial`` files out of
the store into a configurable temporary directory. With the patch installed
unconditionally, Thyra would silently keep writing them beside the
destination, and nothing would fail to say so.

Deferring costs only the exhaustion message logged below. The delays and the
retried error codes are the same on both sides -- upstream's came from the
measurements above. Re-checked against the released v3.4.0 rather than the
merge commit, because #4358 was merged as-is with polishing promised in a
follow-up: same codes (5, 32) and the same 1/5/20/50/200 ms spacing over six
attempts. Upstream sleeps *after* a failure and makes its last attempt outside
the loop; this module sleeps *before* each attempt with a leading zero. Same
schedule, written the other way round.

Removing this module outright is thyra#341, which needs the ``zarr`` ceiling
in pyproject.toml raised past v3.4.0 first.
"""

import contextlib
import logging
import sys
import time
import uuid
from pathlib import Path
from typing import Any, BinaryIO, Callable, Iterator, Literal

logger = logging.getLogger(__name__)

# Codes MoveFileEx reports when the destination cannot be replaced *right
# now*, as opposed to a genuine permission problem. Both are worth a
# retry; anything else propagates untouched.
_TRANSIENT_WINERRORS = frozenset(
    {
        5,  # ERROR_ACCESS_DENIED
        32,  # ERROR_SHARING_VIOLATION
    }
)

# Delay before each successive attempt. The leading 0.0 is the original
# attempt. Every observed failure cleared on the first retry, so the later
# delays are headroom for a busier machine rather than an expected path.
_RETRY_DELAYS = (0.0, 0.001, 0.005, 0.02, 0.05, 0.2)

# The helper zarr-python#4358 added to zarr.storage._local to hold its own
# version of this retry. Its presence is what tells us a Zarr needs nothing
# from us. A name is a weak hinge, so it is chosen to be the sturdiest one
# available: Zarr's own tests import it from there by name, so a rename
# cannot pass upstream CI silently. If it is renamed anyway the check simply
# does not fire and we patch as before, which is where we are today.
_UPSTREAM_RETRY_HELPER = "_move_with_retry"

# Whether the one-time decision below has been taken -- to patch, or to leave
# a Zarr that retries for itself alone. Either way there is nothing to redo.
_installed = False


def _is_transient(exc: OSError) -> bool:
    """Whether the destination was merely busy rather than forbidden."""
    return getattr(exc, "winerror", None) in _TRANSIENT_WINERRORS


def _move_with_retry(
    tmp_path: Path,
    path: Path,
    move: Callable[[Path, Path], Any],
) -> None:
    """Run ``move(tmp_path, path)``, retrying while the target is busy.

    Args:
        tmp_path: The freshly written temporary file.
        path: The destination to move it onto.
        move: The move to perform, as Zarr would have performed it.

    Raises:
        OSError: The last error seen, if every attempt is exhausted, or
            immediately for any error that is not a transient lock.
    """
    last_error: OSError

    for attempt, delay in enumerate(_RETRY_DELAYS):
        if delay:
            time.sleep(delay)
        try:
            move(tmp_path, path)
            if attempt:
                logger.debug(
                    "Zarr metadata rename to %s succeeded on attempt %d "
                    "after a transient Windows lock",
                    path,
                    attempt + 1,
                )
            return
        except OSError as e:
            if not _is_transient(e):
                raise
            last_error = e

    logger.error(
        "Zarr metadata rename to %s still blocked after %d attempts. The "
        "transient form of this clears within a millisecond, so lasting "
        "this long means something really does have the store open: look "
        "for another Python session, napari, or a notebook holding it.",
        path,
        len(_RETRY_DELAYS),
    )
    raise last_error


def make_atomic_write_with_retry(
    safe_move: Callable[[Path, Path], Any],
) -> Callable[..., Any]:
    """Build a replacement for Zarr's ``_atomic_write`` with retries.

    Mirrors ``zarr.storage._local._atomic_write`` (zarr 3.0-3.1, the range
    pinned in pyproject.toml) with the final move retried while the
    destination is busy.

    Args:
        safe_move: Zarr's non-clobbering move, used for exclusive writes.
    """

    @contextlib.contextmanager
    def _atomic_write_with_retry(
        path: Path,
        mode: Literal["r+b", "wb"],
        exclusive: bool = False,
    ) -> Iterator[BinaryIO]:
        tmp_path = path.with_suffix(f".{uuid.uuid4().hex}.partial")
        try:
            with tmp_path.open(mode) as f:
                yield f
            # Zarr uses a non-replacing move when it must not clobber an
            # existing key, and a replacing one otherwise. Keep both.
            if exclusive:
                _move_with_retry(tmp_path, path, safe_move)
            else:
                _move_with_retry(tmp_path, path, lambda src, dst: src.replace(dst))
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise

    _atomic_write_with_retry._thyra_retry_wrapped = True  # type: ignore[attr-defined]
    return _atomic_write_with_retry


def install_windows_atomic_write_retry() -> None:
    """Add a bounded retry to Zarr's atomic metadata writes on Windows.

    A no-op off Windows, where the underlying rename replaces the
    destination atomically and cannot hit this failure, and a no-op if
    called more than once.

    Also a no-op on a Zarr that carries its own retry (zarr-python#4358),
    so that a Zarr new enough to have been fixed keeps its own write path
    instead of this module's copy of the old one. See the module docstring
    for why that matters more than it looks.

    If Zarr's private atomic-write helpers are not shaped as expected -- a
    newer Zarr having reworked them -- nothing is patched and the caller
    keeps Zarr's own behaviour.
    """
    global _installed

    if _installed or sys.platform != "win32":
        return

    try:
        import zarr.storage._local as zarr_local
    except ImportError:  # pragma: no cover - zarr is a hard dependency
        logger.debug("Could not import zarr.storage._local; not patching")
        return

    if hasattr(zarr_local, _UPSTREAM_RETRY_HELPER):
        _installed = True
        logger.debug(
            "This Zarr retries its own atomic metadata writes; leaving its "
            "write path untouched"
        )
        return

    original = getattr(zarr_local, "_atomic_write", None)
    safe_move = getattr(zarr_local, "_safe_move", None)

    if original is None or safe_move is None:
        logger.debug(
            "zarr.storage._local does not expose the expected atomic-write "
            "helpers; leaving Zarr's write path untouched"
        )
        return

    if getattr(original, "_thyra_retry_wrapped", False):
        _installed = True
        return

    zarr_local._atomic_write = make_atomic_write_with_retry(safe_move)
    _installed = True
    logger.debug("Installed bounded retry for Zarr atomic metadata writes")
