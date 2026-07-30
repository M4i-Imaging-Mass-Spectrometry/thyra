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
        "Zarr metadata rename to %s still blocked after %d attempts. "
        "Something is holding the file open -- check for a Python session, "
        "napari, or a notebook with the store loaded.",
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
