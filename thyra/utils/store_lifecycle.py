"""Where a store lives while it is written, and where it goes if that fails.

The output name holds a finished store or nothing, however the run ends.
A conversion writes its store under a temporary name in the output's own
folder, ``<name>.incomplete-<8 hex>``, and renames it to the output name
as its last step. Before, the store was written in place, so a run killed
by the operating system (an out-of-memory kill, a native crash) left a
partial store under the output name: Python had no chance to move it
aside, the retry was refused because the output existed, and ``thyra
validate`` blamed an old Thyra version (issue #525).

A failure Python does see still moves the partial store to
``<name>.failed`` (``.failed2`` and so on), as it did before (issues #245,
#293). A killed run now leaves only the ``.incomplete-<id>`` folder and its
``.thyra_*`` scratch beside the output, for the user to delete, and the
same command can be run again at once.
"""

import logging
import re
import time
from pathlib import Path
from typing import Callable
from uuid import uuid4

logger = logging.getLogger(__name__)

#: What follows the output name in the temporary name.
INCOMPLETE_INFIX = ".incomplete-"

_INCOMPLETE_NAME = re.compile(r"\.incomplete-[0-9a-f]{8}$")

#: How long the final rename keeps trying. On Windows a rename onto a new
#: name now and then fails with ``WinError 5`` while an indexer or a
#: virus scanner still holds a file the write just closed; it clears
#: within a fraction of a second.
_RENAME_ATTEMPTS = 10
_RENAME_PAUSE_S = 0.2


def incomplete_store_path(output: Path) -> Path:
    """A fresh temporary name for ``output``'s store, in its own folder."""
    return output.with_name(f"{output.name}{INCOMPLETE_INFIX}{uuid4().hex[:8]}")


def is_incomplete_store_name(name: str) -> bool:
    """Whether ``name`` is a temporary store name :func:`incomplete_store_path` made."""
    return _INCOMPLETE_NAME.search(name) is not None


def retry_on_os_error(action: Callable[[], object]) -> None:
    """Run ``action``, trying again for a moment when it raises ``OSError``.

    The last error is raised when every attempt failed.
    """
    for attempt in range(_RENAME_ATTEMPTS):
        try:
            action()
            return
        except OSError:
            if attempt == _RENAME_ATTEMPTS - 1:
                raise
            time.sleep(_RENAME_PAUSE_S)


def quarantine_partial_store(partial: Path, output: Path) -> None:
    """Move a partially written store to a sibling ``<output>.failed`` path.

    ``partial`` is where the store was written: the temporary name of a
    conversion, or ``output`` itself for a caller that wrote there. The
    output name stays free for a retry while the partial store remains
    available for diagnosis. A no-op when nothing is at ``partial``.
    """
    if not partial.exists():
        return

    quarantine = output.with_name(f"{output.name}.failed")
    attempt = 1
    while quarantine.exists():
        attempt += 1
        quarantine = output.with_name(f"{output.name}.failed{attempt}")

    try:
        retry_on_os_error(lambda: partial.rename(quarantine))
    except OSError as e:
        # A store under its temporary name does not block a retry; one at
        # the output name does.
        then = "before retrying" if partial == output else "when you can"
        logger.error(
            "The incomplete output was left at %s because it could not be "
            "moved aside (%s). It will not open with "
            "spatialdata.read_zarr(); delete it %s.",
            partial,
            e,
            then,
        )
        return

    logger.error(
        "The incomplete output was moved to %s. It will not open with "
        "spatialdata.read_zarr(); delete it once you no longer need it.",
        quarantine,
    )
