"""Opening a Bruker ``.d``'s sqlite files without writing to them.

A vendor ``.d`` is the user's only copy of an acquisition, and Thyra
reads it. ``sqlite3.connect(path)`` opens read-write, which has three
consequences none of them wanted:

- **It creates side files.** A ``.tdf`` in WAL mode gets a ``-wal`` and a
  ``-shm`` written next to it on the first read. They land inside the
  vendor directory, which may be a read-only share or a directory the
  user backs up byte-for-byte.
- **It takes a lock.** DataAnalysis holding the acquisition open is the
  ordinary case in a lab, not an exotic one, and a read-write open loses
  to it.
- **It can be refused outright.** A share mounted read-only, or a
  ``.d`` on write-protected media, fails to open at all.

``mode=ro`` alone fixes the first and third but not the second: a
read-only connection still waits on the writer's lock. ``immutable=1``
is what steps past it -- it promises sqlite the file will not change
underneath, so no locking is attempted and no journal is consulted. That
promise is the right one here because Thyra never writes to a ``.d`` and
an acquisition being actively rewritten is not a dataset to convert.

The URI is built rather than interpolated because of UNC paths, which is
what a mapped network drive resolves to on Windows and where this lab's
data lives.
"""

import sqlite3
from pathlib import Path
from typing import Union
from urllib.parse import quote

__all__ = ["read_only_uri", "open_read_only"]


def read_only_uri(path: Union[str, Path], immutable: bool = True) -> str:
    """The sqlite URI that opens ``path`` strictly read-only.

    A UNC path -- which is what a mapped network drive resolves to on
    Windows -- starts with ``//server/share``; inside a ``file:`` URI
    that reads as an authority, which sqlite rejects ("invalid uri
    authority"). Doubling the leading slashes leaves the authority empty
    and the path intact, so ``file:////server/share/...`` opens where
    ``file://server/share/...`` does not. Drive-letter and POSIX paths
    are unaffected.

    Args:
        path: The database file.
        immutable: Add ``immutable=1``, which makes the open skip
            locking entirely. Correct for a vendor file Thyra only ever
            reads; leave it off for a database something else may be
            legitimately writing.

    Returns:
        A ``file:`` URI for :func:`sqlite3.connect` with ``uri=True``.
    """
    posix = Path(path).as_posix()
    if posix.startswith("//"):
        posix = "//" + posix
    suffix = "&immutable=1" if immutable else ""
    return f"file:{quote(posix)}?mode=ro{suffix}"


def open_read_only(
    path: Union[str, Path],
    *,
    immutable: bool = True,
    timeout: float = 30.0,
    check_same_thread: bool = True,
) -> sqlite3.Connection:
    """Connect to a vendor sqlite file read-only.

    Args:
        path: The database file.
        immutable: See :func:`read_only_uri`.
        timeout: Seconds to wait for a lock. Only reachable with
            ``immutable=False``; an immutable open takes no locks.
        check_same_thread: Passed through; the TDF reader shares one
            connection across threads and sets this False.

    Returns:
        An open connection. The caller owns it -- wrap it in
        ``contextlib.closing``, because ``with sqlite3.connect(...)``
        commits a transaction and does **not** close the handle.
    """
    return sqlite3.connect(
        read_only_uri(path, immutable=immutable),
        uri=True,
        timeout=timeout,
        check_same_thread=check_same_thread,
    )
