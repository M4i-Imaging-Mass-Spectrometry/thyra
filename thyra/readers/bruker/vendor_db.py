"""Opening a Bruker ``.d``'s sqlite files without writing to them.

A vendor ``.d`` is the user's only copy of an acquisition, and Thyra
reads it. ``sqlite3.connect(path)`` opens read-write, which has three
consequences none of them wanted:

- **It creates side files.** A ``.tdf`` in WAL mode gets a ``-wal`` and a
  ``-shm`` written next to it on the first read. They land inside the
  vendor directory, which may be a read-only share or a directory the
  user backs up byte-for-byte. A *clean* read-write close checkpoints
  and removes them again; several of these opens leaked their handle
  (``with sqlite3.connect(...)`` commits without closing), so they
  stayed.
- **It takes a lock.** DataAnalysis holding the acquisition open is the
  ordinary case in a lab, not an exotic one, and a read-write open loses
  to it.
- **It can be refused outright.** A share mounted read-only, or a
  ``.d`` on write-protected media, fails to open at all.

``mode=ro`` alone fixes the first and third but not the second: a
read-only connection still waits on the writer's lock. ``immutable=1``
is what steps past it -- it promises sqlite the file will not change
underneath, so no locking is attempted and no journal is consulted.

**That promise has to be checked, not assumed.** "No journal is
consulted" includes the write-ahead log, so on a database whose ``-wal``
holds committed rows, ``immutable=1`` returns the state before them --
silently, with no error to notice. Measured, 50 rows checkpointed and 50
more committed into the ``-wal`` of a file a writer still holds open::

    truth                : 100
    mode=ro              : 100
    mode=ro&immutable=1  :  50

A short frame count read as if it were the whole acquisition is the
defect this module exists to remove, not one to introduce, so the
immutable promise is made only when there is a **non-empty** ``-wal``
beside the file to contradict it. (A ``wal_checkpoint(TRUNCATE)`` leaves
a 0-byte one behind; it holds nothing, and downgrading on it would be
self-perpetuating, since a downgraded read leaves its own side files.)

**The downgrade costs something, and it is the lesser cost, not a free
one.** Plain ``mode=ro`` reads the log correctly, but it is not
side-file-free and it is not always able to open. Measured:

    WAL db, no side files   mode=ro    after = [tdf, -shm, -wal]
    -wal present, no -shm   mode=ro    after = [tdf, -shm, -wal]
    ...same, dir denied     mode=ro    OperationalError, cannot open
    ...same, dir denied     immutable  reads, 50 of 50 -- but stale

So on a read-only share holding a database with a live ``-wal``, this
refuses where an immutable open would have answered. That is the right
way round: sqlite cannot read a WAL database without the ``-shm`` it
cannot create there, so the alternative is not a correct answer, it is a
stale one presented as current. Failing is what lets the user close the
writer, or copy the ``.d`` somewhere writable, and get the real number.

That leaves ``immutable=1`` covering the case with no live log, which is
the ordinary one. Measured on a WAL-mode database with no side files
present::

    read-write           during=[tdf, -shm, -wal]  after=[tdf]
    mode=ro              during=[tdf, -shm, -wal]  after=[tdf, -shm, -wal]
    mode=ro&immutable=1  during=[tdf]              after=[tdf]

so of the two read-only spellings, plain ``mode=ro`` is the one that
leaves litter behind.

A rollback journal is not probed for. ``immutable=1`` ignores a hot
``-journal`` the same way it ignores a ``-wal``, and the same reasoning
would apply; it is not handled here because nothing has shown a vendor
``.d`` in that state, and a probe nobody can trigger is a claim nobody
can check.

The URI is built rather than interpolated because of UNC paths, which is
what a mapped network drive resolves to on Windows and where this lab's
data lives.
"""

import logging
import sqlite3
from pathlib import Path
from typing import Union
from urllib.parse import quote

logger = logging.getLogger(__name__)

__all__ = ["read_only_uri", "open_read_only"]


def _has_live_wal(path: Path) -> bool:
    """Whether a write-ahead log beside ``path`` holds anything.

    Existence is not enough: ``PRAGMA wal_checkpoint(TRUNCATE)`` leaves a
    0-byte ``-wal`` behind, which contradicts nothing. Treating that as
    live would also be self-perpetuating, because the downgraded read it
    forces leaves its own ``-wal`` and ``-shm`` for the next one to find.
    """
    try:
        return path.with_name(path.name + "-wal").stat().st_size > 0
    except OSError:
        # Absent, or a directory we cannot stat. Either way there is no
        # log we can show to be live, so the immutable promise stands.
        return False


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
        immutable: Ask for ``immutable=1``, which makes the open skip
            locking entirely. It is **downgraded to plain** ``mode=ro``
            when a ``-wal`` sits beside the file, because immutability
            also skips the write-ahead log and would return the database
            as it stood before the log's committed rows. See the module
            docstring for the measurement.

    Returns:
        A ``file:`` URI for :func:`sqlite3.connect` with ``uri=True``.
    """
    path = Path(path)
    if immutable and _has_live_wal(path):
        logger.debug(
            "%s has a non-empty -wal beside it, so it is opened read-only but "
            "not immutable: an immutable open would skip the log's committed "
            "rows.",
            path,
        )
        immutable = False

    posix = path.as_posix()
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
        immutable: See :func:`read_only_uri`, including when it is
            downgraded.
        timeout: Seconds to wait for a lock. Reached whenever the open
            is not immutable, which includes the downgraded case.
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
