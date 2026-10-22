"""The one rule for every file a Thyra command writes beside the store.

``--log-file``, ``thyra metadata -o``, ``thyra export-metaspace -o`` (and
its default beside the store) and ``thyra-check-ontology --output`` each
opened their file without a check. A path in a folder that did not exist
yet, or a path that was itself a folder, ended in a traceback -- after the
metadata read or the whole ontology check had already run, and for
``--log-file`` before the conversion it was meant to record (issue #525).

They now follow the store path's rule. The path is checked before any work
starts: one that is an existing folder, or that lies under a regular file,
is refused with one line, which each command turns into its own usage
error (exit 2). A folder that does not exist yet is created when the file
is written, as it is for the store, and not before: a run that fails
first leaves no empty folder behind.
"""

from pathlib import Path
from typing import Optional


def output_file_problem(path: Path) -> Optional[str]:
    """Why a file cannot be written at ``path``, or ``None``.

    Only the two shapes no folder creation can fix are refused: ``path``
    is a folder, or the deepest part of it that exists is a regular file.
    """
    if path.is_dir():
        return (
            f"Cannot write {path}: it is a folder. Give the path of a file "
            "to write, for example a name inside that folder."
        )
    ancestor = next((p for p in path.parents if p.exists()), None)
    if ancestor is not None and not ancestor.is_dir():
        return (
            f"Cannot write {path}: {ancestor} is a file, not a folder, so "
            "nothing can be created under it."
        )
    return None


def make_parent_folder(path: Path) -> None:
    """Create the folder ``path`` is written into, if it is missing.

    Called right before the write, once :func:`output_file_problem` has
    passed the path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
