"""Two conventions this suite keeps, read off the source rather than run.

Both exist because the thing they forbid fails *silently*. A ``caplog``
assertion on a Thyra record does not raise when it stops working; it goes
empty, and ``assert "..." not in caplog.text`` on an empty capture passes.
A test that assigns ``sys.argv`` does not fail either -- some later test
does, with an argv it never set.

Deliberately AST-only: no ``import thyra``, no fixtures, no conversion, so
this runs anywhere and reports a filename and a line number rather than a
downstream failure in another file.

Scope is ``test_*.py`` under ``tests/``, which by construction excludes
``tests/conftest.py`` -- the one place both forbidden things belong, since
the fixture that restores ``sys.argv`` has to write it.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Iterator, List, Tuple

TESTS_DIR = Path(__file__).resolve().parent.parent


def _test_files() -> List[Path]:
    return sorted(TESTS_DIR.rglob("test_*.py"))


def _parse(path: Path) -> ast.Module:
    # utf-8-sig, not utf-8. No tracked file carries a BOM any more --
    # tests/unit/metadata/schema/test_cli_metadata.py did until #307 stripped
    # it, and that issue also added the fix-byte-order-marker hook and a test
    # over the index to keep it that way. The tolerant read stays anyway,
    # because the failure it absorbs is out of proportion to its cost: one BOM
    # slipping past both makes ast.parse raise "SyntaxError: invalid
    # non-printable character U+FEFF" here, and this module would report that
    # instead of the convention it exists to check.
    return ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))


def _rel(path: Path) -> str:
    return path.relative_to(TESTS_DIR.parent).as_posix()


def _argv_writes(tree: ast.Module) -> Iterator[Tuple[int, str]]:
    """Yield ``(lineno, spelling)`` for every write to ``sys.argv``.

    Three spellings, because one of them is the idiom the restoring
    fixture itself uses and so the one a reader is most likely to copy:

        sys.argv = [...]      Assign to an Attribute
        sys.argv[:] = [...]   Assign to a Subscript of that Attribute
        sys.argv += [...]     AugAssign to either
    """

    def is_sys_argv(node: ast.expr) -> bool:
        if isinstance(node, ast.Subscript):
            node = node.value
        return (
            isinstance(node, ast.Attribute)
            and node.attr == "argv"
            and isinstance(node.value, ast.Name)
            and node.value.id == "sys"
        )

    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AugAssign):
            targets = [node.target]
        else:
            continue
        for target in targets:
            if is_sys_argv(target):
                yield node.lineno, ast.unparse(target)


def test_no_test_takes_caplog() -> None:
    """Thyra's own log records are captured with ``thyra_logs``.

    ``setup_logging`` sets ``propagate = False`` on the ``thyra`` logger
    and clears its handlers, and that is process-global: after any test in
    the session has invoked the CLI, caplog's root handler never sees
    another Thyra record. Which pytest version is installed decides
    whether that shows up -- 9.x walks ``loggerDict`` and attaches to
    non-propagating loggers, 8.x does not, and neither catches a logger
    that becomes non-propagating inside the capture block.

    ``thyra_logs`` (tests/conftest.py) attaches to the named logger
    instead, so propagation cannot take it away. It carries ``.text``,
    ``.messages`` and the records themselves, so a migration is a rename.
    """
    offenders = []
    for path in _test_files():
        for node in ast.walk(_parse(path)):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if not node.name.startswith("test_"):
                continue
            args = node.args
            names = [a.arg for a in args.args] + [a.arg for a in args.kwonlyargs]
            if "caplog" in names:
                offenders.append(f"{_rel(path)}:{node.lineno} {node.name}")

    assert not offenders, (
        "These tests take pytest's caplog. Use the thyra_logs fixture from "
        "tests/conftest.py instead -- caplog cannot see a Thyra record once "
        "any test in the session has invoked the CLI:\n  " + "\n  ".join(offenders)
    )


def test_no_test_assigns_sys_argv() -> None:
    """Command-line arguments are set with ``monkeypatch``.

    A bare assignment is never undone, so whatever the last CLI test set
    is what every test after it sees. ``monkeypatch.setattr(sys, "argv",
    [...])`` is restored at teardown and is visible to a reader of the one
    test, which the suite-wide fixture is not.
    """
    offenders = []
    for path in _test_files():
        for lineno, spelling in _argv_writes(_parse(path)):
            offenders.append(f"{_rel(path)}:{lineno} {spelling} = ...")

    assert not offenders, (
        'Set argv with monkeypatch.setattr(sys, "argv", [...]), which is '
        "undone at teardown, rather than assigning it:\n  " + "\n  ".join(offenders)
    )


def test_the_guard_recognises_all_three_argv_spellings() -> None:
    """A guard that matched one spelling would read as if it matched all.

    This test cannot fail before the change -- it exercises the helper the
    change adds. It is here because the two tests above pass whether or
    not ``_argv_writes`` understands the subscript and augmented forms,
    and ``sys.argv[:] = ...`` is the form tests/conftest.py itself uses.
    """
    tree = ast.parse(
        "import sys\n"
        "sys.argv = ['a']\n"
        "sys.argv[:] = ['b']\n"
        "sys.argv += ['c']\n"
        "other.argv = ['d']\n"
        "sys.path = ['e']\n"
    )
    assert [lineno for lineno, _ in _argv_writes(tree)] == [2, 3, 4]
