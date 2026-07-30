# tests/__init__.py
#
# This file is load-bearing, not decorative. Without it `tests` is only a
# namespace package, and a *regular* package always wins over namespace
# portions no matter what sys.path order says. Several third-party
# distributions ship their own test suite as a top-level `tests` package
# (msalign 0.2.0 installs tests/__init__.py; choreographer and kaleido dump
# modules into the same directory), so on any machine where one of those is
# installed `from tests.fixtures.mock_msi_generator import ...` resolves to
# site-packages instead of this repo and collection fails. Setting PYTHONPATH
# does not help. Keeping `tests` a regular package -- with __init__.py in every
# subdirectory below it -- makes the import resolve here regardless of what
# else the ambient environment has installed.
