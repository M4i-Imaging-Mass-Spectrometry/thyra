# tests/unit/test_py_typed_marker.py
"""The PEP 561 marker and the ``Typing :: Typed`` classifier must agree.

``pyproject.toml`` declared ``Typing :: Typed`` while ``thyra/py.typed`` did
not exist, so the published wheel asserted typed-ness and shipped no marker.
PEP 561 requires a checker to ignore inline annotations in an installed package
with no marker, so 94% annotation coverage bought a downstream consumer nothing
-- every thyra symbol resolved to ``Any``, and the classifier PyPI shows was
false.

Two honest fixes existed: ship the marker, or drop the classifier.  These tests
pin the one this repository took and keep the two from drifting apart again in
either direction.

None of them reads a BUILT artefact, and none can.  ``uv.lock`` installs thyra
editable for the ``test`` job, and the ``clean-venv-install`` job runs pytest
with ``PYTHONPATH=.``, so ``importlib.resources.files("thyra")`` resolves to
the repository's own ``thyra/`` in every pytest lane -- an assertion here would
pass on the repo file even if the wheel shipped nothing.  The check against the
INSTALLED distribution therefore lives in ``.github/workflows/tests.yml``, and
the last test below guards that it stays there.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]
_MARKER = _REPO_ROOT / "thyra" / "py.typed"
_PYPROJECT = _REPO_ROOT / "pyproject.toml"
_TESTS_WORKFLOW = _REPO_ROOT / ".github" / "workflows" / "tests.yml"

_CLASSIFIER = "Typing :: Typed"


def test_classifier_and_marker_agree() -> None:
    """The declared classifier and the shipped marker must say the same thing.

    REGRESSION GUARD for the defect itself: against the tree before this
    change the classifier was declared and the file did not exist, so this
    fails.  It is the only assertion that fails in BOTH directions of drift --
    adding the classifier without the marker, or deleting the marker while the
    classifier stands -- which is what makes the issue's "two honest fixes"
    mechanical rather than prose.
    """
    with _PYPROJECT.open("rb") as handle:
        pyproject = tomllib.load(handle)

    declares_typed = _CLASSIFIER in pyproject["project"]["classifiers"]
    ships_marker = _MARKER.is_file()

    assert declares_typed == ships_marker, (
        f"pyproject.toml {'declares' if declares_typed else 'omits'} the "
        f"{_CLASSIFIER!r} classifier but thyra/py.typed "
        f"{'exists' if ships_marker else 'does not exist'}. A type checker "
        "reads the marker and not the classifier, so the two must agree: "
        "either ship thyra/py.typed or drop the classifier."
    )


def test_the_marker_is_the_half_that_was_kept() -> None:
    """Of the two honest fixes, this repository ships the marker.

    A REDUNDANT REGRESSION GUARD, disclosed as such: it does fail against the
    tree before this change, but only for the same reason the test above does,
    so it catches nothing there that the agreement test misses.  Its own job is
    the direction that test permits -- the cheapest way to silence a red
    agreement test is to delete the classifier, and this makes that a
    deliberate edit of two tests rather than a one-line pyproject.toml change.
    The annotations are project policy (docs/contributing.md, "Use type hints
    for all public functions"), so the agreement is meant to be reached by
    keeping both, not neither.

    Deliberately does not assert the file is empty: zero bytes is how it was
    created, but PEP 561 also allows a ``partial`` marker, and forbidding that
    here would be an unrelated rule.
    """
    assert _MARKER.is_file(), (
        "thyra/py.typed is missing. The Typing :: Typed classifier in "
        "pyproject.toml is a promise a type checker reads through this file; "
        "without it PEP 561 requires every thyra annotation to be ignored."
    )


def test_the_installed_distribution_is_checked_in_ci() -> None:
    """One CI step must read the marker out of the INSTALLED distribution.

    REGRESSION GUARD: before this change no step in tests.yml mentioned
    py.typed at all.  The marker is a zero-byte file, so no import, test or
    lint failure can notice it going missing from a build -- a change to
    ``[tool.hatch.build.targets]`` that stopped sweeping it would be invisible
    everywhere else, and the pytest lanes structurally cannot see it (see this
    module's docstring).  The placement is the substance of the check, so it is
    asserted rather than left to the step's comment: it must run from outside
    the checkout, or it reads the repo's ``./thyra`` instead, and before the
    test-only installs, so it exercises exactly the closure pyproject.toml
    declares.
    """
    workflow = yaml.safe_load(_TESTS_WORKFLOW.read_text(encoding="utf-8"))
    steps = workflow["jobs"]["clean-venv-install"]["steps"]

    checking = [i for i, step in enumerate(steps) if "py.typed" in step.get("run", "")]
    assert len(checking) == 1, (
        "expected exactly one step of the clean-venv-install job to assert "
        f"the installed distribution ships thyra/py.typed, found {len(checking)}"
    )
    (index,) = checking

    assert "cd /tmp" in steps[index]["run"], (
        "the marker check must leave the checkout before importing thyra, or "
        "it reads the repository's ./thyra rather than the installed "
        "distribution and can never fail"
    )

    test_deps = [
        i for i, step in enumerate(steps) if "pip install pytest" in step.get("run", "")
    ]
    assert test_deps, "the clean-venv-install job no longer installs test deps"
    assert index < min(test_deps), (
        "the marker check must run BEFORE the test-only installs, so it sees "
        "exactly the dependency closure pyproject.toml declares"
    )
