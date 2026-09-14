# tests/unit/test_lane_markers.py
"""Tests that the two lane markers partition the suite by directory.

``pytest -m integration`` is the command both ``tests/README.md`` and the
``integration`` job in ``.github/workflows/tests.yml`` present as "the
integration lane". It used to mean something narrower: the files somebody had
remembered to decorate. 14 of ``tests/integration/``'s 65 tests carried no
marker and were deselected there, including every end-to-end CLI test, while 14
tests under ``tests/unit/`` claimed the lane by decorator. The lane was a false
green in one direction and a misfiling in the other.

``tests/conftest.py`` now stamps the marker from the directory, so the two
agree by construction and the only remaining way to get it wrong is to put a
file in the wrong place. That is what these tests watch for -- and they watch
it through pytest itself rather than by re-implementing the hook, because a
re-implementation would agree with a broken hook.

The assertions read exit codes, which is exact here: ``--collect-only`` exits 5
when everything collected was deselected and 0 when anything survived. A
genuine collection error exits 2, which is why every assertion carries the
child's output -- otherwise a plain ``ImportError`` would be reported as a lane
leak.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import List

_REPO_ROOT = Path(__file__).resolve().parents[2]

# pytest-cov installs a .pth that starts coverage in every child process, which
# would drop .coverage.* files in the repository root for these subprocesses to
# have their collection-time imports combined into the reported figure. This
# file exists to keep the lane honest; it must not quietly inflate coverage
# while doing it.
_COVERAGE_SUBPROCESS_VARS = (
    "COV_CORE_SOURCE",
    "COV_CORE_CONFIG",
    "COV_CORE_DATAFILE",
    "COVERAGE_PROCESS_START",
)

# 5 is pytest's "no tests ran". After --collect-only with a -m expression that
# deselects everything, that is the success we are asserting.
_EXIT_NO_TESTS_COLLECTED = 5


def _collect(*args: str) -> subprocess.CompletedProcess:
    """Run `pytest --collect-only` in a child process from the repository root.

    Args:
        *args: Extra arguments, typically a path and a `-m` expression.

    Returns:
        The completed process, with stdout and stderr captured.
    """
    env = {
        key: value
        for key, value in os.environ.items()
        if key not in _COVERAGE_SUBPROCESS_VARS
    }
    command: List[str] = [
        sys.executable,
        "-m",
        "pytest",
        *args,
        "--collect-only",
        "-q",
        "-p",
        "no:cacheprovider",
    ]
    return subprocess.run(
        command,
        cwd=_REPO_ROOT,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


def _assert_nothing_collected(result: subprocess.CompletedProcess, why: str) -> None:
    """Assert the child deselected everything, showing its output if not."""
    assert result.returncode == _EXIT_NO_TESTS_COLLECTED, (
        f"{why}\nexit {result.returncode} (expected "
        f"{_EXIT_NO_TESTS_COLLECTED}, and 2 means the collection itself "
        f"failed)\n--- stdout ---\n{result.stdout}\n--- stderr ---\n"
        f"{result.stderr}"
    )


def test_nothing_under_tests_integration_escapes_the_integration_lane() -> None:
    """Every test in tests/integration/ must be selected by `-m integration`.

    REGRESSION GUARD for the defect itself. Without the conftest hook this
    collects the 14 undecorated tests -- the six CLI end-to-end tests, three
    imzML conversion tests and the five real-acquisition Bruker tests -- and
    exits 0, which is the false green: `pytest -m integration` would report
    success having never run them.
    """
    result = _collect("tests/integration", "-m", "not integration")
    _assert_nothing_collected(
        result,
        "a test under tests/integration/ is not in the integration lane, so "
        "`pytest -m integration` silently skips it",
    )


def test_no_test_under_tests_unit_claims_the_integration_lane() -> None:
    """No test in tests/unit/ may carry the integration marker.

    REGRESSION GUARD, for the same defect mirrored. Before this change 14
    tests in tests/unit/converters/test_streaming_converter.py were decorated
    `@pytest.mark.integration` and ran in the integration job from outside the
    directory that job claims to be.

    Note this is NOT expressible as `-m "not unit"`: the hook stamps `unit` on
    everything under the directory whatever else it carries, so that
    expression deselects the whole tree and passes even with a stray
    `integration` decorator sitting in it.
    """
    result = _collect("tests/unit", "-m", "integration")
    _assert_nothing_collected(
        result,
        "a test under tests/unit/ claims the integration lane; move the file "
        "to tests/integration/ rather than marking it",
    )


def test_no_test_under_tests_integration_claims_the_unit_lane() -> None:
    """No test in tests/integration/ may carry the unit marker.

    NOT A REGRESSION GUARD -- this passes against the tree before the hook
    existed, because no such decorator has ever been written. It is meaningful
    only because the hook's `elif` never stamps `unit` on an item under
    tests/integration/, so the marker can appear there by hand and no other
    way; this makes that a test failure rather than a quiet double-lane
    membership.
    """
    result = _collect("tests/integration", "-m", "unit")
    _assert_nothing_collected(
        result,
        "a test under tests/integration/ carries the unit marker by hand",
    )


def test_the_two_lanes_cover_the_whole_suite() -> None:
    """Every collected test must land in one lane or the other.

    REGRESSION GUARD, and the broadest of the four: without the hook almost
    nothing carries a marker at all, so this collects 2,583 of 2,690 tests and
    exits 0. Afterwards it also guards the one arrangement the hook cannot
    stamp -- a test file placed directly under tests/, or under a third
    directory beside the two. Such a file would run in the default lane and
    never in `-m integration`, which is the invisible membership this change
    exists to remove.
    """
    result = _collect("-m", "not (unit or integration)")
    _assert_nothing_collected(
        result,
        "a test sits outside both tests/unit/ and tests/integration/, so it "
        "carries no lane marker and belongs to no lane",
    )
