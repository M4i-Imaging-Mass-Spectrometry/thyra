# tests/unit/tools/test_complexity_monitor.py
"""Tests for the CI complexity monitor's ``--files-changed`` fast mode.

Fast mode was dead for the whole of its life: it hardcoded a diff against
``origin/main``, which a shallow ``actions/checkout`` does not provide, and the
resulting failure was swallowed into an empty file list that the caller read as
"could not detect changes, analyze everything". Every PR ran a full-repository
scan while the log announced fast mode.

These tests therefore drive the real thing against real throwaway git
repositories -- a mocked ``subprocess`` would have happily agreed with the
broken version.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

_MONITOR_PATH = (
    Path(__file__).resolve().parents[3]
    / ".github"
    / "scripts"
    / "complexity_monitor.py"
)


def _load_monitor():
    """Import the monitor by path: .github/scripts is not an importable package."""
    spec = importlib.util.spec_from_file_location("complexity_monitor", _MONITOR_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


monitor = _load_monitor()


def _git(repo: Path, *args: str) -> None:
    """Run git in ``repo`` with an identity, so this works on a bare CI runner."""
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=Test",
            "-c",
            "user.email=test@example.com",
            "-c",
            "commit.gpgsign=false",
            *args,
        ],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )


def _commit(repo: Path, relative: str, body: str, message: str) -> None:
    """Write a file and commit it."""
    path = repo / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-m", message)


@pytest.fixture
def stacked_repo(tmp_path):
    """A repo with a PR branch stacked on another feature branch.

    ``main`` -> ``parent`` (adds thyra/parent.py) -> ``child`` (adds
    thyra/child.py), with an ``origin`` remote carrying all three, the way a
    checkout with ``fetch-depth: 0`` sees them.
    """
    upstream = tmp_path / "upstream"
    upstream.mkdir()
    _git(upstream, "init", "-b", "main")
    _commit(upstream, "thyra/base.py", "def base():\n    return 1\n", "base")

    _git(upstream, "checkout", "-b", "parent")
    _commit(upstream, "thyra/parent.py", "def parent():\n    return 2\n", "parent")

    _git(upstream, "checkout", "-b", "child")
    _commit(upstream, "thyra/child.py", "def child():\n    return 3\n", "child")

    clone = tmp_path / "clone"
    _git(tmp_path, "clone", str(upstream), str(clone))
    _git(clone, "checkout", "child")
    return clone


@pytest.mark.unit
def test_stacked_pr_diffs_against_its_own_base_not_main(
    stacked_repo, monkeypatch
) -> None:
    """A PR based on a feature branch must not inherit that branch's changes."""
    monkeypatch.chdir(stacked_repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "parent")

    changed = monitor.get_changed_files()

    assert changed.base_ref == "origin/parent"
    assert [p.as_posix() for p in changed.files] == ["thyra/child.py"], (
        "a stacked PR diffed against main would also report thyra/parent.py, "
        "which belongs to the parent PR"
    )


@pytest.mark.unit
def test_pr_against_main_reports_only_its_own_changes(
    stacked_repo, monkeypatch
) -> None:
    """The ordinary case still works: base main, both branch commits reported."""
    monkeypatch.chdir(stacked_repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "main")

    changed = monitor.get_changed_files()

    assert changed.base_ref == "origin/main"
    assert sorted(p.as_posix() for p in changed.files) == [
        "thyra/child.py",
        "thyra/parent.py",
    ]


@pytest.mark.unit
def test_falls_back_to_a_local_branch_when_there_is_no_remote(
    tmp_path, monkeypatch
) -> None:
    """Without an origin remote, the bare base branch name is still usable."""
    repo = tmp_path / "local"
    repo.mkdir()
    _git(repo, "init", "-b", "main")
    _commit(repo, "thyra/base.py", "def base():\n    return 1\n", "base")
    _git(repo, "checkout", "-b", "feature")
    _commit(repo, "thyra/feature.py", "def feature():\n    return 2\n", "feature")

    monkeypatch.chdir(repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "main")

    changed = monitor.get_changed_files()

    assert changed.base_ref == "main"
    assert [p.as_posix() for p in changed.files] == ["thyra/feature.py"]


@pytest.mark.unit
def test_shallow_checkout_raises_instead_of_returning_nothing(
    tmp_path, monkeypatch
) -> None:
    """The bug's exact shape: no base ref must be an error, not an empty list."""
    repo = tmp_path / "shallow"
    repo.mkdir()
    _git(repo, "init", "-b", "detached")
    _commit(repo, "thyra/only.py", "def only():\n    return 1\n", "only")

    monkeypatch.chdir(repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "main")

    with pytest.raises(monitor.ChangedFilesError) as excinfo:
        monitor.get_changed_files()

    message = str(excinfo.value)
    assert "origin/main" in message, "the message should name what it tried"
    assert "fetch-depth: 0" in message, "the message should name the fix"


@pytest.mark.unit
def test_a_missing_base_ref_never_falls_back_to_main(tmp_path, monkeypatch) -> None:
    """An unresolvable stacked base must fail, not quietly become main.

    Falling back to main here would resurrect the original defect in a subtler
    form: a plausible-looking file set that is wrong.
    """
    repo = tmp_path / "no-parent"
    repo.mkdir()
    _git(repo, "init", "-b", "main")
    _commit(repo, "thyra/base.py", "def base():\n    return 1\n", "base")
    _git(repo, "checkout", "-b", "child")
    _commit(repo, "thyra/child.py", "def child():\n    return 2\n", "child")

    monkeypatch.chdir(repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "parent-that-was-deleted")

    with pytest.raises(monitor.ChangedFilesError):
        monitor.get_changed_files()


@pytest.mark.unit
def test_no_changed_python_files_is_an_empty_list_not_a_failure(
    tmp_path, monkeypatch
) -> None:
    """A docs-only PR analyzes nothing, rather than scanning the repository."""
    repo = tmp_path / "docs-only"
    repo.mkdir()
    _git(repo, "init", "-b", "main")
    _commit(repo, "thyra/base.py", "def base():\n    return 1\n", "base")
    _git(repo, "checkout", "-b", "docs")
    _commit(repo, "README.md", "# docs\n", "docs")

    monkeypatch.chdir(repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "main")

    changed = monitor.get_changed_files()

    assert changed.files == []


@pytest.mark.unit
def test_deleted_files_are_dropped(tmp_path, monkeypatch) -> None:
    """A file removed by the PR is in the diff but cannot be analysed."""
    repo = tmp_path / "deletion"
    repo.mkdir()
    _git(repo, "init", "-b", "main")
    _commit(repo, "thyra/gone.py", "def gone():\n    return 1\n", "gone")
    _git(repo, "checkout", "-b", "remove")
    (repo / "thyra" / "gone.py").unlink()
    _git(repo, "add", "-A")
    _git(repo, "commit", "-m", "remove")

    monkeypatch.chdir(repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "main")

    assert monitor.get_changed_files().files == []


@pytest.mark.unit
def test_main_exits_two_and_analyzes_nothing_when_the_base_is_missing(
    tmp_path, monkeypatch, capsys
) -> None:
    """--files-changed must fail loudly rather than full-scan under a fast-mode label.

    This is the end-to-end version of the original bug: in CI the script printed
    a warning and went on to analyze every file in the repository.
    """
    repo = tmp_path / "shallow-main"
    repo.mkdir()
    _git(repo, "init", "-b", "detached")
    _commit(repo, "thyra/only.py", "def only():\n    return 1\n", "only")

    monkeypatch.chdir(repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "main")
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--files-changed", "--no-save"]
    )

    exit_code = monitor.main()

    assert exit_code == monitor.EXIT_CHANGED_FILES_UNKNOWN
    captured = capsys.readouterr()
    assert "ERROR" in captured.err
    assert "Analyzed" not in captured.out, "it must not have analyzed anything"
    assert "analyzing all files" not in captured.out.lower()


@pytest.mark.unit
def test_main_reports_the_base_ref_it_used(stacked_repo, monkeypatch, capsys) -> None:
    """The log must say what fast mode compared against, so it can be checked."""
    monkeypatch.chdir(stacked_repo)
    monkeypatch.setenv("GITHUB_BASE_REF", "parent")
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--files-changed", "--no-save"]
    )

    exit_code = monitor.main()

    assert exit_code == 0
    out = capsys.readouterr().out
    assert "Comparing against base ref: origin/parent" in out
    assert "Analyzing 1 changed files" in out


@pytest.mark.unit
def test_base_ref_candidates_do_not_include_main_when_a_base_is_known(
    monkeypatch,
) -> None:
    """The candidate list itself encodes the rule, so pin it down."""
    monkeypatch.setenv("GITHUB_BASE_REF", "some/feature-branch")
    assert monitor._base_ref_candidates() == [
        "origin/some/feature-branch",
        "some/feature-branch",
    ]

    monkeypatch.delenv("GITHUB_BASE_REF")
    assert monitor._base_ref_candidates() == ["origin/main", "main", "HEAD~1"]
