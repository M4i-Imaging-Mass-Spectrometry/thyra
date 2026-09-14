# tests/unit/tools/test_complexity_monitor.py
"""Tests for the CI complexity monitor's ``--files-changed`` mode and output.

Fast mode was dead for the whole of its life: it hardcoded a diff against
``origin/main``, which a shallow ``actions/checkout`` does not provide, and the
resulting failure was swallowed into an empty file list that the caller read as
"could not detect changes, analyze everything". Every PR ran a full-repository
scan while the log announced fast mode.

These tests therefore drive the real thing against real throwaway git
repositories -- a mocked ``subprocess`` would have happily agreed with the
broken version.

The second group covers what the script prints. ``--quiet`` documented itself
as "Suppress output except violations" while suppressing the violations too,
which made it a flag with no output at any verbosity.
"""

from __future__ import annotations

import importlib.util
import json
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


@pytest.fixture
def tree_with_one_tangled_function(tmp_path):
    """A working directory whose ``thyra`` package holds one function of complexity 4.

    ``main`` analyzes ``./thyra`` when it exists, so the package name matters;
    the threshold each test passes decides whether that function counts as a
    violation.
    """
    package = tmp_path / "thyra"
    package.mkdir()
    (package / "tangled.py").write_text(
        "def tangled(a, b, c):\n"
        "    if a:\n"
        "        return 1\n"
        "    if b:\n"
        "        return 2\n"
        "    if c:\n"
        "        return 3\n"
        "    return 4\n",
        encoding="utf-8",
    )
    return tmp_path


def test_quiet_still_lists_the_violations(
    tree_with_one_tangled_function, monkeypatch, capsys
) -> None:
    """The whole point of the flag: everything suppressed *except* violations.

    It used to suppress those as well, so ``--quiet`` was silent whatever the
    state of the code -- and the workflow step that used it printed a heading
    over nothing on every run.
    """
    monkeypatch.chdir(tree_with_one_tangled_function)
    monkeypatch.setattr(
        sys,
        "argv",
        ["complexity_monitor.py", "--threshold", "3", "--no-save", "--quiet"],
    )

    exit_code = monitor.main()

    assert exit_code == 1
    out = capsys.readouterr().out
    assert "tangled" in out, "the violation itself must survive --quiet"
    assert "Analyzed" not in out, "the summary is what --quiet is for"
    assert "Complexity threshold" not in out


def test_quiet_prints_nothing_when_there_are_no_violations(
    tree_with_one_tangled_function, monkeypatch, capsys
) -> None:
    """The other half of the contract: nothing to report means no output."""
    monkeypatch.chdir(tree_with_one_tangled_function)
    monkeypatch.setattr(
        sys,
        "argv",
        ["complexity_monitor.py", "--threshold", "15", "--no-save", "--quiet"],
    )

    exit_code = monitor.main()

    assert exit_code == 0
    assert capsys.readouterr().out == ""


def test_the_summary_and_the_listing_both_print_without_quiet(
    tree_with_one_tangled_function, monkeypatch, capsys
) -> None:
    """Moving the listing out of the guard must not thin out the default output."""
    monkeypatch.chdir(tree_with_one_tangled_function)
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--threshold", "3", "--no-save"]
    )

    exit_code = monitor.main()

    assert exit_code == 1
    out = capsys.readouterr().out
    assert "Analyzed 1 files, 1 functions" in out
    assert "Complexity threshold: 3" in out
    assert "Violations found: 1" in out
    assert "Maximum complexity: 4" in out
    assert "Top violations:" in out
    assert "tangled (4)" in out


def test_the_workflow_debug_invocation_says_something_on_a_clean_tree(
    tree_with_one_tangled_function, monkeypatch, capsys
) -> None:
    """The flags of complexity-monitoring.yml's "DIRECT COMPLEXITY CHECK" step.

    That step exists to show a human the whole-repo numbers, so on a repository
    with no violations -- the state this one is in -- it still has to report
    them. This is why the step does not pass ``--quiet``, even now that
    ``--quiet`` behaves.
    """
    monkeypatch.chdir(tree_with_one_tangled_function)
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--threshold", "15", "--no-save"]
    )

    exit_code = monitor.main()

    assert exit_code == 0
    out = capsys.readouterr().out
    assert "Violations found: 0" in out
    assert "Analyzed 1 files" in out


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


@pytest.fixture
def tree_with_a_bom_prefixed_function(tmp_path):
    """The same complexity-4 function, written with a UTF-8 BOM in front of it.

    Written with ``write_bytes`` rather than ``write_text(encoding="utf-8-sig")``
    so that the three bytes under test are visible in the test itself.
    """
    package = tmp_path / "thyra"
    package.mkdir()
    (package / "tangled.py").write_bytes(
        b"\xef\xbb\xbf"
        b"def tangled(a, b, c):\n"
        b"    if a:\n"
        b"        return 1\n"
        b"    if b:\n"
        b"        return 2\n"
        b"    if c:\n"
        b"        return 3\n"
        b"    return 4\n"
    )
    return tmp_path


@pytest.fixture
def tree_with_one_good_and_one_broken_file(tmp_path):
    """A ``thyra`` package holding one analysable file and one that cannot parse.

    The mix is the point: with only a broken file the run would be empty either
    way, and the defect under test is a broken file disappearing into an
    otherwise healthy, passing run.
    """
    package = tmp_path / "thyra"
    package.mkdir()
    (package / "fine.py").write_text(
        "def fine(a, b, c):\n"
        "    if a:\n"
        "        return 1\n"
        "    if b:\n"
        "        return 2\n"
        "    if c:\n"
        "        return 3\n"
        "    return 4\n",
        encoding="utf-8",
    )
    (package / "broken.py").write_text("def f(:\n", encoding="utf-8")
    return tmp_path


def test_a_byte_order_mark_does_not_exempt_a_file(
    tree_with_a_bom_prefixed_function, monkeypatch, capsys
) -> None:
    """A BOM must not remove a file from the gate.

    Regression guard. ``ast.parse`` rejects U+FEFF as a non-printable
    character, so reading with plain ``utf-8`` turned the whole file into a
    warning and an empty result list: before the fix this run reported
    "Analyzed 1 files, 0 functions" and exited 0, the violation never found.
    """
    monkeypatch.chdir(tree_with_a_bom_prefixed_function)
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--threshold", "3", "--no-save"]
    )

    exit_code = monitor.main()

    assert exit_code == 1, "the BOM'd file's violation must be found"
    out = capsys.readouterr().out
    assert "Analyzed 1 files, 1 functions" in out
    assert "tangled (4)" in out


def test_an_unparseable_file_raises_instead_of_returning_no_functions(
    tmp_path,
) -> None:
    """The shape of the defect: "could not read it" must not be an empty result.

    Regression guard, and the counterpart of
    ``test_shallow_checkout_raises_instead_of_returning_nothing`` above.
    ``generate_report`` cannot tell an empty list from a file holding no
    functions, so returning one exempted the file silently.
    """
    broken = tmp_path / "broken.py"
    broken.write_text("def f(:\n", encoding="utf-8")

    with pytest.raises(monitor.FileAnalysisError) as excinfo:
        monitor.analyze_file(broken)

    assert "broken.py" in str(excinfo.value), "the message should name the file"
    assert excinfo.value.__cause__ is not None, "the parse error should be chained"


def test_a_file_that_cannot_be_opened_is_reported_the_same_way(tmp_path) -> None:
    """An unreadable path is a gate that did not run, like an unparseable one.

    This case never even had a warning to fall back on: ``OSError`` sat outside
    the old ``except`` clause entirely, so a missing file escaped
    ``analyze_file`` as a bare traceback.
    """
    with pytest.raises(monitor.FileAnalysisError):
        monitor.analyze_file(tmp_path / "does_not_exist.py")


def test_a_deeply_nested_lambda_chain_does_not_escape_as_a_traceback(
    tmp_path,
) -> None:
    """The visitor's own recursion is a parse failure too, and must be caught.

    A ``lambda:`` chain needs neither parentheses nor indentation, so it clears
    the tokenizer's depth limits and then exhausts the stack inside
    ``ComplexityAnalyzer.visit``. Before the fix that ``RecursionError``
    escaped and exited 1, which the workflow prints as "complexity violations
    found" -- the same conflation of a broken gate with a code-quality result
    that this issue exists to remove.
    """
    deep = tmp_path / "deep.py"
    deep.write_text(
        "def outer():\n    f = " + "lambda: " * 500 + "1\n", encoding="utf-8"
    )

    with pytest.raises(monitor.FileAnalysisError):
        monitor.analyze_file(deep)


def test_main_exits_three_and_names_the_file_it_could_not_parse(
    tree_with_one_good_and_one_broken_file, monkeypatch, capsys
) -> None:
    """A file the gate could not read fails the run, under its own exit code.

    Regression guard. Before the fix this exact tree exited 0 with a warning on
    stdout, which the workflow reports as "SUCCESS: No complexity violations
    found".
    """
    monkeypatch.chdir(tree_with_one_good_and_one_broken_file)
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--threshold", "15", "--no-save"]
    )

    exit_code = monitor.main()

    assert exit_code == monitor.EXIT_UNREADABLE_FILES
    err = capsys.readouterr().err
    assert "ERROR" in err
    assert "broken.py" in err, "the log must name the file that was skipped"
    assert "did not cover" in err


def test_the_analysed_count_excludes_a_file_that_could_not_be_parsed(
    tree_with_one_good_and_one_broken_file, monkeypatch, capsys
) -> None:
    """The count must be files read, not files selected.

    Regression guard for the half of the defect the exit code does not cover:
    the skipped file was not merely absent from the analysis, it was
    affirmatively counted as analysed.
    """
    monkeypatch.chdir(tree_with_one_good_and_one_broken_file)
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--threshold", "15", "--no-save"]
    )

    monitor.main()

    out = capsys.readouterr().out
    assert "Analyzed 1 files" in out
    assert "Analyzed 2 files" not in out, "an unparseable file is not an analysed one"


def test_unreadable_files_survive_quiet(
    tree_with_one_good_and_one_broken_file, monkeypatch, capsys
) -> None:
    """--quiet suppresses findings about the code, not news that the gate broke."""
    monkeypatch.chdir(tree_with_one_good_and_one_broken_file)
    monkeypatch.setattr(
        sys,
        "argv",
        ["complexity_monitor.py", "--threshold", "15", "--no-save", "--quiet"],
    )

    exit_code = monitor.main()

    assert exit_code == monitor.EXIT_UNREADABLE_FILES
    captured = capsys.readouterr()
    assert "broken.py" in captured.err
    assert captured.out == "", "the summary is still what --quiet is for"


def test_the_report_is_still_written_when_a_file_could_not_be_parsed(
    tree_with_one_good_and_one_broken_file, monkeypatch
) -> None:
    """Exiting 3 must not cost the artifact the PR comment is built from.

    Regression guard against this fix's own likeliest mistake: returning as
    soon as the failure is known skips the save block, which leaves the
    workflow uploading an empty artifact and posting no comment at all --
    strictly less feedback than the silent exemption being replaced.
    """
    monkeypatch.chdir(tree_with_one_good_and_one_broken_file)
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--threshold", "15", "--save"]
    )

    exit_code = monitor.main()

    assert exit_code == monitor.EXIT_UNREADABLE_FILES
    reports = sorted(
        (tree_with_one_good_and_one_broken_file / "reports" / "complexity").glob(
            "complexity_report_*.json"
        )
    )
    assert reports, "the report must still be saved on a failing run"

    saved = json.loads(reports[-1].read_text(encoding="utf-8"))
    assert saved["files_analyzed"] == 1
    assert len(saved["unreadable_files"]) == 1
    assert "broken.py" in saved["unreadable_files"][0]


def test_the_report_records_what_was_analysed_and_what_was_not() -> None:
    """Both keys are part of the report's contract: the PR comment reads them.

    Cannot fail before the change -- the keys did not exist, so this pins the
    new contract rather than guarding the old defect.
    """
    report = monitor.generate_report(
        [], 15, files_analyzed=1, unreadable_files=["thyra/broken.py: invalid syntax"]
    )

    assert report["files_analyzed"] == 1
    assert report["unreadable_files"] == ["thyra/broken.py: invalid syntax"]


def test_the_unreadable_exit_code_is_its_own_number() -> None:
    """complexity-monitoring.yml branches on the literal 3, so pin it here.

    Cannot fail before the change -- the constant is new. It exists so a later
    edit cannot renumber it without the workflow's ``elif`` being noticed.
    """
    assert monitor.EXIT_UNREADABLE_FILES == 3
    assert monitor.EXIT_UNREADABLE_FILES != monitor.EXIT_CHANGED_FILES_UNKNOWN


def test_the_headroom_line_names_the_function_holding_the_margin(
    tree_with_one_tangled_function, monkeypatch, capsys
) -> None:
    """The gate reports how much room is left, not only whether it was exceeded.

    Cannot fail before the change -- the line is new. The threshold in this
    repository equals the worst function in it, so the margin is zero and the
    gate is a tripwire under whoever next edits that function; printing the
    margin is what keeps CI from being the first to mention it.
    """
    monkeypatch.chdir(tree_with_one_tangled_function)
    monkeypatch.setattr(
        sys, "argv", ["complexity_monitor.py", "--threshold", "5", "--no-save"]
    )

    assert monitor.main() == 0
    out = capsys.readouterr().out
    assert "Closest to the threshold: tangled" in out
    assert "at 4, 1 below the limit of 5" in out


def test_no_tracked_source_file_carries_a_byte_order_mark() -> None:
    """No BOM may enter the tree, whatever .pre-commit-config.yaml is doing.

    Driven from the index rather than the filesystem for two reasons: an
    untracked local scratch file must not fail the suite, and reading what is
    committed is what makes this a guard against ``git commit --no-verify``
    too. No workflow runs pre-commit today (see #291), so this assertion is the
    only part of the BOM rule that runs in CI at all.
    """
    repo_root = Path(__file__).resolve().parents[3]
    listed = subprocess.run(
        ["git", "ls-files", "-z", "--", "thyra/*.py", "tests/*.py"],
        cwd=repo_root,
        check=True,
        capture_output=True,
    )

    tracked = [name for name in listed.stdout.split(b"\0") if name]
    assert tracked, "the query should match this repository's Python sources"

    offenders = [
        name.decode()
        for name in tracked
        if (repo_root / name.decode()).read_bytes().startswith(b"\xef\xbb\xbf")
    ]
    assert offenders == [], "a UTF-8 BOM makes a file unparseable to ast.parse"
