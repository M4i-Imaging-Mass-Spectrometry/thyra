# tests/unit/tools/test_ci_workflows.py
"""Tests that the CI lint gate is actually wired up, and stays wired up.

Six lint, type and security tools were configured across ``.flake8``,
``.pre-commit-config.yaml`` and ``pyproject.toml`` and ran nowhere but a local
hook a contributor may never install. The fix is a ``lint`` job in
``.github/workflows/tests.yml``; the job is its own proof that the tree is
clean, since it either goes green or it does not.

What a job cannot prove is that it will still be there. The integration lane is
the precedent -- ``tests.yml`` records that it "sat with 15 of its 18 tests
failing precisely because no job anywhere ran it" -- so these tests guard the
WIRING rather than the tools' findings: that a job still runs every hook, that
the triggers it inherits are not narrowed, that the hook set is not quietly
shrunk, and that the pins which make a hook reproducible stay pinned.

Two of these assertions (``test_trigger_block_covers_prs_and_main`` and
``test_the_six_tools_and_the_path_guards_are_all_hooks``) pass against the tree
as it was before the lint job existed. They are not regression guards for that
change and cannot fail without it; they guard a later edit instead, and say so
in their own docstrings.
"""

from __future__ import annotations

import re
import shlex
import tomllib
from pathlib import Path
from typing import Any, Dict, List

import yaml

_REPO_ROOT = Path(__file__).resolve().parents[3]
_TESTS_WORKFLOW = _REPO_ROOT / ".github" / "workflows" / "tests.yml"
_PRE_COMMIT_CONFIG = _REPO_ROOT / ".pre-commit-config.yaml"
_UV_LOCK = _REPO_ROOT / "uv.lock"

# PyYAML parses YAML 1.1, where the bare key `on:` is the boolean True rather
# than the string "on". A workflow's trigger block therefore lives under the
# key True, and `doc["on"]` raises KeyError.
_ON = True


def _load_yaml(path: Path) -> Dict[str, Any]:
    """Parse a YAML file from the repository root."""
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _steps_running_pre_commit() -> List[Dict[str, Any]]:
    """Every step in tests.yml whose `run:` invokes `pre-commit run --all-files`."""
    workflow = _load_yaml(_TESTS_WORKFLOW)
    found = []
    for job in workflow["jobs"].values():
        for step in job.get("steps", []):
            run = step.get("run", "")
            if "pre-commit run --all-files" in run:
                found.append(step)
    return found


def _hooks_with_ids() -> List[Dict[str, Any]]:
    """Flatten .pre-commit-config.yaml to the list of hook mappings it declares."""
    config = _load_yaml(_PRE_COMMIT_CONFIG)
    return [hook for repo in config["repos"] for hook in repo["hooks"]]


def _repo_package_name(repo_url: str) -> str:
    """The distribution a hook repo is expected to provide, from its URL."""
    return repo_url.rstrip("/").rsplit("/", 1)[-1].lower()


def _bare_requirement_name(dependency: str) -> str:
    """`bandit[toml]==1.7.9` -> `bandit`: drop extras and any version specifier."""
    name = re.split(r"[<>=!~\[;]", dependency, maxsplit=1)[0]
    return name.strip().lower()


def test_exactly_one_job_runs_every_pre_commit_hook() -> None:
    """A job must run `pre-commit run --all-files`.

    REGRESSION GUARD for the gap this file exists for: before the lint job was
    added this found nothing, because nothing in .github/workflows mentioned
    pre-commit at all. Deleting or renaming the job away fails here.
    """
    steps = _steps_running_pre_commit()
    assert len(steps) == 1, (
        "expected exactly one step in tests.yml to run "
        f"`pre-commit run --all-files`, found {len(steps)}"
    )


def test_the_lint_step_narrows_the_hook_set_in_no_way() -> None:
    """The gate must run the whole config, not a named subset of it.

    REGRESSION GUARD. A positional hook id (`pre-commit run --all-files black`)
    or a `--hook-stage` would still satisfy the test above while quietly
    enforcing a fraction of the file, and would mean a seventh tool added to
    .pre-commit-config.yaml needs a workflow edit before it is enforced.
    """
    (step,) = _steps_running_pre_commit()
    tokens = shlex.split(step["run"])
    after_run = tokens[tokens.index("run", tokens.index("pre-commit")) + 1 :]

    assert "--hook-stage" not in after_run

    declared_ids = {hook["id"] for hook in _hooks_with_ids()}
    named = [token for token in after_run if token in declared_ids]
    assert not named, f"the lint step narrows the gate to {named}"


def test_trigger_block_covers_prs_and_main() -> None:
    """The lint job inherits tests.yml's triggers; they must stay as they are.

    NOT A REGRESSION GUARD -- this passes against the tree before the lint job
    existed, because these triggers already served the test job. It guards a
    later narrowing: tests.yml's own comment says the `push: [main]` trigger is
    load-bearing for releases, and `pull_request` is deliberately unfiltered so
    a stacked PR is not skipped. Both of those are prose today and enforced
    here.
    """
    workflow = _load_yaml(_TESTS_WORKFLOW)
    triggers = workflow[_ON]

    assert "pull_request" in triggers
    assert triggers["pull_request"] is None, (
        "a `branches:` filter on pull_request matches the BASE branch, which "
        "skips every stacked PR until its parent merges"
    )
    assert triggers["push"]["branches"] == ["main"]


def test_every_additional_dependency_is_version_pinned() -> None:
    """Hook `additional_dependencies` must pin with `==`, bar the self-extra.

    REGRESSION GUARD. An unpinned extra resolves whatever is newest the day a
    hook env happens to be built, so the same commit passes on one machine and
    fails on another -- and once this gate blocks the repository, an upstream
    release of a typeshed stub turns an unrelated PR red.

    The exemption is for an entry naming the hook repo's OWN package, such as
    bandit's `bandit[toml]`. `rev` already pins that, and pinning it as well
    does not install: pre-commit builds the env with `pip install .` from its
    clone, whose metadata reports 0.0.0, and pip refuses the contradiction.
    """
    config = _load_yaml(_PRE_COMMIT_CONFIG)
    unpinned = [
        (hook["id"], dependency)
        for repo in config["repos"]
        for hook in repo["hooks"]
        for dependency in hook.get("additional_dependencies", [])
        if "==" not in dependency
        and _bare_requirement_name(dependency)
        != _repo_package_name(repo.get("repo", ""))
    ]
    assert not unpinned, f"unpinned hook dependencies: {unpinned}"


def test_the_six_tools_and_the_path_guards_are_all_hooks() -> None:
    """The tools the gate exists to run must still be in the config.

    NOT A REGRESSION GUARD -- every one of these ids was already present before
    the lint job existed. It guards the other direction: with CI running the
    hooks, deleting one is now the cheapest way to make a red check green, and
    this makes that a deliberate test edit rather than a one-line config edit.
    """
    declared_ids = {hook["id"] for hook in _hooks_with_ids()}
    expected = {
        "black",
        "isort",
        "flake8",
        "mypy",
        "bandit",
        "pydocstyle",
        "no-lab-share-paths",
        "no-user-home-paths",
    }
    assert expected <= declared_ids, f"missing hooks: {sorted(expected - declared_ids)}"


def test_the_black_hook_runs_the_black_the_lockfile_resolves() -> None:
    """The hook's black and `uv run black .`'s black must be one version.

    REGRESSION GUARD, and the only assertion that makes "CI and a developer
    machine cannot diverge" enforceable rather than aspirational. They had
    diverged: `rev: 24.4.2` against a lockfile resolving 26.5.1, so the
    `uv run black .` docs/contributing.md documents reformatted five files the
    hook considered formatted. Under a CI gate that is a contributor doing
    exactly what CONTRIBUTING says and getting a red PR.
    """
    black_hook_revs = [
        repo["rev"]
        for repo in _load_yaml(_PRE_COMMIT_CONFIG)["repos"]
        if any(hook["id"] == "black" for hook in repo["hooks"])
    ]
    assert len(black_hook_revs) == 1

    with _UV_LOCK.open("rb") as handle:
        lock = tomllib.load(handle)
    locked = [
        package["version"] for package in lock["package"] if package["name"] == "black"
    ]
    assert len(locked) == 1, "uv.lock should resolve exactly one black"

    assert black_hook_revs[0] == locked[0], (
        f"the black hook runs {black_hook_revs[0]} but uv.lock resolves "
        f"{locked[0]}; `uv run black .` and the hook would disagree"
    )
