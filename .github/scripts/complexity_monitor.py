#!/usr/bin/env python3
"""Complexity monitoring script for CI/CD pipeline."""

import argparse
import ast
import json
import os

# nosec B404: CI tooling that shells out to git with a fixed argument list and
# no shell. B603 and B607, which cover the calls themselves, are already
# skipped repository-wide in pyproject.toml.
import subprocess  # nosec B404
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, NamedTuple

# Generous enough that a cold pack on a CI runner cannot make the timeout the
# thing that fails the build: since --files-changed now refuses to guess, a
# timeout here is a hard error rather than a quiet fallback.
GIT_TIMEOUT_SECONDS = 30

# Distinct from 1 ("complexity violations found") so that a broken --files-changed
# mode cannot be mistaken for a code-quality failure.
EXIT_CHANGED_FILES_UNKNOWN = 2

# The console listing is a pointer, not the record: every violation is in the
# saved JSON report, which the workflow uploads as an artifact.
TOP_VIOLATIONS_SHOWN = 5


class ComplexityResult(NamedTuple):
    """Result of complexity analysis for a function."""

    file: str
    line: int
    function: str
    complexity: int


class ComplexityAnalyzer(ast.NodeVisitor):
    """AST visitor to calculate cyclomatic complexity."""

    def __init__(self):
        """Initialize the complexity analyzer."""
        self.complexity = 1  # Base complexity
        self.function_complexities: List[ComplexityResult] = []
        self.current_file = ""

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        """Visit function definition and calculate complexity."""
        old_complexity = self.complexity
        self.complexity = 1  # Reset for this function

        # Visit all child nodes
        self.generic_visit(node)

        # Store result
        result = ComplexityResult(
            file=self.current_file,
            line=node.lineno,
            function=node.name,
            complexity=self.complexity,
        )
        self.function_complexities.append(result)

        # Restore previous complexity
        self.complexity = old_complexity

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        """Visit async function definition."""
        # Handle async functions the same way as regular functions
        old_complexity = self.complexity
        self.complexity = 1  # Reset for this function

        # Visit all child nodes
        self.generic_visit(node)

        # Store result
        result = ComplexityResult(
            file=self.current_file,
            line=node.lineno,
            function=node.name,
            complexity=self.complexity,
        )
        self.function_complexities.append(result)

        # Restore previous complexity
        self.complexity = old_complexity

    def visit_If(self, node: ast.If) -> None:
        """Visit if statement."""
        self.complexity += 1
        self.generic_visit(node)

    def visit_For(self, node: ast.For) -> None:
        """Visit for loop."""
        self.complexity += 1
        self.generic_visit(node)

    def visit_AsyncFor(self, node: ast.AsyncFor) -> None:
        """Visit async for loop."""
        self.complexity += 1
        self.generic_visit(node)

    def visit_While(self, node: ast.While) -> None:
        """Visit while loop."""
        self.complexity += 1
        self.generic_visit(node)

    def visit_Try(self, node: ast.Try) -> None:
        """Visit try statement."""
        # Each except handler adds complexity
        self.complexity += len(node.handlers)
        if node.orelse:
            self.complexity += 1
        if node.finalbody:
            self.complexity += 1
        self.generic_visit(node)

    def visit_With(self, node: ast.With) -> None:
        """Visit with statement."""
        self.complexity += 1
        self.generic_visit(node)

    def visit_AsyncWith(self, node: ast.AsyncWith) -> None:
        """Visit async with statement."""
        self.complexity += 1
        self.generic_visit(node)

    def visit_BoolOp(self, node: ast.BoolOp) -> None:
        """Visit boolean operation (and/or)."""
        # Each additional condition adds complexity
        self.complexity += len(node.values) - 1
        self.generic_visit(node)


def analyze_file(file_path: Path) -> List[ComplexityResult]:
    """Analyze a Python file for cyclomatic complexity."""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()

        tree = ast.parse(content, filename=str(file_path))
        analyzer = ComplexityAnalyzer()
        analyzer.current_file = str(file_path)
        analyzer.visit(tree)

        return analyzer.function_complexities

    except (SyntaxError, UnicodeDecodeError) as e:
        print(f"Warning: Could not analyze {file_path}: {e}")
        return []


def find_python_files(root_dir: Path) -> List[Path]:
    """Find all Python files in the directory."""
    python_files = []

    # Skip common directories that don't need complexity checking
    skip_dirs = {
        ".git",
        ".venv",
        "__pycache__",
        ".pytest_cache",
        "node_modules",
        ".tox",
        "build",
        "dist",
        ".eggs",
    }

    for path in root_dir.rglob("*.py"):
        # Skip if any parent directory is in skip_dirs
        if any(part in skip_dirs for part in path.parts):
            continue
        python_files.append(path)

    return python_files


class ChangedFiles(NamedTuple):
    """The changed Python files, and the git ref they were compared against."""

    base_ref: str
    files: List[Path]


class ChangedFilesError(RuntimeError):
    """Raised when the set of changed files cannot be determined.

    This is deliberately not representable as an empty file list: "nothing
    changed" and "we have no idea what changed" call for opposite responses,
    and conflating them is what let --files-changed scan the whole repository
    while announcing itself as fast mode.
    """


def _run_git(args: List[str]) -> "subprocess.CompletedProcess[str]":
    """Run a git command and return it without raising on a non-zero status.

    Args:
        args: Arguments after ``git`` itself.

    Returns:
        The completed process, whatever its exit status. Callers decide what a
        failure means; nothing here falls through to a silent default.

    Raises:
        ChangedFilesError: if git could not be run at all.
    """
    try:
        return subprocess.run(
            ["git", *args],
            capture_output=True,
            text=True,
            check=False,
            timeout=GIT_TIMEOUT_SECONDS,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ChangedFilesError(f"could not run `git {' '.join(args)}`: {exc}") from exc


def _base_ref_candidates() -> List[str]:
    """Return the refs to compare against, most authoritative first.

    On a GitHub ``pull_request`` run, GITHUB_BASE_REF holds the branch the PR
    targets, which is not necessarily ``main``: this workflow also runs on
    stacked PRs, whose base is another feature branch. A stacked PR diffed
    against ``main`` reports the union of its own changes and its parent's, so
    when the real base is known it is the only acceptable answer -- there is no
    fallback to ``main`` from here.
    """
    base_branch = os.environ.get("GITHUB_BASE_REF", "").strip()
    if base_branch:
        return [f"origin/{base_branch}", base_branch]

    # Not a pull_request run: a developer invoking this by hand, or a push or
    # schedule run that passed --files-changed anyway.
    return ["origin/main", "main", "HEAD~1"]


def _resolve_base_ref() -> str:
    """Return the first candidate base ref that exists in this checkout.

    Returns:
        A ref name that ``git diff`` can resolve.

    Raises:
        ChangedFilesError: if none of the candidates exist, which on CI almost
            always means the checkout is too shallow to contain the base branch.
    """
    candidates = _base_ref_candidates()
    for ref in candidates:
        # --verify --quiet exits non-zero rather than printing when the ref is
        # missing, which is what separates "no such ref" from "diff failed".
        probe = _run_git(["rev-parse", "--verify", "--quiet", f"{ref}^{{commit}}"])
        if probe.returncode == 0:
            return ref

    base_branch = os.environ.get("GITHUB_BASE_REF", "").strip()
    if base_branch:
        hint = (
            f"GITHUB_BASE_REF is {base_branch!r}, so this is a pull_request run "
            "whose checkout does not contain its own base branch. Give "
            "actions/checkout `fetch-depth: 0`."
        )
    else:
        hint = (
            "GITHUB_BASE_REF is unset. Set it to the branch to compare against, "
            "or run this from a checkout that has an origin/main ref."
        )
    raise ChangedFilesError(
        f"no base ref found (tried: {', '.join(candidates)}). {hint}"
    )


def get_changed_files() -> ChangedFiles:
    """Get the Python files changed relative to this branch's base.

    Returns:
        The base ref that was used, and the changed Python files that still
        exist on disk. An empty file list means exactly that: the diff was
        computed and touched no Python file.

    Raises:
        ChangedFilesError: if the base ref cannot be resolved or the diff fails.
    """
    base_ref = _resolve_base_ref()

    result = _run_git(["diff", "--name-only", f"{base_ref}...HEAD"])
    if result.returncode != 0:
        raise ChangedFilesError(
            f"`git diff --name-only {base_ref}...HEAD` exited "
            f"{result.returncode}: {result.stderr.strip()}"
        )

    changed = []
    for line in result.stdout.splitlines():
        name = line.strip()
        if not name.endswith(".py"):
            continue
        path = Path(name)
        # Deleted files are in the diff but cannot be analysed.
        if path.exists():
            changed.append(path)

    return ChangedFiles(base_ref=base_ref, files=changed)


def generate_report(results: List[ComplexityResult], threshold: int) -> Dict:
    """Generate complexity report."""
    violations = [r for r in results if r.complexity > threshold]

    # Calculate statistics
    if results:
        avg_complexity = sum(r.complexity for r in results) / len(results)
        max_complexity = max(r.complexity for r in results)
    else:
        avg_complexity = 0
        max_complexity = 0

    # Complexity distribution
    complexity_ranges = {
        "1-5 (Low)": 0,
        "6-10 (Moderate)": 0,
        "11-15 (High)": 0,
        "16-20 (Very High)": 0,
        "21+ (Critical)": 0,
    }

    for result in results:
        if result.complexity <= 5:
            complexity_ranges["1-5 (Low)"] += 1
        elif result.complexity <= 10:
            complexity_ranges["6-10 (Moderate)"] += 1
        elif result.complexity <= 15:
            complexity_ranges["11-15 (High)"] += 1
        elif result.complexity <= 20:
            complexity_ranges["16-20 (Very High)"] += 1
        else:
            complexity_ranges["21+ (Critical)"] += 1

    return {
        "timestamp": datetime.now().isoformat(),
        "threshold": threshold,
        "total_functions": len(results),
        "total_violations": len(violations),
        "average_complexity": round(avg_complexity, 2),
        "max_complexity": max_complexity,
        "complexity_distribution": complexity_ranges,
        "high_complexity_functions": sorted(
            [
                {
                    "file": r.file,
                    "line": r.line,
                    "function": r.function,
                    "complexity": r.complexity,
                }
                for r in violations
            ],
            key=lambda x: x["complexity"],
            reverse=True,
        ),
    }


def print_violations(report: Dict) -> None:
    """Print the worst offenders in a report.

    Lives outside the ``--quiet`` guard in :func:`main` on purpose. The flag
    advertises "suppress output except violations", so this listing is the one
    thing it must not silence.
    """
    print("Top violations:")
    for i, func in enumerate(
        report["high_complexity_functions"][:TOP_VIOLATIONS_SHOWN], 1
    ):
        print(
            f"  {i}. {func['file']}:{func['line']} - "
            f"{func['function']} ({func['complexity']})"
        )


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Monitor cyclomatic complexity")
    parser.add_argument(
        "--threshold", type=int, default=10, help="Complexity threshold (default: 10)"
    )
    parser.add_argument("--save", action="store_true", help="Save report to file")
    parser.add_argument(
        "--no-save", action="store_true", help="Don't save report to file"
    )
    parser.add_argument(
        "--quiet", action="store_true", help="Suppress output except violations"
    )
    parser.add_argument(
        "--files-changed",
        action="store_true",
        help=(
            "Only analyze files changed against this branch's base -- "
            "$GITHUB_BASE_REF on a pull_request run, otherwise origin/main "
            f"(speeds up PR checks; exits {EXIT_CHANGED_FILES_UNKNOWN} if the "
            "base cannot be resolved)"
        ),
    )

    args = parser.parse_args()

    # Find all Python files
    root_dir = Path(".")
    if (root_dir / "thyra").exists():
        # Focus on the thyra package
        python_files = find_python_files(root_dir / "thyra")
    else:
        python_files = find_python_files(root_dir)

    # Filter to only changed files if requested
    if args.files_changed:
        try:
            changed = get_changed_files()
        except ChangedFilesError as e:
            print(
                f"ERROR: --files-changed cannot tell which files changed: {e}",
                file=sys.stderr,
            )
            print(
                "ERROR: refusing to scan the whole repository and call it fast mode.",
                file=sys.stderr,
            )
            return EXIT_CHANGED_FILES_UNKNOWN

        changed_files_set = set(changed.files)
        python_files = [f for f in python_files if f in changed_files_set]
        if not args.quiet:
            print(f"Comparing against base ref: {changed.base_ref}")
            print(f"Analyzing {len(python_files)} changed files (--files-changed mode)")

    if not python_files:
        if not args.quiet:
            print("No Python files found to analyze")
        return 0

    # Analyze all files
    all_results = []
    for file_path in python_files:
        results = analyze_file(file_path)
        all_results.extend(results)

    # Generate report
    report = generate_report(all_results, args.threshold)

    violations = report["total_violations"]

    # Print summary
    if not args.quiet:
        print(
            f"Analyzed {len(python_files)} files, {report['total_functions']} functions"
        )
        print(f"Complexity threshold: {args.threshold}")
        print(f"Violations found: {violations}")

        if violations > 0:
            print(f"Average complexity: {report['average_complexity']}")
            print(f"Maximum complexity: {report['max_complexity']}")

    # Not guarded by --quiet: the flag reads "Suppress output except
    # violations", and until now it suppressed these too, which left a
    # --quiet run with no way to say anything at all.
    if violations > 0:
        if not args.quiet:
            print()  # separate the listing from the statistics above it
        print_violations(report)

    # Save report if requested
    if args.save and not args.no_save:
        reports_dir = Path("reports/complexity")
        reports_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = reports_dir / f"complexity_report_{timestamp}.json"

        with open(report_file, "w") as f:
            json.dump(report, f, indent=2)

        if not args.quiet:
            print(f"\nReport saved to: {report_file}")

    # Return exit code based on violations
    return 1 if report["total_violations"] > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
