#!/usr/bin/env python3
"""Complexity monitoring script for CI/CD pipeline."""

import argparse
import ast
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, NamedTuple


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


def get_changed_files() -> List[Path]:
    """Get list of changed Python files from git diff.

    Compares current branch against origin/main to find modified files.
    Falls back to comparing against HEAD~1 if origin/main doesn't exist.
    Returns empty list if git operations fail.

    Returns:
        List of Path objects for changed Python files that exist
    """
    # nosec B404: CI tooling that shells out to git with a fixed argument
    # list and no shell. B603 and B607, which cover the calls themselves,
    # are already skipped repository-wide in pyproject.toml.
    import subprocess  # nosec B404

    try:
        # Try comparing against origin/main first
        result = subprocess.run(
            ["git", "diff", "--name-only", "origin/main...HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )

        if result.returncode != 0 or not result.stdout.strip():
            # Fallback: compare against HEAD~1
            result = subprocess.run(
                ["git", "diff", "--name-only", "HEAD~1", "HEAD"],
                capture_output=True,
                text=True,
                check=True,
                timeout=5,
            )

        files = result.stdout.strip().split("\n")
        python_files = []

        for f in files:
            if f and f.endswith(".py"):
                file_path = Path(f)
                if file_path.exists():
                    python_files.append(file_path)

        return python_files

    except (
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
        FileNotFoundError,
    ) as e:
        print(f"Warning: Could not detect changed files via git: {e}")
        print("Falling back to analyzing all files")
        return []


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


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Monitor cyclomatic complexity")
    parser.add_argument(
        "--threshold", type=int, default=10, help="Complexity threshold (default: 10)"
    )
    parser.add_argument("--save", action="store_true", help="Save report to file")
    parser.add_argument(
        "--trends",
        action="store_true",
        help="Generate trend analysis (requires --save)",
    )
    parser.add_argument(
        "--no-save", action="store_true", help="Don't save report to file"
    )
    parser.add_argument(
        "--quiet", action="store_true", help="Suppress output except violations"
    )
    parser.add_argument(
        "--files-changed",
        action="store_true",
        help="Only analyze files changed compared to main branch (speeds up PR checks)",
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
        changed_files = get_changed_files()
        if changed_files:
            # Filter python_files to only include changed files
            changed_files_set = set(changed_files)
            python_files = [f for f in python_files if f in changed_files_set]
            if not args.quiet:
                print(
                    f"Analyzing {len(python_files)} changed files (--files-changed mode)"
                )
        else:
            if not args.quiet:
                print("Warning: Could not detect changed files, analyzing all files")

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

    # Print summary
    if not args.quiet:
        print(
            f"Analyzed {len(python_files)} files, {report['total_functions']} functions"
        )
        print(f"Complexity threshold: {args.threshold}")
        print(f"Violations found: {report['total_violations']}")

        if report["total_violations"] > 0:
            print(f"Average complexity: {report['average_complexity']}")
            print(f"Maximum complexity: {report['max_complexity']}")

            print("\nTop violations:")
            for i, func in enumerate(report["high_complexity_functions"][:5], 1):
                func_info = func["function"]
                print(
                    f"  {i}. {func['file']}:{func['line']} - "
                    f"{func_info} ({func['complexity']})"
                )

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
