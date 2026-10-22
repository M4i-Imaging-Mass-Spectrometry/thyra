# thyra/tools/check_ontology.py
"""Command-line tool to check ontology terms in imzML files."""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import List, Optional

from ..metadata.ontology.cache import ONTOLOGY
from ..utils.output_files import make_parent_folder, output_file_problem


def _create_parser():
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(description="Check ontology terms in imzML files")
    parser.add_argument("input", help="imzML file or directory to check")
    parser.add_argument("--output", help="Save results to JSON file")
    parser.add_argument(
        "--verbose", action="store_true", help="Log progress at INFO level to stderr"
    )
    return parser


def _setup_logging(verbose: bool):
    """Setup logging configuration."""
    if verbose:
        logging.basicConfig(level=logging.INFO)


def _print_file_results(results: dict, args, unknown_list: list):
    """Print results for single file validation.

    The summary is what the validator computed, not a banner and the
    literal string "Test summary for file" -- which is what this printed
    for every run without ``--verbose``, so the term counts, the
    known/unknown split and the most common terms were all discarded
    after being calculated (issue #299).

    ``--verbose`` no longer gates it. The flag turns on INFO logging; it
    is not a request for the tool's only output. ``_generate_summary``
    emits its own "Ontology Validation Summary" banner as its first two
    lines, so nothing needs to be printed around it.
    """
    if not args.output:
        print(results.get("summary", "No summary available."))
        print()

        if unknown_list:
            print(f"Found {len(unknown_list)} unknown terms:")
            for term in unknown_list[:20]:
                print(f"  - {term['accession']}: {term['name']}")
        else:
            print("No unknown terms encountered.")
        print()


def _print_directory_results(results: dict):
    """Print results for directory validation.

    ``all_unknown_terms`` comes most often seen first, so the list under
    "Most common" is that, in the same order on every run. It used to be
    a set, printed in hash order (issue #529).
    """
    print(f"Checked {results['files_checked']} files")
    print(f"Found {len(results['all_unknown_terms'])} unique unknown terms")

    if results["all_unknown_terms"]:
        counts = results.get("unknown_term_counts", {})
        print("\nMost common unknown terms:")
        for term in results["all_unknown_terms"][:20]:
            seen = f" (seen {counts[term]} times)" if term in counts else ""
            print(f"  - {term}{seen}")
        print()

    unreadable = results.get("unreadable_files", {})
    if unreadable:
        print(f"Could not read {len(unreadable)} files:")
        for path, reason in unreadable.items():
            print(f"  - {path}: {reason}")
        print()


def _save_results_if_requested(results: dict, output_path: str):
    """Save results to JSON file if requested."""
    if output_path:
        make_parent_folder(Path(output_path))
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        print(f"Results saved to {output_path}")
        print(ONTOLOGY.report_unknown_terms())
    else:
        print(ONTOLOGY.report_unknown_terms())


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point for ontology validation tool.

    Returns:
        0 when every file was read, 1 when a file could not be read, and
        2 when the input does not exist or ``--output`` is a folder or lies
        under a file. A path that does not exist used to be searched as an
        empty folder: "Checked 0 files", exit 0.
    """
    parser = _create_parser()
    args = parser.parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: path {input_path} does not exist.", file=sys.stderr)
        return 2
    # Before the check runs, so a path no folder creation can fix does not
    # throw a long folder check away at the end (issue #525).
    if args.output:
        problem = output_file_problem(Path(args.output))
        if problem is not None:
            print(f"Error: {problem}", file=sys.stderr)
            return 2
    _setup_logging(args.verbose)

    from ..metadata.validator import (
        UNREADABLE_ERRORS,
        ImzMLOntologyValidator,
        describe_unreadable,
    )

    validator = ImzMLOntologyValidator()

    if input_path.is_file():
        try:
            results = validator.validate_file(input_path)
        except UNREADABLE_ERRORS as error:
            reason = describe_unreadable(error)
            print(f"Error: cannot read {input_path}: {reason}", file=sys.stderr)
            return 1
        # validate_file reports `unknown_terms` as a count; the terms
        # themselves are in `unknown_list`.
        unknown_list = results.get("unknown_list", [])
        _print_file_results(results, args, unknown_list)
    else:
        results = validator.validate_directory(input_path)
        _print_directory_results(results)

    _save_results_if_requested(results, args.output)
    return 1 if results.get("unreadable_files") else 0


if __name__ == "__main__":
    raise SystemExit(main())
