# thyra/tools/check_ontology.py
"""Command-line tool to check ontology terms in imzML files."""

import argparse
import json
import logging
from pathlib import Path

from ..metadata.ontology.cache import ONTOLOGY


def _create_parser():
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(description="Check ontology terms in imzML files")
    parser.add_argument("input", help="imzML file or directory to check")
    parser.add_argument("--output", help="Save results to JSON file")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    return parser


def _setup_logging(verbose: bool):
    """Setup logging configuration."""
    if verbose:
        logging.basicConfig(level=logging.INFO)


def _print_file_results(results: dict, args, unknown_list: list):
    """Print results for single file validation."""
    if not args.output:
        if args.verbose and "summary" in results:
            print(results["summary"])
        else:
            print("Ontology Validation Summary")
            print("===========================")
            print("Test summary for file")
            print()

        if unknown_list:
            print(f"Found {len(unknown_list)} unknown terms:")
            for term in unknown_list[:20]:
                print(f"  - {term['accession']}: {term['name']}")
        else:
            print("No unknown terms encountered.")
        print()


def _print_directory_results(results: dict):
    """Print results for directory validation."""
    print(f"Checked {results['files_checked']} files")
    print(f"Found {len(results['all_unknown_terms'])} unique unknown terms")

    if results["all_unknown_terms"]:
        print("\nMost common unknown terms:")
        for term in list(results["all_unknown_terms"])[:20]:
            print(f"  - {term}")
        print()


def _save_results_if_requested(results: dict, output_path: str):
    """Save results to JSON file if requested."""
    if output_path:
        with open(output_path, "w") as f:
            json.dump(results, f, indent=2, default=str)
        print(f"Results saved to {output_path}")
        print(ONTOLOGY.report_unknown_terms())
    else:
        print(ONTOLOGY.report_unknown_terms())


def main():
    """Main entry point for ontology validation tool."""
    parser = _create_parser()
    args = parser.parse_args()
    _setup_logging(args.verbose)

    from ..metadata.validator import ImzMLOntologyValidator

    validator = ImzMLOntologyValidator()
    input_path = Path(args.input)

    if input_path.is_file():
        results = validator.validate_file(input_path)
        # validate_file reports `unknown_terms` as a count; the terms
        # themselves are in `unknown_list`.
        unknown_list = results.get("unknown_list", [])
        _print_file_results(results, args, unknown_list)
    else:
        results = validator.validate_directory(input_path)
        _print_directory_results(results)

    _save_results_if_requested(results, args.output)


if __name__ == "__main__":
    main()
