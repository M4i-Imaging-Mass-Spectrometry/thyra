# thyra/metadata/validator.py
"""Validate ontology terms in imzML files."""

import logging
import xml.etree.ElementTree as ET  # nosec B405
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

from .ontology.cache import ONTOLOGY

logger = logging.getLogger(__name__)

#: What ``validate_file`` raises for a file it cannot read: an empty, cut
#: short or wrongly encoded file is an ``ET.ParseError`` (an encoding the
#: parser cannot use, a ``ValueError``), and a folder or a locked file an
#: ``OSError``. ``validate_directory`` reports and skips them.
UNREADABLE_ERRORS = (ET.ParseError, OSError, ValueError)


def describe_unreadable(error: BaseException) -> str:
    """One line saying why a file could not be read."""
    if isinstance(error, ET.ParseError):
        return f"not a readable XML file ({error})"
    if isinstance(error, OSError):
        return error.strerror or str(error)
    return str(error) or type(error).__name__


class ImzMLOntologyValidator:
    """Validate ontology terms in imzML files."""

    def __init__(self):
        """Initialize the ontology validator."""
        self.found_terms: Dict[str, int] = {}
        self.unknown_terms: Dict[str, List[str]] = {}

    def validate_file(self, imzml_path: Path) -> Dict[str, Any]:
        """Validate all ontology terms in an imzML file.

        Raises:
            One of ``UNREADABLE_ERRORS`` when the file cannot be read as XML.
        """
        logger.info(f"Validating ontology terms in {imzml_path}")

        # Parse XML
        tree = ET.parse(imzml_path)  # nosec B314
        root = tree.getroot()

        # Find all cvParam elements
        ns = {"mzml": "http://psi.hupo.org/ms/mzml"}
        cv_params = root.findall(".//mzml:cvParam", ns) or root.findall(".//cvParam")

        total_terms = 0
        known_terms = 0
        unknown_terms = 0
        unknown_list: List[Dict[str, Any]] = []
        term_counts: Dict[str, int] = {}

        for cv_param in cv_params:
            accession = cv_param.get("accession", "")
            name = cv_param.get("name", "")
            value = cv_param.get("value", "")

            total_terms += 1

            if accession:
                term_counts[accession] = term_counts.get(accession, 0) + 1

            term = ONTOLOGY.get_term(accession)
            if term:
                known_terms += 1
            else:
                unknown_terms += 1
                unknown_list.append(
                    {
                        "accession": accession,
                        "name": name,
                        "value": value,
                        "validation_url": ONTOLOGY.validate_against_online(accession),
                    }
                )

                self.unknown_terms.setdefault(accession, []).append("unknown")

        results: Dict[str, Any] = {
            "total_terms": total_terms,
            "known_terms": known_terms,
            "unknown_terms": unknown_terms,
            "unknown_list": unknown_list,
            "term_counts": term_counts,
        }
        results["summary"] = self._generate_summary(results)

        return results

    def _generate_summary(self, results: Dict[str, Any]) -> str:
        """Generate a human-readable summary."""
        lines = [
            "Ontology Validation Summary",
            "===========================",
            f"Total CV terms found: {results['total_terms']}",
            # To avoid division by zero if no terms are found
            f"Unique CV terms: {len(results['term_counts'])}",
            (
                (
                    f"Known terms: {results['known_terms']} "
                    f"({100*results['known_terms']/results['total_terms']:.1f}%)"
                )
                if results["total_terms"] > 0
                else "Known terms: 0 (0.0%)"
            ),
            (
                (
                    f"Unknown terms: {results['unknown_terms']} "
                    f"({100*results['unknown_terms']/results['total_terms']:.1f}%)"
                )
                if results["total_terms"] > 0
                else "Unknown terms: 0 (0.0%)"
            ),
        ]

        if results["term_counts"]:
            lines.extend(["", "Most Common Terms:", "------------------"])
            sorted_terms = sorted(
                results["term_counts"].items(),
                key=lambda item: item[1],
                reverse=True,
            )
            for accession, count in sorted_terms[:15]:
                # ``get_term`` returns ``(name, data_type)``, so the name
                # is [0]. Reading [1] printed the XSD type -- "xsd:string"
                # for IMS:1000080, ``None`` for MS:1000040 -- where the
                # term's name belongs. Nothing caught it because this
                # summary was never printed (issue #299).
                term_details = ONTOLOGY.get_term(accession)
                term_name = term_details[0] if term_details else "Unknown Term"
                lines.append(f"- {accession} ({term_name}): {count} times")

        if results["unknown_list"]:
            lines.extend(["", "Unknown Terms:", "--------------"])
            for term in results["unknown_list"][:10]:
                lines.append(f"- {term['accession']}: {term['name']}")
                if term["validation_url"]:
                    lines.append(f"  Check: {term['validation_url']}")

            if len(results["unknown_list"]) > 10:
                lines.append(f"... and {len(results['unknown_list']) - 10} more")

        return "\n".join(lines)

    def validate_directory(self, directory: Path) -> Dict[str, Any]:
        """Validate all imzML files in a directory and its subfolders.

        A file that cannot be read is recorded under ``unreadable_files``
        with the reason, and the other files are still checked. A folder
        whose name ends in ``.imzML`` is not a file and is left out.

        ``all_unknown_terms`` is a list, most often seen first, then by
        accession, so the order is the same on every run;
        ``unknown_term_counts`` holds how often each was seen.
        """
        # Case-insensitive on every platform: ``sample.imzml`` was found on
        # Windows and missed on Linux.
        imzml_files = sorted(
            path
            for path in directory.glob("**/*.imzML", case_sensitive=False)
            if path.is_file()
        )

        per_file_results: Dict[str, Any] = {}
        unreadable_files: Dict[str, str] = {}
        counts: Counter = Counter()

        for imzml_file in imzml_files:
            try:
                results = self.validate_file(imzml_file)
            except UNREADABLE_ERRORS as error:
                reason = describe_unreadable(error)
                logger.warning(f"Skipping {imzml_file}: {reason}")
                unreadable_files[str(imzml_file)] = reason
                continue
            per_file_results[str(imzml_file)] = results
            counts.update(term["accession"] for term in results["unknown_list"])

        ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        return {
            "files_checked": len(per_file_results),
            "all_unknown_terms": [accession for accession, _ in ranked],
            "unknown_term_counts": dict(ranked),
            "unreadable_files": unreadable_files,
            "per_file_results": per_file_results,
        }
