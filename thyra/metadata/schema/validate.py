# thyra/metadata/schema/validate.py
"""Validate MSI metadata documents against the schema.

Validation happens in three layers:

1. ``schema_version`` compatibility (semantic-version rules),
2. structural validation through the pydantic models, and
3. ontology checks: PSI-MS/IMS/UO accessions must exist in the shipped
   tables with a matching label; fields with a designated ontology
   (NCBITaxon, UBERON, CHEBI) warn when the accession's prefix is not
   the designated one.

Errors mean the document does not conform; warnings mean it conforms
but says something suspicious.  ``thyra validate`` exits non-zero only
on errors, so warnings do not break CI pipelines that gate on it.
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from pydantic import ValidationError

from ..ontology.cache import ONTOLOGY
from .models import MSI_METADATA_SCHEMA_VERSION, MSIMetadata, OntologyTerm

logger = logging.getLogger(__name__)

# Prefixes checked against the locally shipped ontology tables.
_LOCAL_PREFIXES = ("MS", "IMS", "UO")

# The designated ontology for each ``*_term`` field, as
# (section, field) -> accepted prefixes.  A term from another ontology
# is a warning, not an error: EFO and friends legitimately cover some
# of the same ground.
EXPECTED_TERM_PREFIXES: Dict[Tuple[str, str], Tuple[str, ...]] = {
    ("sample", "organism_term"): ("NCBITaxon",),
    ("sample", "organism_part_term"): ("UBERON",),
    ("preparation", "matrix_term"): ("CHEBI",),
    ("ms_analysis", "polarity_term"): ("MS",),
    ("ms_analysis", "ionisation_source_term"): ("MS",),
    ("ms_analysis", "analyzer_term"): ("MS",),
}


@dataclass(frozen=True)
class ValidationIssue:
    """One finding from validation.

    Attributes:
        severity: ``"error"`` (does not conform) or ``"warning"``
            (conforms but suspicious).
        location: Dotted path into the document, e.g.
            ``"ms_analysis.pixel_size_um.x"``.
        message: Human-readable description.
    """

    severity: str
    location: str
    message: str


def _check_schema_version(doc: Dict[str, Any]) -> List[ValidationIssue]:
    """Semantic-version compatibility of ``schema_version``."""
    issues: List[ValidationIssue] = []
    version = doc.get("schema_version")
    if version is None:
        issues.append(
            ValidationIssue(
                "error",
                "schema_version",
                "schema_version is missing; a stored document must state "
                "which schema version it conforms to",
            )
        )
        return issues
    if not isinstance(version, str) or version.count(".") != 2:
        issues.append(
            ValidationIssue(
                "error",
                "schema_version",
                f"schema_version {version!r} is not a MAJOR.MINOR.PATCH string",
            )
        )
        return issues

    try:
        major, minor, _ = (int(part) for part in version.split("."))
    except ValueError:
        issues.append(
            ValidationIssue(
                "error",
                "schema_version",
                f"schema_version {version!r} is not a MAJOR.MINOR.PATCH string",
            )
        )
        return issues

    impl_major, impl_minor, _ = (
        int(part) for part in MSI_METADATA_SCHEMA_VERSION.split(".")
    )
    if major != impl_major:
        issues.append(
            ValidationIssue(
                "error",
                "schema_version",
                f"document is schema major version {major}, this Thyra "
                f"implements {MSI_METADATA_SCHEMA_VERSION}",
            )
        )
    elif minor > impl_minor:
        issues.append(
            ValidationIssue(
                "warning",
                "schema_version",
                f"document is schema {version}, newer than the implemented "
                f"{MSI_METADATA_SCHEMA_VERSION}; fields added since may be "
                "rejected as unknown",
            )
        )
    return issues


def _check_term(
    term: Optional[OntologyTerm], location: str, prefixes: Tuple[str, ...]
) -> List[ValidationIssue]:
    """Ontology checks for one ``*_term`` field."""
    if term is None:
        return []

    issues: List[ValidationIssue] = []
    prefix = term.accession.split(":", 1)[0]

    if prefix not in prefixes:
        issues.append(
            ValidationIssue(
                "warning",
                location,
                f"accession {term.accession} is from '{prefix}'; the "
                f"designated ontology here is {' / '.join(prefixes)}",
            )
        )

    if prefix in _LOCAL_PREFIXES:
        entry = ONTOLOGY.terms.get(term.accession)
        if entry is None:
            issues.append(
                ValidationIssue(
                    "error",
                    location,
                    f"accession {term.accession} is not a known {prefix} term",
                )
            )
        elif entry[0].lower() != term.name.lower():
            issues.append(
                ValidationIssue(
                    "warning",
                    location,
                    f"name {term.name!r} does not match the {prefix} label "
                    f"{entry[0]!r} for {term.accession}",
                )
            )
    return issues


def _check_ontology_terms(meta: MSIMetadata) -> List[ValidationIssue]:
    """Ontology checks across every ``*_term`` field of the document."""
    issues: List[ValidationIssue] = []
    for (section, field), prefixes in EXPECTED_TERM_PREFIXES.items():
        term = getattr(getattr(meta, section), field)
        issues.extend(_check_term(term, f"{section}.{field}", prefixes))
    return issues


def validate_document(
    doc: Any,
) -> Tuple[Optional[MSIMetadata], List[ValidationIssue]]:
    """Validate one metadata document.

    Args:
        doc: The parsed document (normally a dict read from a store's
            ``uns["msi_metadata"]`` or from a JSON file).

    Returns:
        ``(model, issues)``.  The model is ``None`` when structural
        validation failed; issues carry everything found, errors first
        is not guaranteed -- filter on ``severity``.
    """
    if not isinstance(doc, dict):
        return None, [
            ValidationIssue(
                "error", "", f"document must be a mapping, got {type(doc).__name__}"
            )
        ]

    issues = _check_schema_version(doc)
    if any(issue.severity == "error" for issue in issues):
        return None, issues

    try:
        meta = MSIMetadata.model_validate(doc)
    except ValidationError as exc:
        for error in exc.errors():
            location = ".".join(str(part) for part in error["loc"])
            issues.append(ValidationIssue("error", location, error["msg"]))
        return None, issues

    issues.extend(_check_ontology_terms(meta))
    return meta, issues
