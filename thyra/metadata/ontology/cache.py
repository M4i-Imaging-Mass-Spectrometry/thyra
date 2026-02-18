# thyra/metadata/ontology/cache.py
"""Local ontology cache for imzML terms."""

import logging
from typing import Dict, List, Optional, Set, Tuple

from ._ims import terms as ims_terms
from ._ms import terms as ms_terms
from ._uo import terms as uo_terms

logger = logging.getLogger(__name__)


class OntologyCache:
    """Local cache for MS/IMS/UO ontology terms."""

    def __init__(self):
        """Initialize the ontology cache by loading all term definitions."""
        self.terms: Dict[str, Tuple[str, Optional[str]]] = {}
        self.unknown_terms: Set[str] = set()

        # Load all terms
        self.terms.update(uo_terms)
        self.terms.update(ms_terms)
        self.terms.update(ims_terms)

        # Create reverse lookup
        self.name_to_accession = {
            name.lower(): acc for acc, (name, _) in self.terms.items()
        }

        logger.info(f"Loaded {len(self.terms)} ontology terms")
        logger.info(
            f"  - MS: {len([k for k in self.terms if k.startswith('MS:')])} terms"
        )
        logger.info(
            f"  - IMS: {len([k for k in self.terms if k.startswith('IMS:')])} terms"
        )
        logger.info(
            f"  - UO: {len([k for k in self.terms if k.startswith('UO:')])} terms"
        )

    def get_term(self, accession: str) -> Optional[Tuple[str, Optional[str]]]:
        """Get term name and data type by accession."""
        term = self.terms.get(accession)
        if not term and accession not in self.unknown_terms:
            self.unknown_terms.add(accession)
            logger.warning(f"Unknown ontology term: {accession}")
        return term

    def check_unknown_terms(self) -> Set[str]:
        """Return set of unknown terms encountered."""
        return self.unknown_terms.copy()

    def report_unknown_terms(self) -> str:
        """Generate a report of unknown terms."""
        if not self.unknown_terms:
            return "No unknown terms encountered."

        report = [f"Found {len(self.unknown_terms)} unknown terms:"]

        # Group by prefix
        by_prefix: Dict[str, List[str]] = {}
        for term in sorted(self.unknown_terms):
            prefix = term.split(":")[0] if ":" in term else "OTHER"
            by_prefix.setdefault(prefix, []).append(term)

        for prefix, terms in sorted(by_prefix.items()):
            report.append(f"\n{prefix} terms ({len(terms)}):")
            for term in terms[:10]:  # Show first 10
                report.append(f"  - {term}")
            if len(terms) > 10:
                report.append(f"  ... and {len(terms) - 10} more")

        return "\n".join(report)

    def validate_against_online(self, accession: str) -> Optional[str]:
        """Get URL to validate a term online."""
        prefix = accession.split(":")[0] if ":" in accession else None

        urls = {
            "MS": (
                f"https://www.ebi.ac.uk/ols/ontologies/ms/terms?iri="
                f"http://purl.obolibrary.org/obo/{accession.replace(':', '_')}"
            ),
            "IMS": (
                f"https://www.ebi.ac.uk/ols/ontologies/ims/terms?iri="
                f"http://purl.obolibrary.org/obo/{accession.replace(':', '_')}"
            ),
            "UO": (
                f"https://www.ebi.ac.uk/ols/ontologies/uo/terms?iri="
                f"http://purl.obolibrary.org/obo/{accession.replace(':', '_')}"
            ),
        }

        return urls.get(prefix) if prefix is not None else None


# Global instance
ONTOLOGY = OntologyCache()
