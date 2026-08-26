# thyra/metadata/schema/vocab.py
"""Controlled-vocabulary normalisation for MSI metadata fields.

Vendor files spell the same fact many ways ("positive", "POS", "+",
"positive scan").  These helpers map the spellings the readers actually
produce onto a canonical label plus the PSI-MS term for it, so that the
``msi_metadata`` block is comparable across formats.

Term names are always resolved from the ontology tables shipped in
:mod:`thyra.metadata.ontology`, never hardcoded here, so a label in a
written block matches the local table by construction.

Unrecognised input normalises to ``None`` rather than guessing: an
absent field is honest, a wrong ontology term is not.
"""

from typing import Dict, Optional, Tuple

from ..ontology.cache import ONTOLOGY
from .models import OntologyTerm

POLARITY_ACCESSIONS: Dict[str, str] = {
    "positive": "MS:1000130",  # positive scan
    "negative": "MS:1000129",  # negative scan
}

_POLARITY_ALIASES: Dict[str, str] = {
    "positive": "positive",
    "pos": "positive",
    "+": "positive",
    "positive scan": "positive",
    "positive polarity": "positive",
    "positive ion mode": "positive",
    "negative": "negative",
    "neg": "negative",
    "-": "negative",
    "negative scan": "negative",
    "negative polarity": "negative",
    "negative ion mode": "negative",
}

IONISATION_SOURCE_ACCESSIONS: Dict[str, str] = {
    "MALDI": "MS:1000075",  # matrix-assisted laser desorption ionization
    "ESI": "MS:1000073",  # electrospray ionization
    "DESI": "MS:1002011",  # desorption electrospray ionization
    "SIMS": "MS:1000402",  # secondary ionization
}

_IONISATION_SOURCE_ALIASES: Dict[str, str] = {
    "maldi": "MALDI",
    "matrix-assisted laser desorption ionization": "MALDI",
    "matrix assisted laser desorption ionization": "MALDI",
    "esi": "ESI",
    "electrospray": "ESI",
    "electrospray ionization": "ESI",
    "electrospray ionisation": "ESI",
    "desi": "DESI",
    "desorption electrospray ionization": "DESI",
    "desorption electrospray ionisation": "DESI",
    "sims": "SIMS",
    "tof-sims": "SIMS",
    "tof sims": "SIMS",
    "secondary ionization": "SIMS",
    "secondary ionisation": "SIMS",
    "secondary ion mass spectrometry": "SIMS",
}

ANALYZER_ACCESSIONS: Dict[str, str] = {
    "TOF": "MS:1000084",  # time-of-flight
    "Orbitrap": "MS:1000484",  # orbitrap
    "FTICR": "MS:1000079",  # fourier transform ion cyclotron resonance MS
    "Ion trap": "MS:1000264",  # ion trap
}

_ANALYZER_ALIASES: Dict[str, str] = {
    "tof": "TOF",
    "time-of-flight": "TOF",
    "time of flight": "TOF",
    "qtof": "TOF",
    "q-tof": "TOF",
    "linear tof": "TOF",
    "reflector tof": "TOF",
    "orbitrap": "Orbitrap",
    "fticr": "FTICR",
    "ft-icr": "FTICR",
    "fourier transform ion cyclotron resonance mass spectrometer": "FTICR",
    "ion trap": "Ion trap",
}


def term_from_accession(accession: str) -> OntologyTerm:
    """Build an :class:`OntologyTerm` from the local ontology tables.

    Args:
        accession: A PSI-MS / IMS / UO accession, e.g. ``"MS:1000075"``.

    Returns:
        The term with its name taken from the shipped table.

    Raises:
        KeyError: If the accession is not in the local tables.  The
            vocabularies in this module only reference accessions that
            are, so hitting this means the tables and this module have
            drifted.
    """
    entry = ONTOLOGY.terms.get(accession)
    if entry is None:
        raise KeyError(f"Accession {accession} is not in the local ontology tables")
    return OntologyTerm(accession=accession, name=entry[0])


def _normalize(value: object) -> Optional[str]:
    """Lower-cased, stripped string form of a raw metadata value."""
    if not isinstance(value, str):
        return None
    stripped = value.strip().lower()
    return stripped or None


def normalize_polarity(value: object) -> Optional[Tuple[str, OntologyTerm]]:
    """Map a raw polarity value to ``("positive"|"negative", term)``.

    Returns ``None`` for anything unrecognised, including mixed-mode
    values -- a dataset alternating polarities has no single truthful
    value for this field.
    """
    key = _normalize(value)
    if key is None:
        return None
    canonical = _POLARITY_ALIASES.get(key)
    if canonical is None:
        return None
    return canonical, term_from_accession(POLARITY_ACCESSIONS[canonical])


def normalize_ionisation_source(value: object) -> Optional[Tuple[str, OntologyTerm]]:
    """Map a raw ionisation source value to ``(label, term)``."""
    key = _normalize(value)
    if key is None:
        return None
    canonical = _IONISATION_SOURCE_ALIASES.get(key)
    if canonical is None:
        return None
    return canonical, term_from_accession(IONISATION_SOURCE_ACCESSIONS[canonical])


def normalize_analyzer(value: object) -> Optional[Tuple[str, OntologyTerm]]:
    """Map a raw mass analyzer value to ``(label, term)``."""
    key = _normalize(value)
    if key is None:
        return None
    canonical = _ANALYZER_ALIASES.get(key)
    if canonical is None:
        return None
    return canonical, term_from_accession(ANALYZER_ACCESSIONS[canonical])
