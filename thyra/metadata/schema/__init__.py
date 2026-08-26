# thyra/metadata/schema/__init__.py
"""The versioned, ontology-mapped MSI metadata schema.

See :mod:`thyra.metadata.schema.models` for the schema itself,
``docs/metadata-schema.md`` for the storage contract, and
``thyra validate`` / ``thyra export-metaspace`` for the CLI.
"""

from .builder import build_msi_metadata
from .metaspace import to_metaspace
from .models import (
    MSI_METADATA_SCHEMA_VERSION,
    MSI_METADATA_UNS_KEY,
    MSAnalysis,
    MSIMetadata,
    OntologyTerm,
    PixelSizeUm,
    Provenance,
    ResolvingPower,
    SampleInformation,
    SamplePreparation,
)
from .store_io import deep_merge, read_msi_metadata_blocks
from .validate import ValidationIssue, validate_document

__all__ = [
    "MSI_METADATA_SCHEMA_VERSION",
    "MSI_METADATA_UNS_KEY",
    "MSAnalysis",
    "MSIMetadata",
    "OntologyTerm",
    "PixelSizeUm",
    "Provenance",
    "ResolvingPower",
    "SampleInformation",
    "SamplePreparation",
    "ValidationIssue",
    "build_msi_metadata",
    "deep_merge",
    "read_msi_metadata_blocks",
    "to_metaspace",
    "validate_document",
]
