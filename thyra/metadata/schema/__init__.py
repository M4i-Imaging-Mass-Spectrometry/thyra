"""The versioned, ontology-mapped MSI metadata schema.

See :mod:`thyra.metadata.schema.models` for the schema itself,
``docs/metadata-schema.md`` for the storage contract, and
``thyra validate`` / ``thyra export-metaspace`` for the CLI.

**What this schema covers, and what it deliberately does not** (issue
#278). ``uns`` holds fourteen keys and only one of them, ``msi_metadata``,
is built here. That is not two standards of rigour applied to one kind of
thing; it is three kinds of thing:

* **The schema block**, ``uns["msi_metadata"]`` -- closed, versioned and
  ontology-mapped. Validation refuses a key the models do not define, so a
  typo fails loudly instead of becoming a field nobody reads.
* **Open provenance** -- ``essential_metadata``, ``format_specific``,
  ``acquisition_params``, ``raw_metadata``, ``regions``. These carry what
  the source said, in the source's own vocabulary, and a closed schema is
  the wrong shape for them by construction: the whole point is to preserve
  vendor keys nobody has enumerated. ``raw_metadata`` exists *because* the
  schema refuses unknown keys -- it is where the refused material goes.
* **Computed data** -- ``average_spectrum``, ``mobility_heatmap``,
  ``feature_axis``, ``mobility_grid``, ``mobility_axis``,
  ``mobility_marginal``, ``msms_schedule``, ``demultiplexed_current``,
  ``average_spectrum_per_region``. Arrays and tables that live in ``uns``
  because anndata puts non-matrix data there. They are outputs, not
  descriptions, and a metadata schema has nothing to say about them.

So the boundary is "can a fixed vocabulary describe this", not "was anyone
rigorous here". Every one of the fourteen is documented in
``docs/output-format.md`` with its storage contract, which is what a
consumer reads to know what it may rely on.
"""

from .builder import build_msi_metadata, forget_resolved_table
from .metaspace import to_metaspace
from .models import (
    MSI_METADATA_SCHEMA_VERSION,
    MSI_METADATA_UNS_KEY,
    MSI_VAR_PRECURSOR_COLUMN,
    MSI_VAR_PRECURSOR_INDEX_COLUMN,
    MSI_VAR_REQUIRED_COLUMNS,
    MSI_VAR_RESERVED_COLUMNS,
    Fragmentation,
    IonMobility,
    IsolationWindow,
    MobilityGrid,
    MSAnalysis,
    MSIMetadata,
    OntologyTerm,
    PixelSizeUm,
    ProcessingStep,
    Provenance,
    ResolvingPower,
    SampleInformation,
    SamplePreparation,
    SoftwareRef,
)
from .store_io import deep_merge, read_msi_metadata_blocks
from .validate import ValidationIssue, check_store_var_conventions, validate_document

__all__ = [
    "MSI_METADATA_SCHEMA_VERSION",
    "MSI_METADATA_UNS_KEY",
    "MSI_VAR_PRECURSOR_COLUMN",
    "MSI_VAR_PRECURSOR_INDEX_COLUMN",
    "MSI_VAR_REQUIRED_COLUMNS",
    "MSI_VAR_RESERVED_COLUMNS",
    "Fragmentation",
    "IonMobility",
    "IsolationWindow",
    "MSAnalysis",
    "MSIMetadata",
    "MobilityGrid",
    "OntologyTerm",
    "PixelSizeUm",
    "ProcessingStep",
    "Provenance",
    "ResolvingPower",
    "SampleInformation",
    "SamplePreparation",
    "SoftwareRef",
    "ValidationIssue",
    "build_msi_metadata",
    "check_store_var_conventions",
    "deep_merge",
    "forget_resolved_table",
    "read_msi_metadata_blocks",
    "to_metaspace",
    "validate_document",
]
