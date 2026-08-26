# thyra/metadata/schema/models.py
"""Pydantic models for the versioned MSI metadata schema.

The schema is the structured, ontology-mapped description of an MSI
dataset that a converted SpatialData store carries in
``table.uns["msi_metadata"]``.  Its base fields mirror the METASPACE
submission form (organism, organ, condition, sample preparation,
matrix, polarity, ionisation source, analyzer, resolving power, pixel
size) so that the metadata a store carries is sufficient for a
METASPACE submission without retyping -- see
:mod:`thyra.metadata.schema.metaspace`.

Ontology mapping: free-text fields are paired with an optional
:class:`OntologyTerm` (``*_term``) carrying a CURIE accession --
NCBITaxon for organism, UBERON for organism part, CHEBI for the MALDI
matrix, and PSI-MS for polarity, ionisation source and analyzer.  The
PSI-MS/IMS/UO accessions are resolvable against the tables shipped in
:mod:`thyra.metadata.ontology`.

Versioning: ``schema_version`` follows semantic-version rules.  A
reader implementing major version N must reject documents with a
different major version, accept documents with an older minor version,
and may warn on a newer minor version.  Purely additive optional
fields bump the minor version; anything else bumps the major version.
The JSON Schema rendering of these models is committed next to this
module and kept in sync by a unit test; regenerate it with
``python -m thyra.metadata.schema.generate``.
"""

import json
from typing import Any, Dict, List, Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator

# The schema version this code implements and writes.
MSI_METADATA_SCHEMA_VERSION = "0.1.0"

# Where the block lives inside a converted store:
# ``table.uns["msi_metadata"]``.  This location is a stable contract
# (like ``uns["essential_metadata"]`` and the ``coordinate_systems``
# root attribute); consumers read it from here and nowhere else.
MSI_METADATA_UNS_KEY = "msi_metadata"

# The committed JSON Schema artifact for this schema version.
SCHEMA_JSON_FILENAME = "msi_metadata_schema_v0_1.json"

# Fixed var column conventions for the MSI table.  ``mz`` is required
# and written by every converter; the remaining names are reserved for
# annotation results so downstream consumers can rely on one spelling
# (see docs/metadata-schema.md).  Nothing may reuse these names with a
# different meaning; ``thyra validate`` checks the ``mz`` contract.
MSI_VAR_REQUIRED_COLUMNS = ("mz",)
MSI_VAR_RESERVED_COLUMNS = (
    "mz",
    "formula",
    "adduct",
    "annotation_source",
    "fdr",
)

_CURIE_PATTERN = r"^[A-Za-z_][A-Za-z0-9_]*:\S+$"
_SEMVER_PATTERN = r"^\d+\.\d+\.\d+$"


class _SchemaModel(BaseModel):
    """Base for all schema models.

    Unknown keys are rejected rather than silently carried, so a typo
    in a hand-authored document fails validation instead of being
    stored as an unread field.
    """

    model_config = ConfigDict(extra="forbid")


class OntologyTerm(_SchemaModel):
    """A controlled-vocabulary reference: CURIE accession plus label."""

    accession: str = Field(
        pattern=_CURIE_PATTERN,
        description=(
            "CURIE identifier, e.g. 'MS:1000075', 'NCBITaxon:10090', "
            "'UBERON:0002107' or 'CHEBI:90695'."
        ),
    )
    name: str = Field(
        min_length=1,
        description="The term's primary label in its ontology.",
    )


class PixelSizeUm(_SchemaModel):
    """In-plane raster pitch in micrometres."""

    x: float = Field(gt=0, description="Pixel size along x in micrometres.")
    y: float = Field(gt=0, description="Pixel size along y in micrometres.")


class ResolvingPower(_SchemaModel):
    """Mass resolving power quoted at a reference m/z."""

    value: float = Field(gt=0, description="Resolving power (FWHM definition).")
    at_mz: float = Field(gt=0, description="The m/z the value is quoted at.")


class SampleInformation(_SchemaModel):
    """What the sample is (METASPACE ``Sample_Information``)."""

    organism: Optional[str] = Field(
        default=None, description="Species, e.g. 'Mus musculus'."
    )
    organism_term: Optional[OntologyTerm] = Field(
        default=None, description="NCBITaxon term for the organism."
    )
    organism_part: Optional[str] = Field(
        default=None, description="Organ or organism part, e.g. 'liver'."
    )
    organism_part_term: Optional[OntologyTerm] = Field(
        default=None, description="UBERON term for the organism part."
    )
    condition: Optional[str] = Field(
        default=None, description="E.g. 'wildtype', 'diseased'."
    )
    sample_growth_conditions: Optional[str] = Field(
        default=None, description="E.g. intervention or treatment."
    )


class SamplePreparation(_SchemaModel):
    """How the sample was prepared (METASPACE ``Sample_Preparation``)."""

    sample_stabilisation: Optional[str] = Field(
        default=None, description="Preservation method, e.g. 'fresh frozen'."
    )
    tissue_modification: Optional[str] = Field(
        default=None, description="Chemical modification, e.g. 'derivatised'."
    )
    matrix: Optional[str] = Field(
        default=None,
        description=(
            "MALDI matrix, e.g. '2,5-dihydroxybenzoic acid (DHB)'. "
            "'none' for matrix-free techniques (DESI, SIMS)."
        ),
    )
    matrix_term: Optional[OntologyTerm] = Field(
        default=None, description="CHEBI term for the matrix compound."
    )
    matrix_application: Optional[str] = Field(
        default=None, description="Matrix application device or protocol."
    )
    solvent: Optional[str] = Field(default=None, description="Solvent used.")


class MSAnalysis(_SchemaModel):
    """How the data was acquired (METASPACE ``MS_Analysis``).

    ``pixel_size_um`` is the one field that is always known at
    conversion time (conversion refuses to run without a pixel size),
    so it is the one required field of the section.
    """

    polarity: Optional[Literal["positive", "negative"]] = Field(
        default=None, description="Ion polarity mode."
    )
    polarity_term: Optional[OntologyTerm] = Field(
        default=None,
        description="PSI-MS scan polarity term (MS:1000130 / MS:1000129).",
    )
    ionisation_source: Optional[str] = Field(
        default=None, description="E.g. 'MALDI', 'DESI', 'SIMS'."
    )
    ionisation_source_term: Optional[OntologyTerm] = Field(
        default=None, description="PSI-MS ionisation type term."
    )
    analyzer: Optional[str] = Field(
        default=None, description="E.g. 'TOF', 'Orbitrap', 'FTICR'."
    )
    analyzer_term: Optional[OntologyTerm] = Field(
        default=None, description="PSI-MS mass analyzer type term."
    )
    instrument_model: Optional[str] = Field(
        default=None, description="Instrument model as reported by the source."
    )
    detector_resolving_power: Optional[ResolvingPower] = Field(
        default=None, description="Resolving power at a reference m/z."
    )
    pixel_size_um: PixelSizeUm = Field(
        description="In-plane raster pitch in micrometres."
    )

    @model_validator(mode="after")
    def _polarity_and_term_agree(self) -> "MSAnalysis":
        """A polarity string and term that contradict are worse than either alone."""
        if self.polarity is not None and self.polarity_term is not None:
            expected = {
                "positive": "MS:1000130",
                "negative": "MS:1000129",
            }[self.polarity]
            if self.polarity_term.accession != expected:
                raise ValueError(
                    f"polarity '{self.polarity}' contradicts polarity_term "
                    f"{self.polarity_term.accession} (expected {expected})"
                )
        return self


class SoftwareRef(_SchemaModel):
    """A software agent, the way mzQC records analysis software."""

    name: str = Field(min_length=1, description="Software name, e.g. 'thyra'.")
    version: str = Field(min_length=1, description="Software version.")
    uri: Optional[str] = Field(default=None, description="Homepage or repository URL.")


class ProcessingStep(_SchemaModel):
    """One processing action between the raw data and this store.

    Modeled on mzQC provenance: an ordered list of steps, each naming
    the software that performed it and the parameters it ran with.  The
    converter records its own steps (conversion, mass axis resampling);
    downstream tools append theirs (normalisation, peak picking,
    annotation) when they modify the store.
    """

    name: str = Field(
        min_length=1,
        description=(
            "What was done, e.g. 'conversion', 'mass axis resampling', "
            "'normalisation', 'peak picking', 'annotation'."
        ),
    )
    software: SoftwareRef = Field(description="The software that did it.")
    parameters: Dict[str, Union[str, int, float, bool]] = Field(
        default_factory=dict,
        description="The parameters the step ran with.",
    )


class Provenance(_SchemaModel):
    """Who wrote the block and from what.

    Deliberately small: the store's ``uns["essential_metadata"]`` and
    root attributes already record dimensions, mass range, spectrum
    type and conversion timestamp.  This section records only what is
    needed to interpret the *metadata block itself*, and contains no
    timestamps so that two conversions of the same input produce an
    identical block (the uns parity tests rely on this).
    """

    thyra_version: str = Field(
        min_length=1, description="Thyra version that wrote the block."
    )
    source_format: Optional[str] = Field(
        default=None,
        description="Detected input format, e.g. 'imzml', 'bruker', 'phi'.",
    )
    source_path: Optional[str] = Field(
        default=None, description="Path of the source data at conversion time."
    )
    pixel_size_source: Optional[Literal["default", "manual", "automatic"]] = Field(
        default=None,
        description=(
            "How the pixel size was determined: read from the source "
            "metadata ('automatic'), supplied by the user ('manual'), or "
            "the 1.0 um fallback ('default')."
        ),
    )


class MSIMetadata(_SchemaModel):
    """The versioned MSI metadata document.

    ``sample`` and ``preparation`` cannot be auto-populated from raw
    files and default to empty; ``ms_analysis`` and ``provenance`` are
    written by the converter for every store.
    """

    schema_version: str = Field(
        default=MSI_METADATA_SCHEMA_VERSION,
        pattern=_SEMVER_PATTERN,
        description="Semantic version of the schema this document conforms to.",
    )
    sample: SampleInformation = Field(
        default_factory=SampleInformation,
        description="What the sample is; user-supplied.",
    )
    preparation: SamplePreparation = Field(
        default_factory=SamplePreparation,
        description="How the sample was prepared; user-supplied.",
    )
    ms_analysis: MSAnalysis = Field(
        description="How the data was acquired; auto-populated where possible."
    )
    processing: List[ProcessingStep] = Field(
        default_factory=list,
        description="Ordered processing history, oldest first (mzQC-style).",
    )
    provenance: Provenance = Field(
        description="Who wrote this block and from what source."
    )

    def to_uns_dict(self) -> Dict[str, Any]:
        """Serialise for storage in ``table.uns``.

        ``None`` fields are dropped, and the ``sample`` / ``preparation``
        / ``processing`` sections are omitted entirely when empty --
        following the store convention that a section the source has
        nothing for is omitted rather than written empty, so consumers
        can tell "not available" from "available and empty".

        ``processing`` is stored as a JSON string: it is a list of
        objects, which AnnData/zarr cannot round-trip (the same reason
        ``uns["regions"]`` is JSON).  ``read_msi_metadata_blocks`` and
        ``validate_document`` both decode it transparently.
        """
        data: Dict[str, Any] = self.model_dump(mode="json", exclude_none=True)
        for section in ("sample", "preparation", "processing"):
            if not data.get(section):
                data.pop(section, None)
        if "processing" in data:
            data["processing"] = json.dumps(data["processing"])
        return data
