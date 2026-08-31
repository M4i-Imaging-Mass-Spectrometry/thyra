# thyra/metadata/schema/metaspace.py
"""Export MSI metadata as a METASPACE submission document.

METASPACE (https://metaspace2020.org) is the public platform for
metabolite annotation of MSI data.  Its submission form requires a
metadata JSON; the schema's base fields mirror that form, so the
export is a straight mapping.

Fields METASPACE requires but the document does not carry are emitted
as empty strings and reported as warnings rather than fabricated --
the output is a truthful starting point the submitter completes, not
a pretend-complete record.  The one inference made: a matrix of
``"none"`` for matrix-free ionisation (DESI, SIMS), which is a fact of
the technique, not a guess.
"""

from typing import Any, Dict, List, Optional, Tuple, Union

from .models import MSIMetadata

METASPACE_DATA_TYPE = "Imaging MS"

_MATRIX_FREE_SOURCES = ("DESI", "SIMS")


def _number(value: float) -> Union[int, float]:
    """Integral floats as ints, matching how the form displays them."""
    return int(value) if float(value).is_integer() else float(value)


def _require(value: Optional[str], label: str, warnings: List[str]) -> str:
    """A required METASPACE field: empty string plus warning when unset."""
    if value is None or not str(value).strip():
        warnings.append(
            f"{label} is required for a METASPACE submission and is not set"
        )
        return ""
    return str(value)


def _matrix_fields(meta: MSIMetadata, warnings: List[str]) -> Tuple[str, str]:
    """``MALDI_Matrix`` and ``MALDI_Matrix_Application``.

    Matrix-free sources truthfully get ``"none"``; otherwise the field
    is required.
    """
    preparation = meta.preparation
    source = meta.ms_analysis.ionisation_source

    if preparation.matrix:
        matrix = preparation.matrix
    elif source in _MATRIX_FREE_SOURCES:
        matrix = "none"
    else:
        matrix = _require(None, "Sample_Preparation.MALDI_Matrix", warnings)

    if preparation.matrix_application:
        application = preparation.matrix_application
    elif matrix == "none":
        application = "none"
    else:
        application = _require(
            None, "Sample_Preparation.MALDI_Matrix_Application", warnings
        )

    return matrix, application


def to_metaspace(meta: MSIMetadata) -> Tuple[Dict[str, Any], List[str]]:
    """Render the METASPACE submission metadata document.

    Args:
        meta: A validated MSI metadata document.

    Returns:
        ``(document, warnings)`` -- the submission JSON structure, and
        one warning per required field that had to be left empty.
    """
    warnings: List[str] = []
    sample = meta.sample
    preparation = meta.preparation
    analysis = meta.ms_analysis

    matrix, matrix_application = _matrix_fields(meta, warnings)

    polarity = ""
    if analysis.polarity is not None:
        polarity = analysis.polarity.capitalize()
    else:
        _require(None, "MS_Analysis.Polarity", warnings)

    document: Dict[str, Any] = {
        "Data_Type": METASPACE_DATA_TYPE,
        "Sample_Information": {
            "Organism": _require(
                sample.organism, "Sample_Information.Organism", warnings
            ),
            "Organism_Part": _require(
                sample.organism_part, "Sample_Information.Organism_Part", warnings
            ),
            "Condition": _require(
                sample.condition, "Sample_Information.Condition", warnings
            ),
            "Sample_Growth_Conditions": sample.sample_growth_conditions or "",
        },
        "Sample_Preparation": {
            "Sample_Stabilisation": preparation.sample_stabilisation or "",
            "Tissue_Modification": preparation.tissue_modification or "",
            "MALDI_Matrix": matrix,
            "MALDI_Matrix_Application": matrix_application,
            "Solvent": preparation.solvent or "none",
        },
        "MS_Analysis": {
            "Polarity": polarity,
            "Ionisation_Source": _require(
                analysis.ionisation_source, "MS_Analysis.Ionisation_Source", warnings
            ),
            "Analyzer": _require(analysis.analyzer, "MS_Analysis.Analyzer", warnings),
            "Pixel_Size": {
                "Xaxis": _number(analysis.pixel_size_um.x),
                "Yaxis": _number(analysis.pixel_size_um.y),
            },
        },
    }

    resolving_power = analysis.detector_resolving_power
    if resolving_power is not None:
        document["MS_Analysis"]["Detector_Resolving_Power"] = {
            "mz": _number(resolving_power.at_mz),
            "Resolving_Power": _number(resolving_power.value),
        }
    else:
        warnings.append(
            "MS_Analysis.Detector_Resolving_Power is required for a METASPACE "
            "submission and is not set"
        )

    return document, warnings
