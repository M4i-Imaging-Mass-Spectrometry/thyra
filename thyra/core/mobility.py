"""Ion mobility as a spectral coordinate.

Mobility is a second coordinate on a *feature*, next to m/z. It never
enters ``obs``, a coordinate system or a transform: those say where a
pixel is, mobility says what was measured there. The reader contract in
:class:`thyra.core.base_reader.BaseMSIReader` carries it as an optional
fourth array per pixel; this module holds the vocabulary the readers, the
converters and the metadata schema share so they name the axis the same
way -- with the PSI-MS terms mzPeak also uses.
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

#: The two mobility quantities an instrument records, as PSI-MS terms.
INVERSE_REDUCED_MOBILITY_ACCESSION = "MS:1002815"
DRIFT_TIME_ACCESSION = "MS:1002476"

MOBILITY_KIND_NAMES: Dict[str, str] = {
    INVERSE_REDUCED_MOBILITY_ACCESSION: "inverse reduced ion mobility",
    DRIFT_TIME_ACCESSION: "ion mobility drift time",
}

#: PSI-MS array terms a source may bind a per-point mobility array to
#: (imzML ``referenceableParamGroup``, mzPeak column metadata), mapped to
#: the quantity they carry -- ``None`` when the term itself does not say
#: and the unit has to decide.
MOBILITY_ARRAY_TERMS: Dict[str, Tuple[str, Optional[str]]] = {
    "MS:1002893": ("ion mobility array", None),
    "MS:1003007": ("raw ion mobility array", None),
    "MS:1002816": ("mean ion mobility array", None),
    "MS:1003154": ("deconvoluted ion mobility array", None),
    "MS:1003006": (
        "mean inverse reduced ion mobility array",
        INVERSE_REDUCED_MOBILITY_ACCESSION,
    ),
    "MS:1003008": (
        "raw inverse reduced ion mobility array",
        INVERSE_REDUCED_MOBILITY_ACCESSION,
    ),
    "MS:1003155": (
        "deconvoluted inverse reduced ion mobility array",
        INVERSE_REDUCED_MOBILITY_ACCESSION,
    ),
    "MS:1002477": ("mean ion mobility drift time array", DRIFT_TIME_ACCESSION),
    "MS:1003153": ("raw ion mobility drift time array", DRIFT_TIME_ACCESSION),
    "MS:1003156": (
        "deconvoluted ion mobility drift time array",
        DRIFT_TIME_ACCESSION,
    ),
}

#: Units that identify the quantity when the array term does not.
_UNIT_TO_KIND: Dict[str, str] = {
    "MS:1002814": INVERSE_REDUCED_MOBILITY_ACCESSION,  # volt-second per square centimeter
    "UO:0000028": DRIFT_TIME_ACCESSION,  # millisecond
    "UO:0000010": DRIFT_TIME_ACCESSION,  # second
}


def classify_mobility_array(
    array_accession: Optional[str], unit_accession: Optional[str]
) -> Tuple[Optional[str], Optional[str]]:
    """The quantity a mobility array carries, as ``(accession, name)``.

    The array term decides when it can (``MS:1003006`` is 1/K0 by
    definition); a generic ``ion mobility array`` is disambiguated by its
    unit; anything else stays ``(None, None)`` rather than guessed.
    """
    kind = None
    if array_accession in MOBILITY_ARRAY_TERMS:
        kind = MOBILITY_ARRAY_TERMS[array_accession][1]
    if kind is None and unit_accession is not None:
        kind = _UNIT_TO_KIND.get(unit_accession)
    if kind is None:
        return None, None
    return kind, MOBILITY_KIND_NAMES[kind]


@dataclass(frozen=True)
class MobilityAxis:
    """What a reader knows about its mobility dimension.

    ``values`` is the native axis when the source shares one across pixels:
    the 1/K0 of every feature of a continuous imzML export, or the 1/K0 of
    every TIMS scan of a Bruker frame. It is ``None`` when each pixel
    carries its own mobility values (processed imzML), in which case only
    :meth:`BaseMSIReader.iter_mobility_spectra` has them.

    The quantity and unit are PSI-MS terms (``MS:1002815`` inverse reduced
    ion mobility, ``MS:1002476`` drift time; ``MS:1002814`` volt-second per
    square centimeter, ``UO:0000028`` millisecond) so a store, an mzPeak
    archive and an imzML export describe the same axis in the same words.
    """

    kind_accession: Optional[str]
    kind_name: Optional[str]
    unit_accession: Optional[str]
    unit_name: Optional[str]
    array_accession: Optional[str] = None
    values: Optional[NDArray[np.float64]] = None
    acq_range: Optional[Tuple[float, float]] = None
    calibration: Optional[Dict[str, Any]] = None
    source: Optional[str] = None

    @property
    def is_shared(self) -> bool:
        """Whether every pixel is described by the same ``values`` array."""
        return self.values is not None

    def to_uns(self) -> Dict[str, Any]:
        """The axis as a plain ``uns`` block: plain-name keys, no colons.

        Zarr writes a dict key as a directory name and a colon is not a
        legal Windows path character, so CV accessions are values here,
        never keys.
        """
        block: Dict[str, Any] = {"present": True}
        if self.kind_accession is not None:
            block["type_accession"] = self.kind_accession
            block["type_name"] = self.kind_name
        if self.unit_accession is not None:
            block["unit_accession"] = self.unit_accession
        if self.unit_name is not None:
            block["unit_name"] = self.unit_name
        if self.array_accession is not None:
            block["array_accession"] = self.array_accession
        if self.values is not None:
            block["values"] = np.asarray(self.values, dtype=np.float64)
            block["n_values"] = int(self.values.size)
        lo_hi = self.acq_range
        if lo_hi is None and self.values is not None and self.values.size:
            lo_hi = (float(np.min(self.values)), float(np.max(self.values)))
        if lo_hi is not None:
            block["range_lower"] = float(lo_hi[0])
            block["range_upper"] = float(lo_hi[1])
        if self.calibration:
            block["calibration"] = dict(self.calibration)
        if self.source is not None:
            block["source"] = self.source
        return block

    def to_extractor_report(self) -> Dict[str, Any]:
        """The shape :func:`thyra.metadata.schema.builder` reads from ``format_specific``."""
        report: Dict[str, Any] = {"present": True}
        if self.kind_name is not None:
            report["separation"] = self.kind_name
            report["separation_accession"] = self.kind_accession
        if self.unit_name is not None:
            report["unit"] = self.unit_name
        if self.unit_accession is not None:
            report["unit_accession"] = self.unit_accession
        if self.array_accession is not None:
            report["array_accession"] = self.array_accession
        lo_hi = self.acq_range
        if lo_hi is None and self.values is not None and self.values.size:
            lo_hi = (float(np.min(self.values)), float(np.max(self.values)))
        if lo_hi is not None:
            report["range"] = [float(lo_hi[0]), float(lo_hi[1])]
        if self.values is not None:
            report["n_values"] = int(self.values.size)
        report["shared_axis"] = self.is_shared
        return report
