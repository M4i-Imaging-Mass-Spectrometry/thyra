# thyra/metadata/extractors/solarix_extractor.py
"""Metadata extractor for Bruker solariX (FT-ICR / MRMS) peaks.sqlite data.

Reads from the reader's already-loaded Properties/AcquisitionKeys tables and
SQL aggregates, so no blobs are decoded for metadata extraction.
"""

import logging
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple

from ...core.base_extractor import MetadataExtractor
from ..types import ComprehensiveMetadata, EssentialMetadata

if TYPE_CHECKING:
    from ...readers.bruker.solarix.solarix_reader import SolarixReader

logger = logging.getLogger(__name__)

# Properties.InstrumentFamily observed on a solariX (serial 1272500.00145);
# the mapping is verified against that one instrument only, so the raw
# family number is always surfaced alongside the derived model name.
_SOLARIX_INSTRUMENT_FAMILY = 513

# Polarity enum in AcquisitionKeys, verified against a paired
# positive/negative acquisition of the same sample.
_POLARITY_NAMES = {0: "positive", 1: "negative"}


def _collapse(pair: Tuple[Any, Any]) -> Any:
    """Collapse a (min, max) aggregate to a scalar when constant."""
    low, high = pair
    if low == high:
        return low
    return [low, high]


class SolarixMetadataExtractor(MetadataExtractor):
    """solariX-specific metadata extractor.

    Args:
        reader: The :class:`SolarixReader` to describe.
    """

    def __init__(self, reader: "SolarixReader"):
        """Initialise the extractor.

        Args:
            reader: The :class:`SolarixReader` to describe.
        """
        super().__init__(reader)
        self._reader = reader

    def _extract_essential_impl(self) -> EssentialMetadata:
        """Extract metadata needed to drive conversion."""
        reader = self._reader
        summary = reader.summary
        n_x, n_y, _ = reader.dimensions

        total_peaks = summary["total_peaks"]
        return EssentialMetadata(
            dimensions=reader.dimensions,
            coordinate_bounds=(0.0, float(n_x - 1), 0.0, float(n_y - 1)),
            mass_range=reader.mass_range,
            pixel_size=reader.pixel_size_um,
            n_spectra=summary["n_spectra"],
            total_peaks=total_peaks,
            estimated_memory_gb=(total_peaks * 2 * 8) / (1024**3),
            source_path=str(reader.data_path),
            # The absolute stage-raster origin subtracted by iter_spectra;
            # XIndexPos starts wherever the stage was, not at 0.
            coordinate_offsets=reader.coordinate_offsets,
            spectrum_type="centroid spectrum",
        )

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        """Extract full metadata including vendor-specific detail."""
        return ComprehensiveMetadata(
            essential=self.get_essential(),
            format_specific=self._extract_format_specific(),
            acquisition_params=self._extract_acquisition_params(),
            instrument_info=self._extract_instrument_info(),
            raw_metadata=self._extract_raw_metadata(),
        )

    def _polarity(self) -> Optional[str]:
        """Polarity name via the verified enum, or None when mixed/unknown.

        Verified against a paired positive/negative acquisition:
        0 = positive, 1 = negative. Anything else stays unnamed and is
        available raw in ``acquisition_keys``.
        """
        raw_values = {
            key["polarity_raw"] for key in self._reader.acquisition_keys.values()
        }
        if len(raw_values) != 1:
            return None
        return _POLARITY_NAMES.get(next(iter(raw_values)))

    def _extract_format_specific(self) -> Dict[str, Any]:
        reader = self._reader
        props = reader.properties
        summary = reader.summary
        return {
            "format": "Bruker solariX peaks.sqlite",
            "schema_type": props.get("SchemaType"),
            "schema_version": (
                f"{props.get('SchemaVersionMajor')}."
                f"{props.get('SchemaVersionMinor')}"
            ),
            "acquisition_software": props.get("AcquisitionSoftware"),
            "acquisition_software_vendor": props.get("AcquisitionSoftwareVendor"),
            "acquisition_software_version": props.get("AcquisitionSoftwareVersion"),
            # Absolute stage-raster extents before the 0-based offset.
            "raster_x_range": [summary["x_min"], summary["x_max"]],
            "raster_y_range": [summary["y_min"], summary["y_max"]],
            "pixel_size_source": reader.pixel_size_source,
            "mis_file": str(reader.mis_path) if reader.mis_path else "",
            "regions": reader.get_region_info() or [],
        }

    def _extract_acquisition_params(self) -> Dict[str, Any]:
        reader = self._reader
        props = reader.properties
        stats = reader.get_acquisition_stats()
        params: Dict[str, Any] = {
            "polarity": self._polarity(),
            "mz_acq_range": list(reader.mass_range),
            "acquisition_datetime": props.get("AcquisitionDateTime"),
            "operator_name": props.get("OperatorName"),
            "laser_power": _collapse(stats["laser_power"]),
            "laser_rep_rate": _collapse(stats["laser_rep_rate"]),
            "num_summations": _collapse(stats["num_summations"]),
        }
        # Raw enum values per acquisition key; the full ScanMode /
        # AcquisitionMode / MsLevel enums are undocumented, so the numbers
        # are stored as found rather than labelled.
        for key_id, key in sorted(reader.acquisition_keys.items()):
            prefix = f"acquisition_key_{key_id}"
            params[f"{prefix}_polarity_raw"] = key["polarity_raw"]
            params[f"{prefix}_scan_mode_raw"] = key["scan_mode_raw"]
            params[f"{prefix}_acquisition_mode_raw"] = key["acquisition_mode_raw"]
            params[f"{prefix}_ms_level_raw"] = key["ms_level_raw"]
        return params

    def _extract_instrument_info(self) -> Dict[str, Any]:
        props = self._reader.properties
        family = props.get("InstrumentFamily")
        try:
            family_number: Optional[int] = int(family)
        except (TypeError, ValueError):
            family_number = None

        info: Dict[str, Any] = {
            # The exact string FTICRDetector matches on -- routes the data
            # to nearest-neighbor resampling on the FTICR axis.
            "instrument_type": "FT-ICR",
            "manufacturer": props.get("InstrumentVendor") or "Bruker",
            "serial_number": props.get("InstrumentSerialNumber"),
            "instrument_family": family_number,
            "instrument_source_type": props.get("InstrumentSourceType"),
        }
        # Family 513 -> solariX is verified against one instrument
        # (serial 1272500.00145); other families keep the raw number only.
        if family_number == _SOLARIX_INSTRUMENT_FAMILY:
            info["instrument_model"] = "solariX"
            info["instrument_name"] = "solariX"
        return info

    def _extract_raw_metadata(self) -> Dict[str, Any]:
        reader = self._reader
        return {
            "properties": dict(reader.properties),
            "acquisition_keys": {
                str(key_id): dict(key)
                for key_id, key in reader.acquisition_keys.items()
            },
            "mis_metadata": dict(reader.mis_metadata),
        }
