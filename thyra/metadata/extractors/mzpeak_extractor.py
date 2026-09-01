# thyra/metadata/extractors/mzpeak_extractor.py
"""Metadata extraction for mzPeak archives.

Everything here comes from two places: the per-spectrum columns of
``spectra_metadata.parquet``, and the file-level JSON blobs the archive
carries (Parquet key-value footer, merged with the index's ``metadata``
object -- see :meth:`MzPeakArchive.file_level_metadata`).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import numpy as np

from ...core.base_extractor import MetadataExtractor
from ...resampling.constants import ImzMLAccessions, SpectrumType
from ..types import ComprehensiveMetadata, EssentialMetadata
from .imzml_extractor import UM_PER_UNIT

if TYPE_CHECKING:  # pragma: no cover - typing only
    # Imported for typing only. A runtime import here would be a cycle:
    # thyra.metadata.extractors is imported while thyra.core.base_reader is
    # still initialising, and reaching into thyra.readers from this module
    # pulls every reader package in on top of that half-built module.
    from ...readers.mzpeak.mzpeak_reader import MzPeakArchive

logger = logging.getLogger(__name__)

#: Grid extent terms, reported by the reference converter alongside pixel
#: size. Not used to size the grid -- the observed positions are, so that a
#: dataset whose declared extent disagrees with its own pixels converts to
#: what it actually contains -- but recorded as provenance.
IMS_MAX_COUNT_PIXELS_X = "IMS:1000042"
IMS_MAX_COUNT_PIXELS_Y = "IMS:1000043"


class MzPeakMetadataExtractor(MetadataExtractor):
    """Metadata extractor for mzPeak imaging archives."""

    def __init__(self, archive: "MzPeakArchive", data_path: Path):
        """Initialise the extractor.

        Args:
            archive: An open archive, already validated as imaging.
            data_path: Path to the ``.mzpeak`` file, reported as the source.
        """
        super().__init__(archive)
        self.archive = archive
        self.data_path = Path(data_path)

    # ------------------------------------------------------------------
    # Essential
    # ------------------------------------------------------------------

    def _extract_essential_impl(self) -> EssentialMetadata:
        """Extract the metadata the converter needs to make decisions."""
        index = self.archive.spatial_index()
        coordinates = index.coordinates
        raw = index.raw_positions

        # Grid is sized from what the file actually contains. mzPeak also
        # declares IMS:1000042/43, but a declared extent that disagrees with
        # the positions would silently pad or clip the output, so the
        # declaration is kept as provenance only.
        dimensions = (
            int(coordinates[:, 0].max()) + 1,
            int(coordinates[:, 1].max()) + 1,
            1,
        )
        coordinate_bounds = (
            float(raw[:, 0].min()),
            float(raw[:, 0].max()),
            float(raw[:, 1].min()),
            float(raw[:, 1].max()),
        )

        n_spectra = int(index.spectrum_indices.size)
        total_peaks = self._total_peaks()

        return EssentialMetadata(
            dimensions=dimensions,
            coordinate_bounds=coordinate_bounds,
            mass_range=self._mass_range(),
            pixel_size=self._pixel_size(),
            n_spectra=n_spectra,
            total_peaks=total_peaks,
            estimated_memory_gb=self._estimated_memory_gb(dimensions, total_peaks),
            source_path=str(self.data_path),
            coordinate_offsets=(index.offsets[0], index.offsets[1], 0),
            spectrum_type=self._spectrum_type(),
            peak_counts_per_pixel=self._peak_counts_per_pixel(dimensions),
        )

    def _mass_range(self) -> Tuple[float, float]:
        """Observed m/z range across the archive.

        Read from the per-spectrum ``lowest_observed_mz`` /
        ``highest_observed_mz`` columns rather than by scanning the point
        data: the metadata member has one row per spectrum, so this is a
        9-row read on a 36k-point file and stays a one-row-per-spectrum read
        at any scale.
        """
        table = self.archive.parquet("spectrum", "metadata").read(
            columns=["lowest_observed_mz", "highest_observed_mz"]
        )
        low = np.asarray(table.column("lowest_observed_mz").to_numpy(), dtype=float)
        high = np.asarray(table.column("highest_observed_mz").to_numpy(), dtype=float)
        low = low[np.isfinite(low)]
        high = high[np.isfinite(high)]
        if low.size == 0 or high.size == 0:
            raise ValueError(
                f"{self.data_path} declares no observed m/z range in its "
                f"spectrum metadata."
            )
        return (float(low.min()), float(high.max()))

    def _scan_settings_parameters(self) -> List[dict]:
        """Flatten every parameter of every scan-settings block."""
        settings = self.archive.file_level_metadata().get("scan_settings_list")
        if not isinstance(settings, list):
            return []
        parameters: List[dict] = []
        for block in settings:
            if not isinstance(block, dict):
                continue
            for parameter in block.get("parameters", []) or []:
                if isinstance(parameter, dict):
                    parameters.append(parameter)
        return parameters

    def _pixel_size(self) -> Optional[Tuple[float, float]]:
        """Pixel size in micrometres, or ``None`` when the file omits it.

        Matched on accession, never on name: the reference archive spells
        IMS:1000046 ``"pixel size (x)"`` with parentheses and IMS:1000047
        ``"pixel size y"`` without, so a name match finds one axis and misses
        the other.

        Absent far more often than present -- real Bruker exports carry no
        scan-settings block at all -- in which case this returns ``None`` and
        the converter falls through to ``--pixel-size`` exactly as it does
        for imzML.
        """
        values: Dict[str, float] = {}
        for parameter in self._scan_settings_parameters():
            accession = parameter.get("accession")
            if accession not in (
                ImzMLAccessions.PIXEL_SIZE_X,
                ImzMLAccessions.PIXEL_SIZE_Y,
            ):
                continue
            converted = self._to_micrometres(parameter, accession)
            if converted is not None:
                values.setdefault(accession, converted)

        x = values.get(ImzMLAccessions.PIXEL_SIZE_X)
        y = values.get(ImzMLAccessions.PIXEL_SIZE_Y)
        if x is None and y is None:
            logger.info("Pixel size not found in metadata of %s", self.data_path.name)
            return None
        # A file declaring only one axis is describing a square pixel; imzML
        # writers do this too, and refusing it would reject valid data.
        if x is None:
            x = y
        if y is None:
            y = x
        return (float(x), float(y))

    def _to_micrometres(self, parameter: dict, accession: str) -> Optional[float]:
        """Convert one pixel-size parameter to micrometres.

        A unit outside the known table is refused rather than passed through:
        everything downstream of the extractor is micrometres, so an
        unrecognised unit would become a silent scale error.
        """
        value = parameter.get("value")
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            logger.warning(
                "%s in %s has non-numeric value %r; ignoring it",
                accession,
                self.data_path.name,
                value,
            )
            return None

        unit = parameter.get("unit")
        if unit is None:
            # imzML files routinely omit the unit and mean micrometres.
            return float(value)
        if unit not in UM_PER_UNIT:
            logger.warning(
                "%s in %s declares unsupported unit %r; ignoring the pixel "
                "size rather than assuming a scale",
                accession,
                self.data_path.name,
                unit,
            )
            return None
        return float(value) * UM_PER_UNIT[unit]

    def _spectrum_type(self) -> Optional[str]:
        """Spectrum representation, from the file description or the columns.

        ``file_description.contents`` carries MS:1000127/MS:1000128 for the
        run as a whole. When it does not, the per-spectrum
        ``spectrum_representation`` column is consulted, which the reference
        writer fills with the CV name.
        """
        description = self.archive.file_level_metadata().get("file_description")
        if isinstance(description, dict):
            for parameter in description.get("contents", []) or []:
                if not isinstance(parameter, dict):
                    continue
                accession = parameter.get("accession")
                if accession == ImzMLAccessions.CENTROID_SPECTRUM:
                    return SpectrumType.CENTROID
                if accession == ImzMLAccessions.PROFILE_SPECTRUM:
                    return SpectrumType.PROFILE

        try:
            table = self.archive.parquet("spectrum", "metadata").read(
                columns=["spectrum_representation"]
            )
        except (KeyError, ValueError):
            return None
        values = {
            str(value).strip().lower()
            for value in table.column("spectrum_representation").to_pylist()
            if value
        }
        if SpectrumType.CENTROID in values:
            return SpectrumType.CENTROID
        if SpectrumType.PROFILE in values:
            return SpectrumType.PROFILE
        return None

    def _total_peaks(self) -> int:
        """Points that actually carry a value, across the whole archive.

        ``number_of_data_points`` counts the null-pair padding too, so on the
        reference imaging archive it over-reports by 36%. The padding count
        comes free from Parquet statistics, so the correction costs nothing
        and keeps the converter's pre-allocation honest.
        """
        index = self.archive.spatial_index()
        declared = int(index.point_counts.sum())
        nulls = self.archive.null_count()
        if not nulls:
            return declared
        return max(0, declared - int(nulls))

    def _peak_counts_per_pixel(
        self, dimensions: Tuple[int, int, int]
    ) -> Optional[np.ndarray]:
        """Per-pixel point counts for the streaming converter's CSR indptr.

        Returns ``None`` when the archive contains null-pair padding. The
        per-spectrum counts the file records include that padding, and
        attributing it back to individual spectra would need a full pass over
        the point data -- which is exactly the pass this method exists to
        avoid. Handing the converter counts that are too high would size the
        CSR ``indptr`` wrongly, so the honest answer is to decline and let it
        take its two-pass path.
        """
        if self.archive.null_count():
            logger.info(
                "%s uses null-pair padding, so its per-spectrum point counts "
                "include points that carry no value; declining to supply "
                "per-pixel counts so the converter measures them itself",
                self.data_path.name,
            )
            return None

        index = self.archive.spatial_index()
        n_x = dimensions[0]
        counts = np.zeros(n_x * dimensions[1], dtype=np.int32)
        flat = index.coordinates[:, 1] * n_x + index.coordinates[:, 0]
        counts[flat] = index.point_counts.astype(np.int32, copy=False)
        return counts

    @staticmethod
    def _estimated_memory_gb(
        dimensions: Tuple[int, int, int], total_peaks: int
    ) -> float:
        """Rough dense footprint, used only to pick a conversion strategy."""
        del dimensions
        # 8 bytes of m/z plus 8 of intensity per stored point.
        return float(total_peaks * 16) / (1024.0**3)

    # ------------------------------------------------------------------
    # Comprehensive
    # ------------------------------------------------------------------

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        """Extract everything the archive records, for provenance."""
        metadata = self.archive.file_level_metadata()
        return ComprehensiveMetadata(
            essential=self.get_essential(),
            format_specific=self._format_specific(metadata),
            acquisition_params=self._acquisition_params(),
            instrument_info=self._instrument_info(metadata),
            raw_metadata=dict(metadata),
        )

    def _format_specific(self, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """Container facts worth keeping beside the data."""
        entries = self.archive.index.get("files")
        members = (
            [
                {
                    "name": entry.get("name"),
                    "entity_type": entry.get("entity_type"),
                    "data_kind": entry.get("data_kind"),
                }
                for entry in entries
                if isinstance(entry, dict)
            ]
            if isinstance(entries, list)
            else []
        )
        return {
            "container_version": metadata.get("version"),
            "layout": self.archive.layout(),
            "members": members,
            "run": metadata.get("run"),
            "sample_list": metadata.get("sample_list"),
        }

    def _acquisition_params(self) -> Dict[str, Any]:
        """Scan-settings terms, including the declared grid extent."""
        parameters = self._scan_settings_parameters()
        declared: Dict[str, Any] = {}
        for parameter in parameters:
            accession = parameter.get("accession")
            if accession in (IMS_MAX_COUNT_PIXELS_X, IMS_MAX_COUNT_PIXELS_Y):
                declared[str(accession)] = parameter.get("value")
        return {
            "scan_settings": parameters,
            "declared_grid_extent": declared or None,
        }

    def _instrument_info(self, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """Instrument and software lists, as the archive records them."""
        return {
            "instrument_configuration_list": metadata.get(
                "instrument_configuration_list"
            ),
            "software_list": metadata.get("software_list"),
            "data_processing_method_list": metadata.get("data_processing_method_list"),
        }
