# thyra/metadata/extractors/mzpeak_extractor.py
"""Metadata extraction for mzPeak archives.

Everything here comes from two places: the per-spectrum columns of
``spectra_metadata.parquet``, and the file-level JSON blobs the archive
carries (Parquet key-value footer, merged with the index's ``metadata``
object -- see :meth:`MzPeakArchive.file_level_metadata`).
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import numpy as np

from ...core.base_extractor import MetadataExtractor
from ..constants import ImzMLAccessions, SpectrumType
from ..types import ComprehensiveMetadata, EssentialMetadata
from .imzml_extractor import UM_PER_UNIT, has_old_pixel_size_name

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
#:
#: Mapped to plain names rather than kept as accessions because these become
#: keys in the store's ``uns``, and zarr turns a dict key into a directory
#: name. A colon is not a legal path character on Windows, so accession-keyed
#: provenance makes the whole store unwritable there.
GRID_EXTENT_TERMS = {
    "IMS:1000042": "max_count_of_pixels_x",
    "IMS:1000043": "max_count_of_pixels_y",
}

#: Per axis, the declared pixel count and the declared extent. A pixel size
#: that is not a declared (x, y) pair is tested against these before it is
#: believed.
AXIS_COUNT_AND_EXTENT = (
    ("IMS:1000042", "IMS:1000044"),
    ("IMS:1000043", "IMS:1000045"),
)

#: Relative tolerance of that test. The two readings of a value differ by
#: its square root, so both can pass only when the value is within 2% of 1,
#: where they give the same side. Every public file that could be tested
#: meets the area reading exactly; the margin is for an extent that was
#: rounded when written.
EXTENT_TOLERANCE = 0.01


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
            source_path=str(self.data_path),
            coordinate_offsets=(index.offsets[0], index.offsets[1], 0),
            spectrum_type=self._spectrum_type(),
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
        """Pixel size in micrometres, or ``None`` when the file does not say.

        Terms are found by accession, never by name: the reference archive
        spells IMS:1000046 ``"pixel size (x)"`` with parentheses and
        IMS:1000047 ``"pixel size y"`` without, so a name match finds one
        axis and misses the other.

        A declared pair is taken as written. Anything else is one number
        whose meaning depends on the age of the vocabulary -- a length now,
        an area before 2017 -- so it is tested against the declared pixel
        count and extent, see :meth:`_tested_pixel_size`.

        Absent far more often than present -- real Bruker exports carry no
        scan-settings block at all -- in which case this returns ``None`` and
        the converter falls through to ``--pixel-size`` exactly as it does
        for imzML.
        """
        parameters = self._first_by_accession()
        x = parameters.get(ImzMLAccessions.PIXEL_SIZE_X)
        y = parameters.get(ImzMLAccessions.PIXEL_SIZE_Y)
        if x is None and y is None:
            logger.info("Pixel size not found in metadata of %s", self.data_path.name)
            return None
        if x is not None and y is not None and not self._has_old_name(x):
            return self._declared_pair(x, y)
        return self._tested_pixel_size(parameters)

    def _first_by_accession(self) -> Dict[str, dict]:
        """The first scan-settings parameter of each accession."""
        first: Dict[str, dict] = {}
        for parameter in self._scan_settings_parameters():
            accession = parameter.get("accession")
            if isinstance(accession, str):
                first.setdefault(accession, parameter)
        return first

    @staticmethod
    def _has_old_name(parameter: dict) -> bool:
        """Whether IMS:1000046 is spelled as it was when it gave an area.

        The reference converter copies the parameter into the archive
        unchanged, name included.
        """
        return has_old_pixel_size_name(parameter.get("name"))

    def _declared_pair(self, x: dict, y: dict) -> Optional[Tuple[float, float]]:
        """Both axes as the file declares them, in micrometres.

        An axis that cannot be converted is not filled in from the other
        one: that would be a square pixel taken on trust.
        """
        x_um = self._to_micrometres(x, ImzMLAccessions.PIXEL_SIZE_X)
        y_um = self._to_micrometres(y, ImzMLAccessions.PIXEL_SIZE_Y)
        if x_um is None or y_um is None:
            logger.warning(
                "%s gives a pixel size that cannot be read on both axes, so "
                "no pixel size is taken from it. Pass --pixel-size.",
                self.data_path.name,
            )
            return None
        return (x_um, y_um)

    def _tested_pixel_size(
        self, parameters: Dict[str, dict]
    ) -> Optional[Tuple[float, float]]:
        """A pixel size that is not a declared pair, believed only if tested.

        IMS:1000046 without IMS:1000047, or under its old name, is a length
        in a file written since 2017 and an area in one written before. The
        file's own pixel count and extent say which: the side of a pixel
        times the count is the extent. When they are missing, or fit
        neither reading, no pixel size is returned and the conversion asks
        for ``--pixel-size``, as it does for the same file as imzML.
        """
        x = parameters.get(ImzMLAccessions.PIXEL_SIZE_X)
        if x is None:
            logger.warning(
                "%s gives %s and no %s. Before 2017 that term was "
                '"image shape", so it is not read as a pixel size. '
                "Pass --pixel-size.",
                self.data_path.name,
                ImzMLAccessions.PIXEL_SIZE_Y,
                ImzMLAccessions.PIXEL_SIZE_X,
            )
            return None
        scaled = self._value_and_factor(x, ImzMLAccessions.PIXEL_SIZE_X)
        if scaled is None:
            return None

        reading = self._reading(scaled, self._declared_axes(parameters))
        if reading is None:
            logger.warning(
                "%s gives one pixel size, %s = %r, which is a length in "
                "files written since 2017 and an area in older ones. Its "
                "pixel count (IMS:1000042/43) and extent (IMS:1000044/45) "
                "are missing or fit neither reading, so no pixel size is "
                "taken from it. Pass --pixel-size.",
                self.data_path.name,
                ImzMLAccessions.PIXEL_SIZE_X,
                x.get("value"),
            )
            return None

        kind, side = reading
        y = parameters.get(ImzMLAccessions.PIXEL_SIZE_Y)
        if kind == "length" and y is not None:
            return self._declared_pair(x, y)
        if kind == "area":
            logger.warning(
                "%s: %s = %r is read as the AREA of a pixel, the meaning "
                "the term had before 2017, because its square root times "
                "the declared pixel count is the declared extent. Pixel "
                "size taken: %g um. Pass --pixel-size to set another.",
                self.data_path.name,
                ImzMLAccessions.PIXEL_SIZE_X,
                x.get("value"),
                side,
            )
        else:
            logger.info(
                "%s: pixel size %g um, confirmed by the declared pixel "
                "count and extent",
                self.data_path.name,
                side,
            )
        return (side, side)

    def _declared_axes(self, parameters: Dict[str, dict]) -> List[Tuple[float, float]]:
        """``(pixel count, extent in micrometres)`` of each axis that has both."""
        axes: List[Tuple[float, float]] = []
        for count_term, extent_term in AXIS_COUNT_AND_EXTENT:
            count = parameters.get(count_term)
            extent = parameters.get(extent_term)
            if count is None or extent is None:
                continue
            pixels = self._value_and_factor(count, count_term)
            extent_um = self._to_micrometres(extent, extent_term)
            if pixels is None or extent_um is None:
                continue
            if pixels[0] > 0 and extent_um > 0:
                axes.append((pixels[0], extent_um))
        return axes

    @staticmethod
    def _reading(
        scaled: Tuple[float, float], axes: List[Tuple[float, float]]
    ) -> Optional[Tuple[str, float]]:
        """Which reading of one pixel size fits the declared grid.

        Args:
            scaled: The value as written, and the micrometres per unit.
            axes: ``(pixel count, extent in micrometres)`` per axis.

        Returns:
            ``("length", side)`` or ``("area", side)`` with the side in
            micrometres, or ``None`` when no axis can be tested or an axis
            disagrees. One number describes both axes, so every axis that
            can be tested has to agree.
        """
        value, factor = scaled
        if not axes or not math.isfinite(value) or value <= 0:
            return None
        candidates = (
            ("length", value * factor),
            ("area", math.sqrt(value) * factor),
        )
        for kind, side in candidates:
            if all(
                math.isclose(side * count, extent, rel_tol=EXTENT_TOLERANCE)
                for count, extent in axes
            ):
                return (kind, side)
        return None

    def _value_and_factor(
        self, parameter: dict, accession: str
    ) -> Optional[Tuple[float, float]]:
        """One parameter's value, and the micrometres per unit it declares.

        A missing unit is read as micrometres (factor 1), as it is for
        imzML, so a file and the archive made from it give the same pixel
        size. A unit outside the known table is refused rather than passed
        through: everything downstream of the extractor is micrometres, so
        an unrecognised unit would become a silent scale error. Centimetre
        is left out of the table on purpose: the public files that declare
        it name it "micrometer" in the same parameter.
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
            return (float(value), 1.0)
        if unit not in UM_PER_UNIT:
            logger.warning(
                "%s in %s declares unsupported unit %r; ignoring the value "
                "rather than assuming a scale",
                accession,
                self.data_path.name,
                unit,
            )
            return None
        return (float(value), UM_PER_UNIT[unit])

    def _to_micrometres(self, parameter: dict, accession: str) -> Optional[float]:
        """Convert one length parameter to micrometres.

        ``None`` when the value or the unit cannot be read, see
        :meth:`_value_and_factor`.
        """
        scaled = self._value_and_factor(parameter, accession)
        if scaled is None:
            return None
        value, factor = scaled
        return value * factor

    def _spectrum_type(self) -> Optional[str]:
        """Spectrum representation of the signal member that is read.

        The member decides where it can. mzPeak puts every centroid
        spectrum in the peaks member, so an archive read from there is
        centroid. An archive that fills both members is read from the
        profile one, so it is profile.

        Otherwise the file's own declaration is used; see
        :meth:`_declared_spectrum_type`.
        """
        if self.archive.signal_kind() == "peaks":
            return SpectrumType.CENTROID
        if self.archive.holds_both_representations():
            return SpectrumType.PROFILE
        return self._declared_spectrum_type()

    def _declared_spectrum_type(self) -> Optional[str]:
        """Spectrum representation as the archive declares it.

        ``file_description.contents`` carries MS:1000127/MS:1000128 for the
        run as a whole. When it does not, the per-spectrum
        ``spectrum_representation`` column is consulted. The reference
        converter fills that column with the accession, and the CV name is
        accepted too.
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
        if values & {SpectrumType.CENTROID, ImzMLAccessions.CENTROID_SPECTRUM.lower()}:
            return SpectrumType.CENTROID
        if values & {SpectrumType.PROFILE, ImzMLAccessions.PROFILE_SPECTRUM.lower()}:
            return SpectrumType.PROFILE
        return None

    def _total_peaks(self) -> int:
        """Points that actually carry a value, in the member that is read.

        The declared counts follow the signal member: ``number_of_peaks``
        for a centroid archive, ``number_of_data_points`` otherwise.

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
            name = GRID_EXTENT_TERMS.get(str(parameter.get("accession")))
            if name is not None:
                declared[name] = parameter.get("value")
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
