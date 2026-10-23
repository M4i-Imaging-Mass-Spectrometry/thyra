# thyra/metadata/extractors/bruker_extractor.py
import logging
import math
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ...core.base_extractor import MetadataExtractor
from ...core.positions import whole_sql
from ...errors import ConversionRefused
from ..types import ComprehensiveMetadata, EssentialMetadata

logger = logging.getLogger(__name__)

#: The ``CalibrationInfo`` keys read, and the only ones: see
#: :meth:`BrukerMetadataExtractor._extract_instrument_calibration`.
_CALIBRATION_INFO_KEYS = (
    "CalibrationDateTime",
    "CalibrationSoftware",
    "CalibrationSoftwareVersion",
    "MzStandardDeviationPPM",
    "ReferencePeakMasses",
    "MassesCorrectedCalibration",
)

#: The text keys among them, under the names Thyra writes them.
_CALIBRATION_INFO_TEXT = (
    ("CalibrationDateTime", "calibration_datetime"),
    ("CalibrationSoftware", "calibration_software"),
    ("CalibrationSoftwareVersion", "calibration_software_version"),
)


def _mismatch_reading(raster: Tuple[float, float], beam: Tuple[float, float]) -> str:
    """What a Raster step that differs from BeamScanSize says, by direction.

    Oversampled only when the beam scans more than one step, on some axis
    and on none less; a step larger than the beam leaves gaps between spots.
    """
    if all(b >= r for r, b in zip(raster, beam)):
        return "The beam scan is larger than the step: the acquisition is oversampled."
    if all(b <= r for r, b in zip(raster, beam)):
        return "The beam scan is smaller than the step: the spots leave gaps."
    return "The beam scan is larger than the step on one axis and smaller on the other."


def _float64_count(blob: Any) -> Optional[int]:
    """How many little-endian float64 values a ``CalibrationInfo`` blob holds."""
    if not isinstance(blob, (bytes, bytearray)) or len(blob) % 8:
        return None
    return len(blob) // 8


def _instrument_calibration(stated: Dict[str, Any]) -> Dict[str, Any]:
    """The ``CalibrationInfo`` values Thyra keeps, under Thyra's key names."""
    info: Dict[str, Any] = {}
    for key, name in _CALIBRATION_INFO_TEXT:
        value = stated.get(key)
        if isinstance(value, str) and value.strip():
            info[name] = value.strip()
    stated_sd = stated.get("MzStandardDeviationPPM")
    try:
        sd = math.nan if stated_sd is None else float(stated_sd)
    except (TypeError, ValueError):
        sd = math.nan
    if math.isfinite(sd):
        info["mz_standard_deviation_ppm"] = sd
    n_reference = _float64_count(stated.get("ReferencePeakMasses"))
    if n_reference and n_reference == _float64_count(
        stated.get("MassesCorrectedCalibration")
    ):
        info["n_reference_peaks"] = n_reference
    return info


#: The ``GlobalMetadata`` keys that state the planned raster, in index units.
IMAGING_AREA_KEYS = (
    "ImagingAreaMinXIndexPos",
    "ImagingAreaMaxXIndexPos",
    "ImagingAreaMinYIndexPos",
    "ImagingAreaMaxYIndexPos",
)


def index_value(value: Any) -> Optional[int]:
    """A ``GlobalMetadata`` raster index as an exact ``int``, or ``None``.

    The values are stored as text. ``int(float(...))`` rounded one above
    2**53 onto its neighbour, which grew a 3x2 raster to 5x2 (#520), so a
    whole number is read as one. An empty, blank or NULL value, or one that
    is not a whole number, is ``None``: the key is treated as absent.
    """
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        pass
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number) or not number.is_integer():
        return None
    return int(number)


def finite_value(value: Any) -> Optional[float]:
    """A ``GlobalMetadata`` number as a finite float, or ``None``."""
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(str(value).strip())
    except ValueError:
        return None
    return number if math.isfinite(number) else None


#: The ``GlobalMetadata`` keys that state the acquired m/z range.
MZ_RANGE_KEYS = ("MzAcqRangeLower", "MzAcqRangeUpper")


def global_metadata_values(cursor: Any, keys: Tuple[str, ...]) -> Dict[str, Any]:
    """The raw ``GlobalMetadata`` values of ``keys`` that the database holds."""
    names = ", ".join(f"'{key}'" for key in keys)
    cursor.execute(
        f"SELECT Key, Value FROM GlobalMetadata WHERE Key IN ({names})"  # nosec B608
    )
    return {key: value for key, value in cursor.fetchall()}


def readable_bounds(values: Dict[str, Any]) -> Dict[str, Any]:
    """The ImagingArea and m/z range values that can be read, parsed.

    ``ImagingArea*IndexPos`` as exact ints (:func:`index_value`), the m/z
    range as finite floats. A value that is empty, NULL or not a number is
    left out, as a key the database does not hold is: the caller falls back
    to the frame bounds, or refuses for a missing m/z range.
    """
    bounds: Dict[str, Any] = {}
    for key, value in values.items():
        number = index_value(value) if key in IMAGING_AREA_KEYS else finite_value(value)
        if number is not None:
            bounds[key] = number
    return bounds


def imaging_area_bounds(cursor: Any) -> Dict[str, int]:
    """The ``ImagingArea*IndexPos`` values a ``.d`` states, as exact ints."""
    return readable_bounds(global_metadata_values(cursor, IMAGING_AREA_KEYS))


def analysis_database(data_path: Path) -> Optional[Tuple[str, Path]]:
    """Which analysis database a timsTOF ``.d`` is read from, and its type.

    The one place that decides TSF or TDF; the reader and this extractor
    both take the answer from here, so a store cannot record one format
    while its spectra were read from the other (issue #518).
    ``analysis.tsf`` wins when both are present, unless it is empty and
    ``analysis.tdf`` is not: a 0-byte database is what a stray sqlite
    probe leaves behind, never an acquisition, and real TSF folders carry
    an empty ``analysis.tdf`` for that reason.

    Args:
        data_path: The ``.d`` directory.

    Returns:
        ``("tsf", path)`` or ``("tdf", path)``; ``None`` when neither
        database is present.
    """
    found = [
        (file_type, path)
        for file_type, path in (
            ("tsf", data_path / "analysis.tsf"),
            ("tdf", data_path / "analysis.tdf"),
        )
        if path.exists()
    ]
    if not found:
        return None

    def empty(path: Path) -> bool:
        try:
            return path.stat().st_size == 0
        except OSError:
            return False

    not_empty = [entry for entry in found if not empty(entry[1])]
    return (not_empty or found)[0]


class BrukerMetadataExtractor(MetadataExtractor):
    """Bruker-specific metadata extractor with optimized single-query extraction."""

    def __init__(
        self,
        conn: sqlite3.Connection,
        data_path: Path,
        calibration_metadata: Optional[Dict[str, Any]] = None,
        region: Optional[int] = None,
        skip_total_peaks: bool = False,
        mis_metadata: Optional[Dict[str, Any]] = None,
        spectrum_type: Optional[str] = None,
        file_type: Optional[str] = None,
    ):
        """Initialize Bruker metadata extractor.

        Args:
            conn: Active SQLite database connection
            data_path: Path to the Bruker .d directory
            calibration_metadata: Optional calibration metadata
            region: Optional region number for multi-region datasets.
                When set, coordinate bounds and frame counts are
                computed from only this region's frames.
            skip_total_peaks: If True, ``EssentialMetadata.total_peaks``
                is reported as 0 without scanning the Frames table.
                Used by metadata-only callers (preview_msi) to avoid
                a SELECT SUM(NumPeaks) FROM Frames that is the slow
                step against network-mounted .d folders.  The wizard
                step 2 preview card does not display total_peaks, so
                the missing value is safe.
            mis_metadata: The flexImaging imaging sequence (``.mis``) the
                reader found and parsed, or ``None`` when it found none.
                Handed on as ``raw_metadata["mis_metadata"]``, where the
                rapiflex and solariX extractors already put theirs, rather
                than parsed a second time here.
            spectrum_type: The representation the reader delivers, as a
                :class:`~thyra.metadata.constants.SpectrumType` name, or
                ``None`` when unknown. The database does not state it (a
                TDF is read as centroid or profile by choice of the reader),
                so the reader says.
            file_type: ``"tsf"`` or ``"tdf"``, as the reader decided it.
                ``None`` asks the same rule the reader uses
                (:func:`analysis_database`),
                so the two can never disagree about a folder that holds
                both databases (issue #518).
        """
        super().__init__(conn)
        self.conn = conn
        self.data_path = data_path
        self.calibration_metadata = calibration_metadata
        self._region = region
        self._skip_total_peaks = bool(skip_total_peaks)
        self._mis_metadata = dict(mis_metadata or {})
        # A reader that found a .mis passes its parse, empty or not, so the
        # file is read once; None leaves the lookup to this extractor.
        self._mis_given = mis_metadata is not None
        self._spectrum_type = spectrum_type
        self._file_type = file_type

    def _query_imaging_bounds(self, cursor) -> Dict[str, Any]:
        """Query imaging area bounds and the m/z range from GlobalMetadata.

        Only values that can be read are returned: the ImagingArea keys as
        exact ints (:func:`index_value`), the m/z range as finite floats. An
        empty or NULL value is left out, as a missing key is.
        """
        return readable_bounds(
            global_metadata_values(cursor, IMAGING_AREA_KEYS + MZ_RANGE_KEYS)
        )

    @staticmethod
    def _has_table(cursor, name: str) -> bool:
        """Whether ``name`` is a table of this database.

        Asked rather than caught, because the two imaging tables are
        absent for a reason a caller acts on -- the acquisition has no
        raster -- while a query that fails for any other reason should
        still surface.
        """
        cursor.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
            (name,),
        )
        return cursor.fetchone() is not None

    def _has_frame_positions(self, cursor) -> bool:
        """Whether the frames were placed on a raster.

        ``MaldiFrameInfo`` is the table that gives every frame an
        ``(XIndexPos, YIndexPos)``.  A ``.d`` written by a run that
        imaged nothing -- an electrospray acquisition on a timsOmni, say
        -- has neither it nor ``MaldiFrameLaserInfo``, and everything
        this extractor reports about the raster is then unavailable
        while everything it reports about the mass spectrometry is not.
        """
        return self._has_table(cursor, "MaldiFrameInfo")

    def _query_laser_info(self, cursor):
        """Query beam scan sizes from laser info, or ``None``.

        ``None`` both when ``MaldiFrameLaserInfo`` has no rows and when
        the acquisition has no such table at all.
        """
        if not self._has_table(cursor, "MaldiFrameLaserInfo"):
            return None
        laser_query = """
        SELECT BeamScanSizeX, BeamScanSizeY, SpotSize
        FROM MaldiFrameLaserInfo
        LIMIT 1
        """
        cursor.execute(laser_query)
        return cursor.fetchone()

    def _query_frame_info(self, cursor):
        """Query coordinate bounds and frame count from frame info.

        When a region is selected, only frames from that region are
        included in the bounds and count.

        The bounds come from the frames whose position is a whole number:
        the reader skips the others, and a text ``MAX`` ended in ``int()``
        while a REAL just below a whole number was cut down (issue #527).
        The count stays every frame, as it was.
        """
        whole = f"{whole_sql('XIndexPos')} AND {whole_sql('YIndexPos')}"
        bounds = (
            f"MIN(CASE WHEN {whole} THEN XIndexPos END), "
            f"MAX(CASE WHEN {whole} THEN XIndexPos END), "
            f"MIN(CASE WHEN {whole} THEN YIndexPos END), "
            f"MAX(CASE WHEN {whole} THEN YIndexPos END), "
            "COUNT(*) as frame_count"
        )
        if self._region is not None:
            where = " WHERE RegionNumber = ?"
            cursor.execute(
                f"SELECT {bounds} FROM MaldiFrameInfo{where}",  # nosec B608
                (self._region,),
            )
        else:
            cursor.execute(f"SELECT {bounds} FROM MaldiFrameInfo")  # nosec B608
        return cursor.fetchone()

    def _query_frame_count(self, cursor) -> int:
        """How many frames the acquisition holds.

        The count off ``Frames``, which every ``.d`` has, for the case
        where there is no ``MaldiFrameInfo`` to count rows of.
        """
        try:
            cursor.execute("SELECT COUNT(*) FROM Frames")
            result = cursor.fetchone()
            return int(result[0]) if result and result[0] else 0
        except sqlite3.OperationalError as e:
            logger.warning(f"Could not count frames: {e}")
            return 0

    def _query_total_peaks(self, cursor):
        """Query total peaks from NumPeaks column.

        When a region is selected, joins with MaldiFrameInfo to count
        only peaks from frames in that region.
        """
        try:
            if self._region is not None:
                cursor.execute(
                    "SELECT SUM(f.NumPeaks) "
                    "FROM Frames f "
                    "JOIN MaldiFrameInfo m ON f.Id = m.Frame "
                    "WHERE m.RegionNumber = ? "
                    "AND f.NumPeaks IS NOT NULL "
                    "AND f.NumPeaks > 0",
                    (self._region,),
                )
            else:
                cursor.execute(
                    "SELECT SUM(NumPeaks) FROM Frames "
                    "WHERE NumPeaks IS NOT NULL "
                    "AND NumPeaks > 0"
                )
            result = cursor.fetchone()
            if result and result[0]:
                return int(result[0])
            logger.warning("No NumPeaks data available, total_peaks will be 0")
            return 0
        except sqlite3.OperationalError as e:
            logger.warning(f"Could not query NumPeaks: {e}, " f"total_peaks will be 0")
            return 0

    def _validate_mass_range(self, bounds_data):
        """Validate that mass range data is available."""
        min_mass = bounds_data.get("MzAcqRangeLower")
        max_mass = bounds_data.get("MzAcqRangeUpper")

        if min_mass is None or max_mass is None:
            missing_keys = [
                key
                for key, value in (
                    ("MzAcqRangeLower", min_mass),
                    ("MzAcqRangeUpper", max_mass),
                )
                if value is None
            ]
            # A refusal, not a traceback (#520): the run cannot be read
            # without its m/z range, and this line says which key.
            raise ConversionRefused(
                f"{self.data_path.name} states no usable "
                f"{' or '.join(missing_keys)} in its GlobalMetadata (missing, "
                "empty or not a number), so its m/z range is unknown."
            )

        return min_mass, max_mass

    def _resolve_imaging_bounds(
        self,
        frame_result: tuple,
        bounds_data: Dict[str, Any],
    ) -> Tuple[int, int, int, int]:
        """Resolve imaging area bounds from frame data and global metadata.

        Args:
            frame_result: (min_x, max_x, min_y, max_y, count) from MaldiFrameInfo.
            bounds_data: Global metadata dictionary with ImagingArea keys.

        Returns:
            Tuple of (min_x, max_x, min_y, max_y) imaging bounds.
        """
        min_x_raw, max_x_raw, min_y_raw, max_y_raw, _ = frame_result

        if self._region is not None:
            # Per-region: use actual coordinate bounds from the filtered query
            return (
                min_x_raw or 0,
                max_x_raw or 0,
                min_y_raw or 0,
                max_y_raw or 0,
            )

        # Global: prefer ImagingArea from GlobalMetadata, fall back to frame bounds
        return (
            bounds_data.get("ImagingAreaMinXIndexPos", min_x_raw or 0),
            bounds_data.get("ImagingAreaMaxXIndexPos", max_x_raw or 0),
            bounds_data.get("ImagingAreaMinYIndexPos", min_y_raw or 0),
            bounds_data.get("ImagingAreaMaxYIndexPos", max_y_raw or 0),
        )

    def _essential_without_a_raster(
        self,
        cursor,
        mass_range: Tuple[float, float],
        total_peaks: int,
    ) -> EssentialMetadata:
        """The record for a ``.d`` whose frames were never placed.

        Everything the mass spectrometry says is here; everything the
        raster would have said is the zero that means "none".  The
        absent ``pixel_size`` is what :func:`thyra.convert.convert_msi`
        refuses on, and :class:`~thyra.readers.bruker.timstof.BrukerReader`
        refuses earlier and more specifically for a conversion, so the
        record only ever reaches a metadata-only caller.
        """
        return EssentialMetadata(
            dimensions=(0, 0, 1),
            coordinate_bounds=(0.0, 0.0, 0.0, 0.0),
            mass_range=mass_range,
            pixel_size=None,
            n_spectra=self._query_frame_count(cursor),
            total_peaks=total_peaks,
            source_path=str(self.data_path),
            spectrum_type=self._spectrum_type,
        )

    def _extract_essential_impl(self) -> EssentialMetadata:
        """Extract essential metadata with proper coordinate normalization.

        ``GlobalMetadata`` and ``Frames`` are read first because every
        Bruker ``.d`` has them.  The two MALDI tables hold the raster and
        are consulted only once they are known to exist: asking for
        ``MaldiFrameLaserInfo`` up front made the whole metadata block
        unreachable for an acquisition that imaged nothing, which failed
        with ``no such table: MaldiFrameLaserInfo`` before even the mass
        range had been read.
        """
        cursor = self.conn.cursor()

        # Query database tables
        bounds_data = self._query_imaging_bounds(cursor)
        if self._skip_total_peaks:
            total_peaks = 0
        else:
            total_peaks = self._query_total_peaks(cursor)

        try:
            min_mass, max_mass = self._validate_mass_range(bounds_data)
            mass_range = (float(min_mass), float(max_mass))

            if not self._has_frame_positions(cursor):
                logger.info(
                    "%s has no MaldiFrameInfo table: the acquisition placed "
                    "no frames on a raster, so the metadata describes no "
                    "image.",
                    self.data_path.name,
                )
                return self._essential_without_a_raster(cursor, mass_range, total_peaks)

            laser_result = self._query_laser_info(cursor)
            frame_result = self._query_frame_info(cursor)
            if not frame_result:
                raise ValueError("No data found in MaldiFrameInfo table")

            imaging_min_x, imaging_max_x, imaging_min_y, imaging_max_y = (
                self._resolve_imaging_bounds(frame_result, bounds_data)
            )

            # Store imaging area offsets for coordinate normalization
            imaging_area_offsets = (int(imaging_min_x), int(imaging_min_y), 0)

            # Normalize coordinates to start from 0
            max_x = float(imaging_max_x - imaging_min_x)
            max_y = float(imaging_max_y - imaging_min_y)

            # Extract beam sizes
            beam_x, beam_y, spot_size = (
                laser_result if laser_result else (None, None, None)
            )

            # Build final metadata objects
            dimensions = self._calculate_dimensions_from_coords(0.0, max_x, 0.0, max_y)
            coordinate_bounds = (0.0, float(max_x), 0.0, float(max_y))
            pixel_size = self._resolve_pixel_size_um(beam_x, beam_y)
            _, _, _, _, frame_count = frame_result
            n_spectra = int(frame_count) if frame_count else 0

            return EssentialMetadata(
                dimensions=dimensions,
                coordinate_bounds=coordinate_bounds,
                mass_range=mass_range,
                pixel_size=pixel_size,
                n_spectra=n_spectra,
                total_peaks=total_peaks,
                source_path=str(self.data_path),
                coordinate_offsets=imaging_area_offsets,
                spectrum_type=self._spectrum_type,
            )

        except sqlite3.OperationalError as e:
            logger.error(f"SQL error extracting essential metadata: {e}")
            raise ValueError(
                f"Failed to extract essential metadata from Bruker database: " f"{e}"
            )
        except ConversionRefused:
            # Nothing unexpected about a refusal, and it is already the
            # whole explanation. _resolve_pixel_size_um below reads the
            # sibling .mis, which refuses an entity-bearing document, so
            # one reaches here; the clause under this one would log it at
            # ERROR as "Unexpected error ..." and then re-raise it into
            # convert_msi, which logs it again. Two lines, the first of
            # them untrue.
            raise
        except Exception as e:
            logger.error(f"Unexpected error extracting essential metadata: {e}")
            raise

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        """Extract comprehensive metadata with additional database queries."""
        essential = self.get_essential()

        raw_metadata = self._extract_global_metadata()
        if self._mis_metadata:
            # The optical image, teaching points and areas of the imaging
            # sequence: the only record of how the raster was registered
            # onto the optical image, which the reader has always used and
            # no tsf/tdf store has carried.
            raw_metadata["mis_metadata"] = dict(self._mis_metadata)

        return ComprehensiveMetadata(
            essential=essential,
            format_specific=self._extract_bruker_specific(),
            acquisition_params=self._extract_acquisition_params(),
            instrument_info=self._extract_instrument_info(),
            raw_metadata=raw_metadata,
        )

    def _sequence_metadata(self) -> Optional[Dict[str, Any]]:
        """The parsed ``.mis`` of this ``.d``, or ``None`` when it has none.

        The one the reader adopted, when it passed it in, so the pitch and
        the alignment read the same file and its warnings are logged once;
        else the shared locator's pick (issue #303).

        Raises:
            ConversionRefused: If the ``.mis`` found is a document defusedxml
                refuses.
        """
        if self._mis_given:
            return self._mis_metadata
        from ...readers.bruker.mis_parser import (
            find_mis_file_for_d_folder,
            parse_mis_file,
        )

        mis_path = find_mis_file_for_d_folder(self.data_path)
        return parse_mis_file(mis_path) if mis_path is not None else None

    def _beam_scan_size(
        self, beam_x: object, beam_y: object
    ) -> Optional[Tuple[float, float]]:
        """The beam scan size as two finite numbers above zero, or ``None``.

        sqlite stores whatever was written, so a column declared REAL can
        hold ``'20 um'``. A value that is not a number is treated like a
        missing one, with a warning naming it, instead of ending the
        conversion in ``float()`` even when ``--pixel-size`` was given.
        """
        if not beam_x or not beam_y:
            return None
        try:
            beam = (float(beam_x), float(beam_y))  # type: ignore[arg-type]
        except (TypeError, ValueError):
            beam = (math.nan, math.nan)
        if all(math.isfinite(v) and v > 0 for v in beam):
            return beam
        logger.warning(
            "MaldiFrameLaserInfo BeamScanSize (%r, %r) in %s is not two numbers "
            "above zero, so it is treated as missing",
            beam_x,
            beam_y,
            self.data_path.name,
        )
        return None

    def _resolve_pixel_size_um(
        self,
        beam_x: Optional[float],
        beam_y: Optional[float],
    ) -> Optional[Tuple[float, float]]:
        """Resolve pixel pitch in micrometers, preferring .mis <Raster>.

        The Bruker SDK exposes ``MaldiFrameLaserInfo.BeamScanSizeX/Y`` which is
        the area each laser shot scans, not the raster step. When the
        acquisition oversamples (BeamScanSize > Raster), these values disagree
        and the canonical pixel pitch is the Raster step from the FlexImaging
        ``.mis`` file. Prefer the Raster step; warn when it disagrees with
        BeamScanSize. Fall back to BeamScanSize when no .mis is found.

        Raises:
            ConversionRefused: If a ``.mis`` is found and is a document
                defusedxml refuses. No .mis at all is the fall-back case
                above; one that is present and refused is not, so it is
                not quietly resolved to BeamScanSize.
        """
        raster: Optional[Tuple[float, float]] = None
        mis_data = self._sequence_metadata()
        if mis_data is not None:
            raster_xy = mis_data.get("raster")
            if (
                isinstance(raster_xy, list)
                and len(raster_xy) == 2
                and all(v > 0 for v in raster_xy)
            ):
                raster = (float(raster_xy[0]), float(raster_xy[1]))

        beam = self._beam_scan_size(beam_x, beam_y)

        if raster is not None and beam is not None:
            if raster != beam:
                logger.warning(
                    "Pixel size mismatch: .mis Raster=(%g, %g) um, "
                    "BeamScanSize=(%g, %g) um. %s "
                    "Using Raster step. Override with --pixel-size if needed.",
                    raster[0],
                    raster[1],
                    beam[0],
                    beam[1],
                    _mismatch_reading(raster, beam),
                )
            return raster

        if mis_data is not None and raster is None and beam is not None:
            # The .mis is there and its step is not (#520): say what the
            # pitch is taken from instead, rather than switching silently.
            logger.warning(
                "The .mis states no usable Raster step, so the pixel size is "
                "BeamScanSize (%g, %g) um. Pass --pixel-size if that is wrong.",
                beam[0],
                beam[1],
            )
        if raster is not None:
            return raster
        return beam

    def _calculate_dimensions_from_coords(
        self,
        min_x: Optional[float],
        max_x: Optional[float],
        min_y: Optional[float],
        max_y: Optional[float],
    ) -> Tuple[int, int, int]:
        """Calculate dataset dimensions from coordinate bounds."""
        if min_x is None or max_x is None or min_y is None or max_y is None:
            return (0, 0, 1)  # Default for problematic data

        # Bruker coordinates are typically in position units
        # Calculate grid dimensions assuming integer grid positions
        x_range = int(max_x - min_x) + 1 if max_x > min_x else 1
        y_range = int(max_y - min_y) + 1 if max_y > min_y else 1

        return (max(1, x_range), max(1, y_range), 1)  # Assume 2D data (z=1)

    def _extract_bruker_specific(self) -> Dict[str, Any]:
        """Extract Bruker format-specific metadata."""
        is_tdf = self._is_tdf_format()
        stem = "analysis.tdf" if is_tdf else "analysis.tsf"
        format_specific: Dict[str, Any] = {
            "bruker_format": ("bruker_tdf" if is_tdf else "bruker_tsf"),
            "data_format": ("bruker_tdf" if is_tdf else "bruker_tsf"),
            "data_path": str(self.data_path),
            "database_path": str(self.data_path / stem),
            "binary_file": str(self.data_path / f"{stem}_bin"),
            "is_maldi": self._is_maldi_dataset(),
            "ion_mobility": self._extract_ion_mobility(),
        }

        source_type = self._instrument_source_type()
        if source_type is not None:
            format_specific["instrument_source_type"] = source_type

        # Add calibration metadata if available
        if self.calibration_metadata:
            format_specific["calibration"] = self.calibration_metadata

        instrument_calibration = self._extract_instrument_calibration()
        if instrument_calibration:
            format_specific["instrument_calibration"] = instrument_calibration

        return format_specific

    def _extract_instrument_calibration(self) -> Dict[str, Any]:
        """The m/z calibration the acquisition ran under, as the database states it.

        Every tsf and tdf acquisition read so far -- five, from timsControl
        4.1 to 7.2, MALDI and electrospray -- carries a ``CalibrationInfo``
        table in its own analysis database: the external calibration the
        instrument had when the run started, with its time, the software
        and version that made it, and its fit. ``calibration.sqlite`` is a
        different thing; on a MALDI run its first state is the online
        lock-mass calibration, written as the acquisition starts, and an
        electrospray run has no such file at all.

        Only the keys the calibration section needs are asked for. The
        table also names who calibrated (``CalibrationUser``,
        ``MobilityCalibrationUser``) and the reference list the lab chose,
        whose name can be a person's; none of those is read (see
        :mod:`thyra.metadata.personal_data`). Neither are the mobility
        calibration's keys: on the one imaging TDF read,
        ``MeasuredTimsVoltages`` is all zeros and
        ``MobilityStandardDeviationPercent`` is 3578, so there they do not
        describe a mobility calibration, and no second imaging TDF was at
        hand to say what they describe elsewhere.

        ``n_reference_peaks`` counts the reference masses that carry a
        corrected mass, i.e. the peaks the calibration was fitted to.
        ``MzStandardDeviationPPM`` is the square root of their summed squared
        ppm errors over ``n - 1``, recomputed from those arrays to six
        decimals on four acquisitions; it is handed on as stated, zero
        included, and the schema builder decides what a zero means.

        Returns:
            The stated facts under Thyra's key names, or an empty dict when
            the table is absent, unreadable, or holds rows for a polarity
            this acquisition did not run in.
        """
        cursor = self.conn.cursor()
        try:
            if not self._has_table(cursor, "CalibrationInfo"):
                return {}
            polarity = self._calibration_polarity(cursor)
            if polarity is None:
                return {}
            cursor.execute(
                "SELECT KeyName, Value FROM CalibrationInfo WHERE KeyPolarity = ? "
                "AND KeyName IN (?, ?, ?, ?, ?, ?)",
                (polarity, *_CALIBRATION_INFO_KEYS),
            )
            stated = dict(cursor.fetchall())
        except sqlite3.OperationalError as e:
            logger.debug("Could not read CalibrationInfo: %s", e)
            return {}
        return _instrument_calibration(stated)

    def _calibration_polarity(self, cursor) -> Optional[str]:
        """Which polarity's ``CalibrationInfo`` rows describe this acquisition.

        The table keys every row by polarity. When the frames all ran in one
        polarity, that polarity's rows are the calibration, and a table with
        none for it says nothing about this acquisition. When the frames
        alternate, no one calibration describes them all, so none is taken.
        When they do not say at all -- no column, or no rows -- a table that
        holds one polarity is taken, and one that holds two is not guessed
        between.
        """
        cursor.execute("SELECT DISTINCT KeyPolarity FROM CalibrationInfo")
        stated = {str(row[0]) for row in cursor.fetchall() if row[0] is not None}
        frames = self._distinct_frame_polarities(cursor)
        if frames:
            polarity = str(frames[0])
            return polarity if len(frames) == 1 and polarity in stated else None
        return next(iter(stated)) if len(stated) == 1 else None

    def _instrument_source_type(self) -> Optional[int]:
        """``GlobalMetadata.InstrumentSourceType``, as the file states it.

        Bruker's own code for what ionised the sample.  Recorded raw and
        deliberately not mapped onto an ionisation term: the MALDI
        imaging acquisitions available carry ``1`` and an electrospray
        timsOmni run carries ``11``, which is two points on an enum
        whose labels Bruker does not document in the file.  Guessing the
        rest would put a CV accession in the store on the strength of
        two samples, where ``ms_analysis.ionisation_source`` is instead
        left unset -- which is what an unset field is for.  A MALDI run
        still fills that field, from the laser tables, as before.
        """
        cursor = self.conn.cursor()
        try:
            cursor.execute(
                "SELECT Value FROM GlobalMetadata WHERE Key = " "'InstrumentSourceType'"
            )
            result = cursor.fetchone()
        except sqlite3.OperationalError:
            return None
        if not result or result[0] is None:
            return None
        try:
            return int(str(result[0]).strip())
        except (TypeError, ValueError):
            logger.debug("InstrumentSourceType is not an integer: %r", result[0])
            return None

    def _extract_ion_mobility(self) -> Dict[str, Any]:
        """Describe the mobility dimension of the acquisition.

        A TDF file is a TIMS acquisition: every frame carries ``NumScans``
        mobility scans, each mapping onto an inverse reduced ion mobility
        (1/K0, PSI-MS ``MS:1002815``) through the file's calibration, and
        ``GlobalMetadata`` records the acquired 1/K0 range. A TSF file has
        no mobility dimension at all. The block says which, so a consumer
        can tell "summed over mobility" from "never had any".
        """
        if not self._is_tdf_format():
            return {"present": False}

        info: Dict[str, Any] = {
            "present": True,
            "separation": "inverse reduced ion mobility",
            "separation_accession": "MS:1002815",
            "unit": "volt-second per square centimeter",
            "unit_accession": "MS:1002814",
        }
        cursor = self.conn.cursor()
        try:
            cursor.execute("SELECT MIN(NumScans), MAX(NumScans) FROM Frames")
            row = cursor.fetchone()
            if row and row[0] is not None:
                info["num_scans_min"] = int(row[0])
                info["num_scans_max"] = int(row[1])
            cursor.execute(
                "SELECT Key, Value FROM GlobalMetadata WHERE Key IN "
                "('OneOverK0AcqRangeLower', 'OneOverK0AcqRangeUpper')"
            )
            bounds = {key: float(value) for key, value in cursor.fetchall()}
            lower = bounds.get("OneOverK0AcqRangeLower")
            upper = bounds.get("OneOverK0AcqRangeUpper")
            if lower is not None and upper is not None:
                info["one_over_k0_range"] = [lower, upper]
        except (sqlite3.OperationalError, TypeError, ValueError) as e:
            logger.debug(f"Could not read the mobility range: {e}")
        return info

    def _extract_acquisition_params(self) -> Dict[str, Any]:
        """Extract acquisition parameters from database."""
        params: Dict[str, Any] = {}
        cursor = self.conn.cursor()

        # Extract laser parameters if available
        self._extract_laser_params(cursor, params)

        # Extract timing parameters
        self._extract_timing_params(cursor, params)

        # Polarity, which is a per-frame column rather than a global key
        self._extract_polarity(cursor, params)

        return params

    def _extract_polarity(self, cursor, params: Dict[str, Any]) -> None:
        """Record the ion polarity when every frame agrees on one.

        ``Frames.Polarity`` is a single character per frame, ``'+'`` or
        ``'-'``, which :func:`~thyra.metadata.schema.vocab.normalize_polarity`
        already understands.  It is the only place a Bruker ``.d`` states
        the polarity, so before this the field was simply never filled
        for any Bruker source.

        A file whose frames disagree is left unset rather than reduced to
        one of its two values: an alternating-polarity acquisition has no
        single truthful answer, and the same rule is what
        :func:`~thyra.metadata.schema.builder._polarity_from_cv_params`
        applies to an imzML that declares both.  Only the selected
        region's frames are asked about when a region is selected.
        """
        polarities = self._distinct_frame_polarities(cursor)
        if polarities is None:
            logger.debug("Frames has no Polarity column")
            return
        if len(polarities) != 1:
            if len(polarities) > 1:
                logger.info(
                    "Frames record more than one polarity; the acquisition "
                    "alternated and no single value is recorded."
                )
            return
        value = polarities[0]
        if isinstance(value, str) and value.strip():
            params["polarity"] = value.strip()

    def _distinct_frame_polarities(self, cursor) -> Optional[List[Any]]:
        """Up to two distinct ``Frames.Polarity`` values, or ``None``.

        Two are enough to know the acquisition alternated. ``None`` means
        the question could not be asked -- no such column -- which is not
        the same answer as an empty list. Only the selected region's frames
        are asked about when a region is selected.
        """
        try:
            if self._region is not None:
                cursor.execute(
                    "SELECT DISTINCT f.Polarity FROM Frames f "
                    "JOIN MaldiFrameInfo m ON f.Id = m.Frame "
                    "WHERE m.RegionNumber = ? AND f.Polarity IS NOT NULL "
                    "LIMIT 2",
                    (self._region,),
                )
            else:
                cursor.execute(
                    "SELECT DISTINCT Polarity FROM Frames "
                    "WHERE Polarity IS NOT NULL LIMIT 2"
                )
            return [row[0] for row in cursor.fetchall()]
        except sqlite3.OperationalError:
            return None

    def _extract_laser_params(self, cursor, params: Dict[str, Any]) -> None:
        """Extract laser parameters: per-frame settings, then beam geometry.

        ``LaserPower``, ``NumLaserShots`` and ``LaserRepRate`` are columns
        of ``MaldiFrameInfo``, one row per frame. ``LaserPower`` is the
        percentage of the laser's range that timsControl shows as
        "Laser Power" (70.0 on the acquisition checked), ``LaserRepRate``
        is in Hz. The three are constant within a normal acquisition, so
        each is reported as a scalar when every frame agrees and as a
        ``[min, max]`` pair when they differ -- the way the solariX
        extractor reports the same three -- and only the frames of the
        selected region are aggregated when there is one. The beam scan
        sizes and spot size live in ``MaldiFrameLaserInfo``.

        An earlier version selected ``LaserPower`` and ``LaserFrequency``
        from ``MaldiFrameLaserInfo``, which has never had those columns
        (TDF/TSF schema 3.x), so no real acquisition ever filled any of
        these keys before.
        """
        try:
            if self._region is not None:
                cursor.execute(
                    "SELECT MIN(LaserPower), MAX(LaserPower), "
                    "MIN(NumLaserShots), MAX(NumLaserShots), "
                    "MIN(LaserRepRate), MAX(LaserRepRate) "
                    "FROM MaldiFrameInfo WHERE RegionNumber = ?",
                    (self._region,),
                )
            else:
                cursor.execute(
                    "SELECT MIN(LaserPower), MAX(LaserPower), "
                    "MIN(NumLaserShots), MAX(NumLaserShots), "
                    "MIN(LaserRepRate), MAX(LaserRepRate) "
                    "FROM MaldiFrameInfo"
                )
            result = cursor.fetchone()
            if result:
                self._process_frame_laser_result(result, params)
        except sqlite3.OperationalError:
            logger.debug("Could not extract per-frame laser parameters")

        try:
            result = self._query_laser_info(cursor)
            if result:
                self._process_beam_result(result, params)
        except sqlite3.OperationalError:
            logger.debug("Could not extract beam scan parameters")

    @staticmethod
    def _process_frame_laser_result(result, params: Dict[str, Any]) -> None:
        """Store the per-frame laser aggregates, collapsed when constant."""
        pairs = (
            ("laser_power", result[0], result[1]),
            ("num_laser_shots", result[2], result[3]),
            ("laser_frequency", result[4], result[5]),
        )
        for key, low, high in pairs:
            if low is None or high is None:
                continue
            params[key] = low if low == high else [low, high]

    @staticmethod
    def _process_beam_result(result, params: Dict[str, Any]) -> None:
        """Store the beam scan sizes and spot size."""
        beam_x, beam_y, spot_size = result
        if beam_x is not None:
            params["beam_scan_size_x"] = beam_x
            params["BeamScanSizeX"] = beam_x  # Add both formats for compatibility
        if beam_y is not None:
            params["beam_scan_size_y"] = beam_y
            params["BeamScanSizeY"] = beam_y  # Add both formats for compatibility
        if spot_size is not None:
            params["laser_spot_size"] = spot_size

    def _extract_timing_params(self, cursor, params: Dict[str, Any]) -> None:
        """Extract the acquisition timestamp and method name from GlobalMetadata.

        ``MethodName`` is recorded as timsControl opened it, which may be
        an absolute path on the acquisition PC; only the file name (the
        ``*.m`` directory's name) is kept, since a path names a machine
        the store will not be opened on. ``raw_metadata["global_metadata"]``
        keeps the same name and not the path: every vendor dictionary is
        reduced the same way (see :mod:`thyra.metadata.personal_data`).
        """
        try:
            cursor.execute(
                "SELECT Value FROM GlobalMetadata WHERE Key = " "'AcquisitionDateTime'"
            )
            result = cursor.fetchone()
            if result:
                params["acquisition_datetime"] = result[0]
        except sqlite3.OperationalError:
            pass
        try:
            cursor.execute("SELECT Value FROM GlobalMetadata WHERE Key = 'MethodName'")
            result = cursor.fetchone()
            if result and isinstance(result[0], str):
                name = result[0].replace("\\", "/").rsplit("/", 1)[-1].strip()
                if name:
                    params["method_name"] = name
        except sqlite3.OperationalError:
            pass

    def _extract_instrument_info(self) -> Dict[str, Any]:
        """Extract instrument information from global metadata."""
        instrument = {}
        cursor = self.conn.cursor()

        # Common instrument metadata keys
        # ``InstrumentVendor`` is a GlobalMetadata key on every tsf and
        # tdf acquisition and was the one instrument fact this query did
        # not ask for, so who built a timsTOF-family instrument reached
        # no store while solariX -- reading the same key out of its own
        # Properties table -- recorded it. Mapped to ``manufacturer``,
        # the spelling the other extractors already use.
        #
        # ``MzCalibrationMode`` used to be asked for here and never found:
        # it is not a GlobalMetadata key on any acquisition read (five,
        # tsf and tdf) but a ``CalibrationInfo`` one, and there an
        # undocumented code (1 to 4 across those five). No store ever
        # carried it, and it is not read anywhere now.
        instrument_keys = [
            ("InstrumentName", "instrument_name"),
            ("InstrumentSerialNumber", "instrument_serial_number"),
            ("InstrumentModel", "instrument_model"),
            ("InstrumentVendor", "manufacturer"),
            ("SoftwareVersion", "software_version"),
        ]

        try:
            for db_key, result_key in instrument_keys:
                cursor.execute(
                    "SELECT Value FROM GlobalMetadata WHERE Key = ?", (db_key,)
                )
                result = cursor.fetchone()
                if result:
                    instrument[result_key] = result[0]
        except sqlite3.OperationalError:
            logger.debug("Could not extract instrument metadata")

        return instrument

    def _extract_global_metadata(self) -> Dict[str, Any]:
        """Extract all global metadata from database.

        The two reads are separate because they can fail separately: a
        ``.d`` with no raster has no ``MaldiFrameLaserInfo`` and a
        perfectly good ``GlobalMetadata``, and one ``try`` around both
        would have dropped the per-frame block silently while keeping
        the global one only because it was assigned first.
        """
        raw_metadata: Dict[str, Any] = {}
        cursor = self.conn.cursor()

        try:
            cursor.execute("SELECT Key, Value FROM GlobalMetadata")
            global_metadata = {}
            for key, value in cursor.fetchall():
                global_metadata[key] = value
            raw_metadata["global_metadata"] = global_metadata
        except sqlite3.OperationalError:
            logger.debug("GlobalMetadata table not found or accessible")

        if not self._has_table(cursor, "MaldiFrameLaserInfo"):
            return raw_metadata

        try:
            # Extract frame info for tests that expect it
            cursor.execute(
                "SELECT Id, SpotXPos, SpotYPos, BeamScanSizeX, BeamScanSizeY "
                "FROM MaldiFrameLaserInfo"
            )
            frame_info = []
            for row in cursor.fetchall():
                frame_info.append(
                    {
                        "id": row[0],
                        "x_pos": row[1],
                        "y_pos": row[2],
                        "beam_x": row[3],
                        "beam_y": row[4],
                    }
                )
            raw_metadata["frame_info"] = frame_info
        except sqlite3.OperationalError:
            logger.debug("MaldiFrameLaserInfo is not readable")

        return raw_metadata

    def _is_tdf_format(self) -> bool:
        """Whether the reader reads ``analysis.tdf`` (TDF) rather than ``analysis.tsf``."""
        if self._file_type is None:
            found = analysis_database(self.data_path)
            self._file_type = "tsf" if found is None else found[0]
        return self._file_type == "tdf"

    def _is_maldi_dataset(self) -> bool:
        """Check if this is a MALDI dataset by checking for laser info."""
        cursor = self.conn.cursor()
        try:
            cursor.execute("SELECT COUNT(*) FROM MaldiFrameLaserInfo")
            result = cursor.fetchone()
            return bool(result and result[0] > 0)
        except sqlite3.OperationalError:
            return False
