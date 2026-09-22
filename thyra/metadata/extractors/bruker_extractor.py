# thyra/metadata/extractors/bruker_extractor.py
import logging
import sqlite3
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from ...core.base_extractor import MetadataExtractor
from ...errors import ConversionRefused
from ..types import ComprehensiveMetadata, EssentialMetadata

logger = logging.getLogger(__name__)


class BrukerMetadataExtractor(MetadataExtractor):
    """Bruker-specific metadata extractor with optimized single-query extraction."""

    def __init__(
        self,
        conn: sqlite3.Connection,
        data_path: Path,
        calibration_metadata: Optional[Dict[str, Any]] = None,
        region: Optional[int] = None,
        skip_total_peaks: bool = False,
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
        """
        super().__init__(conn)
        self.conn = conn
        self.data_path = data_path
        self.calibration_metadata = calibration_metadata
        self._region = region
        self._skip_total_peaks = bool(skip_total_peaks)

    def _query_imaging_bounds(self, cursor):
        """Query imaging area bounds from GlobalMetadata."""
        imaging_bounds_query = """
        SELECT Key, Value FROM GlobalMetadata
        WHERE Key IN ('ImagingAreaMinXIndexPos', 'ImagingAreaMaxXIndexPos',
                      'ImagingAreaMinYIndexPos', 'ImagingAreaMaxYIndexPos',
                      'MzAcqRangeLower', 'MzAcqRangeUpper')
        """
        cursor.execute(imaging_bounds_query)
        return {row[0]: float(row[1]) for row in cursor.fetchall()}

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
        """
        if self._region is not None:
            frame_query = """
            SELECT
                MIN(XIndexPos), MAX(XIndexPos),
                MIN(YIndexPos), MAX(YIndexPos),
                COUNT(*) as frame_count
            FROM MaldiFrameInfo
            WHERE RegionNumber = ?
            """
            cursor.execute(frame_query, (self._region,))
        else:
            frame_query = """
            SELECT
                MIN(XIndexPos), MAX(XIndexPos),
                MIN(YIndexPos), MAX(YIndexPos),
                COUNT(*) as frame_count
            FROM MaldiFrameInfo
            """
            cursor.execute(frame_query)
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

        missing_keys = []
        if min_mass is None:
            missing_keys.append("MzAcqRangeLower")
        if max_mass is None:
            missing_keys.append("MzAcqRangeUpper")

        if missing_keys:
            error_msg = (
                f"Missing critical mass range bounds in GlobalMetadata: "
                f"{', '.join(missing_keys)}. Cannot establish mass range."
            )
            logger.error(error_msg)
            raise ValueError(error_msg)

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

        return ComprehensiveMetadata(
            essential=essential,
            format_specific=self._extract_bruker_specific(),
            acquisition_params=self._extract_acquisition_params(),
            instrument_info=self._extract_instrument_info(),
            raw_metadata=self._extract_global_metadata(),
        )

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
        from ...readers.bruker.mis_parser import (
            find_mis_file_for_d_folder,
            parse_mis_file,
        )

        raster: Optional[Tuple[float, float]] = None
        mis_path = find_mis_file_for_d_folder(self.data_path)
        if mis_path is not None:
            mis_data = parse_mis_file(mis_path)
            raster_xy = mis_data.get("raster")
            if (
                isinstance(raster_xy, list)
                and len(raster_xy) == 2
                and all(v > 0 for v in raster_xy)
            ):
                raster = (float(raster_xy[0]), float(raster_xy[1]))

        beam: Optional[Tuple[float, float]] = (
            (float(beam_x), float(beam_y)) if beam_x and beam_y else None
        )

        if raster is not None and beam is not None:
            if raster != beam:
                logger.warning(
                    "Pixel size mismatch: .mis Raster=(%g, %g) um, "
                    "BeamScanSize=(%g, %g) um. Acquisition is oversampled. "
                    "Using Raster step. Override with --pixel-size if needed.",
                    raster[0],
                    raster[1],
                    beam[0],
                    beam[1],
                )
            return raster

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

        return format_specific

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
            rows = cursor.fetchall()
        except sqlite3.OperationalError:
            logger.debug("Frames has no Polarity column")
            return
        if len(rows) != 1:
            if len(rows) > 1:
                logger.info(
                    "Frames record more than one polarity; the acquisition "
                    "alternated and no single value is recorded."
                )
            return
        value = rows[0][0]
        if isinstance(value, str) and value.strip():
            params["polarity"] = value.strip()

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
        the store will not be opened on. The full value is still in
        ``raw_metadata["global_metadata"]``.
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
        instrument_keys = [
            ("InstrumentName", "instrument_name"),
            ("InstrumentSerialNumber", "instrument_serial_number"),
            ("InstrumentModel", "instrument_model"),
            ("SoftwareVersion", "software_version"),
            ("MzCalibrationMode", "mz_calibration_mode"),
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
        """Check if this is TDF format (vs TSF)."""
        return (self.data_path / "analysis.tdf").exists()

    def _is_maldi_dataset(self) -> bool:
        """Check if this is a MALDI dataset by checking for laser info."""
        cursor = self.conn.cursor()
        try:
            cursor.execute("SELECT COUNT(*) FROM MaldiFrameLaserInfo")
            result = cursor.fetchone()
            return bool(result and result[0] > 0)
        except sqlite3.OperationalError:
            return False
