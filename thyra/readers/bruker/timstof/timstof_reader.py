"""Bruker reader implementation combining best features from all implementations.

This module provides a high-performance, memory-efficient reader for
Bruker TSF/TDF data formats with lazy loading, intelligent caching, and
comprehensive error handling.
"""

import logging
import os
import sqlite3
from pathlib import Path
from typing import Any, Callable, Dict, Generator, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

# Set OpenMP thread limit before any SDK imports to control Bruker DLL
# threading
if "OMP_NUM_THREADS" not in os.environ:
    os.environ["OMP_NUM_THREADS"] = "4"
    logging.getLogger(__name__).info(
        "Set OMP_NUM_THREADS=4 to limit Bruker DLL threading"
    )
else:
    current_setting = os.environ.get("OMP_NUM_THREADS")
    logging.getLogger(__name__).info(
        f"Using existing OMP_NUM_THREADS={current_setting}"
    )

from ....core.base_extractor import MetadataExtractor
from ....core.registry import register_reader
from ....metadata.extractors.bruker_extractor import BrukerMetadataExtractor
from ....utils.bruker_exceptions import DataError, FileFormatError, SDKError
from ..base_bruker_reader import BrukerBaseMSIReader
from ..folder_structure import BrukerFolderStructure, BrukerFormat
from ..mis_parser import parse_mis_file
from .sdk.dll_manager import DLLManager
from .sdk.sdk_functions import SDKFunctions

logger = logging.getLogger(__name__)


def build_raw_mass_axis(
    spectra_iterator: Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ],
    progress_callback: Optional[Callable[[int], None]] = None,
) -> Tuple[NDArray[np.float64], int]:
    """Build raw mass axis from spectra iterator.

    Raw Mass axis in case the user wants the full data. Not recommended for
    normal use.
    Future interpolation module will create optimized mass axis using
    min/max mass + bin width.

    Also counts total peaks for COO matrix pre-allocation.

    Args:
        spectra_iterator: Iterator yielding (coords, mzs, intensities) tuples
        progress_callback: Optional callback for progress updates

    Returns:
        Tuple of (numpy array of unique m/z values in ascending order, total_peaks)
    """
    from tqdm import tqdm

    unique_mzs: set[float] = set()
    count = 0
    total_peaks = 0

    # Create progress bar
    pbar = tqdm(
        desc="Building raw mass axis and counting peaks",
        unit=" spectra",
        dynamic_ncols=True,
        bar_format="{l_bar}{bar}| {n_fmt} spectra [{elapsed}<{remaining}]",
    )

    try:
        for coords, mzs, intensities in spectra_iterator:
            if mzs.size > 0:
                unique_mzs.update(mzs)
                total_peaks += len(mzs)
            count += 1
            pbar.update(1)

            # Log progress periodically
            if count % 10000 == 0:
                pbar.set_postfix(
                    {
                        "unique_mz": len(unique_mzs),
                        "total_peaks": total_peaks,
                        "memory_est_mb": len(unique_mzs) * 8 / 1024 / 1024,
                    }
                )

            if progress_callback and count % 100 == 0:
                progress_callback(count)
    finally:
        pbar.close()

    logger.info(f"Total peaks counted: {total_peaks:,}")
    return (np.array(sorted(unique_mzs)), total_peaks)


def _get_frame_coordinates(
    db_path: Path,
    frame_id: int,
    coordinate_offsets: Optional[Tuple[int, int, int]] = None,
) -> Optional[Tuple[int, int, int]]:
    """Get normalized coordinates for a specific frame directly from database.

    Args:
        db_path: Path to the SQLite database file
        frame_id: Frame ID to look up
        coordinate_offsets: Optional coordinate offsets for normalization
                           (x_offset, y_offset, z_offset)

    Returns:
        Tuple of normalized (x, y, z) coordinates (0-based), or None if not
        found
    """
    try:
        with sqlite3.connect(str(db_path)) as conn:
            cursor = conn.cursor()

            # Check if this is MALDI data
            try:
                cursor.execute(
                    "SELECT XIndexPos, YIndexPos FROM MaldiFrameInfo WHERE "
                    "Frame = ?",
                    (frame_id,),
                )
                result = cursor.fetchone()
                if result:
                    x, y = result
                    # Apply coordinate offsets if provided (Bruker-specific
                    # normalization)
                    if coordinate_offsets:
                        offset_x, offset_y, offset_z = coordinate_offsets
                        return (int(x) - offset_x, int(y) - offset_y, 0)
                    else:
                        return (int(x), int(y), 0)
            except sqlite3.OperationalError:
                # No MALDI table, use generated coordinates
                pass

            # For non-MALDI data, generate coordinates (simple sequential
            # mapping)
            return (frame_id - 1, 0, 0)

    except Exception as e:
        logger.warning(f"Error getting coordinates for frame {frame_id}: {e}")
        return None


def _get_frame_count(db_path: Path) -> int:
    """Get total frame count directly from database.

    Args:
        db_path: Path to the SQLite database file

    Returns:
        Total number of frames
    """
    try:
        with sqlite3.connect(str(db_path)) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM Frames")
            return int(cursor.fetchone()[0])
    except Exception as e:
        logger.error(f"Error getting frame count: {e}")
        return 0


@register_reader("bruker")
class BrukerReader(BrukerBaseMSIReader):
    """Bruker reader for TSF/TDF data formats.

    Features:
    - Sequential spectrum iteration
    - Direct database coordinate access
    - Robust SDK integration with fallback mechanisms
    - Comprehensive error handling and recovery
    - Compatible with spatialdata_converter.py interface
    """

    def __init__(
        self,
        data_path: Path,
        use_recalibrated_state: bool = True,
        cache_coordinates: bool = True,
        memory_limit_gb: Optional[float] = None,
        batch_size: Optional[int] = None,
        progress_callback: Optional[Callable[[int, int], None]] = None,
        region: Optional[int] = None,
        **kwargs,
    ):
        """Initialize the Bruker reader.

        Args:
            data_path: Path to Bruker .d directory, or parent folder containing it.
                If a parent folder is provided, the reader will automatically find
                the .d folder within it.
            use_recalibrated_state: Whether to use recalibrated/active calibration state.
                Defaults to True (use active calibration). Set to False to use original
                calibration from data acquisition.
            cache_coordinates: Ignored, maintained for compatibility
            memory_limit_gb: Ignored, maintained for compatibility
            batch_size: Ignored, maintained for compatibility
            progress_callback: Optional callback for progress updates
            region: Region number to process for multi-region datasets.
                None (default): convert all regions (no filtering).
                int: explicitly select that region number only.
            **kwargs: Additional arguments
        """
        super().__init__(data_path, **kwargs)
        self.use_recalibrated_state = use_recalibrated_state
        self.progress_callback = progress_callback
        self._requested_region = region

        # Validate and setup paths
        self._validate_data_path()
        self._detect_file_type()

        # Read calibration metadata
        self._calibration_metadata = self._read_calibration_metadata()

        # Initialize components
        self._setup_components(cache_coordinates, memory_limit_gb, batch_size)

        # Initialize SDK and connections
        self._initialize_sdk()
        self._initialize_database()

        # Region detection and selection (must be after database init)
        self._region_info: List[Tuple[int, int]] = self._detect_regions()
        self._selected_region: Optional[int]
        self._region_frame_ids: Optional[set]
        self._selected_region, self._region_frame_ids = self._select_region()

        # Optical alignment data (parsed from .mis file + database)
        self._mis_metadata: Dict[str, Any] = self._parse_mis_alignment()
        self._positions: List[Dict[str, Any]] = self._build_positions_from_db()
        self._header: Dict[str, Any] = self._build_header_alignment()

        # Cached properties (lazy loaded)
        self._common_mass_axis: Optional[np.ndarray] = None
        self._frame_count: Optional[int] = None
        self._coordinate_offsets: Optional[Tuple[int, int, int]] = None
        self._closed: bool = False  # Track if resources have been closed

        # Preload NumPeaks cache for buffer size optimization
        self._num_peaks_cache: Dict[int, int] = self._preload_frame_num_peaks()

        cache_status = (
            f"with {len(self._num_peaks_cache)} NumPeaks cached"
            if self._num_peaks_cache
            else "with fallback spectrum reading"
        )
        logger.info(
            f"Initialized BrukerReader for {self.file_type.upper()} data at "
            f"{data_path} ({cache_status})"
        )

    def _validate_data_path(self) -> None:
        """Validate the data path and find .d directory if needed.

        If the path is not a .d directory, uses BrukerFolderStructure to
        find a .d folder within the given path. This allows users to pass
        a parent directory containing the .d folder.
        """
        if not self.data_path.exists():
            raise FileFormatError(f"Data path does not exist: {self.data_path}")

        if not self.data_path.is_dir():
            raise FileFormatError(f"Data path must be a directory: {self.data_path}")

        # If already a .d directory, use it directly
        if self.data_path.suffix.lower() == ".d":
            return

        # Otherwise, use BrukerFolderStructure to find the .d folder
        try:
            folder = BrukerFolderStructure(self.data_path)
            info = folder.analyze()

            if info.format != BrukerFormat.TIMSTOF:
                raise FileFormatError(
                    f"Path does not contain timsTOF data: {self.data_path}. "
                    f"Detected format: {info.format.value}"
                )

            # Update data_path to point to the actual .d folder
            self.data_path = info.data_path
            logger.info(f"Found .d folder at: {self.data_path}")

        except ValueError as e:
            raise FileFormatError(str(e)) from e
        except Exception as e:
            raise FileFormatError(
                f"Could not find .d folder in {self.data_path}: {e}"
            ) from e

    def _detect_file_type(self) -> None:
        """Detect whether this is TSF or TDF data."""
        tsf_path = self.data_path / "analysis.tsf"
        tdf_path = self.data_path / "analysis.tdf"

        if tsf_path.exists():
            self.file_type = "tsf"
            self.db_path = tsf_path
        elif tdf_path.exists():
            self.file_type = "tdf"
            self.db_path = tdf_path
        else:
            raise FileFormatError(
                f"No analysis.tsf or analysis.tdf found in {self.data_path}"
            )

        logger.debug(f"Detected file type: {self.file_type.upper()}")

    def _setup_components(
        self,
        cache_coordinates: bool,
        memory_limit_gb: Optional[float],
        batch_size: Optional[int],
    ) -> None:
        """Setup utility components (now minimal)."""
        # No more coordinate cache - using direct database access
        pass

    def _read_calibration_metadata(self) -> Optional[Dict]:
        """Read calibration metadata from calibration.sqlite.

        Returns:
            Dictionary containing calibration metadata, or None if unavailable.
            Keys include:
            - calibration_id: ID of active calibration state
            - calibration_uuid: Unique identifier for this calibration
            - calibration_datetime: When calibration was performed
            - calibration_source: Software that created calibration
            - num_calibration_versions: Total number of calibration states
            - recalibrated: Whether data has been recalibrated after acquisition
            - original_calibration_datetime: Original calibration datetime (if recalibrated)
        """
        cal_file = self.data_path / "calibration.sqlite"

        if not cal_file.exists():
            logger.warning(f"No calibration.sqlite found in {self.data_path}")
            return None

        try:
            conn = sqlite3.connect(cal_file)
            cursor = conn.cursor()

            # Count total calibration versions
            cursor.execute("SELECT COUNT(*) FROM CalibrationState")
            num_versions = cursor.fetchone()[0]

            # Get ACTIVE calibration (highest ID = most recent)
            cursor.execute(
                """
                SELECT Id, Key, DateTime, Source
                FROM CalibrationState
                ORDER BY Id DESC LIMIT 1
            """
            )
            cal_id, cal_uuid, cal_datetime, cal_source = cursor.fetchone()

            # Get original calibration if recalibrated
            original_datetime = None
            if num_versions > 1:
                cursor.execute(
                    """
                    SELECT DateTime FROM CalibrationState
                    ORDER BY Id ASC LIMIT 1
                """
                )
                original_datetime = cursor.fetchone()[0]

            # Get additional metadata from CalibrationInfo
            cursor.execute(
                """
                SELECT KeyName, Value
                FROM CalibrationInfo
                WHERE CalibrationState = ?
                AND KeyName IN ('CalibrationSoftwareVersion', 'CalibrationUser')
            """,
                (cal_id,),
            )

            extra_info = dict(cursor.fetchall())

            conn.close()

            metadata = {
                "calibration_id": cal_id,
                "calibration_uuid": cal_uuid,
                "calibration_datetime": cal_datetime,
                "calibration_source": cal_source,
                "calibration_software_version": extra_info.get(
                    "CalibrationSoftwareVersion"
                ),
                "calibration_user": extra_info.get("CalibrationUser"),
                "num_calibration_versions": num_versions,
                "recalibrated": num_versions > 1,
                "original_calibration_datetime": original_datetime,
                "calibration_file_size": cal_file.stat().st_size,
            }

            # Log which calibration is being used
            if self.use_recalibrated_state:
                recal_info = (
                    f" (recalibrated {num_versions} times)" if num_versions > 1 else ""
                )
                logger.info(
                    f"Using active calibration state {cal_id} from {cal_datetime}"
                    f"{recal_info}"
                )
            else:
                active_info = f", active state is {cal_id}" if num_versions > 1 else ""
                logger.info(
                    f"Using original calibration (use_recalibrated_state=False)"
                    f"{active_info}"
                )

            return metadata

        except Exception as e:
            logger.error(f"Failed to read calibration metadata: {e}")
            return None

    def _initialize_sdk(self) -> None:
        """Initialize the Bruker SDK with error handling."""
        try:
            # Initialize DLL manager
            self.dll_manager = DLLManager(
                data_directory=self.data_path, force_reload=False
            )

            # Initialize SDK functions
            self.sdk = SDKFunctions(self.dll_manager, self.file_type)

            # Open the data file
            self.handle = self.sdk.open_file(
                str(self.data_path), self.use_recalibrated_state
            )

            logger.debug(f"Successfully initialized {self.file_type.upper()} SDK")

        except Exception as e:
            logger.error(f"Failed to initialize SDK: {e}")
            raise SDKError(f"Failed to initialize Bruker SDK: {e}") from e

    def _initialize_database(self) -> None:
        """Initialize database connection with optimizations."""
        try:
            # Open database in read-only mode to avoid locking issues
            # This allows reading from network drives and concurrent access
            db_uri = f"file:{self.db_path}?mode=ro&immutable=1"
            self.conn = sqlite3.connect(
                db_uri, uri=True, timeout=30.0, check_same_thread=False
            )

            # Apply read-only compatible SQLite optimizations
            # Note: journal_mode and synchronous are not needed for read-only access
            self.conn.execute("PRAGMA cache_size = 10000")
            self.conn.execute("PRAGMA temp_store = MEMORY")

            logger.debug("Initialized database connection in read-only mode")

        except sqlite3.OperationalError as e:
            if "database is locked" in str(e):
                logger.error(
                    f"Database is locked. This may occur if another process "
                    f"(e.g., DataAnalysis) has the file open: {self.db_path}"
                )
                raise DataError(
                    f"Database is locked. Please close any other applications "
                    f"that may have this dataset open: {e}"
                ) from e
            elif "unable to open database file" in str(e):
                logger.error(f"Cannot access database file: {self.db_path}")
                raise DataError(f"Cannot access database file: {e}") from e
            else:
                raise DataError(f"Failed to open database: {e}") from e
        except Exception as e:
            logger.error(f"Failed to initialize database: {e}")
            raise DataError(f"Failed to open database: {e}") from e

    def _create_metadata_extractor(self) -> MetadataExtractor:
        """Create Bruker metadata extractor."""
        if not hasattr(self, "conn") or self.conn is None:
            raise ValueError("Database connection not available")
        # Pass selected region so extractor computes per-region bounds
        region_for_extractor = None
        if self._region_frame_ids is not None:
            region_for_extractor = self._selected_region
        return BrukerMetadataExtractor(
            self.conn,
            self.data_path,
            self._calibration_metadata,
            region=region_for_extractor,
        )

    def _detect_regions(self) -> List[Tuple[int, int]]:
        """Detect available regions in the dataset.

        Returns:
            List of (region_number, frame_count) tuples sorted by frame
            count descending. Empty list if no region info available.
        """
        try:
            cursor = self.conn.cursor()
            cursor.execute(
                "SELECT RegionNumber, COUNT(*) as n_frames "
                "FROM MaldiFrameInfo "
                "GROUP BY RegionNumber "
                "ORDER BY COUNT(*) DESC"
            )
            regions = [(int(r), int(n)) for r, n in cursor.fetchall()]
            if len(regions) > 1:
                region_desc = ", ".join(
                    f"Region {r} ({n:,} frames)" for r, n in regions
                )
                logger.info(f"Detected {len(regions)} regions: {region_desc}")
            return regions
        except sqlite3.OperationalError:
            return []

    def _select_region(self) -> Tuple[Optional[int], Optional[set]]:
        """Select region based on user request or convert all regions.

        By default, all regions are converted (no filtering). The user can
        explicitly select a single region via the region= parameter.

        Returns:
            Tuple of (selected_region_number, set_of_frame_ids).
            Both are None when no filtering is needed (default: all regions).
        """
        if not self._region_info or len(self._region_info) <= 1:
            return (None, None)

        # Multiple regions exist
        valid_regions = [r for r, _ in self._region_info]
        total_spectra = sum(n for _, n in self._region_info)

        if self._requested_region is not None:
            # User explicitly requested a specific region
            if self._requested_region not in valid_regions:
                raise ValueError(
                    f"Region {self._requested_region} not found. "
                    f"Available regions: {valid_regions}"
                )
            selected = self._requested_region

            # Get frame IDs for selected region
            cursor = self.conn.cursor()
            cursor.execute(
                "SELECT Frame FROM MaldiFrameInfo WHERE RegionNumber = ?",
                (selected,),
            )
            frame_ids = {int(row[0]) for row in cursor.fetchall()}
            logger.info(
                f"Region {selected}: {len(frame_ids):,} frames selected "
                f"(use region=None to convert all regions)"
            )
            return (selected, frame_ids)

        # Default: convert all regions, no filtering
        logger.info(
            f"Multi-region dataset: converting all {len(self._region_info)} "
            f"regions ({total_spectra:,} total spectra). "
            f"Use region= parameter to select a specific region."
        )
        return (None, None)

    def get_region_map(self) -> Optional[Dict[tuple, int]]:
        """Get per-pixel region mapping from MaldiFrameInfo.

        Maps each normalized (0-based) (x, y) coordinate to its acquisition
        region number. Uses the same coordinate offsets as iter_spectra() so
        the mapping is consistent with obs indices.

        Returns:
            Dict mapping (x, y) tuples to region numbers, or None if
            region information is not available (single-region dataset).
        """
        if not self._region_info or len(self._region_info) <= 1:
            return None

        coordinate_offsets = self._get_coordinate_offsets()
        region_map: Dict[tuple, int] = {}

        try:
            cursor = self.conn.cursor()
            if self._selected_region is not None:
                cursor.execute(
                    "SELECT XIndexPos, YIndexPos, RegionNumber "
                    "FROM MaldiFrameInfo WHERE RegionNumber = ?",
                    (self._selected_region,),
                )
            else:
                cursor.execute(
                    "SELECT XIndexPos, YIndexPos, RegionNumber " "FROM MaldiFrameInfo"
                )

            for x, y, region_num in cursor.fetchall():
                nx, ny = int(x), int(y)
                if coordinate_offsets:
                    nx -= coordinate_offsets[0]
                    ny -= coordinate_offsets[1]
                region_map[(nx, ny)] = int(region_num)

            logger.info(
                f"Built region map: {len(region_map)} pixels across "
                f"{len(set(region_map.values()))} regions"
            )
            return region_map

        except sqlite3.OperationalError as e:
            logger.warning(f"Could not build region map: {e}")
            return None

    def get_region_info(self) -> Optional[list]:
        """Get summary information about acquisition regions.

        Returns:
            List of region summary dicts with region_number and n_spectra,
            or None if region information is not available.
        """
        if not self._region_info or len(self._region_info) <= 1:
            return None

        return [
            {"region_number": region_num, "n_spectra": n_frames}
            for region_num, n_frames in self._region_info
        ]

    def get_common_mass_axis(self) -> NDArray[np.float64]:
        """Return the common mass axis composed of all unique m/z values.

        Returns:
            Array of unique m/z values in ascending order
        """
        if self._common_mass_axis is None:
            self._common_mass_axis = self._build_common_mass_axis()

        return self._common_mass_axis

    def _build_common_mass_axis(self) -> NDArray[np.float64]:
        """Build the common mass axis and cache total peaks."""
        logger.info("Building raw mass axis")

        # Create iterator for mass axis building
        def mz_iterator():
            for coords, mzs, intensities in self._iter_spectra_raw():
                yield coords, mzs, intensities

        # Build raw mass axis using simplified function (returns mass_axis and total_peaks)
        mass_axis, total_peaks = build_raw_mass_axis(
            mz_iterator(), self.progress_callback  # type: ignore[arg-type]
        )

        # Cache total_peaks for later retrieval (if not already set from NumPeaks cache)
        if not hasattr(self, "_total_peaks_from_mass_axis"):
            self._total_peaks_from_mass_axis = total_peaks

        if len(mass_axis) == 0:
            logger.warning("No m/z values found in dataset")
            return np.array([])

        logger.info(f"Built raw mass axis with {len(mass_axis)} unique m/z values")
        return mass_axis

    def iter_spectra(self, batch_size: Optional[int] = None) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Iterate through all spectra sequentially.

        Args:
            batch_size: Ignored, maintained for compatibility

        Yields:
            Tuples of (coordinates, mz_array, intensity_array)
        """
        # Always use simple sequential iteration
        yield from self._iter_spectra_raw()

    def _iter_spectra_raw(
        self,
    ) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Raw spectrum iteration without batching.

        When region filtering is active, only frames belonging to the
        selected region are iterated.
        """
        # Determine which frames to iterate
        if self._region_frame_ids is not None:
            frame_ids: Any = sorted(self._region_frame_ids)
            total = len(frame_ids)
        else:
            total = self._get_frame_count()
            frame_ids = range(1, total + 1)

        coordinate_offsets = self._get_coordinate_offsets()

        # Setup progress tracking
        with tqdm(
            total=total,
            desc="Reading spectra",
            unit="spectrum",
            disable=True,  # Disable to avoid double progress with converter
        ) as pbar:
            for frame_id in frame_ids:
                try:
                    # Get normalized coordinates using persistent connection
                    coords = self._get_frame_coordinates_cached(
                        frame_id, coordinate_offsets
                    )
                    if coords is None:
                        logger.warning(f"No coordinates found for frame {frame_id}")
                        pbar.update(1)
                        continue

                    # OPTIMIZATION: Get buffer size hint from NumPeaks cache
                    buffer_size_hint = self._num_peaks_cache.get(frame_id)

                    # Read spectrum with optimization (or fallback if no hint)
                    mzs, intensities = self.sdk.read_spectrum(
                        self.handle,
                        frame_id,
                        buffer_size_hint=buffer_size_hint,
                    )

                    # Apply intensity threshold filtering if configured
                    mzs, intensities = self._apply_intensity_filter(mzs, intensities)

                    if mzs.size > 0 and intensities.size > 0:
                        yield coords, mzs, intensities

                    pbar.update(1)

                    # Progress callback
                    if self.progress_callback:
                        self.progress_callback(frame_id, total)

                except Exception as e:
                    logger.warning(f"Error reading spectrum for frame {frame_id}: {e}")
                    pbar.update(1)
                    continue

    def _get_frame_count(self) -> int:
        """Get the total number of frames (respects region filtering)."""
        if self._frame_count is None:
            if self._region_frame_ids is not None:
                self._frame_count = len(self._region_frame_ids)
            else:
                self._frame_count = _get_frame_count(self.db_path)

        return self._frame_count

    def _get_coordinate_offsets(self) -> Optional[Tuple[int, int, int]]:
        """Get coordinate offsets from metadata for normalization."""
        if self._coordinate_offsets is None:
            essential_metadata = self.get_essential_metadata()
            self._coordinate_offsets = essential_metadata.coordinate_offsets

        return self._coordinate_offsets

    def _get_frame_coordinates_cached(
        self,
        frame_id: int,
        coordinate_offsets: Optional[Tuple[int, int, int]] = None,
    ) -> Optional[Tuple[int, int, int]]:
        """Get normalized coordinates for a specific frame using persistent connection.

        This avoids opening new SQLite connections for every frame.

        Args:
            frame_id: Frame ID to look up
            coordinate_offsets: Optional coordinate offsets for normalization

        Returns:
            Tuple of normalized (x, y, z) coordinates (0-based), or None if not
            found
        """
        try:
            cursor = self.conn.cursor()

            # Check if this is MALDI data
            try:
                cursor.execute(
                    "SELECT XIndexPos, YIndexPos FROM MaldiFrameInfo WHERE "
                    "Frame = ?",
                    (frame_id,),
                )
                result = cursor.fetchone()
                if result:
                    x, y = result
                    # Apply coordinate offsets if provided (Bruker-specific
                    # normalization)
                    if coordinate_offsets:
                        offset_x, offset_y, offset_z = coordinate_offsets
                        return (int(x) - offset_x, int(y) - offset_y, 0)
                    else:
                        return (int(x), int(y), 0)
            except sqlite3.OperationalError:
                # No MALDI table, use generated coordinates
                pass

            # For non-MALDI data, generate coordinates (simple sequential
            # mapping)
            return (frame_id - 1, 0, 0)

        except Exception as e:
            logger.warning(f"Error getting coordinates for frame {frame_id}: {e}")
            return None

    def get_frame_id_by_coordinates(self, x: int, y: int, z: int = 0) -> Optional[int]:
        """Get frame ID for given coordinates.

        Args:
            x: X coordinate (0-based)
            y: Y coordinate (0-based)
            z: Z coordinate (0-based, default 0 for 2D data)

        Returns:
            Frame ID if found, None otherwise
        """
        try:
            # Adjust coordinates by offsets if needed
            offsets = self._get_coordinate_offsets()
            if offsets:
                x_offset, y_offset, z_offset = offsets
                x += x_offset
                y += y_offset
                z += z_offset

            with sqlite3.connect(str(self.db_path)) as conn:
                cursor = conn.cursor()
                cursor.execute(
                    "SELECT Frame FROM MaldiFrameInfo WHERE "
                    "XIndexPos = ? AND YIndexPos = ?",
                    (x, y),
                )
                result = cursor.fetchone()
                if result:
                    return int(result[0])
                return None

        except Exception as e:
            logger.warning(
                f"Error getting frame ID for coordinates ({x}, {y}, {z}): {e}"
            )
            return None

    def get_spectrum_by_coordinates(
        self, x: int, y: int, z: int = 0
    ) -> Optional[Tuple[NDArray[np.float64], NDArray[np.float64]]]:
        """Get spectrum data for given coordinates.

        Args:
            x: X coordinate (0-based)
            y: Y coordinate (0-based)
            z: Z coordinate (0-based, default 0 for 2D data)

        Returns:
            Tuple of (mz_array, intensity_array) if found, None otherwise
        """
        frame_id = self.get_frame_id_by_coordinates(x, y, z)
        if frame_id is None:
            return None

        try:
            # Get buffer size hint from cache
            buffer_size_hint = self._num_peaks_cache.get(frame_id)

            # Read spectrum
            mzs, intensities = self.sdk.read_spectrum(
                self.handle, frame_id, buffer_size_hint=buffer_size_hint
            )

            # Apply intensity threshold filtering if configured
            mzs, intensities = self._apply_intensity_filter(mzs, intensities)

            if mzs.size > 0 and intensities.size > 0:
                return mzs, intensities
            return None

        except Exception as e:
            logger.warning(f"Error reading spectrum for frame {frame_id}: {e}")
            return None

    def _preload_frame_num_peaks(self) -> Dict[int, int]:
        """Preload NumPeaks values for relevant frames at initialization.

        This optimization avoids the busy wait loop in SDK by providing
        exact buffer sizes for spectrum reading, reducing CPU usage from 100%
        to normal levels.

        When region filtering is active, only caches peaks for frames in
        the selected region.

        Returns:
            Dictionary mapping frame_id -> num_peaks (validated to be <= 65535)
        """
        try:
            with sqlite3.connect(str(self.db_path)) as conn:
                cursor = conn.cursor()
                cursor.execute("SELECT Id, NumPeaks FROM Frames ORDER BY Id")

                # Use uint16 equivalent for memory efficiency (max 65,535
                # peaks)
                num_peaks_cache = {}
                invalid_count = 0

                for frame_id, num_peaks in cursor.fetchall():
                    # Skip frames not in selected region
                    if (
                        self._region_frame_ids is not None
                        and frame_id not in self._region_frame_ids
                    ):
                        continue

                    if num_peaks is not None and 0 < num_peaks <= 65535:
                        num_peaks_cache[frame_id] = int(num_peaks)
                    else:
                        invalid_count += 1
                        if invalid_count <= 5:  # Log first few invalid values
                            logger.debug(
                                f"Invalid NumPeaks value {num_peaks} for "
                                f"frame {frame_id}"
                            )

                if invalid_count > 5:
                    logger.debug(
                        f"... and {invalid_count - 5} more invalid NumPeaks " f"values"
                    )

                memory_mb = len(num_peaks_cache) * 2 / (1024 * 1024)  # uint16 = 2 bytes
                logger.debug(
                    f"Cached NumPeaks for {len(num_peaks_cache)} frames "
                    f"({memory_mb:.1f}MB)"
                )
                return num_peaks_cache

        except Exception as e:
            logger.warning(f"Failed to preload NumPeaks cache: {e}")
            logger.info("Will use fallback retry logic for spectrum reading")
            return {}  # Empty cache triggers fallback behavior

    def close(self) -> None:
        """Close all resources and connections.

        This method is idempotent - safe to call multiple times.
        """
        # Skip if already closed
        if getattr(self, "_closed", False):
            logger.debug("Resources already closed, skipping")
            return

        logger.debug("Closing Bruker reader")

        try:
            # Close SDK handle
            if hasattr(self, "handle") and self.handle:
                self.sdk.close_file(self.handle)
                self.handle = None

            # Close database connection
            if hasattr(self, "conn") and self.conn:
                self.conn.close()
                self.conn = None

            # Mark as closed
            self._closed = True

            logger.info("Successfully closed all resources")

        except Exception as e:
            logger.error(f"Error during cleanup: {e}")
            # Even if there's an error, mark as closed to prevent repeated attempts
            self._closed = True

    @property
    def shape(self) -> Tuple[int, int, int]:
        """Get spatial dimensions (pixel grid) from metadata extractor.

        Note: This is the spatial pixel grid, not the mass axis
        dimensions.
        Mass axis interpolation to common m/z values is handled during
        conversion.

        Returns:
            Tuple of (x_pixels, y_pixels, z_pixels) spatial dimensions
        """
        essential_metadata = self.get_essential_metadata()
        return essential_metadata.dimensions

    @property
    def mass_range(self) -> Tuple[float, float]:
        """Get mass range from metadata extractor.

        Note: This is the acquisition mass range, not the final
        interpolated axis.
        The actual common mass axis for interpolation is built from all
        unique m/z values.

        Returns:
            Tuple of (min_mass, max_mass) in m/z units
        """
        essential_metadata = self.get_essential_metadata()
        return essential_metadata.mass_range

    def __repr__(self) -> str:
        """String representation of the reader."""
        return (
            f"BrukerReader(path={self.data_path}, "
            f"type={self.file_type.upper()}, "
            f"frames={self._get_frame_count()})"
        )

    @property
    def mis_metadata(self) -> Dict[str, Any]:
        """Get parsed .mis metadata for optical alignment."""
        return self._mis_metadata

    def _parse_mis_alignment(self) -> Dict[str, Any]:
        """Parse .mis file for optical alignment if available.

        Uses BrukerFolderStructure to find the .mis file in the parent
        folder hierarchy, then parses it for area definitions and
        teaching points.

        Returns:
            Dictionary of parsed .mis metadata, empty if no .mis found
        """
        try:
            mis_path = self.get_teaching_points_file()
        except (ValueError, OSError):
            # Folder structure analysis can fail for non-standard paths
            return {}

        if mis_path is None:
            return {}

        logger.info(f"Found .mis file for optical alignment: {mis_path.name}")
        metadata = parse_mis_file(mis_path)

        if metadata.get("areas"):
            logger.info(
                f"Parsed {len(metadata['areas'])} area definitions " f"from .mis file"
            )
        return metadata

    def _build_positions_from_db(self) -> List[Dict[str, Any]]:
        """Build position list from MaldiFrameInfo for alignment.

        Creates position dictionaries compatible with the alignment module
        by querying raster coordinates and region numbers from the database.

        Returns:
            List of position dicts with region, raster_x, raster_y keys.
            Empty list if no .mis areas are available.
        """
        if not self._mis_metadata.get("areas"):
            return []

        positions: List[Dict[str, Any]] = []
        try:
            cursor = self.conn.cursor()
            cursor.execute(
                "SELECT XIndexPos, YIndexPos, RegionNumber " "FROM MaldiFrameInfo"
            )
            for x, y, region in cursor.fetchall():
                positions.append(
                    {
                        "region": int(region),
                        "raster_x": int(x),
                        "raster_y": int(y),
                    }
                )
            logger.info(
                f"Built {len(positions)} positions from MaldiFrameInfo "
                f"for optical alignment"
            )
        except sqlite3.OperationalError as e:
            logger.warning(f"Could not query positions for alignment: {e}")

        return positions

    def _build_header_alignment(self) -> Dict[str, Any]:
        """Build header dict with first_raster offsets for alignment.

        The first_raster_x/y values must match the coordinate offsets used
        by BrukerMetadataExtractor so that transform_point() correctly
        de-normalizes 0-based raster coordinates back to originals.

        When a specific region is selected, uses per-region minimums from
        MaldiFrameInfo. Otherwise uses global ImagingArea from GlobalMetadata.

        Returns:
            Dictionary with first_raster_x and first_raster_y keys.
            Empty dict if no .mis areas are available.
        """
        if not self._mis_metadata.get("areas"):
            return {}

        try:
            cursor = self.conn.cursor()

            if self._selected_region is not None:
                # Per-region: use actual min from the selected region's frames
                # This matches BrukerMetadataExtractor's per-region offsets
                cursor.execute(
                    "SELECT MIN(XIndexPos), MIN(YIndexPos) "
                    "FROM MaldiFrameInfo WHERE RegionNumber = ?",
                    (self._selected_region,),
                )
                result = cursor.fetchone()
                if result and result[0] is not None and result[1] is not None:
                    header = {
                        "first_raster_x": int(result[0]),
                        "first_raster_y": int(result[1]),
                    }
                else:
                    return {}
            else:
                # Global: use ImagingArea from GlobalMetadata
                cursor.execute(
                    "SELECT Value FROM GlobalMetadata "
                    "WHERE Key = 'ImagingAreaMinXIndexPos'"
                )
                result_x = cursor.fetchone()
                cursor.execute(
                    "SELECT Value FROM GlobalMetadata "
                    "WHERE Key = 'ImagingAreaMinYIndexPos'"
                )
                result_y = cursor.fetchone()

                if result_x and result_y:
                    header = {
                        "first_raster_x": int(float(result_x[0])),
                        "first_raster_y": int(float(result_y[0])),
                    }
                else:
                    return {}

            logger.info(
                f"Alignment offsets: first_raster=("
                f"{header['first_raster_x']}, "
                f"{header['first_raster_y']})"
            )
            return header
        except (sqlite3.OperationalError, TypeError) as e:
            logger.warning(f"Could not get alignment offsets: {e}")

        return {}

    @property
    def n_spectra(self) -> int:
        """Return the total number of spectra in the dataset.

        Returns:
            Total number of frames (efficient implementation using cached
            frame count)
        """
        return self._get_frame_count()

    def get_total_peak_count(self) -> int:
        """Get total number of peaks across all spectra from NumPeaks cache.

        This is very fast as NumPeaks data is cached from the database at
        initialization.

        Returns:
            Total number of peaks across all spectra
        """
        if not self._num_peaks_cache:
            logger.warning("NumPeaks cache not available, cannot get exact count")
            return 0

        total = sum(self._num_peaks_cache.values())
        logger.info(f"Total peak count from NumPeaks cache: {total:,}")
        return total

    def get_peak_counts_per_pixel(self) -> Optional[np.ndarray]:
        """Get per-pixel peak counts for CSR indptr construction.

        Converts the frame-indexed NumPeaks cache to pixel-indexed array
        using coordinate mapping. When region filtering is active, only
        includes frames from the selected region.

        Returns:
            Array of size n_pixels where arr[pixel_idx] = peak_count.
            pixel_idx = z * (n_x * n_y) + y * n_x + x
            Returns None if NumPeaks cache not available.
        """
        if not self._num_peaks_cache:
            logger.warning("NumPeaks cache not available")
            return None

        # Get dimensions and coordinate offsets
        metadata = self.get_essential_metadata()
        n_x, n_y, n_z = metadata.dimensions
        n_pixels = n_x * n_y * n_z
        coordinate_offsets = metadata.coordinate_offsets

        # Create output array
        peak_counts = np.zeros(n_pixels, dtype=np.int32)

        # Map frame_id -> pixel_idx using coordinate lookup
        cursor = self.conn.cursor()
        try:
            if self._selected_region is not None:
                cursor.execute(
                    "SELECT Frame, XIndexPos, YIndexPos "
                    "FROM MaldiFrameInfo "
                    "WHERE RegionNumber = ?",
                    (self._selected_region,),
                )
            else:
                cursor.execute(
                    "SELECT Frame, XIndexPos, YIndexPos " "FROM MaldiFrameInfo"
                )
            for frame_id, x, y in cursor.fetchall():
                if frame_id not in self._num_peaks_cache:
                    continue

                # Apply coordinate offsets (normalize to 0-based)
                if coordinate_offsets:
                    x = int(x) - coordinate_offsets[0]
                    y = int(y) - coordinate_offsets[1]
                else:
                    x, y = int(x), int(y)

                # Calculate pixel index
                z = 0  # Bruker MSI is typically 2D
                pixel_idx = z * (n_x * n_y) + y * n_x + x

                if 0 <= pixel_idx < n_pixels:
                    peak_counts[pixel_idx] = self._num_peaks_cache[frame_id]

        except Exception as e:
            logger.warning(f"Error mapping peak counts to pixels: {e}")
            return None

        logger.info(f"Mapped peak counts for {n_pixels:,} pixels")
        return peak_counts

    def __del__(self) -> None:
        """Destructor to ensure cleanup."""
        try:
            self.close()
        except Exception as e:
            logger.debug(
                f"Error during cleanup in destructor: {e}"
            )  # Log but don't raise during destruction
