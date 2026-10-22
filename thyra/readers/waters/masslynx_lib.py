# thyra/readers/waters/masslynx_lib.py
"""Python ctypes bindings to the Waters MassLynx native libraries.

Wraps MLReader.dll/libMLReader.so (API) and MassLynxRaw.dll/libMassLynxRaw.so
(base dependency) to provide access to Waters .raw MSI data.

The C API and struct layouts were reverse-engineered from mzmine's jextract-
generated Java bindings (MassLynxLib.java, ScanInfo.java).
"""

import ctypes
import logging
import os
import platform
import time
from ctypes import (
    POINTER,
    Structure,
    c_char_p,
    c_double,
    c_float,
    c_int32,
    c_uint32,
    c_void_p,
)
from dataclasses import dataclass
from enum import Enum, auto
from math import isfinite
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# Sentinel constants (from MassLynxConstants.java)
NO_POSITION: float = -1.0
NO_LOCKMASS_FUNCTION: int = -1
DEFAULT_FLOAT: float = -1.0


class FunctionType(Enum):
    """Waters acquisition function type (from FunctionType.java)."""

    MS = auto()
    IMS_MS = auto()
    MRM = auto()
    LOCKMASS = auto()
    NOT_MS = auto()


class ScanInfoStruct(Structure):
    """Maps to the C struct ScanInfo (44 bytes, no padding).

    Offsets confirmed from mzmine's ScanInfo.java:
      msLevel: 0, polarity: 4, driftScanCount: 8, isProfile: 12,
      precursorMz: 16, quadIsolationStart: 20, quadIsolationEnd: 24,
      collisionEnergy: 28, rt: 32, laserXPos: 36, laserYPos: 40
    """

    # Python 3.14 warns about `_pack_` without `_layout_` off Windows, and
    # 3.19 will refuse it. "ms" is the layout `_pack_` already implied, and
    # older Pythons ignore the attribute.
    _layout_ = "ms"
    _pack_ = 1
    _fields_ = [
        ("msLevel", c_int32),
        ("polarity", c_int32),
        ("driftScanCount", c_int32),
        ("isProfile", c_int32),
        ("precursorMz", c_float),
        ("quadIsolationStart", c_float),
        ("quadIsolationEnd", c_float),
        ("collisionEnergy", c_float),
        ("rt", c_float),
        ("laserXPos", c_float),
        ("laserYPos", c_float),
    ]


@dataclass(frozen=True)
class ScanInfoData:
    """Python representation of a scan's metadata from the ScanInfo C struct."""

    ms_level: int
    polarity: int
    drift_scan_count: int
    is_profile: int
    precursor_mz: float
    rt: float
    laser_x_pos: float  # in mm from DLL
    laser_y_pos: float  # in mm from DLL
    #: The quadrupole isolation window and collision energy of an MS/MS
    #: scan, as MassLynx reports them; 0.0 when a Synapt scan is MS1. A
    #: Cyclic IMS reports a precursor on full scans too, so the function's
    #: declared ``MSMSMode`` decides there (D72). Defaulted
    #: so a scan record built without them (older callers, tests) is intact.
    quad_isolation_start: float = 0.0
    quad_isolation_end: float = 0.0
    collision_energy: float = 0.0

    @property
    def has_position(self) -> bool:
        """True if this scan has a usable laser position.

        Two ways to fail. The documented one is the -1.0 mm sentinel on
        both axes. The other is a non-finite reading: NaN is not equal to
        anything including itself, so it passed the sentinel test, and
        every NaN then became its own column in the grid -- one of them
        made the reported pixel size NaN, two of them changed the pitch of
        the real axis, and the converter's ``pixel_size_um <= 0`` guards
        are False for NaN and infinity, so such a pitch reached the store
        (issue #231).
        """
        x, y = self.laser_x_pos, self.laser_y_pos
        if not (isfinite(x) and isfinite(y)):
            return False
        return not (x == NO_POSITION and y == NO_POSITION)


class WatersLibError(Exception):
    """Raised when Waters DLL operations fail."""

    pass


class MassLynxLib:
    """Wrapper around Waters MassLynx native libraries.

    Loads MassLynxRaw (base dependency) first, then MLReader (API),
    matching the mzmine loading pattern in MassLynxLib.java.
    """

    _instance: Optional["MassLynxLib"] = None

    def __init__(self, lib_dir: Optional[Path] = None):
        """Initialize MassLynxLib and load native libraries.

        Args:
            lib_dir: Optional directory to search for native libraries.
        """
        self._base_lib: Optional[ctypes.CDLL] = None
        self._api_lib: Optional[ctypes.CDLL] = None
        self._load_libraries(lib_dir)
        self._setup_function_signatures()

    @classmethod
    def get_instance(cls, lib_dir: Optional[Path] = None) -> "MassLynxLib":
        """Singleton accessor. Creates instance on first call."""
        if cls._instance is None:
            cls._instance = cls(lib_dir)
        return cls._instance

    @classmethod
    def reset_instance(cls) -> None:
        """Reset singleton (for testing)."""
        cls._instance = None

    def _get_library_names(self) -> Tuple[str, str]:
        """Return (base_lib_name, api_lib_name) for current platform."""
        if platform.system() == "Windows":
            return ("MassLynxRaw.dll", "MLReader.dll")
        return ("libMassLynxRaw.so", "libMLReader.so")

    def _get_search_paths(self, lib_dir: Optional[Path]) -> list:
        """Get ordered list of directories to search for libraries."""
        paths = []
        # 1. Bundled lib directory (highest priority)
        bundled = Path(__file__).parent / "lib"
        if bundled.is_dir():
            paths.append(bundled)
        # 2. Explicit lib_dir argument
        if lib_dir and lib_dir.is_dir():
            paths.append(lib_dir)
        # 3. Environment variable
        env_path = os.environ.get("WATERS_SDK_PATH")
        if env_path:
            p = Path(env_path)
            if p.is_dir():
                paths.append(p)
        return paths

    def _load_libraries(self, lib_dir: Optional[Path]) -> None:
        """Load both libraries in correct order (base first, then API)."""
        base_name, api_name = self._get_library_names()
        search_paths = self._get_search_paths(lib_dir)

        # Load base library first (dependency of API library)
        self._base_lib = self._find_and_load(base_name, search_paths)
        # Then load API library
        self._api_lib = self._find_and_load(api_name, search_paths)

    def _find_and_load(self, lib_name: str, search_paths: list):
        """Find library in search paths and load it."""
        errors = []
        for search_dir in search_paths:
            lib_path = search_dir / lib_name
            if lib_path.exists() and lib_path.is_file():
                try:
                    lib = ctypes.cdll.LoadLibrary(str(lib_path))
                    logger.info(f"Loaded Waters library: {lib_path}")
                    return lib
                except OSError as e:
                    errors.append(f"{lib_path}: {e}")
                    logger.debug(f"Failed to load {lib_path}: {e}")

        searched = [str(p / lib_name) for p in search_paths]
        raise WatersLibError(
            f"Could not find or load {lib_name}. "
            f"Searched: {searched}. "
            f"Errors: {errors}. "
            f"Place the library in thyra/readers/waters/lib/ or set "
            f"WATERS_SDK_PATH environment variable."
        )

    def _ensure_loaded(self) -> ctypes.CDLL:
        """Assert that the API library is loaded and return it.

        This is a type-narrowing helper that converts Optional[ctypes.CDLL]
        to ctypes.CDLL for mypy.
        """
        if self._api_lib is None:
            raise RuntimeError("API library not loaded")
        return self._api_lib

    def _setup_function_signatures(self) -> None:
        """Define argtypes/restype for all DLL functions."""
        if self._api_lib is None:
            raise RuntimeError("API library not loaded")
        dll = self._api_lib

        # -- File open/close --
        dll.openFile.argtypes = [c_char_p]
        dll.openFile.restype = c_void_p

        dll.closeFile.argtypes = [c_void_p]
        dll.closeFile.restype = c_uint32

        # -- File metadata --
        dll.getAcquisitionDate.argtypes = [c_void_p, c_char_p, c_int32]
        dll.getAcquisitionDate.restype = c_int32

        dll.getNumberOfFunctions.argtypes = [c_void_p]
        dll.getNumberOfFunctions.restype = c_uint32

        dll.getNumberOfScans.argtypes = [c_void_p]
        dll.getNumberOfScans.restype = c_uint32

        dll.getNumberOfScansInFunction.argtypes = [c_void_p, c_uint32]
        dll.getNumberOfScansInFunction.restype = c_uint32

        # -- File type queries --
        dll.isDdaFile.argtypes = [c_void_p]
        dll.isDdaFile.restype = c_uint32

        dll.isIonMobilityFile.argtypes = [c_void_p]
        dll.isIonMobilityFile.restype = c_uint32

        dll.isImagingFile.argtypes = [c_void_p]
        dll.isImagingFile.restype = c_uint32

        dll.isSonarFile.argtypes = [c_void_p]
        dll.isSonarFile.restype = c_uint32

        dll.isMsFunction.argtypes = [c_void_p, c_int32]
        dll.isMsFunction.restype = c_uint32

        dll.isIonMobilityFunction.argtypes = [c_void_p, c_int32]
        dll.isIonMobilityFunction.restype = c_int32

        dll.isRawSpectrumContinuum.argtypes = [c_void_p, c_int32]
        dll.isRawSpectrumContinuum.restype = c_int32

        # -- Function metadata --
        dll.getLockmassFunction.argtypes = [c_void_p]
        dll.getLockmassFunction.restype = c_int32

        dll.getNumberOfMrmsInFunction.argtypes = [c_void_p, c_int32]
        dll.getNumberOfMrmsInFunction.restype = c_int32

        dll.getAcquisitionRangeStart.argtypes = [c_void_p, c_uint32]
        dll.getAcquisitionRangeStart.restype = c_double

        dll.getAcquisitionRangeEnd.argtypes = [c_void_p, c_uint32]
        dll.getAcquisitionRangeEnd.restype = c_double

        # -- Lock mass --
        dll.isLockmassCorrected.argtypes = [c_void_p]
        dll.isLockmassCorrected.restype = c_uint32

        dll.applyAutoLockmassCorrection.argtypes = [c_void_p]
        dll.applyAutoLockmassCorrection.restype = c_uint32

        dll.applyCustomLockmassCorrection.argtypes = [c_void_p, c_float]
        dll.applyCustomLockmassCorrection.restype = c_uint32

        # -- Processing options --
        dll.setCentroid.argtypes = [c_void_p, c_int32]
        dll.setCentroid.restype = c_uint32

        dll.setAbsoluteThreshold.argtypes = [c_void_p, c_float]
        dll.setAbsoluteThreshold.restype = None

        # -- Scan data --
        dll.getScanInfo.argtypes = [c_void_p, c_int32, c_int32, POINTER(ScanInfoStruct)]
        dll.getScanInfo.restype = None

        dll.getDataPoints.argtypes = [
            c_void_p,
            c_int32,
            c_int32,
            POINTER(c_double),
            POINTER(c_double),
            c_uint32,
        ]
        dll.getDataPoints.restype = c_uint32

        # -- Ion mobility (drift bins) --
        dll.getNumberOfMobilityScans.argtypes = [c_void_p, c_int32]
        dll.getNumberOfMobilityScans.restype = c_uint32

        dll.getMobility.argtypes = [c_void_p, c_int32]
        dll.getMobility.restype = c_float

        # Only the raw per-bin read is bound. Its centroiding sibling,
        # getMobilityScanDataPoints, is about ten times slower and raised
        # access violations on SYNAPT G1 files (see read_drift_bin).
        dll.getRawMobilityScanDataPoints.argtypes = [
            c_void_p,
            c_int32,
            c_int32,
            c_int32,
            POINTER(c_double),
            POINTER(c_double),
            c_uint32,
        ]
        dll.getRawMobilityScanDataPoints.restype = c_uint32

    def open_file(self, raw_path: str) -> c_void_p:
        """Open a Waters .raw directory. Returns opaque handle.

        Retries up to 10 times on failure (matching mzmine pattern,
        MassLynxDataAccess.java lines 149-161).
        """
        api = self._ensure_loaded()
        for attempt in range(10):
            handle = api.openFile(raw_path.encode("utf-8"))
            if handle and handle != 0:
                logger.info(f"Opened Waters file: {raw_path}")
                return ctypes.c_void_p(handle)
            logger.debug(f"Open attempt {attempt + 1}/10 failed for {raw_path}")
            time.sleep(0.05)
        raise WatersLibError(
            f"Failed to open Waters .raw file after 10 attempts: {raw_path}"
        )

    def close_file(self, handle) -> None:
        """Close a previously opened file handle."""
        self._ensure_loaded().closeFile(handle)

    def is_imaging_file(self, handle) -> bool:
        """Check if the file contains imaging (MSI) data."""
        return bool(self._ensure_loaded().isImagingFile(handle) > 0)

    def is_ion_mobility_file(self, handle) -> bool:
        """Check if the file contains ion mobility data."""
        return bool(self._ensure_loaded().isIonMobilityFile(handle) > 0)

    def is_ion_mobility_function(self, handle, function: int) -> bool:
        """Whether a function carries drift bins.

        Asked separately from :meth:`classify_function`, which reports a
        function as LOCKMASS before it asks about mobility: the tail chunk
        of a mobility raster MassLynx split across functions is named the
        lockmass function and still carries drift bins.
        """
        return bool(self._ensure_loaded().isIonMobilityFunction(handle, function) == 1)

    def is_dda_file(self, handle) -> bool:
        """Check if the file is a DDA acquisition."""
        return bool(self._ensure_loaded().isDdaFile(handle) > 0)

    def get_number_of_functions(self, handle) -> int:
        """Get the total number of acquisition functions."""
        return int(self._ensure_loaded().getNumberOfFunctions(handle))

    def get_number_of_scans(self, handle) -> int:
        """Get the total number of scans across all functions."""
        return int(self._ensure_loaded().getNumberOfScans(handle))

    def get_number_of_scans_in_function(self, handle, function: int) -> int:
        """Get the number of scans in a specific function."""
        return int(self._ensure_loaded().getNumberOfScansInFunction(handle, function))

    def classify_function(self, handle, function: int) -> FunctionType:
        """Classify a function's type.

        Priority order (from MassLynxDataAccess.java lines 329-348):
        1. LOCKMASS if function == getLockmassFunction
        2. IMS_MS if isIonMobilityFunction == 1
        3. MRM if getNumberOfMrmsInFunction > 0
        4. MS if isMsFunction > 0
        5. NOT_MS otherwise
        """
        api = self._ensure_loaded()
        lockmass_func = api.getLockmassFunction(handle)
        if function == lockmass_func:
            return FunctionType.LOCKMASS
        if api.isIonMobilityFunction(handle, function) == 1:
            return FunctionType.IMS_MS
        if api.getNumberOfMrmsInFunction(handle, function) > 0:
            return FunctionType.MRM
        if api.isMsFunction(handle, function) > 0:
            return FunctionType.MS
        return FunctionType.NOT_MS

    def set_centroid(self, handle, centroid: bool) -> None:
        """Configure whether getDataPoints returns centroided or profile data."""
        self._ensure_loaded().setCentroid(handle, 1 if centroid else 0)

    def get_scan_info(self, handle, function: int, scan: int) -> ScanInfoData:
        """Get scan metadata. Returns a frozen dataclass."""
        info = ScanInfoStruct()
        self._ensure_loaded().getScanInfo(handle, function, scan, ctypes.byref(info))
        return ScanInfoData(
            ms_level=info.msLevel,
            polarity=info.polarity,
            drift_scan_count=info.driftScanCount,
            is_profile=info.isProfile,
            precursor_mz=info.precursorMz,
            rt=info.rt,
            laser_x_pos=info.laserXPos,
            laser_y_pos=info.laserYPos,
            quad_isolation_start=info.quadIsolationStart,
            quad_isolation_end=info.quadIsolationEnd,
            collision_energy=info.collisionEnergy,
        )

    def read_spectrum(
        self,
        handle,
        function: int,
        scan: int,
        initial_buf_size: int = 4096,
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Read a complete spectrum with automatic buffer resizing.

        Implements the mzmine double-call pattern (MassLynxDataAccess.java
        lines 380-387): try with initial buffer, if numDp * 8 > bufSize,
        allocate 2x and retry.

        Returns:
            Tuple of (mzs, intensities) as float64 numpy arrays.
        """
        api = self._ensure_loaded()
        mz_buf = np.empty(initial_buf_size, dtype=np.float64)
        int_buf = np.empty(initial_buf_size, dtype=np.float64)

        num_dp = api.getDataPoints(
            handle,
            function,
            scan,
            mz_buf.ctypes.data_as(POINTER(c_double)),
            int_buf.ctypes.data_as(POINTER(c_double)),
            mz_buf.nbytes,
        )

        if num_dp * 8 > mz_buf.nbytes:
            # Buffer too small -- reallocate with 2x the needed size and retry
            buf_size = num_dp * 2
            mz_buf = np.empty(buf_size, dtype=np.float64)
            int_buf = np.empty(buf_size, dtype=np.float64)
            api.getDataPoints(
                handle,
                function,
                scan,
                mz_buf.ctypes.data_as(POINTER(c_double)),
                int_buf.ctypes.data_as(POINTER(c_double)),
                mz_buf.nbytes,
            )

        if num_dp == 0:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

        return mz_buf[:num_dp].copy(), int_buf[:num_dp].copy()

    def get_number_of_drift_bins(self, handle, function: int) -> int:
        """The number of drift bins (mobility scans) of a function."""
        return int(self._ensure_loaded().getNumberOfMobilityScans(handle, function))

    def get_drift_times(self, handle, n_bins: int) -> NDArray[np.float64]:
        """``getMobility`` for each of the first ``n_bins`` bins, in ms.

        On a Synapt this is each bin's drift time. On a Cyclic IMS it is the
        bin index times one pusher period, which the reader places on the
        run's declared bin timing (``WatersReader._arrival_times``).

        ``getMobility`` takes no function argument, and an index past the
        file's last bin raises a C++ exception inside the library (SEH
        0xE06D7363, measured), so the caller passes the count it read from
        :meth:`get_number_of_drift_bins` and nothing beyond it is asked.
        """
        api = self._ensure_loaded()
        return np.array(
            [api.getMobility(handle, i) for i in range(int(n_bins))],
            dtype=np.float64,
        )

    def read_drift_bin(
        self,
        handle,
        function: int,
        scan: int,
        drift_bin: int,
        initial_buf_size: int = 4096,
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """One drift bin's raw points, as ``(mzs, intensities)``.

        ``getRawMobilityScanDataPoints``: the non-zero points of the bin,
        whatever ``setCentroid`` says. Summed over the bins they give the
        scan's profile spectrum exactly, measured on SYNAPT G2-Si and XS
        runs; centroiding each bin instead costs about ten times as much
        and does not sum to anything the library reports. Same two-call
        buffer pattern as :meth:`read_spectrum`.
        """
        api = self._ensure_loaded()
        mz_buf = np.empty(initial_buf_size, dtype=np.float64)
        int_buf = np.empty(initial_buf_size, dtype=np.float64)

        def call() -> int:
            return int(
                api.getRawMobilityScanDataPoints(
                    handle,
                    function,
                    scan,
                    drift_bin,
                    mz_buf.ctypes.data_as(POINTER(c_double)),
                    int_buf.ctypes.data_as(POINTER(c_double)),
                    mz_buf.nbytes,
                )
            )

        num_dp = call()
        if num_dp * 8 > mz_buf.nbytes:
            mz_buf = np.empty(num_dp * 2, dtype=np.float64)
            int_buf = np.empty(num_dp * 2, dtype=np.float64)
            call()
        if num_dp == 0:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
        return mz_buf[:num_dp].copy(), int_buf[:num_dp].copy()

    def get_acquisition_range(
        self, handle, function: int
    ) -> Optional[Tuple[float, float]]:
        """Get m/z acquisition range for a function. Returns None if unavailable."""
        api = self._ensure_loaded()
        start = api.getAcquisitionRangeStart(handle, function)
        end = api.getAcquisitionRangeEnd(handle, function)
        if start < end and start != DEFAULT_FLOAT:
            return (start, end)
        return None

    def is_raw_spectrum_profile(self, handle, function: int) -> bool:
        """Check if the raw (original) spectrum is profile/continuum data."""
        return bool(self._ensure_loaded().isRawSpectrumContinuum(handle, function) > 0)

    def get_acquisition_date(self, handle) -> str:
        """Get acquisition date string."""
        buf = ctypes.create_string_buffer(64)
        self._ensure_loaded().getAcquisitionDate(handle, buf, 64)
        return buf.value.decode("utf-8", errors="replace")

    def is_lockmass_corrected(self, handle) -> bool:
        """Check if lockmass correction has already been applied."""
        return bool(self._ensure_loaded().isLockmassCorrected(handle) > 0)

    def apply_auto_lockmass_correction(self, handle) -> bool:
        """Apply automatic lockmass correction. Returns True if successful.

        Writes ``lmgt.inf`` into the opened folder when it succeeds, like
        :meth:`apply_custom_lockmass_correction`; never call it on a user's
        run (D83).
        """
        return bool(self._ensure_loaded().applyAutoLockmassCorrection(handle) > 0)

    def apply_custom_lockmass_correction(self, handle, reference_mz: float) -> bool:
        """Correct against ``reference_mz``. Returns True if it was applied.

        The library writes the correction to ``lmgt.inf`` in the folder the
        handle was opened on, and every later open of that folder applies
        it (measured, D83). Never call this on a user's run:
        :func:`thyra.readers.waters.lockmass.measure_lock_mass` calls it on
        a throwaway copy.
        """
        api = self._ensure_loaded()
        return bool(api.applyCustomLockmassCorrection(handle, reference_mz) > 0)

    def get_lockmass_function(self, handle) -> int:
        """Get the lockmass function index, or -1 if none."""
        return int(self._ensure_loaded().getLockmassFunction(handle))
