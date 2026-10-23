# thyra/readers/waters/waters_reader.py
"""Waters .raw MSI reader using MassLynx native libraries.

Reads Waters mass spectrometry imaging data by calling the MassLynxRaw
and MLReader native C libraries via ctypes. The spatial pixel grid is
reconstructed from laser X/Y positions stored in each scan's metadata.
"""

import ctypes
import logging
import math
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Set, Tuple, Union

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from ...core.base_extractor import MetadataExtractor
from ...core.base_reader import BaseMSIReader, summed_peak_mask
from ...core.mass_axis import MassAxisAccumulator
from ...core.mobility import (
    DRIFT_TIME_ACCESSION,
    MILLISECOND_UNIT_ACCESSION,
    MOBILITY_KIND_NAMES,
    MobilityAxis,
)
from ...core.msms import (
    COLLISION_INDUCED_DISSOCIATION_ACCESSION,
    FragmentationSchedule,
    IsolationWindow,
)
from ...core.region_map import (
    build_region_map,
    checked_region_number,
    summarize_regions,
)
from ...core.registry import register_reader
from ...errors import ConversionRefused
from ...metadata.extractors.waters_extractor import WatersMetadataExtractor
from .imaging_grid import ImagingGrid, build_imaging_grid, regrid_for_functions
from .instrument import (
    StagePattern,
    WatersInstrument,
    identify_waters_instrument,
    read_cyclic_sequence,
    read_stage_pattern,
)
from .masslynx_lib import FunctionType, MassLynxLib, ScanInfoData

logger = logging.getLogger(__name__)

#: Bytes per scan record in a ``_FUNCnnn.IDX``: 22 on the Synapt G1 runs
#: measured (``_FUNC003.IDX`` of 28050 bytes for its 1275 scans), 30 on
#: Xevo DESI and Cyclic runs.
_IDX_RECORD_SIZES = (22, 30)

#: What the fragmentation schedule names as its origin.
FRAGMENTATION_SOURCE = "waters_masslynx"

#: How far apart two drift steps may be, relative to the step, and still
#: count as one bin width. The library reports drift times as float32, so
#: steps of one pusher period differ in the seventh digit (0.06625149 vs
#: 0.06625150 ms on a Synapt XS run).
_DRIFT_STEP_TOLERANCE = 1e-4

#: ``MSMSMode`` values of a function that fragments nothing: a full scan,
#: and a quadrupole selection (selected ion monitoring).
_NOT_FRAGMENTED_MODES = frozenset({"MS", "SIM"})


def _equal_step(values: NDArray[np.float64]) -> Optional[float]:
    """The common step of ascending drift times, or ``None`` when they have none.

    Every Waters run measured has one. On a Synapt it is the bin width:
    bins one pusher period apart from 0 ms. On a Cyclic IMS it is the
    library's period, which the run's declared bin timing scales (see
    :meth:`WatersReader._arrival_times`). A run that had none would be
    binned like any continuum.
    """
    if values.size < 2:
        return None
    steps = np.diff(values)
    step: float = float(values[-1] - values[0]) / float(values.size - 1)
    if step <= 0 or not np.all(np.abs(steps - step) <= _DRIFT_STEP_TOLERANCE * step):
        return None
    return step


@register_reader("waters")
class WatersReader(BaseMSIReader):
    """Reader for Waters .raw MSI data using MassLynx native libraries.

    Supports Waters imaging data (.raw directories) containing spatial
    laser coordinates. Uses the same native DLLs as mzmine for data access.

    The reader:
    1. Loads MassLynxRaw + MLReader native libraries via ctypes
    2. Opens the .raw directory and verifies it contains imaging data
    3. Classifies acquisition functions (MS, IMS, MRM, lockmass)
    4. Reconstructs the imaging pixel grid from laser coordinates
    5. Works out which functions hold the image, from the laser positions
    6. Iterates their spectra yielding (coords, mzs, intensities) tuples

    Step 5 exists because a Waters function is not necessarily one
    acquisition function. MassLynx caps a ``_FUNC*.DAT`` file at about
    1.6 GB and opens a new *function* when a long imaging run reaches it,
    so one raster arrives as several functions that tile the stage: they
    never share a pixel. It also names the last of them as the file's
    lockmass function, and reports MS level 1 for the first chunk, 2 for
    the middle ones and 0 for the last. None of that survived contact with
    real files (see :doc:`the D7 entry </design-decisions>`), so the reader
    decides from the laser positions instead, which are the same
    measurement the pixel grid is built from:

    * A function landing on pixels no earlier function covers **extends the
      raster** and is converted, whatever level MassLynx reports for it --
      including the chunk MassLynx calls the lockmass function, except while
      the run is read as centroids, which that one function's spectra cannot
      be (see :meth:`_raster_candidates`).
    * Functions competing for the same pixels are a parallel acquisition --
      MSe low and high energy, a data-dependent run, a co-acquired lockmass
      reference. Summing an MS1 and an MS2 spectrum into one pixel would
      make a spectrum of nothing, so among those the MS1 functions win when
      the file has any, and the rest are recorded under
      ``excluded_functions`` in the Waters-specific metadata block.

    A file whose converted functions carry a precursor m/z holds fragment
    spectra and says so through :meth:`get_fragmentation`, unless the
    function declares an ``MSMSMode`` that fragments nothing (``MS``,
    ``SIM``; D72).

    A run whose method declares regions (a ``MultiRegions`` table) is
    gridded region by region, and ``region=`` converts one of them; see
    :mod:`thyra.readers.waters.imaging_grid`.
    """

    selects_regions = True

    def __init__(
        self,
        data_path: Path,
        use_centroid: Optional[bool] = None,
        intensity_threshold: Optional[float] = None,
        region: Optional[Union[int, str]] = None,
        **kwargs,
    ) -> None:
        """Initialize a Waters MSI reader.

        Args:
            data_path: Path to the Waters .raw directory.
            use_centroid: If True, request vendor centroiding from the DLL.
                If False, read the profile trace. ``None`` (the default)
                decides from the instrument: the profile trace on a SELECT
                SERIES MRT, whose vendor peak picker merges near-isobars the
                trace resolves, and the vendor centroid on every other
                Waters instrument. See :mod:`thyra.readers.waters.instrument`.
            intensity_threshold: Minimum intensity value to include.
            region: Convert only this region of a run whose method declares
                regions (numbered from 1 in the method's table). A run that
                declares none is region 1 as a whole. ``None`` (the
                default) converts everything.
            **kwargs: Additional arguments passed to BaseMSIReader.
                ``metadata_only=True`` also opens a run whose scans all
                report one stage position, which a conversion refuses:
                its metadata is still read, and reading a spectrum is
                refused with the same message.
        """
        super().__init__(data_path, intensity_threshold=intensity_threshold, **kwargs)
        self._metadata_only = bool(kwargs.get("metadata_only", False))

        # Validate the .raw directory structure
        self._validate_raw_directory()

        # Which analyser wrote the run decides the default representation.
        # Read from the side files, so it costs no native-library call and
        # is known before the handle is opened.
        self._instrument: WatersInstrument = identify_waters_instrument(self.data_path)
        # The stage pattern the method declares (pixel steps, regions), also
        # read from a side file. A requested region is checked against it
        # here, before the library opens anything.
        self._pattern: Optional[StagePattern] = read_stage_pattern(self.data_path)
        self._region: Optional[int] = self._checked_region(region)
        if use_centroid is None:
            use_centroid = not self._instrument.is_mrt
            logger.info(
                "%s: reading the %s by default (%s; pass --waters-spectrum "
                "to override)",
                self.data_path.name,
                "vendor centroid" if use_centroid else "profile trace",
                self._instrument.decided_by,
            )
        else:
            logger.info(
                "%s: reading the %s as requested (%s)",
                self.data_path.name,
                "vendor centroid" if use_centroid else "profile trace",
                self._instrument.name,
            )

        # Lazy initialization fields
        self._ml: Optional[MassLynxLib] = None
        self._handle: Optional[ctypes.c_void_p] = None
        self._imaging_grid: Optional[ImagingGrid] = None
        self._function_types: Optional[Dict[int, FunctionType]] = None
        self._ms_functions: Optional[List[int]] = None
        #: Functions not converted: function index -> what they held and
        #: why they stayed out.
        self._excluded_functions: Dict[int, Dict[str, Any]] = {}
        #: Raster chunks MassLynx refuses to centroid: index -> pixels lost.
        self._uncentroidable_functions: Dict[int, int] = {}
        #: Functions covering their own pixels off the raster the MS
        #: functions describe: index -> pixels they held.
        self._off_raster_functions: Dict[int, int] = {}
        self._common_mass_axis_cache: Optional[NDArray[np.float64]] = None
        #: Functions MassLynx reports drift bins for, converted or not.
        self._drift_functions: Set[int] = set()
        #: Drift functions left out because the run is not one the library
        #: reads drift bins from (WatersInstrument.mobility_readable).
        self._unreadable_drift_functions: List[int] = []
        #: The drift axis, read once: ``None`` until asked, then the axis,
        #: or ``False`` when the converted functions give none.
        self._mobility_axis: Any = None
        self._use_centroid = use_centroid
        self._closed = False

    def _checked_region(self, region: object) -> Optional[int]:
        """The ``region=`` request as a region number the run declares.

        Raises:
            ConversionRefused: For a value that is no region number, or a
                number the method does not declare. A run that declares no
                regions is region 1 as a whole.
        """
        try:
            number = checked_region_number(region)
        except ValueError as error:
            raise ConversionRefused(
                f"region takes a region number (a whole number, from 1) or "
                f"None, {error}."
            ) from None
        if number is None:
            return None
        declared = [r.number for r in self._pattern.regions] if self._pattern else []
        if declared and number not in declared:
            raise ConversionRefused(
                f"region={number} is not a region of {self.data_path.name}. "
                f"Its method declares region(s) "
                f"{', '.join(str(n) for n in declared)}."
            )
        if not declared and number != 1:
            raise ConversionRefused(
                f"region={number} is not a region of {self.data_path.name}. "
                f"Its method declares no regions, so the whole run is region 1."
            )
        return number

    def _validate_raw_directory(self) -> None:
        """Validate that data_path is a Waters .raw directory with _FUNC*.DAT files."""
        if not self.data_path.is_dir():
            raise ConversionRefused(
                f"Waters .raw path must be a directory: {self.data_path}"
            )

        # Check for _FUNC*.DAT files (case-insensitive)
        func_files = list(self.data_path.glob("_FUNC[0-9][0-9][0-9].DAT"))
        if not func_files:
            func_files = list(self.data_path.glob("_func[0-9][0-9][0-9].dat"))
        if not func_files:
            # Try broader pattern
            func_files = [
                f
                for f in self.data_path.iterdir()
                if f.name.upper().startswith("_FUNC")
                and f.name.upper().endswith(".DAT")
            ]
        if not func_files:
            raise ConversionRefused(
                f"No _FUNC*.DAT files found in {self.data_path}. "
                "Is this a valid Waters .raw directory?"
            )

    def _function_dat_files(self, suffix: str = ".DAT") -> Dict[int, Path]:
        """The ``_FUNC*<suffix>`` files present, keyed by 0-based function index.

        MassLynx numbers the files from one, so ``_FUNC001.DAT`` holds
        function 0. ``suffix`` is matched without regard to case, as the
        files are written either way.
        """
        found: Dict[int, Path] = {}
        for path in self.data_path.iterdir():
            name = path.name.upper()
            if not (name.startswith("_FUNC") and name.endswith(suffix)):
                continue
            stem = name[len("_FUNC") : -len(suffix)]
            if stem.isdigit():
                found[int(stem) - 1] = path
        return found

    def _validate_function_companions(self) -> None:
        """Refuse a function whose index or status file is gone.

        Each ``_FUNCnnn.DAT`` comes with ``_FUNCnnn.IDX`` (one record per
        scan) and ``_FUNCnnn.STS`` (per-scan status). Without either, or
        with an empty ``.IDX``, the library reports the run as "not a
        Waters imaging file" (measured on a damaged copy of a Cyclic ROI
        export, issue #527). Checked before the file is opened, so the
        refusal names the file that is missing.
        """
        dats = self._function_dat_files()
        indexes = self._function_dat_files(".IDX")
        missing = [
            f"_FUNC{f + 1:03d}{suffix} is missing"
            for suffix in (".IDX", ".STS")
            for f in sorted(set(dats) - set(self._function_dat_files(suffix)))
        ]
        # An empty index beside a data file that holds scans reads as "not
        # imaging" too: the library finds no scans through it.
        missing += [
            f"_FUNC{f + 1:03d}.IDX is empty"
            for f, dat in sorted(dats.items())
            if f in indexes
            and indexes[f].stat().st_size == 0
            and dat.stat().st_size > 0
        ]
        if missing:
            raise ConversionRefused(
                f"{self.data_path.name}: {', '.join(sorted(missing))}. "
                "Each _FUNCnnn.DAT needs its .IDX and .STS beside it; restore "
                "them from the acquisition."
            )

    def _validate_scan_counts(self, n_functions: int) -> None:
        """Refuse a function whose data file holds fewer scans than its index.

        A ``_FUNCnnn.DAT`` cut short or emptied, as an unfinished copy
        leaves it, gave fewer scans and a smaller image at exit 0 (issue
        #527): measured, a 9-scan function cut to 40% returns 3 scans and
        still reads as imaging; cut to nothing it returns 0. The ``.IDX``
        holds one fixed-size record per scan -- 22 bytes on the Synapt G1
        runs measured, 30 on Xevo and Cyclic. A size those records explain
        and a library count that matches none of them is a function that
        lost scans, and refused. A size neither record explains may be a
        MassLynx version not measured: it is warned about, not refused.
        """
        assert self._ml is not None
        indexes = self._function_dat_files(".IDX")
        for f in range(n_functions):
            idx = indexes.get(f)
            if idx is None:
                continue
            size = idx.stat().st_size
            returned = self._ml.get_number_of_scans_in_function(self._handle, f)
            listed = [size // r for r in _IDX_RECORD_SIZES if size % r == 0]
            if not listed:
                logger.warning(
                    "_FUNC%03d.IDX of %s is %d bytes, which no known index "
                    "record size divides; the %d scans the library returned "
                    "for function %d could not be checked against it.",
                    f + 1,
                    self.data_path.name,
                    size,
                    returned,
                    f + 1,
                )
                continue
            if returned in listed:
                continue
            raise ConversionRefused(
                f"{self.data_path.name}: the library returned {returned} scans "
                f"of function {f + 1}, but _FUNC{f + 1:03d}.IDX lists "
                f"{max(listed)}. Is the run complete? Restore "
                f"_FUNC{f + 1:03d}.DAT and _FUNC{f + 1:03d}.IDX from the "
                "acquisition."
            )

    def _validate_function_files(self, n_functions: int) -> None:
        """Refuse a .raw directory that has lost one of its function files.

        The DLL keeps reporting a function whose ``_FUNC*.DAT`` is gone --
        with 0 scans, no error and no warning -- so the only symptom is a
        smaller scan total than the acquisition had. Measured on a 16
        function file with ``_FUNC001.DAT`` deleted: still 16 functions,
        scans [0, 50, 55, ...], 820 positioned scans instead of 848, and
        the reader opened it without a word. On a raster MassLynx split
        across functions (functions 0/1/2 holding rows 0-19 / 19-38 /
        38-45) the same damage drops 20 rows of the image (issue #229).
        """
        present = self._function_dat_files()
        missing = [f for f in range(n_functions) if f not in present]
        if not missing:
            return
        raise ConversionRefused(
            f"{self.data_path.name} declares {n_functions} function(s) but "
            f"{len(missing)} of them have no data file: "
            f"{', '.join(f'_FUNC{f + 1:03d}.DAT' for f in missing)}. "
            f"MassLynx still reports those functions, with 0 scans, so "
            f"converting this directory would silently drop whatever they "
            f"held -- on a raster split across functions, entire rows of "
            f"the image. Restore the missing file(s) from the acquisition."
        )

    def _ensure_initialized(self) -> None:
        """Lazy initialization: load DLL, open file, build imaging grid.

        Called before any data access. Safe to call multiple times.
        """
        if self._handle is not None:
            return

        if self._closed:
            raise RuntimeError("Reader has been closed and cannot be reused")

        # Before the library sees the run: without its .IDX or .STS it
        # answers "not an imaging file" or reads fewer scans.
        self._validate_function_companions()

        # Load the native library (singleton)
        self._ml = MassLynxLib.get_instance(lib_dir=Path(__file__).parent / "lib")

        # Open the .raw directory
        self._handle = self._ml.open_file(str(self.data_path))

        # Verify this is an imaging file. A data file cut to nothing also
        # reads as "not imaging", so the scan counts are checked first: that
        # refusal names the function and its missing scans.
        if not self._ml.is_imaging_file(self._handle):
            try:
                self._validate_scan_counts(
                    self._ml.get_number_of_functions(self._handle)
                )
            except ConversionRefused:
                self._ml.close_file(self._handle)
                self._handle = None
                raise
            self._ml.close_file(self._handle)
            self._handle = None
            raise ConversionRefused(
                f"{self.data_path} is not a Waters imaging file. "
                "The native library reports no imaging data."
            )

        # Set centroid mode
        self._ml.set_centroid(self._handle, self._use_centroid)

        # Classify all functions
        n_funcs = self._ml.get_number_of_functions(self._handle)
        self._function_types = {}
        for f in range(n_funcs):
            ft = self._ml.classify_function(self._handle, f)
            self._function_types[f] = ft
            logger.debug(f"Function {f}: {ft.name}")

        # Which functions carry drift bins, asked per function rather than
        # read off the type: classify_function names the lockmass function
        # before it asks about mobility, and on a mobility raster MassLynx
        # split across functions that is the raster's tail.
        self._drift_functions = {
            f
            for f in range(n_funcs)
            if self._ml.is_ion_mobility_function(self._handle, f)
        }

        # The functions that can hold the image: MS, and MS acquired with
        # ion mobility on, whose spectrum is the one summed over its drift
        # bins. Lockmass, MRM and non-MS functions stay out.
        ms_functions = [
            f
            for f, ft in self._function_types.items()
            if ft in (FunctionType.MS, FunctionType.IMS_MS)
        ]

        # A run the library cannot read drift bins from (no "ADC Mode" in
        # _extern.inf, and not a Cyclic IMS): its mobility functions stay out, as they always did
        # before drift time was read, and the run is refused only when
        # nothing else is left. Decided here, before the grid build, so no
        # spectrum of theirs is ever asked for -- an access violation in
        # the library ends the process with no message.
        self._unreadable_drift_functions = (
            [] if self._instrument.mobility_readable else sorted(self._drift_functions)
        )
        if self._unreadable_drift_functions:
            ms_functions = [f for f in ms_functions if f not in self._drift_functions]
            if not ms_functions:
                self._ml.close_file(self._handle)
                self._handle = None
                raise ConversionRefused(
                    self._unreadable_mobility_message(self._unreadable_drift_functions)
                )
        self._ms_functions = ms_functions

        if not self._ms_functions:
            self._ml.close_file(self._handle)
            self._handle = None
            raise ConversionRefused(
                f"No MS functions found in {self.data_path}. "
                f"Function types: {self._function_types}"
            )

        # Refuse a file whose functions have lost their data. The DLL
        # reports a function whose _FUNC*.DAT is gone, with 0 scans, so
        # nothing downstream notices: on a raster MassLynx split across
        # functions, one missing chunk file silently drops its rows of the
        # image (issue #229).
        try:
            self._validate_function_files(n_funcs)
            self._validate_scan_counts(n_funcs)
        except ConversionRefused:
            self._ml.close_file(self._handle)
            self._handle = None
            raise

        # Build the imaging grid from the MS functions alone. They are what
        # defines the raster; everything else is measured against it. The
        # pitch used to come from every function including the ones
        # selection then excluded, so one reference spot parked off the
        # raster turned a 10x5 image at 30 um into an 11x6 one at
        # 1097 x 440 um (issue #232).
        self._imaging_grid = build_imaging_grid(
            self._ml,
            self._handle,
            self._function_types,
            include_functions=ms_functions,
            allow_single_position=self._metadata_only,
            pattern=self._pattern,
            region=self._region,
        )

        # Which functions hold the image: the raster chunks MassLynx named
        # lockmass belong to it, the parallel functions do not. The first
        # step also records the chunks that had to stay out, so it has to
        # run before the second.
        candidates = self._raster_candidates(
            ms_functions, self._function_types, self._imaging_grid
        )
        self._ms_functions = self._select_converted_functions(
            candidates, self._imaging_grid
        )

        # Re-fit the raster to what is actually converted. A rescued chunk
        # extends the image past the extent of the MS functions, so the
        # grid has to grow to hold it; an excluded function's readings must
        # not leave their pitch behind in a grid it contributes no pixel to.
        self._imaging_grid = regrid_for_functions(
            self._imaging_grid,
            self._ms_functions,
            allow_single_position=self._metadata_only,
            pattern=self._pattern,
            region=self._region,
            fitted_to=ms_functions,
        )

        # Nothing so far has read a spectrum: the grid and the function
        # selection use scan records only. Reading one from a mobility
        # function the library cannot handle ends the process, so this is
        # the last point a refusal can be said.
        self._refuse_unreadable_mobility()

        # A single-pixel grid is refused by build_imaging_grid() above rather
        # than warned about here: "this may not be true imaging data,
        # proceeding anyway" produced a 1x1 store with a 0.0 um pitch and
        # reported success, which is the defect in issue #213. A
        # metadata-only reader gets the grid, and the refusal waits for the
        # first spectrum read (_require_raster).

        logger.info(
            f"Initialized Waters reader: {self.data_path.name}, "
            f"{len(self._ms_functions)} MS functions, "
            f"grid {self._imaging_grid.pixel_count_x}x{self._imaging_grid.pixel_count_y}"
        )

    def _unreadable_mobility_message(self, drift: List[int]) -> str:
        """Why a run's mobility functions are not read (see :meth:`_ensure_initialized`)."""
        return (
            f"{self.data_path.name} was acquired with ion mobility on "
            f"(function(s) {', '.join(str(f) for f in drift)}), and its "
            "_extern.inf does not declare 'Acquisition Algorithm: ADC "
            "Mode', nor is it a Cyclic IMS run ('Acquisition Device: "
            "WatersADC' with a Cyclic-IMS header). On every such run "
            "measured (older Synapt G1 acquisitions) the Waters library "
            "crashed or returned inconsistent data, so Thyra does not read "
            "them."
        )

    def _record_unreadable_drift_functions(self) -> None:
        """List the mobility functions left out of a run the library cannot read."""
        left_out = [
            f
            for f in self._unreadable_drift_functions
            if f not in self._excluded_functions
            and self._function_types is not None
            and self._function_types.get(f)
            in (FunctionType.MS, FunctionType.IMS_MS, FunctionType.LOCKMASS)
        ]
        if not left_out:
            return
        grid = self._imaging_grid
        for f in left_out:
            summary = self._function_summary(f, grid) if grid is not None else {}
            self._excluded_functions[f] = {
                **summary,
                "reason": (
                    "carries drift bins the Waters library cannot read on a run "
                    "without 'Acquisition Algorithm: ADC Mode'"
                ),
            }
        logger.warning(
            "%s Function(s) %s are left out; the others are converted.",
            self._unreadable_mobility_message(left_out),
            ", ".join(str(f) for f in left_out),
        )

    def _refuse_unreadable_mobility(self) -> None:
        """Refuse a run whose converted mobility functions have lost their data.

        A converted mobility function without its ``_FUNCnnn.CDT`` drift
        file is the mobility counterpart of a missing ``.DAT`` (issue
        #229). A run the library cannot read never gets here with a mobility
        function: :meth:`_ensure_initialized` leaves those out before the
        grid is built, and refusing here as well is a safety net.
        """
        self._record_unreadable_drift_functions()
        converted = set(self._ms_functions or ())
        drift = sorted(converted & self._drift_functions)
        if not drift:
            return
        reason = None
        if not self._instrument.mobility_readable:
            reason = self._unreadable_mobility_message(drift)
        else:
            present = {path.name.upper() for path in self.data_path.iterdir()}
            missing = [f for f in drift if f"_FUNC{f + 1:03d}.CDT" not in present]
            if missing:
                reason = (
                    f"{self.data_path.name} was acquired with ion mobility on, "
                    f"but {', '.join(f'_FUNC{f + 1:03d}.CDT' for f in missing)} "
                    "is missing. That file holds the drift bins; restore it "
                    "from the acquisition."
                )
        if reason is None:
            return
        if self._ml is not None and self._handle is not None:
            self._ml.close_file(self._handle)
        self._handle = None
        raise ConversionRefused(reason)

    # ------------------------------------------------------------------
    # Which functions become the stored spectrum
    # ------------------------------------------------------------------

    @staticmethod
    def _positioned_scans(func: int, grid: "ImagingGrid") -> List["ScanInfoData"]:
        """The scan records of one function that carry a laser position.

        Reads the per-function index the grid builds once. This used to
        sort and filter the whole scan map on every call, and function
        selection asks five different questions of each function's scans:
        O(F x N log N) over the file, measured at 87 s of selection against
        a 29 s grid build on a single-function 25 M scan raster, and a
        25-function Synapt G1 run pays it 25 times over (issue #235).
        """
        return grid.scans_by_function.get(func, [])

    @classmethod
    def _function_pixels(
        cls, func: int, grid: "ImagingGrid"
    ) -> Set[Tuple[int, int, int]]:
        """The pixels one function's positioned scans land on."""
        pixels = (
            grid.get_coordinates(info) for info in cls._positioned_scans(func, grid)
        )
        return {c for c in pixels if c is not None}

    @classmethod
    def _comes_back_centroided(cls, func: int, grid: "ImagingGrid") -> bool:
        """Whether MassLynx honours the centroid request for this function.

        ``ScanInfo.isProfile`` is read after ``setCentroid``, so it reports
        what ``getDataPoints`` will hand back rather than how the run was
        acquired. It does not always agree with the request: the library's
        centroider skips the function ``getLockmassFunction`` names, which
        on a chunked imaging file is the last chunk of the raster.
        """
        scans = cls._positioned_scans(func, grid)
        return bool(scans) and not scans[0].is_profile

    def _raster_candidates(
        self,
        ms_functions: List[int],
        function_types: Dict[int, FunctionType],
        grid: "ImagingGrid",
    ) -> List[int]:
        """The MS functions, plus any lockmass function that extends the raster.

        MassLynx names the last function of a chunked imaging file as the
        file's lockmass function even though it is the tail of the raster.
        A function it calls lockmass belongs to the image when it covers
        pixels no MS function covers; a real reference function is acquired
        alongside the image and so covers nothing new. Only ever widens a
        file that already has an MS function, so a run with no MS function
        at all is still refused.

        The catch is that the same library will not centroid that function.
        Rescuing it into a store of centroids would put a band of profile
        rows across the top of the image -- measured at 3.3x the neighbouring
        rows' TIC -- so while the run is being read as centroids it stays
        out, and :attr:`_excluded_functions` records how many pixels that
        costs. Reading the run as the profile trace converts the whole image.
        """
        covered: Set[Tuple[int, int, int]] = set()
        for f in ms_functions:
            covered |= self._function_pixels(f, grid)

        extra, uncentroidable, off_raster = [], [], []
        for f, ft in sorted(function_types.items()):
            if ft is not FunctionType.LOCKMASS:
                continue
            if f in self._unreadable_drift_functions:
                # Its spectra cannot be read safely on this run.
                continue
            # The rescue rule reads "covers pixels no MS function covers"
            # as "is the tail of a split raster". A reference spot parked
            # away from the sample covers a pixel no MS function covers
            # too, and used to be converted as an image pixel on that
            # basis. The tail of a split raster continues the raster, so
            # it lands on the same lattice; a spot does not (issue #232).
            #
            # Tested before the pixel count, not after: an off-raster
            # reading has no pixel in this grid at all, so it would
            # otherwise fall out as "covers nothing new" and never be
            # named.
            scans = self._positioned_scans(f, grid)
            if scans and not all(grid.lies_on_raster(s) for s in scans):
                off_raster.append(f)
                continue
            new_pixels = self._function_pixels(f, grid) - covered
            # Scans in a declared region the grid does not hold yet have no
            # pixel until the regrid, but they are new ground all the same.
            unplaced = sum(1 for s in scans if grid.locate(s) is None)
            if not new_pixels and not unplaced:
                continue
            if self._use_centroid and not self._comes_back_centroided(f, grid):
                uncentroidable.append((f, len(new_pixels) + unplaced))
            else:
                extra.append(f)

        for f in off_raster:
            self._off_raster_functions[f] = len(self._positioned_scans(f, grid))
        if off_raster:
            logger.warning(
                "Function(s) %s cover pixels no MS function covers, but "
                "their stage positions do not lie on the raster the MS "
                "functions describe, so they are not the tail of a split "
                "raster. They are excluded: a reference or calibration "
                "spot converted as an image pixel would stretch the grid "
                "around it.",
                ", ".join(str(f) for f in off_raster),
            )

        if uncentroidable:
            lost = sum(n for _f, n in uncentroidable)
            logger.warning(
                "Function(s) %s hold %d pixels (%.1f%% of the image) that no "
                "other function covers, but MassLynx names them the lockmass "
                "function and will not centroid them. They stay out rather "
                "than put profile rows in a table of centroids: pass "
                "--waters-spectrum profile to convert the whole image.",
                ", ".join(str(f) for f, _n in uncentroidable),
                lost,
                100.0 * lost / max(len(covered) + lost, 1),
            )
        if extra:
            logger.warning(
                "MassLynx names function(s) %s the lockmass function, but they "
                "cover pixels no MS function covers, so they are the tail of a "
                "raster MassLynx split across functions. Converting them too.",
                ", ".join(str(f) for f in extra),
            )
        self._uncentroidable_functions = dict(uncentroidable)
        return sorted(set(ms_functions) | set(extra))

    @classmethod
    def _function_ms_level(cls, func: int, grid: "ImagingGrid") -> int:
        """The MS level MassLynx reports for a function's scans.

        Read off the first positioned scan: MassLynx sets the level per
        function, so one scan speaks for all of them. A function without a
        positioned scan contributes no pixel and is treated as MS1 so that
        it is never the reason an MS1 image is refused.
        """
        scans = cls._positioned_scans(func, grid)
        return max(1, int(scans[0].ms_level)) if scans else 1

    @classmethod
    def _pixel_groups(
        cls, ms_functions: List[int], grid: "ImagingGrid"
    ) -> List[List[int]]:
        """Group the functions that compete for the same pixels.

        Functions in one group were acquired in parallel and only one of
        them can be the pixel's spectrum. Separate groups tile the stage --
        they are the chunks MassLynx makes when a raster outgrows one
        ``_FUNC*.DAT`` file -- and every one of them is part of the image.
        """
        groups: List[Tuple[Set[Tuple[int, int, int]], List[int]]] = []
        for f in ms_functions:
            pixels = cls._function_pixels(f, grid)
            for covered, members in groups:
                if covered & pixels:
                    covered |= pixels
                    members.append(f)
                    break
            else:
                groups.append((pixels, [f]))
        return [members for _covered, members in groups]

    def _select_converted_functions(
        self, ms_functions: List[int], grid: "ImagingGrid"
    ) -> List[int]:
        """Keep every raster chunk, and one MS level per pixel within each.

        Chunks of one raster are all converted. Where functions do compete
        for a pixel, the MS1 ones win when the group has any; the rest are
        recorded in :attr:`_excluded_functions` so the store can say the
        acquisition fragmented something even though the stored spectra are
        intact-ion spectra.
        """
        self._excluded_functions = {
            f: {
                **self._function_summary(f, grid),
                "n_unique_pixels": n_pixels,
                "reason": "MassLynx will not centroid this function",
            }
            for f, n_pixels in self._uncentroidable_functions.items()
        }
        # No n_unique_pixels: an off-raster reading has no pixel in this
        # grid, which is the point. _function_summary's n_scans says how
        # much was left out.
        self._excluded_functions.update(
            {
                f: {
                    **self._function_summary(f, grid),
                    "reason": "stage positions do not lie on the imaging raster",
                }
                for f in self._off_raster_functions
            }
        )
        kept: List[int] = []
        for group in self._pixel_groups(ms_functions, grid):
            kept.extend(self._select_within_group(group, grid))
        return sorted(kept)

    def _function_summary(self, func: int, grid: "ImagingGrid") -> Dict[str, Any]:
        """What a function that was not converted held."""
        precursors = self._function_precursors(func, grid)
        return {
            "ms_level": self._function_ms_level(func, grid),
            "precursor_mz": (
                precursors[0]
                if len(precursors) == 1
                else (precursors if precursors else None)
            ),
            "n_scans": len(self._positioned_scans(func, grid)),
        }

    def _select_within_group(self, group: List[int], grid: "ImagingGrid") -> List[int]:
        """Which of a set of functions competing for one pixel set is stored."""
        group = self._without_parallel_drift_functions(group, grid)
        levels = {f: self._function_ms_level(f, grid) for f in group}
        ms1 = [f for f in group if levels[f] == 1]
        if not ms1 or len(ms1) == len(group):
            if len(group) > 1:
                logger.warning(
                    "%d MS functions (%s) cover the same pixels; their spectra "
                    "are summed per pixel.",
                    len(group),
                    ", ".join(str(f) for f in group),
                )
            return list(group)

        for f in group:
            if levels[f] == 1:
                continue
            self._excluded_functions[f] = {
                **self._function_summary(f, grid),
                "reason": "covers the same pixels as an MS1 function",
            }
        excluded = [f for f in group if levels[f] != 1]
        logger.warning(
            "Functions %s are MS level %s and cover the same pixels as the "
            "MS1 function(s) %s. Only the MS1 spectra are converted; the "
            "others are recorded in the Waters metadata block as "
            "'excluded_functions'.",
            ", ".join(str(f) for f in excluded),
            "/".join(sorted({str(levels[f]) for f in excluded})),
            ", ".join(str(f) for f in ms1),
        )
        return ms1

    def _without_parallel_drift_functions(
        self, group: List[int], grid: "ImagingGrid"
    ) -> List[int]:
        """Leave out the functions with drift bins when others share their pixels.

        Before drift time was read, a function acquired with mobility on was
        never a candidate, so a run recording one in parallel with plain MS
        functions stored the plain ones alone. Summing the two now would
        change those stored spectra, and their drift bins would not add up
        to the sum anyway, so the old store is kept and the mobility
        function is listed under ``excluded_functions``.
        """
        drift = [f for f in group if f in self._drift_functions]
        if not drift or len(drift) == len(group):
            return group
        for f in drift:
            self._excluded_functions[f] = {
                **self._function_summary(f, grid),
                "reason": (
                    "carries drift bins and covers the same pixels as a "
                    "function that does not"
                ),
            }
        logger.warning(
            "Function(s) %s carry drift bins and cover the same pixels as "
            "%s, which do not. Only the latter are converted, as before "
            "drift time was read; the former are recorded in the Waters "
            "metadata block as 'excluded_functions'.",
            ", ".join(str(f) for f in drift),
            ", ".join(str(f) for f in group if f not in drift),
        )
        return [f for f in group if f not in drift]

    def _function_precursors(self, func: int, grid: "ImagingGrid") -> List[float]:
        """The distinct precursor m/z values a function's scans report, ascending.

        Empty for a function whose ``_extern.inf`` section declares
        ``MSMSMode`` as ``MS`` or ``SIM``: nothing was fragmented, whatever
        the scans report. A Cyclic IMS full-scan run reports 500.0 on
        every scan, which is its ``Quad.MSFixedMass`` setting, and an
        HDSIM run reports the quadrupole set mass (D72).
        """
        if self._instrument.msms_mode(func) in _NOT_FRAGMENTED_MODES:
            return []
        return self._reported_precursors(func, grid)

    @classmethod
    def _reported_precursors(cls, func: int, grid: "ImagingGrid") -> List[float]:
        """The distinct positive precursor m/z the scans report, whatever the mode."""
        values = {
            round(float(info.precursor_mz), 4)
            for info in cls._positioned_scans(func, grid)
            # Finite as well as positive: an infinite target is no m/z
            # (issue #498).
            if math.isfinite(info.precursor_mz) and info.precursor_mz > 0
        }
        return sorted(values)

    def function_settings(self) -> Dict[str, Dict[str, Any]]:
        """What each converted function's ``_extern.inf`` section declares.

        Keyed by function index as a string, only for functions that declare
        anything: ``msms_mode``, the quadrupole ``sim_set_mass_mz`` of a
        ``SIM`` function (the m/z its scans report, since nothing else
        carries it), and ``ccs_calibration`` as written (``No`` on every
        Cyclic run measured). Synapt runs declare none of these.
        """
        _ml, _handle, grid, _types, ms_functions = self._require_initialized()
        settings: Dict[str, Dict[str, Any]] = {}
        for func in ms_functions:
            entry: Dict[str, Any] = {}
            mode = self._instrument.msms_mode(func)
            if mode is not None:
                entry["msms_mode"] = mode
            if mode == "SIM":
                set_mass = self._reported_precursors(func, grid)
                if len(set_mass) == 1:
                    entry["sim_set_mass_mz"] = set_mass[0]
                elif set_mass:
                    entry["sim_set_mass_mz"] = set_mass
            ccs = self._instrument.function_parameters.get(func, {}).get(
                "CCS Calibration"
            )
            if ccs:
                entry["ccs_calibration"] = ccs
            if entry:
                settings[str(func)] = entry
        return settings

    @property
    def has_fragmentation(self) -> bool:
        """Always True: MassLynx reports the precursor of every function.

        ``has_precursor_spectra`` stays False -- Waters stores one spectrum
        per function per position, so there is nothing within a pixel to
        separate; the functions were already separated at acquisition.
        """
        return True

    def get_fragmentation(self) -> Optional[FragmentationSchedule]:
        """What the stored spectra are: MS1, or the fragment spectra of what.

        A converted function holds fragment spectra when MassLynx reports a
        precursor m/z for it. The reported MS level alone is not evidence:
        on every real multi-function imaging file measured, the chunks of
        one raster come back as level 1, then 2, then 0, with no precursor
        anywhere (see the D7 entry of the decisions page). So an MSe or
        data-dependent file, whose MS/MS functions were left out, reads as
        MS1 here -- what they were is in the Waters metadata block, since
        this describes the spectra in the store -- and so does a chunked
        MS1 raster. A file whose converted functions do carry precursors
        reports them: one isolation window per distinct precursor, and
        ``constant_across_pixels`` false as soon as a function's precursor
        changes from scan to scan. A precursor is not evidence either when
        the function declares ``MSMSMode`` ``MS`` or ``SIM``: a Cyclic IMS
        reports one on every full scan (:meth:`_function_precursors`).
        """
        _ml, _handle, grid, _types, ms_functions = self._require_initialized()
        precursors = {f: self._function_precursors(f, grid) for f in ms_functions}
        fragment_functions = [f for f in ms_functions if precursors[f]]
        if not fragment_functions:
            return FragmentationSchedule(ms_level=1, source=FRAGMENTATION_SOURCE)

        windows: Dict[float, IsolationWindow] = {}
        constant = True
        for f in fragment_functions:
            if len(precursors[f]) != 1:
                constant = False
                continue
            target = precursors[f][0]
            windows.setdefault(target, self._isolation_window(f, target, grid))
        return FragmentationSchedule(
            ms_level=max(
                [2] + [self._function_ms_level(f, grid) for f in fragment_functions]
            ),
            windows=tuple(windows[target] for target in sorted(windows)),
            constant_across_pixels=constant,
            dissociation_accession=(
                COLLISION_INDUCED_DISSOCIATION_ACCESSION if windows else None
            ),
            source=FRAGMENTATION_SOURCE,
        )

    @classmethod
    def _isolation_window(
        cls, func: int, target: float, grid: "ImagingGrid"
    ) -> IsolationWindow:
        """One function's isolation window from its first positioned scan."""
        info = cls._positioned_scans(func, grid)[0]
        lower = upper = None
        start, end = float(info.quad_isolation_start), float(info.quad_isolation_end)
        if 0.0 < start <= target <= end and math.isfinite(end):
            lower, upper = target - start, end - target
        energy = float(info.collision_energy)
        return IsolationWindow(
            target=target,
            lower_offset=lower,
            upper_offset=upper,
            collision_energy=(energy if math.isfinite(energy) and energy > 0 else None),
        )

    def _require_initialized(
        self,
    ) -> Tuple[
        "MassLynxLib",
        ctypes.c_void_p,
        "ImagingGrid",
        Dict[int, "FunctionType"],
        List[int],
    ]:
        """Ensure the reader is initialized and return required fields."""
        self._ensure_initialized()
        if (
            self._ml is None
            or self._handle is None
            or self._imaging_grid is None
            or self._function_types is None
            or self._ms_functions is None
        ):
            raise RuntimeError("Reader not fully initialized")
        return (
            self._ml,
            self._handle,
            self._imaging_grid,
            self._function_types,
            self._ms_functions,
        )

    @property
    def no_raster_reason(self) -> Optional[str]:
        """Why the run cannot be converted, when its scans share one position.

        ``None`` for a raster, and before the run is opened. Only a
        metadata-only reader opens such a run at all; :func:`thyra.preview_msi`
        reads this to report it as unreadable, as it did before metadata
        could be read from it (issue #494).
        """
        grid = self._imaging_grid
        return grid.no_raster_reason if grid is not None else None

    def _require_raster(self, grid: "ImagingGrid") -> None:
        """Refuse to read spectra from a run that has no raster.

        Only a metadata-only reader gets that far: everywhere else the
        grid build refuses the run (issue #213).
        """
        if grid.no_raster_reason is not None:
            raise ConversionRefused(grid.no_raster_reason)

    @property
    def has_shared_mass_axis(self) -> bool:
        """Waters MSI data is typically processed/centroided with varying m/z per pixel."""
        return False

    @property
    def instrument(self) -> WatersInstrument:
        """What the run's own metadata says about the analyser."""
        return self._instrument

    @property
    def use_centroid(self) -> bool:
        """Whether spectra come from the vendor peak picker (else the profile trace)."""
        return self._use_centroid

    def _create_metadata_extractor(self) -> MetadataExtractor:
        """Create Waters metadata extractor."""
        ml, handle, imaging_grid, function_types, ms_functions = (
            self._require_initialized()
        )
        return WatersMetadataExtractor(
            ml=ml,
            handle=handle,
            data_path=self.data_path,
            imaging_grid=imaging_grid,
            function_types=function_types,
            ms_functions=ms_functions,
            instrument=self._instrument,
            use_centroid=self._use_centroid,
            excluded_functions=self._excluded_functions,
            mobility_report=self.mobility_report(),
            function_settings=self.function_settings(),
            cyclic_sequence=(
                read_cyclic_sequence(self.data_path)
                if self._instrument.is_cyclic_ims
                else None
            ),
            stage_pattern=self._pattern,
        )

    def _placed_scans(self) -> Generator[Tuple[int, int, int], None, None]:
        """``(x, y, region)`` of every placed scan, in acquisition order.

        Only meaningful on a grid that uses regions. A scan of a function
        outside the converted set, or of a region the grid does not hold,
        is left out, as the spectrum iteration leaves it out.
        """
        _ml, _handle, grid, _types, ms_functions = self._require_initialized()
        for func in ms_functions:
            for info in grid.scans_by_function.get(func, []):
                coords = grid.locate(info)
                region = grid.region_of(info)
                if coords is not None and region is not None:
                    yield coords[0], coords[1], region

    def get_region_map(self) -> Optional[Dict[Tuple[int, int], int]]:
        """Map 0-based ``(x, y)`` to the declared region, or ``None``.

        ``None`` when the run's method declares no regions or its scans do
        not lie in them. Under ``region=`` only that region is mapped.
        """
        _ml, _handle, grid, _types, _ms = self._require_initialized()
        if not grid.regions:
            return None
        return build_region_map(self._placed_scans()) or None

    def get_region_info(self) -> Optional[list]:
        """Region number, spectrum count and box of each declared region.

        The same keys and box frame as the timsTOF reader, without a name:
        the method numbers its regions and names none. ``None`` when the
        run uses no regions.
        """
        _ml, _handle, grid, _types, _ms = self._require_initialized()
        if not grid.regions:
            return None
        return summarize_regions(self._placed_scans()) or None

    def get_common_mass_axis(self) -> NDArray[np.float64]:
        """Build common mass axis from all unique m/z values across all spectra.

        Since Waters MSI data is typically processed/centroided, each spectrum
        can have different m/z values. This iterates all MS spectra to collect
        all unique values and returns them sorted.

        Returns:
            Sorted array of unique m/z values across all spectra.
        """
        if self._common_mass_axis_cache is not None:
            return self._common_mass_axis_cache

        ml, handle, imaging_grid, function_types, ms_functions = (
            self._require_initialized()
        )
        self._require_raster(imaging_grid)

        total = sum(ml.get_number_of_scans_in_function(handle, f) for f in ms_functions)

        # Folded rather than collected. This held every scan's m/z array in
        # a list and then paid np.concatenate plus np.unique, so the
        # transient was the whole peak payload rather than its distinct
        # values -- three payload-sized allocations live at once. Worst
        # exactly where Waters usually sits: an MRT defaults to a profile
        # trace, where every scan shares the same grid, so the unique count
        # is small and the concatenate is pure waste (#294).
        #
        # ``max_mass_axis_length`` is honoured if the caller passes one
        # (BaseMSIReader takes it for every reader now), but it has no
        # default here -- only imzML sets one. So no acquisition that
        # converts today starts being refused, and turning the cap on for
        # this reader stays a defaults decision rather than a side effect
        # of removing a duplicated builder.
        accumulator = MassAxisAccumulator(total, max_length=self.max_mass_axis_length)
        n_spectra = 0

        with tqdm(total=total, desc="Building mass axis", unit="scan") as pbar:
            for func in ms_functions:
                n_scans = ml.get_number_of_scans_in_function(handle, func)
                for scan in range(n_scans):
                    pbar.update(1)

                    scan_info = imaging_grid.scan_map.get((func, scan))
                    if scan_info is None or not scan_info.has_position:
                        continue
                    # A scan the grid gives no pixel (another region, under
                    # region=) adds nothing to the store, nor to its axis.
                    if imaging_grid.locate(scan_info) is None:
                        continue

                    try:
                        mzs, _ = ml.read_spectrum(handle, func, scan)
                    except Exception as e:
                        logger.debug(
                            f"Error reading spectrum func={func} scan={scan}: {e}"
                        )
                        continue

                    # Outside the read's try/except, deliberately: a
                    # max_mass_axis_length refusal is a ConversionRefused,
                    # which subclasses ValueError, so folding inside that
                    # handler would swallow it and log it at DEBUG as one
                    # unreadable scan. Same reasoning as the imzML build.
                    if mzs.size > 0:
                        n_spectra += 1
                        accumulator.add(mzs)
                    mzs = None

        self._common_mass_axis_cache = accumulator.finish()

        logger.info(
            f"Built common mass axis with {len(self._common_mass_axis_cache):,} "
            f"unique m/z values from {n_spectra:,} spectra"
        )
        return self._common_mass_axis_cache

    def iter_spectra(self) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Iterate through all MS imaging spectra with spatial coordinates.

        Only yields spectra from the converted functions (MS, and MS with
        ion mobility on, read summed over its drift bins). Skips scans
        without valid laser positions.

        Yields:
            Tuple of ((x, y, z), mzs, intensities) where coordinates are
            0-based pixel indices.
        """
        for coords, _, mzs, intensities in self._iter_spectra_with_order():
            yield coords, mzs, intensities

    @property
    def has_acquisition_order(self) -> bool:
        """Always: every scan of a converted function is numbered."""
        return True

    def iter_spectra_with_acquisition_order(self) -> Generator[
        Tuple[
            Tuple[int, int, int],
            int,
            NDArray[np.float64],
            NDArray[np.float64],
        ],
        None,
        None,
    ]:
        """:meth:`iter_spectra` with each scan's number in the run.

        MassLynx numbers scans from 0 within each function, and a long
        raster continues in the next function once one fills (design
        decision D7), so the scans are numbered on through the converted
        functions in order: function ``f``'s scan ``s`` is ``s`` plus the
        scan counts of the converted functions before ``f``. A scan with no
        stage position keeps its number, so the numbers can have gaps. On
        three real runs, one split across three functions, the scans'
        retention time never decreases along these numbers (D26).
        """
        yield from self._iter_spectra_with_order()

    def _iter_spectra_with_order(
        self,
    ) -> Generator[
        Tuple[
            Tuple[int, int, int],
            int,
            NDArray[np.float64],
            NDArray[np.float64],
        ],
        None,
        None,
    ]:
        """The scan loop behind both spectrum iterators, with run-wide numbers."""
        ml, handle, imaging_grid, function_types, ms_functions = (
            self._require_initialized()
        )
        self._require_raster(imaging_grid)

        total = sum(ml.get_number_of_scans_in_function(handle, f) for f in ms_functions)

        with tqdm(
            total=total,
            desc="Reading spectra",
            unit="spectrum",
            disable=getattr(self, "_quiet_mode", False),
        ) as pbar:
            first_of_function = 0
            for func in ms_functions:
                n_scans = ml.get_number_of_scans_in_function(handle, func)
                for scan in range(n_scans):
                    pbar.update(1)

                    # Look up cached scan info from the grid build pass
                    scan_info = imaging_grid.scan_map.get((func, scan))
                    if scan_info is None or not scan_info.has_position:
                        continue

                    # Map laser position to pixel coordinates
                    coords = imaging_grid.get_coordinates(scan_info)
                    if coords is None:
                        continue

                    try:
                        mzs, intensities = ml.read_spectrum(handle, func, scan)

                        # Apply intensity threshold filtering (from BaseMSIReader)
                        mzs, intensities = self._apply_intensity_filter(
                            mzs, intensities
                        )

                        if mzs.size > 0 and intensities.size > 0:
                            yield coords, first_of_function + scan, mzs, intensities
                    except Exception as e:
                        logger.warning(
                            f"Error reading spectrum func={func} scan={scan}: {e}"
                        )
                first_of_function += n_scans

    # ------------------------------------------------------------------
    # Ion mobility
    #
    # A Waters mobility scan has the shape of a Bruker TDF frame: one
    # pixel, a fixed list of drift bins shared by every scan, a different
    # point cloud in each bin. The bins are equal steps of drift time: from
    # 0 ms one pusher period apart on a Synapt (200 bins on every run
    # measured), and from the declared start delay a declared number of
    # pushes apart on a Cyclic IMS (D72).
    # iter_spectra keeps yielding the spectrum summed over the bins.
    # ------------------------------------------------------------------

    @property
    def has_ion_mobility(self) -> bool:
        """True when every converted function carries drift bins.

        A run that converts a function without drift bins next to one
        with them -- two functions summed into the same pixels -- has no
        point cloud that adds up to its stored spectrum, so it reports no
        mobility rather than half of it.
        """
        return self.get_mobility_axis() is not None

    def get_mobility_axis(self) -> Optional[MobilityAxis]:
        """The drift-time axis, in ms, or ``None`` when there is none.

        ``values[b]`` is the drift time of bin ``b``, ascending: the
        library's value on a Synapt, placed by the declared bin timing on a
        Cyclic IMS (:meth:`_arrival_times`). Read once and cached.
        """
        if self._mobility_axis is None:
            axis = self._build_mobility_axis()
            self._mobility_axis = axis if axis is not None else False
        cached = self._mobility_axis
        return cached if isinstance(cached, MobilityAxis) else None

    def _build_mobility_axis(self) -> Optional[MobilityAxis]:
        ml, handle, _grid, _types, ms_functions = self._require_initialized()
        drift = [f for f in ms_functions if f in self._drift_functions]
        if not drift:
            return None
        if len(drift) != len(ms_functions):
            logger.warning(
                "Function(s) %s carry drift bins and %s do not, and their "
                "spectra are summed into the same image. No mobility is "
                "recorded: the drift bins would not add up to the stored "
                "spectrum.",
                ", ".join(str(f) for f in drift),
                ", ".join(str(f) for f in ms_functions if f not in drift),
            )
            return None
        counts = {f: ml.get_number_of_drift_bins(handle, f) for f in drift}
        if len(set(counts.values())) != 1 or not counts[drift[0]]:
            # getMobility takes no function argument, so functions with
            # different bin counts cannot all be described by its axis.
            logger.warning(
                "The mobility functions report different drift bin counts "
                "(%s); no mobility is recorded.",
                ", ".join(f"{f}: {n}" for f, n in counts.items()),
            )
            return None
        values = self._arrival_times(
            drift, ml.get_drift_times(handle, counts[drift[0]])
        )
        if values is None:
            return None
        return MobilityAxis(
            kind_accession=DRIFT_TIME_ACCESSION,
            kind_name=MOBILITY_KIND_NAMES[DRIFT_TIME_ACCESSION],
            unit_accession=MILLISECOND_UNIT_ACCESSION,
            unit_name="millisecond",
            values=values,
            source="waters_raw",
            bin_width=_equal_step(values),
        )

    def _arrival_times(
        self, drift: List[int], library_times: NDArray[np.float64]
    ) -> Optional[NDArray[np.float64]]:
        """Each bin's drift time, from the library's axis and the run's bin timing.

        ``getMobility`` numbers the bins one pusher period apart from 0
        ms. That is the drift time on a Synapt, whose functions declare no
        bin timing. A Cyclic IMS function declares where its bins start
        and how many pushes each spans (:meth:`WatersInstrument.
        drift_bin_timing`), and bin ``i`` then arrives at ``(delay + i *
        pushes)`` periods, the period being the library's own step. On
        the Cyclic runs measured this matches the drift times MassLynx
        writes into its own drift-time export to float32 precision
        (309.002 to 348.403 ms in 0.198 ms steps, D72).

        Returns:
            The drift times, or ``None`` (logged) when the functions
            declare different timings, a timing that is not two whole
            numbers of pushes, or the library's steps are uneven, since
            none of those can be placed on one axis.
        """
        try:
            timings = {f: self._instrument.drift_bin_timing(f) for f in drift}
        except ValueError as error:
            logger.warning(
                "%s: %s; no mobility is recorded.", self.data_path.name, error
            )
            return None
        if all(timing is None for timing in timings.values()):
            return library_times
        if len(set(timings.values())) != 1:
            logger.warning(
                "The mobility functions declare different drift bin timings "
                "(%s); no mobility is recorded.",
                ", ".join(f"{f}: {t}" for f, t in timings.items()),
            )
            return None
        timing = timings[drift[0]]
        period = _equal_step(library_times)
        if timing is None or period is None:
            logger.warning(
                "The Waters library's drift bins are not one pusher period "
                "apart, so the declared bin timing cannot place them; no "
                "mobility is recorded."
            )
            return None
        delay, pushes = timing
        bins = np.arange(library_times.size, dtype=np.float64)
        return (delay + pushes * bins) * period

    @property
    def lossless_spectrum_option(self) -> Optional[str]:
        """The option that makes the stored spectrum the drift bins' sum, when needed.

        The bins are read raw, so they sum to the profile spectrum exactly
        (measured) and to the vendor centroid not at all: centroiding the
        sum is not summing the centroids. ``None`` when the profile trace
        is what is stored.
        """
        return "--waters-spectrum profile" if self._use_centroid else None

    def mobility_report(self) -> Optional[Dict[str, Any]]:
        """What the metadata extractor reports under ``ion_mobility``."""
        axis = self.get_mobility_axis()
        if axis is None:
            return None
        report = axis.to_extractor_report()
        if axis.values is not None:
            report["num_scans"] = int(axis.values.size)
        drift = [f for f in self._ms_functions or () if f in self._drift_functions]
        timing = self._instrument.drift_bin_timing(drift[0]) if drift else None
        if timing is not None:
            # How the axis was placed (_arrival_times), in the run's units.
            report["start_delay_pushes"], report["pushes_per_bin"] = timing
        return report

    def iter_mobility_spectra(self) -> Generator[
        Tuple[
            Tuple[int, int, int],
            NDArray[np.float64],
            NDArray[np.float64],
            NDArray[np.float64],
        ],
        None,
        None,
    ]:
        """Each scan's drift bins as one ``(coords, mz, drift_ms, intensity)`` cloud.

        The same scans, coordinates and order as :meth:`iter_spectra`.
        Each bin is read raw (:meth:`MassLynxLib.read_drift_bin`), so the
        points summed over the bins are the profile spectrum of the scan,
        whichever representation the summed table was read as. The
        intensity threshold, when set, keeps or drops a peak's points
        across all bins together (D34). No progress bar: the caller's pass
        draws one.

        Raises:
            NotImplementedError: When the run carries no drift bins.
        """
        axis = self.get_mobility_axis()
        if axis is None or axis.values is None:
            raise NotImplementedError(
                f"{self.data_path.name} carries no readable drift bins"
            )
        ml, handle, imaging_grid, _types, ms_functions = self._require_initialized()
        self._require_raster(imaging_grid)
        drift_times = axis.values
        n_bins = int(drift_times.size)
        for func in ms_functions:
            n_scans = ml.get_number_of_scans_in_function(handle, func)
            for scan in range(n_scans):
                scan_info = imaging_grid.scan_map.get((func, scan))
                if scan_info is None or not scan_info.has_position:
                    continue
                coords = imaging_grid.get_coordinates(scan_info)
                if coords is None:
                    continue
                # Never past the function's own count: an out-of-range
                # bin index is the kind of call this library answers
                # with an exception it raises in C++.
                bins = min(n_bins, max(0, int(scan_info.drift_scan_count)))
                mz_parts, mobility_parts, intensity_parts = [], [], []
                for b in range(bins):
                    mzs, intensities = ml.read_drift_bin(handle, func, scan, b)
                    if mzs.size == 0:
                        continue
                    mz_parts.append(mzs)
                    intensity_parts.append(intensities)
                    mobility_parts.append(np.full(mzs.size, drift_times[b]))
                if not mz_parts:
                    continue
                mzs = np.concatenate(mz_parts)
                mobility = np.concatenate(mobility_parts)
                intensities = np.concatenate(intensity_parts)
                if self._intensity_threshold is not None:
                    keep = summed_peak_mask(mzs, intensities, self._intensity_threshold)
                    mzs, mobility, intensities = (
                        mzs[keep],
                        mobility[keep],
                        intensities[keep],
                    )
                if mzs.size > 0:
                    yield coords, mzs, mobility, intensities

    def close(self) -> None:
        """Close the file handle and release native resources."""
        if self._closed:
            return
        if self._handle is not None and self._ml is not None:
            try:
                self._ml.close_file(self._handle)
                logger.debug(f"Closed Waters file: {self.data_path}")
            except Exception as e:
                logger.error(f"Error closing Waters file: {e}")
            self._handle = None
        self._closed = True

    def reset(self) -> None:
        """Reset reader for re-iteration.

        Stateless iteration (no internal cursor), so this is a no-op.
        The scan_map cache persists across iterations for efficiency.
        """
        pass

    def __repr__(self) -> str:
        """Return string representation of the reader."""
        grid_info = ""
        if self._imaging_grid:
            grid_info = (
                f", grid={self._imaging_grid.pixel_count_x}"
                f"x{self._imaging_grid.pixel_count_y}"
            )
        return f"WatersReader(path={self.data_path}{grid_info}, centroid={self._use_centroid})"
