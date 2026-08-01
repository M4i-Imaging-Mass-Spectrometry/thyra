# thyra/readers/imzml/imzml_reader.py
import logging
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Tuple, Union, cast

import numpy as np
from numpy.typing import NDArray
from pyimzml.ImzMLParser import ImzMLParser
from tqdm import tqdm

from ...core.base_extractor import MetadataExtractor
from ...core.base_reader import BaseMSIReader
from ...core.registry import register_reader
from ...metadata.extractors.imzml_extractor import ImzMLMetadataExtractor

logger = logging.getLogger(__name__)

# Scratch buffer floor for the processed-mode mass axis build, in ELEMENTS.
# 1Mi elements is 8 MiB at float64. Measured at a 512 MiB m/z payload on a
# shared grid (peak extra RSS / wall clock): 2^14 -> 7.14 MiB / 3.04 s,
# 2^17 -> 6.90 / 3.06, 2^20 -> 15.30 / 2.30, 2^22 -> 40.29 / 2.11,
# 2^24 -> 148.38 / 2.06. A 1024x sweep of the knob moves wall clock 1.5x with
# no cliff anywhere, so this is a knee rather than a fragile constant.
_MASS_AXIS_BATCH_VALUES = 1 << 20

# Block size for the searchsorted probe in ``_filter_absent``, in ELEMENTS.
# Bounds that function's temporaries to O(block) rather than O(batch).
_MASS_AXIS_PROBE_BLOCK = 1 << 20


def _dedupe_sorted(a: NDArray[Any]) -> NDArray[Any]:
    """Deduplicate an ALREADY-SORTED array, as ``np.unique`` would.

    Reproduces numpy's ``_unique1d`` mask including its ``equal_nan=True``
    behaviour, which collapses a trailing run of NaNs to its first element.
    Unlike ``np.unique`` it does not make the extra full-size ``flatten()``
    copy, which is one of the three payload-sized allocations the old
    collect-everything-then-unique build paid for.

    Args:
        a: A sorted array.

    Returns:
        The distinct values of ``a``, in order.
    """
    n = a.size
    if n == 0:
        return a[:0].copy()

    mask = np.empty(n, dtype=bool)
    mask[0] = True

    if a.dtype.kind == "f" and bool(np.isnan(a[-1])):
        first_nan = int(np.searchsorted(a, a[-1], side="left"))
        # The ``first_nan > 0`` guard is load-bearing: an all-NaN input makes
        # the slices empty and ``np.not_equal`` raises on the shape mismatch.
        if first_nan > 0:
            np.not_equal(a[1:first_nan], a[: first_nan - 1], out=mask[1:first_nan])
        mask[first_nan] = True
        mask[first_nan + 1 :] = False
    else:
        np.not_equal(a[1:], a[:-1], out=mask[1:])

    return a[mask]


def _filter_absent(
    acc: NDArray[Any],
    s: NDArray[Any],
    block: int = _MASS_AXIS_PROBE_BLOCK,
) -> NDArray[Any]:
    """Return the elements of ``s`` that are not already in ``acc``.

    Both arrays must be sorted and free of duplicates. ``s`` is walked in
    blocks so the search temporaries stay O(block) rather than O(len(s)).
    Returns ``s`` itself, uncopied, when none of its values are present --
    the all-distinct fast path, worth half a payload at peak.

    Args:
        acc: The running axis, sorted and unique.
        s: A candidate batch, sorted and unique.
        block: Probe block size in elements.

    Returns:
        The subset of ``s`` absent from ``acc``.
    """
    n = s.size
    # ``acc.size == 0`` must be handled here: the ``np.minimum`` below would
    # otherwise index ``acc[-1]`` on an empty array and raise IndexError.
    if n == 0 or acc.size == 0:
        return s

    keep = np.empty(n, dtype=bool)
    n_keep = 0
    is_float = acc.dtype.kind == "f"
    a_len = acc.size

    for lo in range(0, n, block):
        hi = min(lo + block, n)
        sb = s[lo:hi]
        pos = np.searchsorted(acc, sb, side="left")
        out_of_range = pos >= a_len
        np.minimum(pos, a_len - 1, out=pos)
        va = acc[pos]
        eq = va == sb
        if is_float:
            # searchsorted routes a NaN key onto acc's first NaN, but NaN does
            # not compare equal to itself, so without this the same NaN would
            # be re-inserted on every fold.
            np.logical_or(eq, np.isnan(va) & np.isnan(sb), out=eq)
        np.logical_and(eq, ~out_of_range, out=eq)
        np.logical_not(eq, out=eq)
        keep[lo:hi] = eq
        n_keep += int(np.count_nonzero(eq))
        del pos, out_of_range, va, eq, sb

    if n_keep == n:
        return s
    if n_keep == 0:
        return s[:0]
    return s[keep]


def _merge_disjoint(box_a: List[Any], box_b: List[Any]) -> NDArray[Any]:
    """Merge two DISJOINT sorted arrays, releasing each once it is copied.

    The inputs arrive boxed in one-element lists so this function can drop the
    caller's reference as soon as each side has been copied. Passing them as
    plain parameters would keep both alive for the whole call and add a full
    copy of the larger side to the peak.

    The result needs no second deduplication because ``_filter_absent`` has
    already removed the overlap, which saves another full-size copy.

    Args:
        box_a: One-element list holding the first sorted array.
        box_b: One-element list holding the second sorted array.

    Returns:
        The sorted union of the two inputs.
    """
    a, b = box_a[0], box_b[0]
    out = np.empty(a.size + b.size, dtype=a.dtype)

    n_a = a.size
    out[:n_a] = a
    box_a[0] = None
    del a

    out[n_a:] = b
    box_b[0] = None
    del b

    # Exactly two ascending runs, which is the case quicksort handles best.
    out.sort(kind="quicksort")
    return out


def _read_spectrum_mzs(parser: Any, idx: int) -> Optional[NDArray[Any]]:
    """Read one spectrum's m/z values, or None if missing or unreadable.

    Args:
        parser: An initialized ImzML parser.
        idx: Spectrum index.

    Returns:
        The m/z array, or None when the spectrum is empty or failed to read.
    """
    try:
        spectrum_data = parser.getspectrum(idx)
    except Exception as e:
        logger.warning(f"Error getting spectrum {idx}: {e}")
        return None

    if spectrum_data is None or len(spectrum_data) < 1:
        return None
    mzs = spectrum_data[0]
    return mzs if mzs.size else None


class _MassAxisAccumulator:
    """Builds a sorted, unique m/z axis without holding every spectrum.

    Values are copied into a scratch buffer; when that buffer fills, it is
    sorted, deduplicated, and merged into the running axis. The buffer's
    capacity tracks the axis length, so the number of merges is logarithmic in
    the input rather than linear -- a fixed batch size would re-copy the
    accumulator once per batch, which is the quadratic behaviour that makes a
    naive ``np.union1d`` fold far slower than the code it replaces.

    Peak memory therefore depends on the number of *unique* m/z values rather
    than on the size of the file.
    """

    def __init__(self, total_spectra: int, max_length: Optional[int] = None) -> None:
        """Initialize the accumulator.

        Args:
            total_spectra: Spectrum count, used only in the error message.
            max_length: Cap on unique m/z values, or None for unlimited.
        """
        self._total_spectra = total_spectra
        self._max_length = max_length
        self._acc: Optional[NDArray[Any]] = None
        self._buf: Optional[NDArray[Any]] = None
        self._cap = 0
        self._cap_target = _MASS_AXIS_BATCH_VALUES
        self._n = 0
        self.saw_any = False

    def add(self, mzs: NDArray[Any], index: int) -> None:
        """Buffer one spectrum's m/z values, folding first if the buffer is full.

        Args:
            mzs: The spectrum's m/z values. Not modified.
            index: Spectrum index, used only in the error message.
        """
        m = mzs.size
        if not m:
            return
        self.saw_any = True

        if self._buf is None:
            self._allocate(m, mzs.dtype)
        elif self._n + m > self._cap:
            self._fold(index)
            if self._buf is None:
                self._allocate(m, mzs.dtype)
            elif m > self._cap:
                self._buf = None
                self._allocate(m, mzs.dtype, exact=True)

        assert self._buf is not None
        self._buf[self._n : self._n + m] = mzs
        self._n += m

    def finish(self, index: int) -> NDArray[Any]:
        """Fold whatever is buffered and return the axis.

        Args:
            index: Last spectrum index, used only in the error message.

        Returns:
            Sorted, unique m/z values.

        Raises:
            ValueError: If no spectrum yielded any m/z values.
        """
        self._fold(index)
        self._buf = None

        if not self.saw_any or self._acc is None:
            raise ValueError("No spectra found to build common mass axis")
        if self._acc.size == 0:
            raise ValueError("Failed to extract any m/z values")
        return self._acc

    def _allocate(self, m: int, dtype: Any, exact: bool = False) -> None:
        """Allocate the scratch buffer, at least large enough for one spectrum."""
        self._cap = m if exact else max(self._cap_target, m)
        self._buf = np.empty(self._cap, dtype=dtype)

    def _fold(self, index: int) -> None:
        """Sort, deduplicate and merge the buffered batch into the axis."""
        if self._n == 0 or self._buf is None:
            return

        view = self._buf[: self._n]
        view.sort(kind="quicksort")
        run = _dedupe_sorted(view)
        self._n = 0

        if self._cap > _MASS_AXIS_BATCH_VALUES:
            # The scratch has grown to axis size, so release it before the
            # merge allocates. Gating on the batch floor keeps this from firing
            # on every one of the thousands of folds a shared-grid dataset
            # performs, where re-faulting the buffer each time costs more than
            # the memory it frees.
            view = None
            self._buf = None
            self._cap = 0

        if self._acc is None:
            self._acc = run
        else:
            new = _filter_absent(self._acc, run)
            run = None  # drop the batch before the merge allocates
            if new.size:
                box_a, box_b = [self._acc], [new]
                self._acc = None
                new = None
                self._acc = _merge_disjoint(box_a, box_b)

        self._check_limit(index)

        # Record the next capacity but do NOT allocate it here. On the final
        # fold the stream is already exhausted, so allocating would commit a
        # whole axis-sized buffer that is never written.
        self._cap_target = max(_MASS_AXIS_BATCH_VALUES, self._acc.size)
        if self._cap and self._cap != self._cap_target:
            self._buf = None
            self._cap = 0

    def _check_limit(self, index: int) -> None:
        """Stop if the axis has outgrown ``max_length``."""
        if self._max_length is None or self._acc is None:
            return
        if self._acc.size <= self._max_length:
            return
        raise ValueError(
            f"Common mass axis exceeded {self._max_length:,} unique m/z values "
            f"after {index + 1:,} of {self._total_spectra:,} spectra "
            f"({self._acc.size:,} so far). The peak lists in this dataset do "
            "not share m/z values, so a raw axis grows to roughly one column "
            "per peak, which is not usable downstream. Convert with resampling "
            "instead (it is the default; --no-resample disables it), or raise "
            "max_mass_axis_length."
        )


@register_reader("imzml")
class ImzMLReader(BaseMSIReader):
    """Reader for imzML format files with optimizations for performance."""

    def __init__(
        self,
        data_path: Path,
        batch_size: int = 50,
        cache_coordinates: bool = True,
        **kwargs,
    ) -> None:
        """Initialize an ImzML reader.

        Args:
            data_path: Path to the imzML file
            batch_size: Default batch size for spectrum iteration
            cache_coordinates: Whether to cache coordinates upfront
            **kwargs: Additional arguments. ``max_mass_axis_length`` caps the
                number of unique m/z values the processed-mode raw axis may
                reach before the build gives up; None (the default) is
                unlimited. See ``_extract_continuous_mass_axis``.
        """
        super().__init__(data_path, **kwargs)
        self.filepath: Optional[Union[str, Path]] = data_path
        self.batch_size: int = batch_size
        self.cache_coordinates: bool = cache_coordinates
        self.max_mass_axis_length: Optional[int] = kwargs.get("max_mass_axis_length")
        self.parser: Optional[ImzMLParser] = None
        self.ibd_file: Optional[Any] = None
        self.imzml_path: Optional[Path] = None
        self.ibd_path: Optional[Path] = None
        self.is_continuous: bool = False
        self.is_processed: bool = False

        # Parser initialization flag for lazy loading
        self._parser_initialized: bool = False

        # Cached properties
        self._common_mass_axis: Optional[NDArray[np.float64]] = None
        self._coordinates_array: Optional[NDArray[np.int32]] = (
            None  # Fast numpy array cache
        )

        # Store path but don't initialize parser yet - wait for first use
        if data_path is not None:
            self.filepath = data_path

    def _ensure_parser_initialized(self) -> None:
        """Guarantee parser is initialized exactly once."""
        if not self._parser_initialized:
            if self.filepath is None:
                raise ValueError("No file path provided for parser initialization")
            self._initialize_parser(self.filepath)
            self._parser_initialized = True

    def _initialize_parser(self, imzml_path: Union[str, Path]) -> None:
        """Initialize the ImzML parser with the given path.

        Args:
            imzml_path: Path to the imzML file to parse

        Raises:
            ValueError: If the corresponding .ibd file is not found or metadata
                parsing fails
            Exception: If parser initialization fails
        """
        if isinstance(imzml_path, str):
            imzml_path = Path(imzml_path)

        self.imzml_path = imzml_path
        self.ibd_path = imzml_path.with_suffix(".ibd")

        if not self.ibd_path.exists():
            raise ValueError(f"Corresponding .ibd file not found for {imzml_path}")

        # Open the .ibd file for reading
        self.ibd_file = open(self.ibd_path, mode="rb")

        # Initialize the parser
        logger.info(f"Initializing ImzML parser for {imzml_path}")
        try:
            # ElementTree, NOT lxml. pyimzml prunes each <spectrum> out of the
            # tree (`slist.remove(elem)`) while iterparse is still streaming the
            # document, which invalidates libxml2's text-node coalescing
            # accelerator: ctxt->nodelen/nodemem go on describing a text node
            # inside the subtree that was just removed, so later character data
            # is appended at a stale offset and the buffer is doubled on every
            # miss until xmlRealloc fails. lxml surfaces libxml2's
            # XML_ERR_NO_MEMORY as a *syntax* error, so it reads as a corrupt
            # file when nothing is wrong with it:
            #     XMLSyntaxError: xmlSAX2Characters, line 212575, column 1
            # It only bites when the text between </spectrum> and <spectrum> is
            # exactly "\r\n" -- CRLF with no indentation, how IONTOF SurfaceLab
            # writes imzML. Indented or LF-only files coalesce differently and
            # survive, which is why most files never hit it. ElementTree builds
            # the tree in Python, so there is no C parser state to invalidate;
            # it is also pyimzml's own default, parses byte-identically, and is
            # faster here -- 67s vs 118s on a 2.1 GB imzML at the same peak RSS.
            self.parser = ImzMLParser(
                filename=str(imzml_path),
                parse_lib="ElementTree",
                ibd_file=self.ibd_file,
            )
        except Exception as e:
            if self.ibd_file:
                self.ibd_file.close()
            logger.error(f"Failed to initialize ImzML parser: {e}")
            raise

        if self.parser.metadata is None:
            raise ValueError("Failed to parse metadata from imzML file.")

        # Determine file mode
        # Determine if file is continuous mode
        self.is_continuous = (
            "continuous" in self.parser.metadata.file_description.param_by_name
        )
        # Determine if file is processed mode
        self.is_processed = (
            "processed" in self.parser.metadata.file_description.param_by_name
        )

        if self.is_continuous == self.is_processed:
            raise ValueError(
                "Invalid file mode, expected either 'continuous' or " "'processed'."
            )

        # Cache coordinates if requested
        if self.cache_coordinates:
            self._cache_all_coordinates()

    def _cache_all_coordinates(self) -> None:
        """Cache all coordinates for faster access.

        Converts 1-based coordinates from imzML to 0-based coordinates
        for internal use. Uses vectorized numpy operations for speed.
        Stores as numpy array for O(1) index lookup without dict overhead.
        """
        # Parser should already be initialized when this is called from
        # _initialize_parser

        if self.parser is None:
            raise RuntimeError("Parser is not initialized")
        n_coords = len(self.parser.coordinates)
        logger.info(f"Caching {n_coords:,} coordinates...")

        # Vectorized conversion using numpy (much faster than Python loop)
        # np.array() on the coordinates list is the main cost here
        self._coordinates_array = np.array(self.parser.coordinates, dtype=np.int32)

        # Convert to 0-based in place (subtract 1, but z minimum is 0)
        self._coordinates_array[:, :2] -= 1  # x and y
        self._coordinates_array[:, 2] = np.maximum(
            self._coordinates_array[:, 2] - 1, 0
        )  # z

        logger.info(f"Cached {n_coords:,} coordinates as numpy array")

    def _create_metadata_extractor(self) -> MetadataExtractor:
        """Create ImzML metadata extractor."""
        self._ensure_parser_initialized()

        if not self.imzml_path:
            raise ValueError("ImzML path not available")

        return ImzMLMetadataExtractor(self.parser, self.imzml_path)

    @property
    def has_shared_mass_axis(self) -> bool:
        """Check if all spectra share the same m/z axis.

        Returns True for continuous ImzML (all pixels have same m/z values),
        False for processed ImzML (each pixel has different m/z values).
        """
        return self.is_continuous

    def get_common_mass_axis(self) -> NDArray[np.float64]:
        """Return the common mass axis composed of all unique m/z values.

        For continuous mode, returns the m/z values from the first spectrum.
        For processed mode, collects all unique m/z values across spectra.

        Returns:
            NDArray[np.float64]: Array of m/z values in ascending order

        Raises:
            ValueError: If the common mass axis cannot be created
        """
        self._ensure_parser_initialized()

        if self._common_mass_axis is None:
            # We know parser is not None at this point
            parser = cast(ImzMLParser, self.parser)

            if self.is_continuous:
                logger.info("Using m/z values from first spectrum (continuous mode)")
                spectrum_data = parser.getspectrum(0)
                if spectrum_data is None or len(spectrum_data) < 1:
                    raise ValueError("Could not get first spectrum")

                mzs = spectrum_data[0]
                if mzs.size == 0:
                    raise ValueError("First spectrum contains no m/z values")

                self._common_mass_axis = mzs
            else:
                self._common_mass_axis = self._extract_continuous_mass_axis(parser)

        # Return the common mass axis
        return self._common_mass_axis

    def _extract_continuous_mass_axis(self, parser: ImzMLParser) -> NDArray[np.float64]:
        """Build the common mass axis for processed-mode data.

        Returns exactly what ``np.unique(np.concatenate(all_mzs))`` returned:
        sorted, deduplicated, with the dtype following the file's mzPrecision.

        This streams the file rather than collecting it. The previous
        implementation held every spectrum's m/z array in a list, concatenated
        that into one array, and handed the result to ``np.unique`` -- three
        live copies of the whole m/z payload P, plus the output. Measured peak
        was 3.40-3.51x P on data whose spectra share m/z values and 4.32-4.39x
        P when every value is distinct, which put a hard ceiling around a
        40 GB input on a 128 GB machine.

        Here each spectrum is copied into a scratch buffer, and when that
        buffer fills it is sorted, deduplicated, and merged into the running
        axis. Peak memory becomes a function of the number of *unique* values
        rather than the payload: measured 15.2-15.5 MiB, flat, for payloads
        from 32 MiB to 16 GiB of shared-grid data, and 2.00x P in the
        all-distinct case, which is the structural floor for an out-of-place
        merge.

        Args:
            parser: An initialized ImzML parser.

        Returns:
            Sorted, unique m/z values across all spectra.

        Raises:
            ValueError: If no spectrum yielded any m/z values, or if
                ``max_mass_axis_length`` is set and the axis outgrows it.
        """
        logger.info(
            "Building common mass axis from all unique m/z values " "(processed mode)"
        )

        total_spectra = len(parser.coordinates)
        accumulator = _MassAxisAccumulator(
            total_spectra, getattr(self, "max_mass_axis_length", None)
        )
        idx = -1

        with tqdm(
            total=total_spectra,
            desc="Building common mass axis",
            unit="spectrum",
        ) as pbar:
            for idx in range(total_spectra):
                mzs = _read_spectrum_mzs(parser, idx)
                if mzs is not None:
                    # Deliberately outside the read's try/except: a
                    # max_mass_axis_length failure must not be swallowed and
                    # logged as a per-spectrum warning.
                    accumulator.add(mzs, idx)
                mzs = None
                pbar.update(1)

        axis = accumulator.finish(idx)
        logger.info(f"Created common mass axis with {axis.size} unique m/z values")
        return axis

    def _get_spectrum_coordinates(
        self, parser: ImzMLParser, idx: int
    ) -> Tuple[int, int, int]:
        """Get 0-based coordinates for a spectrum."""
        if self._coordinates_array is not None:
            # Fast O(1) numpy array lookup
            row = self._coordinates_array[idx]
            return (int(row[0]), int(row[1]), int(row[2]))

        # Fallback: compute on the fly
        x, y, z = parser.coordinates[idx]
        return cast(
            Tuple[int, int, int],
            (x - 1, y - 1, z - 1 if z > 0 else 0),
        )

    def _process_single_spectrum(
        self, parser: ImzMLParser, idx: int, pbar
    ) -> Optional[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]]
    ]:
        """Process a single spectrum and return its data."""
        try:
            coords = self._get_spectrum_coordinates(parser, idx)
            mzs, intensities = parser.getspectrum(idx)

            # Apply intensity threshold filtering if configured
            mzs, intensities = self._apply_intensity_filter(mzs, intensities)

            if mzs.size > 0 and intensities.size > 0:
                pbar.update(1)
                return coords, mzs, intensities

            pbar.update(1)
            return None
        except Exception as e:
            logger.warning(f"Error processing spectrum {idx}: {e}")
            pbar.update(1)
            return None

    def _iter_spectra_single(
        self, parser: ImzMLParser, total_spectra: int, pbar
    ) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Process spectra one at a time."""
        for idx in range(total_spectra):
            result = self._process_single_spectrum(parser, idx, pbar)
            if result is not None:
                yield result

    def _iter_spectra_batch(
        self, parser: ImzMLParser, total_spectra: int, batch_size: int, pbar
    ) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Process spectra in batches."""
        for batch_start in range(0, total_spectra, batch_size):
            batch_end = min(batch_start + batch_size, total_spectra)
            batch_size_actual = batch_end - batch_start

            for offset in range(batch_size_actual):
                idx = batch_start + offset
                result = self._process_single_spectrum(parser, idx, pbar)
                if result is not None:
                    yield result

    def iter_spectra(self, batch_size: Optional[int] = None) -> Generator[
        Tuple[Tuple[int, int, int], NDArray[np.float64], NDArray[np.float64]],
        None,
        None,
    ]:
        """Iterate through spectra with progress monitoring and batch processing.

        Maps m/z values to the common mass axis using searchsorted for
        accurate representation in the output data structures.

        Args:
            batch_size: Number of spectra to process in each batch (None for
                default)

        Yields:
            Tuple containing:
                - Tuple[int, int, int]: Coordinates (x, y, z) - 0-based
                - NDArray[np.float64]: m/z values array
                - NDArray[np.float64]: Intensity values array

        Raises:
            ValueError: If parser is not initialized and no filepath is
                available
        """
        self._ensure_parser_initialized()

        if batch_size is None:
            batch_size = self.batch_size

        parser = cast(ImzMLParser, self.parser)
        total_spectra = len(parser.coordinates)
        dimensions = self.get_essential_metadata().dimensions
        total_pixels = dimensions[0] * dimensions[1] * dimensions[2]

        logger.info(
            f"Processing {total_spectra} spectra in a grid of " f"{total_pixels} pixels"
        )

        with tqdm(
            total=total_spectra,
            desc="Reading spectra",
            unit="spectrum",
            disable=getattr(self, "_quiet_mode", False),
        ) as pbar:
            if batch_size <= 1:
                yield from self._iter_spectra_single(parser, total_spectra, pbar)
            else:
                yield from self._iter_spectra_batch(
                    parser, total_spectra, batch_size, pbar
                )

    def read(self) -> Dict[str, Any]:
        """Read the entire imzML file and return a structured data dictionary.

        Returns:
            Dict containing:
                - mzs: NDArray[np.float64] - common m/z values array
                - intensities: NDArray[np.float64] - array of intensity arrays
                - coordinates: List[Tuple[int, int, int]] - list of (x,y,z)
                  coordinates
                - width: int - number of pixels in x dimension
                - height: int - number of pixels in y dimension
                - depth: int - number of pixels in z dimension

        Raises:
            ValueError: If parser is not initialized and no filepath is
                available
        """
        self._ensure_parser_initialized()

        # Get common mass axis
        mzs = self.get_common_mass_axis()

        # Get dimensions
        width, height, depth = self.get_essential_metadata().dimensions

        # Collect all spectra
        coordinates: List[Tuple[int, int, int]] = []
        intensities: List[NDArray[np.float64]] = []

        # Iterate through all spectra
        for coords, spectrum_mzs, spectrum_intensities in self.iter_spectra():
            coordinates.append(coords)

            # Convert sparse representation to full array
            full_spectrum = np.zeros(len(mzs), dtype=np.float64)

            # Find indices in the common mass axis using searchsorted
            indices = np.searchsorted(mzs, spectrum_mzs)

            # Ensure indices are within bounds
            valid_indices = indices < len(mzs)
            indices = indices[valid_indices]
            valid_intensities = spectrum_intensities[valid_indices]

            # Fill spectrum array
            full_spectrum[indices] = valid_intensities
            intensities.append(full_spectrum)

        return {
            "mzs": mzs,
            "intensities": np.array(intensities, dtype=np.float64),
            "coordinates": coordinates,
            "width": width,
            "height": height,
            "depth": depth,
        }

    def close(self) -> None:
        """Close all open file handles."""
        if hasattr(self, "ibd_file") and self.ibd_file is not None:
            self.ibd_file.close()
            self.ibd_file = None

        if hasattr(self, "parser") and self.parser is not None:
            if hasattr(self.parser, "m") and self.parser.m is not None:
                self.parser.m.close()
            self.parser = None

    @property
    def n_spectra(self) -> int:
        """Return the total number of spectra in the dataset.

        Returns:
            Total number of spectra (efficient implementation using parser)
        """
        self._ensure_parser_initialized()

        # Use parser coordinates which is efficient
        parser = cast(ImzMLParser, self.parser)
        return len(parser.coordinates)

    def get_total_peak_count(self) -> int:
        """Get total number of peaks across all spectra.

        For ImzML, this requires iterating through spectra to count peaks.

        Returns:
            Total number of peaks across all spectra
        """
        self._ensure_parser_initialized()
        parser = cast(ImzMLParser, self.parser)
        total_spectra = len(parser.coordinates)

        logger.info("Counting peaks across all spectra for exact allocation...")
        total_peaks = 0

        with tqdm(
            total=total_spectra,
            desc="Counting peaks",
            unit="spectrum",
        ) as pbar:
            for idx in range(total_spectra):
                try:
                    mzs, _ = parser.getspectrum(idx)
                    total_peaks += len(mzs)
                except Exception as e:
                    logger.warning(f"Error getting spectrum {idx}: {e}")
                pbar.update(1)

        logger.info(f"Total peak count: {total_peaks:,}")
        return total_peaks

    @property
    def mass_range(self) -> Tuple[float, float]:
        """Return the mass range (min_mz, max_mz) of the dataset.

        Returns:
            Tuple of (min_mz, max_mz) values
        """
        # Get mass range from essential metadata
        essential_metadata = self.get_essential_metadata()
        return essential_metadata.mass_range

    def get_peak_counts_per_pixel(self) -> Optional[NDArray[np.int32]]:
        """Get per-pixel peak counts for CSR indptr construction.

        Returns peak counts collected during metadata extraction.
        This enables optimized streaming conversion without a separate
        counting pass.

        Returns:
            Array of size n_pixels where arr[pixel_idx] = peak_count.
            pixel_idx = z * (n_x * n_y) + y * n_x + x
            Returns None if not available.
        """
        essential_metadata = self.get_essential_metadata()
        return essential_metadata.peak_counts_per_pixel
