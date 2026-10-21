"""The Python entry point: one call, source path to store.

:func:`convert_msi` is what the package is for, and the only function most
callers need. It detects the format, builds the reader, chooses and builds
the converter, runs it, and cleans up after it. The CLI in
:mod:`thyra.__main__` is a thin wrapper over this, so the two cannot
diverge in behaviour -- only in how arguments are spelled.

Most of this module is the part that runs before any data is read.
Arguments are validated up front, against the paths as well as each other,
because a conversion that fails on a bad ``dataset_id`` after an hour of
reading has wasted the hour: the identifier becomes an element name in the
store, so it is checked for the characters a store will accept rather than
discovered to be unusable at write time.

**A failed conversion must not leave something that looks finished.**
:func:`quarantine_partial_output` renames a partial store to ``.failed``,
which is what ``docs/cli.md`` promises of any failed run and what frees
the output path for a retry. This is the counterpart to ``convert()``
catching ``KeyboardInterrupt``: both exist so that the only thing at the
output path is either a complete store or something plainly marked as not
one (issues #245, #293).
"""

import logging
import numbers
import traceback
import warnings
from math import isfinite
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union
from uuid import uuid4

import numpy as np

from .core.base_converter import (
    BaseMSIConverter,
    PixelSizeSource,
    is_usable_pitch,
    log_interrupt_after_save,
)
from .core.base_reader import BaseMSIReader
from .core.registry import detect_format, get_converter_class, get_reader_class
from .errors import ConversionRefused
from .metadata.constants import sample_value_problem
from .reader_options import READER_OPTION_FORMATS, format_name, reader_option_reaches
from .resampling.plan import MassAxisPlan
from .resampling.types import ResamplingConfig
from .utils.windows_paths import prepare_zarr_output_path

logger = logging.getLogger(__name__)

warnings.filterwarnings(
    "ignore",
    message=r"Accession IMS:1000046.*",  # ignore UserWarning
    category=UserWarning,
    module=r"pyimzml.ontology.ontology",
)


def _validate_paths_parameters(
    input_path: Union[str, Path], output_path: Union[str, Path]
) -> bool:
    """Validate path parameters."""
    if not input_path or not isinstance(input_path, (str, Path)):
        logger.error("Input path must be a valid string or Path object")
        return False

    if not output_path or not isinstance(output_path, (str, Path)):
        logger.error("Output path must be a valid string or Path object")
        return False

    return True


#: Characters SpatialData allows in an element name, besides letters and
#: digits. Its rule (``spatialdata._core.validation.check_valid_name``) is
#: not part of the public API, so it is restated here rather than imported:
#: a private import that moves would take the front-door check down with
#: it, and the only cost of the two drifting is that SpatialData refuses a
#: name Thyra let through -- which is exactly today's behaviour.
_ELEMENT_NAME_EXTRA_CHARS = "_-."

#: Coordinate systems every element is written to besides the per-dataset
#: one, which the dataset id names. An id equal to one of them would merge
#: the two. ``"global"`` merging is harmless (both carry the same
#: transform); ``"acquisition"`` is not, so it is refused (D48). Restated
#: from ``converters.spatialdata.acquisition_frame`` rather than imported:
#: importing that package here would load SpatialData with this module.
_RESERVED_DATASET_IDS = ("acquisition",)

#: The longest suffix any route appends to the dataset id to name an
#: element: an optical image, ``<id>_optical_`` plus a stem cut to 30
#: characters (``optical_image.NAME_BUDGET``). Restated rather than
#: imported for the reason given above ``_RESERVED_DATASET_IDS``.
_LONGEST_ELEMENT_SUFFIX = len("_optical_") + 30

#: The most bytes one file or folder name may hold: NTFS and ext4 both stop
#: at 255. Every element is a folder named after it.
_MAX_NAME_BYTES = 255

#: Names Windows keeps for devices, with or without an extension. Under
#: ``handle_3d`` the table folder is the bare id, and ``nul`` there failed
#: at the save after both passes (issue #519). Refused on every platform,
#: so a store written anywhere can be copied to Windows.
_WINDOWS_DEVICE_NAMES = frozenset(
    {"con", "prn", "aux", "nul"}
    | {f"{port}{n}" for port in ("com", "lpt") for n in "123456789¹²³"}
)


def dataset_id_problem(dataset_id: Any) -> Optional[str]:
    """Why ``dataset_id`` cannot name SpatialData elements, or ``None``.

    Every element key in the store is built from this id -- ``<id>_z0``,
    ``<id>_z0_tic``, ``<id>_pixels`` -- so SpatialData's naming rule
    applies to it. It used to be checked for emptiness only, and an id
    with a space or a slash was refused by SpatialData at ``_save_output``,
    after both passes over the source had run (issue #250). The suffixes
    Thyra appends are all legal characters, so an id that passes here
    yields keys that pass there.

    Three more rules come from what reads the store and what holds it on
    disk rather than from SpatialData's writer (issue #519): a leading
    dot, which SpatialData's reader skips as a hidden file, so the store
    opened empty; an id too long for the longest element name built from
    it; and a Windows device name.
    """
    if not isinstance(dataset_id, str):
        return f"Dataset ID must be a string, not {type(dataset_id).__name__}"
    if not dataset_id.strip():
        return "Dataset ID must be a non-empty string"
    if dataset_id in (".", ".."):
        return f"Dataset ID cannot be {dataset_id!r}"
    if dataset_id.startswith("."):
        return (
            f"Dataset ID {dataset_id!r} cannot start with '.': SpatialData "
            "skips such names when it reads a store, so the store would open "
            "empty"
        )
    if dataset_id.startswith("__"):
        return "Dataset ID cannot start with '__'"
    if dataset_id in _RESERVED_DATASET_IDS:
        return (
            f"Dataset ID cannot be {dataset_id!r}: every store writes a "
            "coordinate system of that name"
        )
    if dataset_id.split(".", 1)[0].lower() in _WINDOWS_DEVICE_NAMES:
        return (
            f"Dataset ID cannot be {dataset_id!r}: Windows keeps that name "
            "for a device, so no folder can carry it"
        )
    longest = len(dataset_id.encode("utf-8")) + _LONGEST_ELEMENT_SUFFIX
    if longest > _MAX_NAME_BYTES:
        return (
            f"Dataset ID is too long: the longest element name built from it "
            f"would be {longest} bytes, and a file name holds at most "
            f"{_MAX_NAME_BYTES}. Use at most "
            f"{_MAX_NAME_BYTES - _LONGEST_ELEMENT_SUFFIX} bytes (one per "
            "letter or digit in plain ASCII)"
        )
    bad = sorted(
        {c for c in dataset_id if not (c.isalnum() or c in _ELEMENT_NAME_EXTRA_CHARS)}
    )
    if bad:
        return (
            "Dataset ID names every element in the store, so it may hold only "
            "letters, digits, underscores, dots and hyphens; "
            f"{dataset_id!r} also holds " + ", ".join(repr(c) for c in bad)
        )
    return None


def _validate_sample_parameters(sample_id: Any, sample_name: Any) -> bool:
    """Validate the user-supplied sample id and name."""
    for name, value in (("sample_id", sample_id), ("sample_name", sample_name)):
        problem = sample_value_problem(name, value)
        if problem is not None:
            logger.error(problem)
            return False
    return True


def _validate_string_parameters(format_type: str, dataset_id: str) -> bool:
    """Validate string parameters."""
    if not isinstance(format_type, str) or not format_type.strip():
        logger.error("Format type must be a non-empty string")
        return False

    problem = dataset_id_problem(dataset_id)
    if problem is not None:
        logger.error(problem)
        return False

    return True


def length_problem(value: Any, label: str) -> Optional[str]:
    """Why ``value`` is not a length in micrometres, or ``None``.

    Any real number counts, numpy scalars included: ``np.float32(25)``
    and ``np.int64(25)`` used to be refused as "not a finite positive
    number" while ``np.float64(25)`` passed, because only the last is a
    ``float`` (issue #519). A bool does not count, though Python calls it
    an int: ``True`` became a 1 um z spacing.

    ``<= 0`` is False for NaN and for +infinity, so both used to pass and
    reach the store as the dataset's pixel size or z spacing (issue #231).
    """
    # numpy's bool is no numbers.Real, so only Python's needs naming.
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        return f"{label} must be a number, not {type(value).__name__}"
    if not isfinite(value) or value <= 0:
        return f"{label} must be a finite positive number, got {value!r}"
    return None


def _validate_numeric_parameters(
    pixel_size_um: Optional[float], z_spacing_um: Optional[float] = None
) -> bool:
    """Validate numeric parameters.

    The API is validated as well as the CLI because a caller can hand
    these in directly (issue #231).
    """
    for value, label in ((pixel_size_um, "Pixel size"), (z_spacing_um, "Z spacing")):
        problem = None if value is None else length_problem(value, label)
        if problem is not None:
            logger.error(problem)
            return False
    return True


def _validate_handle_3d(handle_3d: Any) -> bool:
    """Accept only a bool, numpy's included.

    Any truthy value used to build a volume: ``handle_3d="no"`` did
    (issue #519).
    """
    if type(handle_3d) in (bool, np.bool_):
        return True
    logger.error(
        "handle_3d must be True or False, not %s %r",
        type(handle_3d).__name__,
        handle_3d,
    )
    return False


def _validate_streaming(streaming: object) -> bool:
    """Validate the ``streaming`` argument, which selects nothing.

    A no-op argument still has to say when it was misspelled. ``"yes"``,
    ``None`` and ``1`` were all accepted in silence, and so was ``0`` --
    which a caller writes meaning ``False``, and which the ``is False``
    test in ``_create_converter`` misses, because ``0`` and ``False`` are
    equal but not identical. The one value the code answers was the one a
    plausible spelling could hide (issue #261).

    Checked here rather than in ``_create_converter`` so it fails while the
    caller is still looking at their own arguments, before the source is
    opened -- the same reason ``spectrum_type`` is normalised in
    ``ImzMLReader.__init__``.
    """
    if streaming is True or streaming is False or streaming == "auto":
        return True

    logger.error(
        "streaming must be True, False or 'auto', got %r. It selects nothing "
        "either way -- every conversion streams since v3.23 -- so the "
        "argument can simply be dropped.",
        streaming,
    )
    return False


#: The converter keywords ``convert_msi`` documents under ``**kwargs``.
#: ``use_csc`` is kept for older callers; it selects nothing.
CONVERTER_KEYWORDS = (
    "write_mobility_table",
    "mobility_heatmap",
    "mobility_grid",
    "mobility_bins",
    "mobility_min",
    "mobility_max",
    "msms_table",
    "use_csc",
)

#: Keywords that were removed, with what a caller should do instead.
_RETIRED_KEYWORDS = {
    "sparse_format": (
        "sparse_format was removed: every converter writes CSC, which is the "
        "layout an ion image reads down. Drop the argument; for row-wise "
        "access call X.tocsr() on the matrix you read back"
    ),
}


def keyword_problems(kwargs: Dict[str, Any], reader_options: Any) -> List[str]:
    """Each keyword ``convert_msi`` does not use, and why.

    A keyword the converter did not use was ignored but written into the
    store's ``conversion_options`` as if applied, could overwrite a value
    Thyra computed, and failed the save after both passes when JSON could
    not hold it. A misspelt reader option changed nothing in silence
    (issue #519). Only the documented names pass: a reader option given
    as a keyword is pointed to ``reader_options``.
    """
    problems = [_RETIRED_KEYWORDS[name] for name in kwargs if name in _RETIRED_KEYWORDS]
    misplaced = [name for name in kwargs if name in READER_OPTION_FORMATS]
    if misplaced:
        problems.append(
            _names(misplaced)
            + (" is a reader option" if len(misplaced) == 1 else " are reader options")
            + ": pass "
            + ("it" if len(misplaced) == 1 else "them")
            + " in reader_options"
        )
    unknown = [
        name
        for name in kwargs
        if name not in CONVERTER_KEYWORDS
        and name not in _RETIRED_KEYWORDS
        and name not in READER_OPTION_FORMATS
    ]
    if unknown:
        problems.append(
            f"convert_msi takes no keyword {_names(unknown)}; its converter "
            "keywords are " + ", ".join(CONVERTER_KEYWORDS)
        )
    if reader_options is None:
        return problems
    if not isinstance(reader_options, dict):
        problems.append(
            f"reader_options must be a dict, not {type(reader_options).__name__}"
        )
        return problems
    unknown = [key for key in reader_options if key not in READER_OPTION_FORMATS]
    if unknown:
        problems.append(
            f"reader_options holds {_names(unknown)}, which no reader takes. "
            "Known keys: " + ", ".join(READER_OPTION_FORMATS)
        )
    return problems


def _names(names: List[str]) -> str:
    """``'a'``, ``'a' and 'b'``, ``'a', 'b' and 'c'``."""
    quoted = [repr(name) for name in names]
    if len(quoted) == 1:
        return quoted[0]
    return ", ".join(quoted[:-1]) + " and " + quoted[-1]


def _validate_keywords(
    kwargs: Dict[str, Any], reader_options: Optional[Dict[str, Any]]
) -> bool:
    """Refuse, in one line, every keyword ``convert_msi`` does not use."""
    problems = keyword_problems(kwargs, reader_options)
    if problems:
        logger.error("Nothing was converted. %s.", "; ".join(problems))
        return False
    return True


def _drop_options_the_format_ignores(
    options: Dict[str, Any], input_format: str
) -> Dict[str, Any]:
    """Name each reader option this format does not use, and leave it out.

    The CLI said so for its flags (issue #260); the API said nothing, so a
    caller could believe a filter applied that did not (issue #519). One
    WARNING per option, worded as the CLI words it. An option left at
    ``None`` asks for nothing and is dropped without a word.
    """
    kept = {}
    for key, value in options.items():
        if reader_option_reaches(key, input_format):
            kept[key] = value
            continue
        if value is None:
            continue
        formats = READER_OPTION_FORMATS[key] or ()
        logger.warning(
            "reader_options[%r] is ignored on %s input: it only applies to "
            "%s. The store is written as though it had not been given.",
            key,
            format_name(input_format),
            " or ".join(format_name(f) for f in formats),
        )
    return kept


def _validate_input_parameters(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    format_type: str,
    dataset_id: str,
    pixel_size_um: Optional[float],
    z_spacing_um: Optional[float] = None,
    streaming: object = "auto",
    sample_id: Any = None,
    sample_name: Any = None,
) -> bool:
    """Validate all input parameters for convert_msi function."""
    return (
        _validate_paths_parameters(input_path, output_path)
        and _validate_string_parameters(format_type, dataset_id)
        and _validate_sample_parameters(sample_id, sample_name)
        and _validate_numeric_parameters(pixel_size_um, z_spacing_um)
        and _validate_streaming(streaming)
    )


def output_ancestor_problem(output: Path) -> Optional[str]:
    """Why the store cannot be created at ``output``, or ``None``.

    The existence check alone let an output whose parent is a regular
    file through: the reader was opened and the mass axis built before
    the scratch directory's ``mkdir`` failed with a raw ``WinError 183``
    that named a file the user had not asked about (issue #255). The
    store is a directory tree, so the deepest ancestor that exists has
    to be a directory Thyra may write into; the ones above it are
    created on the way, as they always were.

    Shared by the CLI and :func:`convert_msi`: the check used to live in
    the CLI only, so the API read the whole source before failing the
    same way (issue #519).
    """
    ancestor = next((p for p in output.parents if p.exists()), None)
    if ancestor is None:
        return None
    if not ancestor.is_dir():
        return (
            f"Cannot write {output}: {ancestor} is a file, not a directory. "
            "The output is a Zarr store, which is a directory tree, so every "
            "part of its path above it has to be a directory."
        )
    return _unwritable_directory_problem(output, ancestor)


def _unwritable_directory_problem(output: Path, ancestor: Path) -> Optional[str]:
    """Why the store cannot be created under ``ancestor``, or ``None``.

    The check is an actual ``mkdir``, not ``os.access(W_OK)``. On Windows
    ``os.access`` consults the read-only *attribute* and not the ACL, so
    it answers True for a directory the user has no write permission on
    -- measured here on a directory denied ``(WD,AD)`` via ``icacls``:
    ``os.access(W_OK) -> True``, ``st_mode -> 0o40777``, and both
    ``mkdir`` and ``open(w)`` raise ``PermissionError`` errno 13. So the
    refusal could never fire on the one platform this project is
    developed on, and the failure surfaced instead from deep inside the
    conversion, as "Error during conversion" plus a traceback out of the
    CSC scratch allocation, after the metadata scan (issue #312).

    ``mkdir`` rather than a temporary file because the thing being
    refused is a directory tree: the store and its scratch are both
    directories.

    Not ``tempfile.TemporaryDirectory(dir=ancestor)``, which is what the
    issue proposed: ``mkdtemp`` retries ``PermissionError`` when
    ``os.name == 'nt' and os.path.isdir(dir) and os.access(dir, W_OK)``
    -- gated on the very predicate that is broken here. It would spin
    ``TMP_MAX`` times (2,147,483,647 on this build, extrapolating to
    days) and then raise ``FileExistsError``, which is neither the right
    exception nor a usable wait.
    """
    probe = ancestor / f".thyra-write-probe-{uuid4().hex}"
    try:
        probe.mkdir()
    except OSError as e:
        return f"Cannot write {output}: {ancestor} is not writable ({e.strerror or e})"
    try:
        probe.rmdir()
    except OSError as e:
        # The directory is writable, which is the whole question, so a
        # failure to tidy up must not refuse the conversion. It is still
        # said out loud: an indexer or scanner holding the directory
        # would otherwise leave one of these behind per invocation with
        # nobody told.
        logger.warning(
            "Left the write probe %s behind (%s); it is an empty directory "
            "and safe to delete.",
            probe,
            e,
        )
    return None


def _validate_paths(input_path: Path, output_path: Path) -> bool:
    """Validate that input exists and the store can be created at output."""
    if not input_path.exists():
        logger.error(f"Input path does not exist: {input_path}")
        return False

    if output_path.exists():
        logger.error(f"Destination {output_path} already exists.")
        return False

    problem = output_ancestor_problem(output_path)
    if problem is not None:
        logger.error(problem)
        return False

    return True


def quarantine_partial_output(output: Path) -> None:
    """Move a partially written store aside after a failed conversion.

    A conversion that fails part-way through writing leaves an
    incomplete ``.zarr`` at the destination. That store cannot be opened
    (``spatialdata.read_zarr()`` raises), but it looks like a plausible
    artifact, and it also blocks a retry because
    :func:`_validate_paths` refuses to write to an existing path -- so
    the run that failed is the run that makes the same command fail
    again, with a message about the destination rather than about what
    went wrong. Rename it to a sibling ``.failed`` path so the
    destination is clear while the partial store remains available for
    diagnosis.

    This lives in ``convert_msi`` rather than in the CLI because
    ``convert_msi`` is the public entry point, and a library caller --
    Ousia converts through it, not through the CLI -- got no cleanup at
    all. The CLI's own call is kept and is idempotent, so a store already
    moved aside here is not moved twice (issue #293).

    Only called once the destination has been validated as not existing,
    so anything present is this run's own output and is safe to move.
    """
    if not output.exists():
        return

    quarantine = output.with_name(f"{output.name}.failed")
    attempt = 1
    while quarantine.exists():
        attempt += 1
        quarantine = output.with_name(f"{output.name}.failed{attempt}")

    try:
        output.rename(quarantine)
    except OSError as e:
        logger.error(
            "The incomplete output was left at %s because it could not be "
            "moved aside (%s). It will not open with "
            "spatialdata.read_zarr(); delete it before retrying.",
            output,
            e,
        )
        return

    logger.error(
        "The incomplete output was moved to %s. It will not open with "
        "spatialdata.read_zarr(); delete it once you no longer need it.",
        quarantine,
    )


def _create_reader(
    input_path: Path,
    reader_options: Optional[Dict[str, Any]] = None,
    lossless_tables: Optional[List[str]] = None,
) -> Tuple[BaseMSIReader, str]:
    """Create and return a reader for the input format.

    Args:
        input_path: Path to the input MSI data
        reader_options: Optional format-specific reader options (e.g., calibration settings)
        lossless_tables: Names the sibling tables being written that need
            the summed spectrum to keep all of the ion current
            (``"mobility grid"``, ``"demultiplexed MS/MS"``), so each of
            them adds back up to the summed one exactly; empty when none.

    Returns:
        Tuple of (reader instance, detected format string)
    """
    input_format = detect_format(input_path)
    logger.info(f"Detected format: {input_format}")
    reader_class = get_reader_class(input_format)
    logger.info(f"Using reader: {reader_class.__name__}")

    options = _drop_options_the_format_ignores(dict(reader_options or {}), input_format)
    if lossless_tables:
        _force_scan_sum(input_format, options, lossless_tables)
    return reader_class(input_path, **options), input_format


def _force_scan_sum(
    input_format: str, options: Dict[str, Any], wanted: List[str]
) -> None:
    """Say out loud when a sibling table will not add up to the summed one.

    A sibling table -- the mobility grid, the demultiplexed MS/MS table --
    is built from the raw scans, so it reproduces the summed table (the
    grid's marginal over channels, the split's blocks added back up) only
    when that table was built from the same scans. The default
    ``scan_sum`` is built from them, so nothing needs switching; the
    vendor centroid is a peak-picked spectrum over the same ramp that
    keeps only the current inside the peaks it picks (79-99% on measured
    imaging acquisitions), so an explicit ``--tdf-spectrum vendor_centroid`` is
    the caller saying they want the mismatch, which the sibling's ``uns``
    block then records. It is said at WARNING rather than overridden.

    Said before the reader is opened, because the spectrum mode binds
    into the vendor handle -- which is also why it cannot yet know
    whether the source is a PASEF acquisition at all. So it names the
    tables that *would* carry the mismatch; whether each is written is
    then said by name when the converter plans it.
    """
    mode = options.get("tdf_spectrum")
    if mode is None or mode == "scan_sum":
        return
    if input_format != "bruker":
        # The option reaches no other reader, so nothing about the summed
        # spectrum changed and there is no mismatch to announce. It has
        # already been named as ignored here (issues #260, #519).
        return
    subject = " and ".join(f"a {name} table" for name in wanted)
    logger.warning(
        "%s %s written with --tdf-spectrum %s. Such a table reads raw "
        "scans, so it will not add back up to the summed table; its uns "
        "block records by how much.",
        subject[0].upper() + subject[1:],
        "is" if len(wanted) == 1 else "are",
        mode,
    )


def _lossless_spectrum_for(kwargs: Dict[str, Any]) -> List[str]:
    """Which sibling tables, if any, need the lossless summed spectrum.

    ``msms_table`` defaults to *on* in the converter (design decision D2),
    and the CLI forwards the keyword only when the flag was given, so
    reading it as off-by-default meant the common path -- the default-on
    table -- never reached :func:`_force_scan_sum` and the mismatch it
    warns about went unsaid (issue #253). The default here is the
    converter's own.
    """
    wanted = []
    if kwargs.get("mobility_grid", False):
        wanted.append("mobility grid")
    if kwargs.get("msms_table", True):
        wanted.append("demultiplexed MS/MS")
    return wanted


def _determine_pixel_size(
    reader: BaseMSIReader, pixel_size_um: Optional[float], input_format: str
) -> Tuple[float, PixelSizeSource, Dict[str, Any]]:
    """Determine pixel size either from metadata or user input."""
    if pixel_size_um is not None:
        # Manual pixel size was provided
        pixel_size_detection_info = {
            "method": "manual",
            "source_format": input_format,
            "detection_successful": False,
            "note": "Pixel size manually specified via --pixel-size parameter",
        }
        return pixel_size_um, PixelSizeSource.USER_PROVIDED, pixel_size_detection_info

    # Attempt automatic detection
    logger.info("Attempting automatic pixel size detection...")
    essential_metadata = reader.get_essential_metadata()

    if essential_metadata.pixel_size is None:
        # Said once: the refusal is logged where it is caught. Logging it
        # here as well printed the same line twice.
        raise ConversionRefused(
            "Pixel size not found in metadata. Pass it with --pixel-size "
            "(for example --pixel-size 25)."
        )

    # A pitch that is NaN, infinite or not above 0 is no pitch: it used to
    # pass to the store, where JSON cannot hold it (#503). Refused like a
    # missing one, so the user can state the right value.
    if not all(is_usable_pitch(half) for half in essential_metadata.pixel_size):
        raise ConversionRefused(
            "The pixel size in the metadata is not a positive number "
            f"({essential_metadata.pixel_size[0]!r} x "
            f"{essential_metadata.pixel_size[1]!r} um). Pass it with "
            "--pixel-size (for example --pixel-size 25)."
        )

    final_pixel_size = essential_metadata.pixel_size[0]  # Use X size
    logger.info(f"Detected pixel size: {final_pixel_size:.1f} um")

    pixel_size_detection_info = {
        "method": "automatic",
        "detected_x_um": float(essential_metadata.pixel_size[0]),
        "detected_y_um": float(essential_metadata.pixel_size[1]),
        "source_format": input_format,
        "detection_successful": True,
        "note": "Pixel size automatically detected from source metadata",
    }

    return final_pixel_size, PixelSizeSource.AUTO_DETECTED, pixel_size_detection_info


def _create_converter(
    format_type: str,
    reader: BaseMSIReader,
    output_path: Path,
    dataset_id: str,
    pixel_size_um: float,
    pixel_size_source: PixelSizeSource,
    handle_3d: bool,
    pixel_size_detection_info: Dict[str, Any],
    resampling_config: Optional[
        Union[Dict[str, Any], ResamplingConfig, MassAxisPlan]
    ] = None,
    include_optical: bool = True,
    apply_optical_alignment: bool = True,
    streaming: Union[bool, Literal["auto"]] = "auto",
    z_spacing_um: Optional[float] = None,
    sample_id: Optional[str] = None,
    sample_name: Optional[str] = None,
    **kwargs: Any,
) -> BaseMSIConverter:
    """Create and return a converter for the specified format.

    ``streaming`` used to choose between an in-memory converter and a
    streaming one, on request or on an estimated size. There is one
    converter now and it streams (design decision D11), so the argument
    selects nothing; it is accepted so existing calls keep working, and
    ``False`` is answered with a note rather than silently. The value is
    checked by ``_validate_streaming`` before the reader is opened.
    """
    if streaming is False:
        logger.warning(
            "streaming=False asks for the in-memory converter, which was folded "
            "into the streaming one in v3.23: every conversion makes two passes "
            "over the source and never holds the matrix in memory. The output "
            "is the same; the argument no longer selects anything."
        )

    # Deliberately outside any try. This is often the first call that
    # builds the reader's parser, so it is the call a refused file fails on
    # -- and swallowing it once logged the reason at DEBUG, returned False,
    # and let the converter parse the whole file a second time before
    # failing the same way. On a 2.1 GB imzML that is about two minutes of
    # apparent progress with the real reason invisible. The converter's
    # own constructor extracts this metadata inside a try that logs at
    # DEBUG, which is exactly the swallow.
    reader.get_essential_metadata()

    converter_kwargs = {
        "dataset_id": dataset_id,
        "pixel_size_um": pixel_size_um,
        "pixel_size_source": pixel_size_source,
        "handle_3d": handle_3d,
        "z_spacing_um": z_spacing_um,
        "pixel_size_detection_info": pixel_size_detection_info,
        "resampling_config": resampling_config,
        "include_optical": include_optical,
        "apply_optical_alignment": apply_optical_alignment,
        # Explicit, not left to **kwargs: a sample name is not a converter
        # option, and the converter would take it as an unknown one.
        "sample_id": sample_id,
        "sample_name": sample_name,
        **kwargs,
    }
    # The registry's own refusal names the format asked for and the formats
    # there are. This used to go through a wrapper that, for any format_type
    # containing "spatialdata", replaced that with "SpatialData converter
    # unavailable" plus five lines of zarr-upgrade advice -- written for a
    # missing spatialdata, which registration no longer hides.
    converter_class = get_converter_class(format_type.lower())
    logger.info(f"Using converter: {converter_class.__name__}")
    return converter_class(reader, output_path, **converter_kwargs)


def _perform_conversion(converter: BaseMSIConverter) -> bool:
    """Run the conversion and say how it went.

    Closing the reader is not this function's job: ``convert_msi`` opened
    it and holds it in a ``with`` (issue #279).
    """
    logger.info("Starting conversion...")
    result = converter.convert()
    logger.info(f"Conversion {'completed successfully' if result else 'failed'}")
    return bool(result)


def _plain_values(
    handle_3d: Any, pixel_size_um: Any, z_spacing_um: Any
) -> Tuple[bool, Optional[float], Optional[float]]:
    """The validated arguments as plain Python types, taken once.

    An int pitch made an integer ``obs["spatial_x"]`` (issue #519), and a
    numpy bool is no JSON value.
    """
    return (
        bool(handle_3d),
        None if pixel_size_um is None else float(pixel_size_um),
        None if z_spacing_um is None else float(z_spacing_um),
    )


def _checked_resampling_config(config: Any) -> Tuple[Any, bool]:
    """The config normalised, and whether it passed the CLI's checks.

    Checked before the source is opened, as the CLI checks its flags,
    rather than when the converter is built (issue #519). The converter
    normalises the result again, which a ``ResamplingConfig`` passes
    unchanged and without a second word.
    """
    if config is None:
        return None, True
    from .resampling.axis_planner import normalize_resampling_config

    try:
        return normalize_resampling_config(config), True
    except ConversionRefused as e:
        logger.error("%s", str(e))
        return config, False


def _warn_z_spacing_without_3d(z_spacing_um: Optional[float], handle_3d: bool) -> None:
    """Say that a z spacing without ``handle_3d`` changes nothing.

    The 2D route writes one image per slice and never builds a z axis. A
    caller who set it believes their volume is calibrated.
    """
    if z_spacing_um is not None and not handle_3d:
        logger.warning(
            "z_spacing_um=%g was given without handle_3d=True and will be "
            "ignored: without 3D handling each slice is written as its own "
            "2D image and there is no z axis to space out. Pass "
            "handle_3d=True (--handle-3d) to build a volume.",
            z_spacing_um,
        )


def _with_region(
    reader_options: Optional[Dict[str, Any]], region: Optional[Union[int, str]]
) -> Optional[Dict[str, Any]]:
    """``reader_options`` with ``region`` merged in, when one was given."""
    if region is None:
        return reader_options
    return {**(reader_options or {}), "region": region}


def convert_msi(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    format_type: str = "spatialdata",
    dataset_id: str = "msi_dataset",
    pixel_size_um: Optional[float] = None,
    handle_3d: bool = False,
    z_spacing_um: Optional[float] = None,
    resampling_config: Optional[
        Union[Dict[str, Any], ResamplingConfig, MassAxisPlan]
    ] = None,
    reader_options: Optional[Dict[str, Any]] = None,
    include_optical: bool = True,
    apply_optical_alignment: bool = True,
    streaming: Union[bool, Literal["auto"]] = "auto",
    region: Optional[Union[int, str]] = None,
    sample_id: Optional[str] = None,
    sample_name: Optional[str] = None,
    **kwargs: Any,
) -> bool:
    """Convert MSI data to the specified format.

    Provides automatic pixel size detection from metadata or
    accepts user-specified values.

    Args:
        input_path: Path to input MSI data file or directory
        output_path: Path for output file
        format_type: Output format type (default: "spatialdata")
        dataset_id: Identifier for the dataset
        pixel_size_um: In-plane pixel size in micrometers (None for auto).
            Any finite positive real number, numpy's included; a bool is
            refused.
        handle_3d: Whether to process as 3D data (default: False). Only
            True or False is accepted, numpy's bool included.
        z_spacing_um: Distance between consecutive slices in micrometers.
            Only meaningful together with ``handle_3d=True``; ignored
            (with a warning) otherwise. ``None`` (default) means nothing
            supplied one, in which case the in-plane pitch is reused and
            the store records ``z_spacing_source="assumed_isotropic"``
            so the guess is not mistaken for a measurement. Supply it
            whenever the section thickness is known -- the in-plane
            pitch matches it only by coincidence.
        resampling_config: Optional resampling configuration: a dict, a
            ``ResamplingConfig``, or a ``MassAxisPlan`` from
            :func:`thyra.plan_mass_axis`. With a plan the conversion
            builds the planned axis and refuses (returns ``False``,
            nothing written) when its own axis hashes differently.
        reader_options: Optional format-specific reader options. Only the
            keys below are taken; any other key is refused before the
            source is opened. A key the detected format does not use is
            named in one WARNING and left out.

            - intensity_threshold: float - Minimum intensity
              to include. Default: None (no filtering).
            - use_recalibrated_state: bool - For Bruker data,
              use active/recalibrated calibration (default True).
            - tdf_spectrum: "scan_sum" | "vendor_centroid" - For Bruker
              TDF (TIMS) data, how a frame's mobility scans are collapsed
              into one spectrum per pixel (default "scan_sum").
            - region: the same as the ``region`` argument.
            - use_centroid: bool - For Waters .raw, whether MassLynx
              hands back the vendor centroid (True) or the profile
              trace (False). Default None: the profile trace on a
              SELECT SERIES MRT, the vendor centroid on every other
              Waters instrument. The CLI spells it --waters-spectrum.
            - spectrum_type: str - For imzML, declare the spectrum
              representation explicitly: 'profile' or 'centroid'
              (the full CV names are accepted too). Outranks the
              file's own MS:1000127/MS:1000128, so it can correct a
              file that declares the wrong thing; contradicting a
              declaration is logged as a warning. SCiLS Lab spells
              this --rep_type (2026b User Guide, p.81). Default:
              None, meaning detect. **Changes stored values** for
              files where it disagrees with what detection would
              have chosen, because the representation feeds
              instrument and axis-type selection.
            - max_mass_axis_length: int - For imzML, mzPeak, solariX
              and Waters sources converted with --no-resample, give up
              once the raw mass axis exceeds this many unique m/z
              values. This matters when the peak lists share no m/z
              values, where the raw axis grows to roughly one column
              per peak in the whole dataset. Default on imzML:
              10,000,000, which is the limit SCiLS Lab places on the
              same quantity (2026b User Guide, p.76); no limit
              elsewhere. Pass None for no limit.
            - use_appended_calibration: bool - For PHI data, prefer a
              recalibration appended to the file (default True).
        **kwargs: Forwarded to the converter. ``write_mobility_table``
            (default True) writes the mobility-resolved sibling table when
            the source shares one set of (m/z, mobility) features across
            pixels, e.g. an imzML export with a mobility array.
            ``mobility_heatmap`` (default True) stores the mean
            mass-mobility frame on the summed table as
            ``uns["mobility_heatmap"]`` whenever the source has an ion
            mobility dimension (Bruker TDF, imzML with a mobility array,
            Waters drift bins).
            ``mobility_grid`` (default False) fills that same sibling for
            a source that carries mobility per pixel rather than as a
            shared feature axis (Bruker TDF, processed imzML, Waters) by
            binning the point cloud
            onto a common mobility grid; ``mobility_bins`` (default: the
            heatmap's own channels, 256 or one per Waters drift bin),
            ``mobility_min`` and
            ``mobility_max`` size that grid. On a Bruker TDF its marginal
            reproduces the summed table exactly under the default
            ``tdf_spectrum="scan_sum"``; an explicit ``"vendor_centroid"``
            is kept, with a warning, and the mismatch recorded.
            ``msms_table`` (default True) writes the demultiplexed MS/MS
            sibling table when the source isolates several precursors per
            pixel in disjoint mobility slices (Bruker PASEF); refused with
            a reason, and nothing written, on any other source.
            ``use_csc`` is accepted for older callers and selects nothing.
            Any other keyword is refused before the source is opened,
            with one line naming it; a reader option given here is
            pointed to ``reader_options``.
        include_optical: Include optical images (default: True)
        apply_optical_alignment: If True (default) and the MSI source
            carries FlexImaging Area metadata, MSI elements are placed
            in optical-image pixel space at ``"global"``.  Set to
            False to keep MSI in pure micrometer coordinates -- needed
            when a downstream tool (e.g. Ousia's wizard) computes its
            own MSI-to-target registration.
        streaming: Kept for compatibility; selects nothing. Every
            conversion streams -- two passes over the source into
            memory-mapped arrays, the matrix never held in RAM -- since
            the in-memory converter was folded in (v3.23). ``False`` is
            accepted with a warning. Only ``True``, ``False`` and
            ``"auto"`` are accepted at all: a no-op argument still has to
            say when it was misspelled, and ``0`` -- which reads as
            ``False`` but is not it -- used to convert in silence.
        region: Select one region of a multi-region Bruker timsTOF
            ``.d``. Accepts an int (DB RegionNumber) or a str (matched
            against .mis Area Name, falling back to integer parse).
            None (default) converts all regions. Every other format
            refuses it: no other reader can select a region.
            :func:`thyra.list_regions` lists the regions of any input.
            Passed to the reader as reader_options["region"].
        sample_id: Your identifier for the sample, written to
            ``msi_metadata.sample.sample_id``. Never read from the
            source, because a vendor sample name can identify a person.
            One conversion stamps the same value on every region it
            converts. At most 256 printable characters. Default: None
            (nothing is recorded).
        sample_name: Your name for the sample, written to
            ``msi_metadata.sample.sample_name`` under the same rules as
            ``sample_id``. Default: None.

    Returns:
        True if the conversion completed and the store was written, False
        otherwise. Every failure comes back this way, refusals included:
        a refusal Thyra planned for is logged once at ``ERROR`` with the
        whole explanation and its traceback kept for ``DEBUG`` (issue
        #234), so a caller that wants the exception itself has to use the
        converter classes directly.

        **On failure the destination is cleared.** Anything written to
        ``output_path`` by this call is renamed to a sibling
        ``<name>.failed`` (``.failed2`` and so on if one is already
        there), so the same call can be retried without the caller
        deleting a directory first, and a half-written store is not left
        where a finished one is expected -- it would not open, but it
        looks like one. When the rename itself fails -- a live handle on
        the directory is the usual reason on Windows -- the partial store
        stays put and the failure says so at ``ERROR``, naming the path
        to delete; a retry is then refused until it is gone. An
        ``output_path`` that already existed when the call started is
        refused untouched, never renamed. See
        :func:`quarantine_partial_output` (issue #293).
    """
    # Validate input parameters
    if not _validate_input_parameters(
        input_path,
        output_path,
        format_type,
        dataset_id,
        pixel_size_um,
        z_spacing_um,
        streaming,
        sample_id,
        sample_name,
    ) or not (
        _validate_handle_3d(handle_3d) and _validate_keywords(kwargs, reader_options)
    ):
        return False
    handle_3d, pixel_size_um, z_spacing_um = _plain_values(
        handle_3d, pixel_size_um, z_spacing_um
    )
    resampling_config, ok = _checked_resampling_config(resampling_config)
    if not ok:
        return False
    _warn_z_spacing_without_3d(z_spacing_um, handle_3d)

    # Convert to Path objects and validate
    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    logger.info(f"Processing input file: {input_path}")

    if not _validate_paths(input_path, output_path):
        return False

    # A Zarr store's deepest key sits far below the path the caller named,
    # so a legal-looking output path can still blow the Windows 260
    # character limit part-way through the write. No-op elsewhere.
    output_path = prepare_zarr_output_path(output_path, dataset_id)
    # Past this line the destination has been validated as not existing,
    # so anything that appears there is this run's own output and this
    # run is responsible for clearing it away if it fails.
    succeeded = False
    # Bound once created, so an interrupt that escapes the converter can
    # still ask it whether the store was already complete (issue #500).
    converter: Optional[BaseMSIConverter] = None

    reader_options = _with_region(reader_options, region)

    try:
        # Create reader with format-specific options. A mobility grid
        # table decides the summed spectrum's semantics, so it has to be
        # known before the reader is opened -- the choice is bound into
        # the SDK handle.
        reader, input_format = _create_reader(
            input_path,
            reader_options,
            lossless_tables=_lossless_spectrum_for(kwargs),
        )

        # This function opened the reader, so this function closes it, on
        # every path out -- including the ones that fail before the
        # conversion starts. "Pixel size not found in metadata" is the
        # common one, and it used to leave the source open until the
        # garbage collector happened to reach it, which on Windows holds a
        # lock on the file the user is about to retry with (issue #279).
        with reader:
            # Determine pixel size
            final_pixel_size, pixel_size_source, pixel_size_detection_info = (
                _determine_pixel_size(reader, pixel_size_um, input_format)
            )

            # Create converter
            converter = _create_converter(
                format_type,
                reader,
                output_path,
                dataset_id,
                final_pixel_size,
                pixel_size_source,
                handle_3d,
                pixel_size_detection_info,
                resampling_config,
                include_optical=include_optical,
                apply_optical_alignment=apply_optical_alignment,
                streaming=streaming,
                z_spacing_um=z_spacing_um,
                sample_id=sample_id,
                sample_name=sample_name,
                **kwargs,
            )

            succeeded = _perform_conversion(converter)
            return succeeded

    except ConversionRefused as e:
        # A refusal Thyra wrote: the message is the whole explanation, and
        # a traceback in front of it reads like a crash the code did not
        # plan for. The traceback is still there under --log-level DEBUG,
        # for the cases where the refusal itself is the surprise (#234).
        #
        # ``str(e)``, never ``e``: a log handler that retains records --
        # pytest's capture, Ousia's per-session log -- would otherwise hold
        # the exception, its traceback and every frame in it, which pins
        # the reader and leaves its file open (see issue #249).
        logger.error("%s", str(e))
        logger.debug("Refusal raised at:\n%s", traceback.format_exc())
        return False

    except KeyboardInterrupt:
        # The converter catches its own interrupt (BaseMSIConverter.convert);
        # this covers the steps before it -- opening the reader, detecting
        # the pixel size -- so wherever Ctrl-C lands the caller gets False
        # and the CLI's post-conversion step runs. Without it click printed
        # ``Aborted!`` and no partial store was ever moved aside (#245).
        #
        # It also catches the ones that land after the converter has
        # finished the store -- in its scratch release, in closing the
        # reader on the way out of the ``with``. That store is complete
        # and opens, so it stays where it is (#500). ``is True`` because
        # tests stand a mock in for the converter. ``succeeded`` first, so
        # a second Ctrl-C during the log line cannot reach the quarantine.
        if converter is not None and converter.store_complete is True:
            succeeded = True
            log_interrupt_after_save(output_path)
            return True
        logger.error("Conversion interrupted before it started reading.")
        return False

    except Exception as e:
        logger.error(f"Error during conversion: {e}")
        logger.error(f"Detailed traceback:\n{traceback.format_exc()}")
        return False

    finally:
        # Runs after the ``with`` above has closed the reader, so the move
        # cannot race a handle the source still holds. Covers all five
        # ``return False`` paths past the path validation plus anything
        # that escapes this function.
        if not succeeded:
            quarantine_partial_output(output_path)
