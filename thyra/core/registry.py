"""Which format a path is, and which classes handle it.

Two jobs that look like one. Registration is a decorator
(:func:`register_reader`, :func:`register_converter`) that runs at import
time, so importing a reader's module is what populates the tables; nothing
scans the filesystem for plugins. Detection is :func:`detect_format`, and
it is the harder half.

**Which module to import is a table here, and the import happens on first
use** (issue #381). ``thyra/readers/__init__.py`` used to import every
reader package so the decorators ran, which made ``import thyra`` -- and
therefore ``import thyra.metadata.schema``, since a parent package always
initialises first -- pull in every reader, every converter and
``spatialdata``: 3.5 s and 3,580 modules on the interpreter these numbers
were measured on, before any work. :data:`_READER_MODULES` and
:data:`_CONVERTER_MODULES` name the module that declares each format's
decorator, and a lookup imports it if the format is not in the table yet.

There is still one resolution path, not two: every lookup calls
:meth:`MSIRegistry._ensure_registered` first, whether or not the format is
already there, and a module that has already been imported costs a
``sys.modules`` hit. Nothing registers twice, because nothing imports
twice.

**Extensions are not enough, which is why detection inspects the path.**
``.imzML`` and ``.d`` map straight to a format, but ``.raw`` is claimed by
both Waters and PHI and can only be told apart by looking inside, and a
``.d`` directory may be timsTOF, solariX or a flexImaging pre-scan that is
not an imaging run at all. So the detectors here open files and read
headers rather than trusting the name. Each refusal names the path and
says what was missing, because "unsupported format" on a directory the
user believes is supported is the least actionable message this package
could give.

The ``RLock`` guards the registration tables and nothing else. It is not a
general thread-safety claim: a table gains a format the first time that
format is looked up and is only read afterwards. The lazy import below used to memoise
into an unsynchronised module global -- the one piece of mutable global
state here, and the one piece the lock did not reach -- and that global is
gone rather than locked (issue #284). :func:`_get_bruker_folder_structure`
says why locking it would have been worse than leaving it alone.
"""

import importlib
import logging
import re
import zipfile
from pathlib import Path
from threading import RLock
from typing import Dict, Mapping, NoReturn, Type

from ..errors import ConversionRefused
from .base_converter import BaseMSIConverter
from .base_reader import BaseMSIReader

logger = logging.getLogger(__name__)

#: Which module declares ``@register_reader(format)`` for each format.
#:
#: The decorators remain the only thing that registers; this says where to
#: find them so the import can wait until the format is asked for. Keeping
#: the two in step is not left to care: ``tests/unit/test_registry.py``
#: parses every module under ``thyra/readers/`` for decorator calls and
#: asserts this table names exactly those formats, and for each one the
#: module the decorator is actually written in.
_READER_MODULES: Mapping[str, str] = {
    "bruker": "thyra.readers.bruker.timstof.timstof_reader",
    "imzml": "thyra.readers.imzml.imzml_reader",
    "mzpeak": "thyra.readers.mzpeak.mzpeak_reader",
    "phi": "thyra.readers.phi.phi_reader",
    "rapiflex": "thyra.readers.bruker.rapiflex.rapiflex_reader",
    "solarix": "thyra.readers.bruker.solarix.solarix_reader",
    "waters": "thyra.readers.waters.waters_reader",
}

#: The same table for output formats. One entry, and the import it names is
#: what pulls in ``spatialdata`` -- which is the whole 2.7 s the metadata
#: layer used to pay for nothing (issue #381).
_CONVERTER_MODULES: Mapping[str, str] = {
    "spatialdata": "thyra.converters.spatialdata.converter",
}


def _get_bruker_folder_structure():
    """Import ``BrukerFolderStructure`` and ``BrukerFormat`` late.

    Two reasons, and the first outlived the second. The import used to
    close a cycle -- their package reached this module back through
    ``thyra.readers.__init__``, which imported every reader package so the
    ``@register_reader`` decorators ran. Issue #381 emptied that front
    door, so there is no cycle left to break.

    What is left is that detection must not drag the Bruker reader stack
    in. This function runs on every ``.d`` path and on every directory
    that reached :meth:`MSIRegistry._detect_directory_format`, including
    Waters ones, and ``folder_structure`` costs a ``.mis`` parser and
    nothing more. A module-scope import would put it on the path of
    ``import thyra`` again.

    It is a plain import and deliberately not memoised. It used to assign
    the pair to a module-level global under an unsynchronised
    check-then-set -- the only mutable global state in this module, and the
    only part ``MSIRegistry._lock`` did not cover (issue #284). Two things
    decided deleting the global over locking it:

    * The memo bought 0.07 us per call (measured, 200k calls: 0.035 us
      against 0.105 us), because ``sys.modules`` already caches the import
      and a repeat is two dict lookups. This function runs about three
      times per conversion.
    * Holding the registry lock across an import is a deadlock shape, not
      merely a slow one. The holder waits for the ``thyra.readers`` import
      lock; a thread part-way through importing that package holds it and
      waits for the registry lock inside a ``@register_reader`` decorator.
      Neither can finish. Reproduced on CPython 3.13.3 against this exact
      pair of locks: holding the lock across the import hangs both threads
      indefinitely, while not holding it completes the same interleaving in
      1.1 ms.
    """
    from ..readers.bruker.folder_structure import BrukerFolderStructure, BrukerFormat

    return BrukerFolderStructure, BrukerFormat


class MSIRegistry:
    """Which classes handle a format, and which format a path is.

    Populated at import, read for the rest of the process: the
    ``@register_reader`` and ``@register_converter`` decorators fill the two
    tables while ``thyra.readers`` imports, and nothing registers after
    that.

    ``_lock`` covers those two tables and nothing else. It is not a claim
    that the class is safe to use concurrently in general -- the detection
    half holds no state to protect, and callers share one module-level
    instance. What it buys is that the compound accesses cannot see a
    half-applied registration: :meth:`get_reader_class` tests membership and
    then indexes, and on the failure path iterates the table to list what is
    available. The single dict operations underneath are already atomic in
    CPython; those two sequences of them are not.
    """

    def __init__(self):
        """Initialize the MSI registry."""
        self._lock = RLock()
        self._readers: Dict[str, Type[BaseMSIReader]] = {}
        self._converters: Dict[str, Type[BaseMSIConverter]] = {}
        # Extension mapping for file-based formats. '.raw' is deliberately
        # absent: it is claimed by both Waters and PHI and can only be
        # resolved by inspecting the path (see _detect_raw_format).
        self._extension_to_format = {
            ".imzml": "imzml",
            ".d": "bruker",
        }

    def register_reader(
        self, format_name: str, reader_class: Type[BaseMSIReader]
    ) -> None:
        """Register reader class."""
        with self._lock:
            self._readers[format_name] = reader_class
            logger.info(
                f"Registered reader {reader_class.__name__} for format "
                f"'{format_name}'"
            )

    def register_converter(
        self, format_name: str, converter_class: Type[BaseMSIConverter]
    ) -> None:
        """Register converter class."""
        with self._lock:
            self._converters[format_name] = converter_class
            logger.info(
                f"Registered converter {converter_class.__name__} for format "
                f"'{format_name}'"
            )

    def _detect_bruker_format(self, path: Path) -> str:
        """Detect Bruker data format using BrukerFolderStructure.

        Uses the unified BrukerFolderStructure to detect whether the path
        contains timsTOF, Rapiflex, or solariX data.

        Args:
            path: Path to check

        Returns:
            Format name ('bruker' for timsTOF, 'rapiflex' for Rapiflex,
            'solarix' for solariX, or empty string if not a Bruker format)
        """
        BrukerFolderStructure, BrukerFormat = _get_bruker_folder_structure()

        try:
            detected_format = BrukerFolderStructure.detect_format(path)
        except Exception:  # nosec B110 - intentionally ignore detection errors
            # Format detection failure means this is not a Bruker format
            return ""

        if detected_format == BrukerFormat.TIMSTOF:
            return "bruker"
        elif detected_format == BrukerFormat.RAPIFLEX:
            return "rapiflex"
        elif detected_format == BrukerFormat.SOLARIX:
            return "solarix"

        return ""

    def _detect_waters_format(self, path: Path) -> bool:
        """Check if path contains Waters .raw data.

        Waters .raw directories contain _FUNC[0-9]{3}.DAT files
        for each acquisition function.

        Args:
            path: Directory path to check

        Returns:
            True if Waters _FUNC*.DAT files are found
        """
        func_pattern = re.compile(r"_FUNC\d{3}\.DAT", re.IGNORECASE)
        try:
            for item in path.iterdir():
                if func_pattern.match(item.name):
                    return True
        except PermissionError:
            pass
        return False

    def _detect_phi_format(self, path: Path) -> bool:
        """Check if path is a PHI SmartSoft-TOF .raw file.

        PHI writes a single file whose ASCII header opens with the four-byte
        ``SOFH`` magic. This is what distinguishes it from Waters .raw data,
        which is a directory.

        Args:
            path: File path to check

        Returns:
            True if the file starts with the PHI SOFH magic
        """
        try:
            with path.open("rb") as handle:
                return handle.read(4) == b"SOFH"
        except (OSError, PermissionError):
            return False

    def detect_format(self, input_path: Path) -> str:
        """Detect MSI format from input path.

        Supports:
        - .imzml files (ImzML format)
        - .d directories (Bruker timsTOF)
        - .d directories (Bruker solariX / MRMS, via peaks.sqlite)
        - Folders with .dat + _poslog.txt (Bruker Rapiflex)
        - .raw directories (Waters MassLynx)
        - .raw files (PHI SmartSoft-TOF ToF-SIMS)
        """
        if not input_path.exists():
            raise ConversionRefused(f"Input path does not exist: {input_path}")

        format_name = self._detect_format_name(input_path)
        self._validate_format(format_name, input_path)
        return format_name

    def _detect_format_name(self, input_path: Path) -> str:
        """Detect format name from path extension or directory structure."""
        extension = input_path.suffix.lower()
        if extension == ".imzml":
            return "imzml"
        if extension == ".mzpeak":
            return self._detect_mzpeak_format(input_path)
        if extension == ".d":
            return self._detect_bruker_d_format(input_path)
        if extension == ".raw":
            return self._detect_raw_format(input_path)
        if extension in (".imdx", ".kbd"):
            self._raise_shimadzu_in_development(input_path)
        if input_path.is_dir():
            return self._detect_directory_format(input_path)
        self._raise_unsupported_format(input_path)

    def _detect_bruker_d_format(self, input_path: Path) -> str:
        """Validate and detect Bruker format from .d directory."""
        if not input_path.is_dir():
            raise ConversionRefused(
                "Bruker format requires .d directory, " f"got file: {input_path}"
            )
        bruker_format = self._detect_bruker_format(input_path)
        if bruker_format:
            return bruker_format
        BrukerFolderStructure, _ = _get_bruker_folder_structure()
        if BrukerFolderStructure.is_solarix_without_peaks(input_path):
            raise ConversionRefused(
                f"solariX .d directory without peaks.sqlite: {input_path}. "
                "The acquisition holds raw transients (ser) but no processed "
                "peak store, which is what Thyra reads. Export the dataset "
                "as imzML from the Bruker software (DataAnalysis, SCiLS Lab, "
                "or flexImaging) and convert the imzML file instead."
            )
        raise ConversionRefused(
            "Bruker .d directory missing analysis " f"files: {input_path}"
        )

    def _detect_mzpeak_format(self, input_path: Path) -> str:
        """Validate that a ``.mzpeak`` path really is an mzPeak archive.

        Three checks, cheapest first: the extension (already matched by the
        caller), the ZIP local-file-header magic, and the presence of the
        index member. The extension alone is not enough -- it is a young
        format and the name gets attached to loose Parquet directories -- and
        a file failing here should say so rather than failing later inside
        pyarrow with a message about a corrupt footer.

        Args:
            input_path: Path with a ``.mzpeak`` extension.

        Returns:
            ``"mzpeak"``.

        Raises:
            ValueError: If the file is not a ZIP or carries no index member.
        """
        if input_path.is_dir():
            raise ConversionRefused(
                f"mzPeak format requires a .mzpeak archive file, got "
                f"directory: {input_path}"
            )
        try:
            with input_path.open("rb") as handle:
                magic = handle.read(4)
        except (OSError, PermissionError) as exc:
            raise ConversionRefused(f"Cannot read {input_path}: {exc}") from exc

        if magic != b"PK\x03\x04":
            raise ConversionRefused(
                f"Not an mzPeak archive (missing ZIP signature): {input_path}"
            )

        if not zipfile.is_zipfile(input_path):
            raise ConversionRefused(
                f"Not an mzPeak archive (unreadable ZIP container): " f"{input_path}"
            )
        with zipfile.ZipFile(input_path) as archive:
            if "mzpeak_index.json" not in archive.namelist():
                raise ConversionRefused(
                    f"Not an mzPeak archive (no mzpeak_index.json member): "
                    f"{input_path}"
                )
        return "mzpeak"

    def _detect_raw_format(self, input_path: Path) -> str:
        """Resolve a .raw path to the vendor that wrote it.

        Two vendors claim the extension and they are distinguished by shape:
        Waters .raw is a directory of _FUNC*.DAT files, PHI .raw is a single
        file whose header begins with the SOFH magic.
        """
        if input_path.is_dir():
            if self._detect_waters_format(input_path):
                return "waters"
            raise ConversionRefused(
                "Waters .raw directory missing " f"_FUNC*.DAT files: {input_path}"
            )
        if self._detect_phi_format(input_path):
            return "phi"
        raise ConversionRefused(
            f"Unrecognised .raw file: {input_path}. Expected either a Waters "
            "directory containing _FUNC*.DAT files, or a PHI SmartSoft-TOF "
            "file beginning with the SOFH magic."
        )

    def _detect_directory_format(self, input_path: Path) -> str:
        """Detect vendor format from generic directory structure."""
        bruker_format = self._detect_bruker_format(input_path)
        if bruker_format:
            return bruker_format
        if self._detect_waters_format(input_path):
            return "waters"
        self._raise_unsupported_format(input_path)

    def _raise_shimadzu_in_development(self, input_path: Path) -> NoReturn:
        """Raise a guidance error for recognised Shimadzu imaging data.

        Shimadzu iMScope data (.kbd from LabSolutions, .imdx from
        IMAGEREVEAL MS) is recognised but native reading is still in
        development. Until it lands, IMAGEREVEAL MS can export the data
        as imzML, which Thyra converts today.
        """
        raise ConversionRefused(
            f"Shimadzu imaging data detected: {input_path}. Native support "
            "for Shimadzu formats (.imdx, .kbd) is in development. In the "
            "meantime, export the dataset as imzML from IMAGEREVEAL MS and "
            "convert the imzML file instead."
        )

    def _raise_unsupported_format(self, input_path: Path) -> NoReturn:
        """Raise error for unsupported format."""
        available = [
            ".imzml",
            ".mzpeak (HUPO-PSI, experimental)",
            ".d (timsTOF)",
            ".d (solariX/MRMS)",
            "folder (Rapiflex)",
            ".raw directory (Waters)",
            ".raw file (PHI SmartSoft-TOF)",
        ]
        raise ConversionRefused(
            f"Unsupported format for '{input_path}'. "
            f"Supported: {', '.join(available)}"
        )

    def _validate_format(self, format_name: str, input_path: Path) -> None:
        """Validate format-specific requirements."""
        if format_name == "imzml":
            ibd_path = input_path.with_suffix(".ibd")
            if not ibd_path.exists():
                raise ConversionRefused(
                    f"ImzML file requires corresponding .ibd file: {ibd_path}"
                )
        elif format_name == "bruker":
            self._validate_bruker_format(input_path)

    def _validate_bruker_format(self, input_path: Path) -> None:
        """Validate Bruker .d directory structure using BrukerFolderStructure.

        The detection already validated the format, but we do additional
        checks here for better error messages.
        """
        BrukerFolderStructure, BrukerFormat = _get_bruker_folder_structure()

        if not input_path.is_dir():
            raise ConversionRefused(
                f"Bruker format requires .d directory, got file: {input_path}"
            )

        # Use BrukerFolderStructure for validation
        try:
            folder = BrukerFolderStructure(input_path)
            info = folder.analyze()
            if info.format == BrukerFormat.UNKNOWN:
                raise ConversionRefused(
                    f"Bruker .d directory missing analysis files: {input_path}"
                )
        except Exception as e:
            if "missing analysis" in str(e).lower():
                raise
            # Re-check for required files
            has_tsf = (input_path / "analysis.tsf").exists()
            has_tdf = (input_path / "analysis.tdf").exists()
            if not has_tsf and not has_tdf:
                raise ConversionRefused(
                    f"Bruker .d directory missing analysis files: {input_path}"
                ) from e

    def _ensure_registered(
        self, format_name: str, table: Dict[str, type], modules: Mapping[str, str]
    ) -> None:
        """Import the module that registers ``format_name``, if there is one.

        Called before every lookup, so there is one resolution path rather
        than an eager one and a lazy one. A format already in ``table``
        returns after a dict lookup; an unknown format with no table entry
        returns too, and the caller raises the refusal that names what is
        available.

        **The import is taken outside the lock, deliberately.** Holding it
        across an import is the deadlock shape issue #284 measured on this
        exact pair of locks: the holder waits for the module's import lock
        while a thread part-way through that import waits for the registry
        lock inside a ``@register_reader`` decorator. Releasing first costs
        a re-check -- two threads can import the same module concurrently,
        which CPython's per-module import lock serialises anyway, and the
        decorator writing the same class into the same key twice is
        harmless.

        A module that has already been imported does not register again:
        ``importlib.import_module`` returns the cached module without
        re-running it. That matters only where something empties the tables
        by hand (``TestRegistry`` in ``tests/unit/test_registry.py``, which
        restores them itself); nothing in the package does.
        """
        with self._lock:
            if format_name in table:
                return
            module_path = modules.get(format_name)

        if module_path is None:
            return

        importlib.import_module(module_path)

    def get_reader_class(self, format_name: str) -> Type[BaseMSIReader]:
        """Get reader class."""
        self._ensure_registered(format_name, self._readers, _READER_MODULES)
        with self._lock:
            if format_name not in self._readers:
                available = self._known(self._readers, _READER_MODULES)
                raise ConversionRefused(
                    f"No reader for format "
                    f"'{format_name}'. "
                    f"Available: {available}"
                )
            return self._readers[format_name]

    def get_converter_class(self, format_name: str) -> Type[BaseMSIConverter]:
        """Get converter class."""
        self._ensure_registered(format_name, self._converters, _CONVERTER_MODULES)
        with self._lock:
            if format_name not in self._converters:
                available = self._known(self._converters, _CONVERTER_MODULES)
                raise ConversionRefused(
                    f"No converter for format '{format_name}'. Available: "
                    f"{available}"
                )
            return self._converters[format_name]

    @staticmethod
    def _known(table: Mapping[str, type], modules: Mapping[str, str]) -> list:
        """List the formats this registry can serve, loaded or not.

        The refusal used to list the keys registered so far. That was the
        same set every time while every reader was imported up front; with
        the imports deferred it would be whatever the process happened to
        have touched, so a user who mistyped ``imzml`` before converting
        anything would be told ``Available: []``. The union of the two
        tables is the answer that does not depend on history.
        """
        return sorted(set(table) | set(modules))


# Global registry instance
_registry = MSIRegistry()


# Simple public interface
def detect_format(input_path: Path) -> str:
    """Detect MSI format from input path.

    Args:
        input_path: Path to MSI data file or directory

    Returns:
        Format name ('imzml', 'bruker', 'rapiflex', 'solarix', 'waters',
        or 'phi')
    """
    return _registry.detect_format(input_path)


def get_reader_class(format_name: str) -> Type[BaseMSIReader]:
    """Get reader class for format.

    Args:
        format_name: MSI format name

    Returns:
        Reader class for the format
    """
    return _registry.get_reader_class(format_name)


def get_converter_class(format_name: str) -> Type[BaseMSIConverter]:
    """Get converter class for format.

    Args:
        format_name: MSI format name

    Returns:
        Converter class for the format
    """
    return _registry.get_converter_class(format_name)


def register_reader(format_name: str):
    """Decorator for reader registration."""

    def decorator(cls: Type[BaseMSIReader]):
        _registry.register_reader(format_name, cls)
        return cls

    return decorator


def register_converter(format_name: str):
    """Decorator for converter registration."""

    def decorator(cls: Type[BaseMSIConverter]):
        _registry.register_converter(format_name, cls)
        return cls

    return decorator
