# thyra/core/registry.py
import logging
import re
from pathlib import Path
from threading import RLock
from typing import Dict, NoReturn, Type

from .base_converter import BaseMSIConverter
from .base_reader import BaseMSIReader

logger = logging.getLogger(__name__)

# Import BrukerFolderStructure for unified Bruker format detection
# This avoids circular imports by importing lazily in the method
_bruker_folder_structure_module = None


def _get_bruker_folder_structure():
    """Lazy import of BrukerFolderStructure to avoid circular imports."""
    global _bruker_folder_structure_module
    if _bruker_folder_structure_module is None:
        from ..readers.bruker.folder_structure import (
            BrukerFolderStructure,
            BrukerFormat,
        )

        _bruker_folder_structure_module = (BrukerFolderStructure, BrukerFormat)
    return _bruker_folder_structure_module


class MSIRegistry:
    """Thread-safe registry with format detection for MSI data."""

    def __init__(self):
        """Initialize the MSI registry."""
        self._lock = RLock()
        self._readers: Dict[str, Type[BaseMSIReader]] = {}
        self._converters: Dict[str, Type[BaseMSIConverter]] = {}
        # Extension mapping for file-based formats
        self._extension_to_format = {
            ".imzml": "imzml",
            ".d": "bruker",
            ".raw": "waters",
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
        contains timsTOF or Rapiflex data.

        Args:
            path: Path to check

        Returns:
            Format name ('bruker' for timsTOF, 'rapiflex' for Rapiflex,
            or empty string if not a Bruker format)
        """
        BrukerFolderStructure, BrukerFormat = _get_bruker_folder_structure()

        try:
            detected_format = BrukerFolderStructure.detect_format(path)
            if detected_format == BrukerFormat.TIMSTOF:
                return "bruker"
            elif detected_format == BrukerFormat.RAPIFLEX:
                return "rapiflex"
        except Exception:  # nosec B110 - intentionally ignore detection errors
            pass  # Format detection failure means this is not a Bruker format

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

    def detect_format(self, input_path: Path) -> str:
        """Detect MSI format from input path.

        Supports:
        - .imzml files (ImzML format)
        - .d directories (Bruker timsTOF)
        - Folders with .dat + _poslog.txt (Bruker Rapiflex)
        - .raw directories (Waters MassLynx)
        """
        if not input_path.exists():
            raise ValueError(f"Input path does not exist: {input_path}")

        format_name = self._detect_format_name(input_path)
        self._validate_format(format_name, input_path)
        return format_name

    def _detect_format_name(self, input_path: Path) -> str:
        """Detect format name from path extension or directory structure."""
        extension = input_path.suffix.lower()
        if extension == ".imzml":
            return "imzml"
        if extension == ".d":
            return self._detect_bruker_d_format(input_path)
        if extension == ".raw":
            return self._detect_waters_raw_format(input_path)
        if input_path.is_dir():
            return self._detect_directory_format(input_path)
        self._raise_unsupported_format(input_path)

    def _detect_bruker_d_format(self, input_path: Path) -> str:
        """Validate and detect Bruker format from .d directory."""
        if not input_path.is_dir():
            raise ValueError(
                "Bruker format requires .d directory, " f"got file: {input_path}"
            )
        bruker_format = self._detect_bruker_format(input_path)
        if bruker_format:
            return bruker_format
        raise ValueError("Bruker .d directory missing analysis " f"files: {input_path}")

    def _detect_waters_raw_format(self, input_path: Path) -> str:
        """Validate and detect Waters format from .raw directory."""
        if not input_path.is_dir():
            raise ValueError(
                "Waters format requires .raw directory, " f"got file: {input_path}"
            )
        if self._detect_waters_format(input_path):
            return "waters"
        raise ValueError(
            "Waters .raw directory missing " f"_FUNC*.DAT files: {input_path}"
        )

    def _detect_directory_format(self, input_path: Path) -> str:
        """Detect vendor format from generic directory structure."""
        bruker_format = self._detect_bruker_format(input_path)
        if bruker_format:
            return bruker_format
        if self._detect_waters_format(input_path):
            return "waters"
        self._raise_unsupported_format(input_path)

    def _raise_unsupported_format(self, input_path: Path) -> NoReturn:
        """Raise error for unsupported format."""
        available = [
            ".imzml",
            ".d (timsTOF)",
            "folder (Rapiflex)",
            ".raw (Waters)",
        ]
        raise ValueError(
            f"Unsupported format for '{input_path}'. "
            f"Supported: {', '.join(available)}"
        )

    def _validate_format(self, format_name: str, input_path: Path) -> None:
        """Validate format-specific requirements."""
        if format_name == "imzml":
            ibd_path = input_path.with_suffix(".ibd")
            if not ibd_path.exists():
                raise ValueError(
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
            raise ValueError(
                f"Bruker format requires .d directory, got file: {input_path}"
            )

        # Use BrukerFolderStructure for validation
        try:
            folder = BrukerFolderStructure(input_path)
            info = folder.analyze()
            if info.format == BrukerFormat.UNKNOWN:
                raise ValueError(
                    f"Bruker .d directory missing analysis files: {input_path}"
                )
        except Exception as e:
            if "missing analysis" in str(e).lower():
                raise
            # Re-check for required files
            has_tsf = (input_path / "analysis.tsf").exists()
            has_tdf = (input_path / "analysis.tdf").exists()
            if not has_tsf and not has_tdf:
                raise ValueError(
                    f"Bruker .d directory missing analysis files: {input_path}"
                ) from e

    def get_reader_class(self, format_name: str) -> Type[BaseMSIReader]:
        """Get reader class."""
        with self._lock:
            if format_name not in self._readers:
                available = list(self._readers.keys())
                raise ValueError(
                    f"No reader for format "
                    f"'{format_name}'. "
                    f"Available: {available}"
                )
            return self._readers[format_name]

    def get_converter_class(self, format_name: str) -> Type[BaseMSIConverter]:
        """Get converter class."""
        with self._lock:
            if format_name not in self._converters:
                available = list(self._converters.keys())
                raise ValueError(
                    f"No converter for format '{format_name}'. Available: "
                    f"{available}"
                )
            return self._converters[format_name]


# Global registry instance
_registry = MSIRegistry()


# Simple public interface
def detect_format(input_path: Path) -> str:
    """Detect MSI format from input path.

    Args:
        input_path: Path to MSI data file or directory

    Returns:
        Format name ('imzml', 'bruker', 'rapiflex',
        or 'waters')
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
