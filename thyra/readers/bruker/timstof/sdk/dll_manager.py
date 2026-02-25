"""DLL/Library manager for Bruker SDK integration.

This module provides robust loading and management of the Bruker SDK
libraries with comprehensive error handling and fallback mechanisms.
"""

import logging
import platform
from ctypes import CDLL, cdll
from pathlib import Path
from typing import Any, Optional

# Import windll only on Windows
if platform.system() == "Windows":
    from ctypes import windll
else:
    windll = None

from .....utils.bruker_exceptions import SDKError
from .platform_detector import PlatformDetector, get_dll_paths, validate_library_path

logger = logging.getLogger(__name__)


class DLLManager:
    """Manages loading and access to Bruker SDK libraries.

    This class provides a singleton-like interface for managing the SDK
    library with automatic discovery, validation, and error handling.
    """

    _instance: Optional["DLLManager"] = None
    _dll: Optional[CDLL] = None
    _library_path: Optional[Path] = None

    def __new__(
        cls, data_directory: Optional[Path] = None, force_reload: bool = False
    ) -> "DLLManager":
        """Create or return existing DLLManager instance.

        Args:
            data_directory: Optional data directory for local library discovery
            force_reload: Force reloading of the library even if already loaded
        """
        if cls._instance is None or force_reload:
            cls._instance = super().__new__(cls)
            cls._instance._initialize(data_directory, force_reload)
        return cls._instance

    def _initialize(
        self, data_directory: Optional[Path] = None, force_reload: bool = False
    ) -> None:
        """Initialize the DLL manager and load the library.

        Args:
            data_directory: Optional data directory for local library discovery
            force_reload: Force reloading of the library
        """
        if self._dll is not None and not force_reload:
            return

        self._load_library(data_directory)

    def _load_library(self, data_directory: Optional[Path] = None) -> None:
        """Load the Bruker SDK library with fallback mechanisms.

        Args:
            data_directory: Optional data directory for local library discovery

        Raises:
            SDKError: If library cannot be loaded
        """
        platform_name = PlatformDetector.get_platform()

        if platform_name == "macos":
            raise SDKError("Bruker SDK is not supported on macOS")

        # Get potential library paths
        library_paths = get_dll_paths(data_directory)

        for lib_path in library_paths:
            if not validate_library_path(lib_path):
                continue

            try:
                self._dll = self._load_library_at_path(lib_path)
                self._library_path = lib_path
                logger.debug(f"Successfully loaded Bruker SDK from: {lib_path}")
                return
            except Exception as e:
                logger.debug(f"Failed to load library from {lib_path}: {e}")
                continue

        # If no library found, try generic loading (rely on system PATH)
        try:
            library_name = "timsdata" if platform_name == "windows" else "libtimsdata"
            self._dll = self._load_library_by_name(library_name)
            self._library_path = Path(library_name)
            logger.info(
                f"Successfully loaded Bruker SDK using system PATH: {library_name}"
            )
            return
        except Exception as e:
            logger.debug(f"Failed to load library using system PATH: {e}")

        # If all attempts failed, raise error with helpful message
        lib_name = "timsdata.dll" if platform_name == "windows" else "libtimsdata.so"

        # Get the repository DLL folder path for the error message
        repository_dll_folder = Path(__file__).parent / "dll"
        repository_dll_path = repository_dll_folder / lib_name

        error_msg = (
            f"Failed to load Bruker SDK library '{lib_name}'. "
            f"Checked {len(library_paths)} locations.\n\n"
            f"RECOMMENDED SOLUTION (most reliable):\n"
            f"Place {lib_name} in the repository DLL folder:\n"
            f"  {repository_dll_path}\n\n"
            f"Alternative solutions:\n"
            f"1. Set environment variable: BRUKER_SDK_PATH=<path_to_{lib_name}>\n"
            f"2. Place {lib_name} in the same directory as your data (.d folder)\n"
            f"3. Place {lib_name} in your current working directory\n"
            f"4. Add the SDK directory to your system PATH\n\n"
            f"Where to get the SDK:\n"
            f"- Bruker Daltonics official channels\n"
            f"- Bruker timsTOF software installation\n"
            f"- Contact your Bruker representative\n\n"
            f"Platform: {platform_name}\n"
            f"Repository DLL folder exists: {repository_dll_folder.exists()}"
        )
        raise SDKError(error_msg)

    def _load_library_at_path(self, lib_path: Path) -> CDLL:
        """Load library from a specific path.

        Args:
            lib_path: Path to the library file

        Returns:
            Loaded CDLL object

        Raises:
            OSError: If library cannot be loaded
        """
        if PlatformDetector.is_windows():
            if windll is None:
                raise OSError("windll not available on this platform")
            return windll.LoadLibrary(str(lib_path))
        else:
            return cdll.LoadLibrary(str(lib_path))

    def _load_library_by_name(self, lib_name: str) -> CDLL:
        """Load library by name (using system PATH).

        Args:
            lib_name: Name of the library (without extension)

        Returns:
            Loaded CDLL object

        Raises:
            OSError: If library cannot be loaded
        """
        if PlatformDetector.is_windows():
            if windll is None:
                raise OSError("windll not available on this platform")
            return windll.LoadLibrary(lib_name)
        else:
            return cdll.LoadLibrary(lib_name)

    @property
    def dll(self) -> CDLL:
        """Get the loaded DLL object.

        Returns:
            Loaded CDLL object

        Raises:
            SDKError: If no library is loaded
        """
        if self._dll is None:
            raise SDKError("No Bruker SDK library loaded")
        return self._dll

    @property
    def library_path(self) -> Optional[Path]:
        """Get the path to the loaded library."""
        return self._library_path

    @property
    def is_loaded(self) -> bool:
        """Check if a library is currently loaded."""
        return self._dll is not None

    def reload(self, data_directory: Optional[Path] = None) -> None:
        """Force reload of the library.

        Args:
            data_directory: Optional data directory for local library discovery
        """
        self._dll = None
        self._library_path = None
        self._load_library(data_directory)

    def get_function(self, function_name: str) -> Any:
        """Get a function from the loaded library.

        Args:
            function_name: Name of the function to retrieve

        Returns:
            Function object from the library

        Raises:
            SDKError: If library not loaded or function not found
        """
        if not self.is_loaded:
            raise SDKError("No library loaded")

        try:
            return getattr(self.dll, function_name)
        except AttributeError:
            raise SDKError(f"Function '{function_name}' not found in library")

    def __repr__(self) -> str:
        """String representation of the DLL manager."""
        if self.is_loaded:
            return f"DLLManager(loaded={self.library_path})"
        else:
            return "DLLManager(not loaded)"
