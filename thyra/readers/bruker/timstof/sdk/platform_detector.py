"""Platform detection and SDK library path discovery.

This module provides robust platform detection and automatic discovery
of Bruker SDK libraries across different operating systems and
installation paths.
"""

import logging
import platform
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


class PlatformDetector:
    """Detects platform and provides appropriate SDK paths."""

    @staticmethod
    def get_platform() -> str:
        """Get normalized platform identifier.

        Returns:
            Platform string: 'windows', 'linux', or 'macos'
        """
        system = platform.system().lower()
        if system == "darwin":
            return "macos"
        return system

    @staticmethod
    def is_windows() -> bool:
        """Check if running on Windows."""
        return PlatformDetector.get_platform() == "windows"

    @staticmethod
    def is_linux() -> bool:
        """Check if running on Linux."""
        return PlatformDetector.get_platform() == "linux"

    @staticmethod
    def is_macos() -> bool:
        """Check if running on macOS."""
        return PlatformDetector.get_platform() == "macos"


def _get_windows_paths(lib_name: str, data_directory: Optional[Path]) -> List[Path]:
    """Get Windows-specific DLL search paths."""
    paths = []

    # Repository DLL folder (HIGHEST PRIORITY)
    repository_dll_folder = Path(__file__).parent / "dll" / lib_name
    paths.append(repository_dll_folder)
    logger.debug(f"Checking repository DLL folder: {repository_dll_folder}")

    # Current working directory and script location
    paths.extend(
        [
            Path.cwd() / lib_name,
            Path(__file__).parent / lib_name,
            Path(__file__).parent.parent / lib_name,
            Path(__file__).parent.parent.parent / lib_name,
            Path(__file__).parent.parent.parent.parent / lib_name,
        ]
    )

    # Standard installation paths
    paths.extend(
        [
            Path("C:/Program Files/Bruker/timsTOF/sdk") / lib_name,
            Path("C:/Program Files (x86)/Bruker/timsTOF/sdk") / lib_name,
            Path("C:/Bruker/sdk") / lib_name,
            Path("C:/Bruker/timsdata") / lib_name,
        ]
    )

    # User-specific paths
    paths.extend(
        [
            Path.home() / "Desktop" / "MSIConverter" / lib_name,
            Path.home() / "Downloads" / lib_name,
            Path.home() / "Documents" / "Bruker" / lib_name,
        ]
    )

    # Data directory paths
    if data_directory:
        paths.extend(_get_data_directory_paths(lib_name, data_directory))

    # System PATH fallback
    paths.append(Path(lib_name))

    return paths


def _get_data_directory_paths(lib_name: str, data_directory: Path) -> List[Path]:
    """Get paths relative to data directory and network drives."""
    paths = [
        data_directory / lib_name,
        data_directory.parent / lib_name,
    ]

    # Walk up directory tree
    if data_directory.parent.parent:
        paths.append(data_directory.parent.parent / lib_name)
    if data_directory.parent.parent.parent:
        paths.append(data_directory.parent.parent.parent / lib_name)

    # Network drive locations
    try:
        drive = Path(data_directory.anchor)
        if drive not in (Path("C:/"), Path("C:\\")):
            paths.extend(
                [
                    drive / lib_name,
                    drive / "Bruker" / lib_name,
                    drive / "Bruker" / "sdk" / lib_name,
                    drive / "SDK" / lib_name,
                    drive / "timsTOF" / "sdk" / lib_name,
                ]
            )
    except Exception as e:
        logger.debug(f"Could not check network drive locations: {e}")

    return paths


def _get_linux_paths(lib_name: str, data_directory: Optional[Path]) -> List[Path]:
    """Get Linux-specific SO search paths."""
    paths = []

    # Repository DLL folder (HIGHEST PRIORITY)
    repository_dll_folder = Path(__file__).parent / "dll" / lib_name
    paths.append(repository_dll_folder)
    logger.debug(f"Checking repository DLL folder: {repository_dll_folder}")

    # Standard library paths
    paths.extend(
        [
            Path("/usr/lib") / lib_name,
            Path("/usr/local/lib") / lib_name,
            Path("/opt/bruker/lib") / lib_name,
            Path("/usr/lib/x86_64-linux-gnu") / lib_name,
        ]
    )

    # Local data directory
    if data_directory:
        paths.extend(
            [
                data_directory.parent / lib_name,
                data_directory / lib_name,
            ]
        )

    # Current working directory and LD_LIBRARY_PATH
    paths.extend(
        [
            Path.cwd() / lib_name,
            Path(lib_name),
        ]
    )

    return paths


def _get_macos_paths(lib_name: str, data_directory: Optional[Path]) -> List[Path]:
    """Get macOS-specific dylib search paths."""
    paths = []

    # Repository DLL folder (HIGHEST PRIORITY)
    repository_dll_folder = Path(__file__).parent / "dll" / lib_name
    paths.append(repository_dll_folder)
    logger.debug(f"Checking repository DLL folder: {repository_dll_folder}")

    # Standard library paths
    paths.extend(
        [
            Path("/usr/local/lib") / lib_name,
            Path("/opt/bruker/lib") / lib_name,
        ]
    )

    # Local data directory
    if data_directory:
        paths.extend(
            [
                data_directory.parent / lib_name,
                data_directory / lib_name,
            ]
        )

    # Current working directory
    paths.extend(
        [
            Path.cwd() / lib_name,
            Path(lib_name),
        ]
    )

    return paths


def get_dll_paths(data_directory: Optional[Path] = None) -> List[Path]:
    """Get list of potential DLL/SO library paths for the current platform.

    This combines the best path detection logic from timsconvert and imzy
    implementations to provide comprehensive coverage.

    Args:
        data_directory: Optional data directory to check for local libraries

    Returns:
        List of Path objects representing potential library locations

    Note:
        You can set the BRUKER_SDK_PATH environment variable to specify
        a custom location for the Bruker SDK library.
    """
    import os

    paths = []
    platform_name = PlatformDetector.get_platform()

    # Check environment variable first (highest priority)
    sdk_path_env = os.environ.get("BRUKER_SDK_PATH")
    if sdk_path_env:
        sdk_path = Path(sdk_path_env)
        logger.info(f"Using BRUKER_SDK_PATH environment variable: {sdk_path}")
        paths.append(sdk_path)

    # Get platform-specific paths using helper functions
    lib_name = get_library_name()
    if platform_name == "windows":
        paths.extend(_get_windows_paths(lib_name, data_directory))
    elif platform_name == "linux":
        paths.extend(_get_linux_paths(lib_name, data_directory))
    elif platform_name == "macos":
        paths.extend(_get_macos_paths(lib_name, data_directory))

    # Filter to only existing paths and log findings
    existing_paths = []
    checked_count = 0
    for path in paths:
        checked_count += 1
        if path.exists():
            existing_paths.append(path)
            logger.debug(f"Found SDK library: {path}")
        else:
            logger.debug(f"SDK library not found at: {path}")

    if existing_paths:
        logger.debug(
            f"Found {len(existing_paths)} SDK libraries out of {checked_count} checked paths"
        )
    else:
        logger.warning(
            f"No SDK libraries found in {checked_count} standard locations. "
            f"Checked paths include: {[str(p) for p in paths[:5]]}..."
        )

    return existing_paths


def get_library_name() -> str:
    """Get the appropriate library name for the current platform.

    Returns:
        Library filename (e.g., 'timsdata.dll', 'libtimsdata.so')
    """
    platform_name = PlatformDetector.get_platform()

    if platform_name == "windows":
        return "timsdata.dll"
    elif platform_name == "linux":
        return "libtimsdata.so"
    elif platform_name == "macos":
        return "libtimsdata.dylib"
    else:
        raise RuntimeError(f"Unsupported platform: {platform_name}")


def validate_library_path(library_path: Path) -> bool:
    """Validate that a library path is accessible and loadable.

    Args:
        library_path: Path to the library file

    Returns:
        True if the library appears to be valid, False otherwise
    """
    if not library_path.exists():
        return False

    if not library_path.is_file():
        return False

    # Basic size check (libraries should be substantial)
    try:
        size = library_path.stat().st_size
        if size < 1024:  # Less than 1KB is suspicious
            logger.warning(
                f"Library file seems too small: {library_path} ({size} bytes)"
            )
            return False
    except OSError:
        return False

    return True
