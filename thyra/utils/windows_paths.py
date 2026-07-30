r"""Extended-length path support for Zarr stores on Windows.

Windows caps a normal path at 260 characters including the terminating
NUL, so 259 usable. A Zarr store is a directory tree and the limit applies
to each key inside it, not to the path the user typed. Thyra's deepest key
is a metadata document about 95 characters below the store root::

    <out>.zarr\tables\<dataset_id>_z0\uns\<section>\<key>\zarr.<32-hex>.partial

so an output path of roughly 165 characters is already enough to push the
deepest key past the limit. It surfaces as::

    Error saving SpatialData: [Errno 2] No such file or directory:
      '...\imzml_version\zarr.<hash>.partial'

which points at a file the user never named and gives no hint that path
length is the problem.

Prefixing the store path with ``\\?\`` opts that path out of the 260
character limit entirely. Verified here on Windows 11 with
``LongPathsEnabled=0``: a 207 character output path fails, and the same
path with the prefix converts cleanly.

The prefix is applied only when the projected deepest key would not fit
otherwise, so ordinary paths are handled exactly as before. It is not
applied to input paths: readers reach vendor SDKs that may not accept
extended-length syntax.
"""

import logging
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

#: Longest usable normal Windows path, excluding the terminating NUL.
WINDOWS_MAX_PATH = 259

#: Characters Thyra needs below the store root for its deepest key, not
#: counting the dataset id. Measured at 95 for the longest section/key
#: pair Thyra emits (``raw_metadata\\max count of pixels x``) plus the
#: ``zarr.<32-hex>.partial`` name Zarr writes through. Metadata key names
#: come from the source file and can be longer, so this carries headroom
#: above the observed worst case.
DEEPEST_KEY_RESERVE = 128

_EXTENDED_PREFIX = "\\\\?\\"


def _long_paths_enabled() -> bool:
    """Whether Windows long-path support is switched on system-wide.

    When it is, normal paths already work past 260 characters and there is
    nothing to work around.
    """
    try:
        import winreg

        with winreg.OpenKey(
            winreg.HKEY_LOCAL_MACHINE,
            r"SYSTEM\CurrentControlSet\Control\FileSystem",
        ) as key:
            value, _ = winreg.QueryValueEx(key, "LongPathsEnabled")
        return bool(value)
    except (ImportError, OSError):
        # Missing key, no permission, or not Windows: assume the limit
        # applies, which only ever means using the prefix when it was not
        # strictly needed.
        return False


def to_extended_length_path(path: Path) -> Path:
    r"""Rewrite an absolute Windows path into extended-length form.

    ``\\?\`` disables all path normalisation, so the input must already be
    absolute and free of ``..`` and forward slashes. Callers should pass a
    ``Path.resolve()`` result.

    Already-extended paths are returned unchanged, and UNC shares get the
    ``\\?\UNC\`` spelling rather than an invalid double prefix.
    """
    text = str(path)

    if text.startswith(_EXTENDED_PREFIX):
        return path

    if text.startswith("\\\\"):
        # \\server\share -> \\?\UNC\server\share
        return Path(f"{_EXTENDED_PREFIX}UNC\\{text.lstrip(chr(92))}")

    return Path(f"{_EXTENDED_PREFIX}{text}")


def projected_deepest_key_length(output_path: Path, dataset_id: str) -> int:
    """Length of the longest path Thyra expects to write inside the store."""
    return len(str(output_path)) + 1 + len(dataset_id) + DEEPEST_KEY_RESERVE


def prepare_zarr_output_path(output_path: Path, dataset_id: str) -> Path:
    """Return the path to hand to Zarr, extended if the store needs it.

    A no-op off Windows, and a no-op on Windows when long-path support is
    enabled or the projected deepest key fits inside the normal limit.

    Args:
        output_path: The resolved, absolute output store path.
        dataset_id: The dataset identifier, which appears in the deepest key.
    """
    if sys.platform != "win32":
        return output_path

    projected = projected_deepest_key_length(output_path, dataset_id)
    if projected <= WINDOWS_MAX_PATH:
        return output_path

    if _long_paths_enabled():
        logger.debug(
            "Output path is long (projected deepest key %d characters) but "
            "Windows long-path support is enabled; using it as given",
            projected,
        )
        return output_path

    extended = to_extended_length_path(output_path)
    logger.info(
        "Output path is long enough that the deepest Zarr key would exceed "
        "the %d character Windows limit (projected %d). Writing through an "
        "extended-length path so the conversion is not truncated.",
        WINDOWS_MAX_PATH,
        projected,
    )
    return extended
