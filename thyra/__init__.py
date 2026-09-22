"""MSIConverter - Convert Mass Spectrometry Imaging data to SpatialData/Zarr format.

This package provides tools for converting MSI data from various formats
(ImzML, Bruker) into the modern SpatialData/Zarr format with automatic
pixel size detection.

**Nothing heavy is imported here** (issue #381). This module used to import
the reader and converter packages for their registration side effects and
bind the public API eagerly, which cost 3.5 s and 3,580 modules --
``spatialdata``, ``dask`` and ``anndata`` among them -- before any caller
did any work. A parent package always initialises first, so that was also
the price of ``import thyra.metadata.schema``: the metadata layer depends
on none of it in source (every reference from ``thyra/metadata/`` to the
converter is ``TYPE_CHECKING``-guarded) and could not be used without it.

Two changes replaced it. The registry imports a format's module the first
time that format is looked up, from the table in
:mod:`thyra.core.registry`, so registration no longer needs a front door.
The public names below resolve on first attribute access (:pep:`562`), so
``thyra.convert_msi`` and ``from thyra import preview_msi`` still work and
still cost what they always did -- paid by the caller who asks for them.

``__version__`` stays an ordinary module attribute rather than joining
them. It is what ``python-semantic-release`` rewrites at release
(``version_variables`` in ``pyproject.toml``), and
``thyra.metadata.schema.builder`` reads it through ``from thyra import
__version__`` while building a document -- which is only affordable
because reaching this module is now free.

**Where a broken install now fails.** Issue #310 made ``spatialdata`` and
``defusedxml`` abort ``import thyra``, because both had been treated as
optional and the branch that fired swallowed the cause: the user's first
sign of trouble was ``No converter for format 'spatialdata'. Available:
[]``, naming no package and no reason. The swallowing is still gone and
nothing rewrites either ImportError; what moved is where it surfaces.
``spatialdata`` is imported by :mod:`thyra.converters`, which a conversion
reaches while resolving its output format, before it opens the input;
``defusedxml`` by ``thyra.readers.bruker.mis_parser``, which every Bruker
path reaches during format detection. Each raises at its own import with
its own traceback. Checking for them here instead would have cost the
metadata layer the very thing #381 is for: an implementation that only
writes the document would have to have ``spatialdata`` installed to import
the schema.
"""

# Suppress known warnings from dependencies
import importlib
import warnings
from typing import TYPE_CHECKING, Any, List

if TYPE_CHECKING:  # pragma: no cover - for type checkers, not at runtime
    from .convert import convert_msi
    from .converters.spatialdata.converter import SpatialDataConverter
    from .metadata.document import read_metadata_document
    from .preview import MsiPreview, preview_msi

# Suppress remaining dependency warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="dask")
warnings.filterwarnings("ignore", category=FutureWarning, module="spatialdata")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="numba")
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
warnings.filterwarnings(
    "ignore",
    message="The legacy Dask DataFrame implementation is deprecated",
    category=FutureWarning,
)

__version__ = "4.0.0"

#: Public name -> the module that defines it, imported on first access.
_LAZY_ATTRS = {
    "MsiPreview": ".preview",
    "SpatialDataConverter": ".converters.spatialdata.converter",
    "convert_msi": ".convert",
    "preview_msi": ".preview",
    "read_metadata_document": ".metadata.document",
}

#: Subpackages ``import thyra`` used to bind as attributes on its way to
#: registering their contents. Kept reachable as ``thyra.readers`` and
#: ``thyra.converters`` for code that never imported them by name.
_LAZY_SUBMODULES = ("converters", "readers")

# Expose main API
__all__ = [
    "__version__",
    "convert_msi",
    "preview_msi",
    "read_metadata_document",
    "MsiPreview",
    "SpatialDataConverter",
]


def __getattr__(name: str) -> Any:
    """Resolve a public name by importing the module that defines it."""
    if name in _LAZY_SUBMODULES:
        return importlib.import_module(f".{name}", __name__)
    module_path = _LAZY_ATTRS.get(name)
    if module_path is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module_path, __name__), name)
    globals()[name] = value
    return value


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__) | set(_LAZY_SUBMODULES))
