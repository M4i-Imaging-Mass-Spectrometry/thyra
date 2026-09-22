# thyra/readers/__init__.py
"""MSI data readers for various instrument formats.

This package provides reader implementations for different MSI data formats:
- ImzML: Open format for MSI data
- mzPeak: HUPO-PSI Parquet-in-ZIP archive (experimental, read-only)
- Bruker: timsTOF, Rapiflex and solariX data
- Waters: MassLynx .raw imaging data
- PHI: SmartSoft-TOF .raw ToF-SIMS data

Reader organization:
- bruker/: All Bruker formats (timsTOF, Rapiflex, solariX)
- imzml/: ImzML format reader
- mzpeak/: mzPeak archive reader (pyarrow, lazily imported)
- waters/: Waters .raw format reader (native DLL via ctypes)
- phi/: PHI SmartSoft-TOF .raw reader (pure Python, no vendor SDK)

Note that ``.raw`` is claimed by two vendors: Waters stores a directory,
PHI a single file. See :meth:`thyra.core.registry.MSIRegistry.detect_format`.

**This module used to import all five subpackages** so the
``@register_reader`` decorators ran, which is what made every import of
anything under ``thyra`` -- the metadata schema included -- load every
vendor reader (issue #381). The registry now imports a format's module the
first time that format is looked up or detected, from the table in
:mod:`thyra.core.registry`, so nothing here needs to run at import.

``thyra.readers.imzml`` and friends still resolve as attributes, through
:pep:`562`, for code that reached them that way rather than by importing
the submodule. Naming one is what loads it.
"""

import importlib
from typing import TYPE_CHECKING, Any, List

if TYPE_CHECKING:  # pragma: no cover - for type checkers, not at runtime
    from . import bruker, imzml, mzpeak, phi, waters

__all__ = ["bruker", "imzml", "mzpeak", "phi", "waters"]


def __getattr__(name: str) -> Any:
    """Import a reader subpackage the first time it is named."""
    if name in __all__:
        return importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
