"""Bruker MSI reader implementations.

This package provides readers for Bruker MSI data formats:
- timsTOF: TSF/TDF data via SDK (BrukerReader)
- Rapiflex: MALDI-TOF data via pure Python (RapiflexReader)
- solariX: FT-ICR/MRMS peaks.sqlite data via pure Python (SolarixReader)

Organization:
- timstof/: timsTOF reader and SDK integration
- rapiflex/: Rapiflex reader (pure Python)
- solarix/: solariX reader (pure Python)

Common functionality is provided by BrukerBaseMSIReader and
BrukerFolderStructure for folder analysis.

**The three readers resolve lazily** (:pep:`562`). Format detection reaches
this package for ``BrukerFolderStructure`` on every ``.d`` path -- and on
every plain directory, Waters ones included -- and importing the three
readers there would put the Bruker SDK loader and the timsTOF stack on the
path of a Waters conversion. Naming ``BrukerReader``, ``RapiflexReader`` or
``SolarixReader`` still imports it; nothing else does. See issue #381.
"""

import importlib
from typing import TYPE_CHECKING, Any, List

from ...utils.bruker_exceptions import (
    BrukerReaderError,
    DataError,
    FileFormatError,
    SDKError,
)
from .base_bruker_reader import BrukerBaseMSIReader
from .folder_structure import BrukerFolderInfo, BrukerFolderStructure, BrukerFormat

if TYPE_CHECKING:  # pragma: no cover - for type checkers, not at runtime
    from .rapiflex import RapiflexReader
    from .solarix import SolarixReader
    from .timstof.timstof_reader import BrukerReader

#: Reader name -> the module that defines it, imported on first access.
_LAZY_READERS = {
    "BrukerReader": ".timstof.timstof_reader",
    "RapiflexReader": ".rapiflex",
    "SolarixReader": ".solarix",
}

__all__ = [
    # Base classes
    "BrukerBaseMSIReader",
    "BrukerFolderStructure",
    "BrukerFolderInfo",
    "BrukerFormat",
    # Readers
    "BrukerReader",
    "RapiflexReader",
    "SolarixReader",
    # Exceptions
    "BrukerReaderError",
    "DataError",
    "FileFormatError",
    "SDKError",
]


def __getattr__(name: str) -> Any:
    """Import the module defining a reader the first time it is named."""
    module_path = _LAZY_READERS.get(name)
    if module_path is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module_path, __name__), name)
    globals()[name] = value
    return value


def __dir__() -> List[str]:
    return sorted(set(globals()) | set(__all__))
