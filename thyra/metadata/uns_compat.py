# thyra/metadata/uns_compat.py
"""Read-side defense for ``uns`` blocks written by earlier Thyra versions.

Stores written before the JSON rule in
``BaseSpatialDataConverter._collect_optional_sections`` carry string
lists in ``uns`` (notably ``raw_metadata["cvParams"]``) that materialize
as numpy ``StringDType`` arrays when read back through AnnData/zarr.
``copy.deepcopy`` of such an array segfaults the process on numpy
2.1.x-2.2.x (numpy#28609, fixed in 2.3.0) -- and every table copy
deepcopies ``uns``: ``AnnData.copy``, ``spatialdata.polygon_query``,
``bounding_box_query``, and the join paths behind them.

New stores no longer contain these arrays. This module exists for
already-written stores read in an environment pinned to the affected
numpy range.
"""

from collections.abc import Mapping
from typing import Any, Dict

import numpy as np

__all__ = ["sanitize_uns_string_arrays"]


def sanitize_uns_string_arrays(uns: Mapping[str, Any]) -> Dict[str, Any]:
    """Return ``uns`` with numpy string arrays replaced by plain lists.

    Recurses through nested dicts and lists; every other value is kept
    as-is, so the result is safe to assign straight back::

        import spatialdata as sd
        from thyra.metadata import sanitize_uns_string_arrays

        sdata = sd.read_zarr("store.zarr")
        for table in sdata.tables.values():
            table.uns = sanitize_uns_string_arrays(table.uns)

    Needed only for stores written by Thyra versions that stored string
    lists in ``uns`` as arrays, read under numpy 2.1.x-2.2.x, where
    deepcopying a ``StringDType`` array kills the process (see module
    docstring). Calling it anywhere else is harmless.
    """
    return {key: _sanitize(value) for key, value in uns.items()}


def _sanitize(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _sanitize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_sanitize(item) for item in value]
    # ``kind == "T"`` is StringDType, the one dtype whose deepcopy is
    # broken. On numpy builds without StringDType no array has this
    # kind, so the walk is a no-op there.
    if isinstance(value, np.ndarray) and value.dtype.kind == "T":
        return value.tolist()
    return value
