"""Stores built by hand, so each test states exactly what its store holds.

A converted store is a Zarr group with a ``tables`` group of AnnData
tables and an ``images`` group. Writing one directly takes milliseconds
and lets a test choose the pixel positions, the table kinds and the
metadata, which a conversion would decide for it. The tests that need a
real conversion use a committed fixture instead.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import zarr
from scipy import sparse

#: The mass axis of the default table: ten bins, one m/z apart.
MZ = np.arange(500.0, 510.0)

#: An L-shaped region whose positions start far from zero: the corner at
#: column 12, row 7 was never measured.
L_SHAPE = [(10, 5), (11, 5), (12, 5), (10, 6), (11, 6), (12, 6), (10, 7), (11, 7)]


def dense_matrix(n_rows: int, n_columns: int = MZ.size) -> np.ndarray:
    """Intensities no two cells share: row ``r``, column ``c`` holds ``100 r + c + 1``.

    Column 3 is left empty, so a window over it alone has no signal.
    """
    rows, columns = np.indices((n_rows, n_columns))
    values = 100.0 * rows + columns + 1.0
    values[:, 3] = 0.0
    return values


def make_table(
    positions: Sequence[Sequence[int]] = tuple(L_SHAPE),
    mz: np.ndarray = MZ,
    extra_var: Optional[Dict[str, Any]] = None,
    uns: Optional[Dict[str, Any]] = None,
    matrix: Optional[np.ndarray] = None,
    csr: bool = False,
) -> ad.AnnData:
    """One table: a row per position, a column per ``mz`` value."""
    coordinates = np.asarray(positions, dtype=np.int64)
    dense = dense_matrix(len(coordinates), len(mz)) if matrix is None else matrix
    obs = pd.DataFrame(
        {"x": coordinates[:, 0], "y": coordinates[:, 1]},
        index=[f"pixel_{i}" for i in range(len(coordinates))],
    )
    if coordinates.shape[1] == 3:
        obs["z"] = coordinates[:, 2]
    var = pd.DataFrame({"mz": mz}, index=[f"bin_{i}" for i in range(len(mz))])
    for column, values in (extra_var or {}).items():
        var[column] = values
    stored = dict(uns or {})
    stored.setdefault("average_spectrum", dense.mean(axis=0))
    layout = sparse.csr_matrix if csr else sparse.csc_matrix
    return ad.AnnData(X=layout(dense), obs=obs, var=var, uns=stored)


def write_store(
    path: Path,
    tables: Dict[str, ad.AnnData],
    images: Optional[Dict[str, Sequence[np.ndarray]]] = None,
    attrs: Optional[Dict[str, Any]] = None,
) -> Path:
    """Write ``tables`` and ``images`` as a store at ``path``.

    Each image is a list of levels, finest first, named as its OME
    metadata says and not ``s0``, ``s1``.
    """
    root = zarr.open_group(str(path), mode="w")
    root.attrs.update(attrs or {})
    group = root.create_group("tables")
    for name, table in tables.items():
        ad.io.write_elem(group, name, table)
    image_group = root.create_group("images")
    for name, levels in (images or {}).items():
        element = image_group.create_group(name)
        paths = [f"level{index}" for index in range(len(levels))]
        element.attrs["ome"] = {
            "multiscales": [{"datasets": [{"path": path_} for path_ in paths]}]
        }
        for path_, level in zip(paths, levels):
            element.create_array(path_, data=np.asarray(level))
    return path
