"""Read single pieces of a converted store, never the whole of it.

A table of spectra can hold billions of values, and a picture needs a
sliver of them: one block of columns for an ion image, the pixel positions,
the mass axis. Everything here opens the Zarr arrays directly and slices
them, so the cost follows what is asked for and not the size of the store.

Tables are told apart by their ``var`` columns, as
``docs/for-tool-developers.md`` prescribes, and never by their names: the
dataset id is a user setting.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..metadata.schema.models import MSI_METADATA_UNS_KEY
from ..metadata.schema.store_io import subgroup
from ..utils.windows_paths import prepare_zarr_read_path

if TYPE_CHECKING:
    import zarr

#: The three kinds of table a store can hold, and anything else.
MSI = "msi"
MOBILITY = "mobility"
MSMS = "msms"
OTHER = "other"


class QuicklookError(Exception):
    """A reason no picture can be made, written as one line for the user."""


def table_kind(var_columns: Any) -> str:
    """Which kind of table has these ``var`` columns.

    Args:
        var_columns: The names of a table's ``var`` columns.

    Returns:
        ``"mobility"``, ``"msms"``, ``"msi"``, or ``"other"`` for a table
        with no ``mz`` column, which no Thyra version writes.
    """
    columns = set(var_columns)
    if "mobility" in columns:
        return MOBILITY
    if "precursor_mz" in columns:
        return MSMS
    if "mz" in columns:
        return MSI
    return OTHER


@dataclass(frozen=True)
class TableSummary:
    """One table of a store: its name, kind and shape.

    Attributes:
        name: The table's element name.
        kind: ``"msi"``, ``"mobility"``, ``"msms"`` or ``"other"``.
        shape: ``(rows, columns)``, or ``None`` when the table states none.
    """

    name: str
    kind: str
    shape: Optional[Tuple[int, int]]


@dataclass(frozen=True)
class Positions:
    """Where each row of a table sits on the raster.

    Attributes:
        x: Raster column of each row.
        y: Raster row of each row.
        z: Slice of each row, for a table converted as a volume; else
            ``None``.
    """

    x: NDArray[np.int64]
    y: NDArray[np.int64]
    z: Optional[NDArray[np.int64]]


def _array(group: "zarr.Group", name: str, where: str) -> "zarr.Array[Any]":
    """The array ``name`` of ``group``, or a refusal naming what is missing."""
    import zarr

    node = group.get(name)
    if not isinstance(node, zarr.Array):
        raise QuicklookError(f"{where} has no '{name}' array; the store is incomplete.")
    return node


def _column(group: "zarr.Group", name: str, where: str) -> NDArray[Any]:
    """One column of an ``obs`` or ``var`` group as an array.

    A plain column is a Zarr array. One that anndata wrote as a group, a
    nullable integer for instance, is decoded by anndata.
    """
    import zarr

    node = group.get(name)
    if node is None:
        raise QuicklookError(f"{where} has no '{name}' column.")
    if isinstance(node, zarr.Array):
        return np.asarray(node[:])
    import anndata as ad

    return np.asarray(ad.io.read_elem(node))


def _table_shape(table: "zarr.Group") -> Optional[Tuple[int, int]]:
    """``(rows, columns)`` of a table's matrix, without reading it."""
    import zarr

    matrix = table.get("X")
    shape: Any = None
    if isinstance(matrix, zarr.Array):
        shape = matrix.shape
    elif isinstance(matrix, zarr.Group):
        shape = matrix.attrs.get("shape")
    if not isinstance(shape, (list, tuple)) or len(shape) != 2:
        return None
    return int(shape[0]), int(shape[1])


def _natural_order(name: str) -> List[Tuple[int, int, str]]:
    """A sort key that reads the digits in a name as numbers."""
    parts = re.split(r"(\d+)", name)
    return [(0, int(part), "") if part.isdigit() else (1, 0, part) for part in parts]


class MsiTable:
    """One table of spectra, read a piece at a time.

    Attributes:
        name: The table's element name.
        mz: The mass axis, one value per column.
        n_rows: Number of pixels.
        n_columns: Number of m/z bins.
    """

    def __init__(self, name: str, group: "zarr.Group") -> None:
        """Open the table and read its mass axis.

        Args:
            name: The table's element name.
            group: The table's Zarr group.

        Raises:
            QuicklookError: If the table lacks a part a picture needs, or
                its mass axis is not strictly increasing.
        """
        self.name = name
        self._group = group
        where = f"Table {name}"
        var = subgroup(group, "var")
        if var is None:
            raise QuicklookError(
                f"{where} has no 'var' group; the store is incomplete."
            )
        self.mz: NDArray[np.float64] = _column(var, "mz", where).astype(np.float64)
        shape = _table_shape(group)
        if shape is None:
            raise QuicklookError(f"{where} states no shape; the store is incomplete.")
        self.n_rows, self.n_columns = shape
        if self.mz.size != self.n_columns:
            raise QuicklookError(
                f"{where} has {self.n_columns} columns and {self.mz.size} m/z "
                "values; run 'thyra validate' on the store."
            )
        if self.mz.size > 1 and not bool(np.all(np.diff(self.mz) > 0)):
            raise QuicklookError(
                f"{where} has a mass axis that is not strictly increasing; "
                "run 'thyra validate' on the store."
            )

    def positions(self) -> Positions:
        """The raster position of every row.

        Returns:
            The ``x`` and ``y`` columns of ``obs``, and ``z`` when the table
            is a volume.
        """
        where = f"Table {self.name}"
        obs = subgroup(self._group, "obs")
        if obs is None:
            raise QuicklookError(
                f"{where} has no 'obs' group; the store is incomplete."
            )
        x = _column(obs, "x", where).astype(np.int64)
        y = _column(obs, "y", where).astype(np.int64)
        z = _column(obs, "z", where).astype(np.int64) if "z" in obs else None
        if x.size != self.n_rows or y.size != self.n_rows:
            raise QuicklookError(
                f"{where} has {self.n_rows} rows and {x.size} pixel positions; "
                "the store is incomplete."
            )
        return Positions(x, y, z)

    def window_sum(self, first: int, stop: int) -> NDArray[np.float64]:
        """Per-pixel sum of the columns ``first`` up to, not including, ``stop``.

        Only the values of those columns are read: the matrix is stored
        column by column, so a block of columns is one contiguous run of
        its ``data`` and ``indices`` arrays.

        Args:
            first: First column of the block.
            stop: One past its last column.

        Returns:
            One value per row of the table.

        Raises:
            QuicklookError: If the matrix is not stored column by column,
                which every current Thyra version does.
        """
        where = f"Table {self.name}"
        matrix = subgroup(self._group, "X")
        encoding = matrix.attrs.get("encoding-type") if matrix is not None else None
        if matrix is None or encoding != "csc_matrix":
            raise QuicklookError(
                f"{where} does not store its matrix column by column (CSC), "
                "which this command reads. Convert the source again with a "
                "current Thyra version."
            )
        pointers = np.asarray(_array(matrix, "indptr", where)[first : stop + 1])
        start, end = int(pointers[0]), int(pointers[-1])
        rows = np.asarray(_array(matrix, "indices", where)[start:end])
        values = np.asarray(_array(matrix, "data", where)[start:end], dtype=np.float64)
        if rows.size and (int(rows.min()) < 0 or int(rows.max()) >= self.n_rows):
            raise QuicklookError(
                f"{where} names a row outside the table; the store is damaged."
            )
        return np.bincount(rows, weights=values, minlength=self.n_rows)

    def average_spectrum(self) -> Optional[NDArray[np.float64]]:
        """The stored mean spectrum, or ``None`` when the table has none.

        A value that is not a finite number, which no converter writes, is
        read as zero, so that one damaged bin does not blank the drawing.
        """
        import zarr

        uns = subgroup(self._group, "uns")
        node = uns.get("average_spectrum") if uns is not None else None
        if not isinstance(node, zarr.Array):
            return None
        spectrum = np.asarray(node[:], dtype=np.float64).ravel()
        if spectrum.size != self.n_columns:
            return None
        return np.nan_to_num(spectrum, nan=0.0, posinf=0.0, neginf=0.0)

    def metadata_block(self) -> Dict[str, Any]:
        """The table's ``msi_metadata`` block, or an empty mapping.

        Read from the group that is already open. A block that cannot be
        read gives an empty mapping: the page is still worth having.
        """
        uns = subgroup(self._group, "uns")
        node = uns.get(MSI_METADATA_UNS_KEY) if uns is not None else None
        if node is None:
            return {}
        import anndata as ad

        try:
            block: Any = ad.io.read_elem(node)
        except Exception:  # whatever anndata raises for a damaged block
            return {}
        return dict(block) if isinstance(block, Mapping) else {}


class Store:
    """A converted store, opened for reading.

    Attributes:
        path: The path the user gave.
        tables: Every table of the store, sorted by name, a number inside a
            name by its value: slice 2 before slice 10.
    """

    def __init__(self, path: Path, root: "zarr.Group", tables: "zarr.Group") -> None:
        """Wrap an opened store; use :func:`open_store` to open one."""
        self.path = path
        self._root = root
        self._tables = tables
        self.tables: List[TableSummary] = [
            self._summary(name) for name in sorted(tables.keys(), key=_natural_order)
        ]

    def _summary(self, name: str) -> TableSummary:
        table = subgroup(self._tables, name)
        var = subgroup(table, "var") if table is not None else None
        kind = table_kind(var.keys()) if var is not None else OTHER
        shape = _table_shape(table) if table is not None else None
        return TableSummary(name, kind, shape)

    @property
    def name(self) -> str:
        """The store's folder name without its ``.zarr`` ending."""
        return store_name(self.path)

    @property
    def root_attrs(self) -> Dict[str, Any]:
        """The store's root attributes."""
        return dict(self._root.attrs)

    def msi_tables(self) -> List[TableSummary]:
        """The tables of spectra: one per slice, or one for a volume.

        Raises:
            QuicklookError: If the store holds none.
        """
        found = [table for table in self.tables if table.kind == MSI]
        if not found:
            raise QuicklookError(
                f"{self.path} holds no table of spectra, so there is nothing "
                "to draw. Was it converted by Thyra?"
            )
        return found

    def open_table(self, name: str) -> MsiTable:
        """Open one table of spectra by name."""
        group = subgroup(self._tables, name)
        if group is None:
            raise QuicklookError(f"{self.path} has no table named {name}.")
        return MsiTable(name, group)

    def image_names(self) -> List[str]:
        """The names of the store's image elements, sorted."""
        images = subgroup(self._root, "images")
        return sorted(images.keys()) if images is not None else []

    def _image_levels(self, name: str) -> List["zarr.Array[Any]"]:
        """The arrays of one image, finest first.

        The level names come from the image's OME metadata. A group whose
        metadata cannot be read gives its arrays in name order instead.
        """
        import zarr

        group = subgroup(self._root, f"images/{name}")
        if group is None:
            return []
        try:
            ome: Any = group.attrs.get("ome") or dict(group.attrs)
            paths = [str(level["path"]) for level in ome["multiscales"][0]["datasets"]]
        except (KeyError, IndexError, TypeError):
            paths = sorted(group.keys())
        levels = [group.get(path) for path in paths]
        return [level for level in levels if isinstance(level, zarr.Array)]

    def tic_image(self, table: str) -> Optional[NDArray[np.float64]]:
        """The TIC image of one table, or ``None`` when the store has none.

        Returns:
            ``(y, x)`` for a slice and ``(z, y, x)`` for a volume: the
            image's single channel.
        """
        levels = self._image_levels(f"{table}_tic")
        if not levels or levels[0].ndim not in (3, 4):
            return None
        return np.asarray(levels[0][0], dtype=np.float64)

    def image_preview(self, name: str, long_side: int) -> Optional[NDArray[Any]]:
        """A reduced copy of one image, a little larger than ``long_side``.

        The coarsest stored level that is still large enough is read, and
        every n-th pixel of it, so a slide scan of a gigabyte costs what
        the preview holds.

        Args:
            name: The image element.
            long_side: The size the caller will reduce the result to.

        Returns:
            ``(channels, y, x)``, or ``None`` for an image that is not 2D.
        """
        levels = [level for level in self._image_levels(name) if level.ndim == 3]
        if not levels:
            return None
        chosen = levels[0]
        for level in levels[1:]:
            if max(level.shape[1:]) >= long_side:
                chosen = level
        step = max(1, max(chosen.shape[1:]) // long_side)
        return np.asarray(chosen[:, ::step, ::step])


def store_name(path: Path) -> str:
    """A store's folder name without its ``.zarr`` ending."""
    folder = Path(str(path).rstrip("\\/")).name
    return folder[:-5] if folder.lower().endswith(".zarr") else folder


def open_store(path: Path, read_path: Optional[Path] = None) -> Store:
    """Open a converted store for reading.

    Args:
        path: The store's ``.zarr`` folder, as the user gave it.
        read_path: The same folder as ``prepare_zarr_read_path`` returned
            it, when the caller has prepared it already. Preparing a path
            walks the store on Windows, so it is done once.

    Returns:
        The opened store.

    Raises:
        QuicklookError: If the path holds no Zarr group, or a Zarr group
            that is not a SpatialData store.
    """
    import zarr

    if read_path is None:
        read_path = prepare_zarr_read_path(path)
    try:
        root = zarr.open_group(str(read_path), mode="r")
    except Exception as exc:  # zarr raises several unrelated types here
        raise QuicklookError(
            f"{path} is not a Thyra store: no Zarr group can be opened there."
        ) from exc
    tables = subgroup(root, "tables")
    if tables is None:
        raise QuicklookError(f"{path} is not a Thyra store: it has no 'tables' group.")
    return Store(path, root, tables)
