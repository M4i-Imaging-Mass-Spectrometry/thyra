"""Convert a parity dataset with Thyra and read back what the metrics need.

The conversion asks for SCiLS Lab's axis explicitly: the same m/z range,
the same bin count, and the axis law fitted from SCiLS's own axis values.
The resampling method follows SCiLS Lab's documented rule (2027a User
Guide, section 4.3.3): nearest neighbour for centroid data, TIC preserving
for profile data. Nothing about the conversion is otherwise changed; the
harness measures Thyra, it does not tune it.

Reading back streams the table's CSC arrays from the Zarr store in chunks,
so a large table is never held in memory as a matrix.
"""

from __future__ import annotations

import hashlib
import json
import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from tests.parity.datasets import ParityDataset
from tests.parity.metrics import nearest_law, spacing_exponent
from tests.parity.scils_reference import ScilsReference

# The axis types SCiLS Lab writes into a saved import (.sip), mapped to
# Thyra's names for the same laws.
_SIP_AXIS_TYPES = {
    "CONSTANT": "constant",
    "LINEAR": "linear_tof",
    "AXIAL_TOF": "linear_tof",
    "OTOF": "reflector_tof",
    "ORBITRAP": "orbitrap",
    "FTICR": "fticr",
    "MRMS": "fticr",
}

# Elements read from a CSC table per step when streaming it.
_CHUNK = 20_000_000


def matched_axis_type(reference: ScilsReference) -> str:
    """Thyra's name for the law SCiLS Lab's axis follows.

    Read from the axis values, which are what SCiLS actually produced: a
    command-line import saves no parameters to read the type from, and a
    saved one is only a statement of intent. :func:`import_axis_type`
    gives the saved statement for the report to set beside it.
    """
    return nearest_law(spacing_exponent(reference.mz))


def import_axis_type(reference: ScilsReference) -> Optional[str]:
    """Thyra's name for the axis type a saved SCiLS import names, if any."""
    return _SIP_AXIS_TYPES.get(reference.import_axis_type.upper())


# Each law's coordinate, in which its grid is uniform, and the inverse.
# Thyra's generators for these laws lay ``target_bins + 1`` uniform edges in
# this coordinate from min_mz to max_mz and store the m/z midpoints, where
# SCiLS Lab's axis values lie on the uniform grid itself, ends included; the
# constant law is absent because Thyra lays its values there the same way.
_LAW_COORDINATES: Dict[
    str, Tuple[Callable[[float], float], Callable[[float], float]]
] = {
    "linear_tof": (lambda m: m**0.5, lambda u: u**2),
    "reflector_tof": (math.log, math.exp),
    "orbitrap": (lambda m: m**-0.5, lambda u: u**-2),
    "fticr": (lambda m: 1.0 / m, lambda u: 1.0 / u),
}


def _widen_by_half_a_bin(
    axis_type: str, lo: float, hi: float, bins: int
) -> Tuple[float, float]:
    """Move ``lo`` and ``hi`` out by half a grid step in the law's coordinate.

    With the range widened so, Thyra's bin midpoints fall on SCiLS Lab's
    axis values, to second order in the step.
    """
    if axis_type not in _LAW_COORDINATES or bins < 2:
        return lo, hi
    forward, inverse = _LAW_COORDINATES[axis_type]
    u_lo, u_hi = forward(lo), forward(hi)
    half = 0.5 * (u_hi - u_lo) / (bins - 1)
    return inverse(u_lo - half), inverse(u_hi + half)


def matched_resampling_config(
    dataset: ParityDataset, reference: ScilsReference
) -> Dict[str, Any]:
    """The ``resampling_config`` that reproduces SCiLS Lab's axis.

    The same range, bin count and law as SCiLS's axis. The manifest can add
    ``"axis_alignment": "scils_values"``, which widens the range by half a
    bin so Thyra's bin midpoints land on SCiLS's values rather than half a
    bin from them (see ``_LAW_COORDINATES``), and ``resampling_overrides``,
    a dict merged in last.

    Args:
        dataset: The dataset; its ``spectrum`` picks the method unless
            the manifest names one.
        reference: The SCiLS reference, whose axis is matched.

    Returns:
        A dict for ``convert_msi(resampling_config=...)``.
    """
    method = dataset.manifest.get("method") or (
        "nearest_neighbor" if dataset.spectrum == "centroid" else "tic_preserving"
    )
    axis_type = matched_axis_type(reference)
    lo, hi, bins = (
        float(reference.mz[0]),
        float(reference.mz[-1]),
        int(reference.mz.size),
    )
    if dataset.manifest.get("axis_alignment") == "scils_values":
        lo, hi = _widen_by_half_a_bin(axis_type, lo, hi, bins)
    config: Dict[str, Any] = {
        "method": method,
        "axis_type": axis_type,
        "min_mz": lo,
        "max_mz": hi,
        "target_bins": bins,
    }
    config.update(dataset.manifest.get("resampling_overrides") or {})
    return config


def _code_fingerprint() -> str:
    """A hash of the installed Thyra package's source files.

    Part of the store key, so a store converted before a change to a
    reader or the resampler is not reused as if the current code had
    written it. Path, size and modification time of every ``.py`` file:
    cheap, and any edit moves it.
    """
    import thyra

    root = Path(thyra.__file__).parent
    digest = hashlib.sha1()
    for path in sorted(root.rglob("*.py")):
        stat = path.stat()
        digest.update(
            f"{path.relative_to(root)}|{stat.st_size}|{stat.st_mtime_ns}".encode()
        )
    return digest.hexdigest()


def _settings_key(config: Dict[str, Any], options: Dict[str, Any]) -> str:
    """A short stable hash naming the store a set of settings and code produce."""
    from thyra import __version__ as thyra_version

    blob = json.dumps(
        {"r": config, "o": options, "v": thyra_version, "c": _code_fingerprint()},
        sort_keys=True,
    )
    return hashlib.sha1(blob.encode("utf-8")).hexdigest()[:10]


def convert(
    dataset: ParityDataset,
    resampling_config: Dict[str, Any],
    reconvert: bool = False,
) -> Path:
    """Convert a dataset with Thyra, reusing an earlier store for the same settings.

    The store lands in ``<dataset>/thyra/<hash>.zarr`` with a
    ``<hash>.json`` beside it recording the settings and the Thyra
    version. A store whose JSON exists is reused unless ``reconvert``.

    Args:
        dataset: The dataset to convert.
        resampling_config: Passed to ``convert_msi`` unchanged.
        reconvert: Convert even when a store for these settings exists.

    Returns:
        The store path.
    """
    from thyra import __version__ as thyra_version
    from thyra import convert_msi

    options: Dict[str, Any] = dict(dataset.reader_options)
    options["(source)"] = str(dataset.source)
    if dataset.pixel_size_um is not None:
        options["(pixel_size_um)"] = dataset.pixel_size_um
    key = _settings_key(resampling_config, options)
    dataset.thyra_dir.mkdir(parents=True, exist_ok=True)
    store = dataset.thyra_dir / f"{key}.zarr"
    record = dataset.thyra_dir / f"{key}.json"
    if record.is_file() and store.is_dir() and not reconvert:
        return store
    if store.exists():
        shutil.rmtree(store)
    record.unlink(missing_ok=True)
    ok = convert_msi(
        str(dataset.source),
        str(store),
        resampling_config=resampling_config,
        reader_options=dataset.reader_options or None,
        pixel_size_um=dataset.pixel_size_um,
        include_optical=False,
    )
    if not ok:
        raise RuntimeError(f"{dataset.name}: Thyra conversion failed")
    record.write_text(
        json.dumps(
            {
                "thyra_version": thyra_version,
                "resampling_config": resampling_config,
                "reader_options": dataset.reader_options,
                "pixel_size_um": dataset.pixel_size_um,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return store


@dataclass(frozen=True)
class ThyraOutput:
    """What the metrics read from one Thyra store.

    Attributes:
        mz: The table's mass axis (``var["mz"]``).
        mean: ``uns["average_spectrum"]``.
        pixel_x, pixel_y: Each table row's raster position (``obs``).
        tic: Each row's sum over the whole axis.
        sum_squares: Each row's sum of squared intensities.
        rows: Selected rows as dense spectra, keyed by row position.
        table: The table's name in the store.
    """

    mz: NDArray[np.float64] = field(repr=False)
    mean: NDArray[np.float64] = field(repr=False)
    pixel_x: NDArray[np.float64] = field(repr=False)
    pixel_y: NDArray[np.float64] = field(repr=False)
    tic: NDArray[np.float64] = field(repr=False)
    sum_squares: NDArray[np.float64] = field(repr=False)
    rows: Dict[int, NDArray[np.float64]] = field(repr=False)
    table: str


# Name endings of the sibling tables a store can hold beside the summed one.
_SIBLING_SUFFIXES = ("_mobility", "_msms")


def _summed_table(store: Path) -> str:
    """The name of the store's summed spectral table.

    Chosen by name: the sibling tables (``*_mobility``, ``*_msms``) can be
    wider than the summed table, so width says nothing. A store with more
    than one summed table -- a volume, one table per slice -- is refused,
    because the comparison is per slice and the harness has no SCiLS
    reference to pair slices with.
    """
    import zarr

    tables = zarr.open_group(str(store / "tables"), mode="r")
    names = sorted(
        name for name, _ in tables.groups() if not name.endswith(_SIBLING_SUFFIXES)
    )
    if len(names) != 1:
        raise RuntimeError(f"{store}: expected one summed table, found {names}")
    return names[0]


def read_output(store: Path, rows: Optional[Iterable[int]] = None) -> ThyraOutput:
    """Read a Thyra store's axis, mean, pixel positions, row sums and some rows.

    Args:
        store: The ``.zarr`` store.
        rows: Row positions to return as dense spectra.

    Returns:
        A :class:`ThyraOutput`.
    """
    import anndata.io as aio
    import zarr

    name = _summed_table(store)
    group = zarr.open_group(str(store / "tables" / name), mode="r")
    obs = aio.read_elem(group["obs"])
    var = aio.read_elem(group["var"])
    mean = np.asarray(aio.read_elem(group["uns/average_spectrum"]), dtype=np.float64)

    x_group = group["X"]
    encoding = dict(x_group.attrs).get("encoding-type")
    n_obs, n_var = (int(v) for v in dict(x_group.attrs)["shape"])
    wanted = sorted(set(int(r) for r in (rows or [])))
    dense = {r: np.zeros(n_var, dtype=np.float64) for r in wanted}
    tic = np.zeros(n_obs, dtype=np.float64)
    sumsq = np.zeros(n_obs, dtype=np.float64)

    data, indices, indptr = x_group["data"], x_group["indices"], x_group["indptr"]
    pointer = np.asarray(indptr[:], dtype=np.int64)
    total = int(pointer[-1])
    wanted_arr = np.asarray(wanted, dtype=np.int64)
    for start in range(0, total, _CHUNK):
        stop = min(start + _CHUNK, total)
        values = np.asarray(data[start:stop], dtype=np.float64)
        inner = np.asarray(indices[start:stop], dtype=np.int64)
        # The outer index of each stored element: the column for CSC, the
        # row for CSR.
        outer = np.searchsorted(pointer, np.arange(start, stop), side="right") - 1
        row_of = inner if encoding == "csc_matrix" else outer
        col_of = outer if encoding == "csc_matrix" else inner
        tic += np.bincount(row_of, weights=values, minlength=n_obs)
        sumsq += np.bincount(row_of, weights=values * values, minlength=n_obs)
        if wanted:
            hit = np.isin(row_of, wanted_arr)
            rows_hit, cols_hit, vals_hit = row_of[hit], col_of[hit], values[hit]
            for r in np.unique(rows_hit):
                own = rows_hit == r
                np.add.at(dense[int(r)], cols_hit[own], vals_hit[own])

    return ThyraOutput(
        mz=np.asarray(var["mz"], dtype=np.float64),
        mean=mean,
        pixel_x=np.asarray(obs["x"], dtype=np.float64),
        pixel_y=np.asarray(obs["y"], dtype=np.float64),
        tic=tic,
        sum_squares=sumsq,
        rows=dense,
        table=name,
    )
