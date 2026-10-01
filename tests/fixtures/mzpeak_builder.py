"""Construct minimal valid ``.mzpeak`` archives for tests.

There is no mzPeak writer in Thyra and no Rust toolchain in CI, so the test
corpus is built here, in pure Python, from ``zipfile`` and ``pyarrow``.

Everything this module emits is written out literally -- the index JSON is a
dict spelled out below, the column names are string constants, the member
names are hardcoded -- and nothing is derived from
:mod:`thyra.readers.mzpeak`. That separation is the point. The hand-authored
imzML corpus in ``tests/data/fixtures`` exists because pyimzml's parser and
writer "agree on each other's mistakes"; a builder that asked the reader how
to spell things would reintroduce exactly that loop, and every test on top of
it would pass while real archives failed.

The shapes here were measured against the reference implementation
(HUPO-PSI/mzPeak @ 502c3a4) and its shipped sample archives, not inferred
from the draft prose. The chunked layout was measured against archives
written by mzpeak-convert 0.14.0 (okohlbacher/mzPeakConverter @ 0ed311e),
and its encoders follow the specification's own text
(``docs/layouts/chunked-layout.md`` @ ecc8062).
"""

from __future__ import annotations

import json
import struct
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

# Member names as the reference converter writes them. Real readers resolve
# these through the index rather than by name; the fixtures use the reference
# spelling so a reader that (wrongly) hardcodes them still sees valid input
# and its bug shows up somewhere more informative than a missing file.
DATA_MEMBER = "spectra_data.parquet"
PEAKS_MEMBER = "spectra_peaks.parquet"
METADATA_MEMBER = "spectra_metadata.parquet"
SCANS_MEMBER = "spectra_metadata_scans.parquet"
INDEX_MEMBER = "mzpeak_index.json"

POSITION_X_COLUMN = "opt_IMS_1000050_position_x"
POSITION_Y_COLUMN = "opt_IMS_1000051_position_y"
#: How mzpeak-convert 0.16.0 names the position z column.
POSITION_Z_COLUMN = "position_z"

#: Micrometre unit accession, as the reference archive declares it.
UNIT_MICROMETRE = "UO:0000017"

#: Chunk encodings, by the accession a row names in ``chunk_encoding``.
PLAIN = "MS:1000576"
DELTA = "MS:1003089"
NUMPRESS_LINEAR = "MS:1002312"
GRID = "MS:1003826"

#: Grid models.
LINEAR_GRID = "MS:1003824"
SQUARE_ROOT_GRID = "MS:1003825"

#: Intercept, slope and scale of the grid the fixtures are written on. The
#: slope is a power of two, so every m/z with at most ten binary places
#: comes back from its index to the bit.
GRID_PARAMETERS = (0.0, 2.0**-10, 1.0)

#: A stored point: its m/z and intensity, both ``None`` on padding.
Point = Tuple[Optional[float], Optional[float]]


class Spectrum:
    """One spectrum destined for a fixture archive.

    Attributes:
        x: Position along x, in the file's own (typically 1-based) frame.
            ``None`` is written as a null: a scan that belongs to no pixel.
        y: Position along y, or ``None``.
        mzs: m/z values, ascending.
        intensities: Intensity values, same length as ``mzs``.
    """

    def __init__(
        self,
        x: Optional[int],
        y: Optional[int],
        mzs: Sequence[float],
        intensities: Sequence[float],
    ):
        """Store one spectrum's coordinates and arrays."""
        self.x = None if x is None else int(x)
        self.y = None if y is None else int(y)
        self.mzs = np.asarray(mzs, dtype=np.float64)
        self.intensities = np.asarray(intensities, dtype=np.float64)
        if self.mzs.size != self.intensities.size:
            raise ValueError("mzs and intensities must be the same length")


def _point_table(
    spectra: Sequence[Optional[Spectrum]],
    null_pair_after: Optional[int] = None,
) -> pa.Table:
    """Build the point-layout signal table.

    Args:
        spectra: Spectra in ``spectrum_index`` order. ``None`` stands for a
            spectrum with no points in this member: it keeps its index and
            writes no rows.
        null_pair_after: When given, insert a null pair (two rows whose m/z
            and intensity are both null) after this many real points of every
            spectrum, imitating the padding the reference converter writes
            where it has removed a run of zeros.

    Returns:
        A table with the single struct column ``point``.
    """
    indices: List[int] = []
    mzs: List[Optional[float]] = []
    intensities: List[Optional[float]] = []

    for index, spectrum in enumerate(spectra):
        if spectrum is None:
            continue
        for position, (mz, intensity) in enumerate(
            zip(spectrum.mzs, spectrum.intensities)
        ):
            indices.append(index)
            mzs.append(float(mz))
            intensities.append(float(intensity))
            if null_pair_after is not None and position == null_pair_after - 1:
                # A null *pair* -- the reference writer never emits a lone
                # null, and its reader raises on an unpaired one.
                for _ in range(2):
                    indices.append(index)
                    mzs.append(None)
                    intensities.append(None)

    point = pa.StructArray.from_arrays(
        [
            pa.array(indices, type=pa.uint64()),
            pa.array(mzs, type=pa.float64()),
            pa.array(intensities, type=pa.float32()),
        ],
        fields=[
            pa.field("spectrum_index", pa.uint64()),
            pa.field("mz", pa.float64()),
            pa.field("intensity", pa.float32()),
        ],
    )
    return pa.table({"point": point})


def _stored_points(spectrum: Spectrum, null_pair_after: Optional[int]) -> List[Point]:
    """One spectrum's points as they are stored, padding included."""
    points: List[Point] = []
    for position, (mz, intensity) in enumerate(zip(spectrum.mzs, spectrum.intensities)):
        points.append((float(mz), float(intensity)))
        if null_pair_after is not None and position == null_pair_after - 1:
            points += [(None, None), (None, None)]
    return points


def delta_encode(values: Sequence[Optional[float]]) -> List[Optional[float]]:
    """Delta-encode with nulls, as the specification's example does.

    The first value is left out unless it is null. A value after a null is
    written as it is.
    """
    encoded: List[Optional[float]] = []
    last = values[0]
    if last is None:
        encoded.append(None)
    for value in values[1:]:
        if value is None or last is None:
            encoded.append(value)
        else:
            encoded.append(value - last)
        last = value
    return encoded


def _half_bytes_of(residual: int) -> List[int]:
    """One MS-Numpress residual: a count, then the half bytes it keeps.

    The count is how many of the eight half bytes are left out at the top:
    zeros for a value that fits, ones (count plus 8) for a negative one.
    """
    pattern = residual & 0xFFFFFFFF
    digits = [(pattern >> (4 * place)) & 0xF for place in range(8)]
    top = digits[7]
    if top not in (0x0, 0xF):
        return [0, *digits]
    kept = 8
    while kept > 0 and digits[kept - 1] == top:
        kept -= 1
    if top == 0xF and kept == 0:
        # Minus one keeps one half byte; the count cannot say eight ones.
        kept = 1
    left_out = 8 - kept
    return [left_out + (8 if top == 0xF else 0), *digits[:kept]]


def numpress_linear_encode(
    values: Sequence[float], fixed_point: Optional[float] = None
) -> bytes:
    """Encode values with MS-Numpress linear prediction.

    Written from the published description (Teleman et al., Mol Cell
    Proteomics 2014, 13, 1537): the scale as a big-endian double, the first
    two values as little-endian 32-bit integers, then for each further
    value what the prediction ``2 a - b`` missed.

    Args:
        values: The values to encode.
        fixed_point: The scale. Defaults to the largest power of two that
            keeps every value and every residual inside 31 bits, so that
            values with few binary places come back to the bit.
    """

    def scaled_by(scale: float) -> Tuple[List[int], List[int]]:
        whole = [int(value * scale + 0.5) for value in values]
        missed = [
            whole[place] - (2 * whole[place - 1] - whole[place - 2])
            for place in range(2, len(whole))
        ]
        return whole, missed

    if fixed_point is None:
        largest = max([abs(v) for v in values] + [1.0])
        fixed_point = 2.0 ** int(np.floor(np.log2(0x7FFFFFFF / largest)))
        while any(abs(r) > 0x7FFFFFFF for r in scaled_by(fixed_point)[1]):
            fixed_point /= 2.0
    encoded = bytearray(struct.pack(">d", fixed_point))
    scaled, residuals = scaled_by(fixed_point)
    for value in scaled[:2]:
        encoded += struct.pack("<I", value)
    halves: List[int] = []
    for residual in residuals:
        halves += _half_bytes_of(residual)
    if len(halves) % 2:
        halves.append(0)
    for high, low in zip(halves[0::2], halves[1::2]):
        encoded.append((high << 4) | low)
    return bytes(encoded)


def grid_indices(
    values: Sequence[float], grid_type: str, parameters: Sequence[float]
) -> List[int]:
    """Grid indices of m/z values, each as the step from the one before."""
    intercept, slope, scale = parameters
    scaled = np.asarray(values, dtype=np.float64) * scale
    if grid_type == SQUARE_ROOT_GRID:
        scaled = np.sqrt(scaled)
    indices = np.floor((scaled - intercept) / slope + 0.5).astype(np.int64)
    return [int(v) for v in np.diff(indices, prepend=0)]


def _chunk_row(
    index: int,
    points: Sequence[Point],
    encoding: str,
    grid_type: str,
    fixed_point: Optional[float] = None,
) -> Dict[str, Any]:
    """Write one chunk row, spelled out column by column."""
    mzs = [mz for mz, _ in points]
    real = [mz for mz in mzs if mz is not None]
    row: Dict[str, Any] = {
        "spectrum_index": index,
        # Both bounds are zero on a chunk that holds padding alone.
        "mz_chunk_start": real[0] if real else 0.0,
        "mz_chunk_end": real[-1] if real else 0.0,
        "mz_chunk_values": None,
        "chunk_encoding": encoding,
        "intensity": [intensity for _, intensity in points],
        "mz_numpress_linear_bytes": None,
        "mz_grid": None,
    }
    if encoding == PLAIN:
        row["mz_chunk_values"] = mzs[1:]
    elif encoding == DELTA:
        row["mz_chunk_values"] = delta_encode(mzs)
    elif encoding == NUMPRESS_LINEAR:
        # MS-Numpress has no null; the reference writer stores zero.
        filled = [0.0 if mz is None else mz for mz in mzs]
        row["mz_numpress_linear_bytes"] = list(
            numpress_linear_encode(filled, fixed_point)
        )
    elif encoding == GRID:
        if len(real) != len(mzs):
            raise ValueError("a grid has no null; do not pad a grid fixture")
        row["mz_grid"] = {
            "grid_type": grid_type,
            "parameters": list(GRID_PARAMETERS),
            "indices": grid_indices(mzs, grid_type, GRID_PARAMETERS),
        }
    else:
        # An encoding the builder cannot write is still named in the row,
        # so that a reader can be shown refusing it.
        row["mz_chunk_values"] = mzs[1:]
    return row


def chunk_table(
    content: Sequence[Optional[Spectrum]],
    null_pair_after: Optional[int] = None,
    encoding: Union[str, Sequence[str]] = DELTA,
    chunk_points: int = 4,
    grid_type: str = LINEAR_GRID,
    fixed_point: Optional[float] = None,
) -> pa.Table:
    """Build the chunked-layout signal table.

    Args:
        content: Spectra in ``spectrum_index`` order, ``None`` for one with
            no points in this member.
        null_pair_after: As for the point layout.
        encoding: Accession written into ``chunk_encoding``. A sequence is
            dealt out row by row, for a member that mixes encodings.
        chunk_points: Stored points per row. The specification lets a
            writer cut where it likes, as long as the chunks ascend.
        grid_type: Grid model of a grid-encoded row.
        fixed_point: Scale of an MS-Numpress row; see
            :func:`numpress_linear_encode`.

    Returns:
        A table with the single struct column ``chunk``. It carries the
        byte column and the grid column only when a row needs them, which
        is what the converter does.
    """
    encodings = [encoding] if isinstance(encoding, str) else list(encoding)
    rows: List[Dict[str, Any]] = []
    for index, spectrum in enumerate(content):
        if spectrum is None:
            continue
        points = _stored_points(spectrum, null_pair_after)
        for start in range(0, len(points), chunk_points):
            rows.append(
                _chunk_row(
                    index,
                    points[start : start + chunk_points],
                    encodings[len(rows) % len(encodings)],
                    grid_type,
                    fixed_point,
                )
            )

    fields = [
        pa.field("spectrum_index", pa.uint64()),
        pa.field("mz_chunk_start", pa.float64()),
        pa.field("mz_chunk_end", pa.float64()),
        pa.field("mz_chunk_values", pa.large_list(pa.float64())),
        pa.field("chunk_encoding", pa.string()),
        pa.field("intensity", pa.large_list(pa.float32())),
    ]
    if NUMPRESS_LINEAR in encodings:
        fields.append(pa.field("mz_numpress_linear_bytes", pa.large_list(pa.uint8())))
    if GRID in encodings:
        grid = pa.struct(
            [
                pa.field("grid_type", pa.large_string()),
                pa.field("parameters", pa.large_list(pa.float64())),
                pa.field("indices", pa.large_list(pa.uint32())),
            ]
        )
        fields.append(pa.field("mz_grid", grid))
    chunk = pa.array(
        [{field.name: row[field.name] for field in fields} for row in rows],
        type=pa.struct(fields),
    )
    return pa.table({"chunk": chunk})


def _observed(
    profile: Optional[Spectrum], peaks: Optional[Spectrum]
) -> Optional[Spectrum]:
    """The arrays a spectrum's summary columns describe, if it has any."""
    for candidate in (profile, peaks):
        if candidate is not None and candidate.mzs.size:
            return candidate
    return None


def _metadata_table(
    profile: Sequence[Optional[Spectrum]],
    peaks: Sequence[Optional[Spectrum]],
    null_pair_after: Optional[int],
    spectrum_representation: str,
    count_columns: Sequence[str],
) -> pa.Table:
    """Build the one-row-per-spectrum metadata table.

    ``number_of_data_points`` counts stored rows, padding included -- which is
    what the reference archive does, and the reason the reader has to correct
    it before reporting a peak count.

    Each count is null where the spectrum has no rows in that member. The
    reference converter does the same: on a centroid input every
    ``number_of_data_points`` is null, and on a profile input every
    ``number_of_peaks`` is.
    """
    padding = 0 if null_pair_after is None else 2
    spectra = [_observed(a, b) for a, b in zip(profile, peaks)]
    counts = {
        "number_of_data_points": [
            None if s is None else int(s.mzs.size) + padding for s in profile
        ],
        "number_of_peaks": [None if s is None else int(s.mzs.size) for s in peaks],
    }
    table = _summary_table(spectra, spectrum_representation)
    for name in count_columns:
        table = table.append_column(name, pa.array(counts[name], type=pa.uint64()))
    return table


def _summary_table(
    spectra: Sequence[Optional[Spectrum]], spectrum_representation: str
) -> pa.Table:
    """Build the metadata columns that do not depend on the member."""
    return pa.table(
        {
            "index": pa.array(range(len(spectra)), type=pa.uint64()),
            "id": pa.array(
                [f"Scan={i + 1}" for i in range(len(spectra))],
                type=pa.large_string(),
            ),
            "ms_level": pa.array([1] * len(spectra), type=pa.uint8()),
            "time": pa.array(
                [float(i) for i in range(len(spectra))], type=pa.float64()
            ),
            "spectrum_representation": pa.array(
                [spectrum_representation] * len(spectra), type=pa.string()
            ),
            "lowest_observed_mz": pa.array(
                [None if s is None else float(s.mzs.min()) for s in spectra],
                type=pa.float64(),
            ),
            "highest_observed_mz": pa.array(
                [None if s is None else float(s.mzs.max()) for s in spectra],
                type=pa.float64(),
            ),
            "total_ion_current": pa.array(
                [None if s is None else float(s.intensities.sum()) for s in spectra],
                type=pa.float32(),
            ),
        }
    )


def _scans_table(spectra: Sequence[Spectrum], include_positions: bool) -> pa.Table:
    """Build the scans table, with or without the imaging position columns.

    ``source_index`` is the join key back to the spectrum; ``scan_index``
    numbers scans within a spectrum and is deliberately not the same column,
    because a reader that joins on row order passes until it meets a file
    where they differ.
    """
    columns: Dict[str, pa.Array] = {
        "source_index": pa.array(range(len(spectra)), type=pa.uint64()),
        "scan_index": pa.array([0] * len(spectra), type=pa.uint64()),
        "scan_start_time": pa.array(
            [float(i) for i in range(len(spectra))], type=pa.float32()
        ),
    }
    if include_positions:
        columns[POSITION_X_COLUMN] = pa.array([s.x for s in spectra], type=pa.uint32())
        columns[POSITION_Y_COLUMN] = pa.array([s.y for s in spectra], type=pa.uint32())
    return pa.table(columns)


def _scan_settings(
    pixel_size: Optional[Tuple[float, float]],
    grid: Optional[Tuple[int, int]],
    unit: Optional[str],
) -> List[dict]:
    """Build ``scan_settings_list``, imitating the reference spelling.

    The CV names are copied verbatim from the reference archive, parentheses
    and all: IMS:1000046 is written ``"pixel size (x)"`` and IMS:1000047
    ``"pixel size y"``. A reader matching on name rather than accession finds
    one axis and misses the other, and these fixtures are what makes that
    visible.
    """
    parameters: List[dict] = []
    if grid is not None:
        parameters += [
            {
                "name": "max count of pixels x",
                "accession": "IMS:1000042",
                "value": grid[0],
                "unit": None,
            },
            {
                "name": "max count of pixels y",
                "accession": "IMS:1000043",
                "value": grid[1],
                "unit": None,
            },
        ]
    if pixel_size is not None:
        parameters += [
            {
                "name": "pixel size (x)",
                "accession": "IMS:1000046",
                "value": pixel_size[0],
                "unit": unit,
            },
            {
                "name": "pixel size y",
                "accession": "IMS:1000047",
                "value": pixel_size[1],
                "unit": unit,
            },
        ]
    if not parameters:
        return []
    return [
        {
            "id": "scansettings1",
            "source_file_refs": [],
            "targets": [],
            "parameters": parameters,
        }
    ]


def _index_document(
    include_positions: bool,
    data_kind: str,
    column_mapping_key: Optional[str],
    index_metadata: Optional[Dict[str, Any]],
    signal_members: Sequence[str] = (DATA_MEMBER,),
    position_z_column: Optional[str] = None,
) -> Dict[str, Any]:
    """Spell out ``mzpeak_index.json``.

    Args:
        include_positions: Whether the scans entry binds the position terms.
        data_kind: How to spell the signal member's kind. The reference
            parser accepts ``"data_arrays"`` and ``"data arrays"`` alike, so
            both must resolve.
        column_mapping_key: Which key carries the CV bindings --
            ``"column_mapping"``, its serde alias ``"metadata_mapping"``, or
            ``None`` to omit the bindings entirely, which the reference
            struct permits via ``serde(default)``.
        index_metadata: Contents of the index's ``metadata`` object. The
            reference imaging archive leaves this empty and puts everything
            in the Parquet footer; other archives populate both.
        signal_members: Which signal members the archive holds. The
            reference converter lists both, data member first.
        position_z_column: The scans column holding position z, bound to
            IMS:1000052 like the other two, or ``None`` for no such column.
    """
    scans_bindings = [
        {
            "name": "scan start time",
            "path": "scan_start_time",
            "accession": "MS:1000016",
            "unit": "UO:0000031",
        }
    ]
    if include_positions:
        scans_bindings += [
            {
                "name": "position x",
                "path": POSITION_X_COLUMN,
                "accession": "IMS:1000050",
                "unit": None,
            },
            {
                "name": "position y",
                "path": POSITION_Y_COLUMN,
                "accession": "IMS:1000051",
                "unit": None,
            },
        ]
    if position_z_column is not None:
        scans_bindings.append(
            {
                "name": "position z",
                "path": position_z_column,
                "accession": "IMS:1000052",
                "unit": None,
            }
        )

    def entry(name: str, kind: str, bindings: List[dict]) -> Dict[str, Any]:
        item: Dict[str, Any] = {
            "name": name,
            "entity_type": "spectrum",
            "data_kind": kind,
        }
        if column_mapping_key is not None:
            item[column_mapping_key] = bindings
        item["parameters"] = []
        return item

    kinds = {DATA_MEMBER: data_kind, PEAKS_MEMBER: "peaks"}
    return {
        "files": [
            *(entry(name, kinds[name], []) for name in signal_members),
            entry(METADATA_MEMBER, "metadata", []),
            entry(SCANS_MEMBER, "scans", scans_bindings),
        ],
        "metadata": dict(index_metadata) if index_metadata else {},
    }


def _signal_plan(
    spectra: Sequence[Spectrum],
    signal: str,
    centroids: Optional[Sequence[Optional[Spectrum]]],
    empty_peer: bool,
) -> Dict[str, List[Optional[Spectrum]]]:
    """Decide which signal members to write and what each one holds.

    Returns:
        Member name to its per-spectrum content, data member first. A member
        that is absent from the mapping is not written at all; one whose
        entries are all ``None`` is written with zero rows.
    """
    nothing: List[Optional[Spectrum]] = [None] * len(spectra)
    filled: List[Optional[Spectrum]] = [s if s.mzs.size else None for s in spectra]
    if signal == "profile":
        plan = {DATA_MEMBER: filled}
        if empty_peer:
            plan[PEAKS_MEMBER] = nothing
    elif signal == "centroid":
        plan = {DATA_MEMBER: nothing} if empty_peer else {}
        plan[PEAKS_MEMBER] = filled
    elif signal == "both":
        if centroids is None or len(centroids) != len(spectra):
            raise ValueError("signal='both' needs one centroids entry per spectrum")
        plan = {DATA_MEMBER: filled, PEAKS_MEMBER: list(centroids)}
    else:  # pragma: no cover - programming error in a test
        raise ValueError(f"unknown signal {signal!r}")
    return plan


def build_mzpeak(
    path: Path,
    spectra: Sequence[Spectrum],
    *,
    pixel_size: Optional[Tuple[float, float]] = (25.0, 25.0),
    pixel_size_unit: Optional[str] = UNIT_MICROMETRE,
    grid: Optional[Tuple[int, int]] = None,
    include_positions: bool = True,
    row_group_size: Optional[int] = None,
    layout: str = "point",
    data_kind: str = "data_arrays",
    column_mapping_key: Optional[str] = "column_mapping",
    index_metadata: Optional[Dict[str, Any]] = None,
    null_pair_after: Optional[int] = None,
    spectrum_representation: Optional[str] = None,
    footer_metadata: bool = True,
    signal: str = "profile",
    centroids: Optional[Sequence[Optional[Spectrum]]] = None,
    empty_peer: bool = False,
    declare_counts: bool = True,
    file_contents: Optional[List[dict]] = None,
    chunk_encoding: Union[str, Sequence[str]] = DELTA,
    chunk_points: int = 4,
    grid_type: str = LINEAR_GRID,
    numpress_fixed_point: Optional[float] = None,
    file_metadata: Optional[Dict[str, Any]] = None,
    ms_levels: Optional[Sequence[int]] = None,
    ion_mobility: Optional[Sequence[Optional[float]]] = None,
    position_z: Optional[Sequence[Optional[int]]] = None,
    position_z_column: str = POSITION_Z_COLUMN,
) -> Path:
    """Write one ``.mzpeak`` archive and return its path.

    Args:
        path: Destination file. Parent directories must exist.
        spectra: Spectra, in the order they should be indexed. They carry
            the positions, and the signal of the member ``signal`` names.
        pixel_size: ``(x, y)`` pixel size, or ``None`` to omit the terms and
            exercise the "pixel size not found" path.
        pixel_size_unit: Unit accession for the pixel size terms.
        grid: Declared grid extent (IMS:1000042/43), or ``None`` to omit.
        include_positions: When ``False``, the scans member carries no
            position columns, which is what a non-imaging archive looks like.
        row_group_size: Rows per Parquet row group. Set this below a
            spectrum's point count to force spectra across row-group
            boundaries.
        layout: ``"point"`` or ``"chunk"``, for the member ``signal``
            fills. An empty peer is always written in the point layout.
        data_kind: Spelling for the signal member's ``data_kind``.
        column_mapping_key: Key carrying the CV bindings, or ``None``.
        index_metadata: Contents of the index ``metadata`` object.
        null_pair_after: Insert a null pair after this many points of each
            spectrum. Applies to the data member only: padding marks zero
            runs removed from a profile spectrum.
        spectrum_representation: Value for the per-spectrum column. Defaults
            to the CV name matching ``signal``.
        footer_metadata: Whether to write the file-level JSON blobs into the
            metadata member's Parquet key-value footer. ``False`` leaves the
            index ``metadata`` object as the only source, which is how a
            reader that reads just one of the two gets caught.
        signal: Where ``spectra`` are written. ``"profile"`` puts them in
            the data member, ``"centroid"`` in the peaks member, and
            ``"both"`` puts them in the data member and ``centroids`` in the
            peaks member. A spectrum with no points writes no rows and
            records a null count.
        centroids: Peak lists for ``signal="both"``, one per spectrum, with
            ``None`` for a spectrum that has no peak list.
        empty_peer: Also write the member ``signal`` leaves unused, with
            zero rows and an all-null count column. The reference converter
            always does; a fixture without it has one signal member only.
        declare_counts: Whether the metadata member carries the count
            columns at all.
        file_contents: Replaces ``file_description.contents``. The reference
            converter's imaging archives do not always declare the
            representation there, so this is how to leave it out.
        chunk_encoding: For ``layout="chunk"``, the accession of the
            encoding the m/z are stored under, or one per row to mix them.
        chunk_points: For ``layout="chunk"``, stored points per row.
        grid_type: For the grid encoding, the accession of the model.
        numpress_fixed_point: For MS-Numpress, the scale. The default is a
            power of two, which gives the fixture's m/z back to the bit. A
            scale such as 1e5 gives them back a little off, as a real
            archive does.
        file_metadata: File-level blocks written over the defaults, such
            as ``imaging``, ``bruker_maldi``, ``run`` or
            ``instrument_configuration_list``, spelled as mzpeak-convert
            writes them.
        ms_levels: MS level of each spectrum. All 1 by default.
        ion_mobility: Ion mobility value of each scan, ``None`` for none.
            Without it the scans member has no such column.
        position_z: Position z of each scan, ``None`` for a null. Without
            it the scans member has no z column, as mzpeak-convert writes
            for a source that states no z.
        position_z_column: Name of that column. It is bound to IMS:1000052
            in the index whenever the index carries bindings.

    Returns:
        ``path``, for convenience.
    """
    if layout not in ("point", "chunk"):  # pragma: no cover - test error
        raise ValueError(f"unknown layout {layout!r}")

    plan = _signal_plan(spectra, signal, centroids, empty_peer)
    filled_member = PEAKS_MEMBER if signal == "centroid" else DATA_MEMBER
    tables: Dict[str, pa.Table] = {}
    for member, content in plan.items():
        padding = null_pair_after if member == DATA_MEMBER else None
        if layout == "chunk" and member == filled_member:
            tables[member] = chunk_table(
                content,
                padding,
                chunk_encoding,
                chunk_points,
                grid_type,
                numpress_fixed_point,
            )
        else:
            tables[member] = _point_table(content, padding)

    if spectrum_representation is None:
        spectrum_representation = (
            "centroid spectrum" if signal == "centroid" else "profile spectrum"
        )
    if file_contents is None:
        file_contents = [
            {
                "name": spectrum_representation,
                "accession": (
                    "MS:1000127"
                    if spectrum_representation == "centroid spectrum"
                    else "MS:1000128"
                ),
                "value": None,
                "unit": None,
            }
        ]

    scan_settings = _scan_settings(pixel_size, grid, pixel_size_unit)
    footer = {
        "file_description": {
            "contents": file_contents,
            "source_files": [],
        },
        "scan_settings_list": scan_settings,
        "instrument_configuration_list": [],
        "software_list": [],
        "run": {"id": "fixture", "start_time": None},
        **(file_metadata or {}),
    }

    count_columns = {
        DATA_MEMBER: "number_of_data_points",
        PEAKS_MEMBER: "number_of_peaks",
    }
    nothing: List[Optional[Spectrum]] = [None] * len(spectra)
    metadata_table = _metadata_table(
        plan.get(DATA_MEMBER, nothing),
        plan.get(PEAKS_MEMBER, nothing),
        null_pair_after,
        spectrum_representation,
        [count_columns[member] for member in plan] if declare_counts else [],
    )
    if ms_levels is not None:
        metadata_table = metadata_table.set_column(
            metadata_table.schema.get_field_index("ms_level"),
            "ms_level",
            pa.array(list(ms_levels), type=pa.uint8()),
        )
    if footer_metadata:
        metadata_table = metadata_table.replace_schema_metadata(
            {key: json.dumps(value) for key, value in footer.items()}
        )
    scans_table = _scans_table(spectra, include_positions)
    if ion_mobility is not None:
        scans_table = scans_table.append_column(
            "ion_mobility_value", pa.array(list(ion_mobility), type=pa.float64())
        )
    if position_z is not None:
        scans_table = scans_table.append_column(
            position_z_column, pa.array(list(position_z), type=pa.uint32())
        )

    path = Path(path)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        for member, table, group_size in (
            *((name, tables[name], row_group_size) for name in tables),
            (METADATA_MEMBER, metadata_table, None),
            (SCANS_MEMBER, scans_table, None),
        ):
            sink = pa.BufferOutputStream()
            kwargs: Dict[str, Any] = {"write_statistics": True}
            if group_size is not None:
                kwargs["row_group_size"] = group_size
            pq.write_table(table, sink, **kwargs)
            archive.writestr(member, sink.getvalue().to_pybytes())

        archive.writestr(
            INDEX_MEMBER,
            json.dumps(
                _index_document(
                    include_positions,
                    data_kind,
                    column_mapping_key,
                    index_metadata,
                    list(tables),
                    position_z_column if position_z is not None else None,
                ),
                indent=2,
            ),
        )
    return path


def grid_spectra(
    n_x: int,
    n_y: int,
    n_points: int = 8,
    *,
    base: int = 1,
    skip: Sequence[Tuple[int, int]] = (),
) -> List[Spectrum]:
    """Generate a rectangular acquisition with deterministic arrays.

    Args:
        n_x: Pixels along x.
        n_y: Pixels along y.
        n_points: Points per spectrum.
        base: Coordinate origin as written into the file. The reference
            archives are 1-based; ``0`` exercises the normalisation.
        skip: Positions to leave unacquired, so the grid is sparse.

    Returns:
        Spectra in raster order.
    """
    skipped = {(int(x), int(y)) for x, y in skip}
    spectra: List[Spectrum] = []
    for y in range(n_y):
        for x in range(n_x):
            if (x + base, y + base) in skipped:
                continue
            offset = y * n_x + x
            mzs = 100.0 + np.arange(n_points, dtype=np.float64) * 0.5
            intensities = np.arange(1, n_points + 1, dtype=np.float64) + offset * 10.0
            spectra.append(Spectrum(x + base, y + base, mzs, intensities))
    return spectra
