# thyra/readers/mzpeak/chunk_decoding.py
"""Decode the chunked layout of an mzPeak signal member.

In the chunked layout one row holds a run of points of one spectrum, cut
along m/z. The intensities of the run are a plain list. Its m/z values are
stored under the encoding the row names in ``chunk_encoding``
(specification, ``docs/layouts/chunked-layout.md``):

* ``MS:1000576``, no compression: the values as they are, first one left
  out because ``mz_chunk_start`` holds it.
* ``MS:1003089``, delta encoding: the difference of each value to the one
  before it, first value in ``mz_chunk_start``.
* ``MS:1002312``, MS-Numpress linear prediction: a byte buffer in
  ``mz_numpress_linear_bytes``. Lossy by design.
* ``MS:1003826``, coordinate grid encoding: integer grid indices and the
  parameters of the model that turns an index into m/z, in ``mz_grid``.

A decoded row always has as many m/z values as it has intensities. The
intensity list is a plain list in every encoding, so its length is the one
count that needs no decoding, and every decoder is held to it. A row that
disagrees is refused by name.

Two of the encodings are lossy: a decoded m/z is the stored one to within
1e-6 or so. The two ends of a row are an exception. ``mz_chunk_start`` and
``mz_chunk_end`` hold the first and the last m/z of the row as plain
numbers, so the decoded ends are replaced by them. That keeps every
spectrum inside the m/z range the archive declares for it, which is the
range the resampled axis is built on: a first peak that decodes 1e-7 below
that range would be left out of the store.

Null marking follows the point layout: a point whose m/z or intensity is
null is padding for a removed run of zeros. It comes back as ``NaN`` in the
m/z array, and the reader drops it. MS-Numpress has no null, so the
reference writer stores such an m/z as ``0.0`` and the reference reader
turns ``0.0`` back into a null; this module does the same.
"""

from __future__ import annotations

from typing import Any, Dict, NamedTuple, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from ...errors import ConversionRefused

#: Chunk encodings this module decodes, with their CV names.
ENCODING_PLAIN = "MS:1000576"
ENCODING_DELTA = "MS:1003089"
ENCODING_NUMPRESS_LINEAR = "MS:1002312"
ENCODING_GRID = "MS:1003826"

DECODED_ENCODINGS: Dict[str, str] = {
    ENCODING_PLAIN: "no compression",
    ENCODING_DELTA: "truncation, delta prediction and zlib compression",
    ENCODING_NUMPRESS_LINEAR: "MS-Numpress linear prediction compression",
    ENCODING_GRID: "coordinate grid encoding",
}

#: Encodings that give an m/z back a little off.
LOSSY_ENCODINGS = (ENCODING_NUMPRESS_LINEAR, ENCODING_GRID)

#: Grid models this module decodes. Both are open: the term's definition
#: gives the equation.
GRID_LINEAR = "MS:1003824"
GRID_SQUARE_ROOT = "MS:1003825"

DECODED_GRID_TYPES: Dict[str, str] = {
    GRID_LINEAR: "linear grid interpolation",
    GRID_SQUARE_ROOT: "square root grid interpolation",
}

#: Terms a refusal can name although nothing here decodes them.
REFUSED_TERMS: Dict[str, str] = {
    "MS:1002313": "MS-Numpress positive integer compression",
    "MS:1002314": "MS-Numpress short logged float compression",
    "MS:1003822": "grid coordinate interpolation",
    "MS:9999001": "timsTOF ion mobility grid of the reference implementation",
    "MS:9999002": "timsTOF m/z grid of the reference implementation",
}

#: Transforms of the intensity array, by the column that holds them.
INTENSITY_TRANSFORMS: Dict[str, str] = {
    "intensity_numpress_slof_bytes": "MS:1002314",
    "intensity_numpress_pic_bytes": "MS:1002313",
}

#: Columns a chunk row must have for this module to read it.
REQUIRED_COLUMNS = (
    "spectrum_index",
    "mz_chunk_start",
    "mz_chunk_end",
    "chunk_encoding",
)

#: The column each encoding keeps its values in.
ENCODING_COLUMNS: Dict[str, str] = {
    ENCODING_PLAIN: "mz_chunk_values",
    ENCODING_DELTA: "mz_chunk_values",
    ENCODING_NUMPRESS_LINEAR: "mz_numpress_linear_bytes",
    ENCODING_GRID: "mz_grid",
}

#: How far a decoded end of a row may lie from the bound that replaces it.
#: Ten times the largest error measured on a lossy encoding. A bound
#: further off than this is not the value of that point, and is left alone.
BOUND_TOLERANCE = 1e-5

#: Runs up to this length are summed side by side, longer ones one by one.
SHORT_RUN = 256


class DecodedChunks(NamedTuple):
    """The points of a batch of chunk rows, in stored order.

    Attributes:
        spectrum_index: ``spectrum_index`` of each point.
        mz: m/z of each point. ``NaN`` marks padding.
        intensity: Intensity of each point.
    """

    spectrum_index: NDArray[np.int64]
    mz: NDArray[np.float64]
    intensity: NDArray[np.float64]


def describe_term(accession: Optional[str]) -> str:
    """Give a CV term as ``accession (name)``, or alone when unnamed here."""
    if not accession:
        return "(none given)"
    for names in (DECODED_ENCODINGS, DECODED_GRID_TYPES, REFUSED_TERMS):
        if accession in names:
            return f"{accession} ({names[accession]})"
    return str(accession)


def _decoded(terms: Dict[str, str]) -> str:
    """List the terms that are decoded, for a refusal."""
    return ", ".join(describe_term(term) for term in terms)


def validate_chunk_columns(names: Sequence[str], source: str) -> None:
    """Refuse a chunk row this module cannot read, before any row is read.

    Args:
        names: Children of the ``chunk`` struct.
        source: Name of the archive, for the message.

    Raises:
        ConversionRefused: If the rows are not cut along m/z, or the
            intensity array is stored under a transform.
    """
    missing = [name for name in REQUIRED_COLUMNS if name not in names]
    if missing:
        raise ConversionRefused(
            f"{source} uses the mzPeak chunked layout, but its chunks lack "
            f"{', '.join(missing)}. Thyra reads chunks cut along m/z; this "
            f"member holds {', '.join(names)}."
        )
    if "intensity" in names:
        return
    for column, accession in INTENSITY_TRANSFORMS.items():
        if column in names:
            raise ConversionRefused(
                f"{source} stores its intensity array under "
                f"{describe_term(accession)}, in column '{column}'. Thyra "
                f"does not decode that transform. Convert the source again "
                f"with the intensities stored as they are."
            )
    raise ConversionRefused(
        f"{source} uses the mzPeak chunked layout, but its chunks have no "
        f"'intensity' column. This member holds {', '.join(names)}."
    )


def validate_encodings(
    encodings: Sequence[Optional[str]],
    grid_types: Sequence[Optional[str]],
    names: Sequence[str],
    source: str,
) -> None:
    """Refuse an encoding or a grid model this module does not decode.

    Args:
        encodings: Distinct values of ``chunk_encoding``.
        grid_types: Distinct values of ``mz_grid.grid_type``.
        names: Children of the ``chunk`` struct.
        source: Name of the archive, for the message.

    Raises:
        ConversionRefused: Naming the term that is not decoded.
    """
    for encoding in encodings:
        if encoding not in DECODED_ENCODINGS:
            raise ConversionRefused(
                f"{source} stores m/z under the chunk encoding "
                f"{describe_term(encoding)}, which Thyra does not decode. "
                f"It decodes {_decoded(DECODED_ENCODINGS)}."
            )
        column = ENCODING_COLUMNS[encoding]
        if column not in names:
            raise ConversionRefused(
                f"{source} names the chunk encoding "
                f"{describe_term(encoding)}, but its chunks have no "
                f"'{column}' column to decode."
            )
    if ENCODING_GRID not in encodings:
        return
    for grid_type in grid_types:
        if grid_type not in DECODED_GRID_TYPES:
            raise ConversionRefused(
                f"{source} stores m/z on a grid of type "
                f"{describe_term(grid_type)}, which Thyra does not decode. "
                f"It decodes {_decoded(DECODED_GRID_TYPES)}."
            )


def _list_parts(array: Any) -> Tuple[NDArray[Any], NDArray[np.int64]]:
    """Split a list array into its values and the length of each list.

    A null list counts as empty. Null values come back as ``NaN``.
    """
    import pyarrow.compute as pc  # noqa: WPS433 - deliberate lazy import

    lengths = np.asarray(
        pc.list_value_length(array).fill_null(0).to_numpy(zero_copy_only=False),
        dtype=np.int64,
    )
    values = array.flatten().to_numpy(zero_copy_only=False)
    return values, lengths


def _starts(lengths: NDArray[np.int64]) -> NDArray[np.int64]:
    """Offset of each run in the array that holds them end to end."""
    return np.asarray(np.cumsum(lengths) - lengths, dtype=np.int64)


def _positions(
    starts: NDArray[np.int64], lengths: NDArray[np.int64]
) -> NDArray[np.int64]:
    """Every index of every run, run after run."""
    total = int(lengths.sum())
    shift = np.repeat(starts - _starts(lengths), lengths)
    return np.asarray(shift + np.arange(total, dtype=np.int64), dtype=np.int64)


def _row_sums(
    values: NDArray[np.int64], counts: NDArray[np.int64]
) -> NDArray[np.int64]:
    """Running sum of whole numbers, restarted at every row.

    Whole numbers add exactly, so one sum over all rows, less what had
    gathered before each row began, is each row's own sum.
    """
    total = np.cumsum(values)
    first = _starts(counts)
    holds = counts > 0
    before = np.zeros(counts.size, dtype=np.int64)
    before[holds] = total[first[holds]] - values[first[holds]]
    return np.asarray(total - np.repeat(before, counts), dtype=np.int64)


def _refuse_count(source: str, encoding: str, row: int, found: int, n: int) -> None:
    """Raise for a row whose m/z and intensity counts disagree."""
    raise ConversionRefused(
        f"{source} holds a chunk whose m/z and intensity counts disagree: "
        f"row {row} of a batch, encoded as {describe_term(encoding)}, "
        f"decodes to {found} m/z values for {n} intensities."
    )


def running_sums(
    values: NDArray[np.float64],
    starts: NDArray[np.int64],
    lengths: NDArray[np.int64],
) -> None:
    """Replace each run of ``values`` by its running sum, in place.

    Every sum is taken in stored order, one addition per value, so a run
    decodes to the same bits whichever route it takes here. Short runs are
    summed side by side, one position at a time, and long runs one run at a
    time; the number of Python steps stays small either way.

    Args:
        values: The array holding every run.
        starts: First index of each run.
        lengths: Length of each run.
    """
    long_runs = lengths > SHORT_RUN
    for start, length in zip(starts[long_runs], lengths[long_runs]):
        values[start : start + length] = np.cumsum(values[start : start + length])

    short = ~long_runs & (lengths > 1)
    if not short.any():
        return
    order = np.argsort(-lengths[short], kind="stable")
    run_lengths = lengths[short][order]
    run_starts = starts[short][order]
    # Runs are sorted longest first, so the runs that reach a position are
    # always the first few.
    reach = np.searchsorted(-run_lengths, -np.arange(int(run_lengths[0])), side="left")
    total = values[run_starts].copy()
    for position in range(1, int(run_lengths[0])):
        count = int(reach[position])
        index = run_starts[:count] + position
        total[:count] += values[index]
        values[index] = total[:count]


def decode_plain(
    starts: NDArray[np.float64],
    values: NDArray[np.float64],
    lengths: NDArray[np.int64],
    counts: NDArray[np.int64],
    source: str,
) -> NDArray[np.float64]:
    """Decode rows stored without compression.

    Args:
        starts: ``mz_chunk_start`` of each row.
        values: Every row's ``mz_chunk_values``, end to end.
        lengths: Length of each row's ``mz_chunk_values``.
        counts: Points each row must decode to.
        source: Name of the archive, for the message.

    Returns:
        The m/z of every point, row after row.
    """
    wrong = np.flatnonzero(lengths + 1 != counts)
    if wrong.size:
        row = int(wrong[0])
        _refuse_count(
            source, ENCODING_PLAIN, row, int(lengths[row]) + 1, int(counts[row])
        )
    decoded = np.empty(int(counts.sum()), dtype=np.float64)
    first = _starts(counts)
    decoded[first] = starts
    decoded[_positions(first + 1, lengths)] = values
    return decoded


def decode_delta(
    starts: NDArray[np.float64],
    values: NDArray[np.float64],
    lengths: NDArray[np.int64],
    counts: NDArray[np.int64],
    source: str,
) -> NDArray[np.float64]:
    """Decode delta-encoded rows, nulls included.

    A value after a null is stored as it is, and differences resume after
    it. A row normally leaves its first value to ``mz_chunk_start``. A row
    that began with a null keeps that null in the list instead, and then
    the list is as long as the intensities.

    The reference decoder tells the two apart by looking at the second
    value. The intensity count says it outright, so it is used here.

    Args:
        starts: ``mz_chunk_start`` of each row.
        values: Every row's ``mz_chunk_values``, end to end, ``NaN`` for
            null.
        lengths: Length of each row's ``mz_chunk_values``.
        counts: Points each row must decode to.
        source: Name of the archive, for the message.

    Returns:
        The m/z of every point, row after row, ``NaN`` where it was null.
    """
    first = _starts(counts)
    leads = lengths + 1 == counts
    # A list as long as the intensities must open with the null that took
    # the place of the start.
    opens_null = np.zeros(lengths.size, dtype=bool)
    filled = (lengths == counts) & (lengths > 0)
    opens_null[filled] = np.isnan(values[_starts(lengths)[filled]])
    empty = (lengths == 0) & (counts == 0)
    wrong = np.flatnonzero(~leads & ~opens_null & ~empty)
    if wrong.size:
        row = int(wrong[0])
        _refuse_count(
            source, ENCODING_DELTA, row, int(lengths[row]) + 1, int(counts[row])
        )

    decoded = np.empty(int(counts.sum()), dtype=np.float64)
    decoded[first[leads]] = starts[leads]
    decoded[_positions(first + leads, lengths)] = values

    # A run is a stretch of values with no null in it, inside one row.
    valid = ~np.isnan(decoded)
    opens = np.zeros(decoded.size, dtype=bool)
    opens[first[counts > 0]] = True
    before = np.concatenate(([False], valid[:-1]))
    after = np.concatenate((valid[1:], [False]))
    closes = np.concatenate((opens[1:], [True]))
    run_starts = np.flatnonzero(valid & (opens | ~before))
    run_ends = np.flatnonzero(valid & (closes | ~after))
    running_sums(decoded, run_starts, run_ends - run_starts + 1)
    return decoded


def _half_bytes(data: NDArray[np.uint8]) -> NDArray[np.uint8]:
    """The half bytes of ``data``, high half of each byte first."""
    halves = np.empty(2 * data.size, dtype=np.uint8)
    halves[0::2] = data >> 4
    halves[1::2] = data & 0xF
    return halves


def _leading(head: NDArray[Any]) -> NDArray[Any]:
    """Half bytes a residual leaves out, from the half byte that heads it."""
    return np.where(head <= 8, head, head - 8)


def _residual_heads(
    halves: NDArray[np.uint8],
    seeds: NDArray[np.int64],
    ends: NDArray[np.int64],
    wanted: NDArray[np.int64],
) -> Tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]:
    """Find where each residual of each MS-Numpress buffer begins.

    A residual is a half byte that says how many half bytes follow it, so
    the place of one is known only from the one before. Walking them one at
    a time would cost one step per value. Here every half byte is given the
    place the next residual would have if it were one itself, and the
    places reached from the seeds are found by doubling the stride.

    Args:
        halves: The half bytes of every buffer, end to end.
        seeds: Half-byte position of the first residual of each buffer.
        ends: Half-byte position just past each buffer.
        wanted: Residuals each buffer must hold.

    Returns:
        ``(heads, owner, rank)``: the position of every half byte reached,
        ascending, the buffer it lies in, and its number within that
        buffer.
    """
    size = halves.size
    kind = np.int32 if size < np.iinfo(np.int32).max else np.int64
    jump = np.full(size + 1, size, dtype=kind)
    jump[:size] = np.arange(9, size + 9, dtype=kind) - _leading(halves)
    # A residual may end on the last half byte of its buffer. One that
    # would end past it leads nowhere; ``size`` stands for nowhere.
    limit = np.zeros(size + 1, dtype=kind)
    lengths = ends - seeds
    limit[_positions(seeds, lengths)] = np.repeat(ends, lengths)
    jump[jump >= limit] = size
    del limit

    reached = seeds[(wanted > 0) & (seeds < ends)]
    stride = 1
    longest = int(wanted.max()) + 1 if wanted.size else 0
    while stride < longest:
        further = jump[reached]
        reached = np.concatenate((reached, further[further < size]))
        jump = jump[jump]
        stride *= 2
    heads = np.sort(reached).astype(np.int64)
    owner = np.searchsorted(seeds, heads, side="right") - 1
    found = np.bincount(owner, minlength=seeds.size)
    rank = np.arange(heads.size, dtype=np.int64) - _starts(found)[owner]
    return heads, owner, rank


def _residuals(
    halves: NDArray[np.uint8], heads: NDArray[np.int64]
) -> NDArray[np.int64]:
    """Read the signed 32-bit integer that begins at each head."""
    head = halves[heads].astype(np.int64)
    leading = _leading(head)
    value = np.zeros(heads.size, dtype=np.uint64)
    last = halves.size - 1
    for place in range(8):
        stored = place < 8 - leading
        digit = halves[np.minimum(heads + 1 + place, last)]
        value |= np.where(stored, digit, 0).astype(np.uint64) << np.uint64(4 * place)
    # The half bytes left out are all zero, or all one for a negative value.
    shift = (4 * (8 - leading)).astype(np.uint64)
    fill = (np.uint64(0xFFFFFFFF) << shift) & np.uint64(0xFFFFFFFF)
    value |= np.where(head > 8, fill, 0).astype(np.uint64)
    return value.astype(np.uint32).view(np.int32).astype(np.int64)


def _little_endian_uint32(
    data: NDArray[np.uint8], at: NDArray[np.int64]
) -> NDArray[np.int64]:
    """Read an unsigned 32-bit integer at each byte offset."""
    value = np.zeros(at.size, dtype=np.int64)
    for place in range(4):
        value |= data[at + place].astype(np.int64) << (8 * place)
    return value


def _numpress_residuals(
    data: NDArray[np.uint8],
    begin: NDArray[np.int64],
    lengths: NDArray[np.int64],
    counts: NDArray[np.int64],
    source: str,
) -> NDArray[np.int64]:
    """Read every residual of every buffer, held to the expected count."""
    wanted = np.maximum(counts - 2, 0)
    ends = 2 * (begin + lengths)
    # A buffer with no residual has no seed inside it.
    seeds = np.where(wanted > 0, 2 * (begin + 16), ends)
    halves = _half_bytes(data)
    heads, owner, rank = _residual_heads(halves, seeds, ends, wanted)

    found = np.bincount(owner, minlength=counts.size)
    # One half byte may follow the last residual: the padding that fills
    # the last byte. It is zero and it is the last half byte.
    surplus = rank >= wanted[owner]
    padding = surplus & (heads == ends[owner] - 1) & (halves[heads] == 0)
    found -= np.bincount(owner[padding], minlength=counts.size)
    wrong = np.flatnonzero(found != wanted)
    if wrong.size:
        row = int(wrong[0])
        _refuse_count(
            source,
            ENCODING_NUMPRESS_LINEAR,
            row,
            int(found[row]) + 2,
            int(counts[row]),
        )

    heads = heads[~surplus]
    owner = owner[~surplus]
    past = heads + 9 - _leading(halves[heads].astype(np.int64))
    cut = np.flatnonzero(past > ends[owner])
    if cut.size:
        raise ConversionRefused(
            f"{source} holds an MS-Numpress buffer that ends inside a "
            f"value (row {int(owner[cut[0]])} of a batch)."
        )
    return _residuals(halves, heads)


def decode_numpress_linear(
    data: NDArray[np.uint8],
    lengths: NDArray[np.int64],
    counts: NDArray[np.int64],
    source: str,
) -> NDArray[np.float64]:
    """Decode MS-Numpress linear prediction buffers.

    A buffer is a scale as an 8-byte float, the first two values as 4-byte
    integers, and one residual per further value. A value is predicted from
    the two before it, ``2 a - b``, and the residual is what the prediction
    missed (Teleman et al., Mol Cell Proteomics 2014, 13, 1537).

    Args:
        data: The bytes of every buffer, end to end.
        lengths: Bytes in each buffer.
        counts: Points each buffer must decode to.
        source: Name of the archive, for the message.

    Returns:
        The m/z of every point, buffer after buffer, ``NaN`` where the
        stored value was 0.
    """
    header = np.where(counts == 0, 0, np.where(counts == 1, 12, 16))
    short = np.flatnonzero(lengths < header)
    if short.size:
        row = int(short[0])
        raise ConversionRefused(
            f"{source} holds an MS-Numpress buffer of {int(lengths[row])} "
            f"bytes for {int(counts[row])} values (row {row} of a batch); "
            f"that is too short to hold them."
        )

    begin = _starts(lengths)
    first = _starts(counts)
    holds = counts > 0
    pair = counts > 1
    scale = np.ones(counts.size, dtype=np.float64)
    scale_bytes = data[begin[holds][:, None] + np.arange(8)]
    scale[holds] = np.ascontiguousarray(scale_bytes).view(">f8")[:, 0]

    # value[k] = 2 value[k-1] - value[k-2] + residual[k], so the step from
    # one value to the next grows by the residual. These are whole numbers,
    # so two running sums give every value exactly.
    wanted = np.maximum(counts - 2, 0)
    value0 = _little_endian_uint32(data, begin[holds] + 8)
    value1 = _little_endian_uint32(data, begin[pair] + 12)
    growth = np.zeros(int(counts.sum()), dtype=np.int64)
    growth[first[pair] + 1] = value1 - value0[pair[holds]]
    growth[_positions(first + 2, wanted)] = _numpress_residuals(
        data, begin, lengths, counts, source
    )
    step = _row_sums(growth, counts)
    step[first[holds]] = value0
    value = _row_sums(step, counts)

    decoded = value / np.repeat(scale, counts)
    decoded[decoded == 0.0] = np.nan
    return np.asarray(decoded, dtype=np.float64)


def decode_grid(
    grid_type: str,
    parameters: Tuple[NDArray[np.float64], NDArray[np.int64]],
    indices: Tuple[NDArray[Any], NDArray[np.int64]],
    counts: NDArray[np.int64],
    source: str,
) -> NDArray[np.float64]:
    """Decode grid-encoded rows of one grid model.

    A row holds grid indices, each stored as the difference to the one
    before it, and the intercept, slope and scale of its model. The linear
    model gives ``(intercept + slope * index) / scale``; the square root
    model gives the square of the bracket, over the scale.

    Args:
        grid_type: ``MS:1003824`` or ``MS:1003825``.
        parameters: Every row's parameters end to end, and how many each
            row has.
        indices: Every row's indices end to end, and how many each row
            has.
        counts: Points each row must decode to.
        source: Name of the archive, for the message.

    Returns:
        The m/z of every point, row after row.
    """
    values, n_parameters = parameters
    steps, n_indices = indices
    wrong = np.flatnonzero(n_indices != counts)
    if wrong.size:
        row = int(wrong[0])
        _refuse_count(source, ENCODING_GRID, row, int(n_indices[row]), int(counts[row]))
    few = np.flatnonzero(n_parameters < 2)
    if few.size:
        raise ConversionRefused(
            f"{source} holds a grid of type {describe_term(grid_type)} with "
            f"{int(n_parameters[few[0]])} parameters; the model needs an "
            f"intercept and a slope."
        )

    at = _starts(n_parameters)
    intercept = values[at]
    slope = values[at + 1]
    scale = np.ones(counts.size, dtype=np.float64)
    scaled = n_parameters > 2
    scale[scaled] = values[at[scaled] + 2]

    index = _row_sums(steps.astype(np.int64), counts).astype(np.float64)
    bracket = np.repeat(intercept, counts) + index * np.repeat(slope, counts)
    if grid_type == GRID_SQUARE_ROOT:
        bracket = bracket * bracket
    return np.asarray(bracket / np.repeat(scale, counts), dtype=np.float64)


def _decode_grid_rows(
    column: Any, counts: NDArray[np.int64], source: str
) -> NDArray[np.float64]:
    """Decode grid rows, one grid model at a time."""
    import pyarrow as pa  # noqa: WPS433 - deliberate lazy import

    types = np.array([str(value) for value in column.field("grid_type").to_pylist()])
    validate_encodings([ENCODING_GRID], sorted(set(types)), ["mz_grid"], source)
    decoded = np.empty(int(counts.sum()), dtype=np.float64)
    first = _starts(counts)
    for grid_type in sorted(set(types)):
        rows = np.flatnonzero(types == grid_type)
        part = column.take(pa.array(rows))
        parameters = _list_parts(part.field("parameters"))
        decoded[_positions(first[rows], counts[rows])] = decode_grid(
            grid_type,
            (np.asarray(parameters[0], dtype=np.float64), parameters[1]),
            _list_parts(part.field("indices")),
            counts[rows],
            source,
        )
    return decoded


def take_bounds(
    mz: NDArray[np.float64],
    counts: NDArray[np.int64],
    starts: NDArray[np.float64],
    ends: NDArray[np.float64],
) -> None:
    """Replace the decoded ends of each row by the row's bounds, in place.

    The first and the last m/z of a row that are not padding are set to
    ``mz_chunk_start`` and ``mz_chunk_end``, where the decoded value is
    within :data:`BOUND_TOLERANCE` of the bound.

    Args:
        mz: The decoded m/z of every row, ``NaN`` on padding.
        counts: Points in each row.
        starts: ``mz_chunk_start`` of each row.
        ends: ``mz_chunk_end`` of each row.
    """
    rows = np.flatnonzero(counts > 0)
    if rows.size == 0:
        return
    begin = _starts(counts)[rows]
    place = np.arange(mz.size, dtype=np.int64)
    valid = ~np.isnan(mz)
    lowest = np.minimum.reduceat(np.where(valid, place, mz.size), begin)
    highest = np.maximum.reduceat(np.where(valid, place, -1), begin)
    for at, bound in ((lowest, starts[rows]), (highest, ends[rows])):
        held = np.flatnonzero(lowest <= highest)
        near = np.abs(mz[at[held]] - bound[held]) <= BOUND_TOLERANCE
        mz[at[held][near]] = bound[held][near]


def _decode_rows(
    chunk: Any,
    encoding: str,
    rows: NDArray[np.int64],
    starts: NDArray[np.float64],
    counts: NDArray[np.int64],
    source: str,
) -> NDArray[np.float64]:
    """Decode the rows of one encoding, in row order."""
    import pyarrow as pa  # noqa: WPS433 - deliberate lazy import

    column = chunk.field(ENCODING_COLUMNS[encoding]).take(pa.array(rows))
    if encoding == ENCODING_GRID:
        return _decode_grid_rows(column, counts, source)
    values, lengths = _list_parts(column)
    if encoding == ENCODING_NUMPRESS_LINEAR:
        return decode_numpress_linear(
            np.asarray(values, dtype=np.uint8), lengths, counts, source
        )
    values = np.asarray(values, dtype=np.float64)
    if encoding == ENCODING_DELTA:
        return decode_delta(starts, values, lengths, counts, source)
    return decode_plain(starts, values, lengths, counts, source)


def _numbers(chunk: Any, name: str, kind: Any) -> NDArray[Any]:
    """One child of the ``chunk`` struct, as numpy."""
    return np.asarray(chunk.field(name).to_numpy(zero_copy_only=False), dtype=kind)


def decode_chunks(chunk: Any, source: str) -> DecodedChunks:
    """Turn a batch of chunk rows into points.

    Args:
        chunk: The ``chunk`` struct column of the batch.
        source: Name of the archive, for messages.

    Returns:
        The points of the batch, row after row.

    Raises:
        ConversionRefused: If a row names an encoding that is not decoded,
            or decodes to a different count than its intensities.
    """
    names = [chunk.type.field(i).name for i in range(chunk.type.num_fields)]
    validate_chunk_columns(names, source)

    intensity, counts = _list_parts(chunk.field("intensity"))
    intensity = np.asarray(intensity, dtype=np.float64)
    starts = _numbers(chunk, "mz_chunk_start", np.float64)
    ends = _numbers(chunk, "mz_chunk_end", np.float64)
    encodings = np.array(
        [str(value or "") for value in chunk.field("chunk_encoding").to_pylist()]
    )

    mz = np.full(int(counts.sum()), np.nan, dtype=np.float64)
    first = _starts(counts)
    # A chunk of nulls alone is written with both bounds at zero. It holds
    # no m/z, so its points stay padding.
    holds = ~((starts == 0.0) & (ends == 0.0))
    used = sorted(set(encodings[holds]))
    validate_encodings(used, [], names, source)
    for encoding in used:
        rows = np.flatnonzero(holds & (encodings == encoding))
        mz[_positions(first[rows], counts[rows])] = _decode_rows(
            chunk, encoding, rows, starts[rows], counts[rows], source
        )

    take_bounds(mz, counts, starts, ends)
    # Padding is null in both arrays. A point with no intensity cannot be
    # stored, so it is padding whatever its m/z says.
    mz[np.isnan(intensity)] = np.nan
    return DecodedChunks(
        spectrum_index=np.repeat(_numbers(chunk, "spectrum_index", np.int64), counts),
        mz=mz,
        intensity=intensity,
    )
