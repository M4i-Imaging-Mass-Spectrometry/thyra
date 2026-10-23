"""Position, index and region columns as whole numbers, or a refusal.

Vendor files state pixel positions, spectrum numbers and region numbers in
whatever type their writer chose: ``uint64`` in a Parquet column, ``REAL``
or text in a sqlite one, a float a writer computed as ``0.3 / 0.1``. Each
reader used to cast with ``astype(np.int64)`` or ``int()``, so a float just
below a whole number was cut down onto its neighbour's pixel, a ``uint64``
of ``2**63`` or more wrapped negative, and text or ``NULL`` ended in a
traceback (issue #527).

:func:`whole_numbers` is the one door they go through instead. A value is
accepted when it is a whole number within int64 once floating-point noise
is allowed for: a float within :data:`WHOLE_TOLERANCE` of a whole number is
that number (``2.9999999999999996`` is 3), anything further off is not a
position. Everything else is refused in one line that names the column,
how many values fail and the first of them.

The tolerance is for a column whose type is a float, as mzPeak allows for
positions: there a computed value carries writer noise. A vendor sqlite
column declared INTEGER that holds ``100.9999999`` is damage, not noise,
and its readers pass ``tolerance=0``.
"""

import numbers
from typing import Any, Optional

import numpy as np
from numpy.typing import NDArray

from ..errors import ConversionRefused

__all__ = ["WHOLE_TOLERANCE", "whole_number", "whole_numbers", "whole_sql"]

#: How far a float may sit from a whole number and still be that number.
#: Writer noise (``0.3 / 0.1``) is ~1e-15 of a pixel; a value 1e-6 away
#: is still that pixel by any reading, and 0.5 away is no pixel at all.
WHOLE_TOLERANCE = 1e-6

_INT64_MAX = np.iinfo(np.int64).max
_INT64_MIN = np.iinfo(np.int64).min


def _refuse(what: str, source: Any, n_bad: int, n_all: int, first: Any) -> None:
    where = f" in {source}" if source is not None else ""
    raise ConversionRefused(
        f"{n_bad} of {n_all} {what}{where} are not whole numbers within "
        f"int64 (the first is {first!r}), so no pixel can be placed from them"
    )


def _from_objects(
    values: NDArray[Any], what: str, source: Any, tolerance: float
) -> NDArray[np.int64]:
    """Python objects (a sqlite column): ints, floats, or something else."""
    out = np.empty(values.size, dtype=np.int64)
    bad = []
    for i, value in enumerate(values.tolist()):
        number: Optional[int] = None
        if isinstance(value, bool):
            number = None
        elif isinstance(value, numbers.Integral):
            number = int(value)
        elif isinstance(value, numbers.Real):
            number = _whole_float(float(value), tolerance)
        if number is None or not _INT64_MIN <= number <= _INT64_MAX:
            bad.append(value)
            continue
        out[i] = number
    if bad:
        _refuse(what, source, len(bad), values.size, bad[0])
    return out


def _whole_float(value: float, tolerance: float) -> Optional[int]:
    if not np.isfinite(value):
        return None
    nearest = round(value)
    if abs(value - nearest) > tolerance:
        return None
    return int(nearest)


def whole_number(value: Any) -> Optional[int]:
    """One value of a column declared INTEGER, or ``None`` when it is not whole.

    For a reader that skips a bad record rather than refusing the file: an
    ``int``, or a float with no fraction, within int64. Text, ``None``, a
    bool and a float with any fraction at all give ``None``.
    """
    if isinstance(value, bool):
        return None
    if isinstance(value, numbers.Integral):
        number = int(value)
    elif isinstance(value, numbers.Real):
        found = _whole_float(float(value), 0.0)
        if found is None:
            return None
        number = found
    else:
        return None
    return number if _INT64_MIN <= number <= _INT64_MAX else None


def whole_sql(column: str) -> str:
    """SQL that holds where a sqlite ``column`` stores a whole number.

    An INTEGER, or a REAL with no fraction; text, ``NULL``, a REAL with a
    fraction and a REAL past int64 fail it, the last because ``CAST``
    saturates. The same rule as :func:`whole_number`, for a query that
    should leave such rows out.
    """
    return f"({column} = CAST({column} AS INTEGER))"


def whole_numbers(
    values: Any, what: str, source: Any = None, tolerance: float = WHOLE_TOLERANCE
) -> NDArray[np.int64]:
    """``values`` as int64 whole numbers, or a refusal naming the first bad one.

    Args:
        values: A column of positions, indices or region numbers: any numpy
            array or sequence, of any integer, float or object dtype.
        what: What the values are, plural, as it should read mid-sentence --
            ``"x positions"``, ``"spectrum indices"``.
        source: The file or archive, named in the refusal.
        tolerance: How far a float may sit from a whole number and still
            be read as it; 0 for a column declared as an integer.

    Returns:
        The values as an int64 array of the same length.

    Raises:
        ConversionRefused: If a value is not a whole number within int64:
            ``NULL``, text, ``nan``, ``inf``, a float further than
            :data:`WHOLE_TOLERANCE` from a whole number, or an unsigned
            value of ``2**63`` or more.
    """
    array = np.asarray(values)
    if array.ndim != 1:
        array = array.reshape(-1)
    if array.size == 0:
        return array.astype(np.int64)
    kind = array.dtype.kind
    if kind == "b":
        _refuse(what, source, array.size, array.size, array.flat[0])
    if kind == "i":
        return array.astype(np.int64, copy=False)
    if kind == "u":
        over = np.nonzero(array > np.uint64(_INT64_MAX))[0]
        if over.size:
            _refuse(what, source, over.size, array.size, int(array[over[0]]))
        return array.astype(np.int64)
    if kind == "f":
        as_float = array.astype(np.float64, copy=False)
        nearest = np.rint(as_float)
        with np.errstate(invalid="ignore"):
            # 2**63 is exactly representable and already out of range;
            # every float below it in magnitude fits.
            bad = ~np.isfinite(as_float) | (np.abs(as_float - nearest) > tolerance)
            bad |= (nearest >= 2.0**63) | (nearest < -(2.0**63))
        if bad.any():
            first = np.nonzero(bad)[0][0]
            _refuse(what, source, int(bad.sum()), array.size, float(as_float[first]))
        return nearest.astype(np.int64)
    return _from_objects(array, what, source, tolerance)
