"""Direct .ibd reads against pyimzml's offset tables, m/z array only.

pyimzml's ``getspectrum`` seeks and decodes *both* binary arrays on every
call. Two hot paths -- the processed-mode mass-axis build and the
processed-mode metadata scan -- need only the m/z array, so going through
``getspectrum`` doubles the bytes read and decoded per spectrum for
nothing. The functions here perform the same seek/read/frombuffer sequence
pyimzml itself performs, for the m/z half alone.

The fast path is gated on the parser carrying pyimzml's *actual* data
structures, checked by type rather than by presence: test doubles and
wrappers (including ``unittest.mock.Mock``, which fabricates any attribute
asked of it) fall back to the documented ``getspectrum`` API.
"""

from typing import Any

import numpy as np
from numpy.typing import NDArray


def has_direct_tables(parser: Any) -> bool:
    """Whether ``parser`` carries pyimzml's real offset tables.

    Args:
        parser: A (possibly duck-typed) ImzML parser.

    Returns:
        True when the m/z-only read in :func:`read_mzs_direct` is legal.
    """
    return (
        isinstance(getattr(parser, "mzOffsets", None), (list, np.ndarray))
        and isinstance(getattr(parser, "mzLengths", None), (list, np.ndarray))
        and isinstance(getattr(parser, "sizeDict", None), dict)
        and getattr(parser, "m", None) is not None
    )


def read_mzs_direct(parser: Any, idx: int) -> NDArray[Any]:
    """Read one spectrum's m/z array without touching its intensities.

    Callers must have checked :func:`has_direct_tables` first.

    Args:
        parser: An initialized pyimzml ImzMLParser.
        idx: Spectrum index.

    Returns:
        The m/z array, exactly as ``getspectrum(idx)[0]`` would return it.
    """
    n = int(parser.mzLengths[idx])
    if n <= 0:
        return np.array([], dtype=np.float64)
    parser.m.seek(parser.mzOffsets[idx])
    data = parser.m.read(n * parser.sizeDict[parser.mzPrecision])
    return np.frombuffer(data, dtype=parser.mzPrecision)


def read_spectrum_mzs_only(parser: Any, idx: int) -> NDArray[Any]:
    """Read one spectrum's m/z values by the cheapest legal route.

    Args:
        parser: A (possibly duck-typed) ImzML parser.
        idx: Spectrum index.

    Returns:
        The m/z array. Falls back to ``getspectrum`` for parsers without
        pyimzml's offset tables.
    """
    if has_direct_tables(parser):
        return read_mzs_direct(parser, idx)
    return parser.getspectrum(idx)[0]
