"""Which spectra a pixel keeps when it holds spectra of more than one MS level.

Spectra that share a pixel are summed into it. An MS2 spectrum summed with
the MS1 spectrum of its pixel makes a spectrum the instrument never
measured, so a spectrum of level 2 and up is left out where a spectrum of a
lower level sits on the same pixel. Level 0 means the source stated none,
and is read as MS1 there. An MSn spectrum on a pixel of its own is kept:
leaving it out would lose the pixel (D70).

The mzPeak reader and the imzML route apply the same rule through
:func:`msn_left_out`, so the two routes give the same store for one
acquisition. The imzML route reads its levels from the parser, which the
hook in :mod:`thyra.readers.imzml._pyimzml_compat` fills during pyimzml's
one pass over the spectra; the reader and the metadata extractor both ask
:func:`imzml_left_out`, which works the answer out once per parser.
"""

from typing import Any, Optional, Sequence, cast

import numpy as np
from numpy.typing import NDArray

#: The parser attribute that holds one byte per spectrum: the MS level the
#: spectrum states, on itself or through a parameter group it references,
#: and 0 where it states none. Levels above 255 are stored as 255.
MS_LEVELS_ATTRIBUTE = "thyra_ms_levels"

#: The parser attribute :func:`imzml_left_out` keeps its answer in.
_LEFT_OUT_ATTRIBUTE = "thyra_msn_left_out"


def msn_left_out(
    levels: NDArray[np.integer], positions: NDArray[np.integer]
) -> Optional[NDArray[np.bool_]]:
    """Which spectra the per-pixel MSn rule leaves out, or ``None`` for none.

    Args:
        levels: MS level of each spectrum, 0 where none is given.
        positions: ``(n, k)`` position of each spectrum; spectra with equal
            rows share a pixel.

    Returns:
        A boolean mask, ``True`` for each spectrum left out, or ``None``
        when the rule leaves nothing out.
    """
    higher = np.asarray(levels) >= 2
    if not higher.any() or higher.all():
        return None
    _, pixel = np.unique(positions, axis=0, return_inverse=True)
    pixel = np.asarray(pixel, dtype=np.int64).ravel()
    has_first = np.zeros(int(pixel.max()) + 1, dtype=bool)
    has_first[pixel[~higher]] = True
    left_out = higher & has_first[pixel]
    if not left_out.any():
        return None
    return left_out


def imzml_left_out(parser: Any) -> Optional[NDArray[np.bool_]]:
    """The spectra of an imzML parser that the per-pixel MSn rule leaves out.

    The pixel is the spectrum's ``(x, y, z)`` as the parser holds it; a
    spectrum on another plane is on another pixel.

    Args:
        parser: A pyimzml parser.

    Returns:
        A boolean mask over ``parser.coordinates``, ``True`` for each
        spectrum left out, or ``None`` when nothing is left out. A parser
        built without Thyra's hook records no levels, and keeps every
        spectrum.
    """
    own = getattr(parser, "__dict__", {})
    if _LEFT_OUT_ATTRIBUTE in own:
        return cast(Optional[NDArray[np.bool_]], own[_LEFT_OUT_ATTRIBUTE])
    left_out = _work_out(
        own.get(MS_LEVELS_ATTRIBUTE), getattr(parser, "coordinates", None)
    )
    own[_LEFT_OUT_ATTRIBUTE] = left_out
    return left_out


def _work_out(
    levels: Optional[bytearray], coordinates: Optional[Sequence[Sequence[int]]]
) -> Optional[NDArray[np.bool_]]:
    """:func:`msn_left_out` for the levels and coordinates a parser recorded."""
    if not isinstance(levels, (bytes, bytearray)) or coordinates is None:
        return None
    if len(levels) != len(coordinates) or not levels:
        return None
    array = np.frombuffer(bytes(levels), dtype=np.uint8)
    if not (array >= 2).any():
        # The case of nearly every file: no coordinate array is built.
        return None
    positions = np.asarray(coordinates, dtype=np.int64).reshape(len(coordinates), -1)
    return msn_left_out(array, positions)
