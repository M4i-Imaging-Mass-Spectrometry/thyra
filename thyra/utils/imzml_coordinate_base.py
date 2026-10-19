"""Which coordinate value an imzML file means by "the first pixel".

The imzML specification says x and y count from 1, and until issue #244 the
reader took it at its word and subtracted a constant 1. Files written
0-based exist, and on one of those the constant produced ``x = -1`` and
``y = -1`` for the first row and column. A negative index is a legal
*negative* numpy index, so the converter's grid guard dropped those
spectra: a 3x3 file at coordinates 0..2 stored 4 rows and reported
"5 spectra sat outside the declared 2x2x1 grid" -- blaming a grid the file
never declared, and exiting 0.

Why not rebase on the observed minimum, the way z does
-----------------------------------------------------

The z rule -- subtract the smallest z the file holds -- rests on an
argument that is good for z and does not carry over: imzML guarantees
nothing about the base and pyimzml is inconsistent about it. It stops at
z because **z has no physical origin and x and y do**. An
acquisition cropped to a region of the slide legitimately starts at
``x = 5``; rebasing on the observed minimum would slide it to ``x = 0``,
changing the grid width, every ``obs["spatial_x"]``, the TIC image extent
and the pixel footprint of a file that converts correctly today. Nothing
in the file distinguishes that acquisition from a 0-based export whose
first column happens to be empty.

So the rule folds down a 0 and nothing else: ``base = min(observed, 1)``.
A file starting at 0 is 0-based and keeps all its pixels; a file starting
at 1, or at 5, is 1-based and is offset by exactly 1, which is what
happened before. The only files whose stored coordinates move are the ones
that were losing a row and a column.

The declared ``IMS:1000042`` / ``IMS:1000043`` pixel counts are the third
option and are deliberately not used: nothing in Thyra reads them today
except ``mzpeak_extractor``, which carries a comment about a declared
extent disagreeing with the coordinates it ships with. Trusting a declared
extent over the coordinates would be a larger change with a wider blast
radius than the defect it fixes.

The imzML extractor reports the source's own coordinate of index 0 on
each axis as ``EssentialMetadata.coordinate_offsets``, which the converter
writes to ``coordinate_systems.global.coordinate_offsets_px`` (D30). On x
and y that is what was subtracted. On z it is too, except for a file that
states no z: pyimzml gives each of its spectra z = 1, that 1 is
subtracted, and the store records 0, as every format without z does.

Whether a file states z is judged over every spectrum, during the one pass
pyimzml makes anyway (D63). A file states position z on every spectrum or
on none; one that mixes the two, or takes z from a parameter group, is
refused.
"""

import logging
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple, Union
from xml.etree import ElementTree  # nosec B405

import numpy as np
from numpy.typing import ArrayLike

from ..errors import ConversionRefused

logger = logging.getLogger(__name__)

#: The base the imzML specification declares for x and y. A file whose
#: smallest coordinate is larger than this is a cropped acquisition, not a
#: different convention, so it keeps this base.
SPEC_BASE = 1

#: Position z. A spectrum may leave it out, and pyimzml then gives it z = 1.
POSITION_Z_ACCESSION = "IMS:1000052"

#: mzML's XML namespace, spelled the way pyimzml spells it.
_MZML_NS = "{http://psi.hupo.org/ms/mzml}"

#: What one spectrum's scan says about position z. The parser hook in
#: :mod:`thyra.readers.imzml._pyimzml_compat` records one of these per
#: spectrum, in the parser attribute named by ``Z_STATES_ATTRIBUTE``.
#: No IMS:1000052 on the scan, or one without a value: pyimzml gives the
#: spectrum its stand-in z of 1.
Z_UNSTATED = 0
#: IMS:1000052 on the scan itself, with a value.
Z_STATED = 1
#: No z on the scan, but a parameter group the scan references has one.
#: pyimzml does not read it there and gives the spectrum its stand-in 1.
Z_FROM_GROUP = 2

#: The parser attribute that holds one ``Z_*`` byte per spectrum.
Z_STATES_ATTRIBUTE = "thyra_z_states"
#: The parser attribute that holds ``(spectrum index, group id)`` for the
#: first spectrum whose z sits in a parameter group.
Z_GROUP_ATTRIBUTE = "thyra_z_group"

#: The rule both the imzML and the mzPeak route refuse a mixed file with.
PLANE_RULE = "A plane is stated for every spectrum or for none."


def coordinate_bases(
    coordinates: Iterable[Sequence[int]],
    log: bool = False,
) -> Tuple[int, int, int]:
    """The ``(x, y, z)`` values in this file that map onto index 0.

    One pass over the coordinate list, which is what makes it worth
    memoising at the call site: ~60 ms at 900k spectra, against a 63 s
    parse.

    Args:
        coordinates: The parser's coordinate list -- ``(x, y, z)`` triples.
        log: Say at INFO when x or y is rebased below the specification's
            base. The reader asks for it, once per file; the metadata
            extractor computes the same bases and stays quiet.

    Returns:
        ``(x_base, y_base, z_base)``. x and y fold a 0 down and otherwise
        keep the specification's base of 1; z is the smallest value present,
        because it has no origin to preserve. ``(1, 1, 0)`` for an empty
        list, so an empty file still normalises the way the spec says.
    """
    min_x: Optional[int] = None
    min_y: Optional[int] = None
    min_z: Optional[int] = None
    for coord in coordinates:
        x, y, z = int(coord[0]), int(coord[1]), int(coord[2])
        if min_x is None or x < min_x:
            min_x = x
        if min_y is None or y < min_y:
            min_y = y
        if min_z is None or z < min_z:
            min_z = z

    if min_x is None or min_y is None or min_z is None:
        return SPEC_BASE, SPEC_BASE, 0

    x_base = min(min_x, SPEC_BASE)
    y_base = min(min_y, SPEC_BASE)

    if log and (x_base < SPEC_BASE or y_base < SPEC_BASE):
        logger.info(
            "Rebasing this imzML: its smallest coordinates are x=%d and "
            "y=%d, and the specification counts from 1. x is rebased on %d "
            "and y on %d, so no row or column is lost.",
            min_x,
            min_y,
            x_base,
            y_base,
        )

    return x_base, y_base, int(min_z)


def recorded_offsets(
    bases: Tuple[int, int, int],
    z_values: ArrayLike,
    imzml_path: Path,
    z_states: Optional[Union[bytes, bytearray]] = None,
) -> Tuple[int, int, int]:
    """The offsets a store records for this file.

    The bases :func:`coordinate_bases` subtracts, except for a z the file
    does not state: that one is recorded as 0, as every format without z
    records it (D30). pyimzml's stand-in of 1 is subtracted all the same,
    so the spectra still sit on index 0 of z.

    Args:
        bases: What :func:`coordinate_bases` returned for the file.
        z_values: The z of every spectrum, as pyimzml reports it.
        imzml_path: Path to the imzML file.
        z_states: One ``Z_*`` byte per spectrum, as the parser hook
            recorded it. ``None`` for a parser built without the hook; the
            file is then judged by :func:`states_z`.

    Returns:
        ``(x, y, z)`` for ``EssentialMetadata.coordinate_offsets``.

    Raises:
        ConversionRefused: If ``z_states`` shows a file that states z on
            some spectra and not others, or takes it from a group.
    """
    x_base, y_base, z_base = bases
    if z_states is not None:
        stated = z_is_stated(z_states, imzml_path)
    else:
        stated = states_z(z_values, imzml_path)
    return (x_base, y_base, z_base if stated else 0)


def z_is_stated(
    z_states: Union[bytes, bytearray],
    imzml_path: Path,
    z_group: Optional[Tuple[int, str]] = None,
) -> bool:
    """Whether the file states position z, judged over every spectrum (D63).

    A spectrum states z when its scan carries ``IMS:1000052`` itself, with
    a value; that is what pyimzml reads. A file states z on every spectrum
    or on none. When some spectra state it and others do not, the others
    would sit on a plane made from pyimzml's stand-in of 1, so the file is
    refused with the sentence the mzPeak route uses. A z held in a
    parameter group is not read by pyimzml at all, which would sum the
    planes, so that file is refused too, by name.

    Args:
        z_states: One ``Z_*`` byte per spectrum.
        imzml_path: Path to the imzML file, for the message.
        z_group: ``(spectrum index, group id)`` of the first spectrum whose
            z sits in a group, if the hook recorded one.

    Returns:
        True when every spectrum states z, False when none does.

    Raises:
        ConversionRefused: If some spectra state z and others do not, or if
            any spectrum takes z from a parameter group.
    """
    states = np.frombuffer(bytes(z_states), dtype=np.uint8)
    grouped = np.flatnonzero(states == Z_FROM_GROUP)
    if grouped.size:
        index, group = z_group if z_group is not None else (int(grouped[0]), "")
        named = f"parameter group {group!r}" if group else "a parameter group"
        raise ConversionRefused(
            f"{imzml_path} gives {grouped.size} of {states.size} spectra "
            f"their position z ({POSITION_Z_ACCESSION}) through a referenced "
            f"parameter group (spectrum {index} through {named}). The imzML "
            f"parser reads z only from the scan itself and would put those "
            f"spectra on one plane. State position z on each scan."
        )
    stated = int(np.count_nonzero(states == Z_STATED))
    if stated == 0:
        return False
    if stated == states.size:
        return True
    raise ConversionRefused(
        f"{imzml_path} gives {stated} of {states.size} spectra a position z "
        f"({POSITION_Z_ACCESSION}) and the others none. {PLANE_RULE}"
    )


def states_z(z_values: ArrayLike, imzml_path: Path) -> bool:
    """Whether the file states position z, for a parser built without the hook.

    The reader always builds its parser with the hook, which records every
    spectrum's z (see :func:`z_is_stated`); this is what a parser built
    elsewhere falls back to. pyimzml gives a spectrum without
    ``IMS:1000052`` the z of 1, so z values that are all 1 belong to a file
    that states z = 1 or to one that states none. Any other z settles it.
    Otherwise the first spectrum's scan decides, read again the way pyimzml
    read it.

    Args:
        z_values: The z of every spectrum, as pyimzml reports it.
        imzml_path: Path to the imzML file.

    Returns:
        False only when the first spectrum states no z and no z differs
        from 1. A file that cannot be read again counts as stating z, so
        its offsets stay what was subtracted.
    """
    if np.any(np.asarray(z_values, dtype=np.int64) != 1):
        return True
    stated = _first_spectrum_states_z(imzml_path)
    return True if stated is None else stated


def _first_spectrum_states_z(imzml_path: Path) -> Optional[bool]:
    """Whether the first spectrum's scan carries position z with a value.

    Costs the header and one spectrum, however large the file is.

    Returns:
        None when the document cannot be read this way.
    """
    try:
        events = ElementTree.iterparse(str(imzml_path), events=("end",))  # nosec B314
        for _event, elem in events:
            if elem.tag != _MZML_NS + "spectrum":
                continue
            scan = elem.find(f"{_MZML_NS}scanList/{_MZML_NS}scan")
            if scan is None:
                return False
            node = scan.find(f'{_MZML_NS}cvParam[@accession="{POSITION_Z_ACCESSION}"]')
            return node is not None and node.get("value") is not None
    except (ElementTree.ParseError, OSError) as e:
        logger.debug(f"Could not read the first spectrum of {imzml_path}: {e}")
        return None
    return False
