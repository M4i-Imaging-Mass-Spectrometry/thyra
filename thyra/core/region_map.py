"""The ``(x, y)`` to region map of a run, built from its scan list.

The timsTOF and solariX readers, and an mzPeak archive of a Bruker run,
record a region number per scan. A position
measured in two regions is one row in the store, because the converter sums
the spectra of a repeated position (D40). That row takes the region of its
first scan, the same scan whose number becomes its ``acquisition_order``
(D26), and the reader says which regions shared positions (issue #496).

:func:`summarize_regions` builds the per-region summary the readers report
through ``get_region_info`` (count, box, name) from the same scan list, and
:func:`checked_region_number` is the type rule for a ``region=`` request.
"""

import logging
import numbers
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set, Tuple

logger = logging.getLogger(__name__)


def build_region_map(
    scans: Iterable[Tuple[int, int, int]],
) -> Dict[Tuple[int, int], int]:
    """Map each position to the region of its first scan.

    Args:
        scans: ``(x, y, region)`` per scan, 0-based positions, in
            acquisition order.

    Returns:
        Position to region number.
    """
    region_map: Dict[Tuple[int, int], int] = {}
    shared: Dict[Tuple[int, int], Set[int]] = {}
    for x, y, region in scans:
        key = (int(x), int(y))
        region = int(region)
        first = region_map.setdefault(key, region)
        if first != region:
            shared.setdefault(key, {first}).add(region)
    if shared:
        _warn_shared_positions(shared)
    return region_map


def _warn_shared_positions(shared: Dict[Tuple[int, int], Set[int]]) -> None:
    """Say which regions measured the same positions."""
    regions = sorted(set().union(*shared.values()))
    example = min(shared)
    logger.warning(
        "%d position(s) were measured in more than one region (regions %s; "
        "the first is x=%d, y=%d). Each is one row that sums the spectra of "
        "every region there, labelled with the region of its first scan. "
        "Convert one region at a time (region=) to keep them apart.",
        len(shared),
        ", ".join(str(r) for r in regions),
        example[0],
        example[1],
    )


def summarize_regions(
    scans: Iterable[Tuple[int, int, int]],
    names: Optional[Mapping[int, str]] = None,
) -> List[Dict[str, Any]]:
    """One summary per region present in a scan list.

    Args:
        scans: ``(x, y, region)`` per scan, 0-based positions. A position
            given twice counts twice.
        names: Region number to name; a region without a non-empty name
            gets no ``name`` key.

    Returns:
        Dicts with ``region_number``, ``n_spectra``, ``bounds`` as
        ``(x_min, y_min, x_max, y_max)`` and, when named, ``name``, sorted
        by descending ``n_spectra`` and then ``region_number``.
    """
    by_region: Dict[int, List[Tuple[int, int]]] = {}
    for x, y, region in scans:
        by_region.setdefault(int(region), []).append((int(x), int(y)))
    info: List[Dict[str, Any]] = []
    for region, members in by_region.items():
        xs = [x for x, _ in members]
        ys = [y for _, y in members]
        entry: Dict[str, Any] = {
            "region_number": region,
            "n_spectra": len(members),
            "bounds": (min(xs), min(ys), max(xs), max(ys)),
        }
        if names and names.get(region):
            entry["name"] = names[region]
        info.append(entry)
    info.sort(key=lambda entry: (-entry["n_spectra"], entry["region_number"]))
    return info


def checked_region_number(region: object) -> Optional[int]:
    """A ``region=`` request that must be a region number, as a plain int.

    Any integer type counts (a number read back from ``obs["region_number"]``
    is a numpy integer), and so does a decimal string. ``None`` stays
    ``None``.

    Raises:
        ValueError: For a bool, a float, any other type, or a string that
            is not a whole decimal number.
    """
    if region is None:
        return None
    if isinstance(region, bool) or not isinstance(region, (numbers.Integral, str)):
        raise ValueError(f"got {type(region).__name__}")
    if isinstance(region, numbers.Integral):
        return int(region)
    text = region.strip()
    if not text.isdecimal():
        raise ValueError(f"got the string {region!r}")
    return int(text)
