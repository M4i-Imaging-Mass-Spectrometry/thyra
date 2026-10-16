"""The ``(x, y)`` to region map of a run, built from its scan list.

The timsTOF and solariX readers, and an mzPeak archive of a Bruker run,
record a region number per scan. A position
measured in two regions is one row in the store, because the converter sums
the spectra of a repeated position (D40). That row takes the region of its
first scan, the same scan whose number becomes its ``acquisition_order``
(D26), and the reader says which regions shared positions (issue #496).
"""

import logging
from typing import Dict, Iterable, Set, Tuple

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
