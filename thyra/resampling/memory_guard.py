"""How much memory a mass axis costs, and when a conversion refuses one.

Kept apart from the converter so that :func:`thyra.plan_mass_axis` can ask
the same question before it lays an axis, without importing SpatialData
(issue #473). ``thyra.converters.spatialdata.csc_assembly`` re-exports
every name here.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)

#: Process memory one m/z bin of the summed table costs, everything the
#: conversion builds per bin included. :data:`MAX_COUNT_BYTES` sizes the
#: 4 B/bin count array alone, and that is 2 percent of the real figure:
#: around it sit the axis itself, each table's ``average_spectrum`` and
#: the column sums it is divided out of (8 B/bin each), the ``var`` frame
#: with its index, the ``indptr`` and the write cursor, and anndata's own
#: copies during the write. Measured 2026-09-09 with a six-pixel
#: conversion at a forced bin count, peak working set above the
#: pre-conversion baseline:
#:
#: =============  ==================  ==========
#: bins           peak above baseline bytes/bin
#: =============  ==================  ==========
#: 200,000        56 MB               292
#: 1,000,000      231 MB              243
#: 2,000,000      429 MB              225
#: 5,000,000      966 MB              203
#: 10,000,000     1,802 MB            189
#: =============  ==================  ==========
#:
#: Issue #251 measured 208 B/bin at 20M bins and 337 B/bin at 200M, where
#: a 6-pixel dataset peaked at 67.5 GB. 200 is the figure those agree on.
#:
#: The ``var`` index is over half of it -- 117 B/bin peak at 20M bins --
#: and that is not something a cheaper construction fixes. Building the
#: labels as a numpy ``StringDType`` array rather than a Python list was
#: measured at both 5M and 20M bins: 583 MB against 584 MB, and 2,328 MB
#: against 2,329 MB. On pandas 3 the index becomes an Arrow-backed ``str``
#: either way and the peak is the conversion, not the labels. So the cost
#: is real, it is inherent to a string ``var`` index, and the answer is to
#: refuse an axis that cannot afford it rather than to shave the index.
#:
#: The figure is per bin of *one* table. A multi-slice source converted as
#: 2D writes one table per plane and each carries its own ``var`` copy and
#: its own mean spectrum, so its real per-bin cost is a multiple of this
#: and the guard under-projects it. That predates the per-table mean
#: (#243), which is close to neutral here: it dropped one dataset-wide
#: accumulator and one dataset-wide mean (8 B/bin each) and added 8 B/bin
#: per table.
AXIS_BYTES_PER_BIN = 200

#: Fractions of the machine's free memory the projected mass axis may
#: take before the conversion warns, and before it refuses. Fractions
#: rather than sizes, for the reason the grid's own guard gives: an axis
#: that is routine on a workstation is fatal on a laptop.
AXIS_WARN_FRACTION = 0.25
AXIS_REFUSE_FRACTION = 0.5

#: What to assume is free when the machine will not say. Generous on
#: purpose: a guess must never be the thing that refuses a conversion
#: that would have fitted.
ASSUMED_FREE_GB = 8.0


def available_memory_gb() -> float:
    """The machine's free memory in GB, or :data:`ASSUMED_FREE_GB`."""
    try:
        import psutil

        return float(psutil.virtual_memory().available) / 1024**3
    except Exception as e:  # pragma: no cover - platform-dependent
        logger.debug("Could not read available memory: %s", str(e))
        return ASSUMED_FREE_GB


def projected_axis_gb(n_bins: int) -> float:
    """What a mass axis of this many bins is expected to cost in memory."""
    return float(n_bins) * AXIS_BYTES_PER_BIN / 1024**3


def mass_axis_refusal(
    n_bins: int,
    available_gb: Optional[float] = None,
    levers: Optional[str] = None,
) -> Optional[str]:
    """Why a mass axis this wide cannot be converted here, or ``None``.

    Asked *before* the axis is materialised, from the bin count the
    resampling plan resolved, because building it is already expensive:
    300M bins on a six-pixel dataset were refused in 2.3 s by the count
    array's own ceiling, but the process had reached 7.44 GB getting
    there. 200M bins passed that ceiling entirely -- ``span * 4`` is
    800 MB, under the 1 GiB it allows -- and took 217 s and 67.5 GB
    (issue #251).

    ``available_gb`` overrides the machine's own answer, for tests.
    ``levers`` replaces the closing advice (the CLI options by default), for
    a caller that has other ways to make the axis smaller.
    """
    gb = projected_axis_gb(n_bins)
    free = available_memory_gb() if available_gb is None else float(available_gb)
    if levers is None:
        levers = (
            "resample to fewer mass bins (--resample-bins) or ask for a wider "
            "bin (--resample-width-at-mz)"
        )
    if gb > free * AXIS_REFUSE_FRACTION:
        return (
            f"the mass axis would have {int(n_bins):,} bins, and the "
            f"structures built per bin are projected to need {gb:.1f} GB "
            f"({AXIS_BYTES_PER_BIN} bytes per bin, measured) against "
            f"{free:.1f} GB free, more than the {AXIS_REFUSE_FRACTION:.0%} of "
            f"free memory a conversion may take; {levers}"
        )
    if gb > free * AXIS_WARN_FRACTION:
        logger.warning(
            "The mass axis has %s bins, projected to need %.1f GB of the "
            "%.1f GB free. It will be attempted; %s to make it smaller.",
            f"{int(n_bins):,}",
            gb,
            free,
            levers,
        )
    return None
