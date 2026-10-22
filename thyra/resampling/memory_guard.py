"""How much memory a mass axis costs, and when a conversion refuses one.

Kept apart from the converter so that :func:`thyra.plan_mass_axis` can ask
the same question before it lays an axis, without importing SpatialData
(issue #473). ``thyra.converters.spatialdata.csc_assembly`` re-exports
every name here.
"""

import logging
import sys
from pathlib import Path
from typing import Callable, List, Optional, Tuple

logger = logging.getLogger(__name__)

#: Bytes the per-key count array of a table's first pass may take. It is
#: ``span * 4`` with ``span`` the key space -- ``mass axis bins x mobility
#: channels`` for the grid: 4 bytes x bins x channels. The default timsTOF
#: axis grows with the mass range (599,146 bins on a usual 50-1000 m/z run),
#: so 512 channels there need about 1.2 GB. A raw, unresampled axis of
#: several million bins would push it further, and that is a source to resample rather
#: than a table to build: the same axis would also blow the feature
#: ceiling. Kept here rather than in ``csc_assembly`` (which re-exports
#: it), so :func:`mass_axis_refusal` can apply it before an axis is built.
MAX_COUNT_BYTES = 1 << 30

#: Bytes one key costs in that count array (a ``uint32``).
COUNT_BYTES_PER_KEY = 4

#: The most keys one table can count: 268,435,456. The summed table's
#: keys are its m/z bins, so this is also the widest mass axis any
#: conversion can take, however much memory is free (issue #528, F59).
MAX_COUNT_KEYS = MAX_COUNT_BYTES // COUNT_BYTES_PER_KEY

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
#: 2D writes one table per plane, and each extra one costs
#: :data:`TABLE_BYTES_PER_BIN` more; a source with two regions or more
#: costs :data:`REGION_BYTES_PER_BIN` per region on top (issue #528).
AXIS_BYTES_PER_BIN = 200

#: What each table past the first adds per bin, on the 2D route of a
#: multi-plane source: its own ``var`` copy, mean spectrum, column
#: pointers and write cursor. Measured by the 4.3.3 stress test (F116) as
#: the peak at 2M bins less the peak at 10k: 200 B/bin for one plane, 511
#: for ten. 36 is the slope those give.
TABLE_BYTES_PER_BIN = 36

#: What each region adds per bin when a source has two or more: its total
#: intensity and its mean spectrum, 8 B/bin each, plus what building them
#: costs. Measured at 19.5 B/bin per region on a six-region timsTOF copy at
#: 5M bins (F148). One region builds neither array (D57).
REGION_BYTES_PER_BIN = 20

#: Fractions of the free memory (:func:`available_memory_gb`) the
#: projected mass axis may take before the conversion warns, and before it
#: refuses. Fractions rather than sizes, for the reason the grid's own
#: guard gives: an axis that is routine on a workstation is fatal on a
#: laptop.
AXIS_WARN_FRACTION = 0.25
AXIS_REFUSE_FRACTION = 0.5

#: What to assume is free when the machine will not say. Generous on
#: purpose: a guess must never be the thing that refuses a conversion
#: that would have fitted.
ASSUMED_FREE_GB = 8.0

#: A limit at or above this is "no limit": cgroup v1 writes the largest
#: page-aligned ``int64`` when none is set.
_NO_LIMIT = 1 << 60


def available_memory_gb() -> float:
    """The memory this process can still take, in GB.

    The machine's free memory, or less when the process runs under a
    limit of its own: a Windows job object, a Linux cgroup (a container or
    a cluster scheduler) or an address-space limit. Reading the machine
    alone admitted an axis that a 2 GB job limit then failed with an
    allocation traceback (issue #528, F194). :data:`ASSUMED_FREE_GB` when
    the machine will not say.
    """
    try:
        import psutil

        free = float(psutil.virtual_memory().available)
    except Exception as e:  # pragma: no cover - platform-dependent
        logger.debug("Could not read available memory: %s", str(e))
        free = ASSUMED_FREE_GB * 1024**3
    headroom = process_memory_headroom()
    if headroom is not None:
        free = min(free, float(headroom))
    return max(free, 0.0) / 1024**3


def process_memory_headroom() -> Optional[int]:
    """Bytes this process may still allocate under its own limits, or ``None``.

    ``None`` when no limit applies or none could be read. Where several
    apply, the smallest headroom is the answer.
    """
    probes: Tuple[Callable[[], Optional[int]], ...] = (
        _job_object_headroom,
        _cgroup_headroom,
        _address_space_headroom,
    )
    found: List[int] = []
    for probe in probes:
        try:
            headroom = probe()
        except Exception as e:  # a probe must never stop a conversion
            logger.debug("Could not read a process memory limit: %s", str(e))
            continue
        if headroom is not None:
            found.append(max(int(headroom), 0))
    return min(found) if found else None


def _job_object_headroom() -> Optional[int]:
    """What the Windows job this process runs in still allows, or ``None``.

    A job may cap each process's committed memory, the whole job's, or
    both. Committed memory is what such a limit counts, so the process's
    own use is its ``private`` bytes and the job's is its ``JobMemory``.
    """
    if sys.platform != "win32":
        return None
    import ctypes
    from ctypes import wintypes

    import psutil

    class _Basic(ctypes.Structure):
        _fields_ = [
            ("per_process_time", ctypes.c_int64),
            ("per_job_time", ctypes.c_int64),
            ("flags", wintypes.DWORD),
            ("min_working_set", ctypes.c_size_t),
            ("max_working_set", ctypes.c_size_t),
            ("active_processes", wintypes.DWORD),
            ("affinity", ctypes.c_size_t),
            ("priority_class", wintypes.DWORD),
            ("scheduling_class", wintypes.DWORD),
        ]

    class _Extended(ctypes.Structure):
        _fields_ = [
            ("basic", _Basic),
            ("io", ctypes.c_ulonglong * 6),
            ("process_limit", ctypes.c_size_t),
            ("job_limit", ctypes.c_size_t),
            ("peak_process", ctypes.c_size_t),
            ("peak_job", ctypes.c_size_t),
        ]

    class _Usage(ctypes.Structure):
        _fields_ = [
            ("job_memory", ctypes.c_ulonglong),
            ("peak_job_memory", ctypes.c_ulonglong),
        ]

    kernel32 = ctypes.WinDLL("kernel32")
    kernel32.GetCurrentProcess.restype = wintypes.HANDLE
    kernel32.IsProcessInJob.argtypes = [
        wintypes.HANDLE,
        wintypes.HANDLE,
        ctypes.POINTER(wintypes.BOOL),
    ]
    kernel32.QueryInformationJobObject.argtypes = [
        wintypes.HANDLE,
        ctypes.c_int,
        ctypes.c_void_p,
        wintypes.DWORD,
        ctypes.POINTER(wintypes.DWORD),
    ]
    in_job = wintypes.BOOL(False)
    process = kernel32.GetCurrentProcess()
    if not kernel32.IsProcessInJob(process, None, ctypes.byref(in_job)):
        return None
    if not in_job.value:
        return None
    info = _Extended()
    # 9 is JobObjectExtendedLimitInformation; a NULL handle names the
    # job of the calling process.
    if not kernel32.QueryInformationJobObject(
        None, 9, ctypes.byref(info), ctypes.sizeof(info), None
    ):
        return None
    own = int(psutil.Process().memory_info().private)
    found: List[int] = []
    if info.basic.flags & 0x100:  # JOB_OBJECT_LIMIT_PROCESS_MEMORY
        found.append(int(info.process_limit) - own)
    if info.basic.flags & 0x200:  # JOB_OBJECT_LIMIT_JOB_MEMORY
        usage = _Usage()
        # 28 is JobObjectMemoryUsageInformation, Windows 10 and later;
        # without it, this process is counted alone.
        if kernel32.QueryInformationJobObject(
            None, 28, ctypes.byref(usage), ctypes.sizeof(usage), None
        ):
            used = int(usage.job_memory)
        else:
            used = own
        found.append(int(info.job_limit) - used)
    return min(found) if found else None


def _cgroup_headroom(
    proc_cgroup: Path = Path("/proc/self/cgroup"),
    root: Path = Path("/sys/fs/cgroup"),
) -> Optional[int]:
    """What this process's cgroup, and each cgroup above it, still allows.

    cgroup v2 (``memory.max`` and ``memory.current``) and v1
    (``memory.limit_in_bytes`` and ``memory.usage_in_bytes``). Usage
    counts the page cache, and a conversion reading a large source fills
    it. The kernel takes that back before it fails an allocation, so the
    reclaimable part (``inactive_file``) is not counted as used: without
    that, a cgroup whose cache had filled would read as full and refuse a
    conversion that fits. ``None`` off Linux, or with no limit set.
    """
    if not proc_cgroup.is_file():
        return None
    found: List[int] = []
    for line in proc_cgroup.read_text().splitlines():
        parts = line.strip().split(":", 2)
        if len(parts) != 3:
            continue
        _, controllers, path = parts
        if controllers == "":
            found.extend(
                _cgroup_chain(
                    root, path, "memory.max", "memory.current", "inactive_file"
                )
            )
        elif "memory" in controllers.split(","):
            found.extend(
                _cgroup_chain(
                    root / "memory",
                    path,
                    "memory.limit_in_bytes",
                    "memory.usage_in_bytes",
                    "total_inactive_file",
                )
            )
    return min(found) if found else None


def _cgroup_chain(
    base: Path, path: str, limit_name: str, usage_name: str, cache_name: str
) -> List[int]:
    """The headroom at each level from ``base/path`` up to ``base`` that sets a limit.

    Inside a container the path the process reports may not exist under
    the mount, whose root is then the container's own cgroup; such a level
    simply has no files to read.
    """
    folder = base.joinpath(*[p for p in path.strip("/").split("/") if p])
    found: List[int] = []
    while True:
        limit = _read_int(folder / limit_name)
        usage = _read_int(folder / usage_name)
        if limit is not None and usage is not None and limit < _NO_LIMIT:
            cache = _read_stat(folder / "memory.stat", cache_name) or 0
            found.append(limit - max(usage - cache, 0))
        if folder == base or base not in folder.parents:
            return found
        folder = folder.parent


def _read_int(path: Path) -> Optional[int]:
    """An integer file of the cgroup filesystem; ``None`` when absent or ``max``."""
    try:
        text = path.read_text().strip()
    except OSError:
        return None
    return int(text) if text.isdigit() else None


def _read_stat(path: Path, key: str) -> Optional[int]:
    """One counter of a ``memory.stat`` file, or ``None``."""
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return None
    for line in lines:
        name, _, value = line.partition(" ")
        if name == key and value.strip().isdigit():
            return int(value.strip())
    return None


def _address_space_headroom() -> Optional[int]:
    """What ``RLIMIT_AS`` still lets this process map, or ``None``.

    Off Windows only. Memory maps count against it too, and a conversion
    maps its matrix files.
    """
    if sys.platform != "win32":
        import resource

        import psutil

        soft, _hard = resource.getrlimit(resource.RLIMIT_AS)
        if soft == resource.RLIM_INFINITY or soft < 0:
            return None
        return int(soft) - int(psutil.Process().memory_info().vms)
    return None


def axis_bytes_per_bin(n_tables: int = 1, n_regions: int = 0) -> int:
    """What one m/z bin costs a conversion of this many tables and regions."""
    extra_tables = max(int(n_tables) - 1, 0)
    regions = int(n_regions) if int(n_regions) >= 2 else 0
    return (
        AXIS_BYTES_PER_BIN
        + extra_tables * TABLE_BYTES_PER_BIN
        + regions * REGION_BYTES_PER_BIN
    )


def projected_axis_gb(n_bins: int, n_tables: int = 1, n_regions: int = 0) -> float:
    """What a mass axis of this many bins is expected to cost in memory."""
    return float(n_bins) * axis_bytes_per_bin(n_tables, n_regions) / 1024**3


def _per_bin_text(n_tables: int, n_regions: int) -> str:
    """How the per-bin figure is made up, for a message."""
    total = axis_bytes_per_bin(n_tables, n_regions)
    parts = []
    if int(n_tables) > 1:
        parts.append(
            f"{TABLE_BYTES_PER_BIN} for each of {int(n_tables) - 1} more tables"
        )
    if int(n_regions) >= 2:
        parts.append(f"{REGION_BYTES_PER_BIN} for each of {int(n_regions)} regions")
    if not parts:
        return f"{total} bytes per bin, measured"
    return (
        f"{total} bytes per bin, measured: {AXIS_BYTES_PER_BIN} for one table, "
        + " and ".join(parts)
    )


def mass_axis_refusal(
    n_bins: int,
    available_gb: Optional[float] = None,
    levers: Optional[str] = None,
    n_tables: int = 1,
    n_regions: int = 0,
) -> Optional[str]:
    """Why a mass axis this wide cannot be converted here, or ``None``.

    Asked *before* the axis is materialised, from the bin count the
    resampling plan resolved, because building it is already expensive:
    300M bins on a six-pixel dataset were refused in 2.3 s by the count
    array's own ceiling, but the process had reached 7.44 GB getting
    there. 200M bins passed that ceiling entirely -- ``span * 4`` is
    800 MB, under the 1 GiB it allows -- and took 217 s and 67.5 GB
    (issue #251). That ceiling is checked here first, so an axis past
    :data:`MAX_COUNT_KEYS` is refused before it is built however much
    memory is free; it used to be found by the count array, as a bare
    ``MemoryError`` (issue #528).

    ``n_tables`` and ``n_regions`` are the tables the conversion writes
    and the regions it averages, where the caller knows them; each adds to
    the per-bin cost (:func:`axis_bytes_per_bin`). ``available_gb``
    overrides the process's own answer, for tests and for a caller that
    measured it before the conversion started allocating. ``levers``
    replaces the closing advice (the CLI options by default), for a caller
    that has other ways to make the axis smaller.
    """
    if levers is None:
        levers = (
            "resample to fewer mass bins (--resample-bins) or ask for a wider "
            "bin (--resample-width-at-mz)"
        )
    if int(n_bins) > MAX_COUNT_KEYS:
        return (
            f"the mass axis would have {int(n_bins):,} bins, more than the "
            f"{MAX_COUNT_KEYS:,} a table can count; {levers}"
        )
    gb = projected_axis_gb(n_bins, n_tables, n_regions)
    free = available_memory_gb() if available_gb is None else float(available_gb)
    if gb > free * AXIS_REFUSE_FRACTION:
        return (
            f"the mass axis would have {int(n_bins):,} bins, and the "
            f"structures built per bin are projected to need {gb:.1f} GB "
            f"({_per_bin_text(n_tables, n_regions)}) against "
            f"{free:.1f} GB free, more than the {AXIS_REFUSE_FRACTION:.0%} of "
            f"free memory a conversion may take; {levers}"
        )
    if gb > free * AXIS_WARN_FRACTION:
        logger.warning(
            "The mass axis has %s bins, projected to need %.1f GB (%s) of the "
            "%.1f GB free. It will be attempted; %s to make it smaller.",
            f"{int(n_bins):,}",
            gb,
            _per_bin_text(n_tables, n_regions),
            free,
            levers,
        )
    return None
