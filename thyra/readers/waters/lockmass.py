"""Lock-mass correction of Waters m/z values, without writing to the run (D82).

A Waters run can carry a lock-mass function: scans of a reference compound
taken at intervals during the acquisition. MassLynx corrects the m/z of
every other scan against them, and an uncorrected run drifts. On the
SYNAPT XS MALDI image measured for D82 the stored m/z of known brain
lipids ran +11.6 ppm high at the start and +20.6 ppm at the end of a
three-hour raster, while the corrected values stayed within 1 to 3 ppm.

The Waters library does the correction in two steps, both measured from
outside on that run:

1. It measures the reference in the lock-mass scans and turns each
   measurement into one factor, ``reference / measured``, written with the
   scan's retention time to ``lmgt.inf`` **in the run's own folder**.
2. Every later read of that folder multiplies a scan's m/z (and those of
   its drift bins) by the factor interpolated linearly at the scan's
   retention time, holding the first and last factor outside them.
   Reproduced to 0.03 ppm, the float32 rounding of the stored values.

Step 1 is the library's own peak measurement and is kept: Thyra hands it a
throwaway copy of the run that holds only the side files and the lock-mass
function's scans, with sparse placeholders standing in for the other
functions' data, so the library writes its ``lmgt.inf`` there and never
into the user's run. The factors it writes from that copy are the ones it
writes from the full run, digit for digit. Step 2 is done here, so the
correction a store carries is recorded in it rather than left in the
source folder for every later reader to inherit.
"""

import logging
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from ...errors import ConversionRefused

if TYPE_CHECKING:
    from .masslynx_lib import MassLynxLib

logger = logging.getLogger(__name__)

#: The file the Waters library writes the correction to, in the run folder.
LMGT_FILE_NAME = "lmgt.inf"

#: Function data files: what the copy replaces with placeholders. The
#: index and statistics files are small and are copied whole.
_FUNCTION_DATA_SUFFIXES = (".DAT", ".CDT", ".IND")


@dataclass(frozen=True)
class LockMassCorrection:
    """One run's lock-mass factors, and where they came from.

    Attributes:
        reference_mz: The reference m/z the factors were measured against.
        tolerance_da: The search window the library recorded with them, in Da.
        times_min: Retention time of each lock-mass measurement, ascending.
        factors: The m/z multiplier measured at each of those times.
        source: ``"lmgt.inf"`` when the run already carried the correction,
            ``"computed"`` when Thyra had the library measure it.
    """

    reference_mz: float
    tolerance_da: Optional[float]
    times_min: NDArray[np.float64]
    factors: NDArray[np.float64]
    source: str

    def factor_at(self, rt_min: float) -> float:
        """The multiplier for a scan at ``rt_min``: linear, held at both ends."""
        return float(np.interp(rt_min, self.times_min, self.factors))

    def report(self) -> Dict[str, Any]:
        """What the metadata records about the correction."""
        ppm = (self.factors - 1.0) * 1e6
        return {
            "source": self.source,
            "reference_mz": self.reference_mz,
            "tolerance_da": self.tolerance_da,
            "n_measurements": int(self.factors.size),
            "correction_ppm_min": float(ppm.min()),
            "correction_ppm_max": float(ppm.max()),
        }


def find_lmgt(raw_dir: Path) -> Optional[Path]:
    """The run's ``lmgt.inf``, whatever its case, or ``None``."""
    for path in raw_dir.iterdir():
        if path.name.lower() == LMGT_FILE_NAME and path.is_file():
            return path
    return None


def read_lmgt(path: Path, source: str = "lmgt.inf") -> LockMassCorrection:
    """Parse an ``lmgt.inf``.

    The file is plain text: a first line ``1``, then ``reference,tolerance``,
    then one ``retention time (min),factor`` line per lock-mass scan.

    Raises:
        ConversionRefused: When the file is not that shape, so a damaged
            file never silently becomes "no correction".
    """
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    try:
        head = [float(v) for v in lines[1].split(",")]
        rows = np.array(
            [[float(v) for v in line.split(",")] for line in lines[2:]],
            dtype=np.float64,
        )
    except (IndexError, ValueError) as error:
        raise ConversionRefused(
            f"{path} is not a lock-mass correction file Thyra can read ({error})."
        ) from error
    if rows.ndim != 2 or rows.shape[0] == 0 or rows.shape[1] != 2:
        raise ConversionRefused(
            f"{path} holds no lock-mass measurements, so the correction the "
            "Waters library applies from it cannot be known."
        )
    order = np.argsort(rows[:, 0], kind="stable")
    return LockMassCorrection(
        reference_mz=head[0],
        tolerance_da=head[1] if len(head) > 1 else None,
        times_min=rows[order, 0],
        factors=rows[order, 1],
        source=source,
    )


def _is_function_data(name: str) -> bool:
    upper = name.upper()
    return upper.startswith("_FUNC") and upper.endswith(_FUNCTION_DATA_SUFFIXES)


def _function_number(name: str) -> Optional[int]:
    """``_FUNC002.DAT`` -> 2; ``None`` for any other name."""
    stem = name.upper()[len("_FUNC") :].split(".")[0]
    return int(stem) if stem.isdigit() else None


def _sparse_placeholder(path: Path, size: int) -> None:
    """A file of ``size`` bytes that takes no disk space.

    The library checks the other functions' data files are there and the
    size their index says, but reads nothing from them to measure the lock
    mass. On Linux a file extended past its end is sparse already; NTFS
    has to be told first, or the extension is allocated.
    """
    with open(path, "wb") as handle:
        if sys.platform == "win32":
            _mark_sparse(handle.fileno())
        handle.truncate(size)


def _mark_sparse(fd: int) -> None:
    """Set the NTFS sparse flag on an open file (best effort)."""
    if sys.platform != "win32":
        return
    import ctypes
    import msvcrt
    from ctypes import wintypes

    fsctl_set_sparse = 0x000900C4
    returned = wintypes.DWORD()
    ok = ctypes.WinDLL("kernel32").DeviceIoControl(
        wintypes.HANDLE(msvcrt.get_osfhandle(fd)),
        fsctl_set_sparse,
        None,
        0,
        None,
        0,
        ctypes.byref(returned),
        None,
    )
    if not ok:
        logger.debug("Could not mark a lock-mass placeholder sparse")


def _build_measurement_copy(
    raw_dir: Path, lockmass_function: int, target: Path
) -> None:
    """The run as the library needs it to measure the lock mass, in ``target``.

    Every side file is copied, and the lock-mass function's own scans; the
    other functions' data become placeholders of the same size. An
    ``lmgt.inf`` already in the run is left out, so the measurement starts
    from the stored m/z rather than from a previous correction.
    """
    target.mkdir(parents=True)
    lockmass_number = lockmass_function + 1
    for path in raw_dir.iterdir():
        if not path.is_file() or path.name.lower() == LMGT_FILE_NAME:
            continue
        if not _is_function_data(path.name):
            shutil.copy2(path, target / path.name)
        elif _function_number(
            path.name
        ) == lockmass_number and path.name.upper().endswith(".DAT"):
            shutil.copy2(path, target / path.name)
        else:
            _sparse_placeholder(target / path.name, path.stat().st_size)


def measure_lock_mass(
    ml: "MassLynxLib",
    raw_dir: Path,
    lockmass_function: int,
    reference_mz: float,
) -> LockMassCorrection:
    """Have the Waters library measure the lock mass, on a copy of the run.

    Args:
        ml: The loaded library.
        raw_dir: The run. Nothing is written to it.
        lockmass_function: The 0-based function holding the reference scans.
        reference_mz: The m/z of the reference ion.

    Raises:
        ConversionRefused: When the library finds no usable reference peak.
    """
    with tempfile.TemporaryDirectory(prefix="thyra-lockmass-") as scratch:
        copy = Path(scratch) / raw_dir.name
        _build_measurement_copy(raw_dir, lockmass_function, copy)
        handle = ml.open_file(str(copy))
        try:
            applied = ml.apply_custom_lockmass_correction(handle, reference_mz)
        finally:
            ml.close_file(handle)
        written = find_lmgt(copy)
        if not applied or written is None:
            raise ConversionRefused(
                f"The Waters library found no lock-mass peak at m/z "
                f"{reference_mz:g} in function {lockmass_function} of "
                f"{raw_dir.name}. Check the reference m/z: it is the ion "
                "measured in the lock-mass scans, for example a red "
                "phosphorus cluster anion on a MALDI run or leucine "
                "enkephalin on an electrospray one."
            )
        return read_lmgt(written, source="computed")


# ----------------------------------------------------------------------
# What an uncorrected run's reference scans say (D82)
#
# Thyra does not pick a reference for the user: one that matches the wrong
# peak would shift every m/z by up to the search window, as badly as no
# correction and without a sign of it. It does say what the reference scans
# show, so the warning about an uncorrected run carries facts: how far
# their strongest peak moved over the run, which needs no identification,
# and what the scans look like when that is unambiguous.
# ----------------------------------------------------------------------

#: Monoisotopic mass of phosphorus and of an electron, in Da.
_PHOSPHORUS = 30.97376163
_ELECTRON = 0.00054858

#: Leucine enkephalin, the usual electrospray lock-mass reference:
#: [M-H]- in negative mode (True), [M+H]+ in positive mode (False).
_LEUCINE_ENKEPHALIN = {False: 556.2766, True: 554.2620}

#: How far from a reference ion an uncorrected peak may sit and still be
#: that ion. Uncorrected runs measured for D82 were 8 to 28 ppm off.
_MATCH_PPM = 50.0

#: How tightly the red phosphorus ions must agree with each other to be
#: named: they drift together, so their offsets spread little.
_LADDER_SPREAD_PPM = 3.0

#: Red phosphorus is named only on a ladder of at least this many ions.
_LADDER_MIN_IONS = 5

#: Two lock-mass scans further apart than this (min) are different
#: measurements of the drift.
_BLOCK_GAP_MIN = 1.0

_Scan = Tuple[float, NDArray[np.float64], NDArray[np.float64]]


@dataclass(frozen=True)
class LockMassSurvey:
    """What an uncorrected run's lock-mass scans show.

    Attributes:
        peak_mz: The strongest peak of the reference scans, as stored.
        drift_ppm: How far that peak moved from the first group of
            reference scans to the last, in ppm (last minus first).
        span_ppm: The full range of its position over the groups, in ppm.
        n_groups: How many separate groups of reference scans there were.
        looks_like: ``"red phosphorus"``, ``"leucine enkephalin"``, or
            ``None`` when the scans match neither unambiguously.
        suggested_mz: The reference m/z to pass, when ``looks_like`` names one.
        matched_ions: How many of its ions were found.
        offset_ppm: How far the stored m/z of those ions sit from their
            theoretical masses: the size of the uncorrected error *if* the
            identification is right.
    """

    peak_mz: float
    drift_ppm: float
    span_ppm: float
    n_groups: int
    looks_like: Optional[str]
    suggested_mz: Optional[float]
    matched_ions: int
    offset_ppm: Optional[float]

    def report(self) -> Dict[str, Any]:
        """What the metadata records about the survey."""
        return {
            "reference_peak_mz": self.peak_mz,
            "reference_peak_drift_ppm": self.drift_ppm,
            "reference_peak_span_ppm": self.span_ppm,
            "reference_groups": self.n_groups,
            "looks_like": self.looks_like,
            "suggested_reference_mz": self.suggested_mz,
            "matched_ions": self.matched_ions,
            "offset_from_theory_ppm": self.offset_ppm,
        }


def _summed(scans: Sequence[_Scan]) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The scans summed point by point (uncorrected scans share one m/z grid)."""
    mz = np.concatenate([scan[1] for scan in scans])
    intensity = np.concatenate([scan[2] for scan in scans])
    grid, where = np.unique(np.round(mz, 5), return_inverse=True)
    return grid, np.bincount(where, weights=intensity, minlength=grid.size)


def _apex_near(
    mz: NDArray[np.float64],
    intensity: NDArray[np.float64],
    target: float,
    window_ppm: float,
) -> Optional[Tuple[float, float]]:
    """Centroid and height of the tallest peak within ``window_ppm`` of target."""
    half = target * window_ppm * 1e-6
    lo, hi = np.searchsorted(mz, [target - half, target + half])
    if hi <= lo or float(intensity[lo:hi].max()) <= 0:
        return None
    apex = int(lo + np.argmax(intensity[lo:hi]))
    a, b = max(apex - 3, 0), min(apex + 4, mz.size)
    weights = intensity[a:b]
    return float((mz[a:b] * weights).sum() / weights.sum()), float(intensity[apex])


def _recognise(
    mz: NDArray[np.float64], intensity: NDArray[np.float64], negative: bool
) -> Tuple[Optional[str], Optional[float], int, Optional[float]]:
    """Name the reference only when the summed scans leave no doubt."""
    floor = 0.01 * float(intensity.max())
    sign = 1.0 if negative else -1.0
    ladder = []
    for n in range(3, 61, 2):
        theory = n * _PHOSPHORUS + sign * _ELECTRON
        if theory < mz[0] or theory > mz[-1]:
            continue
        found = _apex_near(mz, intensity, theory, _MATCH_PPM)
        if found is not None and found[1] >= floor:
            ladder.append((theory, (found[0] - theory) / theory * 1e6, found[1]))
    if len(ladder) >= _LADDER_MIN_IONS:
        offsets = np.array([ion[1] for ion in ladder])
        if float(offsets.std()) <= _LADDER_SPREAD_PPM:
            strongest = max(ladder, key=lambda ion: ion[2])
            median = float(np.median(offsets))
            return "red phosphorus", strongest[0], len(ladder), median
    theory = _LEUCINE_ENKEPHALIN[negative]
    found = _apex_near(mz, intensity, theory, _MATCH_PPM)
    if found is not None and found[1] >= 0.5 * float(intensity.max()):
        offset = (found[0] - theory) / theory * 1e6
        return "leucine enkephalin", theory, 1, offset
    return None, None, 0, None


def survey_reference_scans(
    scans: Sequence[_Scan], negative: bool
) -> Optional[LockMassSurvey]:
    """What the lock-mass scans of an uncorrected run show.

    Args:
        scans: ``(retention time in min, m/z, intensity)`` of each reference
            scan, as stored (uncorrected), in acquisition order.
        negative: Whether the run was acquired in negative mode.

    Returns:
        The survey, or ``None`` when the scans hold no signal.
    """
    kept = [scan for scan in scans if scan[1].size and float(np.sum(scan[2])) > 0]
    if not kept:
        return None
    mz, intensity = _summed(kept)
    peak = _apex_near(mz, intensity, float(mz[int(np.argmax(intensity))]), 20.0)
    if peak is None:
        return None
    peak_mz = peak[0]

    groups: List[List[_Scan]] = [[kept[0]]]
    for previous, scan in zip(kept, kept[1:]):
        if scan[0] - previous[0] > _BLOCK_GAP_MIN:
            groups.append([])
        groups[-1].append(scan)
    positions = []
    for group in groups:
        group_mz, group_intensity = _summed(group)
        found = _apex_near(group_mz, group_intensity, peak_mz, _MATCH_PPM)
        if found is not None:
            positions.append(found[0])
    drift = span = 0.0
    if len(positions) > 1:
        drift = (positions[-1] - positions[0]) / peak_mz * 1e6
        span = (max(positions) - min(positions)) / peak_mz * 1e6

    looks_like, suggested, matched, offset = _recognise(mz, intensity, negative)
    return LockMassSurvey(
        peak_mz=peak_mz,
        drift_ppm=drift,
        span_ppm=span,
        n_groups=len(groups),
        looks_like=looks_like,
        suggested_mz=suggested,
        matched_ions=matched,
        offset_ppm=offset,
    )
