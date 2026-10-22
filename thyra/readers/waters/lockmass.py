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
from typing import TYPE_CHECKING, Any, Dict, Optional

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
