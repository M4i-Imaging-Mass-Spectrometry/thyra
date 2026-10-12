"""Load the reference outputs :mod:`tests.parity.scils_extract` wrote.

Plain NumPy and JSON: loading needs neither SCiLS Lab nor its client.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class SampleSpectrum:
    """One spot's spectrum as SCiLS Lab returns it.

    Attributes:
        spot_id: The SCiLS spot id.
        rebinned: The spectrum on the common axis.
        raw_mz: The spectrum's own m/z values, before rebinning.
        raw_intensity: The intensities at ``raw_mz``.
    """

    spot_id: int
    rebinned: NDArray[np.float64] = field(repr=False)
    raw_mz: NDArray[np.float64] = field(repr=False)
    raw_intensity: NDArray[np.float64] = field(repr=False)


@dataclass(frozen=True)
class ScilsReference:
    """Everything extracted from one ``.slx``.

    Attributes:
        mz: The common mass axis.
        mean: The mean spectrum, not normalised.
        mean_tic: The mean spectrum under TIC normalisation, if present.
        spot_id, spot_x, spot_y: The spots and their coordinates.
        spot_sum: Each spot's sum over the whole axis.
        normalizations: Per-spot factor for each normalisation, keyed
            ``"tic"``, ``"rms"``, ``"median"``.
        samples: A handful of individual spectra.
        info: The JSON side: versions, import parameters, sizes.
    """

    mz: NDArray[np.float64] = field(repr=False)
    mean: NDArray[np.float64] = field(repr=False)
    mean_tic: Any = field(repr=False)
    spot_id: NDArray[np.int64] = field(repr=False)
    spot_x: NDArray[np.float64] = field(repr=False)
    spot_y: NDArray[np.float64] = field(repr=False)
    spot_sum: NDArray[np.float64] = field(repr=False)
    normalizations: Dict[str, NDArray[np.float64]] = field(repr=False)
    samples: List[SampleSpectrum] = field(repr=False)
    info: Dict[str, Any] = field(repr=False)

    @property
    def import_axis_type(self) -> str:
        """The axis type the saved import parameters name, or ``""``."""
        return str((self.info.get("import") or {}).get("massaxis", {}).get("type", ""))

    @property
    def import_rebinning(self) -> str:
        """The rebinning the saved import parameters name, or ``""``."""
        return str((self.info.get("import") or {}).get("rebinning") or "")


def load_reference(scils_dir: Path) -> ScilsReference:
    """Read ``reference.npz`` and ``reference.json`` from a ``scils/`` folder."""
    with np.load(scils_dir / "reference.npz") as npz:
        arrays = {key: npz[key] for key in npz.files}
    info = json.loads((scils_dir / "reference.json").read_text(encoding="utf-8"))
    sample_ids = arrays.get("sample_spot_id", np.zeros(0, dtype=np.int64))
    rebinned = arrays.get("sample_rebinned")
    samples = [
        SampleSpectrum(
            spot_id=int(sid),
            rebinned=np.asarray(rebinned[i], dtype=np.float64),
            raw_mz=arrays[f"sample_raw_mz_{i}"],
            raw_intensity=arrays[f"sample_raw_int_{i}"],
        )
        for i, sid in enumerate(sample_ids)
    ]
    return ScilsReference(
        mz=arrays["mz"],
        mean=arrays["mean"],
        mean_tic=arrays.get("mean_tic"),
        spot_id=arrays["spot_id"],
        spot_x=arrays["spot_x"],
        spot_y=arrays["spot_y"],
        spot_sum=arrays["spot_sum"],
        normalizations={
            key[len("norm_") :]: arrays[key]
            for key in arrays
            if key.startswith("norm_")
        },
        samples=samples,
        info=info,
    )
