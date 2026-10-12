"""Extract the reference outputs of one SCiLS Lab ``.slx`` file.

This is the only part of the parity harness that talks to SCiLS Lab. It
runs under an interpreter that has the ``scilslab`` client installed (the
wheels ship with SCiLS Lab under ``APIClients/Python``), which is not one of
Thyra's dependencies, so it imports nothing from Thyra or from the rest of
the harness. Everything it reads goes through the documented public API of
that client; the files it writes are plain NumPy and JSON, which the
harness then loads with :mod:`tests.parity.scils_reference`.

For a dataset folder ``<root>/<name>/`` holding ``scils/<something>.slx``
it writes, next to the ``.slx``:

``reference.npz``
    ``mz`` (the common mass axis, read from the mean spectrum),
    ``mean`` (the mean spectrum, not normalised), ``mean_tic`` (the mean
    spectrum under SCiLS's TIC normalisation, when the file has one),
    ``spot_id``, ``spot_x``, ``spot_y``, ``spot_z`` (spot coordinates),
    ``spot_sum`` (each spot's sum over the whole axis, interval mode
    ``sum``), ``norm_<key>`` (each spot's normalisation factor, as the
    ratio of the plain to the normalised full-axis sum), and
    ``sample_spot_id`` with ``sample_rebinned`` (those spots' spectra on
    the common axis) and ``sample_raw_mz_<i>`` / ``sample_raw_int_<i>``
    (the same spots' spectra as stored, before rebinning).
``reference.json``
    The server version, the normalisation names, the import parameters
    from a saved ``.sip`` next to the ``.slx`` when there is one, and the
    sizes of everything above.

Usage::

    <scilslab-python> tests/parity/scils_extract.py <root>/<name> [--samples 8]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

# The normalisation names SCiLS Lab 15 reports, mapped to the keys the
# reference files use. An unlisted name is kept under a slug of itself.
_NORMALIZATION_KEYS = {
    "Total Ion Count": "tic",
    "Root Mean Square": "rms",
    "Median": "median",
}


def _slug(name: str) -> str:
    """Turn a normalisation name into a key usable in an ``.npz`` file."""
    return "".join(c.lower() if c.isalnum() else "_" for c in name).strip("_")


def _sample_spots(spot_ids: np.ndarray, spot_sum: np.ndarray, n: int) -> List[int]:
    """Pick the spots whose individual spectra are kept.

    Evenly spaced over the spot order, plus the spot with the largest sum,
    so a run of empty spots at the edge of the raster cannot make every
    sample a zero spectrum.
    """
    if spot_ids.size == 0:
        return []
    positions = np.unique(np.linspace(0, spot_ids.size - 1, n).round().astype(int))
    chosen = [int(spot_ids[p]) for p in positions]
    brightest = int(spot_ids[int(np.argmax(spot_sum))])
    if brightest not in chosen:
        chosen.append(brightest)
    # get_spectra_matrix refuses spot ids that are not in increasing order.
    return sorted(chosen)


def _read_sip(scils_dir: Path) -> Dict[str, Any]:
    """Return the mass-axis and rebinning settings of a saved import, if any."""
    for sip in sorted(scils_dir.glob("*.sip")):
        try:
            project = json.loads(sip.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        for block in project.values():
            if isinstance(block, dict) and "massaxis" in block:
                return {
                    "file": sip.name,
                    "project_type": block.get("type"),
                    "massaxis": block.get("massaxis"),
                    "rebinning": block.get("storage", {}).get("rebinningType"),
                    "baseline": block.get("baselineremoval", {}).get("type"),
                }
    return {}


def extract(dataset_dir: Path, samples: int) -> Path:
    """Write ``reference.npz`` and ``reference.json`` for one dataset.

    Args:
        dataset_dir: The dataset folder; its ``scils/`` holds one ``.slx``.
        samples: How many evenly spaced spots to keep spectra for.

    Returns:
        The path of the written ``reference.npz``.
    """
    from scilslab import LocalSession  # type: ignore[import-not-found]

    scils_dir = dataset_dir / "scils"
    slx_files = sorted(scils_dir.glob("*.slx"))
    if len(slx_files) != 1:
        raise SystemExit(f"{scils_dir}: expected one .slx, found {len(slx_files)}")

    arrays: Dict[str, Any] = {}
    info: Dict[str, Any] = {"slx": slx_files[0].name, "import": _read_sip(scils_dir)}
    with LocalSession(filename=str(slx_files[0]), timeout=600) as session:
        ds = session.dataset_proxy
        info["server_version"] = ds.get_server_version()
        norms = {uid: name for uid, name in ds.get_normalizations().items()}
        info["normalizations"] = sorted(norms.values())

        mean = ds.get_mean_spectrum()
        mz = np.asarray(mean["mz"], dtype=np.float64)
        arrays["mz"] = mz
        arrays["mean"] = np.asarray(mean["intensities"], dtype=np.float64)
        tic_ids = [uid for uid, name in norms.items() if name == "Total Ion Count"]
        if tic_ids:
            mean_tic = ds.get_mean_spectrum(normalization_id=tic_ids[0])
            arrays["mean_tic"] = np.asarray(mean_tic["intensities"], dtype=np.float64)

        spots = ds.get_region_spots("Regions")
        spot_ids = np.asarray(spots["spot_id"], dtype=np.int64)
        arrays["spot_id"] = spot_ids
        for axis in ("x", "y", "z"):
            arrays[f"spot_{axis}"] = np.asarray(spots[axis], dtype=np.float64)

        # One query over the whole axis gives every spot's sum; the same query
        # under a normalisation gives the normalised sum, and their ratio is
        # the factor SCiLS divides that spot by. The bounds are the axis ends
        # exactly: an interval reaching past either end is refused with "The
        # specified interval contains no values".
        lo, hi = float(mz[0]), float(mz[-1])
        plain = ds.get_ion_intensities(lo, hi, mode="sum")
        order = np.argsort(np.asarray(plain.spot_ids))
        plain_ids = np.asarray(plain.spot_ids)[order]
        if not np.array_equal(plain_ids, np.sort(spot_ids)):
            raise SystemExit("ion intensities do not cover the spots of 'Regions'")
        by_id = dict(zip(plain_ids.tolist(), np.asarray(plain.values)[order].tolist()))
        spot_sum = np.array([by_id[int(s)] for s in spot_ids], dtype=np.float64)
        arrays["spot_sum"] = spot_sum
        for uid, name in norms.items():
            key = _NORMALIZATION_KEYS.get(name, _slug(name))
            normed = ds.get_ion_intensities(lo, hi, mode="sum", normalization_id=uid)
            n_by_id = dict(
                zip(
                    np.asarray(normed.spot_ids).tolist(),
                    np.asarray(normed.values).tolist(),
                )
            )
            normed_sum = np.array([n_by_id[int(s)] for s in spot_ids], dtype=np.float64)
            with np.errstate(divide="ignore", invalid="ignore"):
                arrays[f"norm_{key}"] = np.where(
                    normed_sum != 0, spot_sum / normed_sum, np.nan
                )

        chosen = _sample_spots(spot_ids, spot_sum, samples)
        arrays["sample_spot_id"] = np.asarray(chosen, dtype=np.int64)
        if chosen:
            matrix = ds.get_spectra_matrix(chosen)
            if not np.allclose(np.asarray(matrix.mz, dtype=np.float64), mz):
                raise SystemExit(
                    "spectra matrix axis differs from the mean spectrum axis"
                )
            arrays["sample_rebinned"] = np.asarray(matrix.intensities, dtype=np.float64)
            for i, sid in enumerate(chosen):
                raw = ds.get_spectrum(int(sid))
                arrays[f"sample_raw_mz_{i}"] = np.asarray(raw["mz"], dtype=np.float64)
                arrays[f"sample_raw_int_{i}"] = np.asarray(
                    raw["intensities"], dtype=np.float64
                )

    info["sizes"] = {k: list(v.shape) for k, v in arrays.items()}
    out = scils_dir / "reference.npz"
    np.savez_compressed(out, **arrays)
    (scils_dir / "reference.json").write_text(
        json.dumps(info, indent=2, sort_keys=True), encoding="utf-8"
    )
    return out


def main(argv: List[str]) -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description="Extract SCiLS Lab reference outputs.")
    parser.add_argument("dataset_dirs", nargs="+", type=Path)
    parser.add_argument("--samples", type=int, default=8)
    args = parser.parse_args(argv)
    for dataset_dir in args.dataset_dirs:
        print(f"{dataset_dir.name}: {extract(dataset_dir, args.samples)}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
