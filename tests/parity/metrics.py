"""The numbers a parity run reports, as plain functions of NumPy arrays.

Nothing here knows about SCiLS Lab, Thyra or a store: every function takes
arrays and returns floats or small dicts, so each one can be checked on
synthetic input (``tests/unit/tools/test_parity_harness.py``) without any
reference data. The thresholds these numbers are held to live in
:mod:`tests.parity.thresholds`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

# The spacing law of an axis is read from how its bin width grows with m/z:
# width ~ m/z ** k. These are the exponents of the laws both programs name.
# SCiLS Lab calls them Constant, Axial TOF, Orthogonal TOF, Orbitrap and
# MRMS; Thyra calls them constant, linear_tof, reflector_tof, orbitrap and
# fticr.
LAW_EXPONENTS: Dict[str, float] = {
    "constant": 0.0,
    "linear_tof": 0.5,
    "reflector_tof": 1.0,
    "orbitrap": 1.5,
    "fticr": 2.0,
}


def spacing_exponent(mz: FloatArray) -> float:
    """Fit ``k`` in ``bin width ~ m/z ** k`` over the whole axis.

    A least-squares line through ``log(width)`` against ``log(m/z)``. The
    widths are the differences between neighbouring axis values, placed at
    their midpoints.

    Args:
        mz: The axis, strictly increasing, at least three values.

    Returns:
        The fitted exponent.
    """
    mz = np.asarray(mz, dtype=np.float64)
    widths = np.diff(mz)
    centres = 0.5 * (mz[1:] + mz[:-1])
    keep = widths > 0
    slope, _ = np.polyfit(np.log(centres[keep]), np.log(widths[keep]), 1)
    return float(slope)


def nearest_law(exponent: float) -> str:
    """Name the law whose exponent is closest to a fitted one."""
    return min(LAW_EXPONENTS, key=lambda name: abs(LAW_EXPONENTS[name] - exponent))


def axis_agreement(reference: FloatArray, candidate: FloatArray) -> Dict[str, float]:
    """Compare two mass axes.

    When the bin counts match, axis values are compared index by index.
    When they do not, each reference value is compared with the nearest
    candidate value instead, which still says how far apart the two grids
    sit but no longer says anything about bin identity.

    Args:
        reference: The reference axis (SCiLS Lab).
        candidate: The axis under test (Thyra).

    Returns:
        ``bins_reference``, ``bins_candidate``, ``bins_equal`` (1.0 or
        0.0), ``max_abs_da``, ``max_ppm``, ``median_ppm``,
        ``first_ppm`` and ``last_ppm`` (signed, candidate minus
        reference), ``max_ppm_in_bins`` (the largest offset as a fraction
        of the reference bin width there), and the two fitted spacing
        exponents.
    """
    ref = np.asarray(reference, dtype=np.float64)
    cand = np.asarray(candidate, dtype=np.float64)
    if ref.size == cand.size:
        paired = cand
    else:
        idx = np.clip(np.searchsorted(cand, ref), 1, cand.size - 1)
        left, right = cand[idx - 1], cand[idx]
        paired = np.where(np.abs(ref - left) <= np.abs(right - ref), left, right)
    diff = paired - ref
    ppm = diff / ref * 1e6
    widths = np.gradient(ref)
    return {
        "bins_reference": float(ref.size),
        "bins_candidate": float(cand.size),
        "bins_equal": float(ref.size == cand.size),
        "max_abs_da": float(np.max(np.abs(diff))),
        "max_ppm": float(np.max(np.abs(ppm))),
        "median_ppm": float(np.median(np.abs(ppm))),
        "first_ppm": float(ppm[0]),
        "last_ppm": float(ppm[-1]),
        "max_ppm_in_bins": float(np.max(np.abs(diff) / widths)),
        "exponent_reference": spacing_exponent(ref),
        "exponent_candidate": spacing_exponent(cand),
    }


def _pearson(a: FloatArray, b: FloatArray) -> float:
    """Pearson r, or NaN when either side is constant."""
    if a.size < 2 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def intensity_agreement(
    reference: FloatArray, candidate: FloatArray
) -> Dict[str, float]:
    """Compare two intensity vectors defined on the same bins.

    Args:
        reference: The reference intensities (SCiLS Lab).
        candidate: The intensities under test (Thyra), same length.

    Returns:
        ``rel_l1`` (``sum |c - r| / sum |r|``), ``pearson_r``,
        ``sum_ratio`` (``sum c / sum r``) and ``max_rel_bin`` (the largest
        per-bin difference relative to the reference maximum).
    """
    ref = np.asarray(reference, dtype=np.float64)
    cand = np.asarray(candidate, dtype=np.float64)
    if ref.shape != cand.shape:
        raise ValueError(
            f"intensity vectors differ in shape: {ref.shape} vs {cand.shape}"
        )
    ref_l1 = float(np.sum(np.abs(ref)))
    ref_max = float(np.max(np.abs(ref))) if ref.size else 0.0
    return {
        "rel_l1": (
            float(np.sum(np.abs(cand - ref)) / ref_l1) if ref_l1 else float("nan")
        ),
        "pearson_r": _pearson(ref, cand),
        "sum_ratio": float(np.sum(cand) / np.sum(ref)) if np.sum(ref) else float("nan"),
        "max_rel_bin": (
            float(np.max(np.abs(cand - ref)) / ref_max) if ref_max else float("nan")
        ),
    }


def coarse_edges(lo: float, hi: float, n: int) -> FloatArray:
    """Equal-width edges over ``[lo, hi]`` for :func:`rebin_sum`."""
    return np.linspace(lo, hi, n + 1).astype(np.float64)


def rebin_sum(mz: FloatArray, intensity: FloatArray, edges: FloatArray) -> FloatArray:
    """Sum intensities into coarse bins.

    Used to compare two spectra whose fine axes differ: summing into bins
    much wider than either axis's spacing keeps every count and makes the
    comparison insensitive to where exactly each fine bin sits. Values
    outside the edges are dropped.

    Args:
        mz: The fine axis.
        intensity: The intensities on it.
        edges: Increasing coarse bin edges.

    Returns:
        One summed intensity per coarse bin.
    """
    idx = np.searchsorted(edges, np.asarray(mz, dtype=np.float64), side="right") - 1
    keep = (idx >= 0) & (idx < edges.size - 1)
    return np.bincount(
        idx[keep],
        weights=np.asarray(intensity, dtype=np.float64)[keep],
        minlength=edges.size - 1,
    ).astype(np.float64)


def ratio_summary(reference: FloatArray, candidate: FloatArray) -> Dict[str, float]:
    """Summarise a per-pixel ratio ``candidate / reference``.

    Pixels where either side is not finite (a SCiLS normalisation factor of
    an empty spot is NaN) are left out of everything and counted. Pixels
    where the reference is zero are left out of the ratio and counted; so
    are pixels where the candidate is zero and the reference is not, which
    are the largest disagreement a ratio can show and would otherwise drop
    out of ``max_abs_log2_ratio`` unseen.

    Returns:
        ``pearson_r`` over the finite pixels, ``median_ratio``,
        ``p05_ratio``, ``p95_ratio``, ``max_abs_log2_ratio`` (over pixels
        where both sides are positive), ``zero_reference``,
        ``zero_candidate`` and ``non_finite``.
    """
    ref = np.asarray(reference, dtype=np.float64)
    cand = np.asarray(candidate, dtype=np.float64)
    finite = np.isfinite(ref) & np.isfinite(cand)
    ref, cand = ref[finite], cand[finite]
    nz = ref != 0
    ratio = cand[nz] / ref[nz]
    positive = ratio[ratio > 0]
    return {
        "pearson_r": _pearson(ref, cand),
        "median_ratio": float(np.median(ratio)) if ratio.size else float("nan"),
        "p05_ratio": float(np.percentile(ratio, 5)) if ratio.size else float("nan"),
        "p95_ratio": float(np.percentile(ratio, 95)) if ratio.size else float("nan"),
        "max_abs_log2_ratio": (
            float(np.max(np.abs(np.log2(positive)))) if positive.size else float("nan")
        ),
        "zero_reference": float(np.count_nonzero(~nz)),
        "zero_candidate": float(np.count_nonzero(ratio == 0)),
        "non_finite": float(np.count_nonzero(~finite)),
    }


def grid_index(values: FloatArray) -> NDArray[np.int64]:
    """Turn coordinates on a regular raster into integer steps from the minimum.

    The step is the most common gap between neighbouring distinct
    coordinates, which for a raster is its pitch: a gap left by a missing
    row is a multiple of it, and one irregular position is outvoted rather
    than taken as the pitch. A coordinate set with one distinct value maps
    to zeros.
    """
    v = np.asarray(values, dtype=np.float64)
    distinct = np.unique(np.round(v, 6))
    if distinct.size < 2:
        return np.zeros(v.size, dtype=np.int64)
    gaps = np.round(np.diff(distinct), 6)
    values_, counts = np.unique(gaps, return_counts=True)
    step = float(values_[np.argmax(counts)])
    return np.rint((v - distinct[0]) / step).astype(np.int64)


@dataclass(frozen=True)
class PixelMatch:
    """How reference spots were paired with candidate pixels.

    Attributes:
        reference_index: Positions into the reference spot arrays.
        candidate_index: The paired positions into the candidate pixels.
        orientation: Which axis flips made the pairing, e.g. ``"identity"``
            or ``"flip_y"``.
        shift: ``(dx, dy)`` in raster steps added to the reference after the
            flip. Non-zero when the two sides' bounding boxes differ, e.g.
            when one side lacks an empty edge row the other keeps.
        n_reference: Spots on the reference side.
        n_candidate: Pixels on the candidate side.
        tie_broken_by_tic: Whether several placements paired the same
            number of pixels and the per-pixel TIC correlation chose among
            them. When it did, a high TIC correlation is partly the result
            of the choice, which the report says.
    """

    reference_index: NDArray[np.int64]
    candidate_index: NDArray[np.int64]
    orientation: str
    shift: Tuple[int, int]
    n_reference: int
    n_candidate: int
    tie_broken_by_tic: bool


_ORIENTATIONS: Tuple[Tuple[str, bool, bool], ...] = (
    ("identity", False, False),
    ("flip_y", False, True),
    ("flip_x", True, False),
    ("flip_xy", True, True),
)

#: Largest shift, in raster steps, tried between the two sides.
MAX_SHIFT = 3


def match_pixels(
    ref_x: FloatArray,
    ref_y: FloatArray,
    cand_x: FloatArray,
    cand_y: FloatArray,
    ref_tic: Optional[FloatArray] = None,
    cand_tic: Optional[FloatArray] = None,
) -> PixelMatch:
    """Pair reference spots with candidate pixels by raster position.

    Both sides are reduced to integer raster steps (:func:`grid_index`),
    each anchored at its own minimum, and joined on them. The reference may
    count an axis the other way, and the two bounding boxes may differ by a
    few empty edge rows, so every flip and every shift up to
    :data:`MAX_SHIFT` steps is tried and the placement pairing the most
    pixels kept. Ties -- a full rectangle pairs equally under every flip --
    are broken by the per-pixel TIC correlation when both TIC vectors are
    given, else by the smallest shift and the order of the flips.

    Args:
        ref_x, ref_y: Reference spot coordinates, any unit.
        cand_x, cand_y: Candidate pixel coordinates, any unit.
        ref_tic, cand_tic: Optional per-spot and per-pixel TIC.

    Returns:
        The pairing; see :class:`PixelMatch`.
    """
    rx, ry = grid_index(ref_x), grid_index(ref_y)
    cx, cy = grid_index(cand_x), grid_index(cand_y)
    pad = MAX_SHIFT
    width = int(max(rx.max(), cx.max())) + 2 * pad + 1
    cand_keys = (cy + pad) * width + (cx + pad)
    order = np.argsort(cand_keys)
    sorted_keys = cand_keys[order]

    placements = []
    shifts = [(dx, dy) for dx in range(-pad, pad + 1) for dy in range(-pad, pad + 1)]
    shifts.sort(key=lambda d: (abs(d[0]) + abs(d[1]), d))
    for name, flip_x, flip_y in _ORIENTATIONS:
        fx = rx.max() - rx if flip_x else rx
        fy = ry.max() - ry if flip_y else ry
        for dx, dy in shifts:
            keys = (fy + dy + pad) * width + (fx + dx + pad)
            pos = np.clip(np.searchsorted(sorted_keys, keys), 0, sorted_keys.size - 1)
            hit = sorted_keys[pos] == keys
            placements.append((name, (dx, dy), np.flatnonzero(hit), order[pos[hit]]))

    best = max(p[2].size for p in placements)
    tied = [p for p in placements if p[2].size == best]
    tie_broken = False
    chosen = tied[0]
    if len(tied) > 1 and ref_tic is not None and cand_tic is not None and best > 1:
        tie_broken = True
        ref_t, cand_t = np.asarray(ref_tic), np.asarray(cand_tic)

        def _score(p: Tuple[str, Tuple[int, int], NDArray, NDArray]) -> float:
            r = _pearson(ref_t[p[2]], cand_t[p[3]])
            return -np.inf if np.isnan(r) else r

        chosen = max(tied, key=_score)
    name, shift, ref_index, cand_index = chosen
    return PixelMatch(
        reference_index=np.asarray(ref_index, dtype=np.int64),
        candidate_index=np.asarray(cand_index, dtype=np.int64),
        orientation=name,
        shift=(int(shift[0]), int(shift[1])),
        n_reference=int(rx.size),
        n_candidate=int(cx.size),
        tie_broken_by_tic=tie_broken,
    )
