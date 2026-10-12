"""Compute every parity metric for a dataset and write it up as markdown.

:func:`evaluate` is what both the pytest tests and the command line use,
so a test failure and a report line always come from the same number.

Command line::

    python -m tests.parity.report [--root DIR] [--out DIR] [--reconvert]
                                  [dataset ...]

``--root`` defaults to ``$THYRA_PARITY_DATA``; ``--out`` to
``<root>/_reports``. One ``<dataset>.md`` and one ``<dataset>.json`` are
written per dataset, plus ``summary.md`` with one row each.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from tests.parity import thresholds as T
from tests.parity.datasets import PARITY_ENV, ParityDataset, discover, parity_root
from tests.parity.metrics import (
    axis_agreement,
    coarse_edges,
    intensity_agreement,
    match_pixels,
    nearest_law,
    ratio_summary,
    rebin_sum,
    spacing_exponent,
)
from tests.parity.scils_reference import ScilsReference, load_reference
from tests.parity.thyra_runner import (
    ThyraOutput,
    convert,
    import_axis_type,
    matched_resampling_config,
    read_output,
)

Result = Dict[str, Any]


def _coarse(
    ref_mz: np.ndarray, ref_i: np.ndarray, cand_mz: np.ndarray, cand_i: np.ndarray
) -> Dict[str, float]:
    """Intensity agreement after both sides are summed into coarse bins.

    An empty side (a SCiLS spot with no stored values) has no range to
    share; it is reported as not measurable rather than raising.
    """
    if ref_mz.size == 0 or cand_mz.size == 0:
        nan = float("nan")
        return {"rel_l1": nan, "pearson_r": nan, "sum_ratio": nan, "max_rel_bin": nan}
    lo = max(float(ref_mz[0]), float(cand_mz[0]))
    hi = min(float(ref_mz[-1]), float(cand_mz[-1]))
    edges = coarse_edges(lo, hi, T.COARSE_BINS)
    return intensity_agreement(
        rebin_sum(ref_mz, ref_i, edges), rebin_sum(cand_mz, cand_i, edges)
    )


def scils_only(reference: ScilsReference) -> Result:
    """The SCiLS-side facts a report states even when nothing is compared."""
    exponent = spacing_exponent(reference.mz)
    tic = reference.normalizations.get("tic")
    tic_check = None
    if tic is not None and reference.spot_sum.mean() > 0:
        expected = reference.spot_sum / reference.spot_sum.mean()
        finite = np.isfinite(tic)
        tic_check = float(np.max(np.abs(tic[finite] / expected[finite] - 1.0)))
    return {
        "bins": int(reference.mz.size),
        "min_mz": float(reference.mz[0]),
        "max_mz": float(reference.mz[-1]),
        "spacing_exponent": exponent,
        "law": nearest_law(exponent),
        "import_axis_type": reference.import_axis_type or None,
        "import_rebinning": reference.import_rebinning or None,
        "spots": int(reference.spot_id.size),
        "normalizations": sorted(reference.normalizations),
        "tic_factor_vs_sum_over_mean": tic_check,
        "server_version": reference.info.get("server_version"),
    }


def compare(reference: ScilsReference, output: ThyraOutput) -> Result:
    """Every comparison between a SCiLS reference and a Thyra store."""
    result: Result = {"axis": axis_agreement(reference.mz, output.mz)}
    same_bins = reference.mz.size == output.mz.size

    result["mean"] = {
        "binwise": (
            intensity_agreement(reference.mean, output.mean) if same_bins else None
        ),
        "coarse": _coarse(reference.mz, reference.mean, output.mz, output.mean),
    }

    match = match_pixels(
        reference.spot_x,
        reference.spot_y,
        output.pixel_x,
        output.pixel_y,
        reference.spot_sum,
        output.tic,
    )
    ri, ci = match.reference_index, match.candidate_index
    result["pixels"] = {
        "scils_spots": match.n_reference,
        "thyra_pixels": match.n_candidate,
        "matched": int(ri.size),
        "matched_fraction": float(ri.size / max(match.n_reference, 1)),
        "orientation": match.orientation,
        "shift_steps": list(match.shift),
        "orientation_chosen_by_tic": match.tie_broken_by_tic,
    }
    result["tic"] = result["tic_factor"] = result["rms"] = None
    if ri.size:
        result["tic"] = ratio_summary(reference.spot_sum[ri], output.tic[ci])
        ref_tic = reference.spot_sum[ri]
        result["tic"]["reference_cv"] = (
            float(np.std(ref_tic) / np.mean(ref_tic))
            if np.mean(ref_tic)
            else float("nan")
        )
        # SCiLS's TIC factor is the spot's sum over the mean sum; the same
        # quantity on Thyra's side compares the per-pixel shape of the TIC
        # with any overall scale between the two removed.
        tic_factor = reference.normalizations.get("tic")
        thyra_tic = output.tic[ci]
        if tic_factor is not None and thyra_tic.mean() > 0:
            result["tic_factor"] = ratio_summary(
                tic_factor[ri], thyra_tic / thyra_tic.mean()
            )
        # SCiLS's RMS factor is a dataset constant times the RMS of the values
        # the spot stores (for centroid data the peaks alone, for profile data
        # every bin), so only the correlation is comparable, not the ratio.
        rms = reference.normalizations.get("rms")
        if rms is not None:
            thyra_rms = np.sqrt(output.sum_squares[ci] / max(output.mz.size, 1))
            result["rms"] = ratio_summary(rms[ri], thyra_rms)

    pair = dict(zip(ri.tolist(), ci.tolist()))
    position = {int(s): i for i, s in enumerate(reference.spot_id)}
    spots: List[Result] = []
    for sample in reference.samples:
        row = pair.get(position.get(sample.spot_id, -1))
        if row is None or row not in output.rows:
            spots.append({"spot_id": sample.spot_id, "matched": False})
            continue
        spectrum = output.rows[row]
        spots.append(
            {
                "spot_id": sample.spot_id,
                "matched": True,
                "thyra_row": int(row),
                "binwise": (
                    intensity_agreement(sample.rebinned, spectrum)
                    if same_bins
                    else None
                ),
                "coarse": _coarse(reference.mz, sample.rebinned, output.mz, spectrum),
                "coarse_vs_raw": _coarse(
                    sample.raw_mz, sample.raw_intensity, output.mz, spectrum
                ),
            }
        )
    result["spots"] = spots
    return result


def _rows_for_samples(reference: ScilsReference, output: ThyraOutput) -> List[int]:
    """The Thyra rows paired with the SCiLS sample spots."""
    match = match_pixels(
        reference.spot_x,
        reference.spot_y,
        output.pixel_x,
        output.pixel_y,
        reference.spot_sum,
        output.tic,
    )
    pair = dict(zip(match.reference_index.tolist(), match.candidate_index.tolist()))
    position = {int(s): i for i, s in enumerate(reference.spot_id)}
    rows = [pair.get(position.get(s.spot_id, -1)) for s in reference.samples]
    return [int(r) for r in rows if r is not None]


def _convert_and_compare(
    dataset: ParityDataset, reference: ScilsReference, reconvert: bool
) -> Result:
    """Convert one dataset (or variant) and compare it with the reference."""
    config = matched_resampling_config(dataset, reference)
    started = time.monotonic()
    store = convert(dataset, config, reconvert=reconvert)
    seconds = round(time.monotonic() - started, 1)
    # Two passes: the first finds which rows the sample spots pair with,
    # the second reads those rows. Row sums come with both.
    first = read_output(store)
    output = read_output(store, rows=_rows_for_samples(reference, first))
    return {
        "settings": config,
        "reader_options": dataset.reader_options,
        "store": str(store),
        "convert_seconds": seconds,
        "comparison": compare(reference, output),
    }


def evaluate(dataset: ParityDataset, reconvert: bool = False) -> Result:
    """Convert (or reuse) and compare one dataset and its variants.

    Args:
        dataset: The dataset.
        reconvert: Convert even when a store for these settings exists.

    Returns:
        A JSON-serialisable result: ``scils``, and when Thyra reads the
        source, ``settings``, ``store``, ``comparison`` and ``variants``
        (each of the manifest's variants, with the same keys).
    """
    reference = load_reference(dataset.scils_dir)
    result: Result = {
        "dataset": dataset.name,
        "manifest": dataset.manifest,
        "scils": scils_only(reference),
    }
    if not dataset.thyra_supported:
        result["skipped"] = "Thyra has no reader for this source"
        return result
    result["scils"]["import_axis_type_as_thyra"] = import_axis_type(reference)
    result.update(_convert_and_compare(dataset, reference, reconvert))
    result["variants"] = {
        label: _convert_and_compare(variant, reference, reconvert)
        for label, variant in dataset.variants().items()
    }
    return result


# --- verdicts ---------------------------------------------------------------


def _ok(value: Optional[float]) -> bool:
    """A measured value usable in a comparison: present and not NaN.

    The JSON round trip of a report turns NaN into None, so both mean "not
    measurable", which fails a check rather than raising.
    """
    return value is not None and not math.isnan(value)


def _le(value: Optional[float], limit: float) -> bool:
    return _ok(value) and value <= limit  # type: ignore[operator]


def _ge(value: Optional[float], limit: float) -> bool:
    return _ok(value) and value >= limit  # type: ignore[operator]


def _within(value: Optional[float], lo: float, hi: float) -> bool:
    return value is not None and not math.isnan(value) and lo <= value <= hi


def verdicts(result: Result) -> List[Tuple[str, str, bool]]:
    """``(check, measured, passed)`` for every threshold in :mod:`thresholds`.

    The same list drives the pytest assertions and the report's verdict
    column, so the two cannot disagree.
    """
    comparison = result.get("comparison")
    if not comparison:
        return []
    out: List[Tuple[str, str, bool]] = []
    axis = comparison["axis"]
    out.append(
        (
            "axis bin count",
            f"{int(axis['bins_candidate'])} vs {int(axis['bins_reference'])}",
            bool(axis["bins_equal"]),
        )
    )
    out.append(
        (
            f"axis max ppm <= {T.AXIS_MAX_PPM}",
            _fmt(axis["max_ppm"]),
            _le(axis["max_ppm"], T.AXIS_MAX_PPM),
        )
    )
    out.append(
        (
            f"axis max offset <= {T.AXIS_MAX_BIN_FRACTION} bin",
            _fmt(axis["max_ppm_in_bins"]),
            _le(axis["max_ppm_in_bins"], T.AXIS_MAX_BIN_FRACTION),
        )
    )
    exp_diff = abs(axis["exponent_candidate"] - axis["exponent_reference"])
    out.append(
        (
            f"axis law exponent diff <= {T.AXIS_MAX_EXPONENT_DIFF}",
            _fmt(exp_diff),
            _le(exp_diff, T.AXIS_MAX_EXPONENT_DIFF),
        )
    )
    binwise = comparison["mean"]["binwise"]
    if binwise is not None:
        out.append(
            (
                f"mean rel L1 <= {T.MEAN_MAX_REL_L1}",
                _fmt(binwise["rel_l1"]),
                _le(binwise["rel_l1"], T.MEAN_MAX_REL_L1),
            )
        )
        out.append(
            (
                f"mean Pearson r >= {T.MEAN_MIN_PEARSON}",
                _fmt(binwise["pearson_r"]),
                _ge(binwise["pearson_r"], T.MEAN_MIN_PEARSON),
            )
        )
        out.append(
            (
                f"mean sum ratio in {T.MEAN_SUM_RATIO_RANGE}",
                _fmt(binwise["sum_ratio"]),
                _within(binwise["sum_ratio"], *T.MEAN_SUM_RATIO_RANGE),
            )
        )
    coarse = comparison["mean"]["coarse"]
    out.append(
        (
            f"coarse mean rel L1 <= {T.COARSE_MEAN_MAX_REL_L1}",
            _fmt(coarse["rel_l1"]),
            _le(coarse["rel_l1"], T.COARSE_MEAN_MAX_REL_L1),
        )
    )
    out.append(
        (
            f"coarse mean Pearson r >= {T.COARSE_MEAN_MIN_PEARSON}",
            _fmt(coarse["pearson_r"]),
            _ge(coarse["pearson_r"], T.COARSE_MEAN_MIN_PEARSON),
        )
    )
    pixels = comparison["pixels"]
    out.append(
        (
            f"matched spots >= {T.MIN_MATCHED_FRACTION:.0%}",
            f"{pixels['matched']}/{pixels['scils_spots']}",
            _ge(pixels["matched_fraction"], T.MIN_MATCHED_FRACTION),
        )
    )
    tic = comparison["tic"]
    if tic is not None:
        cv = tic.get("reference_cv")
        if _ok(cv) and cv < T.TIC_MIN_REFERENCE_CV:
            out.append(
                (
                    f"TIC Pearson r (not informative: SCiLS TIC cv {_fmt(cv)} "
                    f"< {T.TIC_MIN_REFERENCE_CV})",
                    _fmt(tic["pearson_r"]),
                    True,
                )
            )
        else:
            out.append(
                (
                    f"TIC Pearson r >= {T.TIC_MIN_PEARSON}",
                    _fmt(tic["pearson_r"]),
                    _ge(tic["pearson_r"], T.TIC_MIN_PEARSON),
                )
            )
        out.append(
            (
                f"TIC median ratio in {T.TIC_MEDIAN_RATIO_RANGE}",
                _fmt(tic["median_ratio"]),
                _within(tic["median_ratio"], *T.TIC_MEDIAN_RATIO_RANGE),
            )
        )
    for spot in comparison["spots"]:
        if not spot.get("matched"):
            out.append((f"spot {spot['spot_id']} paired", "no", False))
            continue
        if spot["binwise"] is not None:
            r = spot["binwise"]["pearson_r"]
            out.append(
                (
                    f"spot {spot['spot_id']} Pearson r >= {T.SPOT_MIN_PEARSON}",
                    _fmt(r),
                    _ge(r, T.SPOT_MIN_PEARSON),
                )
            )
        r = spot["coarse"]["pearson_r"]
        out.append(
            (
                f"spot {spot['spot_id']} coarse r >= {T.SPOT_COARSE_MIN_PEARSON}",
                _fmt(r),
                _ge(r, T.SPOT_COARSE_MIN_PEARSON),
            )
        )
    return out


# --- markdown ---------------------------------------------------------------


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return "nan" if math.isnan(value) else f"{value:.4g}"
    return str(value)


def _table(rows: Sequence[Tuple[str, Any]]) -> List[str]:
    lines = ["| quantity | value |", "|---|---|"]
    lines += [f"| {k} | {_fmt(v)} |" for k, v in rows]
    return lines


def render(result: Result) -> str:
    """One dataset's report as markdown."""
    manifest = result["manifest"]
    scils = result["scils"]
    lines = [f"# Parity: {result['dataset']}", ""]
    if manifest.get("notes"):
        lines += [manifest["notes"], ""]
    lines += [f"SCiLS import: {manifest.get('scils_import', '-')}", ""]
    lines += ["## SCiLS Lab axis", ""]
    lines += _table(
        [
            ("bins", scils["bins"]),
            ("m/z range", f"{scils['min_mz']:.5f} - {scils['max_mz']:.5f}"),
            ("fitted spacing exponent", scils["spacing_exponent"]),
            ("law (nearest)", scils["law"]),
            ("axis type in saved import", scils["import_axis_type"]),
            ("rebinning in saved import", scils["import_rebinning"]),
            ("spots", scils["spots"]),
            ("normalizations", ", ".join(scils["normalizations"])),
            (
                "TIC factor vs sum / mean sum, max rel",
                scils["tic_factor_vs_sum_over_mean"],
            ),
            ("server", scils["server_version"]),
        ]
    )
    lines.append("")
    if result.get("skipped"):
        lines += [f"Comparison skipped: {result['skipped']}.", ""]
        return "\n".join(lines)

    settings = result["settings"]
    lines += ["## Thyra settings", ""]
    lines += _table(sorted(settings.items()))
    lines += [
        "",
        f"Reader options: `{json.dumps(manifest.get('reader_options') or {})}`. "
        f"Conversion (or reuse): {result['convert_seconds']} s.",
        "",
    ]
    comparison = result["comparison"]
    axis = comparison["axis"]
    lines += ["## Axis", ""]
    lines += _table(sorted(axis.items()))
    lines.append("")
    lines += ["## Mean spectrum", ""]
    for kind in ("binwise", "coarse"):
        block = comparison["mean"][kind]
        label = kind if kind == "binwise" else f"coarse ({T.COARSE_BINS} bins)"
        if block is None:
            lines += [f"{label}: not computed (bin counts differ)", ""]
        else:
            lines += [f"{label}:", ""] + _table(sorted(block.items())) + [""]
    lines += ["## Pixels and TIC", ""]
    lines += _table(sorted(comparison["pixels"].items()))
    lines.append("")
    titles = {
        "tic": "TIC, Thyra row sum / SCiLS spot sum",
        "tic_factor": "TIC factor, Thyra (TIC / mean TIC) / SCiLS TIC normalisation factor",
        "rms": (
            "RMS, Thyra RMS over all bins / SCiLS RMS factor. SCiLS's factor is a "
            "dataset constant times the RMS of the stored values (centroid data: the "
            "peaks alone), so read the correlation, not the ratio"
        ),
    }
    for kind, title in titles.items():
        block = comparison.get(kind)
        if block is not None:
            lines += [f"{title}:", ""] + _table(sorted(block.items())) + [""]
    lines += ["## Spot spectra", ""]
    lines += [
        "| spot | Thyra row | binwise r | binwise rel L1 | sum ratio | coarse r | "
        "coarse r vs SCiLS raw |",
        "|---|---|---|---|---|---|---|",
    ]
    for spot in comparison["spots"]:
        if not spot.get("matched"):
            lines.append(f"| {spot['spot_id']} | not paired | | | | | |")
            continue
        b = spot["binwise"] or {}
        lines.append(
            f"| {spot['spot_id']} | {spot['thyra_row']} | {_fmt(b.get('pearson_r'))} | "
            f"{_fmt(b.get('rel_l1'))} | {_fmt(b.get('sum_ratio'))} | "
            f"{_fmt(spot['coarse']['pearson_r'])} | "
            f"{_fmt(spot['coarse_vs_raw']['pearson_r'])} |"
        )
    lines += ["", "## Verdicts", "", "| check | measured | pass |", "|---|---|---|"]
    for name, measured, passed in verdicts(result):
        lines.append(f"| {name} | {measured} | {'yes' if passed else 'NO'} |")
    lines.append("")
    variants = result.get("variants") or {}
    if variants:
        lines += ["## Variants", "", "Reported, not asserted on.", ""]
        lines += _SUMMARY_HEADER
        for label, block in variants.items():
            lines.append(_summary_row(f"{result['dataset']} [{label}]", scils, block))
        lines.append("")
    return "\n".join(lines)


_SUMMARY_HEADER = [
    "| dataset | reader options | method | SCiLS law / bins | Thyra law / bins "
    "| axis max ppm | max offset (bins) | mean r | mean rel L1 | sum ratio | coarse r "
    "| TIC r | TIC median ratio | TIC factor r | matched | checks failed |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
]


def _summary_row(label: str, scils: Result, block: Result) -> str:
    """One summary line for a dataset's main entry or one of its variants."""
    c = block["comparison"]
    axis, mean = c["axis"], c["mean"]
    tic, factor = c.get("tic") or {}, c.get("tic_factor") or {}
    b = mean["binwise"] or {}
    failed = [name for name, _, ok in verdicts(block) if not ok]
    options = (
        json.dumps(block.get("reader_options") or {})
        if block.get("reader_options")
        else "-"
    )
    return (
        f"| {label} | {options} | {block['settings']['method']} | "
        f"{scils['law']} / {scils['bins']} | "
        f"{nearest_law(axis['exponent_candidate'])} / {int(axis['bins_candidate'])} | "
        f"{_fmt(axis['max_ppm'])} | {_fmt(axis['max_ppm_in_bins'])} | "
        f"{_fmt(b.get('pearson_r'))} | {_fmt(b.get('rel_l1'))} | "
        f"{_fmt(b['sum_ratio']) if b else _fmt(mean['coarse']['sum_ratio']) + ' (coarse)'} | "
        f"{_fmt(mean['coarse']['pearson_r'])} | "
        f"{_fmt(tic.get('pearson_r'))} | {_fmt(tic.get('median_ratio'))} | "
        f"{_fmt(factor.get('pearson_r'))} | "
        f"{c['pixels']['matched']}/{c['pixels']['scils_spots']} | {len(failed)} |"
    )


def render_summary(results: Sequence[Result]) -> str:
    """One row per dataset, and one per variant."""
    lines = ["# Parity summary", ""] + list(_SUMMARY_HEADER)
    for result in results:
        scils = result["scils"]
        if result.get("skipped"):
            lines.append(
                f"| {result['dataset']} | - | - | {scils['law']} / {scils['bins']} | "
                f"skipped: {result['skipped']} | | | | | | | | | | | |"
            )
            continue
        lines.append(_summary_row(result["dataset"], scils, result))
        for label, block in (result.get("variants") or {}).items():
            lines.append(_summary_row(f"{result['dataset']} [{label}]", scils, block))
    lines.append("")
    return "\n".join(lines)


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return None if math.isnan(float(value)) else float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point: write one report per dataset and a summary."""
    parser = argparse.ArgumentParser(
        description="Write Thyra vs SCiLS Lab parity reports."
    )
    parser.add_argument("datasets", nargs="*", help="dataset names (default: all)")
    parser.add_argument(
        "--root", type=Path, default=None, help=f"default: ${PARITY_ENV}"
    )
    parser.add_argument(
        "--out", type=Path, default=None, help="default: <root>/_reports"
    )
    parser.add_argument("--reconvert", action="store_true")
    args = parser.parse_args(argv)

    root = args.root or parity_root()
    if root is None:
        parser.error(f"no data: pass --root or set {PARITY_ENV}")
    out = args.out or root / "_reports"
    out.mkdir(parents=True, exist_ok=True)
    chosen = [d for d in discover(root) if not args.datasets or d.name in args.datasets]
    results = []
    for dataset in chosen:
        if not dataset.has_reference:
            print(
                f"{dataset.name}: no reference.npz, run scils_extract.py first",
                flush=True,
            )
            continue
        print(f"{dataset.name}: evaluating", flush=True)
        try:
            result = _jsonable(evaluate(dataset, reconvert=args.reconvert))
        except Exception as exc:  # one dataset's failure must not end the run
            result = {
                "dataset": dataset.name,
                "manifest": dataset.manifest,
                "scils": _jsonable(scils_only(load_reference(dataset.scils_dir))),
                "skipped": f"error: {type(exc).__name__}: {exc}",
            }
        (out / f"{dataset.name}.json").write_text(
            json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
        )
        (out / f"{dataset.name}.md").write_text(render(result), encoding="utf-8")
        results.append(result)
    if not args.datasets:
        (out / "summary.md").write_text(render_summary(results), encoding="utf-8")
    print(render_summary(results), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
