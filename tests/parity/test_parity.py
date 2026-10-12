"""Thyra against SCiLS Lab, one dataset at a time.

Each test asserts one group of the checks :func:`tests.parity.report.verdicts`
lists, so the failure message carries the same measured numbers the report
prints. The thresholds are in :mod:`tests.parity.thresholds`.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import pytest

from tests.parity.report import verdicts

Verdict = Tuple[str, str, bool]


def _checks(result: Dict[str, Any], prefix: str) -> List[Verdict]:
    chosen = [v for v in verdicts(result) if v[0].startswith(prefix)]
    if not chosen:
        pytest.skip(f"no '{prefix}' checks apply to {result['dataset']}")
    return chosen


def _assert_all(result: Dict[str, Any], checks: List[Verdict]) -> None:
    failed = [f"{name}: measured {measured}" for name, measured, ok in checks if not ok]
    assert not failed, f"{result['dataset']}: " + "; ".join(failed)


def test_axis(parity_result: Dict[str, Any]) -> None:
    """Bin count, offset in ppm and in bins, and spacing law."""
    _assert_all(parity_result, _checks(parity_result, "axis"))


def test_mean_spectrum(parity_result: Dict[str, Any]) -> None:
    """Mean spectrum bin by bin, where the bins correspond."""
    _assert_all(parity_result, _checks(parity_result, "mean"))


def test_mean_spectrum_coarse(parity_result: Dict[str, Any]) -> None:
    """Mean spectrum summed into coarse bins, independent of the fine axis."""
    _assert_all(parity_result, _checks(parity_result, "coarse mean"))


def test_pixels_pair(parity_result: Dict[str, Any]) -> None:
    """Enough SCiLS spots find a Thyra pixel at the same raster position."""
    _assert_all(parity_result, _checks(parity_result, "matched spots"))


def test_tic(parity_result: Dict[str, Any]) -> None:
    """Per-pixel TIC: correlation and median ratio."""
    _assert_all(parity_result, _checks(parity_result, "TIC"))


def test_spot_spectra(parity_result: Dict[str, Any]) -> None:
    """A handful of individual spectra."""
    _assert_all(parity_result, _checks(parity_result, "spot"))
