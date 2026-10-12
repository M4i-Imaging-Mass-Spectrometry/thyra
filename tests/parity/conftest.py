"""Fixtures for the parity lane: one evaluation per dataset per session."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import pytest

from tests.parity.datasets import PARITY_ENV, ParityDataset, discover

#: Set to 1 to convert again even when a store for the same settings exists.
RECONVERT_ENV = "THYRA_PARITY_RECONVERT"


def _datasets() -> List[ParityDataset]:
    return [d for d in discover() if d.has_reference]


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Parametrize ``parity_dataset`` over every dataset with a reference.

    With no data the parameter list holds one placeholder that skips, so a
    run without ``THYRA_PARITY_DATA`` reports the lane as skipped rather
    than collecting nothing.
    """
    if "parity_dataset" not in metafunc.fixturenames:
        return
    datasets = _datasets()
    if datasets:
        metafunc.parametrize("parity_dataset", datasets, ids=[d.name for d in datasets])
    else:
        metafunc.parametrize(
            "parity_dataset",
            [
                pytest.param(
                    None,
                    marks=pytest.mark.skip(
                        reason=f"no parity data: set {PARITY_ENV} to a folder of datasets"
                    ),
                )
            ],
            ids=["no-data"],
        )


# A result, or the exception evaluating it raised: a dataset that fails is
# failed once, not converted and read again by each of its tests.
_CACHE: Dict[str, Any] = {}


@pytest.fixture
def parity_result(parity_dataset: ParityDataset) -> Dict[str, Any]:
    """The evaluated result for one dataset, computed once per session."""
    from tests.parity.report import evaluate

    if not parity_dataset.thyra_supported:
        pytest.skip(f"{parity_dataset.name}: Thyra has no reader for this source")
    if parity_dataset.name not in _CACHE:
        reconvert = os.environ.get(RECONVERT_ENV, "") == "1"
        try:
            _CACHE[parity_dataset.name] = evaluate(parity_dataset, reconvert=reconvert)
        except Exception as exc:  # cached and re-raised for every test
            _CACHE[parity_dataset.name] = exc
    cached = _CACHE[parity_dataset.name]
    if isinstance(cached, Exception):
        raise cached
    return cached
