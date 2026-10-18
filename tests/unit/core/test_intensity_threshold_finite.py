"""An intensity threshold must be a finite number (#503, D62).

No peak's sum reaches NaN or +infinity, so either one emptied the store;
-infinity kept everything. Each was recorded in the processing step, where
JSON cannot hold it. The CLI refused them already; the reader now refuses
them on the API too.
"""

from pathlib import Path

import numpy as np
import pytest

from tests.fixtures.mock_msi_generator import MockMSIReader
from thyra.core.base_reader import BaseMSIReader
from thyra.errors import ConversionRefused


def _init(threshold):
    reader = MockMSIReader.__new__(MockMSIReader)
    BaseMSIReader.__init__(reader, Path("mock"), intensity_threshold=threshold)
    return reader


@pytest.mark.parametrize(
    "threshold",
    [float("nan"), float("inf"), float("-inf"), np.float32("nan")],
    ids=["nan", "inf", "-inf", "np32-nan"],
)
def test_a_non_finite_threshold_is_refused(threshold):
    with pytest.raises(ConversionRefused, match="intensity_threshold"):
        _init(threshold)


@pytest.mark.parametrize("threshold", [None, 0, 50, 12.5, np.int64(3)])
def test_a_finite_threshold_or_none_is_taken(threshold):
    assert _init(threshold).intensity_threshold == threshold
