# tests/unit/converters/test_nearest_neighbor_negatives.py
"""``nearest_neighbor`` must not silently discard negative bins.

``_nearest_neighbor_resample`` accumulates each source peak into its nearest
target bin and then keeps the bins that hold something. That filter used to
be ``accumulated > 0``, which drops a bin whose accumulated value came out
negative -- and dropping it raises the stored total ion current above the
input's, which is the opposite of what the method guarantees.

Baseline-subtracted spectra can carry negative intensities. Every export
examined so far filters them upstream, so this has no observed trigger in
practice, but the invariant is cheap to hold and the failure would be silent.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from thyra.converters.spatialdata.base_spatialdata_converter import (
    BaseSpatialDataConverter,
)


def _resample(axis, mzs, intensities):
    stub = SimpleNamespace(_common_mass_axis=axis)
    return BaseSpatialDataConverter._nearest_neighbor_resample(
        stub, np.asarray(mzs, float), np.asarray(intensities, float)
    )


class TestNegativeIntensities:
    """Negative bins are measurements, not absences."""

    def test_total_is_preserved_with_negatives(self):
        axis = np.linspace(250.0, 1200.0, 5_000)
        mzs = [300.0, 300.02, 500.0, 700.0]
        intensities = [-5.0, 3.0, 10.0, -2.0]

        _, values = _resample(axis, mzs, intensities)

        # Dropping the negative bins gave 13.0 against an input total of 6.0.
        assert values.sum() == pytest.approx(sum(intensities))

    def test_a_wholly_negative_bin_survives(self):
        axis = np.linspace(250.0, 1200.0, 5_000)
        indices, values = _resample(axis, [700.0], [-4.0])

        assert indices.size == 1
        assert values[0] == pytest.approx(-4.0)

    def test_bins_cancelling_to_zero_are_dropped(self):
        """Exactly zero still means nothing to store."""
        axis = np.linspace(250.0, 1200.0, 5_000)
        # Both land in the same bin and cancel.
        indices, values = _resample(axis, [700.0, 700.0001], [5.0, -5.0])

        assert indices.size == 0
        assert values.size == 0

    def test_all_positive_is_unchanged(self):
        """The ordinary case must be untouched by the fix."""
        axis = np.linspace(250.0, 1200.0, 5_000)
        mzs = [300.0, 500.0, 700.0]
        intensities = [1.0, 2.0, 3.0]

        indices, values = _resample(axis, mzs, intensities)

        assert indices.size == 3
        assert values.sum() == pytest.approx(6.0)
        assert (values > 0).all()
