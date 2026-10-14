"""The one rule that turns mobility values into a range to bin over (D46).

The heatmap and the grid both bin over the range a reader states, and
every reader states it through :func:`thyra.core.mobility.value_span`.
"""

import numpy as np
import pytest

from thyra.core.mobility import value_span


class TestValueSpan:
    def test_the_span_of_ordinary_values(self):
        assert value_span(np.array([1.1, 0.8, 0.95])) == (0.8, 1.1)

    def test_non_finite_values_are_not_positions(self):
        values = np.array([np.nan, 0.9, np.inf, 1.2, -np.inf])
        assert value_span(values) == (0.9, 1.2)

    @pytest.mark.parametrize(
        "values",
        [
            np.array([]),
            np.array([np.nan, np.nan]),
            np.array([np.inf, -np.inf]),
        ],
        ids=["empty", "all-nan", "all-infinite"],
    )
    def test_no_finite_value_gives_none(self, values):
        assert value_span(values) is None

    def test_a_single_value_spans_nothing(self):
        assert value_span(np.array([1.0, 1.0, np.nan])) is None

    def test_a_list_is_taken_as_well(self):
        assert value_span([0.7, 1.3]) == (0.7, 1.3)

    def test_the_edges_are_plain_floats(self):
        lower, upper = value_span(np.array([0.8, 1.1], dtype=np.float32))
        assert type(lower) is float and type(upper) is float
