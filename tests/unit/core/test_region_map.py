"""The shared region helpers: the per-region summary and the request check."""

import numpy as np
import pytest

from thyra.core.region_map import checked_region_number, summarize_regions


def test_summary_counts_a_position_given_twice_twice():
    info = summarize_regions([(0, 0, 1), (0, 0, 1), (2, 1, 1)])

    assert info == [{"region_number": 1, "n_spectra": 3, "bounds": (0, 0, 2, 1)}]


def test_summary_is_sorted_by_count_then_number():
    scans = [(0, 0, 5), (1, 0, 5), (4, 4, 2), (5, 4, 2), (9, 9, 7)]

    info = summarize_regions(scans)

    assert [(e["region_number"], e["n_spectra"]) for e in info] == [
        (2, 2),
        (5, 2),
        (7, 1),
    ]


def test_summary_names_only_regions_with_a_name():
    info = summarize_regions([(0, 0, 0), (1, 0, 1)], {0: "left", 1: ""})

    assert info[0]["name"] == "left"
    assert "name" not in info[1]


def test_summary_of_nothing_is_empty():
    assert summarize_regions([]) == []


@pytest.mark.parametrize(
    "value, expected",
    [(None, None), (3, 3), (np.int64(4), 4), ("2", 2), (" 5 ", 5)],
)
def test_checked_region_number_accepts_integers_and_decimal_strings(value, expected):
    result = checked_region_number(value)

    assert result == expected
    assert result is None or type(result) is int


@pytest.mark.parametrize("value", [True, 1.0, "a", "1.5", "", "-1", [1]])
def test_checked_region_number_refuses_the_rest(value):
    with pytest.raises(ValueError):
        checked_region_number(value)
