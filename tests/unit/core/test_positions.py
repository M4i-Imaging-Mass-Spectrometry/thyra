"""The one door for position, index and region columns (issue #527)."""

import re

import numpy as np
import pytest

from thyra.core.positions import whole_numbers
from thyra.errors import ConversionRefused


class TestAccepted:
    def test_integers_pass_as_int64(self):
        out = whole_numbers(np.array([1, 2, 3], dtype=np.int32), "x positions")
        assert out.dtype == np.int64
        assert out.tolist() == [1, 2, 3]

    def test_a_float_a_hair_from_a_whole_number_is_that_number(self):
        out = whole_numbers([1.0, 0.3 / 0.1, 100.9999999], "x positions")
        assert out.tolist() == [1, 3, 101]

    def test_uint64_below_2_63_passes(self):
        values = np.array([0, 2**63 - 1], dtype=np.uint64)
        assert whole_numbers(values, "indices").tolist() == [0, 2**63 - 1]

    def test_sqlite_objects_pass_when_whole(self):
        out = whole_numbers(np.array([500, 501.0, 2**40], dtype=object), "x")
        assert out.tolist() == [500, 501, 2**40]

    def test_empty_is_empty(self):
        assert whole_numbers([], "x").dtype == np.int64


class TestRefused:
    @pytest.mark.parametrize(
        "values, first",
        [
            ([1.0, 2.5], "2.5"),
            ([1.0, float("nan")], "nan"),
            ([1.0, float("inf")], "inf"),
            ([1.0, 1.0e19], "1e+19"),
            ([1.0, 2.0**63], "9.223372036854776e+18"),
        ],
    )
    def test_floats_that_are_no_position(self, values, first):
        with pytest.raises(ConversionRefused, match=r"1 of 2 x positions in a\.mzpeak"):
            whole_numbers(values, "x positions", "a.mzpeak")
        with pytest.raises(ConversionRefused, match=re.escape(f"the first is {first}")):
            whole_numbers(values, "x positions")

    def test_uint64_at_or_above_2_63(self):
        values = np.array([5, 2**63, 2**64 - 1], dtype=np.uint64)
        with pytest.raises(
            ConversionRefused, match=r"2 of 3 .* first is 9223372036854775808"
        ):
            whole_numbers(values, "x positions")

    @pytest.mark.parametrize("bad", ["501a", None, 1.0e19, True, 2.5])
    def test_sqlite_objects_that_are_no_position(self, bad):
        with pytest.raises(ConversionRefused, match=re.escape(f"first is {bad!r}")):
            whole_numbers(np.array([500, bad], dtype=object), "x positions")

    def test_no_tolerance_for_a_column_declared_integer(self):
        assert whole_numbers([100.0], "x", tolerance=0.0).tolist() == [100]
        with pytest.raises(ConversionRefused, match="100.9999999"):
            whole_numbers([100.9999999], "x", tolerance=0.0)
        with pytest.raises(ConversionRefused, match="100.9999999"):
            whole_numbers(np.array([100.9999999], dtype=object), "x", tolerance=0.0)

    def test_a_bool_column(self):
        with pytest.raises(ConversionRefused):
            whole_numbers(np.array([True, False]), "x positions")
