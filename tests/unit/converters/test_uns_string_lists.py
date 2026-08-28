"""The JSON rule for non-numeric lists headed for ``uns``.

AnnData/zarr cannot round-trip a list of dicts (each entry is
stringified into ``repr`` output) and materializes any list of strings
as a numpy string array on read-back -- and deepcopying such an array
segfaults the process outright on numpy 2.1-2.2 (numpy#28609). Every
table copy deepcopies ``uns``, so one stray list of strings written by a
vendor extractor is enough to kill ``AnnData.copy`` or any spatial query
for readers pinned to that numpy range.

``_jsonify_string_lists`` is the write-side guard: lists that are not
purely numeric are stored as JSON strings instead. These tests pin the
rule itself; the round-trip through real stores is asserted in
``test_uns_provenance_parity.py``.
"""

import json

import numpy as np
import pytest

from thyra.converters.spatialdata.base_spatialdata_converter import (
    _jsonify_string_lists,
)


def test_list_of_dicts_becomes_json():
    cv_params = [
        {"name": "MS1 spectrum", "accession": "MS:1000579", "value": True},
        {
            "name": "scan start time",
            "accession": "MS:1000016",
            "value": 0.5,
            "unit_name": "minute",
            "unit_accession": "UO:0000031",
        },
    ]
    result = _jsonify_string_lists({"cvParams": cv_params})
    assert isinstance(result["cvParams"], str)
    assert json.loads(result["cvParams"]) == cv_params


def test_list_of_strings_becomes_json():
    result = _jsonify_string_lists({"files": ["a.imzML", "a.ibd"]})
    assert json.loads(result["files"]) == ["a.imzML", "a.ibd"]


@pytest.mark.parametrize(
    "value",
    [
        [1, 2, 3],
        [1.5, 2.5],
        [True, False],
        [[100.0, 200.0], [300.0, 400.0]],
        [np.float64(1.0), np.int32(2)],
        [],
    ],
)
def test_numeric_lists_stay_lists(value):
    assert _jsonify_string_lists({"key": value})["key"] is value


def test_mixed_list_becomes_json():
    result = _jsonify_string_lists({"key": [1, "a"]})
    assert json.loads(result["key"]) == [1, "a"]


def test_none_bearing_list_becomes_json():
    # np.array([1, None]) is an object array, which the writer cannot
    # store faithfully either -- JSON keeps the null.
    result = _jsonify_string_lists({"key": [1, None]})
    assert json.loads(result["key"]) == [1, None]


def test_recursion_reaches_nested_dicts():
    obj = {"outer": {"inner": [{"a": 1}], "n": [1, 2]}}
    result = _jsonify_string_lists(obj)
    assert json.loads(result["outer"]["inner"]) == [{"a": 1}]
    assert result["outer"]["n"] == [1, 2]


def test_scalars_and_arrays_pass_through():
    arr = np.arange(3)
    obj = {"s": "text", "f": 1.5, "arr": arr}
    result = _jsonify_string_lists(obj)
    assert result["s"] == "text"
    assert result["f"] == 1.5
    assert result["arr"] is arr


def test_exotic_values_inside_lists_do_not_raise():
    # A vendor extractor may report anything; the dump must never take
    # the conversion down. Unknown objects fall back to str().
    class Odd:
        def __repr__(self):
            return "Odd()"

    result = _jsonify_string_lists({"key": ["x", Odd(), np.arange(2)]})
    assert json.loads(result["key"]) == ["x", "Odd()", [0, 1]]
