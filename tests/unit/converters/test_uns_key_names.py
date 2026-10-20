"""Every key ``uns`` holds is a name Zarr can store on every platform.

Zarr writes a dict key as the name of a node, a folder on disk. A colon is
not legal there on Windows, so a CV accession used as a key failed the save
after the whole conversion was done, and a slash split the key into a path.
Keys that come from a source file, an mzPeak archive's free-form metadata
blocks for one, can hold either.
"""

import json

import pytest

from thyra.converters.spatialdata.uns_assembler import (
    _jsonify_string_lists,
    serialize_for_zarr,
    zarr_safe_key,
)


@pytest.mark.parametrize(
    ("key", "safe"),
    [
        ("MS:1000529", "MS_1000529"),
        ("a/b", "a_b"),
        ("a\\b", "a_b"),
        ('what? "*" <x> | y', "what_ ___ _x_ _ y"),
        ("tab\there", "tab_here"),
        ("ends with a dot.", "ends with a dot_"),
        ("..", "__"),
        ("", "_"),
        ("plain_key", "plain_key"),
        ("pixel size (x)", "pixel size (x)"),
        (7, "7"),
    ],
)
def test_a_key_becomes_a_safe_name(key, safe):
    """Each unsafe character becomes ``_``; a safe key is kept."""
    assert zarr_safe_key(key) == safe


def test_keys_are_made_safe_at_every_level():
    """The run block of an mzPeak archive, one level down."""
    result = serialize_for_zarr({"run": {"id": "r", "MS:1000529": "s1"}})

    assert result == {"run": {"id": "r", "MS_1000529": "s1"}}


def test_a_safe_key_keeps_its_name_when_another_becomes_it():
    """Two keys that meet are told apart, and no value is lost."""
    result = serialize_for_zarr({"MS:1": "colon", "MS_1": "underscore"})

    assert result == {"MS_1_2": "colon", "MS_1": "underscore"}


def test_keys_inside_a_list_are_kept():
    """A list of dicts is stored as one JSON string, where any key is legal."""
    params = [{"MS:1000529": "s1", "a/b": 1}]
    result = _jsonify_string_lists(serialize_for_zarr({"params": params}))

    assert json.loads(result["params"]) == params


def test_an_object_is_walked_as_its_dict():
    """A vendor object is stored as its attributes, keys made safe."""

    class Vendor:
        def __init__(self):
            self.__dict__["key:with:colons"] = 1

    assert serialize_for_zarr({"v": Vendor()}) == {"v": {"key_with_colons": 1}}
