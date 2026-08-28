"""``sanitize_uns_string_arrays`` -- the read-side defense for old stores.

Stores written before the JSON rule carry string lists in ``uns`` that
read back as numpy ``StringDType`` arrays, and deepcopying one of those
segfaults on numpy 2.1-2.2 (numpy#28609). The helper swaps them for
plain lists so an old store stays usable in a pinned environment. It
must touch nothing else: the point is that assigning the result back
over ``table.uns`` changes only what would otherwise crash.
"""

import copy

import numpy as np
import pytest

from thyra.metadata import sanitize_uns_string_arrays

_HAS_STRINGDTYPE = hasattr(np.dtypes, "StringDType")

requires_stringdtype = pytest.mark.skipif(
    not _HAS_STRINGDTYPE, reason="numpy build has no StringDType"
)


@requires_stringdtype
def test_stringdtype_arrays_become_lists():
    arr = np.array(["{'name': 'MS1 spectrum'}", "b"], dtype=np.dtypes.StringDType())
    uns = {"raw_metadata": {"cvParams": arr}}

    result = sanitize_uns_string_arrays(uns)

    assert result["raw_metadata"]["cvParams"] == ["{'name': 'MS1 spectrum'}", "b"]
    assert isinstance(result["raw_metadata"]["cvParams"], list)
    # The crash path: this deepcopy is what polygon_query and
    # AnnData.copy do to uns, and what segfaults on numpy 2.1-2.2
    # before sanitizing.
    copy.deepcopy(result)


@requires_stringdtype
def test_recursion_reaches_arrays_inside_lists():
    arr = np.array(["x"], dtype=np.dtypes.StringDType())
    result = sanitize_uns_string_arrays({"outer": [{"inner": arr}]})
    assert result["outer"][0]["inner"] == ["x"]


def test_everything_else_is_untouched():
    numeric = np.arange(3)
    unicode_arr = np.array(["a", "b"])  # kind "U": deepcopy-safe, keep as-is
    uns = {
        "essential_metadata": {"dimensions": numeric, "spectrum_type": "centroid"},
        "regions": '[{"region_number": 1}]',
        "unicode": unicode_arr,
        "average_spectrum": numeric,
    }

    result = sanitize_uns_string_arrays(uns)

    assert result["essential_metadata"]["dimensions"] is numeric
    assert result["essential_metadata"]["spectrum_type"] == "centroid"
    assert result["regions"] == uns["regions"]
    assert result["unicode"] is unicode_arr
    assert result["average_spectrum"] is numeric


def test_accepts_any_mapping_and_returns_plain_dict():
    # AnnData's uns is an OverloadedDict, not a dict -- the helper takes
    # any Mapping and hands back a plain dict fit for assignment.
    class MappingLike(dict):
        pass

    result = sanitize_uns_string_arrays(MappingLike({"a": {"b": 1}}))
    assert type(result) is dict
    assert type(result["a"]) is dict
    assert result == {"a": {"b": 1}}
