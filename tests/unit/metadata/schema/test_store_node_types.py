"""An array where the store layout puts a group.

A converter never writes one, so these stores are fabricated. The readers
treat the array as a missing group, which turns three crashes with
unrelated messages into the answer a missing group already gets.
"""

import numpy as np
import pytest
import zarr

from thyra.metadata.schema import check_store_var_conventions
from thyra.metadata.schema.store_io import read_msi_metadata_blocks, subgroup


def _root(tmp_path):
    return zarr.open_group(str(tmp_path / "store.zarr"), mode="a")


class TestSubgroup:
    def test_a_group_comes_back(self, tmp_path):
        root = _root(tmp_path)
        root.create_group("tables").create_group("t")
        assert isinstance(subgroup(root, "tables/t"), zarr.Group)

    def test_an_array_reads_as_absent(self, tmp_path):
        root = _root(tmp_path)
        root.create_array("tables", shape=(2,), dtype="i4")
        assert subgroup(root, "tables") is None

    def test_nothing_reads_as_absent(self, tmp_path):
        assert subgroup(_root(tmp_path), "tables") is None


class TestTablesIsAnArray:
    def test_reading_says_there_is_no_tables_group(self, tmp_path):
        _root(tmp_path).create_array("tables", shape=(2,), dtype="i4")
        with pytest.raises(ValueError, match="no 'tables' group"):
            read_msi_metadata_blocks(tmp_path / "store.zarr")

    def test_checking_says_there_is_no_tables_group(self, tmp_path):
        _root(tmp_path).create_array("tables", shape=(2,), dtype="i4")
        with pytest.raises(ValueError, match="no 'tables' group"):
            check_store_var_conventions(tmp_path / "store.zarr")


class TestTableIsAnArray:
    def test_reading_skips_it(self, tmp_path):
        root = _root(tmp_path)
        root.create_group("tables").create_array("t", shape=(2,), dtype="i4")
        assert read_msi_metadata_blocks(tmp_path / "store.zarr") == {}

    def test_checking_reports_no_var_group(self, tmp_path):
        root = _root(tmp_path)
        root.create_group("tables").create_array("t", shape=(2,), dtype="i4")
        issues = check_store_var_conventions(tmp_path / "store.zarr")["t"]
        assert [i.message for i in issues] == ["table has no var group"]


def test_a_var_array_reports_no_var_group(tmp_path):
    table = _root(tmp_path).create_group("tables").create_group("t")
    table.create_array("var", data=np.array([1.0, 2.0, 3.0]))
    issues = check_store_var_conventions(tmp_path / "store.zarr")["t"]
    assert [i.message for i in issues] == ["table has no var group"]


def test_a_block_that_is_not_a_mapping_is_named(tmp_path):
    uns = _root(tmp_path).create_group("tables").create_group("t").create_group("uns")
    uns.create_array("msi_metadata", shape=(2,), dtype="i4")
    with pytest.raises(TypeError, match="msi_metadata is a ndarray, not a mapping"):
        read_msi_metadata_blocks(tmp_path / "store.zarr")
