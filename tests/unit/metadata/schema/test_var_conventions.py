"""The fixed var column contract, checked straight off the zarr layout.

These tests fabricate minimal zarr stores rather than converting, so
they run without the SpatialData stack and cover the failure shapes a
conversion can never produce.
"""

import numpy as np
import pytest
import zarr

from thyra.metadata.schema import check_store_var_conventions


def _store_with_var(tmp_path, mz=None, extra=None, skip_var=False):
    root = zarr.open_group(str(tmp_path / "store.zarr"), mode="a")
    table = root.create_group("tables").create_group("t")
    if not skip_var:
        var = table.create_group("var")
        if mz is not None:
            var.create_array("mz", data=np.asarray(mz))
        if extra is not None:
            var.create_array(extra, data=np.zeros(3))
    return tmp_path / "store.zarr"


def _messages(issues):
    return " ".join(i.message for i in issues)


class TestCheckStoreVarConventions:
    def test_conforming_var_has_no_issues(self, tmp_path):
        store = _store_with_var(tmp_path, mz=[100.0, 200.5, 300.0])
        assert check_store_var_conventions(store) == {"t": []}

    def test_missing_var_group_is_an_error(self, tmp_path):
        store = _store_with_var(tmp_path, skip_var=True)
        issues = check_store_var_conventions(store)["t"]
        assert issues and issues[0].severity == "error"

    def test_missing_mz_column_is_an_error(self, tmp_path):
        store = _store_with_var(tmp_path, mz=None, extra="flight_time")
        issues = check_store_var_conventions(store)["t"]
        assert "'mz' is missing" in _messages(issues)

    def test_non_monotonic_mz_is_an_error(self, tmp_path):
        store = _store_with_var(tmp_path, mz=[100.0, 300.0, 200.0])
        issues = check_store_var_conventions(store)["t"]
        assert "strictly increasing" in _messages(issues)

    def test_non_finite_mz_is_an_error(self, tmp_path):
        store = _store_with_var(tmp_path, mz=[100.0, np.nan, 300.0])
        issues = check_store_var_conventions(store)["t"]
        assert "non-finite" in _messages(issues)

    def test_not_a_store_raises(self, tmp_path):
        zarr.open_group(str(tmp_path / "plain.zarr"), mode="a")
        with pytest.raises(ValueError, match="tables"):
            check_store_var_conventions(tmp_path / "plain.zarr")


class TestReservedColumnNames:
    def test_spec_reserves_the_annotation_columns(self):
        from thyra.metadata.schema import (
            MSI_VAR_REQUIRED_COLUMNS,
            MSI_VAR_RESERVED_COLUMNS,
        )

        assert MSI_VAR_REQUIRED_COLUMNS == ("mz",)
        assert set(MSI_VAR_RESERVED_COLUMNS) == {
            "mz",
            "formula",
            "adduct",
            "annotation_source",
            "fdr",
        }
