"""The call pattern of a downstream project that builds on Thyra's store.

A project outside this repository converts a file with one fixed call and
then reads the store back by element name. This file pins that call and
what it reads, so that a change that breaks them is seen here and not by
that project. It sits beside ``test_read_lazy_contract.py``, which pins how
the same store is read lazily.

The call is ``convert_msi`` with the keywords ``input_path``,
``output_path``, ``format_type``, ``dataset_id`` and ``resampling_config``,
and it relies on a truthy return that is ``True`` and not just truthy.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest
import spatialdata

from thyra.convert import convert_msi
from thyra.errors import ConversionRefused

_FIXTURE = (
    Path(__file__).resolve().parent.parent
    / "data"
    / "fixtures"
    / "pixel_size_unit_declared.imzML"
)


def _convert(out: Path, resampling_config, **kwargs):
    return convert_msi(
        input_path=_FIXTURE,
        output_path=out,
        format_type="spatialdata",
        dataset_id="scan_a",
        resampling_config=resampling_config,
        **kwargs,
    )


@pytest.fixture(scope="module")
def resampled_store(tmp_path_factory):
    out = tmp_path_factory.mktemp("downstream") / "scan_a.zarr"
    assert _convert(out, {}) is True
    return out


@pytest.fixture(scope="module")
def unresampled_store(tmp_path_factory):
    out = tmp_path_factory.mktemp("downstream") / "scan_a_raw.zarr"
    assert _convert(out, None) is True
    return out


def test_the_keyword_names_exist():
    parameters = inspect.signature(convert_msi).parameters
    for name in ("input_path", "output_path", "format_type", "dataset_id"):
        assert name in parameters
    assert "resampling_config" in parameters


def test_the_return_value_is_true(tmp_path):
    result = _convert(tmp_path / "out.zarr", {})
    assert result
    assert result is True


class TestWhatTheStoreHolds:
    def test_the_element_names(self, resampled_store):
        sdata = spatialdata.read_zarr(resampled_store)
        assert list(sdata.tables) == ["scan_a_z0"]
        assert list(sdata.images) == ["scan_a_z0_tic"]
        assert list(sdata.shapes) == ["scan_a_z0_pixels"]

    def test_the_root_attributes(self, resampled_store):
        attrs = spatialdata.read_zarr(resampled_store).attrs
        info = attrs["pixel_size_detection_info"]
        assert info["detected_x_um"] == 50.0
        assert info["detected_y_um"] == 50.0
        assert attrs["pixel_size_x_um"] == 50.0
        assert attrs["pixel_size_y_um"] == 50.0

    def test_the_table_names_its_shapes(self, resampled_store):
        table = spatialdata.read_zarr(resampled_store).tables["scan_a_z0"]
        assert table.uns["spatialdata_attrs"]["region"] == "scan_a_z0_pixels"


class TestTheResamplingArgument:
    def test_an_empty_dict_resamples(self, resampled_store, unresampled_store):
        resampled = spatialdata.read_zarr(resampled_store).tables["scan_a_z0"]
        raw = spatialdata.read_zarr(unresampled_store).tables["scan_a_z0"]
        assert resampled.n_vars > raw.n_vars

    def test_none_is_accepted(self, unresampled_store):
        assert unresampled_store.is_dir()


class TestAnUnknownKeyword:
    def test_it_returns_false_by_default(self, tmp_path):
        assert _convert(tmp_path / "out.zarr", {}, bogus_keyword=1) is False

    def test_it_raises_a_refusal_naming_it_on_request(self, tmp_path):
        with pytest.raises(ConversionRefused, match="'bogus_keyword'"):
            _convert(tmp_path / "out.zarr", {}, bogus_keyword=1, raise_on_failure=True)
