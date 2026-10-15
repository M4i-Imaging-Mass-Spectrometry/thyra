"""Reading the pyramid level paths out of an OME-Zarr image's attributes.

The streamed optical image reopens the element SpatialData just wrote and
takes its level paths from ``ome.multiscales[0].datasets``. Anything not
shaped like that comes back as None, which the caller reports as an
image that carries no multiscales metadata.
"""

import pytest

from thyra.converters.spatialdata.optical_image import _multiscale_dataset_paths


def _attrs(datasets):
    return {"ome": {"multiscales": [{"datasets": datasets}]}}


def test_the_paths_come_back_in_order():
    attrs = _attrs([{"path": "0"}, {"path": "1"}, {"path": "2"}])
    assert _multiscale_dataset_paths(attrs) == ["0", "1", "2"]


@pytest.mark.parametrize(
    "attrs",
    [
        pytest.param({}, id="no_ome"),
        pytest.param({"ome": []}, id="ome_not_a_dict"),
        pytest.param({"ome": {}}, id="no_multiscales"),
        pytest.param({"ome": {"multiscales": []}}, id="multiscales_empty"),
        pytest.param({"ome": {"multiscales": ["0"]}}, id="multiscale_not_a_dict"),
        pytest.param({"ome": {"multiscales": [{}]}}, id="no_datasets"),
        pytest.param(_attrs({"path": "0"}), id="datasets_not_a_list"),
        pytest.param(_attrs([{"path": "0"}, "1"]), id="dataset_not_a_dict"),
        pytest.param(_attrs([{"path": "0"}, {}]), id="dataset_without_path"),
        pytest.param(_attrs([{"path": 0}]), id="path_not_a_string"),
    ],
)
def test_anything_else_is_none(attrs):
    assert _multiscale_dataset_paths(attrs) is None
