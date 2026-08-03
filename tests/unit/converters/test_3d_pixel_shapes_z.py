# tests/unit/converters/test_3d_pixel_shapes_z.py
"""A volume's pixel footprints sit at the depth of the slice they came from.

``_create_pixel_shapes`` took an ``is_3d`` argument and never read it. Every
slice's footprints were built from ``spatial_x``/``spatial_y`` alone, so a
volume's polygons all stacked flat at z=0 while its TIC image spanned
``n_z * z_spacing_um`` in depth. ``docs/coordinate-systems.md`` promises that
every element in a Thyra store resolves to the same physical frame at
``"global"``; for a volume that promise was false in z, and PR #160 (which made
the *image* honest) did not change it.

Footprints now carry z, taken from ``obs["spatial_z"]`` -- already micrometres,
and the same column the table stores, which is what makes the two elements
agree rather than merely both look plausible.

**Why not the obvious alternative.** spatialdata's shapes model is 2D, which
invites keeping flat geometry and putting the depth in a transformation (one
element per slice, each with a ``Translation`` in z). Measured, that silently
loses: the element parses without complaint, coexists with a 3D image and
survives a zarr round trip, but the transform **drops z** and every slice comes
back at the same depth -- a depth recorded in the metadata that never reaches
the geometry. ``test_a_translation_would_not_have_worked`` pins that, because
it is the change someone will otherwise propose again.

These are footprints, not voxels: a flat square at one depth, since shapely has
no solid. The pixel's extent in z is the slice spacing and is carried by the
image, not by the polygon.
"""

from __future__ import annotations

from typing import Any, Dict

import numpy as np
import pytest

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.converters.spatialdata.base_spatialdata_converter import (
    SPATIALDATA_AVAILABLE,
)
from thyra.converters.spatialdata.spatialdata_2d_converter import SpatialData2DConverter
from thyra.converters.spatialdata.spatialdata_3d_converter import SpatialData3DConverter

pytestmark = pytest.mark.skipif(
    not SPATIALDATA_AVAILABLE,
    reason="SpatialData dependencies not available",
)

_DATASET_ID = "vol"

# Pairwise different, as elsewhere in the 3D tests: no axis length can stand in
# for another and make a wrong answer look right.
_N_X, _N_Y, _N_Z = 5, 3, 2

# 25x apart, so a footprint placed with the in-plane pitch instead of the slice
# spacing is off by more than any rounding could explain.
_PIXEL_SIZE_UM = 10.0
_Z_SPACING_UM = 250.0


def _config(n_z: int = _N_Z) -> MockMSIConfig:
    """A full grid -- every voxel present, so obs covers the whole volume."""
    return MockMSIConfig(
        n_x=_N_X,
        n_y=_N_Y,
        n_z=n_z,
        n_mz_bins=64,
        peaks_per_spectrum=(4, 8),
        sparsity=0.0,
    )


def _convert_volume(tmp_path, n_z: int = _N_Z, **kwargs: Any):
    output_path = tmp_path / "out.zarr"
    converter = SpatialData3DConverter(
        MockMSIReader(_config(n_z)),
        output_path,
        dataset_id=_DATASET_ID,
        pixel_size_um=_PIXEL_SIZE_UM,
        z_spacing_um=_Z_SPACING_UM,
        **kwargs,
    )
    assert converter.convert() is True
    return output_path


@pytest.fixture(scope="module")
def volume(tmp_path_factory) -> Dict[str, Any]:
    """One converted volume, read back off disk."""
    import spatialdata

    store = _convert_volume(tmp_path_factory.mktemp("shapes_z"))
    sdata = spatialdata.read_zarr(str(store))
    return {
        "shapes": sdata.shapes[f"{_DATASET_ID}_pixels"],
        "obs": sdata.tables[_DATASET_ID].obs,
        "image": sdata.images[f"{_DATASET_ID}_tic"],
    }


def _depth(geometry) -> float:
    """The z ordinate of a footprint. Flat, so any vertex will do."""
    return float(geometry.exterior.coords[0][2])


def test_volume_footprints_carry_z(volume):
    """The reproduction. Before, these were 2D and every depth was absent."""
    shapes = volume["shapes"]

    assert len(shapes) == _N_X * _N_Y * _N_Z
    assert all(geometry.has_z for geometry in shapes.geometry)


def test_each_footprint_sits_at_its_own_slice_depth(volume):
    """Depth must be the slice index times the spacing -- in closed form.

    Checked against the arithmetic rather than against the store, so this
    fails if the converter and the expectation drift apart. The set of depths
    is compared too: a bug that put every footprint on plane 0 would satisfy a
    per-row check that only ever looked at plane 0.
    """
    shapes, obs = volume["shapes"], volume["obs"]

    expected_depths = [z * _Z_SPACING_UM for z in range(_N_Z)]
    assert sorted({_depth(g) for g in shapes.geometry}) == expected_depths

    for label, geometry in zip(shapes.index, shapes.geometry):
        z_index = int(obs.loc[label, "z"])
        assert _depth(geometry) == pytest.approx(z_index * _Z_SPACING_UM)


def test_footprints_and_table_agree_at_global(volume):
    """The two elements a consumer overlays must land at the same depth.

    obs["spatial_z"] and the polygons are built in different places. Either
    alone can be self-consistent while the pair disagrees, which is the class
    of defect this whole area keeps producing.
    """
    shapes, obs = volume["shapes"], volume["obs"]

    for label, geometry in zip(shapes.index, shapes.geometry):
        assert _depth(geometry) == pytest.approx(float(obs.loc[label, "spatial_z"]))


def test_footprints_and_image_span_the_same_depth(volume):
    """The volume's depth extent and the footprints' must be the same range.

    The image spans n_z planes scaled by z_spacing; the deepest footprint sits
    on the last plane. If either the Scale or the footprint depth used the
    in-plane pitch, these part company.
    """
    from spatialdata.transformations import get_transformation

    shapes, image = volume["shapes"], volume["image"]

    n_z_planes = np.asarray(image.data).shape[1]
    assert n_z_planes == _N_Z

    axes = ("c", "z", "y", "x")
    matrix = get_transformation(image, "global").to_affine_matrix(
        input_axes=axes, output_axes=axes
    )
    z_scale = float(matrix[1, 1])

    assert z_scale == pytest.approx(_Z_SPACING_UM)
    assert max(_depth(g) for g in shapes.geometry) == pytest.approx(
        (n_z_planes - 1) * z_scale
    )


def test_the_footprint_is_still_pixel_sized(volume):
    """Adding depth must not disturb the in-plane extent.

    A footprint is the pixel's square at one depth; it is not a voxel, and its
    x/y size is still the in-plane pitch.
    """
    geometry = volume["shapes"].geometry.iloc[0]
    xs = [c[0] for c in geometry.exterior.coords]
    ys = [c[1] for c in geometry.exterior.coords]

    assert max(xs) - min(xs) == pytest.approx(_PIXEL_SIZE_UM)
    assert max(ys) - min(ys) == pytest.approx(_PIXEL_SIZE_UM)


def test_a_single_slice_volume_stays_flat(tmp_path):
    """``handle_3d`` alone is not a volume, so it gets no z.

    ``SpatialData3DConverter`` forces ``handle_3d=True`` regardless of the
    data, and a one-plane acquisition through it takes the 2D branch
    everywhere else. Footprints must follow, or a flat dataset acquires a
    depth axis it has no planes for.
    """
    import spatialdata

    store = _convert_volume(tmp_path, n_z=1)
    sdata = spatialdata.read_zarr(str(store))
    shapes = sdata.shapes[f"{_DATASET_ID}_pixels"]

    assert not any(geometry.has_z for geometry in shapes.geometry)


def test_the_2d_route_is_unchanged(tmp_path):
    """The per-slice route writes flat elements and must keep doing so.

    It emits one table and one image per plane, so each element genuinely has
    no third dimension. A z here would be inventing one.
    """
    import spatialdata

    output_path = tmp_path / "flat.zarr"
    converter = SpatialData2DConverter(
        MockMSIReader(_config()),
        output_path,
        dataset_id=_DATASET_ID,
        pixel_size_um=_PIXEL_SIZE_UM,
    )
    assert converter.convert() is True

    sdata = spatialdata.read_zarr(str(output_path))
    for key, shapes in sdata.shapes.items():
        assert not any(
            geometry.has_z for geometry in shapes.geometry
        ), f"{key} grew a z coordinate"


def test_a_translation_would_not_have_worked():
    """Pins the measured negative result behind the POLYGON Z choice.

    Keeping flat 2D geometry and expressing the depth as a ``Translation`` is
    the natural reading of spatialdata's 2D shapes model, and it is the change
    someone will proposeagain on seeing the UserWarning that POLYGON Z
    provokes. It does not work, and it fails *silently*: the transform drops z
    and returns flat geometry, so the depth lives only in the metadata.

    If a future spatialdata makes this work, this test fails -- and that is the
    signal to revisit the choice, not a regression.
    """
    import geopandas as gpd
    import spatialdata
    from shapely.geometry import box
    from spatialdata.models import ShapesModel
    from spatialdata.transformations import Scale, Sequence, Translation

    gdf = gpd.GeoDataFrame(geometry=[box(0.0, 0.0, 1.0, 1.0)], index=["p0"])
    parsed = ShapesModel.parse(
        gdf,
        transformations={
            "global": Sequence(
                [
                    Scale([_PIXEL_SIZE_UM, _PIXEL_SIZE_UM], axes=("x", "y")),
                    Translation([2 * _Z_SPACING_UM], axes=("z",)),
                ]
            )
        },
    )

    transformed = spatialdata.transform(parsed, to_coordinate_system="global")
    geometry = transformed.geometry.iloc[0]

    # x/y scaling lands; the z translation vanishes without a word.
    assert max(c[0] for c in geometry.exterior.coords) == pytest.approx(_PIXEL_SIZE_UM)
    assert not geometry.has_z, (
        "spatialdata now carries z through a Translation on 2D shapes; "
        "revisit whether POLYGON Z is still the right representation"
    )
