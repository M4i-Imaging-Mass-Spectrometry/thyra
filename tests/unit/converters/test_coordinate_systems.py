# tests/unit/converters/test_coordinate_systems.py
"""Tests for the structured coordinate-system contract emitted by the
SpatialData converters.

The contract guarantees that within one zarr produced by Thyra:

  * Both the TIC image and the pixel-polygon shapes resolve to the same
    space at the ``"global"`` coordinate system. In Mode A (no
    FlexImaging optical alignment) that space is physical micrometers.

  * A structured ``coordinate_systems`` attr is written at the zarr
    top level so consumers (Ousia, registration tooling, etc.) do not
    need to guess what ``"global"`` is in.

These tests exercise the code path end-to-end without mocking out the
transform plumbing, so they catch regressions where image and shapes
silently disagree about what ``"global"`` means.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from thyra.convert import convert_msi

pytestmark = pytest.mark.skipif(
    not pytest.importorskip("spatialdata", reason="SpatialData not installed"),
    reason="SpatialData not installed",
)


def _bbox_at_global(element) -> tuple[float, float, float, float]:
    """Compute the bbox of a SpatialData element in its ``"global"`` CS.

    Works for both shapes (GeoDataFrame) and images (xarray DataArray /
    DataTree). Returns ``(xmin, ymin, xmax, ymax)`` after applying the
    element's transform to ``"global"``.
    """
    import numpy as np
    from spatialdata.transformations import get_transformation

    transform = get_transformation(element, to_coordinate_system="global")
    matrix = transform.to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))

    # Shapes: native bbox is just the GeoDataFrame's total_bounds.
    if hasattr(element, "geometry"):
        xmin_n, ymin_n, xmax_n, ymax_n = element.total_bounds
    else:
        # Image / multiscale image: native frame is integer pixel indices.
        # Pull the highest-resolution data to get x/y sizes.
        if hasattr(element, "shape"):
            data = element
        else:
            # DataTree: pick the first scale's "image" variable.
            scales = list(element.children.keys())
            data = element[scales[0]].ds["image"]
        # Dims convention is ("c", "y", "x") for 2D images.
        y_size = data.sizes["y"]
        x_size = data.sizes["x"]
        xmin_n, ymin_n = 0.0, 0.0
        xmax_n, ymax_n = float(x_size), float(y_size)

    corners = np.array(
        [
            [xmin_n, ymin_n, 1.0],
            [xmax_n, ymin_n, 1.0],
            [xmin_n, ymax_n, 1.0],
            [xmax_n, ymax_n, 1.0],
        ]
    )
    transformed = corners @ matrix.T
    xs = transformed[:, 0]
    ys = transformed[:, 1]
    return float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())


def _read_zarr_attrs(zarr_path: Path) -> dict:
    """Read top-level zarr attrs as a plain dict.

    Supports both zarr v2 (``.zattrs``) and v3 (``zarr.json``) layouts.
    """
    z3 = zarr_path / "zarr.json"
    if z3.exists():
        with open(z3) as f:
            doc = json.load(f)
        return doc.get("attributes", {}) or {}
    z2 = zarr_path / ".zattrs"
    if z2.exists():
        with open(z2) as f:
            return json.load(f) or {}
    raise FileNotFoundError(
        f"No zarr attrs file found in {zarr_path}; expected zarr.json or .zattrs"
    )


class TestCoordinateSystemsContract:
    """Contract-level tests that survive across converter implementations."""

    def test_no_alignment_image_and_shapes_agree_at_global(
        self, create_minimal_imzml, temp_dir
    ):
        """Mode A: image and pixel-polygon shapes must resolve to the
        same bbox at ``"global"`` (within half a pixel) when no optical
        alignment is provided.
        """
        import spatialdata as sd

        imzml_path, _, _, _ = create_minimal_imzml
        output_path = temp_dir / "out.zarr"
        pixel_size_um = 2.5

        result = convert_msi(
            str(imzml_path),
            str(output_path),
            format_type="spatialdata",
            dataset_id="ds",
            pixel_size_um=pixel_size_um,
        )
        assert result is True

        sdata = sd.SpatialData.read(str(output_path))
        assert len(sdata.images) >= 1
        assert len(sdata.shapes) >= 1

        image_name = next(iter(sdata.images))
        shapes_name = next(iter(sdata.shapes))
        image_bbox = _bbox_at_global(sdata.images[image_name])
        shapes_bbox = _bbox_at_global(sdata.shapes[shapes_name])

        # Pixel-polygon shapes are box-padded by half a pixel on every
        # side relative to the image grid (a centroid at (0,0) becomes a
        # box from (-px/2, -px/2) to (+px/2, +px/2)). Allow that tolerance.
        tol = pixel_size_um
        assert abs(image_bbox[0] - shapes_bbox[0]) < tol, (
            f"xmin mismatch at global: image={image_bbox[0]}, "
            f"shapes={shapes_bbox[0]}"
        )
        assert abs(image_bbox[1] - shapes_bbox[1]) < tol, (
            f"ymin mismatch at global: image={image_bbox[1]}, "
            f"shapes={shapes_bbox[1]}"
        )
        assert abs(image_bbox[2] - shapes_bbox[2]) < tol, (
            f"xmax mismatch at global: image={image_bbox[2]}, "
            f"shapes={shapes_bbox[2]}"
        )
        assert abs(image_bbox[3] - shapes_bbox[3]) < tol, (
            f"ymax mismatch at global: image={image_bbox[3]}, "
            f"shapes={shapes_bbox[3]}"
        )

    def test_no_alignment_writes_micrometer_coordinate_systems_attr(
        self, create_minimal_imzml, temp_dir
    ):
        """Mode A: the structured ``coordinate_systems.global`` attr
        must be present and declare ``unit="micrometer"`` with the MSI
        pixel size filled in.
        """
        imzml_path, _, _, _ = create_minimal_imzml
        output_path = temp_dir / "out.zarr"
        pixel_size_um = 3.0

        convert_msi(
            str(imzml_path),
            str(output_path),
            format_type="spatialdata",
            dataset_id="ds",
            pixel_size_um=pixel_size_um,
        )

        attrs = _read_zarr_attrs(output_path)
        assert (
            "coordinate_systems" in attrs
        ), "Thyra zarrs must declare a structured coordinate_systems attr"
        cs_global = attrs["coordinate_systems"]["global"]
        assert cs_global["unit"] == "micrometer"
        assert cs_global["pixel_size_um_x"] == pytest.approx(pixel_size_um)
        assert cs_global["pixel_size_um_y"] == pytest.approx(pixel_size_um)
        assert cs_global["reference_element"] is None
        assert cs_global["convention_version"] >= 1
        assert cs_global["produced_by"].startswith("thyra/")

    def test_image_intrinsic_frame_is_pixel_indices(
        self, create_minimal_imzml, temp_dir
    ):
        """The TIC image must be stored intrinsically in raster pixel
        indices, independent of ``pixel_size_um``. The unit conversion
        is expressed purely through the transform to ``"global"``.

        SpatialData's ``Image2DModel.parse`` stamps half-pixel-centered
        coords (0.5, 1.5, 2.5, ...) when no coords are provided, which
        is still pixel-index space; the key invariant is that the step
        between consecutive pixels is exactly 1.0 and the values do not
        scale with ``pixel_size_um``.
        """
        import spatialdata as sd

        imzml_path, _, _, _ = create_minimal_imzml
        output_path = temp_dir / "out.zarr"
        pixel_size_um = 7.5  # deliberately not 1.0 so any um-leak shows up
        convert_msi(
            str(imzml_path),
            str(output_path),
            format_type="spatialdata",
            dataset_id="ds",
            pixel_size_um=pixel_size_um,
        )

        sdata = sd.SpatialData.read(str(output_path))
        image_name = next(iter(sdata.images))
        image = sdata.images[image_name]

        # Resolve to the underlying DataArray (DataTree if multiscale).
        if hasattr(image, "shape"):
            data = image
        else:
            scales = list(image.children.keys())
            data = image[scales[0]].ds["image"]

        x_coord = data.coords["x"].values
        y_coord = data.coords["y"].values

        # Step must be exactly 1.0 (pixel-index spacing), not the um
        # pixel size. This is the regression guard against the old
        # behavior where xarray coords were stamped in micrometers.
        if len(x_coord) > 1:
            steps_x = x_coord[1:] - x_coord[:-1]
            assert (
                steps_x == 1
            ).all(), f"image x coord must step by 1 (pixel index); got {steps_x}"
        if len(y_coord) > 1:
            steps_y = y_coord[1:] - y_coord[:-1]
            assert (
                steps_y == 1
            ).all(), f"image y coord must step by 1 (pixel index); got {steps_y}"

        # Origin must be small (0 or 0.5 from spatialdata's pixel-center
        # convention), not scaled by pixel_size_um.
        assert x_coord[0] < pixel_size_um, (
            f"image x[0]={x_coord[0]} suggests um leakage "
            f"(pixel_size_um={pixel_size_um})"
        )
        assert y_coord[0] < pixel_size_um, (
            f"image y[0]={y_coord[0]} suggests um leakage "
            f"(pixel_size_um={pixel_size_um})"
        )
