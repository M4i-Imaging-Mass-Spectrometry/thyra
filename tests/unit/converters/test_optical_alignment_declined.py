# tests/unit/converters/test_optical_alignment_declined.py
"""What the store says about ``"global"`` when the alignment is declined.

``apply_optical_alignment=False`` is the option a downstream tool takes
when it computes its own MSI-to-target registration: Thyra's alignment
is not the canonical one, so pre-applying it would only have to be
undone. Every MSI element then stays in micrometers.

The alignment *matrix* is still built, though, and deliberately so --
the opt-out path needs its inverse to carry the optical photo into the
same micrometer frame. So the matrix existing says an alignment is
available, never that it was applied.

``_build_coordinate_systems_attr`` read only the matrix, and so declared
``unit="pixel"``, ``pixel_size_um_x/y = None``, a ``reference_element``
and the optical affine on a store whose table, TIC image and polygons
were all in micrometers -- and it suppressed ``stage_offset_um``, which
is written only in the micrometer variant. The one artefact a consumer
is told to read to find out what ``"global"`` means was the one artefact
that had it wrong (issue #288).

The assertion that matters most is the last one: the affine in the attr
and the transform the TIC image actually carries must agree. That is the
invariant the attr exists to publish, and it held in neither direction
before the fix.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import spatialdata
from spatialdata.transformations import get_transformation

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)

N_X, N_Y = 3, 2
PITCH = 25.0
KEY = "m_z0"
OPTICAL = "slide_overview.tif"

#: A stand-in for what ``_build_tic_to_image_affine`` computes from real
#: FlexImaging landmarks: raster indices to optical-photo pixels.
TIC_TO_IMAGE = np.array(
    [[4.0, 0.0, 10.0], [0.0, 4.0, 20.0], [0.0, 0.0, 1.0]], dtype=np.float64
)


class _AlignedConverter(StreamingSpatialDataConverter):
    """A converter that has alignment data, without needing a ``.mis``.

    Only the affine construction is stubbed. The predicate under test,
    the attr builder and the TIC image's transform are all the real
    ones, which is the whole point -- the defect was in how those three
    read the state, not in how the state is computed.
    """

    def _build_tic_to_image_affine(self) -> None:
        self._tic_to_image_matrix = TIC_TO_IMAGE.copy()
        self._primary_optical_filename = OPTICAL


def _convert(output: Path, *, apply_optical_alignment: bool) -> Path:
    converter = _AlignedConverter(
        reader=MockMSIReader(
            MockMSIConfig(
                n_x=N_X,
                n_y=N_Y,
                n_z=1,
                n_mz_bins=64,
                peaks_per_spectrum=(3, 5),
                seed=7,
                pixel_size_um=PITCH,
            )
        ),
        output_path=output,
        dataset_id="m",
        pixel_size_um=PITCH,
        include_optical=False,
        apply_optical_alignment=apply_optical_alignment,
    )
    assert converter.convert() is True
    return output


def _global_cs(store: Path) -> dict:
    doc = json.loads((store / "zarr.json").read_text())
    return doc["attributes"]["coordinate_systems"]["global"]


@pytest.fixture(scope="module")
def declined(tmp_path_factory) -> Path:
    return _convert(
        tmp_path_factory.mktemp("declined") / "out.zarr",
        apply_optical_alignment=False,
    )


@pytest.fixture(scope="module")
def applied(tmp_path_factory) -> Path:
    return _convert(
        tmp_path_factory.mktemp("applied") / "out.zarr",
        apply_optical_alignment=True,
    )


class TestTheAlignmentWasDeclined:
    """Mode A, even though the alignment data is there."""

    def test_the_unit_is_micrometer(self, declined):
        assert _global_cs(declined)["unit"] == "micrometer"

    def test_the_pitch_is_stated_rather_than_nulled(self, declined):
        cs = _global_cs(declined)
        assert cs["pixel_size_um_x"] == pytest.approx(PITCH)
        assert cs["pixel_size_um_y"] == pytest.approx(PITCH)

    def test_no_optical_image_is_named_as_the_reference(self, declined):
        assert _global_cs(declined)["reference_element"] is None

    def test_the_raster_affine_is_the_pitch_not_the_alignment(self, declined):
        np.testing.assert_allclose(
            _global_cs(declined)["raster_to_global_affine"],
            [[PITCH, 0.0, 0.0], [0.0, PITCH, 0.0], [0.0, 0.0, 1.0]],
        )


class TestTheAlignmentWasApplied:
    """Mode B is unchanged; the fix must not collapse the two."""

    def test_the_unit_is_pixel(self, applied):
        assert _global_cs(applied)["unit"] == "pixel"

    def test_the_pitch_is_null(self, applied):
        cs = _global_cs(applied)
        assert cs["pixel_size_um_x"] is None
        assert cs["pixel_size_um_y"] is None

    def test_the_primary_optical_image_is_the_reference(self, applied):
        assert _global_cs(applied)["reference_element"] == OPTICAL

    def test_the_raster_affine_is_the_alignment(self, applied):
        np.testing.assert_allclose(
            _global_cs(applied)["raster_to_global_affine"], TIC_TO_IMAGE
        )


class TestTheAttrAgreesWithTheElement:
    """``raster_to_global_affine`` is a duplicate of the TIC image's own
    transform, published for consumers that read only attrs. A duplicate
    that disagrees with its original is worse than no duplicate."""

    @pytest.mark.parametrize("store_name", ["declined", "applied"])
    def test_the_published_affine_is_the_transform_the_image_carries(
        self, store_name, request
    ):
        store = request.getfixturevalue(store_name)
        sdata = spatialdata.SpatialData.read(str(store))
        matrix = get_transformation(
            sdata.images[f"{KEY}_tic"], to_coordinate_system="global"
        ).to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))
        np.testing.assert_allclose(
            _global_cs(store)["raster_to_global_affine"], matrix, atol=1e-9
        )


class TestAlignmentDataThatYieldsNoMatrix:
    """An alignment result whose ``region_mappings`` is empty.

    ``TeachingPointAlignment.compute_area_alignment`` skips every area
    with no matching region, so a real ``.mis`` whose areas match no
    region in the data produces an ``AreaAlignmentResult`` with an empty
    mapping list. ``_build_tic_to_image_affine`` then returns early and
    no matrix exists.

    The pixel polygons used to gate on ``_alignment_result`` alone, so in
    that state they took the alignment branch while the attr and the TIC
    image took the micrometer one. ``transform_point`` iterates
    ``region_mappings`` and returns ``None`` for every position, so the
    shapes element came out **empty** -- a store with a table, an image
    and no polygons at all.
    """

    def _convert(self, output: Path) -> Path:
        from thyra.alignment.teaching_points import AreaAlignmentResult

        class _NoMappings(StreamingSpatialDataConverter):
            def _compute_optical_alignment(self) -> None:
                self._alignment_result = AreaAlignmentResult(
                    region_mappings=[], first_raster_x=0, first_raster_y=0
                )

            def _build_tic_to_image_affine(self) -> None:
                # What the real one does with no region mappings.
                return

        converter = _NoMappings(
            reader=MockMSIReader(
                MockMSIConfig(
                    n_x=N_X,
                    n_y=N_Y,
                    n_z=1,
                    n_mz_bins=64,
                    peaks_per_spectrum=(3, 5),
                    seed=11,
                    pixel_size_um=PITCH,
                )
            ),
            output_path=output,
            dataset_id="m",
            pixel_size_um=PITCH,
            include_optical=False,
        )
        assert converter.convert() is True
        return output

    @pytest.fixture(scope="class")
    def store(self, tmp_path_factory) -> Path:
        return self._convert(tmp_path_factory.mktemp("nomap") / "out.zarr")

    def test_the_store_declares_micrometers(self, store):
        assert _global_cs(store)["unit"] == "micrometer"

    def test_the_polygons_are_there(self, store):
        """The assertion the old gate failed: not an empty shapes element."""
        sdata = spatialdata.SpatialData.read(str(store))
        shapes = next(iter(sdata.shapes.values()))
        assert len(shapes) == N_X * N_Y

    def test_the_polygons_are_in_micrometers_like_everything_else(self, store):
        sdata = spatialdata.SpatialData.read(str(store))
        shapes = next(iter(sdata.shapes.values()))
        xmin, ymin, xmax, ymax = shapes.total_bounds
        # The raster is N_X x N_Y pixels of PITCH micrometers.
        assert xmax - xmin == pytest.approx(N_X * PITCH)
        assert ymax - ymin == pytest.approx(N_Y * PITCH)
