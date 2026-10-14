# tests/unit/converters/test_acquisition_coordinate_system.py
"""The ``"acquisition"`` coordinate system puts the raster where the source did.

Readers re-zero the raster, so ``"global"`` starts every store at 0 and two
stores cut from one slide land on top of each other. A store in micrometres
also carries ``"acquisition"``: ``"global"`` followed by a translation to
the source's origin, the index offset times the pitch, or Waters' fitted
stage origin (D48). These assertions are on stored stores, read back with
SpatialData, because the defect was that nothing in the store placed them.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pytest
import spatialdata
from spatialdata.transformations import Identity, Scale, Translation, get_transformation

from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
from thyra.convert import dataset_id_problem
from thyra.converters.spatialdata.acquisition_frame import (
    ACQUISITION_COORDINATE_SYSTEM,
    ORIGIN_FROM_RASTER_INDEX,
    ORIGIN_FROM_STAGE_POSITION,
    AcquisitionOrigin,
    acquisition_origin,
    element_transformations,
)
from thyra.converters.spatialdata.streaming_converter import (
    StreamingSpatialDataConverter,
)
from thyra.core.base_converter import PixelSizeSource
from thyra.errors import ConversionRefused

N_X, N_Y = 3, 2
PITCH_X, PITCH_Y = 30.0, 50.0
ACQ = ACQUISITION_COORDINATE_SYSTEM


class _OffsetReader(MockMSIReader):
    """A mock source that reports where its raster started."""

    def __init__(
        self,
        config: MockMSIConfig,
        offsets: Optional[Tuple[int, int, int]] = None,
        stage_origin_um: Optional[Tuple[float, float]] = None,
    ) -> None:
        super().__init__(config)
        self._offsets = offsets
        self._stage_origin = stage_origin_um

    def get_essential_metadata(self):
        essential = super().get_essential_metadata()
        object.__setattr__(essential, "pixel_size", (PITCH_X, PITCH_Y))
        object.__setattr__(essential, "coordinate_offsets", self._offsets)
        object.__setattr__(essential, "stage_origin_um", self._stage_origin)
        return essential


def _config(n_z: int = 1) -> MockMSIConfig:
    return MockMSIConfig(
        n_x=N_X,
        n_y=N_Y,
        n_z=n_z,
        n_mz_bins=64,
        peaks_per_spectrum=(3, 5),
        seed=5,
    )


def _convert(output: Path, reader: MockMSIReader, **kwargs) -> spatialdata.SpatialData:
    converter = StreamingSpatialDataConverter(
        reader=reader,
        output_path=output,
        dataset_id="m",
        pixel_size_um=PITCH_X,
        pixel_size_source=PixelSizeSource.AUTO_DETECTED,
        pixel_size_detection_info={
            "method": "automatic",
            "detected_x_um": PITCH_X,
            "detected_y_um": PITCH_Y,
            "source_format": "mock",
            "detection_successful": True,
        },
        include_optical=False,
        **kwargs,
    )
    assert converter.convert() is True
    return spatialdata.SpatialData.read(str(output))


def _root_attrs(store: Path) -> dict:
    doc = json.loads((store / "zarr.json").read_text())
    return doc.get("attributes", {}) or {}


def _matrix(element, system: str) -> np.ndarray:
    return get_transformation(element, to_coordinate_system=system).to_affine_matrix(
        input_axes=("x", "y"), output_axes=("x", "y")
    )


def _shapes_bbox(element, system: str) -> np.ndarray:
    """The polygons' bounds in ``system``, as ``[xmin, ymin, xmax, ymax]``."""
    xmin, ymin, xmax, ymax = element.total_bounds
    corners = (
        np.array([[xmin, ymin, 1.0], [xmax, ymax, 1.0]]) @ _matrix(element, system).T
    )
    return np.array([corners[0, 0], corners[0, 1], corners[1, 0], corners[1, 1]])


class TestTheOrigin:
    """Where index 0 sits, from what the reader reported."""

    def test_index_offsets_times_the_pitch_of_each_axis(self):
        origin = acquisition_origin((4, 7, 2), None, (PITCH_X, PITCH_Y), 5.0)
        assert origin == AcquisitionOrigin(120.0, 350.0, 10.0, ORIGIN_FROM_RASTER_INDEX)

    def test_z_is_zero_without_a_slice_pitch(self):
        origin = acquisition_origin((1, 1, 3), None, (10.0, 10.0), None)
        assert origin is not None and origin.z_um == 0.0

    def test_a_stage_origin_is_taken_as_it_is(self):
        """Waters: a measured position, negative and fractional."""
        origin = acquisition_origin(None, (-11400.25, 3050.5), (PITCH_X, PITCH_Y), None)
        assert origin == AcquisitionOrigin(
            -11400.25, 3050.5, 0.0, ORIGIN_FROM_STAGE_POSITION
        )

    def test_nothing_reported_is_no_origin(self):
        assert acquisition_origin(None, None, (PITCH_X, PITCH_Y), None) is None

    def test_zero_offsets_are_an_origin_at_zero(self):
        """PHI states (0, 0, 0): a known origin, not an unknown one."""
        origin = acquisition_origin((0, 0, 0), None, (PITCH_X, PITCH_Y), None)
        assert origin == AcquisitionOrigin(0.0, 0.0, 0.0, ORIGIN_FROM_RASTER_INDEX)


class TestElementTransformations:
    def test_no_origin_writes_the_two_systems_it_always_did(self):
        scale = Scale([2.0, 2.0], axes=("x", "y"))
        assert element_transformations("ds", scale, None) == {
            "ds": scale,
            "global": scale,
        }

    def test_an_identity_becomes_the_translation_alone(self):
        origin = AcquisitionOrigin(30.0, 50.0, 7.0, ORIGIN_FROM_RASTER_INDEX)
        transform = element_transformations("ds", Identity(), origin)[ACQ]
        assert isinstance(transform, Translation)
        assert transform.axes == ("x", "y")
        np.testing.assert_allclose(transform.translation, [30.0, 50.0])

    def test_the_translation_comes_after_the_transform(self):
        origin = AcquisitionOrigin(30.0, 50.0, 7.0, ORIGIN_FROM_RASTER_INDEX)
        scale = Scale([2.0, 3.0], axes=("x", "y"))
        matrix = element_transformations("ds", scale, origin)[ACQ].to_affine_matrix(
            input_axes=("x", "y"), output_axes=("x", "y")
        )
        np.testing.assert_allclose(matrix, [[2, 0, 30], [0, 3, 50], [0, 0, 1]])


class TestAStoreWithAnOffset:
    """Bruker's imaging-area minimum, a --region's minimum, rapifleX."""

    OFFSETS = (4, 7, 0)
    ORIGIN = (4 * PITCH_X, 7 * PITCH_Y)

    @pytest.fixture(scope="class")
    @classmethod
    def store(cls, tmp_path_factory) -> Path:
        out = tmp_path_factory.mktemp("acq") / "out.zarr"
        _convert(out, _OffsetReader(_config(), offsets=cls.OFFSETS))
        return out

    @pytest.fixture(scope="class")
    @classmethod
    def sdata(cls, store):
        return spatialdata.SpatialData.read(str(store))

    def test_every_element_is_in_the_acquisition_system(self, sdata):
        for element in [*sdata.images.values(), *sdata.shapes.values()]:
            systems = get_transformation(element, get_all=True)
            assert set(systems) == {"m", "global", ACQ}

    def test_the_tic_image_moves_by_the_origin(self, sdata):
        image = sdata.images["m_z0_tic"]
        delta = _matrix(image, ACQ) - _matrix(image, "global")
        np.testing.assert_allclose(delta[:2, 2], self.ORIGIN)
        np.testing.assert_allclose(delta[:2, :2], 0.0)

    def test_the_polygons_move_by_the_origin(self, sdata):
        shapes = sdata.shapes["m_z0_pixels"]
        np.testing.assert_allclose(
            _shapes_bbox(shapes, ACQ) - _shapes_bbox(shapes, "global"),
            [*self.ORIGIN, *self.ORIGIN],
        )

    def test_image_and_polygons_still_agree(self, sdata):
        """The D32 cell-centre contract, now in the acquisition system too."""
        image = sdata.images["m_z0_tic"]
        corners = np.array([[0.0, 0.0, 1.0], [N_X, N_Y, 1.0]]) @ _matrix(image, ACQ).T
        image_bbox = [corners[0, 0], corners[0, 1], corners[1, 0], corners[1, 1]]
        np.testing.assert_allclose(
            _shapes_bbox(sdata.shapes["m_z0_pixels"], ACQ), image_bbox
        )

    def test_the_global_system_is_unchanged(self, sdata):
        image = sdata.images["m_z0_tic"]
        np.testing.assert_allclose(
            _matrix(image, "global"),
            [
                [PITCH_X, 0.0, -0.5 * PITCH_X],
                [0.0, PITCH_Y, -0.5 * PITCH_Y],
                [0.0, 0.0, 1.0],
            ],
        )

    def test_the_root_attrs_describe_it(self, store):
        systems = _root_attrs(store)["coordinate_systems"]
        acq = systems[ACQ]
        assert acq["unit"] == "micrometer"
        assert acq["origin_um"] == list(self.ORIGIN)
        assert acq["origin_source"] == ORIGIN_FROM_RASTER_INDEX
        assert acq["convention_version"] == 1
        assert (acq["pixel_size_um_x"], acq["pixel_size_um_y"]) == (PITCH_X, PITCH_Y)
        np.testing.assert_allclose(
            acq["raster_to_acquisition_affine"],
            [
                [PITCH_X, 0.0, -0.5 * PITCH_X + self.ORIGIN[0]],
                [0.0, PITCH_Y, -0.5 * PITCH_Y + self.ORIGIN[1]],
                [0.0, 0.0, 1.0],
            ],
        )

    def test_the_attr_affine_is_the_tic_image_transform(self, store, sdata):
        acq = _root_attrs(store)["coordinate_systems"][ACQ]
        np.testing.assert_allclose(
            acq["raster_to_acquisition_affine"],
            _matrix(sdata.images["m_z0_tic"], ACQ),
        )

    def test_the_stage_offset_is_the_same_origin(self, store):
        global_cs = _root_attrs(store)["coordinate_systems"]["global"]
        assert global_cs["coordinate_offsets_px"] == list(self.OFFSETS)
        assert global_cs["stage_offset_um"] == list(self.ORIGIN)


class TestTwoStoresOfOneSlide:
    """Two --region runs of one .d: overlapping in global, side by side here."""

    @pytest.fixture(scope="class")
    @classmethod
    def pair(cls, tmp_path_factory):
        root = tmp_path_factory.mktemp("pair")
        left = _convert(root / "a.zarr", _OffsetReader(_config(), offsets=(10, 5, 0)))
        right = _convert(
            root / "b.zarr", _OffsetReader(_config(), offsets=(10 + N_X, 5, 0))
        )
        return left.shapes["m_z0_pixels"], right.shapes["m_z0_pixels"]

    def test_they_overlap_in_global(self, pair):
        left, right = pair
        np.testing.assert_allclose(
            _shapes_bbox(left, "global"), _shapes_bbox(right, "global")
        )

    def test_they_sit_side_by_side_in_acquisition(self, pair):
        left, right = pair
        a = _shapes_bbox(left, ACQ)
        b = _shapes_bbox(right, ACQ)
        assert b[0] == pytest.approx(a[2])  # right's left edge is left's right
        np.testing.assert_allclose([a[1], a[3]], [b[1], b[3]])


class TestAStageOrigin:
    """Waters: the fitted stage origin, not an index offset."""

    ORIGIN = (-11400.25, 3050.5)

    @pytest.fixture(scope="class")
    @classmethod
    def store(cls, tmp_path_factory) -> Path:
        out = tmp_path_factory.mktemp("waters") / "out.zarr"
        _convert(out, _OffsetReader(_config(), stage_origin_um=cls.ORIGIN))
        return out

    def test_the_elements_move_by_the_stage_origin(self, store):
        sdata = spatialdata.SpatialData.read(str(store))
        image = sdata.images["m_z0_tic"]
        delta = _matrix(image, ACQ) - _matrix(image, "global")
        np.testing.assert_allclose(delta[:2, 2], self.ORIGIN)

    def test_the_root_attrs_say_where_it_came_from(self, store):
        systems = _root_attrs(store)["coordinate_systems"]
        assert systems[ACQ]["origin_um"] == list(self.ORIGIN)
        assert systems[ACQ]["origin_source"] == ORIGIN_FROM_STAGE_POSITION
        assert systems["global"]["stage_offset_um"] == list(self.ORIGIN)
        # Not an integer index, so not written as one.
        assert "coordinate_offsets_px" not in systems["global"]


class TestNoOriginReported:
    @pytest.fixture(scope="class")
    @classmethod
    def store(cls, tmp_path_factory) -> Path:
        out = tmp_path_factory.mktemp("none") / "out.zarr"
        _convert(out, _OffsetReader(_config()))
        return out

    def test_no_element_gets_the_system(self, store):
        sdata = spatialdata.SpatialData.read(str(store))
        for element in [*sdata.images.values(), *sdata.shapes.values()]:
            assert set(get_transformation(element, get_all=True)) == {"m", "global"}

    def test_the_root_attrs_have_no_entry(self, store):
        systems = _root_attrs(store)["coordinate_systems"]
        assert set(systems) == {"global"}
        assert "stage_offset_um" not in systems["global"]


class TestAVolume:
    """z moves by the z offset times the slice pitch; shapes stay flat."""

    Z_SPACING = 4.0
    OFFSETS = (2, 3, 1)

    @pytest.fixture(scope="class")
    @classmethod
    def store(cls, tmp_path_factory) -> Path:
        out = tmp_path_factory.mktemp("vol") / "out.zarr"
        _convert(
            out,
            _OffsetReader(_config(n_z=2), offsets=cls.OFFSETS),
            handle_3d=True,
            z_spacing_um=cls.Z_SPACING,
        )
        return out

    def test_the_volume_moves_in_z_too(self, store):
        sdata = spatialdata.SpatialData.read(str(store))
        (name,) = [k for k in sdata.images if k.endswith("_tic")]
        image = sdata.images[name]
        axes = ("x", "y", "z")

        def matrix(system):
            return get_transformation(
                image, to_coordinate_system=system
            ).to_affine_matrix(input_axes=axes, output_axes=axes)

        delta = matrix(ACQ) - matrix("global")
        np.testing.assert_allclose(
            delta[:3, 3], [2 * PITCH_X, 3 * PITCH_Y, 1 * self.Z_SPACING]
        )

    def test_the_root_attrs_carry_the_z_origin(self, store):
        acq = _root_attrs(store)["coordinate_systems"][ACQ]
        assert acq["origin_um"] == [2 * PITCH_X, 3 * PITCH_Y, self.Z_SPACING]
        assert acq["z_spacing_um"] == self.Z_SPACING


class TestTheNameIsReserved:
    def test_the_front_door_refuses_it(self):
        assert "coordinate system" in dataset_id_problem(ACQ)

    def test_the_converter_refuses_it(self, tmp_path):
        with pytest.raises(ConversionRefused, match="coordinate system"):
            StreamingSpatialDataConverter(
                reader=MockMSIReader(_config()),
                output_path=tmp_path / "out.zarr",
                dataset_id=ACQ,
                include_optical=False,
            )

    def test_global_stays_allowed(self):
        assert dataset_id_problem("global") is None
