# tests/unit/metadata/test_root_attrs.py
"""The store's root attrs, built from arguments, with no converter anywhere.

The sibling of ``test_uns_assembler.py`` and the point of the same split
(issue #351): the attrs used to be a dozen converter methods reading its
own attributes, so asking what a given input produces meant running a
conversion. Here the builder is handed a stub reader, a stub optical
collaborator and a :class:`RootAttrsContext`.

The end-to-end contract -- that the TIC image and the pixel polygons
resolve to the same space at ``"global"`` -- stays in
``tests/unit/converters/test_coordinate_systems.py``, which drives a real
conversion. What is pinned here is the composition.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np
import pytest

from thyra.core.base_converter import ZSpacingSource
from thyra.metadata.root_attrs import RootAttrsBuilder, RootAttrsContext
from thyra.metadata.types import ComprehensiveMetadata, EssentialMetadata


def _essential(coordinate_offsets=None) -> EssentialMetadata:
    essential = EssentialMetadata(
        dimensions=(3, 3, 1),
        coordinate_bounds=(0.0, 2.0, 0.0, 2.0),
        mass_range=(100.0, 1000.0),
        pixel_size=(10.0, 10.0),
        n_spectra=9,
        total_peaks=90,
        source_path="stub.imzML",
        spectrum_type="centroid spectrum",
    )
    if coordinate_offsets is not None:
        object.__setattr__(essential, "coordinate_offsets", coordinate_offsets)
    return essential


def _comprehensive(**overrides) -> ComprehensiveMetadata:
    fields: Dict[str, Any] = dict(
        essential=_essential(),
        format_specific={"files": ["stub.imzML"]},
        acquisition_params={"polarity": "positive"},
        instrument_info={"model": "stub"},
        raw_metadata={},
    )
    fields.update(overrides)
    return ComprehensiveMetadata(**fields)


class _StubReader:
    """Everything the builder asks a reader for, and nothing else."""

    def __init__(self, *, comprehensive=None, essential=None) -> None:
        self._comprehensive = (
            _comprehensive() if comprehensive is None else comprehensive
        )
        self._essential = essential

    def get_comprehensive_metadata(self) -> ComprehensiveMetadata:
        return self._comprehensive

    def get_essential_metadata(self) -> EssentialMetadata:
        if self._essential is None:
            raise RuntimeError("this stub reports no essential metadata")
        return self._essential


class _StubOptical:
    """The four things the builder reads off the optical collaborator."""

    def __init__(
        self,
        *,
        apply_alignment: bool = False,
        tic_to_image: Optional[np.ndarray] = None,
        alignment_element: Optional[str] = None,
        root: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.apply_alignment = apply_alignment
        self.tic_to_image = tic_to_image
        self.alignment_element = alignment_element
        self._root = root

    @property
    def msi_in_pixel_space(self) -> bool:
        return self.apply_alignment and self.tic_to_image is not None

    def root_attr(self) -> Optional[Dict[str, Any]]:
        return self._root


def _builder(reader=None, **overrides) -> RootAttrsBuilder:
    kwargs: Dict[str, Any] = dict(
        dataset_id="stub",
        pixel_size_detection_info=None,
        handle_3d=False,
        conversion_options={},
    )
    kwargs.update(overrides)
    return RootAttrsBuilder(reader or _StubReader(), **kwargs)


def _context(**overrides) -> RootAttrsContext:
    fields: Dict[str, Any] = dict(
        optical=_StubOptical(),
        pixel_size_xy=(30.0, 50.0),
        dimensions=(3, 3, 1),
        non_empty_pixels=9,
        coordinate_bounds=(0.0, 2.0, 0.0, 2.0),
        is_volume=False,
        z_spacing_um=10.0,
        z_spacing_source=ZSpacingSource.ASSUMED_ISOTROPIC,
    )
    fields.update(overrides)
    return RootAttrsContext(**fields)


class TestTheAttrsFromArgumentsAlone:
    def test_the_pitch_is_the_pair_the_context_carries(self):
        attrs = _builder().build(_context())

        assert attrs["pixel_size_x_um"] == 30.0
        assert attrs["pixel_size_y_um"] == 50.0
        assert attrs["pixel_size_units"] == "micrometers"

    def test_the_dataset_info_counts_the_whole_grid(self):
        attrs = _builder().build(_context(dimensions=(4, 5, 2), non_empty_pixels=7))

        assert attrs["msi_dataset_info"] == {
            "dataset_id": "stub",
            "total_grid_pixels": 40,
            "non_empty_pixels": 7,
            "dimensions_xyz": [4, 5, 2],
        }

    def test_the_comprehensive_sections_are_the_readers(self):
        attrs = _builder().build(_context())

        assert attrs["format_specific_metadata"] == {"files": ["stub.imzML"]}
        assert attrs["acquisition_parameters"] == {"polarity": "positive"}
        assert attrs["instrument_information"] == {"model": "stub"}

    @pytest.mark.parametrize(
        "field,key",
        [
            ("format_specific", "format_specific_metadata"),
            ("acquisition_params", "acquisition_parameters"),
            ("instrument_info", "instrument_information"),
        ],
    )
    def test_a_section_the_reader_left_empty_is_omitted(self, field, key):
        reader = _StubReader(comprehensive=_comprehensive(**{field: {}}))

        assert key not in _builder(reader).build(_context())

    def test_detection_info_is_omitted_when_there_was_none(self):
        assert "pixel_size_detection_info" not in _builder().build(_context())

    def test_detection_info_is_copied_not_aliased(self):
        info = {"method": "automatic"}
        attrs = _builder(pixel_size_detection_info=info).build(_context())

        attrs["pixel_size_detection_info"]["method"] = "mutated"
        assert info == {"method": "automatic"}

    def test_no_optical_images_no_key(self):
        assert "optical_images" not in _builder().build(_context())

    def test_the_caller_can_pass_metadata_it_already_read(self):
        reader = _StubReader(comprehensive=_comprehensive(instrument_info={}))
        already = _comprehensive(instrument_info={"model": "passed in"})

        attrs = _builder(reader).build(_context(), already)

        assert attrs["instrument_information"] == {"model": "passed in"}


class TestTheTwoSectionsRescuedFromTheDeadMethod:
    """Issue #67 item 1: what ``_add_comprehensive_metadata`` never wrote.

    It hung its block on ``SpatialData.metadata``, an attribute no release
    of SpatialData has, so the body never ran. Everything else it built is
    already a root attr or already in ``msi_dataset_info``; these two were
    not written anywhere.
    """

    def test_the_source_bounds_are_recorded(self):
        attrs = _builder().build(_context(coordinate_bounds=(1.0, 4.0, 2.0, 8.0)))

        assert attrs["coordinate_bounds"] == [1.0, 4.0, 2.0, 8.0]

    def test_no_bounds_no_key(self):
        assert "coordinate_bounds" not in _builder().build(
            _context(coordinate_bounds=None)
        )

    def test_the_conversion_options_say_how_it_was_run(self):
        attrs = _builder(handle_3d=True, conversion_options={"use_csc": True}).build(
            _context(z_spacing_um=25.0)
        )

        assert attrs["conversion_options"] == {
            "handle_3d": True,
            "pixel_size_um": 30.0,
            "pixel_size_y_um": 50.0,
            "z_spacing_um": 25.0,
            "z_spacing_source": ZSpacingSource.ASSUMED_ISOTROPIC.value,
            "dataset_id": "stub",
            "use_csc": True,
        }

    def test_handle_3d_records_what_was_asked_not_what_was_written(self):
        # A single-slice source converted with handle_3d=True takes the 2D
        # branch, so is_volume is False while the option stays True.
        attrs = _builder(handle_3d=True).build(_context(is_volume=False))

        assert attrs["conversion_options"]["handle_3d"] is True
        assert "z_spacing_um" not in attrs["coordinate_systems"]["global"]


class TestTheCoordinateContract:
    def test_without_alignment_global_is_micrometers(self):
        global_cs = _builder().build(_context())["coordinate_systems"]["global"]

        assert global_cs["unit"] == "micrometer"
        assert global_cs["pixel_size_um_x"] == 30.0
        assert global_cs["pixel_size_um_y"] == 50.0
        assert global_cs["reference_element"] is None
        assert global_cs["raster_to_global_affine"] == [
            [30.0, 0.0, 0.0],
            [0.0, 50.0, 0.0],
            [0.0, 0.0, 1.0],
        ]

    def test_with_alignment_global_is_optical_pixels(self):
        affine = np.array([[2.0, 0.0, 5.0], [0.0, 2.0, 7.0], [0.0, 0.0, 1.0]])
        optical = _StubOptical(
            apply_alignment=True,
            tic_to_image=affine,
            alignment_element="stub_optical_align",
        )

        global_cs = _builder().build(_context(optical=optical))["coordinate_systems"][
            "global"
        ]

        assert global_cs["unit"] == "pixel"
        assert global_cs["pixel_size_um_x"] is None
        assert global_cs["pixel_size_um_y"] is None
        assert global_cs["reference_element"] == "stub_optical_align"
        assert global_cs["raster_to_global_affine"] == affine.tolist()

    def test_an_available_alignment_that_was_not_applied_stays_micrometers(self):
        # Issue #288: the matrix exists whenever FlexImaging data does,
        # because the opt-out path needs its inverse. It says nothing
        # about whether the raster was moved.
        optical = _StubOptical(
            apply_alignment=False,
            tic_to_image=np.eye(3),
            alignment_element="stub_optical_align",
        )

        global_cs = _builder().build(_context(optical=optical))["coordinate_systems"][
            "global"
        ]

        assert global_cs["unit"] == "micrometer"
        assert global_cs["reference_element"] is None

    def test_a_volume_spaces_its_slices(self):
        global_cs = _builder().build(_context(is_volume=True, z_spacing_um=25.0))[
            "coordinate_systems"
        ]["global"]

        assert global_cs["z_spacing_um"] == 25.0
        assert global_cs["z_spacing_source"] == ZSpacingSource.ASSUMED_ISOTROPIC.value

    def test_a_plane_says_nothing_about_z(self):
        global_cs = _builder().build(_context())["coordinate_systems"]["global"]

        assert "z_spacing_um" not in global_cs
        assert "z_spacing_source" not in global_cs

    def test_the_convention_version_does_not_move_for_additive_keys(self):
        global_cs = _builder().build(_context())["coordinate_systems"]["global"]

        assert global_cs["convention_version"] == 1


class TestTheSourceOffsets:
    def test_they_are_recorded_with_their_physical_equivalent(self):
        reader = _StubReader(essential=_essential(coordinate_offsets=(3, 4, 0)))

        global_cs = _builder(reader).build(_context())["coordinate_systems"]["global"]

        assert global_cs["coordinate_offsets_px"] == [3, 4, 0]
        assert global_cs["stage_offset_um"] == [3 * 30.0, 4 * 50.0]

    def test_in_pixel_space_the_micrometre_offset_is_not_written(self):
        # It would be read as optical pixels, which it is not.
        reader = _StubReader(essential=_essential(coordinate_offsets=(3, 4, 0)))
        optical = _StubOptical(apply_alignment=True, tic_to_image=np.eye(3))

        global_cs = _builder(reader).build(_context(optical=optical))[
            "coordinate_systems"
        ]["global"]

        assert global_cs["coordinate_offsets_px"] == [3, 4, 0]
        assert "stage_offset_um" not in global_cs

    def test_a_reader_that_cannot_say_costs_only_the_offsets(self):
        attrs = _builder(_StubReader(essential=None)).build(_context())

        assert "coordinate_offsets_px" not in attrs["coordinate_systems"]["global"]
        assert attrs["msi_dataset_info"]["dataset_id"] == "stub"
