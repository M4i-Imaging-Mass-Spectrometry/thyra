# tests/unit/converters/test_optical_images_collaborator.py
"""``OpticalImages`` on its own, with no converter anywhere near it.

The class owns every piece of optical state a conversion has: the
FlexImaging alignment, the TIC-to-image affine, which file became which
element, and the placeholders waiting for their pixels. That it can be
built from a reader, an output path and a dataset id -- and exercised
end to end from there -- is what "the state has one owner" means; the
converter's own optical tests go through ``convert()`` and so cannot
tell a collaborator from a mixin.

The alignment numbers here are derived from
``TeachingPointAlignment.compute_area_alignment``'s own rule (a region's
raster bounds stretched onto its Area's image bounds), not copied out of
a run, so a change to the affine has to be a deliberate one.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pytest
import tifffile
from spatialdata.transformations import Affine, Identity

from thyra.converters.spatialdata.optical_image import OpticalImages

#: One Area, in optical-photo pixels: (100, 200) to (500, 600).
AREA = {"name": "Area0", "p1": (100, 200), "p2": (500, 600)}
#: A 4 x 3 raster, all of it region 0 -- the index the single Area matches.
N_X, N_Y = 4, 3
IMAGE_FILE = "SlideA_0000.tif"
PIXEL_SIZE = (10.0, 20.0)


class _StubReader:
    """Only what ``OpticalImages`` asks a reader for, and nothing else."""

    def __init__(
        self,
        optical_paths: Optional[List[Path]] = None,
        *,
        areas: Optional[List[Dict[str, Any]]] = None,
        image_file: str = IMAGE_FILE,
    ) -> None:
        self._optical_paths = list(optical_paths or [])
        self.mis_metadata: Dict[str, Any] = {
            "areas": list(areas if areas is not None else [AREA]),
            "ImageFile": image_file,
        }
        self._positions = [
            {"region": 0, "raster_x": x, "raster_y": y}
            for y in range(N_Y)
            for x in range(N_X)
        ]
        self._header = {"first_raster_x": 0, "first_raster_y": 0}

    def get_optical_image_paths(self) -> List[Path]:
        return list(self._optical_paths)

    def get_primary_optical_image_path(self) -> Optional[Path]:
        return None


def _optical(tmp_path: Path, reader: _StubReader, **kwargs: Any) -> OpticalImages:
    kwargs.setdefault("include", True)
    kwargs.setdefault("apply_alignment", True)
    return OpticalImages(
        reader,
        tmp_path / "out.zarr",
        "ds",
        pixel_size_xy=lambda: PIXEL_SIZE,
        **kwargs,
    )


@pytest.fixture
def tiff(tmp_path: Path) -> Path:
    """A real, small, readable RGB TIFF named the way the .mis names it."""
    rng = np.random.default_rng(17)
    path = tmp_path / IMAGE_FILE
    tifffile.imwrite(str(path), rng.integers(0, 256, size=(24, 32, 3), dtype=np.uint8))
    return path


def test_the_alignment_comes_out_of_the_mis_metadata_alone(tmp_path: Path) -> None:
    optical = _optical(tmp_path, _StubReader())
    optical.compute_alignment()

    assert optical.alignment is not None
    assert len(optical.alignment.region_mappings) == 1
    # The .mis names the alignment image; the stem, lowercased, is what
    # every later "is this the primary?" comparison uses.
    assert optical.primary_filename == "slidea_0000"


def test_the_affine_stretches_the_raster_bounds_onto_the_area(tmp_path: Path) -> None:
    optical = _optical(tmp_path, _StubReader())
    optical.compute_alignment()
    optical.build_tic_to_image_affine()

    matrix = optical.tic_to_image
    assert matrix is not None
    assert matrix.shape == (3, 3)
    assert matrix.dtype == np.float64
    # N_X raster columns over (500 - 100) image pixels, N_Y rows over
    # (600 - 200), and the offset puts the first pixel's centre half a
    # raster step inside the Area's own corner.
    scale_x = 400 / N_X
    scale_y = 400 / N_Y
    np.testing.assert_allclose(
        matrix,
        [
            [scale_x, 0.0, 100 + scale_x / 2],
            [0.0, scale_y, 200 + scale_y / 2],
            [0.0, 0.0, 1.0],
        ],
    )


def test_an_area_that_matches_no_region_leaves_no_affine(tmp_path: Path) -> None:
    """The state the pixel polygons are gated against: a result, no matrix."""
    reader = _StubReader()
    reader._positions = [{"region": 3, "raster_x": 0, "raster_y": 0}]
    optical = _optical(tmp_path, reader)
    optical.compute_alignment()
    optical.build_tic_to_image_affine()

    assert optical.alignment is not None
    assert optical.alignment.region_mappings == []
    assert optical.tic_to_image is None


def test_declaring_an_image_records_it_without_reading_its_pixels(
    tmp_path: Path, tiff: Path
) -> None:
    optical = _optical(tmp_path, _StubReader([tiff]))
    optical.compute_alignment()
    images: Dict[str, Any] = {}
    optical.add_images(images)

    assert list(images) == ["ds_optical_highres"]
    assert list(optical.pending) == ["ds_optical_highres"]
    assert optical.sources == {"ds_optical_highres": IMAGE_FILE}
    # The .mis named this file, so it is the alignment image.
    assert optical.alignment_element == "ds_optical_highres"
    assert optical.pending["ds_optical_highres"].source.shape == (3, 24, 32)


def test_the_root_attr_names_the_file_and_the_alignment_element(
    tmp_path: Path, tiff: Path
) -> None:
    """The shape a store's ``optical_images`` root attr has to keep."""
    optical = _optical(tmp_path, _StubReader([tiff]))
    optical.compute_alignment()
    optical.add_images({})

    assert optical.root_attr() == {
        "alignment_element": "ds_optical_highres",
        "elements": {"ds_optical_highres": {"source_file": IMAGE_FILE}},
    }


def test_the_alignment_image_is_an_identity_when_the_msi_moved_to_it(
    tmp_path: Path, tiff: Path
) -> None:
    optical = _optical(tmp_path, _StubReader([tiff]))
    optical.compute_alignment()
    optical.build_tic_to_image_affine()
    optical.add_images({})

    transform = optical.pending["ds_optical_highres"].transformations["global"]
    assert isinstance(transform, Identity)


def test_declining_the_alignment_carries_the_photo_into_micrometers(
    tmp_path: Path, tiff: Path
) -> None:
    """The opt-out route: the photo moves instead, by the affine's inverse."""
    optical = _optical(tmp_path, _StubReader([tiff]), apply_alignment=False)
    optical.compute_alignment()
    optical.build_tic_to_image_affine()
    optical.add_images({})

    transform = optical.pending["ds_optical_highres"].transformations["global"]
    assert isinstance(transform, Affine)
    ps_x, ps_y = PIXEL_SIZE
    expected = np.array(
        [[ps_x, 0.0, 0.0], [0.0, ps_y, 0.0], [0.0, 0.0, 1.0]]
    ) @ np.linalg.inv(optical.tic_to_image)
    np.testing.assert_allclose(
        transform.to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y")),
        expected,
    )
    # Both coordinate systems get the same transform object, as ever.
    assert optical.pending["ds_optical_highres"].transformations["ds"] is transform


def test_excluded_optical_images_leave_every_piece_of_state_empty(
    tmp_path: Path, tiff: Path
) -> None:
    optical = _optical(tmp_path, _StubReader([tiff]), include=False)
    optical.compute_alignment()
    images: Dict[str, Any] = {}
    optical.add_images(images)

    assert images == {}
    assert optical.pending == {}
    assert optical.sources == {}
    assert optical.alignment_element is None
    assert optical.root_attr() is None
