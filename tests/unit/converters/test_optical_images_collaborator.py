# tests/unit/converters/test_optical_images_collaborator.py
"""``OpticalImages`` on its own, with no converter anywhere near it.

The class owns every piece of optical state a conversion has: the
FlexImaging alignment, the TIC-to-image affine, which file became which
element, and the placeholders waiting for their pixels. That it can be
built from a reader, an accessor for the output path and a dataset id
-- and exercised end to end from there -- is what "the state has one
owner" means; the converter's own optical tests go through ``convert()``
and so cannot tell a collaborator from a mixin.

The output path is an accessor rather than a path because the pixels are
streamed long after this object is built, and the converter's
``output_path`` is not fixed at construction. See the last three tests.

The alignment numbers here are derived from
``TeachingPointAlignment.compute_area_alignment``'s own rule (a region's
raster bounds stretched onto its Area's image bounds), not copied out of
a run, so a change to the affine has to be a deliberate one.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional

import numpy as np
import pytest
import tifffile
import zarr
from spatialdata.transformations import Affine, Identity

from thyra.converters.spatialdata.optical_image import (
    OpticalImages,
    StreamedOpticalImage,
)

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
    kwargs.setdefault("output_path", lambda: tmp_path / "out.zarr")
    return OpticalImages(
        reader,
        kwargs.pop("output_path"),
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


# --- The store is read when it is used, not when this object is built -------
#
# Everything above stops at the declaration. The three reads of the output
# path all happen after it: the pixels stream once the SpatialData write
# that carried the placeholders has returned, and a dropped image is
# discarded and unrecorded after that again.
#
# ``convert_msi`` shortens the output path for Windows
# (``prepare_zarr_output_path``) *before* it builds the converter, so the
# production route never moves it underneath this object. A caller that
# stands a converter up itself has to shorten afterwards, and does --
# ``test_internal_failures_are_not_warnings`` does it three times, and the
# CLI's ``_quarantine_partial_output`` records that ``convert_msi``
# "resolves and may extend the output path". This object used to capture
# the path at construction and so kept the pre-move spelling, while the
# converter wrote the store, made its scratch directories and consolidated
# from a path it read at call time.


def _moving_store(tmp_path: Path, tiff: Path):
    """A collaborator with its images declared, and the converter it reads."""
    converter = SimpleNamespace(output_path=tmp_path / "before.zarr")
    optical = _optical(
        tmp_path,
        _StubReader([tiff]),
        output_path=lambda: converter.output_path,
    )
    optical.compute_alignment()
    optical.add_images({})
    assert list(optical.pending) == ["ds_optical_highres"]
    return converter, optical


def test_the_pixels_stream_into_the_store_being_written_now(
    tmp_path: Path, tiff: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The store the converter holds when the pixels stream, not before."""
    converter, optical = _moving_store(tmp_path, tiff)

    handed: List[Path] = []
    monkeypatch.setattr(
        StreamedOpticalImage,
        "stream_pixels",
        lambda self, store_path: handed.append(Path(store_path)),
    )

    converter.output_path = tmp_path / "after.zarr"
    assert optical.stream_pending_pixels() == 1

    assert handed == [tmp_path / "after.zarr"]


def test_a_dropped_image_is_discarded_from_that_same_store(
    tmp_path: Path, tiff: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The tolerant branch reads the path once, and it is the current one.

    An image whose pixels will not decode is taken back out of the store,
    which is only ever the store the placeholder was written into.
    """
    converter, optical = _moving_store(tmp_path, tiff)

    discarded: List[Path] = []

    def _refuse(self, store_path: Path) -> None:
        raise RuntimeError("strips are truncated")

    monkeypatch.setattr(StreamedOpticalImage, "stream_pixels", _refuse)
    monkeypatch.setattr(
        StreamedOpticalImage,
        "discard",
        lambda self, store_path: discarded.append(Path(store_path)),
    )

    converter.output_path = tmp_path / "after.zarr"
    assert optical.stream_pending_pixels() == 0

    assert discarded == [tmp_path / "after.zarr"]


def test_unrecording_a_dropped_image_corrects_the_store_being_written_now(
    tmp_path: Path, tiff: Path
) -> None:
    """The third read, and the one whose failure is silent.

    ``forget_image`` swallows everything it cannot do into a warning,
    because it runs on a store the conversion has just written. Pointed at
    a store that was never written it would warn and return, leaving the
    real store's root attrs naming an element with no pixels -- exactly
    the state the method exists to prevent.
    """
    converter, optical = _moving_store(tmp_path, tiff)
    declared = optical.root_attr()
    assert declared is not None

    # Only the store the conversion is writing now exists on disk, with the
    # attrs the SpatialData write would have carried into it.
    converter.output_path = tmp_path / "after.zarr"
    root = zarr.open_group(str(converter.output_path), mode="w")
    root.attrs["optical_images"] = declared

    optical.forget_image("ds_optical_highres")

    reread = zarr.open_group(str(converter.output_path), mode="r")
    assert "optical_images" not in reread.attrs.asdict()
    assert optical.root_attr() is None
    assert not (tmp_path / "before.zarr").exists()
