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
``TeachingPointAlignment.compute_area_alignment``'s own rules, not copied
out of a run, so a change to the affine has to be a deliberate one. The
stub states no teaching points, so most tests get the approximate fallback
(a region's raster bounds stretched onto its Area's image bounds); one
test gives it teaching points and gets the lattice (D31).
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

from thyra.alignment.teaching_points import (
    AreaAlignmentResult,
    LatticeFit,
    RegionMapping,
)
from thyra.converters.spatialdata.optical_image import (
    OpticalImages,
    StreamedOpticalImage,
)
from thyra.core.base_reader import OpticalImageLabel

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
    # Without teaching points this is the approximate fallback: N_X raster
    # columns over (500 - 100) image pixels, N_Y rows over (600 - 200). The
    # TIC cell i spans [i, i + 1), so coordinate 0 is the Area's corner and
    # the first pixel's centre, at 0.5, is half a raster step inside it.
    scale_x = 400 / N_X
    scale_y = 400 / N_Y
    np.testing.assert_allclose(
        matrix,
        [
            [scale_x, 0.0, 100.0],
            [0.0, scale_y, 200.0],
            [0.0, 0.0, 1.0],
        ],
    )


def test_with_teaching_points_the_affine_is_flexImagings_lattice(
    tmp_path: Path,
) -> None:
    """Teaching points and a reference point give the lattice, not the stretch.

    10 um per photo pixel, stage y up, a 20 um raster: one raster step is 2
    photo pixels whatever the Area's size. The Area here is drawn larger
    than the spots, which the bounding-box stretch would have turned into
    bigger pixels.
    """
    reader = _StubReader(areas=[{"name": "Area0", "p1": (99, 199), "p2": (109, 207)}])
    reader.mis_metadata.update(
        teaching_points=[
            {"image": [0, 0], "stage": [0, 0]},
            {"image": [100, 0], "stage": [1000, 0]},
            {"image": [0, 100], "stage": [0, -1000]},
        ],
        reference_point=[0.0, 0.0],
        raster=[20, 20],
    )
    optical = _optical(tmp_path, reader)
    optical.compute_alignment()
    optical.build_tic_to_image_affine()

    assert optical.alignment is not None
    lattice = optical.alignment.lattice
    assert lattice is not None
    assert lattice.spots_outside == 0
    matrix = optical.tic_to_image
    assert matrix is not None
    np.testing.assert_allclose(matrix, lattice.cell_to_image)
    # One step is 2 photo pixels, and cell 0's centre (0.5, 0.5) is a node.
    np.testing.assert_allclose(matrix[:2, :2], [[2.0, 0.0], [0.0, 2.0]], atol=1e-9)
    centre = matrix @ np.array([0.5, 0.5, 1.0])
    off_node = np.abs((centre[:2] + 1.0) % 2.0 - 1.0)
    np.testing.assert_allclose(off_node, [0.0, 0.0], atol=1e-9)
    corners = np.array(optical.alignment.cell_corners(0, 0))
    np.testing.assert_allclose(corners.mean(axis=0), centre[:2])


class _StatingReader(_StubReader):
    """A source that states its own registration, as an mzPeak archive can."""

    def __init__(self, alignment: AreaAlignmentResult) -> None:
        super().__init__()
        self._alignment = alignment

    def get_image_alignment(self) -> Optional[AreaAlignmentResult]:
        return self._alignment


def test_an_alignment_the_source_states_is_taken_whole(tmp_path: Path) -> None:
    """Before the .mis route, and the TIC is placed by its lattice."""
    cell_to_image = np.array([[50.0, 1.0, 7.0], [-1.0, 50.0, 9.0], [0.0, 0.0, 1.0]])
    stated = AreaAlignmentResult(
        region_mappings=[RegionMapping(0, "0", 5, 8, 2, 4, 0, 210, 0, 160)],
        first_raster_x=5,
        first_raster_y=2,
        pos_to_region={(5 + x, 2 + y): 0 for x in range(4) for y in range(3)},
        lattice=LatticeFit(cell_to_image, (5, 2), (1, 1), (10.0, 10.0), 0, 0, "stated"),
    )
    optical = _optical(tmp_path, _StatingReader(stated))
    optical.compute_alignment()
    optical.build_tic_to_image_affine()

    assert optical.alignment is stated
    assert optical.msi_in_pixel_space
    np.testing.assert_array_equal(optical.tic_to_image, cell_to_image)


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


class _LabellingReader(_StubReader):
    """A reader that copied its images out of a container and named them."""

    def __init__(self, labels: Dict[Path, Optional[OpticalImageLabel]]) -> None:
        super().__init__(list(labels), areas=[])
        self._labels = labels

    def get_optical_image_label(self, path: Path) -> Optional[OpticalImageLabel]:
        return self._labels.get(path)


def test_a_labelled_image_is_named_by_its_label(tmp_path: Path) -> None:
    """Not by the ``_0000`` / ``_0001`` rule, which reads vendor folders.

    mzpeak-convert numbers its members ``image_0000``, ``image_0001``: the
    rule called them the high resolution scan and the derived image.
    """
    rng = np.random.default_rng(3)
    labels: Dict[Path, Optional[OpticalImageLabel]] = {}
    for number, member in enumerate(("images/image_0000.svs", "images/image_0001.tif")):
        path = tmp_path / f"image_000{number}.tif"
        tifffile.imwrite(str(path), rng.integers(0, 256, (8, 8, 3), dtype=np.uint8))
        labels[path] = OpticalImageLabel(path.stem, member)
    optical = _optical(tmp_path, _LabellingReader(labels))
    optical.compute_alignment()
    images: Dict[str, Any] = {}
    optical.add_images(images)

    assert list(images) == ["ds_optical_image_0000", "ds_optical_image_0001"]
    assert optical.root_attr() == {
        "alignment_element": None,
        "elements": {
            "ds_optical_image_0000": {"source_file": "images/image_0000.svs"},
            "ds_optical_image_0001": {"source_file": "images/image_0001.tif"},
        },
    }


def test_an_unlabelled_image_keeps_the_vendor_rule(tmp_path: Path, tiff: Path) -> None:
    """A reader that gives no label for a path leaves it to the file name."""
    optical = _optical(tmp_path, _LabellingReader({tiff: None}))
    optical.compute_alignment()
    images: Dict[str, Any] = {}
    optical.add_images(images)

    assert list(images) == ["ds_optical_highres"]
    assert optical.sources == {"ds_optical_highres": IMAGE_FILE}


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
    # Back to the TIC coordinate, less the half cell its centre sits at,
    # then to micrometres: the spot of index r lands on r * pixel_size.
    to_index = np.array([[1.0, 0.0, -0.5], [0.0, 1.0, -0.5], [0.0, 0.0, 1.0]])
    expected = (
        np.array([[ps_x, 0.0, 0.0], [0.0, ps_y, 0.0], [0.0, 0.0, 1.0]])
        @ to_index
        @ np.linalg.inv(optical.tic_to_image)
    )
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
