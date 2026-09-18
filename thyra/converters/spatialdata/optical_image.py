"""Bounded-memory loading of optical images into a SpatialData store.

The optical images Thyra bundles with a conversion (FlexImaging brightfield,
typically a 10k x 40k RGB TIFF that decodes to 1-2 GB) used to travel the
same road as any other raster: decode the whole page into numpy, hand it to
``Image2DModel.parse`` with ``scale_factors``, and let ``SpatialData.write``
compute the pyramid. Measured on a 1.29 GB decode that road costs ~6.6 GB
over the conversion's baseline, for three reasons, none of them the image:

* dask's ``from_array`` copies every numpy array it is given (``x.copy()``
  in both ``dask.array.core.from_array`` and its expression-based twin, so
  no dask configuration avoids it), and ``Image2DModel.parse`` calls it for
  numpy input: the decoded page is in memory twice before anything is
  written;
* ``multiscale_spatial_image`` builds each pyramid level as
  ``coarsen(...).mean().astype(dtype)`` -- the mean promotes every 4096 x 4096
  block to float64 (128 MB per 16 MB of uint8) and the threaded scheduler
  runs one such block per core;
* all levels are computed in a single ``dask.compute``, so the scheduler is
  free to hold as many of those intermediates as it has threads.

This module keeps the store what that road produced while bounding memory
by one *band* of the source, not by the image:

1. **Declare, then stream.** The converter hands ``SpatialData``
   :meth:`StreamedOpticalImage.placeholder` -- the element with its final
   shapes, dtype, chunks and transformations, whose pixels are lazy zeros.
   spatialdata and ome-zarr create the arrays and write every piece of
   metadata exactly as before; zarr skips the all-zero chunks (its
   ``write_empty_chunks`` default), so declaring stores no chunk.
2. **Level 0 from the TIFF in bands.** :meth:`StreamedOpticalImage.stream_pixels`
   reads row bands through tifffile's zarr adapter, which decodes only the
   strips (or tile rows) a band needs, and writes them into the level-0
   array under :data:`BAND_BUDGET_BYTES`. Bands are whole strips: a page
   stored as one strip cannot be read in pieces, so it is decoded once,
   whole, as before.
3. **Each pyramid level from the one below, on disk.** Every level-*k* write
   unit is assembled one level-*k-1* chunk at a time, each reduced with the
   exact :func:`block_mean`, so a worker holds one source chunk, one
   accumulator and the chunk it is building whatever the image size.

**Formats other than TIFF.** FlexImaging does not always export a TIFF: an
acquisition's ``.mis`` can name a ``.jpg``, and ``.png`` and ``.bmp`` turn up
too. Those are single entropy-coded streams with no addressable row range --
a baseline JPEG cannot be decoded in pieces at all -- so
:class:`OpticalRasterSource` reads them through Pillow and yields the page
once, whole, where :class:`OpticalTiffSource` yields bands. Only step 2
changes: the decoded page still goes straight into the store's level-0 array
instead of through ``Image2DModel.parse`` (no dask copy), and every pyramid
level is still reduced from the level below on disk, which is where the
whole-page route spent most of its memory. :func:`probe_optical_source`
picks between the two by suffix, on the same list discovery filters on.

The pixel values are those ``xarray``'s ``coarsen`` produces: same
``boundary="trim", side="right"`` rule (an odd-sized axis drops its FIRST
row or column), same float64 mean, same truncating cast back to the source
dtype. :func:`block_mean` is unit-tested against xarray for exactly that.

**Failure policy.** An optical image that cannot be decoded used to be
skipped with a warning before anything was written; now the header is read
up front and the pixels only after the store exists, so
:meth:`OpticalImages.stream_pending_pixels` keeps that
contract by dropping the declared element from the store and warning. A
store never keeps an image whose pixels were not written.

**The alignment image is cropped to the section it aligns.** A FlexImaging
``.mis`` names one whole-slide scan and every ``.d`` cut from that slide
names the same one, so a project of eleven sections used to hold eleven
copies of a 28670 x 10606 photo, each placed so that *its own* Area lands on
*its own* MSI -- which is where the other ten sections' pixels are wrong by
centimetres. Measured on a real store: five copies of one photo, composed
placements up to 54,873 um apart. When a conversion covers exactly one
region -- a single-region acquisition, or one region of a multi-region file
selected with ``region=`` -- :class:`OpticalImages` crops the alignment
image to that region's Area box plus :data:`CROP_MARGIN_FRACTION`, streams
and pyramids the crop alone, and composes a :class:`Translation` by the
crop's origin in front of whatever transform the full image would have
carried. Nothing on the MSI side moves: the raster-to-image affine, the
pixel polygons and the ``coordinate_systems`` attr all keep describing the
full image's pixel grid, and the crop's own transform is what places its
pixels back where they were. The crop's origin and size, and the full
image's, are recorded in the store's ``optical_images`` root attr so the
whole slide can be found again without re-reading the ``.d``. A
whole-file conversion of a multi-region acquisition is unchanged.
"""

from __future__ import annotations

import itertools
import logging
import math
import os
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterator,
    List,
    Mapping,
    MutableMapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import dask.array as da
import numpy as np
import tifffile
import xarray as xr
import zarr
from numpy.typing import NDArray
from PIL import Image as PILImage
from spatialdata.models import Image2DModel
from spatialdata.transformations import Affine, Identity, Scale
from spatialdata.transformations import Sequence as SequenceTransform
from spatialdata.transformations import Translation, set_transformation

from ...alignment import AreaAlignmentResult, TeachingPointAlignment
from ._chunking import image_chunks

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ...alignment.teaching_points import RegionMapping
    from ...core.base_reader import BaseMSIReader

logger = logging.getLogger(__name__)

#: How much of a region's Area box is added on each side when the alignment
#: image is cropped to it, as a fraction of the box's own width and height.
#: Enough that the tissue edge and the slide around it are both in the crop;
#: on the real sections measured (a 1028 x 724 px Area on a 28670 x 10606
#: scan) it keeps the crop under one percent of the slide.
CROP_MARGIN_FRACTION = 0.10

#: Decoded bytes one band of the source TIFF may occupy while level 0 is
#: streamed in. Measured on a 36736 px wide RGB brightfield: a whole
#: 4096-row chunk row (451 MB) fits, so every level-0 chunk is encoded once;
#: a 256 MiB budget split it in two and cost a read-modify-write per chunk
#: (+0.7 s per chunk row) for a stage peak only 29 MB lower.
BAND_BUDGET_BYTES = 512 * 1024 * 1024

#: Output rows reduced at a time inside one source chunk. Bounds the
#: accumulator of the block mean to (rows x chunk width) values.
REDUCE_STRIP_ROWS = 512

#: Pyramid write units reduced concurrently. Each worker holds one source
#: chunk, one strip accumulator and the unit it is assembling.
REDUCE_WORKERS = 4

#: The layouts of a TIFF page this module reads: samples interleaved per
#: pixel (what FlexImaging exports and practically every RGB TIFF use),
#: a single sample, or one plane per sample.
PAGE_LAYOUTS = ("YXS", "YX", "SYX")

#: Suffixes :func:`probe_optical_source` routes to :class:`OpticalTiffSource`.
#: Everything else goes to :class:`OpticalRasterSource`.
TIFF_SUFFIXES = (".tif", ".tiff")

#: Pillow modes :class:`OpticalRasterSource` reads as they are, each mapped
#: to the ``(channels, dtype)`` the decoded array has. Anything not listed
#: (palette, CMYK, YCbCr, ...) is converted to RGB, or to RGBA when the file
#: carries transparency: those two are the conversions Pillow supports from
#: every mode, and converting is the only way to know the channel count at
#: probe time without decoding.
_PILLOW_MODES: Dict[str, Tuple[int, str]] = {
    "L": (1, "uint8"),
    "LA": (2, "uint8"),
    "RGB": (3, "uint8"),
    "RGBA": (4, "uint8"),
    "I;16": (1, "uint16"),
    "I;16B": (1, "uint16"),
    "I;16L": (1, "uint16"),
    "I": (1, "int32"),
    "F": (1, "float32"),
}

Shape = Tuple[int, int, int]


def block_mean(block: np.ndarray, factor: int) -> np.ndarray:
    """Mean over ``factor`` x ``factor`` windows of a ``(c, y, x)`` block.

    Reproduces ``xarray``'s ``coarsen(window).mean().astype(block.dtype)``
    value for value. Unsigned integers are summed exactly and floor-divided,
    which equals the truncated float64 mean for values that are never
    negative; signed integers take numpy's float64 mean (the accumulator
    dask picks for the same reduction); floating input keeps its own
    precision and skips NaNs the way xarray's ``nanmean`` does when any are
    present. The cast truncates the way ``ndarray.astype`` does.

    Args:
        block: ``(c, y, x)`` array whose ``y`` and ``x`` are multiples of
            ``factor``. Trim before calling; see :func:`trim_excess`.
        factor: Window edge along ``y`` and ``x``.

    Returns:
        Array of shape ``(c, y // factor, x // factor)`` and ``block.dtype``.
    """
    if block.ndim != 3:
        raise ValueError(f"block_mean expects a (c, y, x) block, got {block.ndim}-d")
    if factor < 1:
        raise ValueError(f"window factor must be positive, got {factor}")
    channels, rows, cols = block.shape
    if rows % factor or cols % factor:
        raise ValueError(
            f"block of shape {block.shape} is not a multiple of window {factor}"
        )
    if block.dtype.kind in "ub":
        # Exact integer path: the float64 sum of these values is exact and
        # /n truncates to the floor, so floor(sum / n) is the same number.
        n = factor * factor
        acc_dtype = (
            np.uint32 if int(np.iinfo(block.dtype).max) * n < 2**32 else np.uint64
        )
        if block.dtype.kind == "b":
            acc_dtype = np.uint32
        acc = np.zeros((channels, rows // factor, cols // factor), dtype=acc_dtype)
        for dy in range(factor):
            for dx in range(factor):
                acc += block[:, dy::factor, dx::factor]
        return (acc // n).astype(block.dtype)
    windows = block.reshape(channels, rows // factor, factor, cols // factor, factor)
    if block.dtype.kind in "fc" and np.isnan(windows).any():
        return np.nanmean(windows, axis=(2, 4)).astype(block.dtype)
    return np.mean(windows, axis=(2, 4)).astype(block.dtype)


def trim_excess(size: int, factor: int) -> int:
    """Leading rows/columns xarray's ``coarsen`` drops for ``boundary="trim"``.

    With ``side="right"`` the excess comes off the FRONT of the axis, so the
    last window always ends on the last element. This is the rule
    ``multiscale_spatial_image`` uses for every pyramid level.
    """
    return size - (size // factor) * factor


def level_shapes(shape: Shape, scale_factors: Sequence[int]) -> List[Shape]:
    """``(c, y, x)`` of every pyramid level, level 0 first.

    Each factor divides the spatial axes of the level before it, flooring
    the way ``coarsen(boundary="trim")`` does. The channel axis is never
    reduced.
    """
    shapes = [(int(shape[0]), int(shape[1]), int(shape[2]))]
    for factor in scale_factors:
        c, y, x = shapes[-1]
        shapes.append((c, y // factor, x // factor))
    return shapes


def band_rows(
    row_bytes: int,
    chunk_rows: int,
    strip_rows: int = 1,
    budget: int = BAND_BUDGET_BYTES,
) -> int:
    """Rows per band: whole strips, within ``budget``, tiling the chunk rows.

    A strip (or tile row) is the smallest unit the TIFF decoder produces, so
    a band is always a whole number of them and never shorter than one --
    a single-strip page yields one band, the whole page, decoded once.
    Within that, the largest divisor of ``chunk_rows`` that fits is
    preferred so bands complete whole level-0 chunks instead of straddling
    them; when the strip height does not divide the chunk height no such
    divisor exists and the budget alone decides, which only means the
    chunks are completed by more than one band.
    """
    if row_bytes <= 0 or chunk_rows <= 0 or strip_rows <= 0:
        raise ValueError("row_bytes, chunk_rows and strip_rows must be positive")
    rows = max(strip_rows, (budget // row_bytes) // strip_rows * strip_rows)
    if rows >= chunk_rows and chunk_rows % strip_rows == 0:
        return chunk_rows
    for candidate in range(min(rows, chunk_rows), 0, -1):
        if chunk_rows % candidate == 0 and candidate % strip_rows == 0:
            return candidate
    return rows


def _write_unit(array: zarr.Array) -> Tuple[int, ...]:
    """The block one write of ``array`` should cover: its shard, else its chunk.

    ``zarr.Array.chunks`` is the inner chunk once an array is sharded, and
    concurrent writes inside one shard clobber each other, so the streaming
    and the reducer key on the shard whenever there is one.
    """
    return tuple(int(n) for n in (array.shards or array.chunks))


def _write_regions(array: zarr.Array) -> Iterator[Tuple[slice, slice, slice]]:
    """Every write unit of a ``(c, y, x)`` array as slices, C order."""
    unit = _write_unit(array)
    ranges = [range(0, size, step) for size, step in zip(array.shape, unit)]
    for starts in itertools.product(*ranges):
        c, y, x = (
            slice(start, min(start + step, size))
            for start, step, size in zip(starts, unit, array.shape)
        )
        yield c, y, x


@dataclass(frozen=True)
class OpticalTiffSource:
    """What the first page of an optical TIFF looks like, without decoding it.

    ``shape`` is always ``(c, y, x)``. ``page_axes`` is one of
    :data:`PAGE_LAYOUTS` and says how a decoded array maps onto that;
    ``strip_rows`` is the height of the page's strips or tiles, the smallest
    row range the decoder can produce on its own.
    """

    path: Path
    shape: Shape
    dtype: np.dtype
    page_axes: str
    strip_rows: int

    @classmethod
    def probe(cls, path: Union[str, Path]) -> "OpticalTiffSource":
        """Read the page header.

        Raises:
            ValueError: for a page layout other than :data:`PAGE_LAYOUTS`
                or a sample format tifffile cannot map to a dtype -- the
                cases the whole-page decode used to reject, at the same
                point, before anything is written.
        """
        path = Path(path)
        with tifffile.TiffFile(path) as tif:
            page = tif.pages[0]
            page_shape = tuple(int(n) for n in page.shape)
            dtype = page.dtype
            axes = str(page.axes)
            chunks = tuple(int(n) for n in page.chunks)
        if dtype is None:
            raise ValueError(f"{path.name}: sample format not supported by tifffile")
        if axes == "YX":
            shape: Shape = (1, page_shape[0], page_shape[1])
            strip_rows = chunks[0]
        elif axes == "YXS":
            shape = (page_shape[2], page_shape[0], page_shape[1])
            strip_rows = chunks[0]
        elif axes == "SYX":
            shape = (page_shape[0], page_shape[1], page_shape[2])
            # tifffile reports planar strips as (rows, width).
            strip_rows = chunks[-2]
        else:
            raise ValueError(
                f"{path.name}: unsupported TIFF page layout {axes} {page_shape}"
            )
        return cls(
            path=path,
            shape=shape,
            dtype=np.dtype(dtype),
            page_axes=axes,
            strip_rows=max(1, min(strip_rows, shape[1])),
        )

    @property
    def row_bytes(self) -> int:
        """Decoded bytes of one full-width row across all channels."""
        return int(self.shape[0]) * int(self.shape[2]) * int(self.dtype.itemsize)

    def bands(
        self, rows: int, window: Optional[Window] = None
    ) -> Iterator[Tuple[int, np.ndarray]]:
        """Yield ``(first_row, band)`` with ``band`` shaped ``(c, rows, x)``.

        Row ranges are decoded through tifffile's zarr adapter, so only the
        strips of a band are ever decoded. When ``rows`` covers the page the
        page is decoded once, whole, with tifffile's own parallel decoder.

        ``window`` is ``(y_start, y_stop, x_start, x_stop)`` in this page's
        pixels: only those rows are decoded, only those columns are kept,
        and ``first_row`` stays in the page's own row numbering. With a
        window the adapter is always used, so a single-strip page decodes
        its one strip and keeps the window of it.
        """
        n_rows = self.shape[1]
        with tifffile.TiffFile(self.path) as tif:
            page = tif.pages[0]
            if window is None and rows >= n_rows:
                yield 0, self._to_cyx(page.asarray())
                return
            y_start, y_stop, x_start, x_stop = window or (0, n_rows, 0, self.shape[2])
            source = zarr.open_array(store=page.aszarr(), mode="r")
            for start in range(y_start, y_stop, rows):
                stop = min(start + rows, y_stop)
                if self.page_axes == "SYX":
                    decoded = source[:, start:stop, x_start:x_stop]
                else:
                    decoded = source[start:stop, x_start:x_stop]
                yield start, self._to_cyx(decoded)

    def _to_cyx(self, decoded: np.ndarray) -> np.ndarray:
        if self.page_axes == "YX":
            return decoded[np.newaxis, :, :]
        if self.page_axes == "YXS":
            return np.moveaxis(decoded, -1, 0)
        return decoded


@contextmanager
def _unlimited_pixels() -> Iterator[None]:
    """Pillow's decompression-bomb ceiling lifted for one open or decode.

    ``Image.MAX_IMAGE_PIXELS`` warns above ~89 megapixels and raises above
    twice that. It guards a process that opens images it did not ask for;
    the file here is the acquisition's own optical scan, sitting next to the
    raw data the user pointed Thyra at, and a whole-slide brightfield export
    passes both thresholds legitimately -- a 36736 x 21000 FlexImaging scan
    is 771 Mpx. The TIFF path has no such ceiling, and a refusal here would
    not be visible as one: the per-image guard turns it into a warning and
    the store silently loses the image. Restored on the way out rather than
    cleared once at import, because the attribute is process-wide and shared
    with every other library in the process.
    """
    previous = PILImage.MAX_IMAGE_PIXELS
    PILImage.MAX_IMAGE_PIXELS = None
    try:
        yield
    finally:
        PILImage.MAX_IMAGE_PIXELS = previous


@dataclass(frozen=True)
class OpticalRasterSource:
    """An optical image in a format that is not TIFF: JPEG, PNG or BMP.

    Same surface as :class:`OpticalTiffSource` -- ``shape`` is ``(c, y, x)``,
    :meth:`bands` yields ``(first_row, band)`` -- so
    :class:`StreamedOpticalImage` never has to know which of the two it
    holds. What differs is that these are single entropy-coded streams with
    no addressable row range, so ``strip_rows`` is the whole height and
    :meth:`bands` decodes the page once, whole, however small the band
    budget is. That is the same thing :class:`OpticalTiffSource` already
    does for a single-strip TIFF.

    ``convert_to`` is the Pillow mode the decode goes through, or ``None``
    when the file's own mode is one :data:`_PILLOW_MODES` reads directly.
    It is decided from the header at :meth:`probe` time so the declared
    channel count and dtype are the ones the pixels will actually have.
    """

    path: Path
    shape: Shape
    dtype: np.dtype
    mode: str
    convert_to: Optional[str]
    strip_rows: int

    @classmethod
    def probe(cls, path: Union[str, Path]) -> "OpticalRasterSource":
        """Read the header.

        Pillow's ``open`` parses the header and stops, so nothing is decoded
        here.

        Raises:
            OSError: if the file is not an image Pillow recognises
                (``UnidentifiedImageError`` is one), or is truncated before
                the header ends -- the same point, and the same per-image
                skip, as :meth:`OpticalTiffSource.probe`.
        """
        path = Path(path)
        with _unlimited_pixels(), PILImage.open(path) as page:
            mode = str(page.mode)
            width, height = (int(n) for n in page.size)
            # Pillow >= 9.5. Says whether the file carries alpha at all,
            # including the palette-with-transparency case that plain
            # ``"A" in mode`` misses.
            transparent = bool(getattr(page, "has_transparency_data", False))
        spec = _PILLOW_MODES.get(mode)
        convert_to: Optional[str] = None
        if spec is None:
            convert_to = "RGBA" if transparent else "RGB"
            spec = _PILLOW_MODES[convert_to]
        channels, dtype = spec
        if width < 1 or height < 1:
            raise ValueError(f"{path.name}: image is {width}x{height}")
        return cls(
            path=path,
            shape=(channels, height, width),
            dtype=np.dtype(dtype),
            mode=mode,
            convert_to=convert_to,
            strip_rows=height,
        )

    @property
    def row_bytes(self) -> int:
        """Decoded bytes of one full-width row across all channels."""
        return int(self.shape[0]) * int(self.shape[2]) * int(self.dtype.itemsize)

    def bands(
        self, rows: int, window: Optional[Window] = None
    ) -> Iterator[Tuple[int, np.ndarray]]:
        """Yield ``(first_row, page)``: these formats decode whole or not at all.

        ``rows`` is accepted for the shared surface and ignored;
        :func:`band_rows` returns the full height for a ``strip_rows`` this
        large anyway, so the caller asks for the whole page regardless.
        With a ``window`` the page is still decoded whole -- there is no
        other way -- and the window of it is what comes out, with
        ``first_row`` at the window's top.
        """
        del rows
        with _unlimited_pixels(), PILImage.open(self.path) as page:
            frame = page if self.convert_to is None else page.convert(self.convert_to)
            decoded = np.asarray(frame)
        cyx = self._to_cyx(decoded)
        if window is None:
            yield 0, cyx
            return
        y_start, y_stop, x_start, x_stop = window
        yield y_start, cyx[:, y_start:y_stop, x_start:x_stop]

    def _to_cyx(self, decoded: np.ndarray) -> np.ndarray:
        """``(y, x)`` or ``(y, x, s)`` as Pillow returns it, to ``(c, y, x)``."""
        if decoded.ndim == 2:
            decoded = decoded[:, :, np.newaxis]
        if decoded.ndim != 3:
            raise ValueError(
                f"{self.path.name}: decoded to a {decoded.ndim}-d array, "
                f"expected 2 or 3 dimensions"
            )
        moved = np.moveaxis(decoded, -1, 0)
        if moved.shape != self.shape:
            raise ValueError(
                f"{self.path.name}: decoded to {moved.shape}, "
                f"but its header said {self.shape}"
            )
        return moved


#: ``(y_start, y_stop, x_start, x_stop)`` in a page's own pixels: the part of
#: it a source is asked to decode.
Window = Tuple[int, int, int, int]

#: What :func:`probe_optical_source` returns. The two classes share the
#: surface :class:`StreamedOpticalImage` uses, not a base class: they have
#: nothing else in common, and a Protocol would only restate this line.
OpticalSource = Union[OpticalTiffSource, OpticalRasterSource]


@dataclass(frozen=True)
class OpticalCrop:
    """Where a cropped alignment image sits in the whole-slide scan it was cut from.

    Everything is in the full image's pixels, integers, origin at its top
    left: ``x0``/``y0`` are the crop's first column and row, ``width`` and
    ``height`` its size, ``full_width``/``full_height`` the scan's. The
    crop's own pixel ``(u, v)`` is the scan's pixel ``(x0 + u, y0 + v)``,
    which is exactly the :class:`Translation` the cropped element carries.
    """

    x0: int
    y0: int
    width: int
    height: int
    full_width: int
    full_height: int
    region_id: int
    region_name: str

    @classmethod
    def around(
        cls,
        box: Tuple[float, float, float, float],
        full_width: int,
        full_height: int,
        *,
        mapping: "RegionMapping",
        margin: float = CROP_MARGIN_FRACTION,
    ) -> Optional["OpticalCrop"]:
        """The crop for a region whose pixels land in ``box``, or ``None``.

        ``box`` is ``(x_min, x_max, y_min, y_max)`` in full-image pixels --
        where the region's raster lands under the alignment affine. It is
        padded by ``margin`` of its own width and height on every side,
        rounded outward to whole pixels and clamped to the image. ``None``
        when there is nothing to cut: the padded box covers the whole scan,
        or does not meet it at all (an Area drawn off the image, which the
        affine can produce from a ``.mis`` whose areas do not match the
        data).
        """
        x_min, x_max, y_min, y_max = box
        pad_x = (x_max - x_min) * margin
        pad_y = (y_max - y_min) * margin
        x0 = max(0, int(math.floor(x_min - pad_x)))
        y0 = max(0, int(math.floor(y_min - pad_y)))
        x1 = min(full_width, int(math.ceil(x_max + pad_x)))
        y1 = min(full_height, int(math.ceil(y_max + pad_y)))
        if x1 <= x0 or y1 <= y0:
            return None
        if x0 == 0 and y0 == 0 and x1 == full_width and y1 == full_height:
            return None
        return cls(
            x0=x0,
            y0=y0,
            width=x1 - x0,
            height=y1 - y0,
            full_width=full_width,
            full_height=full_height,
            region_id=mapping.region_id,
            region_name=mapping.name,
        )

    @property
    def window(self) -> Window:
        """The crop as the row and column range a source decodes."""
        return (self.y0, self.y0 + self.height, self.x0, self.x0 + self.width)

    def to_attr(self) -> Dict[str, Any]:
        """The crop as the store's ``optical_images`` root attr records it."""
        return {
            "origin": [self.x0, self.y0],
            "size": [self.width, self.height],
            "full_size": [self.full_width, self.full_height],
            "region_id": self.region_id,
            "region_name": self.region_name,
            "margin_fraction": CROP_MARGIN_FRACTION,
        }


@dataclass(frozen=True)
class CroppedOpticalSource:
    """A window of an optical source, with the surface of a whole one.

    :class:`StreamedOpticalImage` never learns it holds a crop: ``shape`` is
    the window's, ``bands`` yields the window's rows numbered from the
    window's top, and the pyramid is built from the window's size. The
    decode itself is the inner source's, asked for the window, so a TIFF
    still decodes only the strips the window touches.
    """

    inner: OpticalSource
    crop: OpticalCrop

    @property
    def path(self) -> Path:
        """The file the inner source reads."""
        return self.inner.path

    @property
    def shape(self) -> Shape:
        """``(c, y, x)`` of the window."""
        return (self.inner.shape[0], self.crop.height, self.crop.width)

    @property
    def dtype(self) -> np.dtype:
        """The inner source's dtype."""
        return self.inner.dtype

    @property
    def strip_rows(self) -> int:
        """The inner source's strip height: the decode unit is unchanged."""
        return self.inner.strip_rows

    @property
    def row_bytes(self) -> int:
        """Decoded bytes of one window-wide row across all channels."""
        return int(self.shape[0]) * int(self.crop.width) * int(self.dtype.itemsize)

    def bands(self, rows: int) -> Iterator[Tuple[int, np.ndarray]]:
        """Yield ``(first_row, band)`` in the window's own row numbering."""
        for start, band in self.inner.bands(rows, window=self.crop.window):
            yield start - self.crop.y0, band


#: What :class:`StreamedOpticalImage` streams from: a whole source, or a
#: window of one.
AnyOpticalSource = Union[OpticalSource, CroppedOpticalSource]


def probe_optical_source(path: Union[str, Path]) -> OpticalSource:
    """Read the header of an optical image, whatever raster format it is.

    TIFF goes to :class:`OpticalTiffSource`, which decodes row bands; every
    other format goes to :class:`OpticalRasterSource`, which decodes the
    page whole. The split is by suffix, from the same list discovery filters
    on (``BrukerFolderStructure.OPTICAL_IMAGE_SUFFIXES``), so the two agree
    on what a ``.jpg`` is.

    Raises:
        Whatever the chosen probe raises for a file it cannot read. Callers
        skip that image with a warning; see the failure policy above.
    """
    path = Path(path)
    if path.suffix.lower() in TIFF_SUFFIXES:
        return OpticalTiffSource.probe(path)
    return OpticalRasterSource.probe(path)


@dataclass
class StreamedOpticalImage:
    """One optical image: declared to SpatialData up front, pixels streamed after.

    Build it, put :meth:`placeholder` in the images the SpatialData write
    carries, and call :meth:`stream_pixels` on the store once that write
    has returned. Between the two the element on disk has all its metadata
    and no pixels; :meth:`discard` removes it again if the pixels cannot
    be streamed.

    There is deliberately no hook for per-element attributes. A store's
    ``images/*`` group carries exactly ``ome`` and ``spatialdata_attrs``:
    the ome-zarr writer spatialdata hands the element to composes the group
    metadata itself and does not carry anything else across, so whatever a
    caller puts in the element's ``.attrs`` is dropped by the write -- and
    on the multiscale path ``xr.DataTree.from_dict`` drops it before that,
    so it never even reaches the writer. This class used to take an
    ``attrs`` mapping and the converter used to fill it with the image's
    source filename; nothing of it was ever in a store. Provenance that has
    to survive belongs in the store's root attrs, which is where
    :meth:`OpticalImages.root_attr` now puts
    it, under ``optical_images``.
    """

    source: AnyOpticalSource
    name: str
    chunks: Tuple[int, ...]
    scale_factors: Sequence[int]
    transformations: Dict[str, Any]

    def __post_init__(self) -> None:
        if len(self.chunks) != 3:
            raise ValueError(f"chunks must be (c, y, x), got {self.chunks}")
        self.scale_factors = [int(f) for f in self.scale_factors]

    @property
    def level_shapes(self) -> List[Shape]:
        """``(c, y, x)`` of every pyramid level of this image, level 0 first."""
        return level_shapes(self.source.shape, self.scale_factors)

    def placeholder(self) -> Union[xr.DataArray, xr.DataTree]:
        """The element as ``Image2DModel.parse`` would build it, pixels lazy zeros.

        With scale factors this is the multiscale DataTree ``parse`` builds
        for ``scale_factors``: same levels, chunks, per-level transformations
        and coordinates. Without, it is the single-scale DataArray ``parse``
        returns, so the store gets the single-scale writer it always had.
        Either way the write creates the arrays and metadata and stores no
        chunk.
        """
        dims = ("c", "y", "x")
        if not self.scale_factors:
            image = Image2DModel.parse(
                da.zeros(
                    self.source.shape, chunks=self.chunks, dtype=self.source.dtype
                ),
                dims=dims,
                transformations=dict(self.transformations),
                chunks=self.chunks,
            )
            return image
        base = self.level_shapes[0]
        levels: Dict[str, xr.Dataset] = {}
        for index, shape in enumerate(self.level_shapes):
            coords: Dict[str, Any] = {"c": np.arange(shape[0])}
            for axis, full, size in zip(dims[1:], base[1:], shape[1:]):
                # Pixel centres in level-0 units; spatialdata's own
                # compute_coordinates for a DataTree, verbatim.
                coords[axis] = np.linspace(0, full, size + 1)[:-1] + full / size / 2
            image = xr.DataArray(
                da.zeros(shape, chunks=self.chunks, dtype=self.source.dtype),
                dims=dims,
                coords=coords,
            )
            levels[f"scale{index}"] = xr.Dataset({"image": image})
        tree = xr.DataTree.from_dict(levels)
        set_transformation(tree, dict(self.transformations), set_all=True)
        Image2DModel.validate(tree)
        return tree

    def stream_pixels(self, store_path: Union[str, Path]) -> None:
        """Fill the element's arrays in ``store_path`` with the real pixels.

        Requires the store to hold the element :meth:`placeholder` declared.
        Level 0 comes from the source in bands (one whole-page band for a
        format that cannot be decoded in pieces); each further level from
        the level below, one write unit at a time.
        """
        arrays = self._open_level_arrays(store_path)
        self._stream_level0(arrays[0])
        for index, factor in enumerate(self.scale_factors, start=1):
            self._reduce_level(arrays[index - 1], arrays[index], factor)
        logger.info(
            f"Streamed optical image '{self.name}': level 0 plus "
            f"{len(self.scale_factors)} pyramid level"
            f"{'s' if len(self.scale_factors) != 1 else ''}"
        )

    def discard(self, store_path: Union[str, Path]) -> None:
        """Remove the declared element from ``store_path``, if it is there."""
        root = zarr.open_group(str(store_path), mode="r+", use_consolidated=False)
        images = root["images"]
        if self.name in images:
            del images[self.name]

    def _open_level_arrays(self, store_path: Union[str, Path]) -> List[zarr.Array]:
        root = zarr.open_group(str(store_path), mode="r+", use_consolidated=False)
        element = root[f"images/{self.name}"]
        try:
            datasets = element.attrs["ome"]["multiscales"][0]["datasets"]
        except (KeyError, IndexError, TypeError) as exc:
            raise RuntimeError(
                f"images/{self.name} in {store_path} carries no multiscales metadata"
            ) from exc
        paths = [dataset["path"] for dataset in datasets]
        expected = self.level_shapes
        if len(paths) != len(expected):
            raise RuntimeError(
                f"images/{self.name} has {len(paths)} pyramid levels on disk, "
                f"expected {len(expected)}"
            )
        arrays: List[zarr.Array] = []
        for path, shape in zip(paths, expected):
            array = element[path]
            if not isinstance(array, zarr.Array) or tuple(array.shape) != tuple(shape):
                raise RuntimeError(
                    f"images/{self.name}/{path} is {getattr(array, 'shape', None)}, "
                    f"expected {shape}"
                )
            arrays.append(array)
        return arrays

    def _stream_level0(self, target: zarr.Array) -> None:
        chunk_rows = _write_unit(target)[1]
        # The budget is looked up here, not bound as a default, so a test
        # (or a caller) can narrow it.
        rows = band_rows(
            self.source.row_bytes,
            chunk_rows,
            self.source.strip_rows,
            budget=BAND_BUDGET_BYTES,
        )
        n_bands = -(-self.source.shape[1] // rows)
        logger.debug(
            f"  level 0 of '{self.name}': {n_bands} band(s) of {rows} rows "
            f"({rows * self.source.row_bytes / 1e6:.0f} MB decoded each)"
        )
        for start, band in self.source.bands(rows):
            target[:, start : start + band.shape[1], :] = band
            # The next band is decoded before this name is rebound; drop
            # it now so only one band is ever alive.
            del band

    def _reduce_level(
        self, source: zarr.Array, target: zarr.Array, factor: int
    ) -> None:
        """Write every unit of ``target`` from its footprint in ``source``.

        A unit's footprint is walked one source write unit at a time (an
        odd trim shifts the walk by the excess, so the last step of a row or
        column may touch a second source chunk), and each piece is reduced
        in strips of :data:`REDUCE_STRIP_ROWS` output rows.
        """
        excess_y = trim_excess(source.shape[1], factor)
        excess_x = trim_excess(source.shape[2], factor)
        _, src_rows, src_cols = _write_unit(source)
        out_rows_per_read = max(1, src_rows // factor)
        out_cols_per_read = max(1, src_cols // factor)
        regions = list(_write_regions(target))
        workers = max(1, min(REDUCE_WORKERS, os.cpu_count() or 1, len(regions)))

        def reduce_region(region: Tuple[slice, slice, slice]) -> None:
            c_sel, y_sel, x_sel = region
            unit = np.empty(
                (
                    c_sel.stop - c_sel.start,
                    y_sel.stop - y_sel.start,
                    x_sel.stop - x_sel.start,
                ),
                dtype=target.dtype,
            )
            for y0 in range(0, unit.shape[1], out_rows_per_read):
                y1 = min(y0 + out_rows_per_read, unit.shape[1])
                for x0 in range(0, unit.shape[2], out_cols_per_read):
                    x1 = min(x0 + out_cols_per_read, unit.shape[2])
                    piece = source[
                        c_sel,
                        excess_y
                        + (y_sel.start + y0) * factor : excess_y
                        + (y_sel.start + y1) * factor,
                        excess_x
                        + (x_sel.start + x0) * factor : excess_x
                        + (x_sel.start + x1) * factor,
                    ]
                    for r0 in range(0, y1 - y0, REDUCE_STRIP_ROWS):
                        r1 = min(r0 + REDUCE_STRIP_ROWS, y1 - y0)
                        unit[:, y0 + r0 : y0 + r1, x0:x1] = block_mean(
                            piece[:, r0 * factor : r1 * factor, :], factor
                        )
                    del piece
            target[region] = unit

        with ThreadPoolExecutor(max_workers=workers) as pool:
            # list() so a failing unit raises here, not silently.
            list(pool.map(reduce_region, regions))


def _calc_optical_scale_factors(
    smallest_dim: int,
    min_coarsest_size: int = 1000,
    factor: int = 2,
) -> list:
    """Pick pyramid scale factors for an optical image of given size.

    Decides how many cumulative-doubling downsample levels to generate
    based on the smallest spatial dimension; stops once the next
    halving would drop the short side below ``min_coarsest_size``.

    Mirrors :func:`spatialdata_io.readers._utils._utils.calc_scale_factors`
    so wizard-converted microscopy and Thyra-bundled FlexImaging
    brightfield share the same pyramid shape Xenium's own morphology
    image gets out of spatialdata-io.

    Returns a list of (cumulative) downsample factors to feed to
    ``Image2DModel.parse(scale_factors=...)``.  Empty list means
    "no pyramid needed" (image is already at or below the coarsest
    target size).

    Args:
        smallest_dim: ``min(width, height)`` of the source image.
        min_coarsest_size: stop adding levels once the next halving
            would drop the short side below this value.  Default
            ~1000 px matches spatialdata-io's convention and gives a
            coarsest level that comfortably fits a single viewport
            paint in a few hundred KB.
        factor: downsample factor per step.  Default 2 (mip-map style).
    """
    if smallest_dim <= 0:
        return []
    factors: list = []
    cur = smallest_dim / factor
    while cur >= min_coarsest_size:
        factors.append(factor)
        cur /= factor
    return factors


class OpticalImages:
    """The optical images of one conversion: alignment, declaration, streaming.

    One conversion owns one of these. It holds every piece of state the
    optical images need -- the FlexImaging alignment, the TIC-to-image
    affine, which file became which element, which of those elements is
    the alignment image, and the placeholders still waiting for their
    pixels -- and it is the only thing that writes any of them.

    The three steps, in the order a conversion takes them:

    1. :meth:`compute_alignment` and :meth:`build_tic_to_image_affine`,
       while the reader's metadata is being read. Both run regardless of
       ``apply_alignment``: the opt-out path needs the affine's *inverse*
       to carry the optical photo into micrometers.
    2. :meth:`add_images`, while the elements are assembled -- each image
       is declared as a placeholder with its final shape, chunks and
       transformations, and nothing is decoded.
    3. :meth:`stream_pending_pixels`, once the store exists, which fills
       each placeholder from its file in bands.

    What the converter reads back is the four read-only properties
    (:attr:`alignment`, :attr:`tic_to_image`, :attr:`sources`,
    :attr:`alignment_element`) plus :attr:`apply_alignment`; it never
    writes any of them.
    """

    def __init__(
        self,
        reader: "BaseMSIReader",
        output_path: Path,
        dataset_id: str,
        *,
        include: bool,
        apply_alignment: bool,
        pixel_size_xy: Callable[[], Tuple[float, float]],
    ) -> None:
        """Build the optical images of one conversion.

        Args:
            reader: The MSI reader. Only ``get_optical_image_paths``,
                ``get_primary_optical_image_path``, ``mis_metadata``,
                ``_positions`` and ``_header`` are read, and the last four
                through ``getattr``/``hasattr`` -- a reader here is whatever
                satisfies the interface, not necessarily a
                :class:`~thyra.core.base_reader.BaseMSIReader` subclass.
            output_path: The store being written. Read only after it exists,
                by :meth:`stream_pending_pixels` and :meth:`forget_image`.
            dataset_id: Names the per-dataset coordinate system each image is
                placed in, and prefixes every element name.
            include: Whether to put optical images in the store at all.
            apply_alignment: Whether the MSI raster itself is placed in
                optical-photo pixels. See :meth:`_load_single_image` for what
                the two modes do to an image's transform.
            pixel_size_xy: The conversion's in-plane pitch as ``(x_um, y_um)``,
                called when an image is loaded rather than read here: the
                converter settles its pitch from the reader's metadata after
                this object is built.
        """
        self.reader = reader
        self.output_path = output_path
        self.dataset_id = dataset_id
        self._include = bool(include)
        self._apply_alignment = bool(apply_alignment)
        self._pixel_size_xy = pixel_size_xy
        # Optical-MSI alignment (computed from FlexImaging Area definitions)
        self._alignment: Optional[AreaAlignmentResult] = None
        # Affine matrix mapping TIC raster indices to optical image pixels
        self._tic_to_image: Optional[NDArray[np.float64]] = None
        # Primary optical image filename from .mis <ImageFile> and its
        # dimensions. The filename is public because a test that has no
        # .mis to read plants it; the dimensions are only ever derived
        # here, from the primary image's own header.
        self.primary_filename: Optional[str] = None
        self._primary_dims: Optional[Tuple[int, int]] = None  # (width, height)
        # What the store will say about its optical images: which element
        # each file became, and which of those elements is the alignment
        # image. Filled as the images are declared, read when the root
        # attrs are composed. See :meth:`root_attr`.
        self._sources: Dict[str, str] = {}
        self._alignment_element: Optional[str] = None
        # The window of the alignment image this conversion embeds, when it
        # covers exactly one region. See the module docstring and
        # :meth:`_crop_for_alignment_image`.
        self._crop: Optional[OpticalCrop] = None
        # Optical images declared to SpatialData as placeholders whose pixels
        # still have to be streamed into the store once it is written. See
        # :class:`StreamedOpticalImage` and :meth:`stream_pending_pixels`.
        self._pending: Dict[str, StreamedOpticalImage] = {}

    @property
    def alignment(self) -> Optional[AreaAlignmentResult]:
        """The FlexImaging area alignment, or ``None`` if there is none."""
        return self._alignment

    @property
    def tic_to_image(self) -> Optional[NDArray[np.float64]]:
        """The 3x3 affine from TIC raster indices to optical-photo pixels."""
        return self._tic_to_image

    @property
    def sources(self) -> Mapping[str, str]:
        """Element name to the name of the file it was read from."""
        return self._sources

    @property
    def alignment_element(self) -> Optional[str]:
        """The element the alignment is stated against, if it is in the store."""
        return self._alignment_element

    @property
    def apply_alignment(self) -> bool:
        """Whether the caller asked for the raster in optical-photo pixels."""
        return self._apply_alignment

    @property
    def msi_in_pixel_space(self) -> bool:
        """Whether the MSI raster actually lands in optical-photo pixels.

        Both halves matter. The affine is built whenever FlexImaging
        alignment data exists, ``apply_alignment`` or not -- it is needed
        either way, because opting out uses its *inverse* to carry the
        optical photo into micrometers. So the matrix being present says
        only that an alignment is available, never that it was applied to
        the raster; and the flag being set says only that it was asked
        for, never that there was anything to apply.

        Reading only the matrix is what let the store's own
        ``coordinate_systems`` attr declare ``unit="pixel"`` on a store
        whose every element was in micrometers (issue #288). The TIC
        image's transform, the pixel polygons and that attr each spelled
        the condition out separately and one of the three spelled it
        differently, so it is written here, once, on the object that
        holds both halves.

        The polygons additionally require :attr:`alignment`, which they
        need for the transform itself rather than for the decision. They
        are not gated on it *alone*: the matrix is only built when
        ``region_mappings`` is non-empty, so a .mis whose areas match no
        region left the polygons on the alignment branch while everything
        else took the micrometer one, and every position failed to
        transform.
        """
        return self._apply_alignment and self._tic_to_image is not None

    @property
    def pending(self) -> Mapping[str, StreamedOpticalImage]:
        """The placeholders whose pixels have not been streamed yet."""
        return self._pending

    @property
    def crop(self) -> Optional[OpticalCrop]:
        """The window of the alignment image the store holds, or ``None`` for all of it."""
        return self._crop

    def _converted_region_mapping(self) -> Optional["RegionMapping"]:
        """The one region this conversion covers, or ``None`` if it covers several.

        A single-region acquisition has one mapping. A multi-region file
        converted with ``region=`` still maps every Area -- the reader's
        positions are not filtered -- but the reader says which region it
        selected, and that region's mapping is the one. A multi-region
        file converted whole has several, and ``None`` here is what keeps
        that conversion unchanged.
        """
        if self._alignment is None:
            return None
        mappings = self._alignment.region_mappings
        selected = getattr(self.reader, "_selected_region", None)
        if selected is not None:
            for mapping in mappings:
                if mapping.region_id == selected:
                    return mapping
            return None
        if len(mappings) == 1:
            return mappings[0]
        return None

    def _crop_for_alignment_image(
        self, full_width: int, full_height: int
    ) -> Optional[OpticalCrop]:
        """Where to cut the alignment image, or ``None`` to keep it whole.

        The box is where the converted region's raster actually lands under
        :attr:`tic_to_image` -- the same affine the TIC and the polygons are
        placed by -- not the Area's own corners. For a single region the two
        coincide by construction (the affine stretches the raster bounds
        onto the Area); for a selected region of a multi-region file the
        affine is the global stretch and the Area corners are only where
        the ``.mis`` drew them, so the affine's answer is the one that
        agrees with the MSI. Half a raster step is added on each side so the
        outermost pixels' full footprint is inside the box before the
        margin is.
        """
        mapping = self._converted_region_mapping()
        if mapping is None or self._alignment is None or self._tic_to_image is None:
            return None
        matrix = self._tic_to_image
        scale_x, scale_y = float(matrix[0, 0]), float(matrix[1, 1])
        # Raster bounds in the TIC grid's own indices, i.e. relative to the
        # position the affine calls (0, 0).
        rx0 = mapping.raster_min_x - self._alignment.first_raster_x
        rx1 = mapping.raster_max_x - self._alignment.first_raster_x
        ry0 = mapping.raster_min_y - self._alignment.first_raster_y
        ry1 = mapping.raster_max_y - self._alignment.first_raster_y
        x_lo = float(matrix[0, 2]) + rx0 * scale_x - abs(scale_x) / 2.0
        x_hi = float(matrix[0, 2]) + rx1 * scale_x + abs(scale_x) / 2.0
        y_lo = float(matrix[1, 2]) + ry0 * scale_y - abs(scale_y) / 2.0
        y_hi = float(matrix[1, 2]) + ry1 * scale_y + abs(scale_y) / 2.0
        box = (min(x_lo, x_hi), max(x_lo, x_hi), min(y_lo, y_hi), max(y_lo, y_hi))
        return OpticalCrop.around(box, full_width, full_height, mapping=mapping)

    def compute_alignment(self) -> None:
        """Compute optical-MSI alignment from reader metadata.

        For FlexImaging data with Area definitions, this computes the
        transformation that maps MSI raster coordinates to optical image
        pixel coordinates.
        """
        # Check if reader has FlexImaging-specific metadata
        if not hasattr(self.reader, "mis_metadata"):
            logger.debug("Reader does not have mis_metadata, skipping alignment")
            return

        mis_metadata = getattr(self.reader, "mis_metadata", {})
        areas = mis_metadata.get("areas", [])

        if not areas:
            logger.debug("No Area definitions found, skipping alignment")
            return

        # Get required data for alignment
        positions = getattr(self.reader, "_positions", [])
        header = getattr(self.reader, "_header", {})

        if not positions:
            logger.warning("No position data available for alignment")
            return

        first_raster_x = header.get("first_raster_x", 0)
        first_raster_y = header.get("first_raster_y", 0)

        # Store the primary optical image filename from <ImageFile>
        image_file = mis_metadata.get("ImageFile", "")
        if image_file:
            self.primary_filename = Path(image_file).stem.lower()
            logger.info(f"Primary alignment image from .mis: {image_file}")

        # Compute area-based alignment
        try:
            aligner = TeachingPointAlignment()
            self._alignment = aligner.compute_area_alignment(
                areas=areas,
                poslog_positions=positions,
                first_raster_x=first_raster_x,
                first_raster_y=first_raster_y,
            )
            logger.info(
                f"Computed optical alignment with "
                f"{len(self._alignment.region_mappings)} region mappings"
            )
        except (KeyError, TypeError, ValueError) as e:
            # Thyra arithmetic over Thyra's own parse of the .mis and the
            # poslog (issue #280). What it can legitimately meet is a
            # record those files did not fill: an Area without ``p1``, a
            # position without ``region``, a corner that is not a pair.
            # Not AttributeError -- nothing here reads an attribute, so one
            # would be a defect in the aligner and must keep its traceback.
            logger.warning(f"Failed to compute optical alignment: {e}")
            self._alignment = None

    def build_tic_to_image_affine(self) -> None:
        """Build affine matrix mapping TIC raster-index coords to image pixels.

        When optical alignment is available, this creates a 3x3 affine matrix
        that transforms TIC image coordinates (integer raster indices) into
        optical image pixel coordinates, so the TIC overlays correctly on the
        optical image in SpatialData.

        For single-region data, uses that region's mapping directly.
        For multi-region data, computes a global affine from the overall
        raster bounds and overall image bounds across all regions.

        The matrix encodes: image_pixel = scale * raster_index + offset
        where offset places the first pixel center at image_min + half_pixel.
        """
        if self._alignment is None:
            return
        if not self._alignment.region_mappings:
            return

        mappings = self._alignment.region_mappings

        if len(mappings) == 1:
            # Single region: use its mapping directly
            rm = mappings[0]
            n_raster_x = rm.raster_max_x - rm.raster_min_x + 1
            image_width = rm.image_max_x - rm.image_min_x
            scale_x = image_width / max(1, n_raster_x)
            n_raster_y = rm.raster_max_y - rm.raster_min_y + 1
            image_height = rm.image_max_y - rm.image_min_y
            scale_y = image_height / max(1, n_raster_y)
            half_x = scale_x / 2.0
            half_y = scale_y / 2.0
            tx = rm.image_min_x + half_x
            ty = rm.image_min_y + half_y
        else:
            # Multi-region: compute global affine from overall bounds.
            # The TIC grid covers the full normalized raster space
            # (0..n_x-1, 0..n_y-1). We map this to the bounding box
            # of all region image areas.
            first_rx = self._alignment.first_raster_x
            first_ry = self._alignment.first_raster_y

            # Global raster bounds (original coords)
            global_raster_min_x = min(rm.raster_min_x for rm in mappings)
            global_raster_max_x = max(rm.raster_max_x for rm in mappings)
            global_raster_min_y = min(rm.raster_min_y for rm in mappings)
            global_raster_max_y = max(rm.raster_max_y for rm in mappings)

            # Global image bounds
            global_img_min_x = min(rm.image_min_x for rm in mappings)
            global_img_max_x = max(rm.image_max_x for rm in mappings)
            global_img_min_y = min(rm.image_min_y for rm in mappings)
            global_img_max_y = max(rm.image_max_y for rm in mappings)

            n_raster_x = global_raster_max_x - global_raster_min_x + 1
            n_raster_y = global_raster_max_y - global_raster_min_y + 1
            image_width = global_img_max_x - global_img_min_x
            image_height = global_img_max_y - global_img_min_y

            scale_x = image_width / max(1, n_raster_x)
            scale_y = image_height / max(1, n_raster_y)
            half_x = scale_x / 2.0
            half_y = scale_y / 2.0

            # The TIC grid index (0,0) corresponds to original raster
            # position (first_rx, first_ry). We need to account for
            # any gap between first_rx and global_raster_min_x.
            offset_raster_x = first_rx - global_raster_min_x
            offset_raster_y = first_ry - global_raster_min_y

            tx = global_img_min_x + half_x + offset_raster_x * scale_x
            ty = global_img_min_y + half_y + offset_raster_y * scale_y

        # 3x3 affine: [[sx, 0, tx], [0, sy, ty], [0, 0, 1]]
        self._tic_to_image = np.array(
            [
                [scale_x, 0, tx],
                [0, scale_y, ty],
                [0, 0, 1],
            ],
            dtype=np.float64,
        )
        logger.info(
            f"Built TIC-to-image affine: "
            f"scale=({scale_x:.2f}, {scale_y:.2f}), "
            f"offset=({tx:.1f}, {ty:.1f})"
        )

    def add_images(self, images: MutableMapping[str, Any]) -> None:
        """Load and add optical images from the reader to data structures.

        Finds the optical images associated with the MSI data and adds them
        as image layers in the SpatialData output. The primary alignment image
        (from .mis <ImageFile>) is loaded first so its dimensions are known
        when computing Scale transforms for the other images.

        Args:
            images: The store's images, to add each placeholder to
        """
        if not self._include:
            return

        optical_paths = self.reader.get_optical_image_paths()
        if not optical_paths:
            logger.debug("No optical images found")
            return

        self._adopt_reader_primary_optical()
        logger.info(f"Found {len(optical_paths)} optical image(s)")

        # Load primary image first so we know its dimensions for scaling others
        primary_paths = [p for p in optical_paths if self._is_primary_optical(p)]
        other_paths = [p for p in optical_paths if not self._is_primary_optical(p)]

        for image_path in primary_paths + other_paths:
            try:
                self._load_single_image(image_path, images)
            except Exception as e:
                logger.warning(f"Failed to load optical image {image_path.name}: {e}")

    def _adopt_reader_primary_optical(self) -> None:
        """Take the .mis alignment image from the reader if nothing else set it.

        :meth:`compute_alignment` normally records it while it is
        reading the .mis, but it returns early whenever there is nothing to
        align -- no Area definitions, no positions -- and then the primary
        image is never named even though the .mis names it. The Bruker
        readers resolve it during folder discovery regardless, so ask them.

        Asked for with ``getattr`` even though ``BaseMSIReader`` defines it,
        because a reader here is whatever satisfies the interface and not
        necessarily a subclass: this project's own
        ``tests/unit/converters/test_streaming_converter.py`` passes in a
        plain class that implements the methods and inherits nothing. This
        method is optional and arrived after those readers were written, so
        not having it has to mean "no designated image", not a crash in the
        middle of a conversion.
        """
        if self.primary_filename:
            return
        resolve = getattr(self.reader, "get_primary_optical_image_path", None)
        primary = resolve() if callable(resolve) else None
        if primary is not None:
            self.primary_filename = Path(primary).stem.lower()
            logger.info(f"Primary alignment image from .mis: {primary.name}")

    def _is_primary_optical(self, image_path: Path) -> bool:
        """Check if an image file is the primary alignment image from .mis.

        Matched on the stem, not the whole filename: what the .mis names and
        what is on disk agree on the name but not always on the spelling of
        the extension, and a folder does not hold the same stem twice.
        """
        if not self.primary_filename:
            return False
        return image_path.stem.lower() == self.primary_filename

    def _scale_transform(self, x_size: int, y_size: int) -> Any:
        """Compute a Scale transform for a non-primary optical image.

        Maps the image's pixel coordinates to the primary alignment image's
        coordinate space using the dimension ratio.

        Args:
            x_size: Width of the non-primary image
            y_size: Height of the non-primary image

        Returns:
            Scale transform, or Identity if no primary dimensions available
        """
        if self._primary_dims is None:
            return Identity()

        primary_w, primary_h = self._primary_dims
        scale_x = primary_w / x_size
        scale_y = primary_h / y_size

        logger.info(f"  Scale to primary: ({scale_x:.4f}, {scale_y:.4f})")
        return Scale([scale_x, scale_y], axes=("x", "y"))

    def _to_um_transform(self) -> Any:
        """Build an Affine mapping primary-optical pixels to MSI um.

        Composes the inverse of the tic-to-image affine (so optical
        pixel -> MSI raster index) with the MSI pixel size (so raster
        index -> um).  Used only when ``apply_optical_alignment=False``
        and FlexImaging metadata is available -- it places the optical
        image into the same micrometer "global" frame as the MSI so a
        downstream registration step (e.g. Ousia's EscDat wizard) can
        map both elements together with a single composed affine.

        The math: ``tic_to_image_matrix`` is a 3x3 affine encoding
        ``image_pixel = scale * raster_index + offset``.  Inverting
        and composing with scale-by-pixel-size yields::

            um = pixel_size_um * inv(tic_to_image) @ optical_pixel

        Returns an :class:`Affine` over ``(x, y)`` input + output axes.
        Caller should not invoke when ``tic_to_image`` is None.
        """
        if self._tic_to_image is None:
            raise RuntimeError(
                "_to_um_transform called without "
                "a tic_to_image matrix; check call-site guard."
            )
        inv = np.linalg.inv(self._tic_to_image)
        # Scale matrix: [[ps_x, 0, 0], [0, ps_y, 0], [0, 0, 1]]
        ps_x, ps_y = self._pixel_size_xy()
        scale_mat = np.array(
            [[ps_x, 0.0, 0.0], [0.0, ps_y, 0.0], [0.0, 0.0, 1.0]],
            dtype=np.float64,
        )
        matrix = scale_mat @ inv
        return Affine(
            matrix,
            input_axes=("x", "y"),
            output_axes=("x", "y"),
        )

    def _load_single_image(
        self, image_path: Path, images: MutableMapping[str, Any]
    ) -> None:
        """Load a single optical image and add it to data structures.

        The primary image (identified by .mis <ImageFile>) gets an Identity
        transform. Other images get a Scale transform mapping their pixel
        coordinates to the primary image's coordinate space.

        Args:
            image_path: Path to the optical image (TIFF, JPEG, PNG or BMP)
            images: The store's images, to add the placeholder to
        """
        # Generate a clean name for the image layer
        image_name = self._element_name(image_path)

        logger.info(f"Loading optical image: {image_path.name} as '{image_name}'")

        # Only the page header is read here. The pixels never enter this
        # process whole: the element is declared to SpatialData as a lazy
        # placeholder with the final shape, dtype, chunking and pyramid,
        # and stream_pending_pixels() fills it in bands once the store
        # exists. See this module's own docstring for why (the whole-page
        # route cost ~5x the decoded image in transient memory).
        # probe raises for a layout or sample format it cannot read; the
        # per-image guard in add_images turns that into the same
        # "skip with a warning" the whole-page decode used to give.
        probed = probe_optical_source(image_path)
        source: AnyOpticalSource = probed
        n_channels, y_size, x_size = source.shape

        # Determine transform.  Two cases:
        #
        # 1. apply_optical_alignment=True (default): "global" is
        #    optical-image pixel space.  Primary image is Identity;
        #    non-primary images Scale to match primary dims.
        #
        # 2. apply_optical_alignment=False (e.g. Ousia wizard):
        #    "global" is MSI micrometer space.  Map the primary
        #    image's pixel coordinates into MSI um using the inverse
        #    of the tic-to-image affine, then scale by pixel_size_um.
        #    This way the optical image lands alongside the MSI in
        #    the same um frame and downstream registration steps
        #    map both together.
        um_mode = not self._apply_alignment and self._tic_to_image is not None
        is_primary = self._is_primary_optical(image_path)
        crop: Optional[OpticalCrop] = None
        if is_primary:
            # The FULL dimensions, crop or not: the other images scale into
            # the whole scan's pixel grid, which is the frame everything on
            # the MSI side is stated in.
            self._primary_dims = (x_size, y_size)
            if um_mode:
                transform = self._to_um_transform()
                logger.info(
                    f"  Primary image -> um via inverse alignment: {x_size}x{y_size}"
                )
            else:
                transform = Identity()
                logger.info(f"  Primary alignment image: {x_size}x{y_size}")
            crop = self._crop_for_alignment_image(x_size, y_size)
            if crop is not None:
                # The crop's pixel (u, v) is the scan's (x0 + u, y0 + v):
                # that translation goes FIRST, and then whatever carried the
                # whole scan into the store's frame carries the crop. So
                # the raster-to-image affine, the pixel polygons and the
                # coordinate_systems attr keep describing the scan's grid.
                shift = Translation([crop.x0, crop.y0], axes=("x", "y"))
                transform = (
                    shift
                    if isinstance(transform, Identity)
                    else SequenceTransform([shift, transform])
                )
                source = CroppedOpticalSource(probed, crop)
                n_channels, y_size, x_size = source.shape
                self._crop = crop
                logger.info(
                    f"  Cropped to region {crop.region_id} ('{crop.region_name}'): "
                    f"{crop.width}x{crop.height} at ({crop.x0}, {crop.y0}) of "
                    f"{crop.full_width}x{crop.full_height}"
                )
        elif self._primary_dims is not None:
            # Non-primary: first scale to match primary, then if
            # we're in um-mode, chain through the same um affine.
            base = self._scale_transform(x_size, y_size)
            if um_mode:
                transform = SequenceTransform([base, self._to_um_transform()])
            else:
                transform = base
        else:
            transform = Identity()

        # Multi-scale pyramid + chunked layout.
        #
        # Without scale_factors a single-scale image is written and any
        # downstream viewer has to read full-resolution tiles at every
        # zoom level.  For a typical FlexImaging brightfield (10k x 10k+
        # pixels) that is the difference between an instant first paint
        # and a multi-second stall every time the user pans or zooms.
        #
        # We mirror what spatialdata-io's xenium reader does for its
        # morphology images: scale_factors=[2, 2, 2, 2] gives the viewer
        # five pyramid levels.  Here we adapt the level count to the
        # image's smallest spatial dimension so tiny images don't waste
        # levels and huge ones get enough to keep the coarsest level
        # fast (< ~1000 px short side).
        #
        # chunks=(1, 4096, 4096) stores each channel as 4k x 4k blocks
        # so a viewer's 512 x 512 tile read decompresses at most one
        # chunk per request.
        smallest = min(y_size, x_size)
        scale_factors = _calc_optical_scale_factors(smallest)
        streamed = StreamedOpticalImage(
            source=source,
            name=image_name,
            chunks=image_chunks(2),  # (1, 4096, 4096); sharding seam, see _chunking
            scale_factors=scale_factors,
            transformations={
                self.dataset_id: transform,
                "global": transform,
            },
        )
        # Keyed by element name, as the images dict is: a second file that
        # maps to the same name replaces the first, the way the dict
        # assignment always did, only now with a warning.
        earlier = self._pending.get(image_name)
        if earlier is not None:
            logger.warning(
                f"Optical image '{image_name}' from {earlier.source.path.name} "
                f"is replaced by {image_path.name}, which maps to the same name"
            )
        images[image_name] = streamed.placeholder()
        self._pending[image_name] = streamed
        # Which file this element came from, and whether it is the one the
        # alignment is stated against. Neither is recoverable from the
        # store otherwise: the element name drops the extension and
        # rewrites the stem (_0000 -> highres, and a stem over 30
        # characters is truncated), and the alignment image is only
        # distinguishable by its transform, and only when the alignment
        # was applied. Recorded here, written by the root attrs.
        self._sources[image_name] = image_path.name
        if is_primary:
            self._alignment_element = image_name

        pyramid_desc = (
            f", {len(scale_factors)} pyramid level{'s' if len(scale_factors) != 1 else ''}"
            if scale_factors
            else " (no pyramid; image small enough)"
        )
        logger.info(
            f"Added optical image '{image_name}': {x_size}x{y_size} "
            f"({n_channels} channel{'s' if n_channels > 1 else ''}){pyramid_desc}"
            "; pixels stream in once the store is written"
        )

    def stream_pending_pixels(self) -> int:
        """Fill every optical image declared so far with its pixels.

        Call once the SpatialData write that carried the placeholders has
        returned and before metadata is consolidated. Each image streams
        from its TIFF in bands and builds its pyramid level by level on
        disk, so memory stays bounded by one band, not by the image.

        A TIFF whose pixels cannot be read is dropped from the store with a
        warning and the conversion goes on without it -- the tolerance the
        whole-page decode had, when the same failure happened before
        anything was written. Only a failure to drop the element propagates,
        because an image with metadata and no pixels is a corrupt store.

        Returns:
            The number of images whose pixels are now in the store.
        """
        pending, self._pending = self._pending, {}
        streamed = 0
        for image in pending.values():
            logger.info(f"Streaming optical image pixels: '{image.name}'")
            try:
                image.stream_pixels(self.output_path)
            except Exception as e:  # mirrors the per-image guard in add_images
                logger.warning(
                    f"Failed to load optical image {image.source.path.name}: {e}; "
                    f"dropping '{image.name}' from the store"
                )
                image.discard(self.output_path)
                self.forget_image(image.name)
                continue
            streamed += 1
        return streamed

    def forget_image(self, name: str) -> None:
        """Take a dropped optical image back out of the store's root attrs.

        The root attrs were composed and written by the ``SpatialData``
        write that carried the placeholders, so an image
        :meth:`stream_pending_pixels` then drops is still named in
        them. Metadata naming an element that is not in the store is worse
        than none -- a consumer that trusts it gets a KeyError where it
        would otherwise have fallen back -- so the store is corrected
        before ``zarr.consolidate_metadata`` runs, which copies whatever is
        here into the consolidated document.

        Args:
            name: Element name of the image whose pixels could not be read.
        """
        self._sources.pop(name, None)
        was_alignment = self._alignment_element == name
        if was_alignment:
            self._alignment_element = None
        try:
            root = zarr.open_group(
                str(self.output_path), mode="r+", use_consolidated=False
            )
            optical = root.attrs.get("optical_images")
            if isinstance(optical, dict):
                remaining = self.root_attr()
                if remaining is None:
                    del root.attrs["optical_images"]
                else:
                    root.attrs["optical_images"] = remaining
            if was_alignment:
                self._clear_reference_element(root)
        except Exception as e:  # pragma: no cover - a store we just wrote
            logger.warning(
                f"Could not unrecord the dropped optical image '{name}' from "
                f"the store's attrs: {e}"
            )

    @staticmethod
    def _clear_reference_element(root: Any) -> None:
        """Null ``coordinate_systems.global.reference_element`` in ``root``.

        Reassigns the whole attr: a zarr attribute is a value, so mutating
        the dict a read returns changes nothing on disk.
        """
        systems = root.attrs.get("coordinate_systems")
        if not isinstance(systems, dict) or not isinstance(systems.get("global"), dict):
            return
        root.attrs["coordinate_systems"] = {
            **systems,
            "global": {**systems["global"], "reference_element": None},
        }

    def _element_name(self, image_path: Path) -> str:
        """Generate a clean name for an optical image layer.

        The suffix is dropped, so the same acquisition exported as a .tif or
        a .jpg lands under the same element name.

        Args:
            image_path: Path to the optical image

        Returns:
            Clean name for the image layer (e.g., 'optical_0000', 'optical_deriv')
        """
        stem = image_path.stem.lower()

        # Extract meaningful suffix from filename
        if "_0000" in stem:
            suffix = "highres"
        elif "_0001" in stem:
            suffix = "derived"
        elif "deriv" in stem:
            suffix = "overview"
        else:
            # Use stem with special chars replaced
            suffix = stem.replace(" ", "_").replace("-", "_")
            # Truncate if too long
            if len(suffix) > 30:
                suffix = suffix[:30]

        return f"{self.dataset_id}_optical_{suffix}"

    def root_attr(self) -> Optional[Dict[str, Any]]:
        """What each optical element in this store came from, and which aligns.

        Two facts the store could not state before:

        * **Which file.** The element name is derived, not the filename:
          :meth:`_element_name` drops the extension, maps
          ``_0000``/``_0001``/``deriv`` onto ``highres``/``derived``/
          ``overview`` and truncates anything else at 30 characters, so
          ``sample_0000.tif`` and ``sample_0000.jpg`` both land under
          ``<dataset_id>_optical_highres`` and neither name survives.
        * **Which one the alignment is stated against.** The .mis names it
          in ``<ImageFile>``; the store only ever implied it, through the
          transform (the alignment image gets ``Identity``, the others a
          ``Scale`` into its pixel grid) and only when
          ``apply_optical_alignment=True``. A consumer that wanted the
          alignment image had to guess -- Ousia guesses alphabetically.

        Written into the store's own root attrs rather than onto the
        elements, because per-element attributes do not survive the write:
        see :class:`StreamedOpticalImage`.

        * **Whether it is a crop, and of what.** The alignment element of a
          single-region conversion is a window of the whole-slide scan,
          and its element carries the translation that puts the window
          back. ``crop`` records the window's origin and size and the
          scan's full size, in the scan's pixels, so the whole slide can
          be located -- or a "show the whole slide" built -- without
          re-reading the ``.d``. Absent on an element that is the whole
          file.

        Returns:
            ``{"alignment_element": str | None, "elements": {name:
            {"source_file": str, "crop"?: {...}}}}``, or ``None`` when the
            conversion put no optical image in the store -- the section is
            omitted rather than written empty, as every other optional
            section here is.
        """
        if not self._sources:
            return None
        elements: Dict[str, Dict[str, Any]] = {
            name: {"source_file": source_file}
            for name, source_file in self._sources.items()
        }
        if self._crop is not None and self._alignment_element in elements:
            elements[self._alignment_element]["crop"] = self._crop.to_attr()
        return {
            "alignment_element": self._alignment_element,
            "elements": elements,
        }
