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
picks between the two by the file's first bytes: a TIFF is read by
tifffile whatever it is named.

**What is stored is what the file shows** (D69). A TIFF's colour model is
applied (palette, inverted grey, CMYK and YCbCr become grey or RGB), a TIFF
whose channels span pages is read as one image, a 16-bit colour PNG keeps
its 16 bits, and an orientation tag is applied by
:class:`OrientedOpticalSource`, band by band. The alignment image and its
resampled copies are the exception to the last rule: they stay as encoded.

The pixel values are those ``xarray``'s ``coarsen`` produces: same
``boundary="trim", side="right"`` rule (an odd-sized axis drops its FIRST
row or column), same float64 mean, same truncating cast back to the source
dtype. :func:`block_mean` is unit-tested against xarray for exactly that.

**Failure policy.** An optical image that cannot be decoded used to be
skipped with a warning before anything was written; now the header is read
up front and the pixels only after the store exists, so
:meth:`OpticalImages.stream_pending_pixels` keeps that
contract by dropping the declared element from the store and warning. A
store never keeps an image whose pixels were not written. A PNG whose
header declares more pixel data than the file can hold, and a whole-page
decode larger than the free memory, are refused at the header, before the
store declares anything.

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
from dataclasses import dataclass, field, replace
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
import imagecodecs
import numpy as np
import psutil
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
from .acquisition_frame import AcquisitionOrigin, element_transformations

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ...alignment.teaching_points import RegionMapping
    from ...core.base_reader import BaseMSIReader, OpticalImageLabel

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

#: Suffixes :func:`probe_optical_source` routes to :class:`OpticalTiffSource`
#: even when the file does not start like a TIFF, so tifffile names what is
#: wrong with it. Everything else goes by its first bytes: see
#: :data:`TIFF_SIGNATURES`.
TIFF_SUFFIXES = (".tif", ".tiff")

#: The first four bytes of a TIFF and a BigTIFF, in either byte order. A file
#: that starts with one is read by tifffile whatever its suffix, so the same
#: bytes give the same store under ``.tif``, ``.png`` or an archive's declared
#: media type. Pillow would read such a file too, but turns it by its
#: orientation tag on some compressions and not on others.
TIFF_SIGNATURES = (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")

#: The EXIF and TIFF tag that says how the stored pixels are turned.
ORIENTATION_TAG = 274

#: Orientation tag values that read the stored rows bottom to top, and those
#: that swap rows and columns. See :func:`orient`.
_ROWS_REVERSED = (3, 4, 6, 7)
_TRANSPOSED = (5, 6, 7, 8)

#: Leading axes of a TIFF series that are channels of one picture: tifffile's
#: channel and sample axes, and its codes for a run of pages it has no other
#: name for. Any other leading axis (a z stack, a time series) keeps its first
#: plane only, with a warning.
_CHANNEL_AXES = "CSQI"

#: The most a deflate stream can expand: 258 bytes from a 2-bit code. A PNG
#: whose header declares more decoded bytes than this times its file size
#: cannot hold them, and is refused before anything is allocated.
DEFLATE_MAX_RATIO = 1032

#: Samples per pixel of each PNG colour type (grey, RGB, palette, grey with
#: alpha, RGB with alpha).
_PNG_SAMPLES = {0: 1, 2: 3, 3: 1, 4: 2, 6: 4}

#: Pillow's raw modes for PNGs with 16 bits a sample that Pillow can only
#: narrow to 8 bits, each with the channel count the file holds. These are
#: decoded by imagecodecs at full depth instead. A 16-bit grey PNG has a
#: 16-bit Pillow mode and stays with Pillow.
_PNG_WIDE_RAWMODES = {"RGB;16B": 3, "RGBA;16B": 4, "LA;16B": 2}

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


def orient(cyx: np.ndarray, tag: int) -> np.ndarray:
    """Turn a ``(c, y, x)`` array upright as its orientation tag says.

    ``tag`` is the EXIF/TIFF orientation (tag 274): 1 is upright, 2 to 4
    mirror or turn by a half turn, 5 to 8 swap rows and columns. The result
    equals Pillow's ``ImageOps.exif_transpose`` of the same pixels. The
    operations are separable, so a row band of the stored image comes out
    as the matching row or column band of the upright one; see
    :class:`OrientedOpticalSource`.
    """
    if tag == 2:
        turned = cyx[:, :, ::-1]
    elif tag == 3:
        turned = cyx[:, ::-1, ::-1]
    elif tag == 4:
        turned = cyx[:, ::-1, :]
    elif tag == 5:
        turned = cyx.transpose(0, 2, 1)
    elif tag == 6:
        turned = np.rot90(cyx, -1, axes=(1, 2))
    elif tag == 7:
        turned = cyx[:, ::-1, ::-1].transpose(0, 2, 1)
    elif tag == 8:
        turned = np.rot90(cyx, 1, axes=(1, 2))
    else:
        return cyx
    return np.ascontiguousarray(turned)


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


def _multiscale_dataset_paths(attrs: Mapping[str, Any]) -> Optional[List[str]]:
    """The ``path`` of each level in an OME-Zarr image's ``multiscales[0]``.

    Returns None when the metadata is missing or not shaped as OME-Zarr
    writes it: a dict, a non-empty list, a dict, a list of dicts that
    each name a string path.
    """
    ome = attrs.get("ome")
    if not isinstance(ome, dict):
        return None
    multiscales = ome.get("multiscales")
    if not isinstance(multiscales, list) or not multiscales:
        return None
    first = multiscales[0]
    if not isinstance(first, dict):
        return None
    datasets = first.get("datasets")
    if not isinstance(datasets, list):
        return None
    paths: List[str] = []
    for dataset in datasets:
        path = dataset.get("path") if isinstance(dataset, dict) else None
        if not isinstance(path, str):
            return None
        paths.append(path)
    return paths


def _tag_numbers(value: Any, count: int) -> Optional[Tuple[float, ...]]:
    """``count`` numbers from a TIFF tag value, rationals divided out.

    tifffile returns a RATIONAL array as numerator, denominator pairs; this
    takes either that or the plain numbers. ``None`` when the value has
    neither length or holds a zero denominator.
    """
    numbers = tuple(float(n) for n in np.ravel(np.asarray(value, dtype=np.float64)))
    if len(numbers) == count:
        return numbers
    if len(numbers) == 2 * count and all(numbers[1::2]):
        return tuple(n / d for n, d in zip(numbers[0::2], numbers[1::2]))
    return None


def _rows_at_a_time(rows: int, row_values: int) -> int:
    """Rows whose float64 working copy stays near 64 MB, at least one."""
    return max(1, (64 * 1024 * 1024) // max(1, 8 * row_values))


@dataclass(frozen=True)
class TiffColour:
    """How a TIFF page's decoded samples become the colour that is stored.

    tifffile decodes samples as the file stores them and leaves the colour
    model (the photometric tag) to the caller. ``model`` is what this does
    with them: ``"as_is"`` for grey and RGB, and for JPEG-compressed YCbCr,
    which tifffile already decodes to RGB; ``"palette"`` looks each index up
    in the colour map, giving RGB; ``"miniswhite"`` inverts grey, so 0 is
    black as everywhere else in the store; ``"cmyk"`` turns ink into RGB the
    way Pillow does; ``"ycbcr"`` turns uncompressed luma and chroma into RGB
    with the file's own coefficients. Extra samples (alpha) after the colour
    samples are kept as they are.
    """

    model: str = "as_is"
    colormap: Optional[np.ndarray] = field(default=None, compare=False)
    white: int = 0
    colour_samples: int = 0
    ycbcr: Tuple[float, ...] = ()

    @classmethod
    def of(cls, page: Any, name: str) -> "TiffColour":
        """Read the colour model of a TIFF page.

        Raises:
            ValueError: for a colour model this cannot turn into grey or
                RGB (CIELab, a transparency mask, multi-ink separations,
                subsampled uncompressed YCbCr, inverted signed or floating
                samples). The image is then skipped with a warning, like
                any other image Thyra cannot read.
        """
        photometric = tifffile.PHOTOMETRIC(int(page.photometric))
        if photometric in (tifffile.PHOTOMETRIC.MINISBLACK, tifffile.PHOTOMETRIC.RGB):
            return cls()
        readers = {
            tifffile.PHOTOMETRIC.PALETTE: cls._palette,
            tifffile.PHOTOMETRIC.MINISWHITE: cls._inverted_grey,
            tifffile.PHOTOMETRIC.SEPARATED: cls._separated,
            tifffile.PHOTOMETRIC.YCBCR: cls._luma_chroma,
        }
        read = readers.get(photometric)
        if read is None:
            raise ValueError(
                f"{name}: colour model {photometric.name} is not read; "
                "only grey, RGB, palette, CMYK and YCbCr are"
            )
        return read(page, name)

    @staticmethod
    def _layout(page: Any) -> Tuple[int, int, np.dtype, int]:
        """``(samples, colour samples, dtype, bits)`` of a TIFF page."""
        samples = int(page.samplesperpixel)
        colour = samples - len(page.extrasamples or ())
        return samples, colour, np.dtype(page.dtype), int(page.bitspersample)

    @classmethod
    def _palette(cls, page: Any, name: str) -> "TiffColour":
        samples, _, dtype, _ = cls._layout(page)
        colormap = page.colormap
        if colormap is None or samples != 1 or dtype.kind != "u":
            raise ValueError(f"{name}: palette image without a usable colour map")
        return cls(model="palette", colormap=cls._narrow_colormap(colormap))

    @classmethod
    def _inverted_grey(cls, page: Any, name: str) -> "TiffColour":
        _, colour, dtype, bits = cls._layout(page)
        if dtype.kind == "b":
            return cls(model="miniswhite", white=1, colour_samples=colour)
        if dtype.kind != "u":
            raise ValueError(f"{name}: inverted grey with {dtype} samples")
        return cls(model="miniswhite", white=2**bits - 1, colour_samples=colour)

    @classmethod
    def _separated(cls, page: Any, name: str) -> "TiffColour":
        _, colour, dtype, bits = cls._layout(page)
        inkset = page.tags.valueof(332, 1)
        if inkset != 1 or colour != 4 or dtype.kind != "u":
            raise ValueError(
                f"{name}: separated colour with {colour} inks (inkset {inkset}); "
                "only CMYK is read"
            )
        return cls(model="cmyk", white=2**bits - 1, colour_samples=4)

    @classmethod
    def _luma_chroma(cls, page: Any, name: str) -> "TiffColour":
        if int(page.compression) in (6, 7):
            # OJPEG and JPEG: tifffile's decoder already returns RGB.
            return cls()
        _, colour, dtype, bits = cls._layout(page)
        # TIFF 6.0's default for an absent YCbCrSubSampling is (2, 2).
        if page.subsampling is None or tuple(page.subsampling) != (1, 1):
            raise ValueError(f"{name}: subsampled uncompressed YCbCr is not read")
        if colour != 3 or dtype.kind != "u":
            raise ValueError(f"{name}: YCbCr with {colour} {dtype} samples")
        white = float(2**bits - 1)
        half = float(2 ** (bits - 1))
        luma = _tag_numbers(page.tags.valueof(529), 3) or (0.299, 0.587, 0.114)
        reference = _tag_numbers(page.tags.valueof(532), 6) or (
            0.0,
            white,
            half,
            white,
            half,
            white,
        )
        return cls(
            model="ycbcr", white=int(white), colour_samples=3, ycbcr=luma + reference
        )

    @staticmethod
    def _narrow_colormap(colormap: Any) -> np.ndarray:
        """The colour map as 8-bit when it holds 8-bit colours, else 16-bit.

        TIFF stores colour-map entries as 16-bit numbers. Writers fill them
        from 8-bit colours as ``v * 257`` or ``v * 256``, and some write the
        8-bit value itself; each of those is stored as 8-bit. A map that
        uses the full 16 bits keeps them.
        """
        table = np.asarray(colormap).astype(np.uint16)
        if int(table.max(initial=0)) < 256:
            return table.astype(np.uint8)
        low, high = table & 0xFF, table >> 8
        if not low.any() or np.array_equal(low, high):
            return high.astype(np.uint8)
        return table

    def channels(self, samples: int) -> int:
        """Stored channels for a page with ``samples`` samples per pixel."""
        if self.model == "palette":
            return 3
        if self.model == "cmyk":
            return samples - 1
        return samples

    def dtype(self, decoded: np.dtype) -> np.dtype:
        """Stored dtype for samples decoded as ``decoded``."""
        if self.model == "palette" and self.colormap is not None:
            return self.colormap.dtype
        return decoded

    def apply(self, cyx: np.ndarray) -> np.ndarray:
        """Decoded ``(c, y, x)`` samples to the stored colour."""
        if self.model == "palette" and self.colormap is not None:
            return np.take(self.colormap, cyx[0], axis=1, mode="clip")
        if self.model == "miniswhite":
            out = np.array(cyx, copy=True)
            colour = out[: self.colour_samples]
            if out.dtype.kind == "b":
                np.logical_not(colour, out=colour)
            else:
                np.subtract(self.white, colour, out=colour, casting="unsafe")
            return out
        if self.model in ("cmyk", "ycbcr"):
            n_out = 3 + cyx.shape[0] - self.colour_samples
            out = np.empty((n_out,) + cyx.shape[1:], dtype=cyx.dtype)
            out[3:] = cyx[self.colour_samples :]
            step = _rows_at_a_time(cyx.shape[1], cyx.shape[0] * cyx.shape[2])
            convert = self._cmyk if self.model == "cmyk" else self._ycbcr
            for start in range(0, cyx.shape[1], step):
                rows = slice(start, start + step)
                out[:3, rows] = convert(cyx[:, rows].astype(np.float64))
            return out
        return cyx

    def _cmyk(self, ink: np.ndarray) -> np.ndarray:
        """Pillow's CMYK to RGB: each colour is ``(W - ink) * (W - K) / W``."""
        white = float(self.white)
        rgb = (white - ink[:3]) * (white - ink[3]) / white
        return np.floor(rgb + 0.5)

    def _ycbcr(self, samples: np.ndarray) -> np.ndarray:
        """TIFF 6.0 YCbCr to RGB, from the file's coefficients and reference."""
        red, green, blue = self.ycbcr[:3]
        ref = self.ycbcr[3:]
        white = float(self.white)
        coding = float(2 ** (int(round(math.log2(self.white + 1))) - 1) - 1)
        luma = (samples[0] - ref[0]) * white / (ref[1] - ref[0])
        cb = (samples[1] - ref[2]) * coding / (ref[3] - ref[2])
        cr = (samples[2] - ref[4]) * coding / (ref[5] - ref[4])
        r = cr * (2 - 2 * red) + luma
        b = cb * (2 - 2 * blue) + luma
        g = (luma - blue * b - red * r) / green
        return np.clip(np.floor(np.stack([r, g, b]) + 0.5), 0, white)


def _tiff_orientation(page: Any) -> Optional[int]:
    """The page's orientation tag, or ``None`` when it carries none."""
    tag = page.tags.get(ORIENTATION_TAG)
    if tag is None:
        return None
    try:
        return int(np.ravel(tag.value)[0])
    except (TypeError, ValueError, IndexError):
        return None


@dataclass(frozen=True)
class OpticalTiffSource:
    """What an optical TIFF looks like, without decoding it.

    ``shape`` is always ``(c, y, x)``. ``page_axes`` is one of
    :data:`PAGE_LAYOUTS`, or ``"series"`` for an image whose channels are
    spread over several pages (an OME-TIFF), and says how a decoded array
    maps onto that; ``strip_rows`` is the height of the page's strips or
    tiles, the smallest row range the decoder can produce on its own.
    ``colour`` turns the decoded samples into the stored colour.
    ``orientation`` is the first page's orientation tag; this source yields
    the pixels as stored, and :class:`OrientedOpticalSource` turns them.

    For a series, ``series_index`` holds ``None`` for each leading axis that
    is kept as channels and ``0`` for each that keeps its first plane only;
    ``dropped_planes`` counts the pages left out that way, along
    ``dropped_axes``, so the caller can say so.
    """

    path: Path
    shape: Shape
    dtype: np.dtype
    page_axes: str
    strip_rows: int
    colour: TiffColour = field(default_factory=TiffColour)
    orientation: Optional[int] = None
    series_index: Tuple[Optional[int], ...] = ()
    series_samples_last: bool = False
    dropped_planes: int = 0
    dropped_axes: str = ""

    @classmethod
    def probe(cls, path: Union[str, Path]) -> "OpticalTiffSource":
        """Read the page header, and the first series when it spans pages.

        Raises:
            ValueError: for a page layout other than :data:`PAGE_LAYOUTS`,
                a sample format tifffile cannot map to a dtype or a colour
                model :class:`TiffColour` cannot read -- the cases the
                whole-page decode used to reject, at the same point, before
                anything is written.
        """
        path = Path(path)
        with tifffile.TiffFile(path) as tif:
            page: Any = tif.pages[0]
            page_shape = tuple(int(n) for n in page.shape)
            dtype = page.dtype
            axes = str(page.axes)
            chunks = tuple(int(n) for n in page.chunks)
            if dtype is None:
                raise ValueError(
                    f"{path.name}: sample format not supported by tifffile"
                )
            colour = TiffColour.of(page, path.name)
            orientation = _tiff_orientation(page)
            series = tif.series[0] if tif.series else None
            pages = len(series.pages) if series is not None else 1
            spread = series is not None and pages > 1
            series_axes = str(series.axes) if spread and series is not None else ""
            series_shape = (
                tuple(int(n) for n in series.shape)
                if spread and series is not None
                else ()
            )
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
        source = cls(
            path=path,
            shape=(colour.channels(shape[0]), shape[1], shape[2]),
            dtype=colour.dtype(np.dtype(dtype)),
            page_axes=axes,
            strip_rows=max(1, min(strip_rows, shape[1])),
            colour=colour,
            orientation=orientation,
        )
        if spread:
            return source._with_series(series_axes, series_shape, pages)
        return source

    def _with_series(
        self, axes: str, shape: Tuple[int, ...], pages: int
    ) -> "OpticalTiffSource":
        """This source over the whole first series, when its pages are channels.

        Each page must be one plane of the same picture: grey or RGB, laid
        out ``YX`` or ``YXS``. The leading axes in :data:`_CHANNEL_AXES`
        become channels, page after page; any other leading axis keeps its
        first plane. When only one plane is kept, or the pages cannot be
        stacked, the first page is read on its own and the rest of the
        series' ``pages`` are counted as left out.
        """
        dims = [int(n) for n in shape]
        position = axes.find("YX")
        trailing = axes[position:] if position >= 0 else ""
        spatial = dims[position : position + 2] if position >= 0 else []
        if (
            self.colour.model != "as_is"
            or self.page_axes not in ("YX", "YXS")
            or trailing not in ("YX", "YXS")
            or tuple(spatial) != tuple(self.shape[1:])
        ):
            return replace(self, dropped_planes=pages - 1, dropped_axes=axes)
        leading = axes[:position]
        sizes = dims[:position]
        index = tuple(None if a in _CHANNEL_AXES else 0 for a in leading)
        kept = int(np.prod([n for a, n in zip(leading, sizes) if a in _CHANNEL_AXES]))
        dropped = int(np.prod(sizes)) - kept
        dropped_axes = "".join(a for a in leading if a not in _CHANNEL_AXES)
        if kept <= 1:
            return replace(self, dropped_planes=dropped, dropped_axes=dropped_axes)
        samples = dims[-1] if trailing == "YXS" else 1
        return replace(
            self,
            shape=(kept * samples, spatial[0], spatial[1]),
            page_axes="series",
            series_index=index,
            series_samples_last=trailing == "YXS",
            dropped_planes=dropped,
            dropped_axes=dropped_axes,
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
        A series is always read through the adapter, over all its pages.

        ``window`` is ``(y_start, y_stop, x_start, x_stop)`` in this page's
        pixels: only those rows are decoded, only those columns are kept,
        and ``first_row`` stays in the page's own row numbering. With a
        window the adapter is always used, so a single-strip page decodes
        its one strip and keeps the window of it.
        """
        n_rows = self.shape[1]
        with tifffile.TiffFile(self.path) as tif:
            if self.page_axes == "series":
                store = tif.series[0].aszarr(level=0)
            else:
                page: Any = tif.pages[0]
                if window is None and rows >= n_rows:
                    yield 0, self.colour.apply(self._to_cyx(page.asarray()))
                    return
                store = page.aszarr()
            y_start, y_stop, x_start, x_stop = window or (0, n_rows, 0, self.shape[2])
            source = zarr.open_array(store=store, mode="r")
            for start in range(y_start, y_stop, rows):
                stop = min(start + rows, y_stop)
                # zarr types reads as array-or-scalar; a slice of a
                # numpy-buffered array is an ndarray, so asarray copies nothing.
                if self.page_axes == "series":
                    decoded = np.asarray(
                        source[self._series_selection(start, stop, x_start, x_stop)]
                    )
                elif self.page_axes == "SYX":
                    decoded = np.asarray(source[:, start:stop, x_start:x_stop])
                else:
                    decoded = np.asarray(source[start:stop, x_start:x_stop])
                yield start, self.colour.apply(self._to_cyx(decoded))

    def _series_selection(
        self, start: int, stop: int, x_start: int, x_stop: int
    ) -> Tuple[Union[int, slice], ...]:
        """The zarr selection of a row band across the kept planes."""
        leading: Tuple[Union[int, slice], ...] = tuple(
            slice(None) if i is None else i for i in self.series_index
        )
        selection = leading + (slice(start, stop), slice(x_start, x_stop))
        if self.series_samples_last:
            selection += (slice(None),)
        return selection

    def _to_cyx(self, decoded: np.ndarray) -> np.ndarray:
        if self.page_axes == "YX":
            return decoded[np.newaxis, :, :]
        if self.page_axes == "YXS":
            return np.moveaxis(decoded, -1, 0)
        if self.page_axes == "series":
            if self.series_samples_last:
                decoded = np.moveaxis(decoded, -1, -3)
            return decoded.reshape((-1,) + decoded.shape[-2:])
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


def _png_decoded_bytes(path: Path) -> Optional[int]:
    """Bytes the IDAT stream of a PNG inflates to, read from its header.

    Each row is one filter byte plus its packed samples. ``None`` when the
    file does not start with a PNG signature and an IHDR chunk.
    """
    with path.open("rb") as handle:
        head = handle.read(33)
    if len(head) < 33 or head[:8] != b"\x89PNG\r\n\x1a\n" or head[12:16] != b"IHDR":
        return None
    width = int.from_bytes(head[16:20], "big")
    height = int.from_bytes(head[20:24], "big")
    depth, colour_type = head[24], head[25]
    bits = depth * _PNG_SAMPLES.get(colour_type, 4)
    return height * (1 + (width * bits + 7) // 8)


def _check_decodable(path: Path, decoded_bytes: int, copies: int) -> None:
    """Refuse an image its own file, or this machine, cannot decode.

    A PNG that declares more decoded bytes than :data:`DEFLATE_MAX_RATIO`
    times its file size is damaged or forged: no deflate stream expands
    that far. Any whole-page decode that needs more memory than is free is
    refused too, before the store declares the image, rather than failing
    part-way with an empty reason after minutes of swapping.

    Raises:
        ValueError: naming the declared size and what it does not fit.
    """
    inflated = _png_decoded_bytes(path)
    size = path.stat().st_size
    if inflated is not None and inflated > DEFLATE_MAX_RATIO * size:
        raise ValueError(
            f"{path.name}: its header declares {inflated / 1e9:.2f} GB of "
            f"pixel data, more than its {size} bytes can hold"
        )
    needed = decoded_bytes * copies
    available = int(psutil.virtual_memory().available)
    if needed > available:
        raise ValueError(
            f"{path.name}: decoding it needs about {needed / 1e9:.1f} GB of "
            f"memory and {available / 1e9:.1f} GB is free"
        )


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
    ``wide_png`` marks a PNG with 16 bits a colour sample, which Pillow can
    only narrow to 8 bits; imagecodecs decodes it at full depth instead.
    ``orientation`` is the EXIF orientation tag, read but not applied here:
    see :class:`OrientedOpticalSource`.
    """

    path: Path
    shape: Shape
    dtype: np.dtype
    mode: str
    convert_to: Optional[str]
    strip_rows: int
    wide_png: bool = False
    orientation: Optional[int] = None

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
            ValueError: if the declared size cannot be decoded; see
                :func:`_check_decodable`.
        """
        path = Path(path)
        with _unlimited_pixels(), PILImage.open(path) as page:
            mode = str(page.mode)
            width, height = (int(n) for n in page.size)
            # Pillow >= 9.5. Says whether the file carries alpha at all,
            # including the palette-with-transparency case that plain
            # ``"A" in mode`` misses.
            transparent = bool(getattr(page, "has_transparency_data", False))
            rawmode = cls._png_rawmode(page)
            orientation = cls._exif_orientation(page)
        if width < 1 or height < 1:
            raise ValueError(f"{path.name}: image is {width}x{height}")
        convert_to: Optional[str] = None
        if rawmode in _PNG_WIDE_RAWMODES:
            channels, dtype = _PNG_WIDE_RAWMODES[rawmode], "uint16"
        else:
            spec = _PILLOW_MODES.get(mode)
            if spec is None:
                convert_to = "RGBA" if transparent else "RGB"
                spec = _PILLOW_MODES[convert_to]
            channels, dtype = spec
        decoded = channels * height * width * np.dtype(dtype).itemsize
        _check_decodable(path, decoded, 3 if convert_to else 2)
        return cls(
            path=path,
            shape=(channels, height, width),
            dtype=np.dtype(dtype),
            mode=mode,
            convert_to=convert_to,
            strip_rows=height,
            wide_png=rawmode in _PNG_WIDE_RAWMODES,
            orientation=int(orientation) if isinstance(orientation, int) else None,
        )

    @staticmethod
    def _exif_orientation(page: Any) -> Any:
        """The EXIF orientation from the header alone, or ``None``.

        Read from the EXIF block the header parse already holds (a JPEG's
        APP1, a PNG's eXIf chunk), not through ``getexif()``: on a PNG that
        decodes the whole image to look for a chunk after the pixels.
        """
        raw = page.info.get("exif")
        if not isinstance(raw, (bytes, bytearray)) or not raw:
            return None
        exif = PILImage.Exif()
        try:
            exif.load(bytes(raw))
        except Exception:  # a damaged EXIF block only loses the tag
            return None
        return exif.get(ORIENTATION_TAG)

    @staticmethod
    def _png_rawmode(page: Any) -> Optional[str]:
        """The raw mode Pillow decodes a PNG's samples from, else ``None``."""
        if getattr(page, "format", None) != "PNG" or not page.tile:
            return None
        # (codec, extents, offset, args): a plain tuple before Pillow 11,
        # a named one since, so it is read by position.
        args = page.tile[0][3]
        rawmode = args[0] if isinstance(args, tuple) and args else args
        return rawmode if isinstance(rawmode, str) else None

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
        if self.wide_png:
            decoded = np.asarray(imagecodecs.png_decode(self.path.read_bytes()))
        else:
            with _unlimited_pixels(), PILImage.open(self.path) as page:
                frame = (
                    page if self.convert_to is None else page.convert(self.convert_to)
                )
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
        # Rounded outward, but a value a hair from a whole pixel is that
        # pixel: an affine fitted from teaching points gives 4.9999999 for 5.
        eps = 1e-6
        x0 = max(0, int(math.floor(x_min - pad_x + eps)))
        y0 = max(0, int(math.floor(y_min - pad_y + eps)))
        x1 = min(full_width, int(math.ceil(x_max + pad_x - eps)))
        y1 = min(full_height, int(math.ceil(y_max + pad_y - eps)))
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


@dataclass(frozen=True)
class OrientedOpticalSource:
    """A whole optical source, turned upright as its orientation tag says.

    ``shape`` is the upright one: a quarter turn (tags 5 to 8) swaps height
    and width. The inner source still decodes stored row bands, and each
    band is turned with :func:`orient` on its own, which gives the matching
    band of the upright image: rows for tags 2 to 4, columns for 5 to 8,
    counted from the far edge when the tag reads the stored rows bottom to
    top. :meth:`blocks` says where each one goes, so memory stays bounded
    by one band as before.
    """

    inner: OpticalSource
    tag: int

    @property
    def path(self) -> Path:
        """The file the inner source reads."""
        return self.inner.path

    @property
    def shape(self) -> Shape:
        """``(c, y, x)`` of the upright image."""
        c, y, x = self.inner.shape
        return (c, x, y) if self.tag in _TRANSPOSED else (c, y, x)

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
        """Decoded bytes of one stored row: a band is still stored rows."""
        return self.inner.row_bytes

    def blocks(self, rows: int) -> Iterator[Tuple[int, int, np.ndarray]]:
        """Yield ``(first_row, first_column, block)`` of the upright image."""
        stored_rows = self.inner.shape[1]
        for start, band in self.inner.bands(rows):
            stop = start + band.shape[1]
            first = stored_rows - stop if self.tag in _ROWS_REVERSED else start
            block = orient(band, self.tag)
            del band
            if self.tag in _TRANSPOSED:
                yield 0, first, block
            else:
                yield first, 0, block


#: What :class:`StreamedOpticalImage` streams from: a whole source, a
#: window of one, or one turned upright.
AnyOpticalSource = Union[OpticalSource, CroppedOpticalSource, OrientedOpticalSource]


def _source_blocks(
    source: AnyOpticalSource, rows: int
) -> Iterator[Tuple[int, int, np.ndarray]]:
    """``(first_row, first_column, block)`` of any source, in its own frame."""
    if isinstance(source, OrientedOpticalSource):
        yield from source.blocks(rows)
        return
    for start, band in source.bands(rows):
        yield start, 0, band


def _starts_like_tiff(path: Path) -> bool:
    """Whether the file's first bytes are a TIFF or BigTIFF signature."""
    try:
        with path.open("rb") as handle:
            return handle.read(4) in TIFF_SIGNATURES
    except OSError:
        return False


def probe_optical_source(path: Union[str, Path]) -> OpticalSource:
    """Read the header of an optical image, whatever raster format it is.

    TIFF goes to :class:`OpticalTiffSource`, which decodes row bands; every
    other format goes to :class:`OpticalRasterSource`, which decodes the
    page whole. A TIFF is told by its first bytes, not its suffix, so the
    same file gives the same store whatever it is named (an archive member
    copied out under its declared media type, a scan saved as ``.png``). A
    file with a TIFF suffix goes to tifffile even when it does not start
    like one, so the warning says what tifffile found wrong with it.

    Raises:
        Whatever the chosen probe raises for a file it cannot read. Callers
        skip that image with a warning; see the failure policy above.
    """
    path = Path(path)
    if _starts_like_tiff(path) or path.suffix.lower() in TIFF_SUFFIXES:
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
    #: What a warning calls the image: the file, or the archive member a
    #: reader copied it out of. None names the file the source reads.
    source_file: Optional[str] = None

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
        if isinstance(images, zarr.Group) and self.name in images:
            del images[self.name]

    def _open_level_arrays(self, store_path: Union[str, Path]) -> List[zarr.Array]:
        root = zarr.open_group(str(store_path), mode="r+", use_consolidated=False)
        element = root[f"images/{self.name}"]
        paths = _multiscale_dataset_paths(element.attrs)
        if not isinstance(element, zarr.Group) or paths is None:
            raise RuntimeError(
                f"images/{self.name} in {store_path} carries no multiscales metadata"
            )
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
        # Stored row bands fill target columns once the image is turned a
        # quarter, so the bands then tile the column unit instead.
        turned = (
            isinstance(self.source, OrientedOpticalSource)
            and self.source.tag in _TRANSPOSED
        )
        chunk_rows = _write_unit(target)[2 if turned else 1]
        # The budget is looked up here, not bound as a default, so a test
        # (or a caller) can narrow it.
        rows = band_rows(
            self.source.row_bytes,
            chunk_rows,
            self.source.strip_rows,
            budget=BAND_BUDGET_BYTES,
        )
        stored = (
            self.source.inner
            if isinstance(self.source, OrientedOpticalSource)
            else self.source
        )
        n_bands = -(-stored.shape[1] // rows)
        logger.debug(
            f"  level 0 of '{self.name}': {n_bands} band(s) of {rows} rows "
            f"({rows * self.source.row_bytes / 1e6:.0f} MB decoded each)"
        )
        for y0, x0, block in _source_blocks(self.source, rows):
            # A band of stored rows is a row band of the target, or a
            # column band when the image is turned a quarter.
            target[:, y0 : y0 + block.shape[1], x0 : x0 + block.shape[2]] = block
            # The next band is decoded before this name is rebound; drop
            # it now so only one band is ever alive.
            del block

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
                    # asarray only narrows zarr's array-or-scalar read type:
                    # the slice is already an ndarray, so nothing is copied.
                    piece = np.asarray(
                        source[
                            c_sel,
                            excess_y
                            + (y_sel.start + y0) * factor : excess_y
                            + (y_sel.start + y1) * factor,
                            excess_x
                            + (x_sel.start + x0) * factor : excess_x
                            + (x_sel.start + x1) * factor,
                        ]
                    )
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


#: Characters an optical element name keeps after ``<dataset_id>_optical_``,
#: and how many of them.
_NAME_CHARACTERS = frozenset("0123456789abcdefghijklmnopqrstuvwxyz_.")
NAME_BUDGET = 30


def _clean_stem(stem: str, budget: int = NAME_BUDGET) -> str:
    """A file stem made fit for an element name.

    Lower case; every character other than an ASCII letter, a digit, ``_``
    or ``.`` becomes ``_``, so the name is one spatialdata accepts and one
    every file system can hold. Dots at either end are stripped before and
    after the cut to ``budget`` characters: Windows drops a trailing dot
    from a folder name, which left the root attrs naming an element that
    was not there. A stem with nothing left is called ``image``.
    """
    cleaned = "".join(c if c in _NAME_CHARACTERS else "_" for c in stem.lower())
    return cleaned.strip(".")[:budget].strip(".") or "image"


def _numbered(suffix: str, number: int, budget: int = NAME_BUDGET) -> str:
    """``suffix`` with ``_<number>`` appended, inside the same budget."""
    tail = f"_{number}"
    return suffix[: max(0, budget - len(tail))].rstrip(".") + tail


def _reason(
    error: BaseException, path: Optional[Path] = None, display: str = ""
) -> str:
    """An exception as a warning gives it: its text, else its type.

    With ``path``, every spelling of it in the text (decoders quote the
    file they opened, some with doubled backslashes) becomes ``display``:
    for an image a reader copied out of an archive, the warning names the
    member, not a temporary copy.
    """
    text = str(error).strip() or type(error).__name__
    if path is not None and display:
        for spelling in (repr(str(path))[1:-1], str(path)):
            text = text.replace(spelling, display)
    return text


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
        output_path: Callable[[], Path],
        dataset_id: str,
        *,
        include: bool,
        apply_alignment: bool,
        pixel_size_xy: Callable[[], Tuple[float, float]],
        acquisition_origin: Optional[Callable[[], Optional[AcquisitionOrigin]]] = None,
    ) -> None:
        """Build the optical images of one conversion.

        Args:
            reader: The MSI reader. Only ``get_optical_image_paths``,
                ``get_primary_optical_image_path``, ``mis_metadata``,
                ``_positions`` and ``_header`` are read, and the last four
                through ``getattr``/``hasattr`` -- a reader here is whatever
                satisfies the interface, not necessarily a
                :class:`~thyra.core.base_reader.BaseMSIReader` subclass.
            output_path: The store being written, called when it is needed
                rather than read here: the converter's output path is not
                fixed at construction. ``convert_msi`` shortens it for
                Windows before building the converter, but a caller that
                stands one up itself has to do that afterwards, and the
                converter reads its own ``output_path`` at call time
                everywhere else. Capturing the path here streamed the
                pixels into whatever the store was called before the
                reassignment, which left the element the SpatialData write
                declared holding nothing but its fill value. Called only
                after the store exists, by :meth:`stream_pending_pixels`
                and :meth:`forget_image`.
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
            acquisition_origin: Where the source put grid index 0, called
                when an image is loaded, for the same reason. An image
                carried into the MSI's micrometres gets the
                ``"acquisition"`` system from it (D48); ``None`` writes
                none.
        """
        self.reader = reader
        self._output_path = output_path
        self.dataset_id = dataset_id
        self._include = bool(include)
        self._apply_alignment = bool(apply_alignment)
        self._pixel_size_xy = pixel_size_xy
        self._acquisition_origin = acquisition_origin
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
        # Which file stem holds each element name, so a second file that
        # cleans to a taken name gets a numbered one instead of replacing
        # the first; see :meth:`_unique_name`. Lower-cased keys, because
        # spatialdata refuses names that differ only in case.
        self._name_owners: Dict[str, str] = {}
        # The orientation tag of each element whose file carries one other
        # than 1, and whether the pixels were turned by it.
        self._orientations: Dict[str, Dict[str, Any]] = {}

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
        placed by -- not the Area's own corners, which are only where the
        ``.mis`` drew the outline. It is the extent of the region's cells,
        so the outermost pixels' full footprint is inside the box before
        the margin is, whatever the affine's rotation.
        """
        mapping = self._converted_region_mapping()
        if mapping is None or self._alignment is None or self._tic_to_image is None:
            return None
        matrix = self._tic_to_image
        # The region's cells in the TIC grid's own coordinates: cell i spans
        # [i, i + 1), so the block runs from the first index to last + 1.
        rx0 = mapping.raster_min_x - self._alignment.first_raster_x
        rx1 = mapping.raster_max_x - self._alignment.first_raster_x + 1
        ry0 = mapping.raster_min_y - self._alignment.first_raster_y
        ry1 = mapping.raster_max_y - self._alignment.first_raster_y + 1
        corners = (
            np.array(
                [[rx0, ry0, 1.0], [rx1, ry0, 1.0], [rx1, ry1, 1.0], [rx0, ry1, 1.0]]
            )
            @ np.asarray(matrix, dtype=np.float64).T
        )
        box = (
            float(corners[:, 0].min()),
            float(corners[:, 0].max()),
            float(corners[:, 1].min()),
            float(corners[:, 1].max()),
        )
        return OpticalCrop.around(box, full_width, full_height, mapping=mapping)

    def compute_alignment(self) -> None:
        """Compute optical-MSI alignment from reader metadata.

        For FlexImaging data with Area definitions, this computes the
        transformation that maps MSI raster coordinates to optical image
        pixel coordinates.
        """
        # A source that states its own registration (an mzPeak archive with
        # a teach_points affine) hands it over whole.
        stated: Optional[Callable[[], Optional[AreaAlignmentResult]]] = getattr(
            self.reader, "get_image_alignment", None
        )
        alignment = stated() if callable(stated) else None
        if alignment is not None:
            self._alignment = alignment
            logger.info(
                "Optical alignment stated by the source, with "
                f"{len(alignment.region_mappings)} region mappings"
            )
            return

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

        raster = mis_metadata.get("raster")
        reference = mis_metadata.get("reference_point")
        # The raster the source planned over its Areas, where it states
        # one (timsTOF GlobalMetadata, the rapifleX .dat header): it
        # places a run that stopped before filling its Area (D59).
        width, height = header.get("raster_width"), header.get("raster_height")
        planned = header.get("planned_grid") or (
            (first_raster_x, first_raster_y, width, height)
            if width and height
            else None
        )

        # Compute area-based alignment
        try:
            aligner = TeachingPointAlignment()
            self._alignment = aligner.compute_area_alignment(
                areas=areas,
                poslog_positions=positions,
                first_raster_x=first_raster_x,
                first_raster_y=first_raster_y,
                teaching_points=mis_metadata.get("teaching_points"),
                reference_point=tuple(reference) if reference else None,
                raster_step=tuple(raster) if raster else None,
                planned_grid=planned,
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
        """Build the affine from the TIC image's coordinates to image pixels.

        The matrix maps the TIC image's own coordinates, where cell
        ``(i, j)`` spans ``[i, i + 1) x [j, j + 1)`` as SpatialData draws
        it, to optical-image pixels. A cell's centre therefore lands on its
        spot, and the polygons, the crop and ``raster_to_global_affine`` use
        the same matrix.

        With the teaching-point lattice (D31) the matrix is the lattice's:
        exact, one for every region, rotation included. Without it, each
        region is stretched over its Area's bounding box (one region), or
        all regions over the box of all Areas, which is approximate.
        """
        if self._alignment is None:
            return
        if not self._alignment.region_mappings:
            return

        lattice = getattr(self._alignment, "lattice", None)
        if lattice is not None:
            self._tic_to_image = np.array(lattice.cell_to_image, dtype=np.float64)
            logger.info(
                "Built TIC-to-image affine from the teaching-point lattice: "
                f"{np.round(self._tic_to_image[:2], 4).tolist()}"
            )
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

        # tx/ty above place a cell's centre for an integer index; the TIC
        # image puts that centre at index + 0.5, so move back half a cell.
        tx -= half_x
        ty -= half_y

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
        when computing Scale transforms for the other images, or telling
        that a Scale cannot place one (D59).

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
                shown = self._display_name(image_path)
                logger.warning(
                    f"Failed to load optical image {shown}: "
                    f"{_reason(e, image_path, shown)}"
                )

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
        resolve: Optional[Callable[[], Optional[Path]]] = getattr(
            self.reader, "get_primary_optical_image_path", None
        )
        primary = resolve() if callable(resolve) else None
        if primary is not None:
            self.primary_filename = Path(primary).stem.lower()
            logger.info(f"Primary alignment image named by the reader: {primary.name}")

    def _is_primary_optical(self, image_path: Path) -> bool:
        """Check if an image file is the primary alignment image from .mis.

        Matched on the stem, not the whole filename: what the .mis names and
        what is on disk agree on the name but not always on the spelling of
        the extension, and a folder does not hold the same stem twice.
        """
        if not self.primary_filename:
            return False
        return image_path.stem.lower() == self.primary_filename

    def _scale_transform(self, x_size: int, y_size: int) -> Optional[Any]:
        """Compute a Scale transform for a non-primary optical image.

        Maps the image's pixel coordinates to the primary alignment image's
        coordinate space using the dimension ratio. That is right only for
        a resampled copy of the alignment image, such as FlexImaging's
        quarter-size ``_0001`` scan: the same ratio on both axes, give or
        take the rounding of one pixel. FlexImaging's overview photo is the
        source the alignment image was made from, turned by up to a degree
        or so; a turn widens the two sides by different ratios, and a Scale
        cannot express it (#497, D59).

        Args:
            x_size: Width of the non-primary image
            y_size: Height of the non-primary image

        Returns:
            Scale transform, or ``None`` when there is no alignment image
            to scale onto or the image is not a resampled copy of it.
        """
        if self._primary_dims is None:
            return None

        primary_w, primary_h = self._primary_dims
        scale_x = primary_w / x_size
        scale_y = primary_h / y_size
        # The same ratio on both axes, to within the rounding of the copy's
        # sides: half a copy pixel on each moves the cross product by up to
        # half the sum of the primary's sides.
        if abs(x_size * primary_h - y_size * primary_w) > max(primary_w, primary_h):
            logger.warning(
                f"  Not placed: {x_size}x{y_size} is not a resampled copy of the "
                f"{primary_w}x{primary_h} alignment image (scales {scale_x:.4f} "
                f"and {scale_y:.4f} differ). It is kept in its own coordinate "
                "system, not in 'global'."
            )
            return None

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

        The math: ``tic_to_image`` maps TIC coordinates, in which raster
        index ``i`` has its centre at ``i + 0.5``, to optical pixels. In
        micrometres a pixel's centre is ``obs["spatial_x"] = i *
        pixel_size_um``, where its polygon is centred. So the spot under
        an optical pixel goes back to its index and then to that centre::

            um = pixel_size_um * (inv(tic_to_image) @ optical_pixel - 0.5)

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
        to_index = np.array(
            [[1.0, 0.0, -0.5], [0.0, 1.0, -0.5], [0.0, 0.0, 1.0]], dtype=np.float64
        )
        matrix = scale_mat @ to_index @ inv
        return Affine(
            matrix,
            input_axes=("x", "y"),
            output_axes=("x", "y"),
        )

    def _transformations(self, transform: Any, placed_in_um: bool) -> Dict[str, Any]:
        """An image's transformations, ``"acquisition"`` included when placed.

        Only an image carried into the MSI's micrometres has a place in the
        source's frame. An optical-pixel store has no acquisition system at
        all (D31). An image with no alignment never comes here: it keeps a
        coordinate system of its own and nothing else (D69).
        """
        origin = (
            self._acquisition_origin()
            if placed_in_um and self._acquisition_origin is not None
            else None
        )
        return element_transformations(self.dataset_id, transform, origin)

    def _load_single_image(
        self, image_path: Path, images: MutableMapping[str, Any]
    ) -> None:
        """Load a single optical image and add it to data structures.

        An image is placed in ``"global"`` and the dataset's system only
        when the store holds an alignment for it: the primary image (the
        one the .mis <ImageFile> or the archive's teach points name) once
        the TIC-to-image affine exists, with an Identity transform (or its
        crop's translation, or the inverse affine into micrometres), and a
        resampled copy of it, with a Scale into the primary's pixel grid.
        Every other image keeps one coordinate system of its own, named
        after its element, with Identity and nothing else (D59, D69).

        The primary image and its resampled copies are stored as encoded,
        because the teaching points may refer to the encoded pixels; every
        other image is turned upright by its orientation tag (D69).

        Args:
            image_path: Path to the optical image (TIFF, JPEG, PNG or BMP)
            images: The store's images, to add the placeholder to
        """
        label = self._reader_label(image_path)
        source_file = label.source_file if label is not None else image_path.name
        owner = (label.name if label is not None else image_path.stem).lower()
        image_name = self._unique_name(self._element_suffix(image_path, label), owner)

        logger.info(f"Loading optical image: {source_file} as '{image_name}'")

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
        if isinstance(probed, OpticalTiffSource) and probed.dropped_planes:
            logger.warning(
                f"  {source_file}: only the first plane along "
                f"{probed.dropped_axes or 'its pages'} is kept; "
                f"{probed.dropped_planes} page(s) are left out"
            )

        is_primary = self._is_primary_optical(image_path)
        um_mode = not self._apply_alignment and self._tic_to_image is not None
        transform, placed, keep_encoded, placed_source = self._place(probed, is_primary)
        source = self._oriented(
            probed, placed_source, source_file, image_name, keep_encoded
        )
        n_channels, y_size, x_size = source.shape

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
        # An unplaced image gets a coordinate system of its own, named after
        # it: in "global" any transform would claim a placement.
        transformations = (
            self._transformations(transform, um_mode)
            if placed
            else {image_name: Identity()}
        )
        streamed = StreamedOpticalImage(
            source=source,
            name=image_name,
            chunks=image_chunks(2),  # (1, 4096, 4096); sharding seam, see _chunking
            scale_factors=scale_factors,
            transformations=transformations,
            source_file=source_file,
        )
        # Keyed by element name, as the images dict is: a second file with
        # the same stem (the same picture saved as .tif and .jpg) replaces
        # the first, the way the dict assignment always did, with a
        # warning. A different stem never lands here; see _unique_name.
        earlier = self._pending.get(image_name)
        if earlier is not None:
            logger.warning(
                f"Optical image '{image_name}' from "
                f"{earlier.source_file or earlier.source.path.name} "
                f"is replaced by {source_file}, which maps to the same name"
            )
            self._orientations.pop(image_name, None)
        images[image_name] = streamed.placeholder()
        self._pending[image_name] = streamed
        self._name_owners[image_name.lower()] = owner
        if probed.orientation not in (None, 1) and 1 <= int(probed.orientation) <= 8:
            self._orientations[image_name] = {
                "tag": int(probed.orientation),
                "applied": not keep_encoded,
            }
        # Which file this element came from, and whether it is the one the
        # alignment is stated against. Neither is recoverable from the
        # store otherwise: the element name drops the extension and
        # rewrites the stem (_0000 -> highres, and a stem is cleaned, cut
        # at 30 characters and numbered on a clash), and the alignment
        # image is only distinguishable by its transform, and only when
        # the alignment was applied. Recorded here, written by the root
        # attrs. For an image a reader copied out of a container, the
        # member it came from, not the copy.
        self._sources[image_name] = source_file
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

    def _place(
        self, probed: OpticalSource, is_primary: bool
    ) -> Tuple[Any, bool, bool, AnyOpticalSource]:
        """Where an image goes: ``(transform, placed, keep_encoded, source)``.

        Two cases, both only with an alignment:

        1. apply_optical_alignment=True (default): "global" is optical-image
           pixel space. The primary image is Identity; a resampled copy of
           it Scales to match the primary's dimensions.
        2. apply_optical_alignment=False (e.g. Ousia wizard): "global" is
           MSI micrometer space. The primary image's pixel coordinates map
           into MSI um through the inverse of the tic-to-image affine, then
           the pixel size, so the optical image lands alongside the MSI in
           the same um frame and downstream registration steps map both
           together.

        Without the affine (no .mis, an alignment that failed, an archive
        image with no teach points) "global" says nothing about where the
        image is, so it is not placed there at all (D69). ``keep_encoded``
        is true for the primary image and its resampled copies, which are
        not turned by their orientation tag. ``source`` is the crop of the
        primary image when the conversion covers one region.
        """
        aligned = self._tic_to_image is not None
        if is_primary:
            transform, source = self._place_primary(probed)
            placed, keep_encoded = aligned, True
        else:
            base = self._scale_transform(probed.shape[2], probed.shape[1])
            transform, source = self._place_copy(base), probed
            placed, keep_encoded = aligned and base is not None, base is not None
        if not aligned:
            logger.info(
                "  No alignment for this image in the store: it is kept in a "
                "coordinate system of its own, not in 'global'"
            )
        elif not placed and self._primary_dims is None:
            logger.info(
                "  The alignment image is not in the store, so nothing places "
                "this one: it is kept in a coordinate system of its own"
            )
        return transform, placed, keep_encoded, source

    def _place_primary(self, probed: OpticalSource) -> Tuple[Any, AnyOpticalSource]:
        """The primary image's transform, and its crop when it has one."""
        _, y_size, x_size = probed.shape
        # The FULL dimensions, crop or not: the other images scale into the
        # whole scan's pixel grid, which is the frame everything on the MSI
        # side is stated in.
        self._primary_dims = (x_size, y_size)
        transform: Any = Identity()
        if not self._apply_alignment and self._tic_to_image is not None:
            transform = self._to_um_transform()
            logger.info(
                f"  Primary image -> um via inverse alignment: {x_size}x{y_size}"
            )
        else:
            logger.info(f"  Primary alignment image: {x_size}x{y_size}")
        crop = self._crop_for_alignment_image(x_size, y_size)
        if crop is None:
            return transform, probed
        # The crop's pixel (u, v) is the scan's (x0 + u, y0 + v): that
        # translation goes FIRST, and then whatever carried the whole scan
        # into the store's frame carries the crop. So the raster-to-image
        # affine, the pixel polygons and the coordinate_systems attr keep
        # describing the scan's grid.
        shift = Translation([crop.x0, crop.y0], axes=("x", "y"))
        transform = (
            shift
            if isinstance(transform, Identity)
            else SequenceTransform([shift, transform])
        )
        self._crop = crop
        logger.info(
            f"  Cropped to region {crop.region_id} ('{crop.region_name}'): "
            f"{crop.width}x{crop.height} at ({crop.x0}, {crop.y0}) of "
            f"{crop.full_width}x{crop.full_height}"
        )
        return transform, CroppedOpticalSource(probed, crop)

    def _place_copy(self, base: Optional[Any]) -> Any:
        """A non-primary image's transform from its Scale onto the primary.

        First the Scale to match the primary, then, in micrometre mode, the
        same um affine. An image a Scale cannot place (``base`` is
        ``None``) stays out of both frames (D59); the stored dimensions
        decide, before any turn.
        """
        if base is None:
            return Identity()
        if not self._apply_alignment and self._tic_to_image is not None:
            return SequenceTransform([base, self._to_um_transform()])
        return base

    @staticmethod
    def _oriented(
        probed: OpticalSource,
        source: AnyOpticalSource,
        source_file: str,
        image_name: str,
        keep_encoded: bool,
    ) -> AnyOpticalSource:
        """``source`` turned upright by its orientation tag, where that is right.

        The alignment image and its resampled copies stay as encoded, with
        a warning, because the teaching points may be stated against the
        encoded pixels and no tagged FlexImaging file exists to check. A
        tag outside 1 to 8 is not applied either.
        """
        tag = probed.orientation
        if tag in (None, 1):
            return source
        if not 1 <= int(tag) <= 8:
            logger.warning(
                f"  {source_file}: orientation tag {tag} is not one of 1 to 8; "
                "the image is stored as encoded"
            )
            return source
        if keep_encoded:
            logger.warning(
                f"  {source_file} carries orientation tag {tag} but is the "
                "alignment image or a resampled copy of it, so it is stored as "
                f"encoded in '{image_name}': the teaching points may refer to "
                "the encoded pixels"
            )
            return source
        logger.info(f"  Turned upright by its orientation tag {tag}")
        return OrientedOpticalSource(probed, int(tag))

    def _unique_name(self, suffix: str, owner: str) -> str:
        """The element name for ``suffix``, numbered when another file holds it.

        A name already held by a file with the same stem is returned as it
        is, and that file is replaced (the same picture exported twice).
        A name held by a file with another stem -- two stems that differ
        only past the 30-character cut, or by characters the cleaning
        turns into ``_`` -- gets ``_1``, ``_2`` and so on, inside the same
        budget. Names are compared ignoring case, as spatialdata does.
        """
        number = 0
        candidate = suffix
        while True:
            name = f"{self.dataset_id}_optical_{candidate}"
            holder = self._name_owners.get(name.lower())
            if holder is None or holder == owner:
                return name
            number += 1
            candidate = _numbered(suffix, number)

    def _display_name(self, image_path: Path) -> str:
        """What a warning calls an image: the archive member, else the file."""
        try:
            label = self._reader_label(image_path)
        except Exception:  # a reader's naming must not hide the real failure
            label = None
        return label.source_file if label is not None else image_path.name

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
        store = self._output_path()
        streamed = 0
        for image in pending.values():
            logger.info(f"Streaming optical image pixels: '{image.name}'")
            try:
                image.stream_pixels(store)
            except Exception as e:  # mirrors the per-image guard in add_images
                shown = image.source_file or image.source.path.name
                logger.warning(
                    f"Failed to load optical image {shown}: "
                    f"{_reason(e, image.source.path, shown)}; "
                    f"dropping '{image.name}' from the store"
                )
                image.discard(store)
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
        self._orientations.pop(name, None)
        was_alignment = self._alignment_element == name
        if was_alignment:
            self._alignment_element = None
        try:
            root = zarr.open_group(
                str(self._output_path()), mode="r+", use_consolidated=False
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

    def _reader_label(self, image_path: Path) -> Optional[OpticalImageLabel]:
        """The reader's name for an image it copied out of its source.

        Asked for with ``getattr`` for the reason
        :meth:`_adopt_reader_primary_optical` gives: a reader here need not
        inherit from ``BaseMSIReader``.
        """
        resolve: Optional[Callable[[Path], Optional[OpticalImageLabel]]] = getattr(
            self.reader, "get_optical_image_label", None
        )
        if not callable(resolve):
            return None
        return resolve(image_path)

    def _element_suffix(
        self, image_path: Path, label: Optional[OpticalImageLabel] = None
    ) -> str:
        """The part of an optical element's name after ``<dataset_id>_optical_``.

        The file suffix is dropped, so the same acquisition exported as a
        .tif or a .jpg lands under the same element name. The ``_0000`` /
        ``_0001`` rule is for files named by the acquisition software; an
        image a reader copied out of a container is named by the reader's
        label. Anything else is the file stem, cleaned by
        :func:`_clean_stem`. :meth:`_unique_name` numbers it when another
        file already holds the name.

        Args:
            image_path: Path to the optical image
            label: The reader's name for the image, if it gave one.

        Returns:
            The suffix (e.g. 'highres', 'derived', 'slide_1_')
        """
        if label is not None:
            return _clean_stem(label.name)

        stem = image_path.stem.lower()

        # Extract meaningful suffix from filename
        if "_0000" in stem:
            suffix = "highres"
        elif "_0001" in stem:
            suffix = "derived"
        elif "deriv" in stem:
            suffix = "overview"
        else:
            suffix = _clean_stem(stem)

        return suffix

    def root_attr(self) -> Optional[Dict[str, Any]]:
        """What each optical element in this store came from, and which aligns.

        Two facts the store could not state before:

        * **Which file.** The element name is derived, not the filename:
          :meth:`_element_suffix` drops the extension, maps
          ``_0000``/``_0001``/``deriv`` onto ``highres``/``derived``/
          ``overview`` and cleans anything else (see :func:`_clean_stem`)
          and cuts it at 30 characters, so ``sample_0000.tif`` and
          ``sample_0000.jpg`` both land under
          ``<dataset_id>_optical_highres`` and neither name survives.
        * **Which one the alignment is stated against.** The .mis names it
          in ``<ImageFile>``; the store only ever implied it, through the
          transform (the alignment image gets ``Identity``, a resampled
          copy a ``Scale`` into its pixel grid) and only when
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
        * **How it was turned.** ``orientation`` records the file's
          orientation tag, where it is one other than 1, and whether the
          pixels were turned upright by it (``applied``). The alignment
          image and its resampled copies are stored as encoded (D69).

        Returns:
            ``{"alignment_element": str | None, "elements": {name:
            {"source_file": str, "crop"?: {...}, "orientation"?: {"tag":
            int, "applied": bool}}}}``, or ``None`` when the
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
        for name, orientation in self._orientations.items():
            if name in elements:
                elements[name]["orientation"] = dict(orientation)
        return {
            "alignment_element": self._alignment_element,
            "elements": elements,
        }
