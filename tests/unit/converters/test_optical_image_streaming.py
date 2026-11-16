"""The bounded optical-image path writes what the whole-page path wrote.

Three contracts, each pinned against the library code it replaces rather
than against a stored fixture:

* :func:`block_mean` / :func:`trim_excess` reproduce xarray's
  ``coarsen(boundary="trim", side="right").mean().astype(dtype)`` value for
  value, odd sizes and all dtypes an optical TIFF can carry included;
* :meth:`StreamedOpticalImage.placeholder` is structurally what
  ``Image2DModel.parse`` builds, with and without ``scale_factors``: same
  levels, shapes, chunks, dtype, per-level transformations and coordinates;
* declaring the placeholder through ``SpatialData.write`` and then streaming
  the pixels leaves a store whose optical element -- metadata, chunk grid,
  every pixel of every level -- equals the one the parse-and-write route
  produces, and the store never sees the page whole.

Plus the tolerance the old route had: a TIFF whose pixels cannot be read
is dropped with a warning, not fatal, on both converter routes. A TIFF cut
short is found from its header, so that holds whatever the decoder does
with the bytes that are missing.

And what the store says about the images it holds: which file each optical
element came from and which one the alignment is stated against, in the
root attrs -- the only place that survives, since the ome-zarr writer keeps
nothing a caller puts on the element itself.
"""

from __future__ import annotations

import json
import logging
import shutil
import struct
import threading
from pathlib import Path
from types import SimpleNamespace
from typing import List

import numpy as np
import pytest
import xarray as xr
from PIL import Image as PILImage

pytest.importorskip("spatialdata")
tifffile = pytest.importorskip("tifffile")
zarr = pytest.importorskip("zarr")

from spatialdata import SpatialData  # noqa: E402
from spatialdata.models import Image2DModel  # noqa: E402
from spatialdata.transformations import Scale, get_transformation  # noqa: E402

from thyra.converters.spatialdata import optical_image  # noqa: E402
from thyra.converters.spatialdata.optical_image import (  # noqa: E402
    OpticalTiffSource,
    StreamedOpticalImage,
    band_rows,
    block_mean,
    level_shapes,
    trim_excess,
)

CHUNKS = (1, 16, 16)
TRANSFORM = Scale([2.0, 0.5], axes=("x", "y"))
TRANSFORMATIONS = {"slide": TRANSFORM, "global": TRANSFORM}


def _xarray_reference(block: np.ndarray, factor: int) -> np.ndarray:
    """What multiscale_spatial_image computes for one pyramid step."""
    image = xr.DataArray(block, dims=("c", "y", "x"))
    return (
        image.coarsen(dim={"y": factor, "x": factor}, boundary="trim", side="right")
        .mean()
        .astype(block.dtype)
        .values
    )


def _stream_reference(block: np.ndarray, factor: int) -> np.ndarray:
    """The same step the way stream_pixels computes it: trim, then block mean."""
    ey = trim_excess(block.shape[1], factor)
    ex = trim_excess(block.shape[2], factor)
    return block_mean(block[:, ey:, ex:], factor)


def _parse_reference(whole: np.ndarray, scale_factors):
    """The route this module replaces: whole page -> Image2DModel.parse."""
    kwargs = {"transformations": dict(TRANSFORMATIONS), "chunks": CHUNKS}
    if scale_factors:
        kwargs["scale_factors"] = scale_factors
    return Image2DModel.parse(
        xr.DataArray(
            whole,
            dims=("c", "y", "x"),
            coords={
                "c": np.arange(whole.shape[0]),
                "y": np.arange(whole.shape[1]),
                "x": np.arange(whole.shape[2]),
            },
        ),
        **kwargs,
    )


def _streamed(source: OpticalTiffSource, scale_factors, name="optical"):
    return StreamedOpticalImage(
        source=source,
        name=name,
        chunks=CHUNKS,
        scale_factors=scale_factors,
        transformations=dict(TRANSFORMATIONS),
    )


def _cyx(path: Path) -> np.ndarray:
    decoded = tifffile.imread(str(path))
    return decoded[np.newaxis] if decoded.ndim == 2 else np.moveaxis(decoded, -1, 0)


@pytest.mark.parametrize(
    "dtype",
    [np.uint8, np.uint16, np.int16, np.uint32, np.int32, np.float32, np.float64],
)
@pytest.mark.parametrize("shape", [(3, 8, 8), (1, 9, 7), (2, 15, 33), (3, 5, 6)])
@pytest.mark.parametrize("factor", [2, 3])
def test_block_mean_matches_xarray_coarsen(dtype, shape, factor):
    rng = np.random.default_rng(int(np.prod(shape)) * factor)
    if np.issubdtype(dtype, np.integer):
        info = np.iinfo(dtype)
        block = rng.integers(info.min, info.max, size=shape, dtype=dtype)
    else:
        block = rng.random(size=shape).astype(dtype) * 300 - 100
    expected = _xarray_reference(block, factor)
    actual = _stream_reference(block, factor)
    assert actual.dtype == expected.dtype
    np.testing.assert_array_equal(actual, expected)


def test_block_mean_extremes_and_nans():
    # Saturated unsigned windows: the exact accumulator must not overflow.
    full = np.full((2, 6, 4), 255, dtype=np.uint8)
    np.testing.assert_array_equal(
        _stream_reference(full, 2), _xarray_reference(full, 2)
    )
    full32 = np.full((1, 4, 6), np.iinfo(np.uint32).max, dtype=np.uint32)
    np.testing.assert_array_equal(
        _stream_reference(full32, 3), _xarray_reference(full32, 3)
    )
    # NaNs are skipped the way xarray's nanmean skips them.
    floats = np.arange(2 * 4 * 4, dtype=np.float32).reshape(2, 4, 4)
    floats[0, 0, 0] = np.nan
    floats[1, 2:4, 2:4] = np.nan
    expected = _xarray_reference(floats, 2)
    actual = _stream_reference(floats, 2)
    np.testing.assert_array_equal(np.isnan(actual), np.isnan(expected))
    np.testing.assert_array_equal(
        actual[~np.isnan(actual)], expected[~np.isnan(expected)]
    )


def test_trim_drops_the_leading_rows():
    """side="right" keeps the LAST full windows; the excess comes off the front."""
    block = np.arange(3 * 5 * 4, dtype=np.uint8).reshape(3, 5, 4)
    reference = _xarray_reference(block, 2)
    assert trim_excess(5, 2) == 1
    assert trim_excess(4, 2) == 0
    # Row 0 is dropped: the first output row averages rows 1 and 2.
    np.testing.assert_array_equal(
        reference[0, 0],
        block[0, 1:3, :].reshape(2, 2, 2).mean(axis=(0, 2)).astype(np.uint8),
    )
    np.testing.assert_array_equal(_stream_reference(block, 2), reference)


def test_block_mean_rejects_untrimmed_input():
    with pytest.raises(ValueError, match="not a multiple"):
        block_mean(np.zeros((1, 5, 4), dtype=np.uint8), 2)
    with pytest.raises(ValueError, match="\\(c, y, x\\)"):
        block_mean(np.zeros((4, 4), dtype=np.uint8), 2)


def test_level_shapes_floor_like_coarsen():
    assert level_shapes((3, 11748, 36736), [2, 2, 2]) == [
        (3, 11748, 36736),
        (3, 5874, 18368),
        (3, 2937, 9184),
        (3, 1468, 4592),
    ]
    assert level_shapes((1, 10, 10), []) == [(1, 10, 10)]


def test_band_rows_are_whole_strips_that_tile_the_chunk_rows():
    row = 36736 * 3  # the 110 KB rows of the real brightfield
    # A whole chunk row fits the budget: one band per chunk row.
    assert band_rows(1000, 4096, budget=10_000_000) == 4096
    assert band_rows(row, 4096, strip_rows=1, budget=512 << 20) == 4096
    # Half a chunk row fits: the largest divisor, still a whole number of strips.
    assert band_rows(row, 4096, strip_rows=1, budget=256 << 20) == 2048
    assert band_rows(row, 4096, strip_rows=8, budget=256 << 20) == 2048
    # Strips that do not divide the chunk: whole strips within the budget.
    assert band_rows(row, 4096, strip_rows=1000, budget=256 << 20) == 2000
    # A single-strip page (rowsperstrip == height): one band, the whole page.
    assert band_rows(row, 4096, strip_rows=11748, budget=256 << 20) == 11748
    # Never less than one strip, even when a strip is over budget.
    assert band_rows(1 << 30, 4096, strip_rows=64, budget=1) == 64
    assert band_rows(1 << 30, 4096, strip_rows=1, budget=1) == 1
    with pytest.raises(ValueError):
        band_rows(0, 4096)


@pytest.fixture
def rgb_tiff(tmp_path: Path) -> Path:
    """An odd-sized RGB TIFF with strips of one row, like a FlexImaging export."""
    rng = np.random.default_rng(7)
    image = rng.integers(0, 256, size=(37, 53, 3), dtype=np.uint8)
    # A fully black band, so an all-fill chunk exercises write_empty_chunks.
    image[10:20, :, :] = 0
    path = tmp_path / "brightfield_0000.tiff"
    tifffile.imwrite(str(path), image, rowsperstrip=1, compression="lzw")
    return path


def test_probe_reads_only_the_header(rgb_tiff: Path):
    source = OpticalTiffSource.probe(rgb_tiff)
    assert source.shape == (3, 37, 53)
    assert source.dtype == np.uint8
    assert source.page_axes == "YXS"
    assert source.strip_rows == 1
    assert source.row_bytes == 53 * 3


def test_bands_decode_row_ranges(rgb_tiff: Path):
    source = OpticalTiffSource.probe(rgb_tiff)
    whole = _cyx(rgb_tiff)
    starts = []
    for start, band in source.bands(rows=16):
        starts.append(start)
        assert band.shape[0] == 3
        np.testing.assert_array_equal(band, whole[:, start : start + band.shape[1], :])
    assert starts == [0, 16, 32]


def test_probe_reports_strip_and_tile_heights(tmp_path: Path):
    rng = np.random.default_rng(1)
    image = rng.integers(0, 256, size=(70, 90, 3), dtype=np.uint8)
    default = tmp_path / "default.tif"
    tifffile.imwrite(str(default), image)  # tifffile's default: one strip
    assert OpticalTiffSource.probe(default).strip_rows == 70
    strips = tmp_path / "strips.tif"
    tifffile.imwrite(str(strips), image, rowsperstrip=8)
    assert OpticalTiffSource.probe(strips).strip_rows == 8
    tiles = tmp_path / "tiles.tif"
    tifffile.imwrite(str(tiles), image, tile=(16, 16))
    assert OpticalTiffSource.probe(tiles).strip_rows == 16


def test_single_strip_page_is_decoded_once_whole(tmp_path: Path, monkeypatch):
    """A page stored as one strip cannot be read in pieces: one band, one decode."""
    rng = np.random.default_rng(2)
    image = rng.integers(0, 256, size=(64, 40, 3), dtype=np.uint8)
    path = tmp_path / "single_strip.tif"
    tifffile.imwrite(str(path), image)
    source = OpticalTiffSource.probe(path)
    assert source.strip_rows == 64
    rows = band_rows(source.row_bytes, 16, source.strip_rows, budget=40 * 3 * 4)
    assert rows == 64
    calls = []
    original = tifffile.TiffPage.asarray

    def counting_asarray(self, *args, **kwargs):
        calls.append(1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(tifffile.TiffPage, "asarray", counting_asarray)
    bands = list(source.bands(rows))
    assert [start for start, _ in bands] == [0]
    np.testing.assert_array_equal(bands[0][1], np.moveaxis(image, -1, 0))
    assert len(calls) == 1


def test_probe_rejects_what_the_decoder_would_have(tmp_path: Path, monkeypatch):
    rng = np.random.default_rng(3)
    volume = rng.integers(0, 256, size=(4, 16, 16), dtype=np.uint8)
    path = tmp_path / "volume.tif"
    # One ZYX page (ImageDepth), the kind of 3-d page that is not an image.
    tifffile.imwrite(str(path), volume, photometric="minisblack", volumetric=True)
    with pytest.raises(ValueError, match="layout"):
        OpticalTiffSource.probe(path)

    # tifffile reports a (SampleFormat, BitsPerSample) pair it cannot map as
    # dtype None; np.dtype(None) would silently be float64.
    class _UntypedPage:
        shape = (16, 16, 3)
        dtype = None
        axes = "YXS"
        chunks = (1, 16, 3)

    class _UntypedTiff:
        pages = [_UntypedPage()]

        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_exc):
            return False

    monkeypatch.setattr(optical_image.tifffile, "TiffFile", _UntypedTiff)
    with pytest.raises(ValueError, match="sample format"):
        OpticalTiffSource.probe(path)


def test_grayscale_and_planar_layouts(tmp_path: Path):
    rng = np.random.default_rng(3)
    gray = rng.integers(0, 65535, size=(21, 30), dtype=np.uint16)
    gray_path = tmp_path / "gray.tif"
    tifffile.imwrite(str(gray_path), gray, rowsperstrip=4)
    source = OpticalTiffSource.probe(gray_path)
    assert source.shape == (1, 21, 30)
    assert source.dtype == np.uint16
    bands = list(source.bands(rows=8))
    assert [start for start, _ in bands] == [0, 8, 16]
    np.testing.assert_array_equal(
        np.concatenate([b for _, b in bands], axis=1), gray[np.newaxis]
    )

    # One plane per sample: already (c, y, x), and streamed by rows too.
    planar = rng.integers(0, 256, size=(3, 40, 13), dtype=np.uint8)
    planar_path = tmp_path / "planar.tif"
    tifffile.imwrite(
        str(planar_path),
        planar,
        photometric="rgb",
        planarconfig="separate",
        rowsperstrip=8,
    )
    source = OpticalTiffSource.probe(planar_path)
    assert source.page_axes == "SYX"
    assert source.shape == (3, 40, 13)
    assert source.strip_rows == 8
    bands = list(source.bands(rows=16))
    assert [start for start, _ in bands] == [0, 16, 32]
    np.testing.assert_array_equal(np.concatenate([b for _, b in bands], axis=1), planar)


@pytest.mark.parametrize("scale_factors", [[2, 2], []])
def test_placeholder_mirrors_image2dmodel_parse(rgb_tiff: Path, scale_factors):
    placeholder = _streamed(
        OpticalTiffSource.probe(rgb_tiff), scale_factors
    ).placeholder()
    reference = _parse_reference(_cyx(rgb_tiff), scale_factors)
    assert type(placeholder) is type(reference)

    if scale_factors:
        assert list(placeholder.children) == list(reference.children)
        pairs = [
            (placeholder[level]["image"], reference[level]["image"])
            for level in reference.children
        ]
    else:
        pairs = [(placeholder, reference)]
    for ours, theirs in pairs:
        assert ours.dims == theirs.dims
        assert ours.shape == theirs.shape
        assert ours.dtype == theirs.dtype
        assert ours.data.chunks == theirs.data.chunks
        for axis in ("c", "y", "x"):
            np.testing.assert_array_equal(
                ours.coords[axis].values, theirs.coords[axis].values
            )
        for system in TRANSFORMATIONS:
            assert get_transformation(ours, system) == get_transformation(
                theirs, system
            )


def _read_element(store: Path, name: str):
    group = zarr.open_group(str(store), mode="r", use_consolidated=False)[
        f"images/{name}"
    ]
    attrs = group.attrs.asdict()
    datasets = attrs["ome"]["multiscales"][0]["datasets"]
    levels = {}
    for dataset in datasets:
        array = group[dataset["path"]]
        levels[dataset["path"]] = (array.metadata.to_dict(), array[:])
    files = sorted(
        str(f.relative_to(store / "images" / name))
        for f in (store / "images" / name).rglob("*")
        if f.is_file()
    )
    return attrs, levels, files


@pytest.mark.parametrize("scale_factors", [[2, 2], []])
def test_declared_then_streamed_store_equals_parse_and_write(
    rgb_tiff: Path, tmp_path: Path, scale_factors
):
    name = "optical"
    whole = _cyx(rgb_tiff)
    expected_store = tmp_path / "expected.zarr"
    SpatialData(images={name: _parse_reference(whole, scale_factors)}).write(
        expected_store
    )

    streamed = _streamed(OpticalTiffSource.probe(rgb_tiff), scale_factors, name)
    actual_store = tmp_path / "actual.zarr"
    SpatialData(images={name: streamed.placeholder()}).write(actual_store)
    declared_files = [
        f for f in (actual_store / "images" / name).rglob("*") if f.is_file()
    ]
    # Declaring stores metadata only: one zarr.json per level plus the group's.
    assert all(f.name == "zarr.json" for f in declared_files)
    assert len(declared_files) == len(scale_factors) + 2

    streamed.stream_pixels(actual_store)

    expected_attrs, expected_levels, expected_files = _read_element(
        expected_store, name
    )
    actual_attrs, actual_levels, actual_files = _read_element(actual_store, name)
    assert json.dumps(actual_attrs, sort_keys=True) == json.dumps(
        expected_attrs, sort_keys=True
    )
    assert list(actual_levels) == list(expected_levels)
    assert list(actual_levels) == [f"s{i}" for i in range(len(scale_factors) + 1)]
    for path in expected_levels:
        expected_meta, expected_pixels = expected_levels[path]
        actual_meta, actual_pixels = actual_levels[path]
        assert json.dumps(actual_meta, sort_keys=True, default=str) == json.dumps(
            expected_meta, sort_keys=True, default=str
        )
        np.testing.assert_array_equal(actual_pixels, expected_pixels)
    # Same chunk files, so all-fill chunks are skipped the same way too.
    assert actual_files == expected_files

    # And the pixels are the TIFF's.
    np.testing.assert_array_equal(actual_levels["s0"][1], whole)
    # Reading back through spatialdata still sees a valid image.
    sdata = SpatialData.read(str(actual_store))
    element = sdata.images[name]
    if scale_factors:
        assert list(element.children) == ["scale0", "scale1", "scale2"]
    else:
        assert isinstance(element, xr.DataArray)
    axes = ("x", "y")
    np.testing.assert_array_equal(
        get_transformation(element, "slide").to_affine_matrix(axes, axes),
        TRANSFORM.to_affine_matrix(axes, axes),
    )


def test_streaming_a_small_band_budget_still_matches(
    rgb_tiff: Path, tmp_path: Path, monkeypatch
):
    """Bands narrower than a chunk row (partial-chunk writes) change nothing."""
    monkeypatch.setattr(optical_image, "BAND_BUDGET_BYTES", 53 * 3 * 3)  # three rows
    monkeypatch.setattr(optical_image, "REDUCE_STRIP_ROWS", 2)
    source = OpticalTiffSource.probe(rgb_tiff)
    assert band_rows(source.row_bytes, 16, source.strip_rows, budget=53 * 3 * 3) == 2
    streamed = _streamed(source, [2, 2])
    store = tmp_path / "narrow.zarr"
    SpatialData(images={"optical": streamed.placeholder()}).write(store)
    bands_seen = []
    original = OpticalTiffSource.bands

    def counting_bands(self, rows):
        bands_seen.append(rows)
        return original(self, rows)

    monkeypatch.setattr(OpticalTiffSource, "bands", counting_bands)
    streamed.stream_pixels(store)
    assert bands_seen == [2]

    whole = _cyx(rgb_tiff)
    _, levels, _ = _read_element(store, "optical")
    np.testing.assert_array_equal(levels["s0"][1], whole)
    np.testing.assert_array_equal(levels["s1"][1], _stream_reference(whole, 2))
    np.testing.assert_array_equal(
        levels["s2"][1], _stream_reference(_stream_reference(whole, 2), 2)
    )


def test_the_reducer_never_reads_through_the_tiff_decoder(
    rgb_tiff: Path, tmp_path: Path, monkeypatch
):
    """Issue #284 item 3, measured, and it does not hold.

    The issue proposed that the reducer's worker threads read concurrently
    from one ``zarr.Array`` opened over tifffile's zarr adapter, backed by a
    single ``TiffFile`` handle, and that this needed a lock or a handle per
    worker if the adapter turned out not to be safe for concurrent reads.

    No worker ever touches that adapter. ``page.aszarr()`` has exactly one
    caller, :meth:`OpticalTiffSource.bands`, which level 0 drives from the
    calling thread inside its own ``with tifffile.TiffFile(...)``. What
    :meth:`_reduce_level` hands the pool is ``arrays[k-1]`` from
    :meth:`_open_level_arrays` -- the pyramid level already written to the
    output store, a plain zarr array over a ``LocalStore``.

    So this pins three things: the TIFF decoder is only ever entered from
    the calling thread, the pool nonetheless really does run wide (or the
    first assertion would be vacuous), and what it reads is the store on
    disk.
    """
    decoder_threads: List[threading.Thread] = []
    worker_threads: List[threading.Thread] = []
    reduced_from: List[object] = []

    real_bands = OpticalTiffSource.bands
    real_block_mean = optical_image.block_mean
    real_reduce = StreamedOpticalImage._reduce_level

    def spy_bands(self, rows):
        decoder_threads.append(threading.current_thread())
        yield from real_bands(self, rows)

    def spy_block_mean(*args, **kwargs):
        # Only _reduce_level's workers call this.
        worker_threads.append(threading.current_thread())
        return real_block_mean(*args, **kwargs)

    def spy_reduce(self, source, target, factor):
        reduced_from.append(source.store)
        return real_reduce(self, source, target, factor)

    monkeypatch.setattr(OpticalTiffSource, "bands", spy_bands)
    monkeypatch.setattr(optical_image, "block_mean", spy_block_mean)
    monkeypatch.setattr(StreamedOpticalImage, "_reduce_level", spy_reduce)

    streamed = _streamed(OpticalTiffSource.probe(rgb_tiff), [2, 2])
    store = tmp_path / "threads.zarr"
    SpatialData(images={"optical": streamed.placeholder()}).write(store)
    streamed.stream_pixels(store)

    assert decoder_threads == [threading.current_thread()]
    assert len({thread.name for thread in worker_threads}) > 1
    assert reduced_from and all(
        isinstance(store_, zarr.storage.LocalStore) for store_ in reduced_from
    )


def test_stream_pixels_refuses_a_store_without_the_element(
    rgb_tiff: Path, tmp_path: Path
):
    source = OpticalTiffSource.probe(rgb_tiff)
    streamed = _streamed(source, [2])
    store = tmp_path / "wrong.zarr"
    # A different pyramid was declared under that name.
    SpatialData(images={"optical": _streamed(source, [2, 2]).placeholder()}).write(
        store
    )
    with pytest.raises(RuntimeError, match="pyramid levels"):
        streamed.stream_pixels(store)
    shutil.rmtree(store)
    with pytest.raises((KeyError, FileNotFoundError)):
        streamed.stream_pixels(store)


def test_discard_removes_the_declared_element(rgb_tiff: Path, tmp_path: Path):
    streamed = _streamed(OpticalTiffSource.probe(rgb_tiff), [2])
    store = tmp_path / "declared.zarr"
    SpatialData(images={"optical": streamed.placeholder()}).write(store)
    assert (store / "images" / "optical").exists()
    streamed.discard(store)
    assert not (store / "images" / "optical").exists()
    streamed.discard(store)  # idempotent
    zarr.consolidate_metadata(str(store))
    assert list(SpatialData.read(str(store)).images) == []


# ---- the converter routes ----------------------------------------------


def _convert(reader, output_path: Path):
    from thyra.converters.spatialdata.streaming_converter import (
        StreamingSpatialDataConverter,
    )

    converter = StreamingSpatialDataConverter(
        reader, output_path, dataset_id="ds", pixel_size_um=10.0, use_csc=True
    )
    return converter, converter.convert()


def _mock_reader(*optical_paths: Path):
    from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader

    return MockMSIReader(
        MockMSIConfig(n_x=4, n_y=3, n_mz_bins=32, peaks_per_spectrum=(2, 4)),
        optical_image_paths=list(optical_paths),
    )


def test_converters_stream_the_pixels_after_the_write(rgb_tiff: Path, tmp_path: Path):
    """Both write paths declare the placeholder, then fill it: the store has the pixels."""
    output_path = tmp_path / "rgb.zarr"
    converter, success = _convert(_mock_reader(rgb_tiff), output_path)
    assert success is True
    # Nothing is left waiting: every declared image was streamed.
    assert converter.optical.pending == {}

    whole = _cyx(rgb_tiff)
    _, levels, _ = _read_element(output_path, "ds_optical_highres")
    # 37 x 53 is far below the 1000 px coarsest-level floor: a single level.
    assert list(levels) == ["s0"]
    np.testing.assert_array_equal(levels["s0"][1], whole)
    sdata = SpatialData.read(str(output_path))
    np.testing.assert_array_equal(sdata.images["ds_optical_highres"].values, whole)


def test_the_pixels_follow_an_output_path_reassigned_after_construction(
    rgb_tiff: Path, tmp_path: Path
):
    """A converter whose output path moved after it was built still fills its store.

    ``convert_msi`` shortens the output path for Windows
    (``prepare_zarr_output_path``) *before* it builds the converter, so the
    production route never moves it. A caller that stands a converter up
    itself has to shorten afterwards: ``test_internal_failures_are_not_warnings``
    does it three times, and the CLI's ``_quarantine_partial_output`` records
    that ``convert_msi`` "resolves and may extend the output path".

    ``OpticalImages`` captured the path at construction, so it streamed into
    whatever the store had been called before the move, while the converter
    wrote, scratched and consolidated at a path it read at call time. What
    came out was a store whose optical element had its metadata, its chunk
    grid and its name in the root attrs, and not one pixel -- every chunk
    still the fill value -- and a conversion reported as failed, because
    dropping the element again went to the same absent path and that failure
    is the one ``stream_pending_pixels`` deliberately propagates.
    """
    from thyra.converters.spatialdata.streaming_converter import (
        StreamingSpatialDataConverter,
    )

    declared = tmp_path / "declared.zarr"
    written = tmp_path / "written.zarr"

    converter = StreamingSpatialDataConverter(
        _mock_reader(rgb_tiff),
        declared,
        dataset_id="ds",
        pixel_size_um=10.0,
        use_csc=True,
    )
    converter.output_path = written

    assert converter.convert() is True
    assert converter.optical.pending == {}
    assert not declared.exists()

    whole = _cyx(rgb_tiff)
    _, levels, _ = _read_element(written, "ds_optical_highres")
    np.testing.assert_array_equal(levels["s0"][1], whole)
    sdata = SpatialData.read(str(written))
    np.testing.assert_array_equal(sdata.images["ds_optical_highres"].values, whole)


@pytest.fixture
def truncated_tiff(tmp_path: Path) -> Path:
    """Header intact, strips cut: probe succeeds, decoding fails."""
    rng = np.random.default_rng(11)
    image = rng.integers(0, 256, size=(40, 30, 3), dtype=np.uint8)
    path = tmp_path / "truncated_0000.tiff"
    tifffile.imwrite(str(path), image, rowsperstrip=1, compression="lzw")
    data = path.read_bytes()
    path.write_bytes(data[: len(data) // 2])
    source = OpticalTiffSource.probe(path)
    assert source.shape == (3, 40, 30)
    return path


def test_unreadable_pixels_drop_the_image_not_the_conversion(
    truncated_tiff: Path, tmp_path: Path, thyra_logs
):
    """A TIFF that cannot be decoded is skipped with a warning, as it always was."""
    output_path = tmp_path / "truncated.zarr"
    with thyra_logs("thyra.converters", logging.WARNING) as records:
        converter, success = _convert(_mock_reader(truncated_tiff), output_path)
    assert success is True
    assert converter.optical.pending == {}
    assert any(
        "Failed to load optical image truncated_0000.tiff" in record.message
        for record in records
    )
    sdata = SpatialData.read(str(output_path))
    assert [k for k in sdata.images if "optical" in k] == []
    assert not (output_path / "images" / "ds_optical_highres").exists()
    # The rest of the store is intact and consolidated.
    assert "ds_z0_tic" in sdata.images
    assert (output_path / "zarr.json").exists()


# ---- a TIFF cut short ---------------------------------------------------
#
# tifffile raised on a strip or tile the file does not hold until 2026.11.12,
# which fills it with zeros and returns the page. These tests hold on both:
# the cut is found from the header, before anything is decoded.

#: ``tifffile.imwrite`` arguments for a 64 x 48 RGB page, the segment the
#: file is cut at, and the first row that segment holds.
_CUT_LAYOUTS = {
    "lzw strips": ({"rowsperstrip": 8, "compression": "lzw"}, 5, 40),
    "raw strips": ({"rowsperstrip": 8}, 5, 40),
    "deflate tiles": ({"tile": (16, 16), "compression": "zlib"}, 6, 32),
    "raw tiles": ({"tile": (16, 16)}, 6, 32),
    # One plane per sample: the cut is in the last plane, so the first two
    # are whole and the read still has to be refused.
    "planar strips": (
        {"rowsperstrip": 8, "planarconfig": "separate", "photometric": "rgb"},
        21,
        40,
    ),
}


def _cut_at(path: Path, offset: int, kept: int = 0) -> None:
    """End the file ``kept`` bytes after ``offset``."""
    path.write_bytes(path.read_bytes()[: offset + kept])


def _rows(source: OpticalTiffSource, rows: int, window=None) -> np.ndarray:
    return np.concatenate([band for _, band in source.bands(rows, window)], axis=1)


@pytest.mark.parametrize("kept", [0, 3])
@pytest.mark.parametrize("layout", sorted(_CUT_LAYOUTS))
def test_a_cut_tiff_is_refused_where_it_is_cut_and_read_where_it_is_not(
    tmp_path: Path, layout: str, kept: int
):
    """No read returns pixels the file does not hold; a read above the cut is whole.

    ``kept`` bytes of the first lost segment are left in the file: none, or
    a few, which is the case a decoder can pad instead of refusing.
    """
    arguments, segment, lost_row = _CUT_LAYOUTS[layout]
    rng = np.random.default_rng(13)
    image = rng.integers(1, 256, size=(64, 48, 3), dtype=np.uint8)
    whole = np.moveaxis(image, -1, 0)
    path = tmp_path / "cut.tif"
    planar = "planarconfig" in arguments
    tifffile.imwrite(str(path), whole if planar else image, **arguments)
    with tifffile.TiffFile(path) as tif:
        offsets = [int(n) for n in tif.pages[0].dataoffsets]
    assert offsets == sorted(offsets)
    _cut_at(path, offsets[segment], kept)

    source = OpticalTiffSource.probe(path)
    assert source.shape == (3, 64, 48)
    # Decoded whole, and in bands through the zarr adapter.
    for rows in (64, 8):
        with pytest.raises(ValueError, match="the file is cut short"):
            list(source.bands(rows))
    # Everything above the cut is the picture.
    np.testing.assert_array_equal(
        _rows(source, 8, (0, lost_row, 3, 40)), whole[:, :lost_row, 3:40]
    )
    with pytest.raises(ValueError, match="the file is cut short"):
        list(source.bands(8, (0, lost_row + 1, 3, 40)))


def test_a_cut_tile_is_only_refused_by_a_read_that_needs_it(tmp_path: Path):
    """Tiles are tested by column too: the adapter never reads the others."""
    rng = np.random.default_rng(14)
    image = rng.integers(1, 256, size=(64, 48, 3), dtype=np.uint8)
    whole = np.moveaxis(image, -1, 0)
    path = tmp_path / "cut_tile.tif"
    tifffile.imwrite(str(path), image, tile=(16, 16))
    with tifffile.TiffFile(path) as tif:
        offsets = [int(n) for n in tif.pages[0].dataoffsets]
    # Four rows of three tiles. The last tile is rows 48-63, columns 32-47.
    _cut_at(path, offsets[11], 5)
    source = OpticalTiffSource.probe(path)
    np.testing.assert_array_equal(_rows(source, 16, (0, 64, 0, 32)), whole[:, :, :32])
    np.testing.assert_array_equal(_rows(source, 16, (0, 48, 0, 48)), whole[:, :48])
    with pytest.raises(ValueError, match="tile 12 of 12"):
        list(source.bands(16, (0, 64, 0, 33)))


def _fake_page(offsets, counts, **layout) -> SimpleNamespace:
    """A page as ``_refuse_cut_short`` reads one: two tables and a layout.

    Four rows of four uint8 samples, one row a strip, unless ``layout``
    says otherwise: every segment needs four bytes.
    """
    fields = {
        "keyframe": None,
        "dataoffsets": tuple(offsets),
        "databytecounts": tuple(counts),
        "imagelength": 4,
        "imagewidth": 4,
        "is_tiled": False,
        "rowsperstrip": 1,
        "tilelength": 0,
        "tilewidth": 0,
        "samplesperpixel": 1,
        "planarconfig": 1,
        "dtype": np.dtype(np.uint8),
        "bitspersample": 8,
        "compression": 1,
        "fillorder": 1,
    }
    fields.update(layout)
    return SimpleNamespace(**fields)


def test_only_pixels_the_file_does_not_hold_count_as_cut():
    refuse = optical_image._refuse_cut_short
    name = Path("slide.tif")
    everything = (0, 4, 0, 4)

    # A segment ending on the last byte is whole; one byte fewer is not.
    page = _fake_page((100, 104, 108, 112), (4, 4, 4, 4))
    refuse(name, page, everything, 116)
    with pytest.raises(ValueError) as refused:
        refuse(name, page, everything, 115)
    assert str(refused.value) == (
        "slide.tif: the file is cut short: strip 4 of 4 ends at byte 116 "
        "and the file has 115 bytes"
    )
    # Rows the read does not touch are not tested.
    refuse(name, page, (0, 3, 0, 4), 115)

    # A byte count that overshoots the file is not a cut while the bytes
    # that are there hold the segment's samples.
    generous = _fake_page((100, 104, 108, 112), (4, 4, 4, 40))
    refuse(name, generous, everything, 116)
    with pytest.raises(ValueError, match="strip 4 of 4 ends at byte 152"):
        refuse(name, generous, everything, 115)

    # Offset 0 or count 0 is a segment the file never held: fill, not a cut.
    sparse = _fake_page((0, 104, 5000, 112), (4, 4, 0, 4))
    refuse(name, sparse, everything, 116)

    # A table that stops before the grid does (a last row of tiles left
    # out) has nothing to test past its end.
    tiles = {"is_tiled": True, "tilelength": 2, "tilewidth": 2}
    short_table = _fake_page((100, 104), (4, 4), **tiles)
    refuse(name, short_table, everything, 108)
    with pytest.raises(ValueError, match="tile 2 of 2"):
        refuse(name, short_table, everything, 107)

    # A tile on the edge may hold only its part of the page: here the last
    # of a 3 x 3 page in 2 x 2 tiles is whole at four, two or one sample.
    edge = _fake_page(
        (100, 104, 108, 112), (4, 4, 4, 4), imagelength=3, imagewidth=3, **tiles
    )
    for held in (4, 2, 1):
        refuse(name, edge, (0, 3, 0, 3), 112 + held)
    with pytest.raises(ValueError, match="tile 4 of 4"):
        refuse(name, edge, (0, 3, 0, 3), 115)

    # An offset at the top of the 64-bit range must not wrap into the file.
    wild = _fake_page((2**64 - 8, 104, 108, 112), (16, 4, 4, 4))
    with pytest.raises(ValueError, match="strip 1 of 4"):
        refuse(name, wild, everything, 1000)

    # A page of a series read as one run sits whole pages after the first.
    refuse(name, page, everything, 216, shift=100)
    with pytest.raises(ValueError, match="strip 4 of 4 ends at byte 216"):
        refuse(name, page, everything, 215, shift=100)

    # Samples that are not whole bytes are not counted here: a segment
    # with some of its bytes is left to the decoder, one with none is not.
    bilevel = _fake_page((100, 101, 102, 103), (1, 1, 1, 9), bitspersample=1)
    refuse(name, bilevel, everything, 104)
    with pytest.raises(ValueError, match="strip 4 of 4"):
        refuse(name, bilevel, everything, 103)


@pytest.mark.parametrize("compression", [None, "lzw", "zlib"])
def test_a_byte_count_past_the_end_does_not_drop_a_whole_picture(
    tmp_path: Path, compression
):
    """Some writers declare more bytes than a strip has; the pixels are all there."""
    rng = np.random.default_rng(17)
    image = rng.integers(1, 256, size=(64, 48, 3), dtype=np.uint8)
    path = tmp_path / "generous.tif"
    tifffile.imwrite(str(path), image, rowsperstrip=8, compression=compression)
    with tifffile.TiffFile(path) as tif:
        page = tif.pages[0]
        tag = page.tags[279]  # StripByteCounts
        code = {3: "H", 4: "I", 16: "Q"}[int(tag.dtype)]
        last = int(page.databytecounts[-1])
        # The premise: the last strip is the last thing in the file.
        assert int(page.dataoffsets[-1]) + last == path.stat().st_size
        where = tag.valueoffset + (tag.count - 1) * struct.calcsize(code)
        order = tif.byteorder
    data = bytearray(path.read_bytes())
    struct.pack_into(order + code, data, where, last + 100)
    path.write_bytes(bytes(data))

    source = OpticalTiffSource.probe(path)
    for rows in (64, 8):
        np.testing.assert_array_equal(_rows(source, rows), np.moveaxis(image, -1, 0))


def _pages_on_disk(path: Path) -> List[List[int]]:
    with tifffile.TiffFile(path) as tif:
        return [[int(n) for n in page.dataoffsets] for page in tif.series[0].pages]


def test_a_cut_series_is_refused_by_the_pages_it_reads(tmp_path: Path):
    """Channels across pages: a cut in a page that is read refuses the read.

    Compressed, so that every page's header sits in front of its own
    pixels and survives the cut. (An uncompressed multi-page file keeps its
    later headers behind the pixels; cut, it is a shorter file with whole
    pages, not a page cut short.)
    """
    rng = np.random.default_rng(15)
    cyx = rng.integers(1, 65535, size=(3, 32, 24), dtype=np.uint16)
    path = tmp_path / "channels.tif"
    tifffile.imwrite(
        str(path),
        cyx,
        photometric="minisblack",
        metadata={"axes": "CYX"},
        rowsperstrip=8,
        compression="zlib",
    )
    source = OpticalTiffSource.probe(path)
    assert source.page_axes == "series" and source.shape == (3, 32, 24)
    np.testing.assert_array_equal(_rows(source, 8), cyx)

    pages = _pages_on_disk(path)
    assert len(pages) == 3 and all(len(strips) == 4 for strips in pages)
    # The third channel loses its rows from 16 on.
    _cut_at(path, pages[2][2], 2)
    source = OpticalTiffSource.probe(path)
    assert source.page_axes == "series" and source.shape == (3, 32, 24)
    with pytest.raises(ValueError, match="the file is cut short: strip 3 of 4"):
        list(source.bands(8))
    np.testing.assert_array_equal(_rows(source, 8, (0, 16, 0, 24)), cyx[:, :16])
    with pytest.raises(ValueError, match="the file is cut short"):
        list(source.bands(8, (0, 17, 0, 24)))


def test_a_cut_in_a_plane_that_is_not_read_does_not_refuse_the_image(tmp_path: Path):
    """Only the first z plane is stored, so only its pages are tested."""
    rng = np.random.default_rng(16)
    zcyx = rng.integers(1, 65535, size=(2, 3, 32, 24), dtype=np.uint16)
    path = tmp_path / "stack.tif"
    tifffile.imwrite(
        str(path),
        zcyx,
        photometric="minisblack",
        metadata={"axes": "ZCYX"},
        rowsperstrip=8,
        compression="zlib",
    )
    pages = _pages_on_disk(path)
    assert len(pages) == 6
    # The last page is the third channel of the second z plane.
    _cut_at(path, pages[5][1], 2)
    source = OpticalTiffSource.probe(path)
    assert source.page_axes == "series" and source.shape == (3, 32, 24)
    assert source.dropped_planes == 3
    np.testing.assert_array_equal(_rows(source, 8), zcyx[0])


def test_a_series_read_as_one_run_is_tested_where_its_pixels_are(tmp_path: Path):
    """Uncompressed pages that follow one another: tifffile reads them as one run.

    The later pages' headers sit behind the pixels, so a cut file has lost
    them. What matters is whether the pixels of the run are there.
    """
    rng = np.random.default_rng(18)
    cyx = rng.integers(1, 65535, size=(3, 32, 24), dtype=np.uint16)
    intact = tmp_path / "run.tif"
    tifffile.imwrite(
        str(intact), cyx, photometric="minisblack", metadata={"axes": "CYX"}
    )
    with tifffile.TiffFile(intact) as tif:
        first = tif.series[0].dataoffset
    # The premise: one run, three pages of 32 x 24 uint16 from ``first``.
    assert first is not None
    page_bytes = 32 * 24 * 2
    data = intact.read_bytes()
    assert len(data) > first + 3 * page_bytes

    # Every pixel is there, the headers behind them are not.
    whole = tmp_path / "pixels_only.tif"
    whole.write_bytes(data[: first + 3 * page_bytes + 4])
    source = OpticalTiffSource.probe(whole)
    assert source.page_axes == "series" and source.shape == (3, 32, 24)
    np.testing.assert_array_equal(_rows(source, 8), cyx)

    # The third channel ends five bytes into row 16.
    cut = tmp_path / "cut.tif"
    cut.write_bytes(data[: first + 2 * page_bytes + 16 * 24 * 2 + 5])
    source = OpticalTiffSource.probe(cut)
    assert source.page_axes == "series" and source.shape == (3, 32, 24)
    with pytest.raises(ValueError, match="the file is cut short"):
        list(source.bands(8))


def test_the_warning_says_the_file_is_cut_short(
    truncated_tiff: Path, tmp_path: Path, thyra_logs
):
    """The reason is Thyra's own, the same whatever the decoder would have done."""
    with thyra_logs("thyra.converters", logging.WARNING) as records:
        _, success = _convert(_mock_reader(truncated_tiff), tmp_path / "cut.zarr")
    assert success is True
    size = truncated_tiff.stat().st_size
    reasons = [r.message for r in records if "Failed to load optical" in r.message]
    assert len(reasons) == 1
    assert "truncated_0000.tiff: the file is cut short: strip " in reasons[0]
    assert "of 40 ends at byte " in reasons[0]
    assert f"and the file has {size} bytes" in reasons[0]
    assert reasons[0].endswith("dropping 'ds_optical_highres' from the store")


def test_same_name_from_different_stems_keeps_both_files(tmp_path: Path, thyra_logs):
    """a_0000 and b_0000 both clean to "highres": the second is numbered."""
    rng = np.random.default_rng(5)
    first = tmp_path / "a_0000.tif"
    second = tmp_path / "b_0000.tif"
    one = rng.integers(0, 256, size=(20, 24, 3), dtype=np.uint8)
    two = rng.integers(0, 256, size=(18, 22, 3), dtype=np.uint8)
    tifffile.imwrite(str(first), one)
    tifffile.imwrite(str(second), two)
    output_path = tmp_path / "out.zarr"
    with thyra_logs("thyra.converters", logging.WARNING) as records:
        converter, success = _convert(_mock_reader(first, second), output_path)
    assert success is True
    assert converter.optical.pending == {}
    assert not any("is replaced by" in r.message for r in records)
    _, levels, _ = _read_element(output_path, "ds_optical_highres")
    np.testing.assert_array_equal(levels["s0"][1], np.moveaxis(one, -1, 0))
    _, levels, _ = _read_element(output_path, "ds_optical_highres_1")
    np.testing.assert_array_equal(levels["s0"][1], np.moveaxis(two, -1, 0))


def test_same_stem_keeps_the_last_file(tmp_path: Path, thyra_logs):
    """The same picture saved as .tif and .jpg: last wins, with a warning."""
    rng = np.random.default_rng(5)
    first = tmp_path / "sample_0000.tif"
    second = tmp_path / "sample_0000.jpg"
    tifffile.imwrite(str(first), rng.integers(0, 256, size=(20, 24, 3), dtype=np.uint8))
    PILImage.fromarray(rng.integers(0, 256, size=(18, 22, 3), dtype=np.uint8)).save(
        second, quality=95
    )
    output_path = tmp_path / "out.zarr"
    with thyra_logs("thyra.converters", logging.WARNING) as records:
        converter, success = _convert(_mock_reader(first, second), output_path)
    assert success is True
    assert converter.optical.pending == {}
    assert any(
        "is replaced by sample_0000.jpg" in r.message
        and "maps to the same name" in r.message
        for r in records
    )
    group = zarr.open_group(str(output_path), mode="r", use_consolidated=False)
    assert [n for n in group["images"] if "optical" in n] == ["ds_optical_highres"]
    _, levels, _ = _read_element(output_path, "ds_optical_highres")
    assert levels["s0"][1].shape == (3, 18, 22)


# ---- what the store says about its optical images ----------------------


def _root_attrs(store: Path) -> dict:
    """The store's own top-level attributes, as a consumer reads them."""
    with open(store / "zarr.json") as f:
        return json.load(f)["attributes"]


def _convert_designating(primary_stem: str, reader, output_path: Path, **fields):
    """Convert with the .mis alignment image already resolved.

    On a real Bruker acquisition ``OpticalImages.primary_filename`` is set
    while the .mis is read; the mock reader has no .mis, so it is planted
    here. Extra ``fields`` are set on the same object -- that object owns
    this state, so planting it on the converter would only be ignored.
    """
    from thyra.converters.spatialdata.streaming_converter import (
        StreamingSpatialDataConverter,
    )

    converter = StreamingSpatialDataConverter(
        reader, output_path, dataset_id="ds", pixel_size_um=10.0, use_csc=True
    )
    converter.optical.primary_filename = primary_stem
    for name, value in fields.items():
        setattr(converter.optical, name, value)
    return converter, converter.convert()


def test_the_store_says_which_file_each_optical_element_came_from(
    tmp_path: Path,
) -> None:
    """The element name is derived; the filename is only in the root attrs."""
    rng = np.random.default_rng(3)
    highres = tmp_path / "SlideA_0000.tif"
    overview = tmp_path / "SlideA_deriv.tif"
    tifffile.imwrite(
        str(highres), rng.integers(0, 256, size=(20, 24, 3), dtype=np.uint8)
    )
    tifffile.imwrite(
        str(overview), rng.integers(0, 256, size=(10, 12, 3), dtype=np.uint8)
    )

    output_path = tmp_path / "sources.zarr"
    _, success = _convert(_mock_reader(highres, overview), output_path)
    assert success is True

    optical = _root_attrs(output_path)["optical_images"]
    assert optical["elements"] == {
        "ds_optical_highres": {"source_file": "SlideA_0000.tif"},
        "ds_optical_overview": {"source_file": "SlideA_deriv.tif"},
    }
    # Nothing designated one: the block still carries the provenance.
    assert optical["alignment_element"] is None
    # And it is in the root attrs because it cannot be on the element: the
    # ome-zarr writer keeps exactly these two keys, whatever it was handed.
    attrs, _, _ = _read_element(output_path, "ds_optical_highres")
    assert sorted(attrs) == ["ome", "spatialdata_attrs"]


def test_a_store_with_no_optical_image_omits_the_section(
    tmp_path: Path,
) -> None:
    """Omitted rather than written empty, like every other optional section."""
    output_path = tmp_path / "none.zarr"
    _, success = _convert(_mock_reader(), output_path)
    assert success is True
    assert "optical_images" not in _root_attrs(output_path)


def test_the_store_names_the_alignment_element_not_the_filename(
    tmp_path: Path,
) -> None:
    """Both places that state it name an element key that resolves."""
    import spatialdata

    rng = np.random.default_rng(4)
    primary = tmp_path / "SlideB_0000.tif"
    other = tmp_path / "SlideB_deriv.tif"
    tifffile.imwrite(
        str(primary), rng.integers(0, 256, size=(20, 24, 3), dtype=np.uint8)
    )
    tifffile.imwrite(str(other), rng.integers(0, 256, size=(10, 12, 3), dtype=np.uint8))

    output_path = tmp_path / "aligned.zarr"
    # An alignment matrix is what puts "global" in optical pixels, which is
    # the variant that fills reference_element in.
    _, success = _convert_designating(
        "slideb_0000",
        _mock_reader(primary, other),
        output_path,
        _tic_to_image=np.eye(3, dtype=np.float64),
    )
    assert success is True

    attrs = _root_attrs(output_path)
    assert attrs["optical_images"]["alignment_element"] == "ds_optical_highres"
    cs_global = attrs["coordinate_systems"]["global"]
    assert cs_global["unit"] == "pixel"
    assert cs_global["reference_element"] == "ds_optical_highres"

    # The point of naming the element: it resolves.
    sdata = spatialdata.SpatialData.read(str(output_path))
    assert cs_global["reference_element"] in sdata.images


def test_a_dropped_optical_image_leaves_no_trace_in_the_attrs(
    truncated_tiff: Path, tmp_path: Path, thyra_logs
) -> None:
    """Attrs naming an element the store does not hold are worse than none."""
    output_path = tmp_path / "dropped.zarr"
    with thyra_logs("thyra.converters", logging.WARNING) as records:
        _, success = _convert_designating(
            "truncated_0000",
            _mock_reader(truncated_tiff),
            output_path,
            _tic_to_image=np.eye(3, dtype=np.float64),
        )
    assert success is True
    # The premise: it was declared and then dropped, not never added.
    assert any(
        "dropping 'ds_optical_highres' from the store" in r.message for r in records
    )
    assert not (output_path / "images" / "ds_optical_highres").exists()

    attrs = _root_attrs(output_path)
    assert "optical_images" not in attrs
    # The alignment image went with it; "global" is still that image's
    # pixel grid, there is just no element here that is it.
    assert attrs["coordinate_systems"]["global"]["unit"] == "pixel"
    assert attrs["coordinate_systems"]["global"]["reference_element"] is None
    # The consolidated copy the reader actually uses agrees.
    consolidated = SpatialData.read(str(output_path))
    assert [k for k in consolidated.images if "optical" in k] == []
    assert (
        consolidated.attrs["coordinate_systems"]["global"]["reference_element"] is None
    )
