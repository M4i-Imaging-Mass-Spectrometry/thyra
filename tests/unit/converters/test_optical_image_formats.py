"""Optical images that are not TIFF reach the store.

Bruker discovery used to glob ``*.tif``/``*.tiff`` only, so an acquisition
whose ``.mis`` names a ``.jpg`` -- which real Rapiflex slides do -- converted
to a store holding nothing but the generated TIC image, with no warning,
because there was nothing to warn about: no file had been found.

Three things are pinned here:

* :func:`probe_optical_source` routes by suffix and reads a JPEG, PNG or BMP
  header without decoding it, reporting the ``(c, y, x)`` shape and dtype the
  pixels will actually have -- including the modes that have to be converted
  to get a channel count at all (palette) and the ones that are not 8-bit;
* a declared-then-streamed store built from one of those equals what
  ``Image2DModel.parse`` and ``SpatialData.write`` produce from the decoded
  array, pyramid levels included, exactly as the TIFF path does;
* end to end: a Rapiflex acquisition whose ``.mis`` names a ``.jpg`` converts
  to a store with that image in it.
"""

from __future__ import annotations

import json
import logging
import struct
import zlib
from pathlib import Path
from typing import List

import numpy as np
import pytest

pytest.importorskip("spatialdata")
zarr = pytest.importorskip("zarr")
PILImage = pytest.importorskip("PIL.Image")
tifffile = pytest.importorskip("tifffile")

from PIL import ImageOps  # noqa: E402
from spatialdata import SpatialData  # noqa: E402
from spatialdata.models import Image2DModel  # noqa: E402
from spatialdata.transformations import Scale  # noqa: E402

from thyra.converters.spatialdata import optical_image  # noqa: E402
from thyra.converters.spatialdata.optical_image import (  # noqa: E402
    OpticalImages,
    OpticalRasterSource,
    OpticalTiffSource,
    OrientedOpticalSource,
    StreamedOpticalImage,
    orient,
    probe_optical_source,
)

CHUNKS = (1, 16, 16)
TRANSFORM = Scale([2.0, 0.5], axes=("x", "y"))
TRANSFORMATIONS = {"slide": TRANSFORM, "global": TRANSFORM}


def _rgb(height: int = 37, width: int = 53, seed: int = 3) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(0, 256, size=(height, width, 3), dtype=np.uint8)


def _save(array: np.ndarray, path: Path, mode: str | None = None) -> Path:
    image = PILImage.fromarray(array)
    if mode is not None:
        image = image.convert(mode)
    image.save(path)
    return path


@pytest.fixture
def rgb_png(tmp_path: Path) -> Path:
    """Lossless, so the decoded pixels are exactly what was written."""
    return _save(_rgb(), tmp_path / "scan_0000.png")


# ---- probing -----------------------------------------------------------


def test_probe_routes_by_suffix(tmp_path: Path, rgb_png: Path):
    tiff = tmp_path / "scan.tif"
    pytest.importorskip("tifffile").imwrite(str(tiff), _rgb())

    assert isinstance(probe_optical_source(tiff), OpticalTiffSource)
    assert isinstance(probe_optical_source(rgb_png), OpticalRasterSource)
    # Case of the suffix does not decide which decoder is used.
    upper = tmp_path / "scan_upper.TIF"
    upper.write_bytes(tiff.read_bytes())
    assert isinstance(probe_optical_source(upper), OpticalTiffSource)


@pytest.mark.parametrize("suffix", [".jpg", ".jpeg", ".png", ".bmp"])
def test_probe_reads_the_header_of_every_accepted_format(tmp_path: Path, suffix: str):
    source = probe_optical_source(_save(_rgb(), tmp_path / f"scan_0000{suffix}"))

    assert source.shape == (3, 37, 53)
    assert source.dtype == np.dtype("uint8")
    # One band, the whole page: these streams have no addressable row range.
    assert source.strip_rows == 37
    assert source.row_bytes == 3 * 53


@pytest.mark.parametrize(
    "mode, expected_channels, expected_dtype",
    [
        ("L", 1, "uint8"),
        ("RGB", 3, "uint8"),
        ("RGBA", 4, "uint8"),
        ("P", 3, "uint8"),  # palette: no channel count until it is converted
        ("I;16", 1, "uint16"),
    ],
)
def test_probe_reports_the_shape_the_pixels_will_have(
    tmp_path: Path, mode: str, expected_channels: int, expected_dtype: str
):
    """Whatever Pillow mode the file is in, the declaration matches the decode."""
    if mode == "I;16":
        base = (np.arange(37 * 53).reshape(37, 53) % 65535).astype(np.uint16)
        path = _save(base, tmp_path / "scan_0000.png")
    else:
        path = _save(_rgb(), tmp_path / "scan_0000.png", mode=mode)

    source = probe_optical_source(path)
    assert source.shape == (expected_channels, 37, 53)
    assert source.dtype == np.dtype(expected_dtype)

    ((_, decoded),) = source.bands(rows=source.shape[1])
    assert decoded.shape == source.shape
    assert decoded.dtype == source.dtype


def test_probe_refuses_a_file_that_is_not_an_image(tmp_path: Path):
    path = tmp_path / "scan_0000.jpg"
    path.write_bytes(b"not a jpeg, just bytes")

    # OSError, which UnidentifiedImageError is: the same class of failure
    # the TIFF probe raises, and the same per-image skip downstream.
    with pytest.raises(OSError):
        probe_optical_source(path)


def test_probe_does_not_refuse_a_large_image(tmp_path: Path, monkeypatch):
    """Pillow's decompression-bomb ceiling does not drop a big slide scan."""
    monkeypatch.setattr(PILImage, "MAX_IMAGE_PIXELS", 4)
    source = probe_optical_source(_save(_rgb(), tmp_path / "scan_0000.png"))

    assert source.shape == (3, 37, 53)
    ((_, decoded),) = source.bands(rows=37)
    assert decoded.shape == (3, 37, 53)
    # And the process-wide setting is left exactly as it was found.
    assert PILImage.MAX_IMAGE_PIXELS == 4


def test_decoded_pixels_are_the_file(rgb_png: Path):
    source = probe_optical_source(rgb_png)
    bands = list(source.bands(rows=1))  # asking for one row still gets the page

    assert [start for start, _ in bands] == [0]
    np.testing.assert_array_equal(
        bands[0][1], np.moveaxis(np.asarray(PILImage.open(rgb_png)), -1, 0)
    )


# ---- the store ---------------------------------------------------------


def _read_element(store: Path, name: str):
    group = zarr.open_group(str(store), mode="r", use_consolidated=False)[
        f"images/{name}"
    ]
    attrs = group.attrs.asdict()
    levels = {}
    for dataset in attrs["ome"]["multiscales"][0]["datasets"]:
        levels[dataset["path"]] = group[dataset["path"]][:]
    return attrs, levels


@pytest.mark.parametrize("scale_factors", [[2, 2], []])
def test_streamed_png_store_equals_parse_and_write(
    rgb_png: Path, tmp_path: Path, scale_factors
):
    """The non-TIFF source produces the store the whole-page route produced."""
    name = "optical"
    whole = np.moveaxis(np.asarray(PILImage.open(rgb_png)), -1, 0)

    kwargs = {"transformations": dict(TRANSFORMATIONS), "chunks": CHUNKS}
    if scale_factors:
        kwargs["scale_factors"] = scale_factors
    expected_store = tmp_path / "expected.zarr"
    SpatialData(
        images={name: Image2DModel.parse(whole, dims=("c", "y", "x"), **kwargs)}
    ).write(expected_store)

    streamed = StreamedOpticalImage(
        source=probe_optical_source(rgb_png),
        name=name,
        chunks=CHUNKS,
        scale_factors=scale_factors,
        transformations=dict(TRANSFORMATIONS),
    )
    actual_store = tmp_path / "actual.zarr"
    SpatialData(images={name: streamed.placeholder()}).write(actual_store)
    streamed.stream_pixels(actual_store)

    expected_attrs, expected_levels = _read_element(expected_store, name)
    actual_attrs, actual_levels = _read_element(actual_store, name)
    assert json.dumps(actual_attrs, sort_keys=True) == json.dumps(
        expected_attrs, sort_keys=True
    )
    assert list(actual_levels) == [f"s{i}" for i in range(len(scale_factors) + 1)]
    for path, expected_pixels in expected_levels.items():
        np.testing.assert_array_equal(actual_levels[path], expected_pixels)
    np.testing.assert_array_equal(actual_levels["s0"], whole)


def _convert_with_optical(optical: Path, output_path: Path):
    from tests.fixtures.mock_msi_generator import MockMSIConfig, MockMSIReader
    from thyra.converters.spatialdata.streaming_converter import (
        StreamingSpatialDataConverter,
    )

    reader = MockMSIReader(
        MockMSIConfig(n_x=4, n_y=3, n_mz_bins=32, peaks_per_spectrum=(2, 4)),
        optical_image_paths=[optical],
    )
    converter = StreamingSpatialDataConverter(
        reader, output_path, dataset_id="ds", pixel_size_um=10.0, use_csc=True
    )
    return converter, converter.convert()


def test_conversion_carries_a_jpeg_into_the_store(tmp_path: Path):
    jpeg = _save(_rgb(), tmp_path / "scan_0000.jpg")
    output_path = tmp_path / "jpeg.zarr"

    converter, success = _convert_with_optical(jpeg, output_path)

    assert success is True
    assert converter.optical.pending == {}
    sdata = SpatialData.read(str(output_path))
    stored = sdata.images["ds_optical_highres"]
    assert stored.shape == (3, 37, 53)
    assert stored.dtype == np.dtype("uint8")
    np.testing.assert_array_equal(
        stored.values, np.moveaxis(np.asarray(PILImage.open(jpeg)), -1, 0)
    )


def test_an_undecodable_jpeg_drops_the_image_not_the_conversion(
    tmp_path: Path, thyra_logs
):
    """Header intact, entropy-coded data cut: the tolerance TIFF already had."""
    jpeg = _save(_rgb(), tmp_path / "scan_0000.jpg")
    data = jpeg.read_bytes()
    jpeg.write_bytes(data[: len(data) // 3])
    assert probe_optical_source(jpeg).shape == (3, 37, 53)

    output_path = tmp_path / "truncated.zarr"
    # Not caplog: this test converts, and a conversion can reach setup_logging,
    # which sets propagate = False on the `thyra` logger process-wide. After
    # that caplog's root handler never sees another Thyra record, so a caplog
    # assertion here passes alone and fails in the full suite.
    with thyra_logs("thyra", logging.WARNING) as records:
        converter, success = _convert_with_optical(jpeg, output_path)

    assert success is True
    assert converter.optical.pending == {}
    assert any(
        "Failed to load optical image scan_0000.jpg" in record.getMessage()
        for record in records
    )
    sdata = SpatialData.read(str(output_path))
    assert [key for key in sdata.images if "optical" in key] == []
    assert "ds_z0_tic" in sdata.images


# ---- the whole way through, from a Rapiflex folder ---------------------


def _rapiflex_acquisition(folder: Path, optical_name: str) -> Path:
    """A 3x3 Rapiflex acquisition whose .mis names ``optical_name``.

    The smallest folder RapiflexReader accepts: a .dat holding a header, an
    offset table and nine spectra, the position log and info file that go
    with it, and a .mis naming the optical image and the Area it covers.
    """
    folder.mkdir(parents=True, exist_ok=True)
    n_spots, n_datapoints = 9, 64

    with open(folder / "sample_poslog.txt", "w", encoding="utf-8") as handle:
        handle.write("#Timestamp Pos X Y Z\n")
        for y in range(3):
            for x in range(3):
                handle.write(
                    f"2026-01-01 12:00:00.000 R00X{x}Y{y} "
                    f"{1000.0 + x * 20:.1f} {2000.0 + y * 20:.1f} 0.0\n"
                )

    (folder / "sample_info.txt").write_text(
        "Rapiflex Info File\n"
        "Name of Sample: Optical Format Fixture\n"
        f"Number of Spots: {n_spots}\n"
        "Number of Shots: 100\n"
        "Spectrum Size: 64\n"
        "Detector Gain: 2.0\n"
        "Mass Start: 100.0\n"
        "Mass End: 500.0\n"
        "Acquisition Mode: REFLECTOR\n"
        "Instrument Serial Number: TEST123\n"
        "Laser Power: 50\n"
        "Sample Rate: 1.0\n"
        f"DataPoints: {n_datapoints}\n"
        "Method: TestMethod.par\n"
        "flexImaging Version: 5.0.0\n"
        "flexControl Version: 4.0.0\n"
        "Raster: 20,20\n"
        "Start Time: Mon, 01.01.2026 12:00:00\n"
        "End Time: Mon, 01.01.2026 12:30:00\n",
        encoding="utf-8",
    )

    (folder / "sample.mis").write_text(
        '<?xml version="1.0"?>\n'
        '<ImagingSequence flexImagingVersion="5.0.0">\n'
        "  <Method>TestMethod.par</Method>\n"
        f"  <ImageFile>{optical_name}</ImageFile>\n"
        # Absolute, from the machine that acquired it: provenance only.
        f"  <OriginalImage>D:\\FlexImaging\\2026\\{optical_name}</OriginalImage>\n"
        "  <Raster>20,20</Raster>\n"
        # 2 um per photo pixel, stage y up the photo as FlexImaging has it,
        # so the 20 um raster is 10 px and spot (x, y) sits at photo pixel
        # (10 + 10x, 10 + 10y), on the lattice through the first teaching
        # point. The Area is drawn half a step outside the spots (D31).
        "  <TeachPoint>10,10;1000,2000</TeachPoint>\n"
        "  <TeachPoint>40,10;1060,2000</TeachPoint>\n"
        "  <TeachPoint>10,40;1000,1940</TeachPoint>\n"
        '  <Area Name="region0" Type="0">\n'
        "    <Point>5,5</Point>\n"
        "    <Point>35,35</Point>\n"
        "  </Area>\n"
        "</ImagingSequence>\n",
        encoding="utf-8",
    )

    header_size = 48
    data_start = header_size + n_spots * 4
    with open(folder / "sample.dat", "wb") as handle:
        handle.write(
            struct.pack(
                "<12I", header_size, 256, 0, 0, 3, 3, n_datapoints, 0, 0, 0, 0, 0
            )
        )
        for index in range(n_spots):
            handle.write(struct.pack("<I", data_start + index * n_datapoints * 4))
        for index in range(n_spots):
            spectrum = np.zeros(n_datapoints, dtype=np.float32)
            spectrum[10 + index] = 100.0 + index * 10
            spectrum[30] = 50.0
            handle.write(spectrum.tobytes())

    return folder


def test_rapiflex_acquisition_with_a_jpeg_optical_image(tmp_path: Path):
    """The reported defect: .mis names a .jpg, and the store now has it.

    On the glob this replaces the store came back with the generated TIC
    image alone.
    """
    from thyra.convert import convert_msi

    folder = _rapiflex_acquisition(tmp_path / "acquisition", "sample_0000.jpg")
    expected = _save(_rgb(seed=8), folder / "sample_0000.jpg")
    # Two more JPEGs in the same folder, as a real slide has: neither is the
    # one the .mis names.
    _save(_rgb(29, 31, seed=9), folder / "slide_overview.jpg")
    _save(_rgb(23, 19, seed=10), folder / "sample_0001.jpg")

    output_path = tmp_path / "rapiflex.zarr"
    assert convert_msi(
        str(folder), str(output_path), dataset_id="rapiflex", pixel_size_um=20.0
    )

    sdata = SpatialData.read(str(output_path))
    assert "rapiflex_optical_highres" in sdata.images
    # The .mis draws one Area, so the alignment image in the store is the
    # window its pixels cover, 5 to 35 on both axes (plus the margin, clamped
    # to the 37 x 53 page), not the whole page; the store says where the
    # window sits. The bytes are still the JPEG's own, which is what this
    # test is for.
    root = json.loads((output_path / "zarr.json").read_text(encoding="utf-8"))
    crop = root["attributes"]["optical_images"]["elements"]["rapiflex_optical_highres"][
        "crop"
    ]
    assert crop["full_size"] == [53, 37]
    x0, y0 = crop["origin"]
    width, height = crop["size"]
    assert (x0, y0, width, height) == (2, 2, 36, 35)
    whole = np.moveaxis(np.asarray(PILImage.open(expected)), -1, 0)
    np.testing.assert_array_equal(
        sdata.images["rapiflex_optical_highres"].values,
        whole[:, y0 : y0 + height, x0 : x0 + width],
    )
    # The image the .mis named is the alignment image; the other two are
    # carried along beside it, scaled into its pixel space when they are a
    # resampled copy of it (D59).
    assert "rapiflex_optical_derived" in sdata.images
    assert "rapiflex_optical_slide_overview" in sdata.images


def test_a_declined_alignment_puts_the_photo_on_the_tic(tmp_path: Path):
    """With the alignment declined, the photo, the TIC and ``obs`` agree.

    The fixture puts spot ``(x, y)`` at photo pixel ``(10 + 10x, 10 + 10y)``.
    Carried into micrometres by the photo's own transform, that point must be
    the row's ``obs`` position, and so must the centre of the TIC cell
    ``(x, y)`` (D32). Before D32 the TIC cell sat half a step off both.
    """
    from spatialdata.transformations import get_transformation

    from thyra.convert import convert_msi

    folder = _rapiflex_acquisition(tmp_path / "acquisition", "sample_0000.jpg")
    _save(_rgb(seed=8), folder / "sample_0000.jpg")
    output_path = tmp_path / "declined.zarr"
    assert convert_msi(
        str(folder),
        str(output_path),
        dataset_id="rapiflex",
        pixel_size_um=20.0,
        apply_optical_alignment=False,
    )

    sdata = SpatialData.read(str(output_path))
    root = json.loads((output_path / "zarr.json").read_text(encoding="utf-8"))
    assert root["attributes"]["coordinate_systems"]["global"]["unit"] == "micrometer"
    photo = "rapiflex_optical_highres"
    x0, y0 = root["attributes"]["optical_images"]["elements"][photo]["crop"]["origin"]

    def matrix(element):
        return np.asarray(
            get_transformation(element, "global").to_affine_matrix(
                input_axes=("x", "y"), output_axes=("x", "y")
            )
        )

    photo_to_um = matrix(sdata.images[photo])
    tic_to_um = matrix(sdata.images["rapiflex_z0_tic"])
    obs = sdata.tables["rapiflex_z0"].obs
    assert len(obs) == 9
    for _, row in obs.iterrows():
        x, y = int(row["x"]), int(row["y"])
        spot = photo_to_um @ [10 + 10 * x - x0, 10 + 10 * y - y0, 1.0]
        cell = tic_to_um @ [x + 0.5, y + 0.5, 1.0]
        position = [row["spatial_x"], row["spatial_y"]]
        np.testing.assert_allclose(spot[:2], position, atol=1e-9)
        np.testing.assert_allclose(cell[:2], position, atol=1e-9)


# ---- orientation, colour, series, 16-bit PNG and size (D69) -------------


class _PhotoReader:
    """A reader that holds photos and states nothing about them."""

    def __init__(self, paths: List[Path]) -> None:
        self._paths = list(paths)

    def get_optical_image_paths(self) -> List[Path]:
        return list(self._paths)

    def get_primary_optical_image_path(self):
        return None


def _declare(tmp_path: Path, *paths: Path) -> OpticalImages:
    """An ``OpticalImages`` over ``paths`` with every image declared."""
    optical = OpticalImages(
        _PhotoReader(list(paths)),
        lambda: tmp_path / "photos.zarr",
        "ds",
        pixel_size_xy=lambda: (10.0, 10.0),
        include=True,
        apply_alignment=True,
    )
    optical.compute_alignment()
    optical.add_images({})
    return optical


def _stored(tmp_path: Path, *paths: Path):
    """Declare, write and stream ``paths``: the state, the store, what was declared."""
    optical = _declare(tmp_path, *paths)
    declared = dict(optical.pending)  # streaming empties the pending images
    images = {name: image.placeholder() for name, image in declared.items()}
    SpatialData(images=images).write(tmp_path / "photos.zarr")
    optical.stream_pending_pixels()
    return optical, tmp_path / "photos.zarr", declared


def _level0(store: Path, name: str) -> np.ndarray:
    return _read_element(store, name)[1]["s0"]


def _cyx_of(image) -> np.ndarray:
    array = np.asarray(image)
    return array[np.newaxis] if array.ndim == 2 else np.moveaxis(array, -1, 0)


def _upright(raw: np.ndarray, tag: int) -> np.ndarray:
    """Pillow's own upright version of ``raw`` carrying orientation ``tag``."""
    image = PILImage.fromarray(raw)
    image.getexif()[274] = tag
    return _cyx_of(ImageOps.exif_transpose(image))


def _tagged_jpeg(path: Path, raw: np.ndarray, tag: int) -> Path:
    exif = PILImage.Exif()
    exif[274] = tag
    PILImage.fromarray(raw).save(path, quality=95, exif=exif.tobytes())
    return path


def _tagged_tiff(path: Path, raw: np.ndarray, tag: int) -> Path:
    tifffile.imwrite(str(path), raw, extratags=[(274, 3, 1, tag, True)])
    return path


def _photo(tmp_path: Path, kind: str, raw: np.ndarray, tag: int) -> Path:
    if kind == "jpeg":
        return _tagged_jpeg(tmp_path / "photo.jpg", raw, tag)
    path = _tagged_tiff(tmp_path / "photo.tif", raw, tag)
    if kind == "tiff_named_png":
        return path.rename(tmp_path / "photo.png")
    return path


@pytest.mark.parametrize("kind", ["jpeg", "tiff", "tiff_named_png"])
@pytest.mark.parametrize("tag", range(1, 9))
def test_a_photo_is_stored_upright_by_its_orientation_tag(
    tmp_path: Path, kind: str, tag: int
):
    raw = _rgb(6, 10, seed=tag)
    path = _photo(tmp_path, kind, raw, tag)
    if kind == "jpeg":
        # The decoded JPEG is the reference: lossy, but the same pixels.
        expected = _cyx_of(ImageOps.exif_transpose(PILImage.open(path)))
    else:
        expected = _upright(raw, tag)

    assert probe_optical_source(path).orientation == tag
    optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    streamed = declared[name].source
    # The declared shape is the upright one: a quarter turn swaps the sides.
    assert streamed.shape == expected.shape
    np.testing.assert_array_equal(_level0(store, name), expected)
    recorded = optical.root_attr()["elements"][name].get("orientation")
    if tag == 1:
        assert recorded is None
        assert not isinstance(streamed, OrientedOpticalSource)
    else:
        assert recorded == {"tag": tag, "applied": True}
        assert isinstance(streamed, OrientedOpticalSource)


def test_a_tiff_is_read_by_its_bytes_not_its_suffix(tmp_path: Path):
    named_tif = _tagged_tiff(tmp_path / "a.tif", _rgb(6, 10), 1)
    named_png = tmp_path / "a.png"
    named_png.write_bytes(named_tif.read_bytes())
    named_dat = tmp_path / "a.dat"
    named_dat.write_bytes(named_tif.read_bytes())

    for path in (named_tif, named_png, named_dat):
        assert isinstance(probe_optical_source(path), OpticalTiffSource)


@pytest.mark.parametrize("tag", range(1, 9))
def test_orient_equals_pillows_exif_transpose(tag: int):
    raw = _rgb(5, 9, seed=20 + tag)

    np.testing.assert_array_equal(
        orient(np.moveaxis(raw, -1, 0), tag), _upright(raw, tag)
    )
    # A grey page, which Pillow turns with the same rule.
    grey = raw[:, :, 0]
    np.testing.assert_array_equal(orient(grey[np.newaxis], tag), _upright(grey, tag))


def test_a_tag_outside_one_to_eight_is_left_as_encoded(tmp_path: Path, thyra_logs):
    raw = _rgb(6, 10)
    path = _tagged_tiff(tmp_path / "photo.tif", raw, 12)

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), np.moveaxis(raw, -1, 0))
    assert any("orientation tag 12" in r.getMessage() for r in records)
    assert "orientation" not in optical.root_attr()["elements"][name]


@pytest.mark.parametrize("tag", [3, 6])
def test_a_tagged_tiff_streams_its_bands_into_place(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, tag: int
):
    """Row bands of the stored image land as column or reversed row bands."""
    raw = _rgb(40, 24, seed=tag)
    path = tmp_path / "tagged.tif"
    tifffile.imwrite(str(path), raw, rowsperstrip=4, extratags=[(274, 3, 1, tag, True)])
    source = probe_optical_source(path)
    assert source.orientation == tag
    # One strip (4 rows of 72 bytes) a band: ten bands for 40 rows.
    monkeypatch.setattr(optical_image, "BAND_BUDGET_BYTES", 4 * 24 * 3)
    blocks: List[tuple] = []
    original = optical_image._source_blocks

    def spy(source, rows):
        for block in original(source, rows):
            blocks.append(block[:2])
            yield block

    monkeypatch.setattr(optical_image, "_source_blocks", spy)
    streamed = StreamedOpticalImage(
        source=OrientedOpticalSource(source, tag),
        name="optical",
        chunks=CHUNKS,
        scale_factors=[],
        transformations=dict(TRANSFORMATIONS),
    )
    store = tmp_path / "tagged.zarr"
    SpatialData(images={"optical": streamed.placeholder()}).write(store)
    streamed.stream_pixels(store)

    assert len(blocks) > 3
    expected = np.rot90(raw, -1) if tag == 6 else raw[::-1, ::-1]
    np.testing.assert_array_equal(
        _level0(store, "optical"), np.moveaxis(expected, -1, 0)
    )


def test_the_alignment_image_stays_as_encoded_with_a_warning(
    tmp_path: Path, thyra_logs
):
    """The teaching points may refer to the encoded pixels: nothing is turned."""
    from thyra.convert import convert_msi

    folder = _rapiflex_acquisition(tmp_path / "acquisition", "sample_0000.jpg")
    path = _tagged_jpeg(folder / "sample_0000.jpg", _rgb(seed=8), 6)
    output_path = tmp_path / "tagged.zarr"

    with thyra_logs("thyra", logging.WARNING) as records:
        assert convert_msi(
            str(folder), str(output_path), dataset_id="rapiflex", pixel_size_um=20.0
        )

    sdata = SpatialData.read(str(output_path))
    root = json.loads((output_path / "zarr.json").read_text(encoding="utf-8"))
    element = root["attributes"]["optical_images"]["elements"][
        "rapiflex_optical_highres"
    ]
    assert element["orientation"] == {"tag": 6, "applied": False}
    x0, y0 = element["crop"]["origin"]
    width, height = element["crop"]["size"]
    whole = _cyx_of(PILImage.open(path))
    np.testing.assert_array_equal(
        sdata.images["rapiflex_optical_highres"].values,
        whole[:, y0 : y0 + height, x0 : x0 + width],
    )
    named = [
        r.getMessage()
        for r in records
        if r.levelno == logging.WARNING and "orientation tag" in r.getMessage()
    ]
    assert len(named) == 1
    assert "sample_0000.jpg" in named[0] and "tag 6" in named[0]


def test_a_dotted_ome_name_resolves_and_does_not_clash(tmp_path: Path):
    """F138: ``x.ome.png`` is the alignment image; ``x.png`` is a second one."""
    from thyra.convert import convert_msi

    ome = "HE_section_04_brightfield_40x.ome.png"
    folder = _rapiflex_acquisition(tmp_path / "acquisition", ome)
    _save(_rgb(seed=8), folder / ome)
    output_path = tmp_path / "ome.zarr"
    assert convert_msi(
        str(folder), str(output_path), dataset_id="rapiflex", pixel_size_um=20.0
    )
    sdata = SpatialData.read(str(output_path))
    root = json.loads((output_path / "zarr.json").read_text(encoding="utf-8"))
    attrs = root["attributes"]
    reference = attrs["coordinate_systems"]["global"]["reference_element"]
    alignment = attrs["optical_images"]["alignment_element"]
    assert reference == alignment
    assert alignment in sdata.images

    plain = tmp_path / "plain.zarr"
    _save(_rgb(seed=9), folder / "HE_section_04_brightfield_40x.png")
    assert convert_msi(
        str(folder), str(plain), dataset_id="rapiflex", pixel_size_um=20.0
    )
    sdata = SpatialData.read(str(plain))
    attrs = json.loads((plain / "zarr.json").read_text(encoding="utf-8"))["attributes"]
    optical_names = [n for n in sdata.images if "optical" in n]
    assert len(optical_names) == 2
    assert attrs["optical_images"]["alignment_element"] in optical_names
    assert attrs["coordinate_systems"]["global"]["reference_element"] in optical_names


# ---- TIFF colour models (F20) -------------------------------------------


def test_a_palette_tiff_is_stored_as_rgb(tmp_path: Path):
    colours = _rgb(1, 16, seed=1)[0]  # 16 8-bit colours
    index = (np.arange(6 * 10).reshape(6, 10) % 16).astype(np.uint8)
    colormap = np.zeros((3, 256), np.uint16)
    colormap[:, :16] = colours.T.astype(np.uint16) * 257
    path = tmp_path / "palette.tif"
    tifffile.imwrite(str(path), index, photometric="palette", colormap=colormap)

    source = probe_optical_source(path)
    assert source.shape == (3, 6, 10)
    assert source.dtype == np.dtype("uint8")
    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    np.testing.assert_array_equal(
        _level0(store, name), np.moveaxis(colours[index], -1, 0)
    )


def test_a_palette_with_sixteen_bit_colours_keeps_them(tmp_path: Path):
    index = (np.arange(6 * 10).reshape(6, 10) % 4).astype(np.uint8)
    colormap = np.zeros((3, 256), np.uint16)
    colormap[:, :4] = [[1000, 2000, 3000, 65535], [5, 40000, 777, 0], [9, 8, 7, 6]]
    path = tmp_path / "palette16.tif"
    tifffile.imwrite(str(path), index, photometric="palette", colormap=colormap)

    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    stored = _level0(store, name)
    assert stored.dtype == np.dtype("uint16")
    np.testing.assert_array_equal(stored, colormap[:, index])


def test_a_miniswhite_tiff_is_stored_inverted(tmp_path: Path):
    grey = _rgb(6, 10)[:, :, 0]
    path = tmp_path / "white.tif"
    tifffile.imwrite(str(path), grey, photometric="miniswhite")

    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), (255 - grey)[np.newaxis])


def test_a_bilevel_miniswhite_tiff_is_stored_inverted(tmp_path: Path):
    bits = (np.arange(6 * 10).reshape(6, 10) % 3 == 0).astype(bool)
    path = tmp_path / "bilevel.tif"
    tifffile.imwrite(str(path), bits, photometric="miniswhite")

    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    stored = _level0(store, name)
    assert stored.dtype == np.dtype("bool")
    np.testing.assert_array_equal(stored, ~bits[np.newaxis])


def test_a_cmyk_tiff_is_stored_as_rgb_the_way_pillow_converts_it(tmp_path: Path):
    rng = np.random.default_rng(4)
    ink = rng.integers(0, 256, size=(6, 10, 4), dtype=np.uint8)
    path = tmp_path / "cmyk.tif"
    tifffile.imwrite(str(path), ink, photometric="separated")

    assert probe_optical_source(path).shape == (3, 6, 10)
    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    expected = _cyx_of(PILImage.fromarray(ink, "CMYK").convert("RGB"))
    np.testing.assert_array_equal(_level0(store, name), expected)


def test_an_uncompressed_ycbcr_tiff_is_stored_as_rgb(tmp_path: Path):
    rng = np.random.default_rng(6)
    samples = rng.integers(16, 240, size=(6, 10, 3), dtype=np.uint8)
    samples[0, :4] = [100, 128, 128]  # grey: chroma at its centre
    path = tmp_path / "ycbcr.tif"
    tifffile.imwrite(str(path), samples, photometric="ycbcr", subsampling=(1, 1))

    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    stored = _level0(store, name)
    np.testing.assert_array_equal(stored[:, 0, :4], np.full((3, 4), 100))
    expected = _cyx_of(PILImage.fromarray(samples, "YCbCr").convert("RGB"))
    assert np.abs(stored.astype(int) - expected.astype(int)).max() <= 1


@pytest.mark.parametrize("subsampling", [None, (2, 2)])
def test_uncompressed_ycbcr_without_four_four_four_is_refused(subsampling):
    """An absent YCbCrSubSampling tag means (2, 2) in TIFF 6.0, not (1, 1)."""

    class _Tags:
        def valueof(self, _key, default=None):
            return default

    class _Page:
        photometric = 6
        compression = 1
        samplesperpixel = 3
        extrasamples = ()
        dtype = np.dtype(np.uint8)
        bitspersample = 8
        tags = _Tags()

    page = _Page()
    page.subsampling = subsampling  # type: ignore[attr-defined]
    with pytest.raises(ValueError, match="subsampled"):
        optical_image.TiffColour.of(page, "ycbcr.tif")


def test_a_cielab_tiff_is_skipped_with_a_warning(tmp_path: Path, thyra_logs):
    path = tmp_path / "lab.tif"
    tifffile.imwrite(str(path), _rgb(6, 10), photometric="cielab")
    with pytest.raises(ValueError, match="CIELAB"):
        probe_optical_source(path)

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical = _declare(tmp_path, path)

    assert optical.pending == {}
    assert any(
        "Failed to load optical image lab.tif" in r.getMessage()
        and "CIELAB" in r.getMessage()
        for r in records
    )


# ---- multi-page series (F137) -------------------------------------------


def test_an_ome_tiff_with_channels_across_pages_is_stored_whole(tmp_path: Path):
    rng = np.random.default_rng(8)
    cyx = rng.integers(0, 256, size=(3, 7, 9), dtype=np.uint8)
    path = tmp_path / "channels.ome.tif"
    tifffile.imwrite(str(path), cyx, ome=True, metadata={"axes": "CYX"})

    assert probe_optical_source(path).shape == (3, 7, 9)
    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), cyx)


def test_a_z_stack_keeps_its_first_plane_with_a_warning(tmp_path: Path, thyra_logs):
    rng = np.random.default_rng(9)
    zyx = rng.integers(0, 256, size=(3, 7, 9), dtype=np.uint8)
    path = tmp_path / "stack.ome.tif"
    tifffile.imwrite(str(path), zyx, ome=True, metadata={"axes": "ZYX"})

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), zyx[:1])
    assert any(
        "only the first plane along Z is kept" in r.getMessage() for r in records
    )


# ---- every plane behind one page header (D88) ----------------------------
#
# tifffile's ``truncate=True`` writes one page header and every plane behind
# it, and so does ImageJ above 4 GB. A file that lost its later headers looks
# the same: tifffile writes them behind the pixels, so a cut removes them
# first. Each test checks what tifffile sees in its file before it asserts.


def _planes(shape, seed: int, dtype=np.uint8) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(1, 250, size=shape).astype(dtype)


def _seen_by_tifffile(path: Path):
    """``(kind, axes, shape, page headers)`` of the first series."""
    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        return series.kind, series.axes, tuple(series.shape), len(series.pages)


def _cut_after_planes(path: Path, planes: float, extra: int = 0) -> Path:
    """Cut ``path`` so it holds ``planes`` planes of pixels, plus ``extra`` bytes."""
    with tifffile.TiffFile(str(path)) as tif:
        page = tif.pages[0]
        first, plane = int(page.dataoffsets[0]), int(page.nbytes)
    path.write_bytes(path.read_bytes()[: first + int(planes * plane) + extra])
    return path


def _warnings_of(records) -> List[str]:
    return [r.getMessage() for r in records if r.levelno >= logging.WARNING]


@pytest.mark.parametrize(
    "shape, dtype, options",
    [
        ((3, 9, 11), np.uint8, {"imagej": True}),
        ((3, 9, 11), np.uint8, {"photometric": "minisblack"}),
        ((3, 9, 11), np.uint16, {"photometric": "minisblack"}),
        (
            (3, 9, 11),
            np.uint8,
            {"photometric": "minisblack", "metadata": {"axes": "CYX"}},
        ),
        ((3, 9, 11), np.float32, {"imagej": True, "byteorder": ">"}),
        ((3, 9, 11), np.uint8, {"photometric": "minisblack", "bigtiff": True}),
        ((2, 9, 11, 3), np.uint8, {"photometric": "rgb"}),
    ],
    ids=["imagej", "shaped", "shaped16", "shaped-axes", "imagej-be", "big", "rgb"],
)
def test_planes_behind_one_page_header_are_stored_like_planes_with_their_own(
    tmp_path: Path, thyra_logs, shape, dtype, options
):
    planes = _planes(shape, 21, dtype)
    one = tmp_path / "one_header.tif"
    each = tmp_path / "each_header.tif"
    tifffile.imwrite(str(one), planes, truncate=True, **options)
    tifffile.imwrite(str(each), planes, **options)
    assert _seen_by_tifffile(one)[2:] == (shape, 1)
    assert _seen_by_tifffile(each)[2:] == (shape, shape[0])

    expected = planes if planes.ndim == 3 else np.moveaxis(planes, -1, 1)
    expected = expected.reshape((-1,) + expected.shape[-2:])
    probed = probe_optical_source(one)
    reference = probe_optical_source(each)
    assert probed.shape == reference.shape == expected.shape
    assert (probed.page_axes, probed.series_index, probed.dropped_planes) == (
        reference.page_axes,
        reference.series_index,
        reference.dropped_planes,
    )
    assert probed.cut_short == reference.cut_short == ""

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, one)
    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), expected)
    assert _warnings_of(records) == []


def test_planes_behind_one_page_header_are_read_in_bands_and_windows(tmp_path: Path):
    planes = _planes((3, 40, 17), 22, np.uint16)
    path = tmp_path / "strips.tif"
    tifffile.imwrite(
        str(path), planes, photometric="minisblack", truncate=True, rowsperstrip=8
    )
    assert _seen_by_tifffile(path)[2:] == ((3, 40, 17), 1)

    source = OpticalTiffSource.probe(path)
    assert source.strip_rows == 8
    bands = list(source.bands(16))
    assert [start for start, _ in bands] == [0, 16, 32]
    np.testing.assert_array_equal(np.concatenate([b for _, b in bands], axis=1), planes)
    (window,) = list(source.bands(40, window=(5, 29, 3, 12)))
    np.testing.assert_array_equal(window[1], planes[:, 5:29, 3:12])


def test_a_z_stack_behind_one_page_header_keeps_its_first_plane_with_a_warning(
    tmp_path: Path, thyra_logs
):
    zyx = _planes((3, 7, 9), 23)
    path = tmp_path / "stack.tif"
    tifffile.imwrite(
        str(path), zyx, imagej=True, truncate=True, metadata={"axes": "ZYX"}
    )
    assert _seen_by_tifffile(path) == ("imagej", "ZYX", (3, 7, 9), 1)

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), zyx[:1])
    assert any(
        "only the first plane along Z is kept; 2 page(s) are left out" in message
        for message in _warnings_of(records)
    )


def test_channels_of_the_first_slice_are_kept_behind_one_page_header(tmp_path: Path):
    zcyx = _planes((2, 3, 7, 9), 24)
    path = tmp_path / "hyperstack.tif"
    tifffile.imwrite(
        str(path), zcyx, imagej=True, truncate=True, metadata={"axes": "ZCYX"}
    )
    assert _seen_by_tifffile(path) == ("imagej", "ZCYX", (2, 3, 7, 9), 1)

    probed = probe_optical_source(path)
    assert (probed.shape, probed.dropped_planes, probed.dropped_axes) == (
        (3, 7, 9),
        3,
        "Z",
    )
    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), zcyx[0])


@pytest.mark.parametrize("options", [{"imagej": True}, {"photometric": "minisblack"}])
def test_a_file_cut_behind_its_pixels_keeps_every_plane_and_says_so(
    tmp_path: Path, thyra_logs, options
):
    """The later page headers are gone and every pixel is still there."""
    planes = _planes((3, 9, 11), 25)
    path = tmp_path / "headers_lost.tif"
    tifffile.imwrite(str(path), planes, **options)
    assert _seen_by_tifffile(path)[2:] == ((3, 9, 11), 3)
    _cut_after_planes(path, 3)
    assert _seen_by_tifffile(path)[2:] == ((3, 9, 11), 1)

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), planes)
    (warning,) = _warnings_of(records)
    assert "headers_lost.tif: it points to a page header at byte" in warning
    assert "may have been cut short" in warning


def test_a_file_cut_inside_a_plane_it_declares_keeps_the_first_with_warnings(
    tmp_path: Path, thyra_logs
):
    planes = _planes((3, 9, 11), 26)
    path = tmp_path / "cut.tif"
    tifffile.imwrite(str(path), planes, photometric="minisblack")
    _cut_after_planes(path, 2.5)
    # tifffile still declares three planes, behind the one header it finds.
    assert _seen_by_tifffile(path) == ("shaped", "QYX", (3, 9, 11), 1)

    probed = probe_optical_source(path)
    assert (probed.shape, probed.page_axes, probed.dropped_planes) == (
        (1, 9, 11),
        "YX",
        2,
    )
    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), planes[:1])
    messages = _warnings_of(records)
    assert any(
        "cut.tif: its header declares 3 planes and the file holds the pixels of 2"
        in message
        for message in messages
    )
    assert any("2 page(s) are left out" in message for message in messages)


@pytest.mark.parametrize(
    "axes, shape, held, kept, dropped",
    [
        # Slices outermost: the first slice's channels are the first planes.
        ("ZCYX", (2, 3, 7, 9), 4, (slice(0, 1), slice(None)), 3),
        # Channels outermost: the last channel's first slice is plane 4.
        ("CZYX", (3, 2, 7, 9), 5, (slice(None), slice(0, 1)), 3),
        ("CZYX", (3, 2, 7, 9), 4, None, 5),
    ],
    ids=["first-planes-held", "last-plane-held", "last-plane-lost"],
)
def test_a_cut_file_keeps_its_channels_when_the_planes_read_are_held(
    tmp_path: Path, axes, shape, held, kept, dropped
):
    planes = _planes(shape, 27)
    path = tmp_path / "cut_stack.tif"
    tifffile.imwrite(
        str(path),
        planes,
        photometric="minisblack",
        truncate=True,
        metadata={"axes": axes},
    )
    _cut_after_planes(path, held)
    assert _seen_by_tifffile(path) == ("shaped", axes, shape, 1)

    probed = probe_optical_source(path)
    assert f"declares 6 planes and the file holds the pixels of {held}" in (
        probed.cut_short
    )
    assert probed.dropped_planes == dropped
    read = np.concatenate([band for _, band in probed.bands(4)], axis=1)
    expected = planes[0, :1] if kept is None else planes[kept].reshape(-1, 7, 9)
    np.testing.assert_array_equal(read, expected)


def test_an_imagej_file_cut_inside_its_pixels_keeps_the_first_page_and_says_so(
    tmp_path: Path, thyra_logs
):
    planes = _planes((3, 9, 11), 28)
    path = tmp_path / "hyperstack_cut.tif"
    tifffile.imwrite(str(path), planes, imagej=True)
    _cut_after_planes(path, 2.5)
    # tifffile sets the ImageJ header aside: a whole one-page TIFF is left.
    assert _seen_by_tifffile(path) == ("generic", "YX", (9, 11), 1)

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), planes[:1])
    (warning,) = _warnings_of(records)
    assert (
        "its ImageJ header declares 3 images and the file holds the pixels of 2"
        in warning
    )


@pytest.mark.parametrize("bigtiff", [False, True], ids=["classic", "bigtiff"])
def test_pages_lost_from_a_file_with_no_description_are_told_by_the_header_chain(
    tmp_path: Path, thyra_logs, bigtiff
):
    planes = _planes((3, 9, 11), 29)
    path = tmp_path / "plain_cut.tif"
    tifffile.imwrite(
        str(path), planes, photometric="minisblack", metadata=None, bigtiff=bigtiff
    )
    assert _seen_by_tifffile(path)[2:] == ((3, 9, 11), 3)
    assert probe_optical_source(path).cut_short == ""
    _cut_after_planes(path, 1.5)
    assert _seen_by_tifffile(path)[2:] == ((9, 11), 1)

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), planes[:1])
    (warning,) = _warnings_of(records)
    assert "plain_cut.tif: it points to a page header at byte" in warning


@pytest.mark.parametrize(
    "options",
    [
        {"photometric": "rgb"},
        {"photometric": "rgb", "bigtiff": True},
        {"photometric": "rgb", "compression": "zlib", "tile": (16, 16)},
    ],
    ids=["strips", "bigtiff", "tiles"],
)
def test_a_whole_file_is_not_called_cut_short(tmp_path: Path, thyra_logs, options):
    path = tmp_path / "whole.tif"
    tifffile.imwrite(str(path), _rgb(40, 48), **options)
    pages = tmp_path / "pages.ome.tif"
    tifffile.imwrite(
        str(pages), _planes((3, 7, 9), 30), ome=True, metadata={"axes": "CYX"}
    )

    assert probe_optical_source(path).cut_short == ""
    assert probe_optical_source(pages).cut_short == ""
    with thyra_logs("thyra.converters", logging.WARNING) as records:
        _declare(tmp_path, path, pages)
    assert _warnings_of(records) == []


def test_a_buffer_error_from_tifffile_becomes_a_reason_in_words(
    tmp_path: Path, monkeypatch
):
    path = tmp_path / "photo.tif"
    tifffile.imwrite(str(path), _rgb(6, 10))

    def short_read(cls, path):
        raise struct.error("unpack requires a buffer of 12 bytes")

    monkeypatch.setattr(OpticalTiffSource, "_read_header", classmethod(short_read))
    with pytest.raises(ValueError, match="photo.tif: the file ends inside a page"):
        probe_optical_source(path)


def test_a_file_that_ends_inside_a_page_header_is_skipped_in_words(
    tmp_path: Path, thyra_logs
):
    path = tmp_path / "header_cut.tif"
    tifffile.imwrite(str(path), _planes((3, 9, 11), 31), imagej=True)
    _cut_after_planes(path, 3, extra=4)
    try:
        _seen_by_tifffile(path)
    except struct.error:
        pass
    else:
        pytest.skip("this tifffile reads a page header the file ends inside")

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical = _declare(tmp_path, path)

    assert optical.pending == {}
    (warning,) = _warnings_of(records)
    assert "Failed to load optical image header_cut.tif" in warning
    assert "the file ends inside a page header" in warning


def _first_page_header(path: Path):
    """``(offset, entries, {tag code: entry offset})`` of the first page header."""
    with tifffile.TiffFile(str(path)) as tif:
        page = tif.pages[0]
        entries = {tag.code: tag.offset for tag in page.tags.values()}
        return int(page.offset), len(page.tags), entries


def _point_past_the_last_header(path: Path, target: int) -> Path:
    """Make a classic TIFF's only page header name a next one at ``target``."""
    offset, entries, _ = _first_page_header(path)
    data = bytearray(path.read_bytes())
    pointer = offset + 2 + entries * 12
    data[pointer : pointer + 4] = struct.pack("<I", target)
    path.write_bytes(bytes(data))
    return path


def test_an_ome_tiff_that_also_carries_an_imagej_description_is_whole(
    tmp_path: Path, thyra_logs
):
    """Two descriptions in one page header, as Micro-Manager writes them."""
    planes = _planes((3, 9, 11), 32)
    path = tmp_path / "two_descriptions.ome.tif"
    note = "ImageJ=1.53\nimages=3\nchannels=3\nhyperstack=true\n"
    tifffile.imwrite(
        str(path),
        planes,
        ome=True,
        metadata={"axes": "CYX"},
        extratags=[(271, "s", 0, note, False)],
    )
    # Written as tag 271 so tifffile accepts it, then renumbered to 270.
    _, _, entries = _first_page_header(path)
    data = bytearray(path.read_bytes())
    data[entries[271] : entries[271] + 2] = struct.pack("<H", 270)
    path.write_bytes(bytes(data))
    with tifffile.TiffFile(str(path)) as tif:
        assert tif.is_imagej and tif.series[0].kind == "ome"

    assert probe_optical_source(path).cut_short == ""
    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), planes)
    assert _warnings_of(records) == []


def test_one_page_with_an_imagej_count_it_cannot_hold_is_told_as_it_is(
    tmp_path: Path, thyra_logs
):
    """Nothing says the file was cut: the warning states the two numbers."""
    plane = _planes((9, 11), 33)
    path = tmp_path / "one_slice.tif"
    tifffile.imwrite(
        str(path), plane, description="ImageJ=1.53\nimages=3\nslices=3\n", metadata=None
    )
    assert _seen_by_tifffile(path)[1:] == ("YX", (9, 11), 1)

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical, store, declared = _stored(tmp_path, path)

    (name,) = declared
    np.testing.assert_array_equal(_level0(store, name), plane[np.newaxis])
    (warning,) = _warnings_of(records)
    assert (
        "its ImageJ header declares 3 images and the file holds the pixels of 1"
        in warning
    )
    assert "may have been cut short" in warning


@pytest.mark.parametrize(
    "target, told",
    [
        (1_000_000, "it points to a page header at byte 1000000 and ends at byte"),
        # Inside the pixels: the bytes there are no page header, so they
        # say nothing about what the file holds. tifffile drops them too.
        (300, ""),
    ],
    ids=["past-the-end", "into-the-pixels"],
)
def test_a_stray_pointer_behind_the_last_page_header(tmp_path: Path, target, told):
    path = tmp_path / "stray.tif"
    rgb = np.full((40, 32, 3), 238, dtype=np.uint8)
    tifffile.imwrite(str(path), rgb, photometric="rgb")
    whole = probe_optical_source(path)
    _point_past_the_last_header(path, target)

    probed = probe_optical_source(path)
    assert told in probed.cut_short and bool(told) == bool(probed.cut_short)
    assert probed.shape == whole.shape == (3, 40, 32)
    (band,) = list(probed.bands(40))
    np.testing.assert_array_equal(band[1], np.moveaxis(rgb, -1, 0))


def test_a_trailing_axis_of_one_is_counted_like_planes_with_their_own_headers(
    tmp_path: Path,
):
    """tifffile keeps that axis in some versions and drops it in others."""
    planes = _planes((3, 9, 11, 1), 34)
    one = tmp_path / "unit_one.tif"
    each = tmp_path / "unit_each.tif"
    tifffile.imwrite(str(one), planes, photometric="minisblack", truncate=True)
    tifffile.imwrite(str(each), planes, photometric="minisblack")
    assert _seen_by_tifffile(one)[3] == 1
    assert _seen_by_tifffile(each)[3] == 3

    probed = probe_optical_source(one)
    reference = probe_optical_source(each)
    assert (probed.shape, probed.dropped_planes, probed.dropped_axes) == (
        reference.shape,
        reference.dropped_planes,
        reference.dropped_axes,
    )
    # Every plane is stored, or the ones left out are counted: never neither.
    assert probed.shape[0] == 3 or probed.dropped_planes == 2


def test_a_metamorph_stack_keeps_the_planes_behind_its_one_page_header(
    tmp_path: Path,
):
    """Built by hand: one page with the two UIC tags, three planes behind it."""
    planes = _planes((3, 9, 11), 35, np.uint16)
    path = tmp_path / "stack.stk"
    per_plane = [0, 1, 2450000, 1000, 2450000, 1000]
    tifffile.imwrite(
        str(path),
        planes[0],
        photometric="minisblack",
        metadata=None,
        extratags=[
            (33628, "I", 2, [0, 1], False),
            (33629, "I", 18, per_plane * 3, False),
        ],
    )
    _, _, entries = _first_page_header(path)
    data = bytearray(path.read_bytes())
    uic1, uic2 = entries[33628], entries[33629]
    data[uic1 + 4 : uic1 + 8] = struct.pack("<I", 1)
    # UIC2 as six rationals a plane: its count is the number of planes.
    data[uic2 + 2 : uic2 + 4] = struct.pack("<H", 5)
    data[uic2 + 4 : uic2 + 8] = struct.pack("<I", 3)
    path.write_bytes(bytes(data) + planes[1].tobytes() + planes[2].tobytes())
    assert _seen_by_tifffile(path) == ("stk", "IYX", (3, 9, 11), 1)

    probed = probe_optical_source(path)
    assert (probed.shape, probed.cut_short) == ((3, 9, 11), "")
    read = np.concatenate([band for _, band in probed.bands(4)], axis=1)
    np.testing.assert_array_equal(read, planes)


def test_a_one_header_file_cut_after_its_probe_is_refused_at_the_read(tmp_path: Path):
    """The read tests every plane it needs, not only the one with a header."""
    planes = _planes((3, 9, 11), 36)
    path = tmp_path / "cut_later.tif"
    tifffile.imwrite(str(path), planes, photometric="minisblack", truncate=True)
    source = OpticalTiffSource.probe(path)
    assert (source.shape, source.cut_short) == ((3, 9, 11), "")

    _cut_after_planes(path, 2.5)
    with pytest.raises(ValueError, match="cut short"):
        list(source.bands(9))


# ---- 16-bit PNGs (F136) --------------------------------------------------


def _png_chunk(kind: bytes, body: bytes) -> bytes:
    return (
        struct.pack(">I", len(body))
        + kind
        + body
        + struct.pack(">I", zlib.crc32(kind + body))
    )


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _write_png16(path: Path, samples: np.ndarray, colour_type: int) -> Path:
    """A 16-bit PNG of ``samples`` ``(h, w, c)``, written by hand."""
    height, width = samples.shape[:2]
    rows = b"".join(b"\x00" + samples[y].astype(">u2").tobytes() for y in range(height))
    header = struct.pack(">IIBBBBB", width, height, 16, colour_type, 0, 0, 0)
    path.write_bytes(
        PNG_SIGNATURE
        + _png_chunk(b"IHDR", header)
        + _png_chunk(b"IDAT", zlib.compress(rows))
        + _png_chunk(b"IEND", b"")
    )
    return path


@pytest.mark.parametrize(
    "colour_type, channels",
    [(0, 1), (2, 3), (4, 2), (6, 4)],
    ids=["L", "RGB", "LA", "RGBA"],
)
def test_a_sixteen_bit_png_keeps_its_depth_and_channels(
    tmp_path: Path, colour_type: int, channels: int
):
    rng = np.random.default_rng(colour_type)
    samples = rng.integers(0, 65536, size=(6, 10, channels), dtype=np.uint16)
    samples[0, 0] = 65535
    samples[0, 1] = 257  # distinct high and low bytes
    path = _write_png16(tmp_path / "deep.png", samples, colour_type)

    source = probe_optical_source(path)
    assert source.shape == (channels, 6, 10)
    assert source.dtype == np.dtype("uint16")
    optical, store, declared = _stored(tmp_path, path)
    (name,) = declared
    stored = _level0(store, name)
    assert stored.dtype == np.dtype("uint16")
    np.testing.assert_array_equal(stored, np.moveaxis(samples, -1, 0))


# ---- a header that promises more than the file can hold (F164) -----------


def _forged_png(path: Path, width: int, height: int) -> Path:
    header = struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0)
    path.write_bytes(
        PNG_SIGNATURE
        + _png_chunk(b"IHDR", header)
        + _png_chunk(b"IDAT", zlib.compress(b"\x00" * 4, 9))
        + _png_chunk(b"IEND", b"")
    )
    return path


def test_a_png_that_declares_more_than_it_can_hold_is_refused(tmp_path: Path):
    path = _forged_png(tmp_path / "forged.png", 40000, 40000)
    assert path.stat().st_size < 100

    with pytest.raises(ValueError, match="can hold"):
        probe_optical_source(path)


def test_the_warning_names_the_forged_png(tmp_path: Path, thyra_logs):
    path = _forged_png(tmp_path / "forged.png", 40000, 40000)

    with thyra_logs("thyra.converters", logging.WARNING) as records:
        optical = _declare(tmp_path, path)

    assert optical.pending == {}
    assert any(
        "forged.png" in r.getMessage() and "can hold" in r.getMessage() for r in records
    )


def test_a_small_honest_png_is_not_refused(tmp_path: Path):
    path = _save(np.zeros((64, 64, 3), np.uint8), tmp_path / "flat.png")

    assert probe_optical_source(path).shape == (3, 64, 64)
