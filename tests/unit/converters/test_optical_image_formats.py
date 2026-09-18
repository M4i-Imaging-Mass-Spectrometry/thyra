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
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("spatialdata")
zarr = pytest.importorskip("zarr")
PILImage = pytest.importorskip("PIL.Image")

from spatialdata import SpatialData  # noqa: E402
from spatialdata.models import Image2DModel  # noqa: E402
from spatialdata.transformations import Scale  # noqa: E402

from thyra.converters.spatialdata.optical_image import (  # noqa: E402
    OpticalRasterSource,
    OpticalTiffSource,
    StreamedOpticalImage,
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
        "  <TeachPoint>10,10;1000,2000</TeachPoint>\n"
        "  <TeachPoint>40,10;1040,2000</TeachPoint>\n"
        "  <TeachPoint>10,40;1000,2040</TeachPoint>\n"
        '  <Area Name="region0" Type="0">\n'
        "    <Point>10,10</Point>\n"
        "    <Point>40,40</Point>\n"
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
    # The .mis draws one Area, so the alignment image in the store is that
    # Area's window of the JPEG (plus the margin, clamped to the 37 x 53
    # page), not the whole page; the store says where the window sits. The
    # bytes are still the JPEG's own, which is what this test is for.
    root = json.loads((output_path / "zarr.json").read_text(encoding="utf-8"))
    crop = root["attributes"]["optical_images"]["elements"]["rapiflex_optical_highres"][
        "crop"
    ]
    assert crop["full_size"] == [53, 37]
    x0, y0 = crop["origin"]
    width, height = crop["size"]
    assert (x0, y0, width, height) == (7, 7, 36, 30)
    whole = np.moveaxis(np.asarray(PILImage.open(expected)), -1, 0)
    np.testing.assert_array_equal(
        sdata.images["rapiflex_optical_highres"].values,
        whole[:, y0 : y0 + height, x0 : x0 + width],
    )
    # The image the .mis named is the alignment image; the other two are
    # carried along beside it, scaled into its pixel space.
    assert "rapiflex_optical_derived" in sdata.images
    assert "rapiflex_optical_slide_overview" in sdata.images
