"""Numbers to pixels: the colour map, the scale and the figure."""

from __future__ import annotations

import io

import numpy as np
import pytest
from PIL import Image

from thyra.quicklook.ion_image import IonImage, MzWindow
from thyra.quicklook.render import (
    FIGURE_LONG_SIDE,
    NOT_MEASURED,
    VIRIDIS,
    colour_bar,
    colourise,
    format_intensity,
    integer_scale,
    ion_image_figure,
    png_bytes,
    raster_image,
    scale_top,
)

NAN = np.nan


def _ion(image: np.ndarray) -> IonImage:
    window = MzWindow(760.513, 10.0, 4, 7, 760.5077, 760.5177)
    measured = image[np.isfinite(image)]
    return IonImage(
        "msi_dataset_z0",
        window,
        image,
        float(measured.min()),
        float(measured.max()),
        int(measured.size),
    )


class TestTheColourMap:
    def test_it_is_viridis_from_end_to_end(self):
        assert VIRIDIS.shape == (256, 3)
        assert VIRIDIS.dtype == np.uint8
        assert tuple(VIRIDIS[0]) == (0x44, 0x01, 0x54)
        assert tuple(VIRIDIS[255]) == (0xFD, 0xE7, 0x25)

    def test_it_only_ever_gets_lighter(self):
        """What makes it readable in grey, and what a rainbow map fails."""
        luminance = VIRIDIS.astype(float) @ np.array([0.2126, 0.7152, 0.0722])

        assert np.all(np.diff(luminance) > -0.5)
        assert luminance[-1] - luminance[0] > 150

    def test_no_colour_of_it_can_be_taken_for_not_measured(self):
        distance = np.abs(VIRIDIS.astype(int) - np.array(NOT_MEASURED)).sum(axis=1)

        assert distance.min() > 60


class TestScaleTop:
    def test_it_is_the_99th_percentile_of_the_pixels_with_signal(self):
        image = np.concatenate([np.zeros(500), np.arange(1.0, 101.0), [NAN] * 7])

        assert scale_top(image) == pytest.approx(99.01)

    def test_one_very_bright_pixel_does_not_set_the_scale(self):
        image = np.append(np.full(999, 10.0), 1e6)

        assert scale_top(image) == pytest.approx(10.0)

    def test_no_signal_is_zero(self):
        assert scale_top(np.array([0.0, 0.0, NAN])) == 0.0


class TestColourise:
    def test_zero_top_and_above_top_take_the_ends_of_the_map(self):
        pixels = colourise(np.array([[0.0, 50.0, 100.0, 400.0]]), top=100.0)

        assert tuple(pixels[0, 0]) == tuple(VIRIDIS[0])
        assert tuple(pixels[0, 1]) == tuple(VIRIDIS[128])
        assert tuple(pixels[0, 2]) == tuple(VIRIDIS[255])
        assert tuple(pixels[0, 3]) == tuple(VIRIDIS[255])

    def test_not_measured_is_drawn_apart_from_zero(self):
        pixels = colourise(np.array([[0.0, NAN]]), top=100.0)

        assert tuple(pixels[0, 1]) == NOT_MEASURED
        assert tuple(pixels[0, 0]) != NOT_MEASURED

    def test_an_image_without_signal_is_the_lowest_colour_throughout(self):
        pixels = colourise(np.zeros((2, 2)), top=0.0)

        assert {tuple(pixel) for pixel in pixels.reshape(-1, 3)} == {tuple(VIRIDIS[0])}


class TestRasterImage:
    @pytest.mark.parametrize(
        "height, width, factor",
        [(36, 48, 10), (480, 100, 1), (481, 10, 1), (3, 4, 120), (300, 300, 2)],
    )
    def test_the_factor_is_whole_and_reaches_the_long_side(self, height, width, factor):
        assert integer_scale(height, width, 480) == factor

    def test_a_small_raster_is_not_a_postage_stamp(self):
        pixels = colourise(np.arange(36.0 * 48).reshape(36, 48), top=1727.0)

        image = raster_image(pixels, FIGURE_LONG_SIDE)

        assert image.mode == "RGB"
        assert image.size == (480, 360)

    def test_each_pixel_becomes_one_block_of_its_own_colour(self):
        """Nearest neighbour: no colour appears that was not in the raster."""
        pixels = colourise(np.array([[0.0, 100.0], [NAN, 50.0]]), top=100.0)

        enlarged = np.asarray(raster_image(pixels, 8))

        assert enlarged.shape == (8, 8, 3)
        for row in range(2):
            for column in range(2):
                block = enlarged[row * 4 : row * 4 + 4, column * 4 : column * 4 + 4]
                assert (block == pixels[row, column]).all()

    def test_a_figure_keeps_every_pixel_of_a_large_raster(self):
        pixels = np.zeros((1000, 500, 3), dtype=np.uint8)

        assert raster_image(pixels, 480).size == (500, 1000)

    def test_the_report_reduces_a_large_raster(self):
        pixels = np.zeros((1000, 500, 3), dtype=np.uint8)

        assert raster_image(pixels, 480, reduce=True).size == (240, 480)


class TestTheFigure:
    def test_it_is_an_rgb_image_that_holds_the_enlarged_raster(self):
        figure = ion_image_figure(_ion(np.ones((36, 48))), "sample, msi_dataset_z0")

        assert figure.mode == "RGB"
        assert figure.width > 480 and figure.height > 360

    def test_the_raster_is_in_it_pixel_for_pixel(self):
        image = np.arange(12.0).reshape(3, 4)
        figure = np.asarray(ion_image_figure(_ion(image), "sample"))
        raster = np.asarray(raster_image(colourise(image, scale_top(image)), 480))

        height, width = raster.shape[:2]
        found = any(
            (figure[top : top + height, 16 : 16 + width] == raster).all()
            for top in range(figure.shape[0] - height + 1)
        )
        assert found

    def test_the_colour_bar_carries_the_whole_map_top_colour_first(self):
        bar = np.asarray(colour_bar(256))

        assert tuple(bar[0, 0]) == tuple(VIRIDIS[255])
        assert tuple(bar[-1, 0]) == tuple(VIRIDIS[0])
        figure = np.asarray(ion_image_figure(_ion(np.ones((36, 48))), "sample"))
        colours = {tuple(pixel) for pixel in figure.reshape(-1, 3)}
        assert {tuple(VIRIDIS[0]), tuple(VIRIDIS[128]), tuple(VIRIDIS[255])} <= colours

    def test_the_legend_for_not_measured_appears_only_when_needed(self):
        whole = np.ones((36, 48))
        holed = whole.copy()
        holed[0, 0] = NAN

        with_legend = ion_image_figure(_ion(holed), "sample")
        without = ion_image_figure(_ion(whole), "sample")

        assert with_legend.height > without.height
        below_raster = np.asarray(with_legend)[-30:]
        assert (below_raster == NOT_MEASURED).all(axis=2).any()

    def test_a_line_scan_keeps_a_colour_bar_that_can_be_read(self):
        """One row of pixels: the bar and its two labels still get their room."""
        figure = np.asarray(ion_image_figure(_ion(np.ones((1, 600))), "sample"))

        colours = {tuple(pixel) for pixel in figure.reshape(-1, 3)}
        assert {tuple(VIRIDIS[0]), tuple(VIRIDIS[255])} <= colours
        bar = (figure == VIRIDIS[0]).all(axis=2).any(axis=0)
        assert np.count_nonzero((figure[:, bar] != 255).any(axis=(1, 2))) >= 40

    def test_a_raster_with_no_pixels_has_a_factor_all_the_same(self):
        assert integer_scale(0, 0, 480) >= 1

    def test_it_is_written_as_a_png(self):
        data = png_bytes(ion_image_figure(_ion(np.ones((36, 48))), "sample"))

        assert data[:8] == b"\x89PNG\r\n\x1a\n"
        reopened = Image.open(io.BytesIO(data))
        assert reopened.format == "PNG" and reopened.mode == "RGB"


@pytest.mark.parametrize(
    "value, text",
    [
        (0.0, "0"),
        (5384.2, "5380"),
        (12.0, "12"),
        (0.01234, "0.0123"),
        (1.2e7, "1.2e+07"),
    ],
)
def test_an_intensity_is_three_digits_without_an_exponent_where_it_fits(value, text):
    assert format_intensity(value) == text
