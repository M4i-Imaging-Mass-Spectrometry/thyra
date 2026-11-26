"""Numbers to pixels: the colour map, the scale and the figure."""

from __future__ import annotations

import io

import numpy as np
import pytest
from PIL import Image

from thyra.quicklook.ion_image import IonImage, MzWindow
from thyra.quicklook.render import (
    FIGURE_LONG_SIDE,
    MAX_PIXEL_RATIO,
    NOT_MEASURED,
    SQUARE,
    VIRIDIS,
    colour_bar,
    colourise,
    figure_lines,
    format_intensity,
    integer_scale,
    ion_image_figure,
    pixel_size_note,
    pixel_stretch,
    png_bytes,
    raster_image,
    scale_top,
    stretched_size,
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


def _blocks(image: np.ndarray, axis: int) -> np.ndarray:
    """The length of each run of equal lines along ``axis``."""
    lines = np.moveaxis(image, axis, 0)
    edges = np.flatnonzero((lines[1:] != lines[:-1]).any(axis=(1, 2))) + 1
    return np.diff(np.concatenate(([0], edges, [lines.shape[0]])))


class TestPixelStretch:
    @pytest.mark.parametrize(
        "pixel_size, stretch",
        [
            ((70.0, 150.0), (1.0, 150 / 70)),
            ((150.0, 70.0), (150 / 70, 1.0)),
            ((1.0, 5000.0), (1.0, MAX_PIXEL_RATIO)),
            ((5000.0, 1.0), (MAX_PIXEL_RATIO, 1.0)),
            ((1e-300, 1e300), (1.0, MAX_PIXEL_RATIO)),
        ],
    )
    def test_the_shorter_side_is_the_unit_and_the_ratio_is_bounded(
        self, pixel_size, stretch
    ):
        assert pixel_stretch(pixel_size) == pytest.approx(stretch)

    @pytest.mark.parametrize(
        "pixel_size",
        [None, (20.0, 20.0), (20.0, 20.0 * (1 + 1e-12)), (10.0, 10.000001)],
    )
    def test_square_unknown_and_float_noise_are_square(self, pixel_size):
        assert pixel_stretch(pixel_size) == SQUARE
        assert pixel_size_note(pixel_size) is None

    @pytest.mark.parametrize(
        "pixel_size",
        [
            (10.0, 10.00011),
            (1.000004999, 1.000015),
            (0.9999995, 1.00001),
            (99998.0, 1e5),
        ],
    )
    def test_two_sizes_that_count_as_different_never_print_alike(self, pixel_size):
        """Six digits are printed; a smaller difference is square instead."""
        note = pixel_size_note(pixel_size)

        assert note is not None
        x, y = note.split(" um")[0].split(" x ")
        assert x != y

    def test_the_note_states_both_sides(self):
        assert pixel_size_note((70.0, 150.0)) == "70 x 150 um, drawn to scale"
        assert pixel_size_note((12.5, 33.3), "µm") == "12.5 x 33.3 µm, drawn to scale"

    def test_past_the_bound_the_note_says_what_is_drawn(self):
        assert pixel_size_note((1.0, 5000.0)) == (
            "1 x 5000 um, drawn 1 to 20: not to scale"
        )
        assert pixel_size_note((5000.0, 1.0)) == (
            "5000 x 1 um, drawn 20 to 1: not to scale"
        )
        assert pixel_size_note((1.0, 20.0)) == "1 x 20 um, drawn to scale"

    def test_a_ratio_that_prints_as_the_bound_is_not_called_out_of_scale(self):
        assert pixel_size_note((1.0, 20.00001)) == "1 x 20 um, drawn to scale"
        assert pixel_size_note((1.0, 20.0003)) == (
            "1 x 20.0003 um, drawn 1 to 20: not to scale"
        )

    def test_the_note_of_a_figure_has_no_letter_the_built_in_font_lacks(self):
        """Pillow's built-in font has no micro sign; it would draw a box."""
        note = pixel_size_note((70.0, 150.0))

        assert note is not None and note.isascii()


class TestARasterOfPixelsThatAreNotSquare:
    @pytest.mark.parametrize(
        "pixel_size", [(70.0, 150.0), (150.0, 70.0), (10.0, 33.0), (1.0, 20.0)]
    )
    def test_each_pixel_is_one_rectangle_of_its_own_colour(self, pixel_size):
        """Whole pixels, no colour that was not in the raster."""
        pixels = np.random.default_rng(7).integers(0, 256, (7, 9, 3), dtype=np.uint8)
        wide, tall = pixel_stretch(pixel_size)

        drawn = np.asarray(raster_image(pixels, 480, stretch=(wide, tall)))

        heights, widths = _blocks(drawn, 0), _blocks(drawn, 1)
        assert (len(heights), len(widths)) == (7, 9)
        rebuilt = np.repeat(np.repeat(pixels, heights, axis=0), widths, axis=1)
        assert (drawn == rebuilt).all()
        # Along the shorter side of a pixel every block is the same whole
        # number; along the longer one, the two whole numbers around it.
        short, long = (widths, heights) if tall > wide else (heights, widths)
        assert len(set(short)) == 1
        exact = short[0] * max(wide, tall)
        assert set(long) <= {int(np.floor(exact)), int(np.ceil(exact))}

    @pytest.mark.parametrize("shape", [(12, 30), (40, 180), (3, 4), (1, 600), (97, 5)])
    @pytest.mark.parametrize("pixel_size", [(70.0, 150.0), (150.0, 70.0), (10.0, 33.0)])
    @pytest.mark.parametrize("reduce", [False, True])
    def test_the_size_is_within_a_pixel_of_the_measured_area(
        self, shape, pixel_size, reduce
    ):
        height, width = shape
        area = (width * pixel_size[0]) / (height * pixel_size[1])

        drawn_width, drawn_height = stretched_size(
            height, width, 480, pixel_stretch(pixel_size), reduce
        )

        if area >= 1:
            assert abs(drawn_width / area - drawn_height) <= 1
        else:
            assert abs(drawn_height * area - drawn_width) <= 1
        assert max(drawn_width, drawn_height) >= 480

    def test_a_square_stretch_is_the_drawing_it_always_was(self):
        pixels = np.random.default_rng(7).integers(0, 256, (36, 48, 3), dtype=np.uint8)

        for long_side in (20, 480):
            for reduce in (False, True):
                plain = raster_image(pixels, long_side, reduce=reduce)
                square = raster_image(pixels, long_side, reduce=reduce, stretch=SQUARE)
                assert plain.size == square.size
                assert (np.asarray(plain) == np.asarray(square)).all()

    def test_a_figure_keeps_every_pixel_of_a_large_raster(self):
        pixels = np.zeros((1000, 500, 3), dtype=np.uint8)

        image = raster_image(pixels, 480, stretch=pixel_stretch((70.0, 150.0)))

        assert image.size == (500, 2143)

    def test_the_report_reduces_a_large_raster_to_the_long_side(self):
        pixels = np.zeros((1000, 500, 3), dtype=np.uint8)
        stretch = pixel_stretch((70.0, 150.0))

        image = raster_image(pixels, 480, reduce=True, stretch=stretch)

        assert image.size == (112, 480)

    def test_the_report_does_not_average_a_raster_it_would_not_reduce(self):
        """Long pixels take 200 rows past 320 pixels; the raster is not past it."""
        pixels = np.random.default_rng(7).integers(
            0, 256, (200, 200, 3), dtype=np.uint8
        )
        stretch = pixel_stretch((70.0, 150.0))

        drawn = np.asarray(raster_image(pixels, 320, reduce=True, stretch=stretch))

        heights, widths = _blocks(drawn, 0), _blocks(drawn, 1)
        assert drawn.shape[:2] == (429, 200)
        assert set(widths) == {1} and set(heights) <= {2, 3}
        rebuilt = np.repeat(np.repeat(pixels, heights, axis=0), widths, axis=1)
        assert (drawn == rebuilt).all()

    def test_a_raster_is_reduced_when_a_square_one_of_its_shape_would_be(self):
        stretch = pixel_stretch((70.0, 150.0))

        assert stretched_size(320, 320, 320, stretch, reduce=True) == (320, 686)
        assert stretched_size(321, 321, 320, stretch, reduce=True) == (149, 320)

    def test_the_largest_image_of_a_report_is_twenty_times_its_long_side(self):
        stretch = pixel_stretch((1.0, 1e9))

        assert stretched_size(480, 480, 480, stretch, reduce=True) == (480, 9600)
        assert stretched_size(481, 480, 480, stretch, reduce=True) == (24, 480)

    def test_a_reduced_raster_never_loses_a_side(self):
        pixels = np.zeros((1, 5000, 3), dtype=np.uint8)
        stretch = pixel_stretch((5000.0, 1.0))

        assert raster_image(pixels, 480, reduce=True, stretch=stretch).size == (480, 1)

    def test_float_noise_does_not_add_one_to_a_whole_factor(self):
        """Five rows of 7 by 8 um reach 480 pixels at 84 times, not 85."""
        assert stretched_size(5, 3, 480, pixel_stretch((7.0, 8.0))) == (3 * 84, 480)

    def test_a_raster_with_no_pixels_is_left_alone(self):
        empty = np.zeros((0, 0, 3), dtype=np.uint8)

        assert raster_image(empty, 480, stretch=(1.0, 2.0)).size == (0, 0)


class TestTheFigure:
    def test_pixels_that_are_not_square_get_a_line_that_states_both_sides(self):
        ion = _ion(np.ones((3, 4)))

        assert len(figure_lines(ion, "sample", 1.0)) == 2
        assert len(figure_lines(ion, "sample", 1.0, (20.0, 20.0))) == 2
        assert figure_lines(ion, "sample", 1.0, (70.0, 150.0))[2] == (
            "Pixel size: 70 x 150 um, drawn to scale."
        )

    def test_the_stretched_raster_is_in_it_pixel_for_pixel(self):
        image = np.arange(12.0).reshape(3, 4)
        stretch = pixel_stretch((70.0, 150.0))
        figure = np.asarray(ion_image_figure(_ion(image), "sample", (70.0, 150.0)))
        coloured = colourise(image, scale_top(image))
        raster = np.asarray(raster_image(coloured, 480, stretch=stretch))

        height, width = raster.shape[:2]
        assert (width, height) == (4 * 75, 482)
        found = any(
            (figure[top : top + height, 16 : 16 + width] == raster).all()
            for top in range(figure.shape[0] - height + 1)
        )
        assert found

    def test_square_pixels_give_the_figure_they_always_did(self):
        ion = _ion(np.arange(12.0).reshape(3, 4))

        plain = png_bytes(ion_image_figure(ion, "sample"))

        assert png_bytes(ion_image_figure(ion, "sample", (20.0, 20.0))) == plain
        assert png_bytes(ion_image_figure(ion, "sample", None)) == plain

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
