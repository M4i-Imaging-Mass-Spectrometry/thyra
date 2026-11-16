"""Turn a raster of numbers into pixels, with numpy and Pillow only.

matplotlib is not a dependency of Thyra and this module is the reason it
does not have to become one: a colour map is a lookup table, a colour bar
is that table drawn as a strip, and the labels are a few lines of text.

Three choices are made here, once, for both commands:

- The colour map is viridis, which is perceptually uniform and readable
  with the common colour-vision deficiencies. It is in the public domain.
- A position that was never measured is drawn in a neutral grey, which
  viridis does not contain, so it cannot be read as a low or high value.
- The colour scale runs from zero to the 99th percentile of the pixels
  that carry signal. A handful of very bright pixels would otherwise set
  the scale and leave the rest of the image dark. The pixels above that
  value all take the top colour, and the highest value is printed.
"""

from __future__ import annotations

import io
from typing import TYPE_CHECKING, Any, Tuple

import numpy as np
from numpy.typing import NDArray

from .ion_image import IonImage

if TYPE_CHECKING:
    from PIL import Image, ImageDraw, ImageFont

#: Viridis at 33 evenly spaced positions, as ``rrggbb``. The 256 colours are
#: interpolated between them; no channel is then more than 2 of 255 away from
#: the published map.
_VIRIDIS_ANCHORS = (
    "440154470d6048186a482374472d7b4537814240863e49893b528b375b8d33638d"
    "2f6b8e2c728e297a8e26828e23898e21918c1f978b1f9f8821a68527ad8131b57b"
    "3dbc744cc26c5cc8636ece5881d34d95d840aadc32c0df25d5e21aeae51afde725"
)

#: The colour of a position that was never measured.
NOT_MEASURED = (190, 190, 190)

#: The share of the pixels with signal that the colour scale spans.
SCALE_PERCENTILE = 99.0

#: A raster is enlarged until its longer side is at least this many pixels.
FIGURE_LONG_SIDE = 480

_MARGIN = 16
_BAR_WIDTH = 18
_GAP = 10
_INK = (20, 20, 20)
_PAPER = (255, 255, 255)


def _build_lut() -> NDArray[np.uint8]:
    """The 256 viridis colours, ``(256, 3)``."""
    raw = bytes.fromhex(_VIRIDIS_ANCHORS)
    anchors = np.frombuffer(raw, dtype=np.uint8).reshape(-1, 3).astype(np.float64)
    at = np.linspace(0.0, 255.0, len(anchors))
    channels = [np.interp(np.arange(256), at, anchors[:, c]) for c in range(3)]
    return np.rint(np.stack(channels, axis=1)).astype(np.uint8)


#: Colour of each of the 256 levels, lowest first.
VIRIDIS: NDArray[np.uint8] = _build_lut()


def scale_top(image: NDArray[np.float64]) -> float:
    """The value the top colour stands for.

    Args:
        image: A raster, ``NaN`` where nothing was measured.

    Returns:
        The 99th percentile of the values above zero, or 0.0 when no pixel
        carries signal.
    """
    signal = image[np.isfinite(image) & (image > 0)]
    if signal.size == 0:
        return 0.0
    return float(np.percentile(signal, SCALE_PERCENTILE))


def colourise(image: NDArray[np.float64], top: float) -> NDArray[np.uint8]:
    """Colour a raster: viridis from zero to ``top``, grey where not measured.

    Args:
        image: A raster, ``NaN`` where nothing was measured.
        top: The value that takes the top colour; larger values take it too.

    Returns:
        ``(height, width, 3)`` RGB.
    """
    measured = np.isfinite(image)
    scaled = np.zeros(image.shape, dtype=np.float64)
    if top > 0:
        np.divide(image, top, out=scaled, where=measured)
    levels = np.rint(np.clip(scaled, 0.0, 1.0) * 255).astype(np.intp)
    pixels = VIRIDIS[levels]
    pixels[~measured] = NOT_MEASURED
    return pixels


def integer_scale(height: int, width: int, long_side: int) -> int:
    """The whole factor that brings the longer side to at least ``long_side``."""
    return max(1, -(-long_side // max(height, width, 1)))


def raster_image(
    pixels: NDArray[np.uint8], long_side: int, reduce: bool = False
) -> "Image.Image":
    """An RGB raster as an image whose longer side is about ``long_side``.

    A small raster is enlarged by a whole factor, each pixel becoming a
    square block, so no value is blurred into its neighbour.

    Args:
        pixels: ``(height, width, 3)`` RGB, or ``(height, width)`` grey.
        long_side: The size to aim for.
        reduce: Whether a larger raster is reduced to ``long_side`` by
            averaging. The report does, to stay small; a figure keeps
            every pixel.

    Returns:
        The image.
    """
    from PIL import Image

    image = Image.fromarray(pixels)
    height, width = pixels.shape[:2]
    factor = integer_scale(height, width, long_side)
    if factor > 1:
        return image.resize((width * factor, height * factor), Image.Resampling.NEAREST)
    if reduce and max(height, width) > long_side:
        image.thumbnail((long_side, long_side), Image.Resampling.BOX)
    return image


def colour_bar(height: int, width: int = _BAR_WIDTH) -> "Image.Image":
    """The colour map as a vertical strip, the top colour at the top."""
    from PIL import Image

    levels = np.rint(np.linspace(255, 0, height)).astype(np.intp)
    strip = np.repeat(VIRIDIS[levels][:, np.newaxis, :], width, axis=1)
    return Image.fromarray(strip)


def format_intensity(value: float) -> str:
    """An intensity in three significant digits, as ``5380`` or ``1.2e+07``."""
    return format(float(f"{value:.3g}"), "g")


def _font(size: int) -> Any:
    """Pillow's built-in font at ``size``, or its fixed-size one."""
    from PIL import ImageFont

    try:
        return ImageFont.load_default(size=size)
    except (ImportError, OSError):  # a Pillow built without FreeType
        return ImageFont.load_default()


def _text_size(
    draw: "ImageDraw.ImageDraw",
    text: str,
    font: "ImageFont.FreeTypeFont | ImageFont.ImageFont",
) -> Tuple[int, int]:
    left, top, right, bottom = draw.textbbox((0, 0), text, font=font)
    return int(right - left), int(bottom)


def ion_image_figure(ion: IonImage, source: str) -> "Image.Image":
    """One ion image as a figure: title, raster, colour bar and legend.

    Args:
        ion: The ion image.
        source: What the image was made from, printed under the title.

    Returns:
        An RGB image on a white ground.
    """
    from PIL import Image, ImageDraw

    top = scale_top(ion.image)
    raster = raster_image(colourise(ion.image, top), FIGURE_LONG_SIDE)
    window = ion.window
    title = f"m/z {window.target:.4f} +/- {window.ppm:g} ppm"
    lines = [
        f"{source}: {window.n_bins} m/z bins, "
        f"{window.low:.4f} to {window.high:.4f}",
        f"Colour: 0 to {format_intensity(top)}. "
        f"Highest pixel: {format_intensity(ion.highest)}.",
    ]
    unmeasured = bool(np.isnan(ion.image).any())
    legend = "not measured"
    labels = (format_intensity(top), "0")

    title_font, text_font = _font(20), _font(13)
    probe = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    title_width, title_height = _text_size(probe, title, title_font)
    line_height = _text_size(probe, "Ag", text_font)[1] + 4
    label_width = max(_text_size(probe, label, text_font)[0] for label in labels)
    legend_width = line_height + 2 + _text_size(probe, legend, text_font)[0]

    body_width = raster.width + _GAP + _BAR_WIDTH + 6 + label_width
    text_width = max(
        [title_width, legend_width]
        + [_text_size(probe, line, text_font)[0] for line in lines]
    )
    width = 2 * _MARGIN + max(body_width, text_width)
    header = _MARGIN + title_height + 6 + line_height * len(lines) + 6
    footer = (line_height + 6 if unmeasured else 0) + _MARGIN
    # A line scan is a few pixels high: the colour bar keeps room for its
    # two labels however flat the raster is.
    body_height = max(raster.height, 4 * line_height)
    canvas = Image.new("RGB", (width, header + body_height + footer), _PAPER)
    draw = ImageDraw.Draw(canvas)

    draw.text((_MARGIN, _MARGIN), title, fill=_INK, font=title_font)
    for index, line in enumerate(lines):
        y = _MARGIN + title_height + 6 + index * line_height
        draw.text((_MARGIN, y), line, fill=_INK, font=text_font)

    canvas.paste(raster, (_MARGIN, header))
    bar_left = _MARGIN + raster.width + _GAP
    canvas.paste(colour_bar(body_height), (bar_left, header))
    label_left = bar_left + _BAR_WIDTH + 6
    low_y = header + body_height - line_height
    draw.text((label_left, header), labels[0], fill=_INK, font=text_font)
    draw.text((label_left, low_y), labels[1], fill=_INK, font=text_font)
    if unmeasured:
        y = header + body_height + 6
        swatch = (_MARGIN, y, _MARGIN + line_height - 4, y + line_height - 4)
        draw.rectangle(swatch, fill=NOT_MEASURED)
        draw.text((_MARGIN + line_height + 2, y), legend, fill=_INK, font=text_font)
    return canvas


def png_bytes(image: "Image.Image") -> bytes:
    """An image encoded as PNG."""
    buffer = io.BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    return buffer.getvalue()
