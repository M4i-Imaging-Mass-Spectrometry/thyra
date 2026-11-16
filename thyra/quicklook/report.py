"""One HTML file that shows what a converted store holds.

The page is for somebody who has just converted their first dataset and
wants to see that it worked. It is one file: the images are embedded, the
spectrum is drawn inline, and nothing is fetched, so it opens without a
network connection and can be sent to a colleague as it is.

**What is printed from the store is an allow-list.** A store keeps the
path it was converted from, and its source metadata can hold a sample
name that identifies a person. The page prints the instrument, the pixel
size, the Thyra version and the *name* of the source file, and reads no
other metadata field. A path is cut down to its last part whichever
separator it was written with.
"""

from __future__ import annotations

import base64
import html
import math
from pathlib import PureWindowsPath
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from .ion_image import DEFAULT_PPM, make_ion_image, place_pixels
from .render import (
    NOT_MEASURED,
    VIRIDIS,
    colourise,
    format_intensity,
    png_bytes,
    raster_image,
    scale_top,
)
from .store import (
    MOBILITY,
    MSI,
    MSMS,
    MsiTable,
    Positions,
    QuicklookError,
    Store,
    TableSummary,
)

if TYPE_CHECKING:
    from PIL import Image

#: How many ion images the page shows for each table of spectra.
PEAK_COUNT = 6

#: Two peaks shown side by side are at least this far apart in m/z, so an
#: isotope peak is not shown beside the peak it belongs to.
PEAK_SEPARATION = 1.5

_TIC_LONG_SIDE = 480
_ION_LONG_SIDE = 320
_OPTICAL_LONG_SIDE = 640
_SPECTRUM_POINTS = 850

_KIND_LABELS = {MSI: "spectra", MOBILITY: "ion mobility", MSMS: "MS/MS"}

_STYLE = """
:root { color-scheme: light dark; }
body { font-family: system-ui, sans-serif; line-height: 1.45;
  max-width: 62rem; margin: 2rem auto; padding: 0 1rem; }
h1 { margin-bottom: 0.2rem; }
h2, h3 { margin-top: 2.2rem; }
table { border-collapse: collapse; }
th, td { text-align: left; padding: 0.3rem 1.5rem 0.3rem 0;
  border-bottom: 1px solid rgba(128, 128, 128, 0.4); }
figure { margin: 0; }
figcaption, .note { font-size: 0.9rem; opacity: 0.85; }
img { max-width: 100%; height: auto; image-rendering: pixelated; }
img.photo { image-rendering: auto; }
svg { width: 100%; height: auto; }
.grid { display: grid; gap: 1.2rem;
  grid-template-columns: repeat(auto-fill, minmax(15rem, 1fr)); }
.bar, .swatch { display: inline-block; height: 0.8rem; vertical-align: middle;
  border: 1px solid rgba(128, 128, 128, 0.6); }
.bar { width: 8rem; background: linear-gradient(to right, BAR_STOPS); }
.swatch { width: 0.8rem; background: SWATCH; }
code { font-size: 0.95em; }
"""


def _css() -> str:
    stops = ", ".join("#%02x%02x%02x" % tuple(VIRIDIS[i]) for i in range(0, 256, 15))
    swatch = "#%02x%02x%02x" % NOT_MEASURED
    return _STYLE.replace("BAR_STOPS", stops).replace("SWATCH", swatch)


def file_name_only(path: Any) -> Optional[str]:
    """The last part of a stored path, whichever separator it uses.

    Args:
        path: A path as a store records it, or anything else.

    Returns:
        The file name, or ``None`` when ``path`` is not a non-empty string.
    """
    if not isinstance(path, str) or not path.strip():
        return None
    return PureWindowsPath(path.strip().rstrip("\\/")).name or None


def pick_peaks(
    mz: NDArray[np.float64],
    spectrum: NDArray[np.float64],
    count: int = PEAK_COUNT,
    separation: float = PEAK_SEPARATION,
) -> List[int]:
    """The tallest peaks of a spectrum that stand apart from each other.

    A peak is a bin taller than the bin before it and at least as tall as
    the one after. They are taken tallest first, skipping any closer than
    ``separation`` to one already taken.

    Args:
        mz: The mass axis.
        spectrum: One intensity per bin.
        count: How many peaks to return at most.
        separation: The smallest m/z distance between two of them.

    Returns:
        The column of each peak, tallest first.
    """
    if spectrum.size == 0:
        return []
    before = np.concatenate(([-np.inf], spectrum[:-1]))
    after = np.concatenate((spectrum[1:], [-np.inf]))
    peaks = np.flatnonzero((spectrum > 0) & (spectrum > before) & (spectrum >= after))
    # Enough candidates that the separation rule cannot run out of them on
    # any real spectrum, few enough that the loop below stays short.
    tallest = peaks[np.argsort(spectrum[peaks], kind="stable")[::-1][:5000]]
    chosen: List[int] = []
    for column in tallest:
        if all(abs(mz[column] - mz[other]) >= separation for other in chosen):
            chosen.append(int(column))
            if len(chosen) == count:
                break
    return chosen


def _data_uri(image: "Image.Image") -> str:
    return "data:image/png;base64," + base64.b64encode(png_bytes(image)).decode("ascii")


def _figure(image: "Image.Image", caption: str, css_class: str = "") -> str:
    attribute = f' class="{css_class}"' if css_class else ""
    return (
        f'<figure><img{attribute} src="{_data_uri(image)}" '
        f'width="{image.width}" height="{image.height}" alt="{html.escape(caption)}">'
        f"<figcaption>{html.escape(caption)}</figcaption></figure>"
    )


def _nice_ticks(low: float, high: float, wanted: int = 6) -> List[float]:
    """Round tick positions between ``low`` and ``high``."""
    span = high - low
    if span <= 0:
        return [low]
    raw = span / wanted
    magnitude = 10 ** math.floor(math.log10(raw))
    step = next(m * magnitude for m in (1, 2, 5, 10) if m * magnitude >= raw)
    first = math.ceil(low / step) * step
    return [first + i * step for i in range(int((high - first) / step + 1e-9) + 1)]


def spectrum_svg(mz: NDArray[np.float64], spectrum: NDArray[np.float64]) -> str:
    """A spectrum as an inline SVG, reduced to what a screen can show.

    The mass range is cut into at most 850 equal steps and each step is
    drawn at the tallest value inside it, so a narrow peak keeps its
    height however many bins the axis has.

    Args:
        mz: The mass axis, increasing.
        spectrum: One intensity per bin.

    Returns:
        The ``<svg>`` element. It takes its colour from the text around it.
    """
    width, height = 900, 300
    left, right, top, bottom = 70, 15, 15, 40
    plot_width, plot_height = width - left - right, height - top - bottom
    low, high = float(mz[0]), float(mz[-1])
    span = (high - low) or 1.0
    peak = float(spectrum.max()) if spectrum.size else 0.0

    steps = min(_SPECTRUM_POINTS, spectrum.size)
    edges = np.searchsorted(mz, low + span * np.arange(steps) / steps)
    filled = np.diff(np.append(edges, spectrum.size)) > 0
    tallest = np.maximum.reduceat(spectrum, np.minimum(edges, spectrum.size - 1))
    tallest = np.where(filled, tallest, 0.0)

    xs = left + plot_width * (np.arange(steps) + 0.5) / steps
    ys = top + plot_height * (1.0 - tallest / (peak or 1.0))
    floor = top + plot_height
    points = " ".join(f"L{x:.1f},{y:.1f}" for x, y in zip(xs, ys))
    path = f"M{left},{floor} {points} L{left + plot_width},{floor} Z"

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        'aria-label="Mean spectrum" fill="currentColor" font-size="13">',
        f'<path d="{path}" stroke="currentColor" stroke-width="1" '
        'fill-opacity="0.2" stroke-linejoin="round"/>',
        f'<path d="M{left},{top} V{floor} H{left + plot_width}" fill="none" '
        'stroke="currentColor" stroke-opacity="0.6"/>',
        f'<text x="{left - 8}" y="{top + 5}" text-anchor="end">'
        f"{format_intensity(peak)}</text>",
        f'<text x="{left - 8}" y="{floor + 4}" text-anchor="end">0</text>',
    ]
    for tick in _nice_ticks(low, high):
        x = left + plot_width * (tick - low) / span
        parts.append(
            f'<path d="M{x:.1f},{floor} v5" stroke="currentColor" '
            'stroke-opacity="0.6"/>'
            f'<text x="{x:.1f}" y="{floor + 20}" text-anchor="middle">{tick:g}</text>'
        )
    parts.append(
        f'<text x="{left + plot_width / 2:.0f}" y="{height - 4}" '
        'text-anchor="middle">m/z</text></svg>'
    )
    return "".join(parts)


def _text(value: Any) -> Optional[str]:
    return value.strip() if isinstance(value, str) and value.strip() else None


def _pixel_size(block: Mapping[str, Any]) -> str:
    analysis = block.get("ms_analysis") or {}
    size = analysis.get("pixel_size_um") or {}
    try:
        x, y = float(size["x"]), float(size["y"])
    except (KeyError, TypeError, ValueError):
        return "not stated"
    stated = f"{x:g} µm" if x == y else f"{x:g} x {y:g} µm"
    if (block.get("provenance") or {}).get("pixel_size_source") == "default":
        stated += " (assumed: the source states none)"
    return stated


def table_facts(
    table: MsiTable, positions: Positions, block: Mapping[str, Any]
) -> List[Tuple[str, str]]:
    """The short table of facts about one table of spectra.

    Args:
        table: The table.
        positions: Its pixel positions.
        block: Its ``msi_metadata`` block, or an empty mapping.

    Returns:
        ``(label, value)`` pairs, in the order they are shown.
    """
    analysis = block.get("ms_analysis") or {}
    provenance = block.get("provenance") or {}
    maker, model = _text(analysis.get("manufacturer")), _text(
        analysis.get("instrument_model")
    )
    if maker and model and maker.lower() in model.lower():
        maker = None
    instrument = " ".join(part for part in (maker, model) if part)
    raster = (
        f"{int(np.ptp(positions.x)) + 1} x {int(np.ptp(positions.y)) + 1} pixels"
        if table.n_rows
        else "no pixels"
    )
    if positions.z is not None and table.n_rows:
        raster += f", {np.unique(positions.z).size} slices"
    mz_range = (
        f"{table.mz[0]:.4f} to {table.mz[-1]:.4f}" if table.mz.size else "no bins"
    )
    version = _text(provenance.get("thyra_version"))
    return [
        ("Instrument", instrument or "not stated in the source"),
        ("Pixel size", _pixel_size(block)),
        ("Raster", raster),
        ("Measured pixels", f"{table.n_rows:,}"),
        ("m/z range", mz_range),
        ("m/z bins", f"{table.n_columns:,}"),
        ("Converted with", f"Thyra {version}" if version else "not recorded"),
        (
            "Source file",
            file_name_only(provenance.get("source_path")) or "not recorded",
        ),
    ]


def _rows(rows: Sequence[Sequence[str]], header: Sequence[str] = ()) -> str:
    """An HTML table; the first cell of each row is its heading."""
    out = ["<table>"]
    if header:
        cells = "".join(f"<th>{html.escape(cell)}</th>" for cell in header)
        out.append(f"<tr>{cells}</tr>")
    for row in rows:
        cells = f"<th>{html.escape(row[0])}</th>" + "".join(
            f"<td>{html.escape(cell)}</td>" for cell in row[1:]
        )
        out.append(f"<tr>{cells}</tr>")
    out.append("</table>")
    return "".join(out)


def _tic_raster(
    tic: NDArray[np.float64], positions: Positions
) -> Optional[NDArray[np.float64]]:
    """The TIC image on the same raster as the ion images.

    The stored image is read at each measured pixel and laid out again, so
    that a position without a measurement is told apart from a measured
    pixel with no signal. An image the positions do not index is shown as
    stored, when it is a single plane.

    A volume of one plane is the case to mind: its table carries ``z``, all
    of one value, and its TIC image is still a single plane.
    """
    stored = tic if tic.ndim == 2 and tic.size else None
    index: Tuple[NDArray[np.int64], ...] = (positions.y, positions.x)
    if tic.ndim == 3 and positions.z is not None:
        index = (positions.z,) + index
    elif positions.z is not None and np.unique(positions.z).size > 1:
        return stored
    inside = tic.ndim == len(index) and all(
        axis.size and int(axis.min()) >= 0 and int(axis.max()) < size
        for axis, size in zip(index, tic.shape)
    )
    return place_pixels(positions, tic[index]) if inside else stored


def _command_example(store: Store, mz: float) -> str:
    """The ``thyra ion-image`` command for this store, ready to paste."""
    name = store.path.name
    if any(character.isspace() for character in name):
        name = f'"{name}"'
    return f"thyra ion-image {name} {mz:.4f}"


def _ion_image_section(
    store: Store,
    table: MsiTable,
    positions: Positions,
    spectrum: NDArray[np.float64],
    heading: str,
) -> List[str]:
    """Ion images of the tallest peaks of one table's mean spectrum.

    A table whose matrix cannot be read a block of columns at a time, one
    written by an old Thyra version, gets the reason in place of the
    images: the rest of the page does not depend on the matrix.
    """
    peaks = pick_peaks(table.mz, spectrum)
    if not peaks:
        return []
    out = [f"<{heading}>Ion images of the tallest peaks</{heading}>"]
    figures = []
    try:
        for column in peaks:
            target = float(table.mz[column])
            ion = make_ion_image(table, target, DEFAULT_PPM, positions)
            top = scale_top(ion.image)
            pixels = colourise(ion.image, top)
            caption = f"m/z {target:.4f}, 0 to {format_intensity(top)}"
            figures.append(
                _figure(raster_image(pixels, _ION_LONG_SIDE, reduce=True), caption)
            )
    except QuicklookError as exc:
        return out + [f'<p class="note">{html.escape(str(exc))}</p>']
    example = _command_example(store, float(table.mz[peaks[0]]))
    return out + [
        f'<p class="note">The {len(figures)} tallest peaks of the mean '
        f"spectrum that are at least {PEAK_SEPARATION:g} m/z apart, each "
        f"within {DEFAULT_PPM:g} ppm. To draw one as a file, with any m/z value: "
        f"<code>{html.escape(example)}</code></p>",
        f'<div class="grid">{"".join(figures)}</div>',
    ]


def _spectra_sections(store: Store, summary: TableSummary, level: int) -> List[str]:
    """Facts, TIC image, mean spectrum and ion images of one table."""
    table = store.open_table(summary.name)
    positions = table.positions()
    heading = f"h{level}"
    facts = table_facts(table, positions, table.metadata_block())
    out = [_rows([list(pair) for pair in facts])]

    tic = store.tic_image(summary.name)
    raster = _tic_raster(tic, positions) if tic is not None else None
    if raster is not None:
        top = scale_top(raster)
        image = raster_image(colourise(raster, top), _TIC_LONG_SIDE, reduce=True)
        caption = f"Total ion current of every pixel, 0 to {format_intensity(top)}."
        out += [f"<{heading}>TIC image</{heading}>", _figure(image, caption)]

    spectrum = table.average_spectrum()
    if spectrum is None or not table.mz.size:
        return out
    out += [
        f"<{heading}>Mean spectrum</{heading}>",
        spectrum_svg(table.mz, spectrum),
    ]
    return out + _ion_image_section(store, table, positions, spectrum, heading)


def _as_photo(pixels: NDArray[Any]) -> NDArray[np.uint8]:
    """A ``(channels, y, x)`` image as 8-bit grey or RGB, ``(y, x[, 3])``.

    Three or more channels are read as red, green and blue; fewer, as the
    first channel in grey. Anything but 8-bit values is stretched from its
    lowest value, or zero when that is lower, to its highest.
    """
    if pixels.dtype != np.uint8:
        values = np.nan_to_num(pixels.astype(np.float64), posinf=0.0, neginf=0.0)
        low = min(float(values.min()), 0.0) if values.size else 0.0
        span = (float(values.max()) - low) if values.size else 0.0
        scaled = (values - low) * (255.0 / span if span > 0 else 0.0)
        pixels = np.rint(scaled).astype(np.uint8)
    if pixels.shape[0] >= 3:
        return np.ascontiguousarray(np.moveaxis(pixels[:3], 0, -1))
    return np.ascontiguousarray(pixels[0])


def _optical_section(store: Store, level: int) -> List[str]:
    """The store's optical images, reduced; nothing when it has none."""
    tic_names = {f"{table.name}_tic" for table in store.tables}
    recorded: Any = store.root_attrs.get("optical_images")
    elements: Dict[str, Any] = {}
    if isinstance(recorded, dict) and isinstance(recorded.get("elements"), dict):
        elements = recorded["elements"]
    figures = []
    for name in store.image_names():
        if name in tic_names:
            continue
        preview = store.image_preview(name, _OPTICAL_LONG_SIDE)
        if preview is None or preview.size == 0:
            continue
        image = raster_image(_as_photo(preview), 1)
        image.thumbnail((_OPTICAL_LONG_SIDE, _OPTICAL_LONG_SIDE))
        info = elements.get(name)
        source = (
            file_name_only(info.get("source_file")) if isinstance(info, dict) else None
        )
        caption = f"{name}, from {source}" if source else name
        figures.append(_figure(image, caption, "photo"))
    if not figures:
        return []
    title = "Optical image" if len(figures) == 1 else "Optical images"
    return [
        f"<h{level}>{title}</h{level}>",
        f'<div class="grid">{"".join(figures)}</div>',
    ]


def _table_list(store: Store) -> str:
    rows = []
    for table in store.tables:
        rows_, columns = (
            (f"{table.shape[0]:,}", f"{table.shape[1]:,}")
            if table.shape
            else ("?", "?")
        )
        rows.append([table.name, _KIND_LABELS.get(table.kind, "other"), rows_, columns])
    return _rows(rows, ("Table", "Holds", "Pixels", "Columns"))


def build_report(store: Store) -> str:
    """The whole report page for one store.

    Args:
        store: The opened store.

    Returns:
        A complete HTML document that refers to nothing outside itself.

    Raises:
        QuicklookError: If the store holds no table of spectra, or one of
            them lacks a part the page needs.
    """
    tables = store.msi_tables()
    name = html.escape(store.name)
    body = [
        f"<h1>{name}</h1>",
        f'<p class="note">What the store {html.escape(store.path.name)} holds. '
        'Colours run from 0 <span class="bar"></span> to the value under each '
        'image; <span class="swatch"></span> was not measured.</p>',
    ]
    for summary in tables:
        if len(tables) > 1:
            body.append(f"<h2>{html.escape(summary.name)}</h2>")
        level = 3 if len(tables) > 1 else 2
        body += _spectra_sections(store, summary, level)
    body += _optical_section(store, 2)
    body += ["<h2>Tables</h2>", _table_list(store)]

    return (
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="color-scheme" content="light dark">\n'
        f"<title>{name}: Thyra report</title>\n<style>{_css()}</style>\n</head>\n"
        "<body>\n" + "\n".join(body) + "\n</body>\n</html>\n"
    )
