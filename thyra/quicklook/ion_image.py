"""From an m/z value to a raster: the window of bins, then the pixels.

An ion image is the per-pixel sum of the m/z bins inside a window around
the value asked for. The window is stated in parts per million of that
value, because that is how a mass tolerance is stated in the field, and a
bin belongs to it when its m/z lies inside, edges included.

The raster covers the bounding box of the measured pixels. A position the
instrument never measured has no row in the table, and is ``NaN`` here so
that it can be drawn apart from a measured pixel with no signal.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from numpy.typing import NDArray

from .store import MsiTable, Positions, QuicklookError

#: Half-width of the default window, in parts per million of the m/z value.
DEFAULT_PPM = 10.0


@dataclass(frozen=True)
class MzWindow:
    """The m/z bins one ion image adds up.

    Attributes:
        target: The m/z value asked for.
        ppm: The half-width asked for, in parts per million of ``target``.
        first: The first column inside the window.
        stop: One past the last column inside it.
        low: The m/z of the first column.
        high: The m/z of the last column.
    """

    target: float
    ppm: float
    first: int
    stop: int
    low: float
    high: float

    @property
    def n_bins(self) -> int:
        """How many bins the window holds."""
        return self.stop - self.first


def _ppm_that_reaches(target: float, value: float) -> int:
    """The smallest whole ppm half-width around ``target`` that holds ``value``."""
    needed = max(1, math.ceil(abs(value - target) / target * 1e6))
    while not _edges(target, needed)[0] <= value <= _edges(target, needed)[1]:
        needed += 1
    return needed


def _edges(target: float, ppm: float) -> tuple[float, float]:
    """The window's lower and upper m/z."""
    tolerance = target * ppm / 1e6
    return target - tolerance, target + tolerance


def select_bins(mz: NDArray[np.float64], target: float, ppm: float) -> MzWindow:
    """The columns whose m/z lies within ``ppm`` of ``target``.

    Args:
        mz: The mass axis, strictly increasing.
        target: The m/z value asked for.
        ppm: Half-width of the window, in parts per million of ``target``.

    Returns:
        The window. It holds at least one bin.

    Raises:
        QuicklookError: If the window lies outside the mass axis, or inside
            it with no bin in it. The second names the nearest bin and the
            ``--ppm`` that would reach it.
    """
    if mz.size == 0:
        raise QuicklookError("The table has no m/z bins.")
    low, high = _edges(target, ppm)
    if high < mz[0] or low > mz[-1]:
        raise QuicklookError(
            f"m/z {target:.4f} is outside the mass axis of this store, which "
            f"runs from {mz[0]:.4f} to {mz[-1]:.4f}."
        )
    first = int(np.searchsorted(mz, low, side="left"))
    stop = int(np.searchsorted(mz, high, side="right"))
    if first == stop:
        neighbours = mz[max(first - 1, 0) : first + 1]
        nearest = float(neighbours[np.argmin(np.abs(neighbours - target))])
        raise QuicklookError(
            f"No m/z bin lies within {ppm:g} ppm of {target:.4f}. The nearest "
            f"bin is at {nearest:.4f}; --ppm {_ppm_that_reaches(target, nearest)} "
            "would reach it."
        )
    return MzWindow(target, ppm, first, stop, float(mz[first]), float(mz[stop - 1]))


def _plane(
    x: NDArray[np.int64],
    y: NDArray[np.int64],
    values: NDArray[np.float64],
    shape: tuple[int, int],
) -> NDArray[np.float64]:
    image = np.full(shape, np.nan, dtype=np.float64)
    image[y, x] = values
    return image


def place_pixels(
    positions: Positions, values: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Lay one value per row out on the raster.

    The raster is the bounding box of the measured pixels, so a region
    whose positions start far from zero is not drawn in a corner of an
    empty frame. Positions with no row stay ``NaN``.

    A volume is drawn as its slices side by side, in order, left to right
    and then top to bottom, one ``NaN`` column or row apart.

    Args:
        positions: Where each row sits.
        values: One value per row.

    Returns:
        The raster, ``NaN`` where nothing was measured.

    Raises:
        QuicklookError: If there are no rows to place.
    """
    if values.size == 0:
        raise QuicklookError("The table holds no pixels, so there is nothing to draw.")
    x = positions.x - positions.x.min()
    y = positions.y - positions.y.min()
    shape = (int(y.max()) + 1, int(x.max()) + 1)
    if positions.z is None:
        return _plane(x, y, values, shape)

    slices: List[NDArray[np.float64]] = []
    for z in np.unique(positions.z):
        rows = positions.z == z
        slices.append(_plane(x[rows], y[rows], values[rows], shape))
    columns = math.ceil(math.sqrt(len(slices)))
    lines = math.ceil(len(slices) / columns)
    height, width = shape
    sheet = np.full(
        (lines * (height + 1) - 1, columns * (width + 1) - 1), np.nan, dtype=np.float64
    )
    for index, tile in enumerate(slices):
        top = (index // columns) * (height + 1)
        left = (index % columns) * (width + 1)
        sheet[top : top + height, left : left + width] = tile
    return sheet


@dataclass(frozen=True)
class IonImage:
    """One ion image and what it was made from.

    Attributes:
        table: The table the image was read from.
        window: The m/z bins that were added up.
        image: The raster, ``NaN`` where nothing was measured.
        lowest: The smallest per-pixel sum.
        highest: The largest per-pixel sum.
        n_pixels: The number of measured pixels.
    """

    table: str
    window: MzWindow
    image: NDArray[np.float64]
    lowest: float
    highest: float
    n_pixels: int


def make_ion_image(
    table: MsiTable,
    target: float,
    ppm: float = DEFAULT_PPM,
    positions: Optional[Positions] = None,
) -> IonImage:
    """Read one ion image from a table of spectra.

    Args:
        table: The table to read.
        target: The m/z value.
        ppm: Half-width of the window, in parts per million of ``target``.
        positions: The table's pixel positions, when the caller has read
            them already; read here otherwise.

    Returns:
        The image, with the window it covers.

    Raises:
        QuicklookError: If the window holds no bin, or the table no pixel.
    """
    window = select_bins(table.mz, target, ppm)
    values = table.window_sum(window.first, window.stop)
    image = place_pixels(positions or table.positions(), values)
    return IonImage(
        table=table.name,
        window=window,
        image=image,
        lowest=float(values.min()),
        highest=float(values.max()),
        n_pixels=int(values.size),
    )


def tallest_peak_near(
    table: MsiTable, target: float, reach: float = 1.0
) -> Optional[float]:
    """The m/z of the tallest bin of the mean spectrum within ``reach``.

    Used to point somewhere useful when an ion image comes out blank: a
    value typed to one decimal often falls between two peaks.

    Args:
        table: The table the image was read from.
        target: The m/z value asked for.
        reach: How far from ``target`` to look, in m/z.

    Returns:
        The m/z, or ``None`` when the table has no mean spectrum or the
        mean spectrum has no signal that close.
    """
    spectrum = table.average_spectrum()
    if spectrum is None:
        return None
    first = int(np.searchsorted(table.mz, target - reach, side="left"))
    stop = int(np.searchsorted(table.mz, target + reach, side="right"))
    if first == stop or float(spectrum[first:stop].max()) <= 0:
        return None
    return float(table.mz[first + int(np.argmax(spectrum[first:stop]))])
