"""List the acquisition regions of an MSI input (issue #448).

A multi-region input converts into one raster that holds every region, and
only the timsTOF and Waters readers can select one region (``region=``). To write one
store per region a caller has to know the regions first. :func:`list_regions`
answers that from the same metadata-only read :func:`thyra.preview_msi`
does: no spectrum is decoded and no store is written.

Unlike :attr:`thyra.MsiPreview.regions`, which is ``None`` for one region so
a wizard can ask "is this one section?", the list here always has at least
one entry. A loop over it then works on every input.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Union

from .errors import ConversionRefused
from .preview import _build_reader, _close_quietly


def _whole_run(reader: Any) -> Dict[str, Any]:
    """One entry for an input with no region concept.

    Numbered 1, which is what ``obs["region_number"]`` holds when a reader
    gives no region map. A reader that states one region lists it itself,
    with the number ``region=`` accepts (D57).
    """
    essential = reader.get_essential_metadata()
    width, height = essential.dimensions[0], essential.dimensions[1]
    return {
        "region_number": 1,
        "name": None,
        "n_spectra": (
            int(essential.n_spectra) if essential.n_spectra_counted else None
        ),
        "bounds": (0, 0, int(width) - 1, int(height) - 1),
    }


def _normalise(entry: Dict[str, Any]) -> Dict[str, Any]:
    """Give a reader's region entry the four keys every entry has."""
    bounds = entry.get("bounds")
    return {
        "region_number": int(entry["region_number"]),
        "name": entry.get("name"),
        "n_spectra": int(entry["n_spectra"]),
        "bounds": tuple(int(v) for v in bounds) if bounds is not None else None,
    }


def list_regions(path: Union[str, Path]) -> List[Dict[str, Any]]:
    """List the acquisition regions of an MSI input without converting it.

    Each entry has four keys:

    * ``region_number``: the value the store writes to
      ``obs["region_number"]`` for that region. On a timsTOF ``.d`` it is
      the DB RegionNumber, which ``region=`` accepts, on a single-region
      run too. On a Waters ``.raw`` it is the region's place in the
      method file's ``MultiRegions`` table, counted from 1.
    * ``name``: the ``.mis`` Area name of the region, or ``None`` when the
      source names none. The whole-run entry of a format with no regions
      has no name.
    * ``n_spectra``: how many spectra the region holds. A position
      measured twice counts twice. ``None`` when the format cannot count
      without reading spectra.
    * ``bounds``: ``(x_min, y_min, x_max, y_max)``, 0-based and inclusive,
      in the pixel grid of the converted store. ``None`` when the reader
      could not read a region's box.

    A format without regions gives one entry, numbered 1, that covers the
    whole raster.

    Args:
        path: The input, as :func:`thyra.convert_msi` takes it.

    Returns:
        One dict per region, sorted by ``region_number``.

    Raises:
        ConversionRefused: When the input cannot be read.
    """
    path = Path(path)
    reader, error = _build_reader(path)
    if reader is None:
        raise ConversionRefused(f"{path.name}: could not be read: {error}")
    try:
        found = reader.get_region_info()
        if not found:
            return [_whole_run(reader)]
        entries = [_normalise(entry) for entry in found]
        return sorted(entries, key=lambda entry: entry["region_number"])
    finally:
        _close_quietly(reader)
