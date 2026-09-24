"""Read a source's metadata document without converting it.

The metadata-only counterpart of :func:`thyra.convert.convert_msi`: it
builds the same versioned ``msi_metadata`` block a conversion would write
into ``table.uns``, straight from the raw file, and writes nothing.

Two kinds of caller want this.  One is anybody deciding what a folder
holds before spending hours converting it -- the same audience as
:func:`thyra.preview_msi`, which answers a shorter question with a
smaller record.  The other is a source Thyra will never convert: a
metadata document describes an acquisition, and plenty of acquisitions
have no image in them.  A Bruker ``.d`` from an electrospray run states
its instrument, polarity, mass range, precursor schedule and acquisition
time exactly as an imaging run does; what it does not state is a pixel
size, because there is no raster.  Such a document is built here and
does not validate against schema 0.6.0, which requires the pitch.  That
is the point rather than a failure: see design decision D23 in
docs/design-decisions.md.

Every reader is built with ``metadata_only=True``, the same way
:mod:`thyra.preview` builds one, so for every format but Waters no
spectra are read and no vendor SDK is loaded. The Waters extractor still
loads MassLynx and passes over every MS scan for the mass range and the
spectrum count.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from ..core.registry import detect_format, get_reader_class
from .schema.builder import build_metadata_document

logger = logging.getLogger(__name__)


def _fragmentation_report(reader: Any) -> Optional[Dict[str, Any]]:
    """What the reader says about fragmentation, in the builder's shape.

    Asked through the base reader contract, so a format that learns to
    report a schedule later needs no change here.  A reader that raises
    costs the document one section rather than the whole document: the
    schedule is the most format-specific thing in the block and the
    least load-bearing.
    """
    try:
        if not reader.has_fragmentation:
            return None
        schedule = reader.get_fragmentation()
    except Exception as e:  # pragma: no cover - reader-defined
        logger.warning("Could not describe the fragmentation: %s", str(e))
        return None
    if schedule is None:
        return None
    report: Dict[str, Any] = schedule.to_extractor_report()
    return report


def _pixel_size(essential: Any) -> Optional[Tuple[float, float]]:
    """The in-plane pitch the source states, or ``None`` when it states none."""
    pitch = getattr(essential, "pixel_size", None)
    if pitch is None:
        return None
    return (float(pitch[0]), float(pitch[1]))


def read_metadata_document(path: Path | str) -> Dict[str, Any]:
    """Build the ``msi_metadata`` document for a raw MSI source.

    Args:
        path: Path to an input in any format Thyra reads -- the same
            paths :func:`thyra.convert.convert_msi` accepts.

    Returns:
        The document as a plain dict, in the shape
        :func:`~thyra.metadata.schema.store_io.read_msi_metadata_blocks`
        returns for a converted store.  ``processing`` is empty: nothing
        has been done to the data.  ``ms_analysis.pixel_size_um`` is
        absent when the source states no raster pitch, which makes the
        document invalid against schema 0.6.0 -- the caller reports that
        rather than this function hiding it.

    Raises:
        ValueError: If the format cannot be detected or has no reader.
            The registry's own message names what Thyra does read.
        Exception: Whatever the reader raises for a file it cannot open.
            Deliberately not swallowed, unlike :func:`thyra.preview_msi`,
            whose contract is to never fail: a preview that cannot read a
            folder still has a card to draw, while a metadata document
            that cannot read one has nothing to write.
    """
    path = Path(path)
    source_format = detect_format(path)
    reader_class = get_reader_class(source_format)

    with reader_class(path, metadata_only=True) as reader:
        comprehensive = reader.get_comprehensive_metadata()
        essential = getattr(comprehensive, "essential", None)
        pitch = _pixel_size(essential)
        if pitch is None:
            logger.info(
                "%s states no pixel size; the document describes an "
                "acquisition rather than an image.",
                path.name,
            )
        return build_metadata_document(
            comprehensive,
            pixel_size_um=pitch,
            pixel_size_source="automatic" if pitch is not None else None,
            source_format=source_format,
            fragmentation=_fragmentation_report(reader),
        )
