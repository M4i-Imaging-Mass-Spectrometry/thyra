"""Essential metadata for an imzML, from the document head alone.

The metadata-only counterpart to :class:`ImzMLMetadataExtractor`, used by
``preview_msi`` (issue #360).  It subclasses it rather than restating it:
of the twenty-odd things that extractor reads, only three come from the
spectrum list -- the grid, the spectrum count and the mass range -- and
everything else (pixel size and its unit, the instrument configurations
the axis-type guess matches on, the declared mobility array, the UUID,
the declared representation) is read from ``parser.metadata``, which
:class:`~thyra.readers.imzml.header.ImzMLHeaderParser` provides in full.

So this class overrides exactly the methods that would have reached for a
spectrum, and inherits the rest unchanged.  That is the point: a second
implementation of "what is this file" would be a second answer to the
question, and the two would drift.

What it reports differently, and why that is not a drift:

- ``dimensions`` is the raster the file declares (``IMS:1000042`` /
  ``IMS:1000043``), not the extent of the coordinates present.  A
  conversion must keep using the coordinates -- it has to place every
  pixel it writes -- but a preview is describing an acquisition, and the
  two differ only when one was interrupted.
- ``total_peaks`` is 0, not counted, exactly as ``skip_total_peaks``
  does for Bruker: the preview card does not show it, and nothing in the
  resampling decision consults it (``avg_peaks_per_spectrum`` feeds only
  ``is_high_density_profile``, which no detector calls).
- ``coordinate_bounds`` is zeros and ``coordinate_offsets`` is ``None``.
  Both describe where the spectra sit, which is what was not read.
- ``dimensions[2]`` is 1, because imzML has no header term for the number
  of planes -- z is a property of the coordinates. So ``is_3d`` is False
  here for a 3D acquisition that the coordinate path would report as 3D.
  Nothing in the preview reads it (the card shows x and y), and nothing
  on the conversion path reaches this class.
- The spectrum type and the instrument configuration follow what the
  first spectrum states, where a conversion follows every converted
  spectrum (D78). They differ only for a file whose later spectra say
  otherwise, as for polarity (D35).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional, Sequence, Set, Tuple, cast

from ...readers.imzml.header import ImzMLHeaderParser, scan_observed_mz_range
from ..types import EssentialMetadata
from .imzml_extractor import ImzMLMetadataExtractor

logger = logging.getLogger(__name__)


class ImzMLHeaderExtractor(ImzMLMetadataExtractor):
    """Read an imzML's metadata without reading its spectra."""

    def __init__(
        self,
        parser: ImzMLHeaderParser,
        imzml_path: Path,
        spectrum_type: Optional[str] = None,
    ):
        """Initialize a header-only ImzML metadata extractor.

        Args:
            parser: The parsed document head.  Stands in for
                :class:`ImzMLParser`, which is why the base class's type
                is widened here rather than satisfied.
            imzml_path: Path to the ``.imzML`` file.
            spectrum_type: Optional override for the spectrum
                representation, as on the full extractor.
        """
        super().__init__(cast(Any, parser), imzml_path, spectrum_type=spectrum_type)
        self.header = parser

    def _extract_essential_impl(self) -> EssentialMetadata:
        """Build essential metadata from the head of the document."""
        n_x, n_y = self._declared_raster()
        mass_range = self._recorded_mass_range()
        n_spectra = self.header.n_spectra
        spectrum_type = self._detect_centroid_spectrum()

        return EssentialMetadata(
            dimensions=(n_x, n_y, 1),
            # Not read: see the module docstring.  Zeros are what an
            # unread extent reports everywhere in this package.
            coordinate_bounds=(0.0, 0.0, 0.0, 0.0),
            mass_range=mass_range if mass_range is not None else (0.0, 0.0),
            mass_range_known=mass_range is not None,
            pixel_size=self._extract_pixel_size_fast(),
            n_spectra=int(n_spectra) if n_spectra is not None else 0,
            n_spectra_counted=n_spectra is not None,
            total_peaks=0,
            source_path=str(self.imzml_path),
            coordinate_offsets=None,
            spectrum_type=spectrum_type,
            spectrum_type_guessed=self._spectrum_type_guessed,
        )

    def _declared_raster(self) -> Tuple[int, int]:
        """The raster the file declares, which the caller has checked for.

        :func:`head_shortfall` gates this path, so reaching it without a
        usable declaration is a bug rather than a bad file, and raises
        accordingly.
        """
        raster = declared_raster(self.header)
        if raster is None:
            raise ValueError(
                f"{self.imzml_path} declares no raster geometry "
                "(IMS:1000042/IMS:1000043); a metadata-only read cannot "
                "report its dimensions."
            )
        return raster

    def _get_mass_range_complete(self) -> Tuple[Tuple[float, float], int]:
        """The recorded mass range, never a scan of the binary.

        The inherited implementation walks the spectrum list and decodes
        the ``.ibd``, which is the cost this path exists to avoid.

        Raises:
            ValueError: If the file records no observed m/z range.  A
                caller asking for the range specifically -- rather than
                for the essential metadata, which reports it as unknown
                -- gets told, because the alternative is an invented
                range reaching an axis (see the base implementation).
        """
        mass_range = self._recorded_mass_range()
        if mass_range is None:
            raise ValueError(
                f"{self.imzml_path} records no observed m/z range "
                "(MS:1000528/MS:1000527), and a metadata-only read does "
                "not decode spectra to find one."
            )
        return mass_range, 0

    def _recorded_mass_range(self) -> Optional[Tuple[float, float]]:
        """The observed-m/z extrema the document records, or ``None``.

        The one call site for the scan on this path, so the mode reaches
        it once: a continuous file (``IMS:1000030``, every spectrum one
        m/z array) states its range on its first spectrum and the scan
        stops there; a processed file is scanned in full.  The mode comes
        off ``parser.metadata`` exactly as the full extractor reads it,
        which is what keeps the two paths agreeing about which files are
        which.
        """
        return scan_observed_mz_range(
            self.imzml_path, shared_axis=self._is_continuous_mode()
        )

    def _spectrum_count(self) -> int:
        """The count the file declares, not one taken over coordinates."""
        return int(self.header.n_spectra or 0)

    def _spectrum_polarity_accessions(self) -> Set[str]:
        """The polarity terms the first spectrum states, as pyimzml reads it.

        The head ends at the second spectrum. A conversion reads every
        spectrum's terms, so the two differ only for a file whose later
        spectra state the other polarity (D35).
        """
        return set(cast(ImzMLHeaderParser, self.parser).first_spectrum_polarity)

    def _converted_spectra(self) -> Optional[Sequence[int]]:
        """Every spectrum the head recorded, which is the first one.

        The spectrum type and the instrument configuration follow the
        first spectrum here; a conversion follows every converted
        spectrum, so the two differ only for a file whose later spectra
        say otherwise (D78, as D35 for polarity).
        """
        return None


def head_shortfall(parser: ImzMLHeaderParser) -> Optional[str]:
    """Why the head cannot describe this file, or ``None`` if it can.

    The gate on the metadata-only path, and the single place the reasons
    live -- so the line logged when a file falls back to the full parse
    says which one it was, rather than naming whichever was checked first.
    The head answers only when its declarations fit together (D63); any
    other file is read the full way, as a conversion reads it.

    - **It declares no spectra.**  ``<spectrumList count="0">``, a
      negative count, or no count at all.  The coordinate path refuses a
      file that really has none by name ("No coordinates found"), and it
      should keep refusing it rather than be described as an empty raster.

    - **It is numbered from 0**, which makes the raster declaration
      ambiguous.  ``max count of pixels x`` is defined as a count, but
      pyimzml's own writer fills it with the largest coordinate it saw --
      the same number for a 1-based file, one short for a 0-based one, so
      a 3x3 raster written from origin 0 declares ``2``.  Nothing in the
      head says which was meant, so the coordinate path settles the base
      exactly, as it has since issue #244.  A file that states no
      position at all on its first spectrum is ambiguous the same way and
      goes the same route, where it will fail honestly if it really has
      no coordinates.

    - **It declares no raster**, or a side below 1.  The spec puts
      ``IMS:1000042`` / ``IMS:1000043`` in ``<scanSettings>`` and every
      writer met so far fills them in, but nothing enforces it, and a side
      of 0 or less describes no acquisition.

    - **It declares more spectra than its raster holds.**  No acquisition
      fits that, and it is also how a 0-based file shows itself when its
      first spectrum is not at 0: pyimzml's writer declares the largest
      coordinate, one short per axis, so a 4 x 3 raster numbered from 0 and
      listed from its last pixel declares 3 x 2 for 12 spectra.  A count
      equal to the raster is the normal full raster and is trusted.

    Known limit: a wrong count that still fits its raster, and a 0-based
    file whose spectra fit inside its one-short declaration with the first
    spectrum away from 0, are still described from the head.  Both touch
    only the preview and the metadata document, never a store.
    """
    if parser.n_spectra is None:
        return "its spectrum list declares no count"
    if parser.n_spectra < 1:
        return f"its spectrum list declares {parser.n_spectra} spectra"

    position = parser.first_position
    if position is None or position[0] < 1 or position[1] < 1:
        return (
            f"its first spectrum is at {position}, so a declared raster "
            f"could be a count or a largest index"
        )

    raster = declared_raster(parser)
    if raster is None:
        return "it declares no raster geometry (IMS:1000042/IMS:1000043)"
    n_x, n_y = raster
    if n_x < 1 or n_y < 1:
        return f"it declares a raster of {n_x} x {n_y}"
    if parser.n_spectra > n_x * n_y:
        return (
            f"it declares {parser.n_spectra:,} spectra for a raster of "
            f"{n_x} x {n_y}, which holds {n_x * n_y:,}"
        )

    return None


def declared_raster(parser: ImzMLHeaderParser) -> Optional[Tuple[int, int]]:
    """``(n_x, n_y)`` as the document head declares them, or ``None``.

    Reads the declaration and nothing else.  Whether the declaration can
    be *trusted* is :func:`head_shortfall`'s question, and it is asked
    before this path is taken at all.
    """
    n_x = parser.imzmldict.get("max count of pixels x")
    n_y = parser.imzmldict.get("max count of pixels y")
    try:
        if n_x is None or n_y is None:
            return None
        return (int(n_x), int(n_y))
    except (TypeError, ValueError):
        logger.debug("Unusable declared raster: %r x %r", n_x, n_y)
        return None
