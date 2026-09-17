"""The two metadata records every format is reduced to.

:class:`EssentialMetadata` is what a conversion needs before it can start:
the grid, the extent, the mass range, the pitch, and the counts. Seven
formats fill in the same frozen dataclass, so everything downstream is
written against one shape rather than against whichever vendor produced
the file. :class:`ComprehensiveMetadata` wraps it and adds what the source
records for provenance, which is format-shaped by nature and stays in
dictionaries.

**A field here is a promise that survives to disk**, which is what makes
the absences deliberate. ``n_spectra_counted`` exists because 0 spectra
and "not counted" are different facts and a store that conflates them
misreports a PHI preview (issue #240); ``mass_range_known`` is the same
distinction for a range, and exists because an imzML that records no
observed-m/z cvParams leaves one genuinely unread (issue #360).
``pixel_size`` is the in-plane
raster pitch and says nothing about z, which is why ``z_spacing_um`` is
separate and usually ``None``. ``coordinate_offsets`` records what
rebasing the file's coordinates subtracted, so a store can say where its
origin came from.

Two fields were removed rather than kept for symmetry:
``peak_counts_per_pixel`` built a CSR ``indptr`` that went with D10, and
``estimated_memory_gb`` fed a size-based router that went with D11. Both
went on being computed and written into every store long after the thing
that read them was deleted (issues #271, #272). A field nothing reads is
not free here: it reaches disk and a consumer believes it.
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple


@dataclass(frozen=True)
class EssentialMetadata:
    """Critical metadata for processing decisions and interpolation setup.

    Attributes:
        dimensions: Grid dimensions as ``(x, y, z)``.
        coordinate_bounds: Spatial extent as ``(min_x, max_x, min_y, max_y)``.
        mass_range: Mass-to-charge range as ``(min_mz, max_mz)``.
        pixel_size: In-plane pixel dimensions as ``(x_um, y_um)`` in
            micrometres, or ``None`` when not detected.  This is the
            raster pitch and says nothing about z; see ``z_spacing_um``.
        n_spectra: Total number of spectra in the dataset -- the ones
            actually present, not the positions the raster covers.
            Meaningless unless ``n_spectra_counted`` is True.
        n_spectra_counted: Whether ``n_spectra`` was counted at all. False
            only on the metadata-only path of a format that cannot count
            spectra without decoding them: PHI stores a stream of ion
            events rather than a list of spectra, so which pixels carry one
            is a property of the data, not of the header (issue #240). When
            False, ``n_spectra`` and ``total_peaks`` are 0 meaning "not
            counted", which is not the same as "none present" -- callers
            that display a count must say so rather than print the zero.
        total_peaks: Total number of peaks across all spectra (used for
            sparse matrix pre-allocation).
        mass_range_known: Whether ``mass_range`` was read at all. False
            only on the metadata-only path of a file that records no
            range of its own: an imzML written without the observed-m/z
            cvParams (``MS:1000528``/``MS:1000527``) states its extrema
            nowhere but in the binary, and a preview does not decode
            spectra to find them (issue #360). When False, ``mass_range``
            is ``(0.0, 0.0)`` meaning "not read", which is not a range
            any acquisition has -- callers must say "unknown" rather
            than print it.
        source_path: Absolute path to the source data.
        coordinate_offsets: Raw coordinate offsets ``(x, y, z)`` used to
            normalise coordinates to 0-based indexing.
        spectrum_type: Spectrum type string (e.g. ``"centroid spectrum"``),
            used to guide resampling decisions.
        z_spacing_um: Distance between consecutive slices in micrometres,
            or ``None`` when the source cannot report it.  This is a
            *physical* distance set by how the sections were cut, and is
            unrelated to the in-plane pitch in ``pixel_size`` -- the two
            agree only by coincidence.  No reader populates this today:
            imzML has no standard term for it, and 3D acquisitions are
            frequently non-consecutive sections, so the value usually has
            to be supplied by the user.  The field exists so a format that
            *does* record it has somewhere to report it, and so the
            converter's precedence chain has something to consult.
    """

    dimensions: Tuple[int, int, int]
    coordinate_bounds: Tuple[float, float, float, float]
    mass_range: Tuple[float, float]
    pixel_size: Optional[Tuple[float, float]]
    n_spectra: int
    total_peaks: int
    source_path: str
    coordinate_offsets: Optional[Tuple[int, int, int]] = None
    spectrum_type: Optional[str] = None
    z_spacing_um: Optional[float] = None
    n_spectra_counted: bool = True
    mass_range_known: bool = True

    @property
    def has_pixel_size(self) -> bool:
        """Check if pixel size information is available."""
        return self.pixel_size is not None

    @property
    def is_3d(self) -> bool:
        """Check if dataset is 3D (z > 1)."""
        return self.dimensions[2] > 1


@dataclass
class ComprehensiveMetadata:
    """Complete metadata including format-specific details.

    Wraps :class:`EssentialMetadata` and adds vendor-specific information
    that is not needed for conversion but useful for provenance and QC.

    Attributes:
        essential: Core metadata required for conversion.
        format_specific: Vendor-specific metadata (e.g. ImzML CV params,
            Bruker property tables).
        acquisition_params: Acquisition parameters such as polarity,
            scan range, and laser settings.
        instrument_info: Instrument model, serial number, and software
            version.
        raw_metadata: Unprocessed metadata exactly as read from the
            source file, preserved for round-trip fidelity.
    """

    essential: EssentialMetadata
    format_specific: Dict[str, Any]
    acquisition_params: Dict[str, Any]
    instrument_info: Dict[str, Any]
    raw_metadata: Dict[str, Any]

    @property
    def dimensions(self) -> Tuple[int, int, int]:
        """Convenience access to dimensions from essential metadata."""
        return self.essential.dimensions

    @property
    def pixel_size(self) -> Optional[Tuple[float, float]]:
        """Convenience access to pixel size from essential metadata."""
        return self.essential.pixel_size

    @property
    def coordinate_bounds(self) -> Tuple[float, float, float, float]:
        """Convenience access to coordinate bounds from essential metadata."""
        return self.essential.coordinate_bounds
