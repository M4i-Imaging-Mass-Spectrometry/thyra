# thyra/converters/spatialdata/base_spatialdata_converter.py

"""The SpatialData write path every converter shares.

**Exceptions are logged as text here, never as objects.** Every
``logger.warning("...: %s", str(e))`` in this file could read ``e`` and
render the same line -- but a log record keeps its arguments, a handler
that keeps records keeps the record, and an exception object drags along
its traceback, the frames reachable through ``tb_frame.f_back`` and every
memmap those frames hold. The finalize path builds each table's AnnData
over the CSC scratch memmaps, so one retained record is enough to keep the
scratch directory mapped, and Windows will not delete a mapped file: the
directory survives the conversion and the user is told to remove it by
hand. pytest's ``caplog`` retains records, and so does Ousia's per-session
log capture, which is the consumer that met it (issue #249).
"""

import logging
import warnings
from abc import ABC, abstractmethod
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import geopandas as gpd
import numpy as np
import pandas as pd
import zarr
from anndata import AnnData
from numpy.typing import NDArray
from shapely.geometry import box
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity

from ...core.base_converter import BaseMSIConverter, PixelSizeSource
from ...core.base_reader import BaseMSIReader
from ...core.conversion_state import ConversionState
from ...errors import ConversionRefused
from ...metadata.root_attrs import RootAttrsBuilder, RootAttrsContext
from ...metadata.schema import MSI_VAR_RESERVED_COLUMNS
from ...metadata.uns_assembler import UnsAssembler, UnsContext
from ...resampling.axis_planner import AxisPlanner, normalize_resampling_config
from ...resampling.mobility_grid import MOBILITY_CHANNELS
from ...resampling.strategies import ResamplingStrategy, build_strategy
from ...resampling.types import AxisLinearisation, ResamplingConfig
from ._chunking import table_write_config
from .optical_image import OpticalImages
from .sibling_tables import SiblingContext, SiblingTables

logger = logging.getLogger(__name__)


@contextmanager
def _suppress_upstream_warnings():
    """Suppress known upstream warnings from ome_zarr, zarr v3, and spatialdata."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Passing storage-related arguments",
            category=FutureWarning,
        )
        warnings.filterwarnings(
            "ignore", message="Object at.*is not recognized", category=UserWarning
        )
        warnings.filterwarnings(
            "ignore",
            message="Consolidated metadata is currently not",
            category=UserWarning,
        )
        yield


class BaseSpatialDataConverter(BaseMSIConverter, ABC):
    """Base converter for MSI data to SpatialData format with shared functionality."""

    def __init__(
        self,
        reader: BaseMSIReader,
        output_path: Path,
        dataset_id: str = "msi_dataset",
        pixel_size_um: float = 1.0,
        pixel_size_source: PixelSizeSource = PixelSizeSource.DEFAULT,
        handle_3d: bool = False,
        z_spacing_um: Optional[float] = None,
        pixel_size_detection_info: Optional[Dict[str, Any]] = None,
        resampling_config: Optional[Union[Dict[str, Any], ResamplingConfig]] = None,
        include_optical: bool = True,
        apply_optical_alignment: bool = True,
        write_mobility_table: bool = True,
        mobility_heatmap: bool = True,
        mobility_grid: bool = False,
        mobility_bins: int = MOBILITY_CHANNELS,
        mobility_min: Optional[float] = None,
        mobility_max: Optional[float] = None,
        msms_table: bool = True,
        **kwargs: Any,
    ) -> None:
        """Initialize the base SpatialData converter.

        Args:
            reader: MSI data reader
            output_path: Path for output file
            dataset_id: Identifier for the dataset
            pixel_size_um: In-plane size of each pixel in micrometers
            pixel_size_source: How pixel size was determined
            handle_3d: Whether to process as 3D data (True) or 2D slices
                (False)
            z_spacing_um: Distance between consecutive slices in
                micrometers.  Only meaningful when a true volume is
                written (``handle_3d=True`` and more than one slice).
                ``None`` (default) falls back to the in-plane pitch and
                records that as an assumption rather than a measurement;
                see :meth:`BaseMSIConverter._resolve_z_spacing`.
            pixel_size_detection_info: Optional metadata about pixel size
                detection
            resampling_config: Optional resampling configuration dict
            include_optical: Whether to include optical images in output
                (default: True)
            write_mobility_table: When the reader shares one set of
                (m/z, mobility) feature pairs across pixels, also write
                them as a mobility-resolved sibling table
                (``{table}_mobility``) beside the summed MSI table
                (default: True). Never changes the MSI table itself.
            mobility_heatmap: When the reader has an ion mobility
                dimension, accumulate the mean mass-mobility frame from
                the raw scan read and store it on the summed table as
                ``uns["mobility_heatmap"]`` (default: True). Costs one
                extra pass over the source; see ``mobility_heatmap.py``.
            mobility_grid: When the reader carries mobility per pixel
                rather than as a shared feature axis (Bruker TDF), bin
                the point cloud onto a common mobility grid and write the
                result as the same ``{table}_mobility`` sibling
                (default: False -- opt in, it costs an extra pass and a
                far larger table). Ignored by a source that already
                shares a feature axis, which needs no grid.
            mobility_bins: Channels the grid divides the mobility range
                into (default: 256). The default is the mass-mobility
                heatmap's own channel count over the same edges, which is
                what lets a box drawn on the heatmap index grid channels
                directly; changing it gives that up.
            mobility_min: Lower edge of the grid, in the axis unit.
                ``None`` (default) takes the smallest mobility value the
                source's axis actually holds, which is also where the
                heatmap starts.
            mobility_max: Upper edge of the grid; ``None`` takes the
                largest value the axis holds.
            msms_table: When the source isolates several precursors per
                pixel in disjoint mobility slices (Bruker PASEF), also
                write them split apart as a demultiplexed sibling table
                (``{table}_msms``) beside the summed MSI table (default:
                True -- the summed spectrum of such a pixel is a mixture
                of unrelated fragment spectra, so the split is the
                accurate representation; design decision D2). Refused
                with a reason rather than approximated when the schedule
                is not separable, which is also what happens on every
                source that is not MS/MS, so the default costs a source
                that has nothing to split nothing; see ``msms_table.py``.
                Never changes the MSI table itself.
            apply_optical_alignment: If True (default) and the MSI source
                has FlexImaging Area metadata, compute an alignment that
                places MSI raster coordinates in optical-image pixel
                space.  Set to False when a downstream tool owns the
                alignment (e.g. Ousia's wizard registers MSI to Xenium
                via EscDat and does not want Thyra to pre-rotate the
                MSI into FlexImaging's optical frame).  When False the
                MSI elements land in pure micrometer coordinates at
                ``"global"``.
            **kwargs: Additional keyword arguments

        Raises:
            ConversionRefused: If pixel_size_um is not positive, dataset_id
                is empty, or ``sparse_format`` is passed
        """
        # ``sparse_format`` chose between CSC and CSR until v3.22. CSC is now
        # the only layout written, so the keyword has nothing left to select --
        # but an unknown keyword lands in ``self.options`` without a word, and a
        # caller who asked for CSR would get CSC and no signal. That is the
        # failure this removal was meant to end, so say it instead.
        if "sparse_format" in kwargs:
            raise ConversionRefused(
                "sparse_format was removed: every converter writes CSC, which "
                "is the layout an ion image reads down. Drop the argument; for "
                "row-wise access call X.tocsr() on the matrix you read back."
            )

        # Validate inputs
        if pixel_size_um <= 0:
            raise ConversionRefused(
                f"pixel_size_um must be positive, got {pixel_size_um}"
            )
        if not dataset_id.strip():
            raise ConversionRefused("dataset_id cannot be empty")

        # Extract pixel_size_detection_info from kwargs if provided
        kwargs_filtered = dict(kwargs)
        if (
            pixel_size_detection_info is None
            and "pixel_size_detection_info" in kwargs_filtered
        ):
            pixel_size_detection_info = kwargs_filtered.pop("pixel_size_detection_info")

        super().__init__(
            reader,
            output_path,
            dataset_id=dataset_id,
            pixel_size_um=pixel_size_um,
            pixel_size_source=pixel_size_source,
            handle_3d=handle_3d,
            z_spacing_um=z_spacing_um,
            **kwargs_filtered,
        )

        self._non_empty_pixel_count: int = 0
        # The m/z range a peak has to be inside to be kept, as opposed to
        # the span of the axis points themselves. Adopted from the planner
        # in _setup_mass_axis(); ``None`` means "no resampled axis was
        # built", and the span is then the axis's own. See
        # :func:`thyra.resampling.binning.kept_mz_range`.
        self._axis_range: Optional[Tuple[float, float]] = None
        # The coordinate the common mass axis may be indexed through, once
        # the planner has built the axis and usable_linearisation() has
        # checked it; None means every placement onto it is a search
        # (D21). It belongs to the axis, not to the
        # strategy that bins with it, which is why the sibling sinks read
        # it here rather than off ``_resampler``: an interpolated
        # conversion still writes sibling tables, and they still place
        # peaks onto this axis.
        self._axis_linearisation: Optional[AxisLinearisation] = None
        # The operator that places one spectrum onto the common mass axis:
        # built beside the axis in _setup_mass_axis(), from the method the
        # config or the detector chain settled on. ``None`` means no
        # resampled axis was built, and the spectrum is mapped onto the
        # reader's own axis instead (--no-resample).
        self._resampler: Optional[ResamplingStrategy] = None
        # Intensities dropped for being non-finite or negative, and
        # whether that has been said yet. See _count_unusable_intensities().
        self._unusable_intensities: int = 0
        self._unusable_intensities_warned: bool = False
        self._pixel_size_detection_info = pixel_size_detection_info
        # ``convert.py`` detects an (x, y) pair, hands the converter the x
        # half as ``pixel_size_um`` and puts both in the detection info.
        # Reading the y half back here is what stops the store describing
        # an anisotropic raster as square (issue #228); everything the
        # converter writes takes its y pitch from this attribute.
        #
        # Only for a pitch that was actually detected. A caller who states
        # one gets it on both axes -- that is what stating it means -- and
        # a source ``convert.py`` could not detect leaves the placeholder,
        # which :meth:`_adopt_detected_pixel_size` settles from the reader
        # once its metadata is loaded.
        detected_y = (pixel_size_detection_info or {}).get("detected_y_um")
        if (
            detected_y is not None
            and pixel_size_source is PixelSizeSource.AUTO_DETECTED
        ):
            self.pixel_size_y_um = float(detected_y)
            self._log_anisotropic_raster()
        self._resampling_config = (
            normalize_resampling_config(resampling_config)
            if resampling_config is not None
            else None
        )
        # What the planner resolved, adopted in _setup_mass_axis();
        # consumed by the uns assembler's processing provenance, through
        # _uns_context().
        self._resolved_resampling_plan: Optional[Dict[str, Any]] = None

        # The tables written beside the summed one -- the mobility-resolved
        # sibling and the demultiplexed MS/MS sibling -- and everything that
        # decides them: the mobility grid one of them bins onto, the raw
        # pass they share with the mass-mobility heatmap, and the scratch
        # directories every table's memmaps live in until it is written.
        # Everything it reads that a conversion decides later travels per
        # call in a :class:`~thyra.converters.spatialdata.sibling_tables.
        # SiblingContext` (see :meth:`_sibling_context`).
        self.siblings = SiblingTables(
            self.reader,
            self._store_path,
            write_mobility_table=write_mobility_table,
            mobility_heatmap=mobility_heatmap,
            mobility_grid=mobility_grid,
            mobility_bins=mobility_bins,
            mobility_min=mobility_min,
            mobility_max=mobility_max,
            write_msms_table=msms_table,
        )

        # Everything that decides which m/z values the store will carry:
        # the method, the axis law, the bin width and the range, and the
        # caches the detectors are read through. Built here because an
        # explicit ``--resample-method`` is checked against the detector
        # chain while there is still someone to warn (issue #246); the axis
        # itself is laid in _setup_mass_axis(), once the reader's metadata
        # is loaded.
        self.axis_planner = AxisPlanner(self.reader, self._resampling_config)

        # The store's optical images and everything they own: the
        # FlexImaging alignment, the TIC-to-image affine, which file became
        # which element and the placeholders waiting for their pixels. Built
        # last because it is handed this converter's pitch accessor, which
        # only reads a settled pitch once the reader's metadata is in.
        # Handed the store accessor for the same reason: the pixels stream
        # after the write, and ``output_path`` can be reassigned in between.
        self.optical = OpticalImages(
            self.reader,
            self._store_path,
            self.dataset_id,
            include=include_optical,
            apply_alignment=apply_optical_alignment,
            pixel_size_xy=self._resolved_pixel_size_xy,
        )

        # The table's ``uns`` block and everything that composes it. Handed
        # the same pitch accessor for the same reason, and the three other
        # things it needs that are settled here; everything it reads that a
        # conversion decides later travels per call in an
        # :class:`~thyra.metadata.uns_assembler.UnsContext` (see
        # :meth:`_uns_context`).
        self.uns = UnsAssembler(
            self.reader,
            pixel_size_xy=self._resolved_pixel_size_xy,
            pixel_size_detection_info=self._pixel_size_detection_info,
            resampling_config=self._resampling_config,
        )

        # The store's own attrs, the assembler's sibling. The optical
        # collaborator reaches it per call rather than here: whether the
        # raster landed in optical-photo pixel space decides the
        # coordinate contract, and the converter may replace the
        # collaborator after this point.
        self.root_attrs = RootAttrsBuilder(
            self.reader,
            dataset_id=self.dataset_id,
            pixel_size_detection_info=self._pixel_size_detection_info,
            handle_3d=self.handle_3d,
            conversion_options=self.options,
        )

    def convert(self) -> bool:
        """Run the base workflow; release the tables' scratch on every exit path.

        The scratch cleanup must run on success, on an exception and on a
        KeyboardInterrupt alike: a route that released its temp directory
        on the success path only leaked 79.5 GiB into one user's system
        temp before a manual sweep.
        """
        try:
            return super().convert()
        finally:
            self.siblings.release_scratch()

    def build_uns_metadata(self) -> Dict[str, Any]:
        """The provenance block every table of this conversion carries.

        The converter's orchestration point: it packs what it currently
        knows into an :class:`~thyra.metadata.uns_assembler.UnsContext`
        and hands it to :meth:`~thyra.metadata.uns_assembler.UnsAssembler.build`,
        which composes the mapping. See that method for what lands in it
        and why there is exactly one place it is composed.

        Returns:
            Mapping of ``uns`` key to the value to store. Empty if the
            reader cannot produce comprehensive metadata at all.
        """
        return self.uns.build(self._uns_context())

    def _uns_context(self) -> UnsContext:
        """What the assembler needs that this conversion decided after setup.

        Packed fresh per call, never cached: the sibling keys are decided
        per slice and taken back when a builder declines, the resolved
        plan and the mobility grid are filled while the conversion runs,
        and ``pixel_size_source`` is reassigned when a pitch is detected
        from the reader's metadata.

        ``_region_info`` is read through ``getattr`` because a converter
        that never ran :meth:`_initialize_conversion` -- which is how
        several tests ask for the block -- does not have it yet, and the
        block is written from whatever regions are known, or omitted.
        """
        return UnsContext(
            mobility_table_key=self.siblings.mobility_table_key,
            msms_table_key=self.siblings.msms_table_key,
            mobility_grid=self.siblings.mobility_grid,
            resolved_resampling_plan=self._resolved_resampling_plan,
            region_info=getattr(self, "_region_info", None),
            pixel_size_source=self.pixel_size_source,
            mobility_heatmap=self._ensure_mobility_heatmap,
        )

    def _sibling_context(self) -> SiblingContext:
        """What the sibling tables need that this conversion decided after setup.

        Packed fresh per call for the same reason :meth:`_uns_context` is:
        the common mass axis and its linearisation are laid in
        :meth:`_setup_mass_axis`, the grid dimensions in
        :meth:`_initialize_conversion`, and the three delegates are read
        off collaborators a test may have replaced since.
        """
        return SiblingContext(
            mass_axis=self._common_mass_axis,
            linearisation=self._axis_linearisation,
            dimensions=self._dimensions,
            n_spectra=self._get_total_spectra_count,
            build_uns=self.build_uns_metadata,
            fragmentation=self.uns.fragmentation,
            unname_declined=self.uns.unname_declined_siblings,
        )

    def _ensure_mobility_heatmap(self) -> Optional[Dict[str, Any]]:
        """The conversion's mass-mobility heatmap block, built on first demand.

        The converter's orchestration point, as
        :meth:`build_uns_metadata` is the ``uns`` assembler's: it packs
        what it knows and hands it over. The block is a sibling concern --
        built from the raw mobility pass the sibling tables share -- that
        lands in the summed table's own ``uns``, which is why the
        assembler asks for it through a callable.
        """
        return self.siblings.heatmap(self._sibling_context())

    def _resolved_pixel_size_xy(self) -> Tuple[float, float]:
        """The in-plane pixel pitch as ``(x_um, y_um)``.

        One pair, and every block of the store is written from it: the
        root attrs, ``coordinate_systems.global`` and its affine, the
        image and shapes transformations, the pixel footprints,
        ``obs["spatial_x"]``/``["spatial_y"]`` and the ``msi_metadata``
        block.

        It used to be only the last of those. The converter carried a
        single float -- detection kept the source's x pitch and discarded
        the y one -- while this method read the true pair back out of the
        detection info for the metadata block alone. A raster acquired at
        30 x 50 um was then rendered at 30 x 30 and the store's own blocks
        contradicted each other (issue #228).
        """
        return (float(self.pixel_size_um), float(self.pixel_size_y_um))

    def _store_path(self) -> Path:
        """Where this conversion's store is being written, right now.

        ``output_path`` is not fixed at construction. ``convert_msi``
        shortens it for Windows (``prepare_zarr_output_path``) *before*
        building the converter, so the production route is safe, but a
        caller that stands a converter up itself has to shorten afterwards
        -- the CLI's ``_quarantine_partial_output`` says as much, and the
        suite does it in ``test_internal_failures_are_not_warnings``.

        The one accessor both collaborators that outlive construction are
        handed, rather than a snapshot each: :class:`SiblingTables` makes
        every table's scratch in its parent, and :class:`OpticalImages`
        streams its pixels into it once the write has returned. What the
        converter does with the path itself -- the write, the consolidate
        -- reads ``output_path`` directly, at call time, which is the same
        thing.
        """
        return self.output_path

    def _add_metadata_to_uns(self, adata) -> None:
        """Apply :meth:`build_uns_metadata` to an AnnData about to be written."""
        self.uns.apply(adata, self._uns_context())

    def _initialize_conversion(self) -> None:
        """Override parent initialization to preserve resampled mass axis."""
        logger.info("Loading essential dataset information...")
        try:
            # Load essential metadata first (fast, single query for Bruker)
            essential = self.reader.get_essential_metadata()
            # Hand it to the planner rather than let it ask again: on a
            # Bruker source that is a query, not an attribute.
            self.axis_planner.adopt_essential_metadata(essential)

            self._dimensions = essential.dimensions
            if any(d <= 0 for d in self._dimensions):
                raise ConversionRefused(
                    f"Invalid dimensions: {self._dimensions}. All dimensions "
                    f"must be positive."
                )

            # Store essential metadata for use throughout conversion
            self._coordinate_bounds = essential.coordinate_bounds
            self._n_spectra = essential.n_spectra

            # Override pixel size only if using default and metadata is available
            self._adopt_detected_pixel_size(essential)

            # After pixel size, because the fallback is the pixel size.
            self._resolve_z_spacing(essential)

            # Handle mass axis setup
            self._setup_mass_axis()

            # Only load comprehensive metadata if needed (lazy loading)
            self._metadata = None  # Will be loaded on demand

            # Fetch region data from reader (available for multi-region datasets)
            # Always populate region metadata for a consistent schema:
            # - obs["region_number"] exists on every dataset
            # - uns["regions"] always has at least one entry
            # Single-region or unknown-region datasets default to region 1.
            self._region_map = self.reader.get_region_map()
            self._region_info = self.reader.get_region_info()
            if self._region_info:
                logger.info(
                    f"Region metadata available: {len(self._region_info)} regions"
                )
            else:
                self._region_info = [{"region_number": 1, "n_spectra": self._n_spectra}]

            # Compute optical alignment for FlexImaging data when
            # available.  The alignment data is computed REGARDLESS of
            # apply_optical_alignment so it can be used to place the
            # FlexImaging optical image relative to the MSI -- even
            # when a downstream tool (e.g. Ousia's wizard) is doing
            # its own MSI-to-target registration and does not want
            # Thyra to pre-rotate the MSI itself.  The flag only
            # controls whether MSI elements use the alignment; the
            # optical image always uses it (if available) to land in
            # the same "global" frame as the MSI.
            self.optical.compute_alignment()
            self.optical.build_tic_to_image_affine()
            if not self.optical.apply_alignment:
                logger.info(
                    "apply_optical_alignment=False: MSI elements will "
                    "use micrometer coordinates; optical image (if "
                    "any) will use inverse-alignment to land in the "
                    "same micrometer frame."
                )

            logger.info(f"Dataset dimensions: {self._dimensions}")
            logger.info(f"Coordinate bounds: {self._coordinate_bounds}")
            logger.info(f"Total spectra: {self._n_spectra}")
            if self._common_mass_axis is None:
                raise RuntimeError("Common mass axis is None after initialization")
            logger.info(f"Common mass axis length: {len(self._common_mass_axis)}")
        except ConversionRefused:
            # Said once, by whoever catches it. Re-prefixing a refusal
            # with the stage it came from adds nothing a user can act on.
            raise
        except Exception as e:
            logger.error(f"Error during initialization: {e}")
            raise

    def _setup_mass_axis(self) -> None:
        """Adopt the common mass axis the planner settles, and check it.

        The planner decides which m/z values the store carries -- built
        from the resampling plan, or the source's own under
        ``--no-resample`` -- and hands back the axis together with the
        three things the rest of the conversion has to know about it: the
        range the bins were laid across, the coordinate the axis may be
        indexed through, and the method the peaks are placed with. They
        are adopted together, because the four must never describe
        different axes.

        The memory budget travels the other way: the planner asks
        :meth:`_refuse_wide_mass_axis` at the moment a bin count becomes
        knowable on each route, which is before a resampled axis is
        materialised and after a raw one has been handed over (issue
        #251). It stays here because what it projects is the per-bin cost
        of the CSC assembly this converter writes, not anything about the
        axis.
        """
        settled = self.axis_planner.settle(self._refuse_wide_mass_axis)
        self._common_mass_axis = settled.axis
        self._axis_range = settled.axis_range
        self._axis_linearisation = settled.linearisation
        self._resolved_resampling_plan = settled.provenance

        # The axis and the operator that places spectra onto it are
        # adopted together and never separately: the triple below is
        # everything a strategy is allowed to know about the axis, and the
        # method and the gap tolerance are everything it is allowed to
        # know about the request. Nothing else in the converter decides
        # what resampling does to a spectrum.
        config = self._resampling_config
        if settled.method is None or config is None:
            self._resampler = None
        else:
            self._resampler = build_strategy(
                settled.method,
                settled.axis,
                settled.axis_range,
                settled.linearisation,
                config.gap_tolerance_da,
            )

        self._refuse_non_finite_axis()

    def _refuse_wide_mass_axis(self, n_bins: int) -> None:
        """Refuse a mass axis whose per-bin cost will not fit in memory.

        Raises:
            ConversionRefused: When the projection exceeds the fraction of
                free memory a conversion may take.
        """
        from .csc_assembly import mass_axis_refusal

        refusal = mass_axis_refusal(int(n_bins))
        if refusal is not None:
            raise ConversionRefused(refusal)

    def _refuse_non_finite_axis(self) -> None:
        """Refuse a mass axis carrying NaN or infinity.

        ``var.mz`` is this axis. Nothing on the way here tested it for
        finiteness, so a source with NaN m/z converted at exit 0 and
        ``thyra validate`` then refused the result with
        ``ERROR var.mz: 'mz' contains non-finite values``: the converter
        wrote a store its own validator rejects (issue #248). Dropping
        the values instead is not an option on the raw path -- the axis
        *is* the source's m/z values, and every spectrum's indices are
        paired with a full intensity array -- so this is a refusal.

        One pass over the axis, which is at most a few million float64s
        and is walked several times either way.
        """
        axis = self._common_mass_axis
        if axis is None or len(axis) == 0:
            return
        finite = np.isfinite(axis)
        if bool(finite.all()):
            return
        n_bad = int(axis.size - finite.sum())
        raise ConversionRefused(
            f"The common mass axis holds {n_bad:,} non-finite value(s) out of "
            f"{axis.size:,} (NaN or infinity). It becomes var['mz'], where "
            "every value has to be a real mass, so this store would be "
            "written and then refused by `thyra validate`. On a resampled "
            "conversion, check --resample-min-mz / --resample-max-mz; on "
            "--no-resample, the axis is the source's own m/z values and the "
            "source is what carries them."
        )

    def _drop_unusable_intensities(
        self, mzs: NDArray[np.float64], intensities: NDArray[np.float64]
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Drop intensities that are not a measurement, before resampling.

        Two kinds, one rule (issue #248):

        *Non-finite.* A NaN intensity reached ``X.data`` unchallenged,
        turned ``uns["average_spectrum"]`` into NaN, and passed
        validation -- one unusable value poisoning every aggregate over
        the column it sits in.

        *Negative.* An intensity is a count of ions, or an ADC reading
        proportional to one. A negative value is not a smaller
        measurement; it is a baseline subtraction that overshot, and no
        consumer of a converted store reads it as signal -- a TIC, an ion
        image and a mean spectrum all take it as removing current that
        was never there. Keeping it also made the two resampling methods
        disagree in kind rather than in detail: nearest_neighbor stored
        it (giving a pixel a TIC of 0), while tic_preserving found a
        non-positive total for the spectrum and zeroed the whole thing.

        Dropped here, before either method runs, so both see the same
        spectrum and produce the same store. Dropping rather than
        refusing, because baseline-subtracted profile data really does
        carry small negatives and refusing it would leave no way to
        convert those files at all; the count says how much went.

        The arrays are returned unchanged -- the same objects -- when
        there is nothing to drop, which is every ordinary spectrum, so
        the shared-axis identity fast path downstream is untouched.
        """
        if intensities.size == 0:
            return mzs, intensities
        keep = np.isfinite(intensities) & (intensities >= 0)
        if bool(keep.all()):
            return mzs, intensities
        self._count_unusable_intensities(
            int(intensities.size - keep.sum()),
            int(intensities.size),
            intensities[~keep],
        )
        return mzs[keep], intensities[keep]

    def _count_unusable_intensities(
        self, n_dropped: int, n_total: int, dropped: NDArray[np.float64]
    ) -> None:
        """Record intensities dropped as unusable, warning once.

        Once per conversion rather than once per spectrum, for the reason
        ``ResamplingStrategy``'s out-of-range counter gives: a source with
        negatives usually has them in every spectrum.
        ``_unusable_intensities`` keeps the running total, and like that
        counter it counts *calls*, so it double-counts a two-pass
        conversion; read the warning for the per-spectrum figure.
        """
        self._unusable_intensities += n_dropped
        if self._unusable_intensities_warned:
            return
        self._unusable_intensities_warned = True
        n_nan = int(np.count_nonzero(~np.isfinite(dropped)))
        logger.warning(
            "Dropping intensities that are not a measurement -- %d of %d in "
            "the first spectrum affected (%d non-finite, %d negative). A NaN "
            "would make every aggregate over its column NaN, and a negative "
            "value is a baseline subtraction that overshot, not signal. Both "
            "resampling methods see the spectrum without them, so they agree "
            "on what is stored.",
            n_dropped,
            n_total,
            n_nan,
            n_dropped - n_nan,
        )

    @staticmethod
    def _coalesce_duplicate_bins(
        indices: NDArray[np.int_], intensities: NDArray[np.float64]
    ) -> Tuple[NDArray[np.int_], NDArray[np.float64]]:
        """Sum intensities that landed on the same axis bin within one spectrum.

        Without resampling every m/z maps to its own axis entry, so this is
        the identity for ordinary data. A spectrum that repeats an m/z value
        -- an ion mobility export lists a feature once per mobility -- maps
        two entries to one column, and writing both would leave a
        duplicate ``(row, col)`` in the sparse matrix (the scatter also
        requires a row's bins to be unique). Summing here is what the old
        in-memory route did through scipy's COO conversion, so a store
        reads back the same.
        """
        if indices.size < 2:
            return indices, intensities
        # Compared pairwise rather than through ``np.diff``. Measured at
        # 331k bins, which is a realistic axis: 0.521 ms -> 0.060 ms hot
        # (8.8x), 0.630 -> 0.147 cold (4.3x). This runs once per spectrum
        # on the whole axis, for every spectrum of a --no-resample
        # conversion (issue #311).
        #
        # The reason is the temporary, not short-circuiting. ``np.diff``
        # materialises an 8-byte intp array the length of the axis where
        # the comparison writes a 1-byte bool -- measured, 2,647,992 bytes
        # against 330,999. Stopping early buys almost nothing by
        # comparison (1.06-1.13x between a duplicate at the front and a
        # fully ascending run), and on the memoised arange this mostly
        # runs against there is no False to stop at anyway.
        #
        # The answer is identical: same predicate, same pairs.
        if bool((indices[1:] > indices[:-1]).all()):
            return indices, intensities
        unique, inverse = np.unique(indices, return_inverse=True)
        summed = np.bincount(
            np.asarray(inverse).ravel(),
            weights=np.asarray(intensities, dtype=np.float64),
            minlength=unique.size,
        )
        return unique.astype(indices.dtype, copy=False), summed

    def build_region_numbers(self, x_values, y_values) -> NDArray[np.int32]:
        """``obs["region_number"]`` for the given pixel positions, in row order.

        Written on every dataset for a consistent schema, so a consumer
        does not have to branch on whether the acquisition had regions:
        without a region map every pixel is region 1, which is also what
        ``uns["regions"]`` reports for that case. Only the Bruker timsTOF
        reader produces a map today; a position missing from it gets -1.

        Kept as one method because the three write paths that used to
        build an obs table each had their own copy of this rule, and the
        hand-written streaming layout had no copy at all and simply
        omitted the column. One obs builder remains; the rule stays here
        so the sibling tables' builders can share it too.

        Args:
            x_values: Pixel x index per obs row.
            y_values: Pixel y index per obs row.

        Returns:
            Region number per obs row.
        """
        n = len(x_values)
        region_map = getattr(self, "_region_map", None)
        if region_map is None:
            return np.ones(n, dtype=np.int32)

        keys = zip(x_values.tolist(), y_values.tolist())
        return np.fromiter(
            (region_map.get(key, -1) for key in keys), dtype=np.int32, count=n
        )

    def _create_mass_dataframe(self) -> pd.DataFrame:
        """Create m/z dataframe for variable metadata.

        Returns:
            DataFrame with m/z values

        Raises:
            ValueError: If common mass axis is not initialized
        """
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized")

        n_channels = len(self._common_mass_axis)
        columns: Dict[str, Any] = {"mz": self._common_mass_axis}
        columns.update(self._validated_mass_axis_annotations())

        return pd.DataFrame(
            columns,
            index=[f"mz_{i}" for i in range(n_channels)],
        )

    def _validated_mass_axis_annotations(self) -> Dict[str, Any]:
        """Reader-supplied per-channel columns that fit the axis being written.

        A reader whose native axis is not m/z (flight time, drift time) can
        keep that axis alongside ``mz`` so the conversion stays reversible.
        Annotations whose length does not match are dropped rather than
        raising: that is the expected outcome when resampling rebuilds the
        axis, and it must not fail an otherwise good conversion.
        """
        if self._common_mass_axis is None:
            return {}
        n_channels = len(self._common_mass_axis)

        # getattr rather than a bare call: readers predating this hook, and
        # test doubles that do not subclass BaseMSIReader, simply have nothing
        # to contribute and should not have to raise to say so.
        getter = getattr(self.reader, "get_mass_axis_annotations", None)
        if getter is None:
            return {}

        try:
            annotations = getter()
        except Exception as exc:  # pragma: no cover - reader-defined
            # Format eagerly. Passing the exception itself to the logger keeps
            # it alive inside the log record, and with it its traceback, the
            # caller frames reachable through tb_frame.f_back, and any memmap
            # those frames hold -- which on Windows leaves the backing file
            # locked long after the conversion has finished.
            detail = f"{type(exc).__name__}: {exc}"
            logger.warning("Reader failed to supply mass axis annotations: %s", detail)
            return {}

        validated: Dict[str, Any] = {}
        for name, values in (annotations or {}).items():
            if name == "mz":
                logger.warning("Ignoring mass axis annotation named 'mz'")
                continue
            if name in MSI_VAR_RESERVED_COLUMNS:
                logger.info(
                    "Mass axis annotation %r uses a var column name the "
                    "metadata spec reserves for annotation results; it must "
                    "carry that meaning (see docs/metadata-schema.md)",
                    name,
                )
            if len(values) != n_channels:
                logger.info(
                    "Dropping mass axis annotation %r: %d values for a "
                    "%d channel axis (expected when resampling rebuilds it)",
                    name,
                    len(values),
                    n_channels,
                )
                continue
            validated[name] = values
        return validated

    def _create_pixel_shapes(self, adata: AnnData) -> "ShapesModel":
        """Create geometric shapes for pixels with proper transformations.

        When optical alignment is available (FlexImaging with Area definitions),
        shapes are created in optical image pixel coordinates for proper overlay.
        Otherwise, shapes use physical (micrometer) coordinates.

        **Footprints are two-dimensional, including on a multi-slice
        volume.** A slice's depth is carried by the TIC image's ``Scale``
        and by ``obs["spatial_z"]``; it is deliberately not also put on
        the polygon geometry.

        Three options were measured. All three are imperfect, so what
        follows is the reasoning rather than a claim that this one is
        free:

        * **``POLYGON Z``** -- geometrically the most honest, and what
          7317792 shipped in v3.2.0. It breaks ``spatialdata``'s spatial
          queries: measured on a 5x3x2 volume, a bounding box enclosing
          the whole dataset with 1000um of margin returned 26 of 30
          footprints and 26 of 30 table rows, with no exception and no
          warning. ``shapely.force_2d`` on the same geometry restored
          30 of 30. A z-restricted query returned the same rows whether
          z was inside or far outside the data, so z is ignored by the
          query path rather than honoured.
        * **Flat 2D shapes carrying the depth in a transformation**
          (one element per slice, each with a ``Translation`` in z).
          The element parses, coexists with a 3D image and round-trips,
          but the transform **silently drops z** and every slice comes
          back at the same depth -- a depth in the metadata that never
          reaches the geometry. ``test_3d_pixel_shapes_z.py`` pins this,
          because it is the change someone will otherwise propose.
        * **Flat, with the depth on the image and in ``obs`` only** --
          what this does.

        The deciding argument is that ``spatialdata`` itself asks for
        this. ``ShapesModel.validate`` warns that a 3-dimensional
        geometry column "could led to unexpected behaviors" and names
        ``force_2d()`` as the remedy; the query result above is that
        behaviour. 3D shapes are not on the upstream roadmap
        (scverse/spatialdata#109 has been idle since June 2023 and
        covers images, labels and transformations only), and the live
        2.5D discussion (#961) scopes itself to points, images and
        labels. Serial-section MSI is 2.5D in that taxonomy.

        What this costs: ``docs/coordinate-systems.md``'s promise that
        every element agrees at ``"global"`` is exact in x and y and
        silent in z for the shapes element. That is a documented gap
        rather than a wrong answer, which a truncated query is not. See
        ``docs/output-format.md`` for the consumer-facing statement.

        Reinstating z means restoring the ``is_3d`` parameter at all
        three call sites as well; it was removed with the geometry so it
        could not sit unread, which is how the original defect arose.

        Args:
            adata: AnnData object containing coordinates

        Returns:
            SpatialData shapes model
        """
        geometries = []

        # Track valid indices (for alignment mode where we skip empty positions)
        valid_indices: Optional[List[int]] = None

        # Use FlexImaging alignment for MSI shapes only when both
        # data is available AND the caller hasn't opted out via
        # apply_optical_alignment=False.  Opt-out leaves MSI in pure
        # micrometer coordinates so a downstream alignment step
        # (e.g. Ousia's EscDat registration) is the canonical mapping.
        #
        # The matrix is part of the condition, not just the alignment
        # result, so this agrees with the TIC image and with the
        # coordinate_systems attr. Gating on `optical.alignment` alone was
        # not equivalent: `build_tic_to_image_affine` returns early when
        # `region_mappings` is empty, which
        # `TeachingPointAlignment.compute_area_alignment` produces from a
        # real .mis whose areas match no region. In that
        # state the attr and the TIC image took the micrometer branch
        # while this took the alignment branch, `transform_point` returned
        # None for every position, and the shapes element came out with
        # zero polygons.
        use_msi_alignment = (
            self.optical.msi_in_pixel_space and self.optical.alignment is not None
        )

        if use_msi_alignment:
            # Use optical alignment - transform raster coords to image pixels
            # Only create shapes for positions that have actual spectra (in pos_to_region)
            raster_x: NDArray[np.int_] = adata.obs["x"].values
            raster_y: NDArray[np.int_] = adata.obs["y"].values

            # Default half-pixel size (fallback if region lookup fails)
            default_half_pixel = (10.0, 10.0)
            if self.optical.alignment.region_mappings:
                # Use first region as default
                rm = self.optical.alignment.region_mappings[0]
                default_half_pixel = rm.get_half_pixel_size()

            valid_indices = []
            for i in range(len(adata)):
                rx, ry = int(raster_x[i]), int(raster_y[i])
                img_coords = self.optical.alignment.transform_point(rx, ry)

                if img_coords is not None:
                    ix, iy = img_coords
                    # Get region-specific half-pixel size (may be non-square)
                    half_pixel = self.optical.alignment.get_half_pixel_size(rx, ry)
                    if half_pixel is None:
                        half_pixel = default_half_pixel

                    half_x, half_y = half_pixel
                    pixel_box = box(
                        ix - half_x,
                        iy - half_y,
                        ix + half_x,
                        iy + half_y,
                    )
                    geometries.append(pixel_box)
                    valid_indices.append(i)
                # Skip positions without spectra - empty grid cells

            n_skipped = len(adata) - len(valid_indices)
            if n_skipped > 0:
                logger.info(
                    f"Created {len(geometries)} shapes using optical alignment "
                    f"(skipped {n_skipped} empty grid positions)"
                )
            else:
                logger.info(
                    f"Created {len(geometries)} shapes using optical alignment "
                    f"(image pixel coordinates)"
                )
        else:
            # Standard physical coordinates (micrometers)
            from shapely import box as shapely_box_vectorized

            x_coords: NDArray[np.float64] = adata.obs["spatial_x"].values
            y_coords: NDArray[np.float64] = adata.obs["spatial_y"].values
            half_x_um = self.pixel_size_um / 2
            half_y_um = self.pixel_size_y_um / 2

            # Footprints are flat, on every route including volumes. A
            # slice's depth lives on the TIC image's Scale and in
            # obs["spatial_z"]; see the docstring for why it is not also
            # put on the geometry. Built through shapely's vectorised box
            # constructor -- one C call for the whole table instead of one
            # Python-level geometry per pixel.
            geometries = shapely_box_vectorized(
                x_coords - half_x_um,
                y_coords - half_y_um,
                x_coords + half_x_um,
                y_coords + half_y_um,
            )

        # Create GeoDataFrame with appropriate index
        if valid_indices is not None:
            # Use filtered indices for alignment mode
            filtered_index = adata.obs.index[valid_indices]
            gdf = gpd.GeoDataFrame(geometry=geometries, index=filtered_index)
        else:
            gdf = gpd.GeoDataFrame(geometry=geometries, index=adata.obs.index)

        # Set up transform
        transform = Identity()
        transformations = {self.dataset_id: transform, "global": transform}

        # Parse shapes
        shapes = ShapesModel.parse(gdf, transformations=transformations)
        return shapes

    def _save_output(self, state: ConversionState) -> bool:
        """Save the data to SpatialData format.

        Args:
            state: Data structures to save

        Returns:
            True if saving was successful, False otherwise
        """
        try:
            # Create SpatialData object with images included
            sdata = SpatialData(
                tables=state.tables,
                shapes=state.shapes,
                images=state.images,
            )

            # Add metadata
            self.add_metadata(sdata)

            # Write to disk. table_write_config() must wrap the write itself,
            # not just the SpatialData construction: anndata creates the table's
            # zarr arrays lazily inside sdata.write, and that is where the shard
            # budget is read. Without it the table lands in ~1 KB shards -- one
            # file per KB -- and the write dominates conversion wall-clock.
            with _suppress_upstream_warnings(), table_write_config():
                sdata.write(str(self.output_path))
                # The optical images above are placeholders; their pixels
                # stream into the store now that it exists.
                self.optical.stream_pending_pixels()
                zarr.consolidate_metadata(str(self.output_path))
            logger.info(f"Successfully saved SpatialData to {self.output_path}")
            # Every table was written from its memmaps; nothing holds them
            # now but this object and the caller's mapping.
            del sdata
            self.siblings.release_scratch(state.tables)
            return True
        except Exception as e:
            logger.error(f"Error saving SpatialData: {e}")
            import traceback

            logger.debug(f"Detailed traceback:\n{traceback.format_exc()}")
            return False

    def add_metadata(self, metadata: "SpatialData") -> None:
        """Add comprehensive metadata to the SpatialData object.

        Args:
            metadata: SpatialData object to add metadata to
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")

        # Call parent to prepare structured metadata
        super().add_metadata(metadata)

        # Get comprehensive metadata object for detailed access
        comprehensive_metadata_obj = self.reader.get_comprehensive_metadata()

        if not hasattr(metadata, "attrs") or metadata.attrs is None:
            metadata.attrs = {}

        logger.info("Adding comprehensive metadata to SpatialData.attrs")
        metadata.attrs.update(self.build_root_attrs(comprehensive_metadata_obj))

    def build_root_attrs(
        self, comprehensive_metadata_obj: Any = None
    ) -> Dict[str, Any]:
        """The store's root attrs, composed by :class:`RootAttrsBuilder`."""
        return self.root_attrs.build(
            self._root_attrs_context(), comprehensive_metadata_obj
        )

    def _root_attrs_context(self) -> RootAttrsContext:
        """What the builder needs that this conversion settled after setup.

        Packed fresh per call, never cached: the pitch is adopted from
        reader metadata, the dimensions and bounds arrive with the
        essential metadata, and the non-empty pixel count is only final
        once both passes have run.
        """
        if self._dimensions is None:
            raise RuntimeError("Dimensions are not initialized")
        return RootAttrsContext(
            optical=self.optical,
            pixel_size_xy=self._resolved_pixel_size_xy(),
            dimensions=self._dimensions,
            non_empty_pixels=self._non_empty_pixel_count,
            coordinate_bounds=self._coordinate_bounds,
            is_volume=self._is_volume,
            z_spacing_um=self.z_spacing_um,
            z_spacing_source=self.z_spacing_source,
        )

    @abstractmethod
    def _create_data_structures(self) -> ConversionState:
        """Create data structures for the specific converter type."""
        pass

    @abstractmethod
    def _finalize_data(self, state: ConversionState) -> None:
        """Finalize data structures for the specific converter type."""
        pass
