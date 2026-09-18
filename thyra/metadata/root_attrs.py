# thyra/metadata/root_attrs.py

"""The store's own attrs: what one conversion writes at the Zarr root.

Sibling of :mod:`thyra.metadata.uns_assembler`, and here for the same
reason. ``uns`` is the table's provenance; these are the store's, and the
two used to be composed by twenty-odd methods of one converter reading
its own attributes, so the only way to ask what a given input produces
was to run a conversion. A :class:`RootAttrsBuilder` is handed a reader,
the optical collaborator and a :class:`RootAttrsContext` instead.

**Root attrs are the only place provenance survives.** SpatialData drops
an element's ``.attrs`` on write -- only ``ome`` and ``spatialdata_attrs``
reach the store -- so anything that has to be readable back about the
conversion as a whole belongs here and nowhere else.

**Sections the reader has nothing for are omitted** rather than written
empty, matching the ``uns`` block, so a consumer can tell "not available
from this format" from "available and empty".
"""

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Mapping, Optional, Tuple

import pandas as pd

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..converters.spatialdata.optical_image import OpticalImages
    from ..core.base_reader import BaseMSIReader

logger = logging.getLogger(__name__)

#: Schema version for the structured ``coordinate_systems`` attr.
#: Bump when the schema shape changes in a way consumers need to notice.
COORDINATE_SYSTEMS_SCHEMA_VERSION: int = 1


@dataclass(frozen=True)
class RootAttrsContext:
    """What the converter knows when the root attrs are composed.

    Every field is state the converter settles after it is constructed:
    the pitch is adopted from the reader's metadata, the dimensions and
    the coordinate bounds arrive with the essential metadata, the
    non-empty pixel count is only final once both passes have run, and
    the z spacing is resolved against the detected pitch. Passing them
    per call is what keeps the builder from reaching back into the
    converter for any of it.

    Attributes:
        optical: The conversion's optical images, read for the
            ``optical_images`` attr and for whether the raster landed in
            optical-photo pixel space. Passed per call rather than held:
            the converter owns the collaborator and may replace it after
            the builder exists, and a captured reference then composes
            the coordinate contract from the wrong one.
        pixel_size_xy: The resolved in-plane pitch as ``(x_um, y_um)``.
            One pair, and every block of the store is written from it
            (issue #228).
        dimensions: The grid as ``(n_x, n_y, n_z)``.
        non_empty_pixels: Grid positions that carry a spectrum.
        coordinate_bounds: The source's own bounds, or ``None`` when it
            reported none.
        is_volume: Whether this conversion writes a real multi-slice
            volume -- ``handle_3d`` and a z axis to space out, not the
            flag alone.
        z_spacing_um: Slice pitch, meaningful only for a volume.
        z_spacing_source: Where that pitch came from; read as ``.value``.
    """

    optical: "OpticalImages"
    pixel_size_xy: Tuple[float, float]
    dimensions: Tuple[int, int, int]
    non_empty_pixels: int
    coordinate_bounds: Optional[Any]
    is_volume: bool
    z_spacing_um: float
    z_spacing_source: Any


class RootAttrsBuilder:
    """The root attrs of one conversion, composed from arguments.

    One entry point, :meth:`build`. The converter calls it once, when the
    store is about to be written, and copies the result onto
    ``SpatialData.attrs``.
    """

    def __init__(
        self,
        reader: "BaseMSIReader",
        *,
        dataset_id: str,
        pixel_size_detection_info: Optional[Dict[str, Any]],
        handle_3d: bool,
        conversion_options: Optional[Mapping[str, Any]] = None,
    ) -> None:
        """Build the root-attr builder of one conversion.

        Args:
            reader: The MSI reader. Only ``get_comprehensive_metadata``
                and ``get_essential_metadata`` are read.
            dataset_id: The store's dataset identifier.
            pixel_size_detection_info: What pitch detection reported, or
                ``None``.
            handle_3d: Whether the caller asked for a volume. Recorded as
                asked, which is why it is not
                ``RootAttrsContext.is_volume`` -- that one says whether a
                volume was actually written.
            conversion_options: Any further options the conversion was
                given, recorded verbatim.
        """
        self.reader = reader
        self._dataset_id = dataset_id
        self._pixel_size_detection_info = pixel_size_detection_info
        self._handle_3d = handle_3d
        self._conversion_options = dict(conversion_options or {})

    def build(
        self, ctx: RootAttrsContext, comprehensive_metadata_obj: Any = None
    ) -> Dict[str, Any]:
        """The root-level attributes every write path must persist, identically.

        The store's own attrs, as opposed to the table's ``uns`` block the
        assembler owns. Sibling of that method and here for the same
        reason: the streaming path used to hand-write its Zarr layout and
        composed its own, shorter set -- 7 attributes against 10 on real
        ``pea.imzML``, missing ``coordinate_systems``,
        ``format_specific_metadata`` and ``msi_dataset_info``.

        ``coordinate_systems`` is the one that matters most in practice:
        it is the structured contract saying what unit ``"global"`` is in,
        and Ousia and the registration tooling read it rather than
        guessing. A streaming store simply did not have it, and at the
        time the route was chosen by size, so the datasets that lost it
        were the largest ones. Every store is written through one path
        now, so this is the one place the attrs are composed.

        Args:
            ctx: What the converter knows about the store being written.
            comprehensive_metadata_obj: Already-read comprehensive
                metadata, when the caller has it. ``None`` reads it from
                the reader.

        Returns:
            Mapping of root attribute name to value.
        """
        if comprehensive_metadata_obj is None:
            comprehensive_metadata_obj = self.reader.get_comprehensive_metadata()

        attrs = self._pixel_size_attrs(ctx)
        self._add_comprehensive_sections(attrs, comprehensive_metadata_obj)
        return attrs

    def _pixel_size_attrs(self, ctx: RootAttrsContext) -> Dict[str, Any]:
        """Pixel size, conversion provenance and the coordinate contract."""
        try:
            from .. import __version__

            version = __version__
        except ImportError:
            version = "unknown"

        px, py = ctx.pixel_size_xy
        attrs: Dict[str, Any] = {
            "pixel_size_x_um": float(px),
            "pixel_size_y_um": float(py),
            "pixel_size_units": "micrometers",
            "coordinate_system": "physical_micrometers",
            "msi_converter_version": version,
            "conversion_timestamp": pd.Timestamp.now().isoformat(),
        }

        # Structured coordinate-system contract describing what "global"
        # actually means in this zarr. Consumers (e.g. Ousia, registration
        # tooling) read this to know what unit "global" is in and how to
        # convert to micrometers without guessing.
        attrs["coordinate_systems"] = self._coordinate_systems_attr(ctx, version)

        if self._pixel_size_detection_info is not None:
            attrs["pixel_size_detection_info"] = dict(self._pixel_size_detection_info)
            logger.info(
                "Added pixel size detection info: %s", self._pixel_size_detection_info
            )

        n_x, n_y, n_z = ctx.dimensions
        attrs["msi_dataset_info"] = {
            "dataset_id": self._dataset_id,
            "total_grid_pixels": n_x * n_y * n_z,
            "non_empty_pixels": ctx.non_empty_pixels,
            "dimensions_xyz": list(ctx.dimensions),
        }

        # The two sections worth keeping from the dead
        # _add_comprehensive_metadata (issue #67 item 1). That method
        # assembled a block and hung it on SpatialData.metadata, an
        # attribute no release of SpatialData has, so its body never ran
        # and nothing it built ever reached a store. Most of what it
        # assembled is already here -- the detection info as
        # pixel_size_detection_info, the pixel count in msi_dataset_info
        # -- and the element counts it took off the SpatialData object
        # describe a store the reader can already see. These two do not
        # survive anywhere else.
        if ctx.coordinate_bounds is not None:
            attrs["coordinate_bounds"] = list(ctx.coordinate_bounds)
        attrs["conversion_options"] = {
            "handle_3d": self._handle_3d,
            "pixel_size_um": float(px),
            "pixel_size_y_um": float(py),
            "z_spacing_um": ctx.z_spacing_um,
            "z_spacing_source": ctx.z_spacing_source.value,
            "dataset_id": self._dataset_id,
            **self._conversion_options,
        }

        optical = ctx.optical.root_attr()
        if optical is not None:
            attrs["optical_images"] = optical

        return attrs

    def _coordinate_systems_attr(
        self, ctx: RootAttrsContext, thyra_version: str
    ) -> Dict[str, Any]:
        """Build the structured coordinate-system contract attr.

        This describes what ``"global"`` means in the produced zarr so
        that downstream consumers can render and convert without
        guessing.

        Two variants are emitted depending on whether FlexImaging optical
        alignment was applied during conversion:

        - No alignment (``global = micrometer``): the TIC image carries a
          ``Scale(pixel_size_um)`` and pixel-polygon shapes are stored in
          micrometers with ``Identity``. Both elements agree at
          ``global``. ``pixel_size_um_x/y`` are filled with the MSI grid
          pixel size, since "global" is in physical micrometers and there
          is no canonical raster image other than the MSI itself.

        - With alignment (``global = pixel``): the TIC image carries an
          ``Affine`` mapping raster indices to optical-image pixels and
          shapes are stored directly in optical-image pixels with
          ``Identity``. The optical image is the canonical reference.
          ``pixel_size_um_x/y`` are typically unknown at conversion time
          (FlexImaging does not generally calibrate the optical photo to
          um); leave them null and let the consumer fill in.

        - With alignment data but ``apply_optical_alignment=False``: the
          micrometer variant, exactly as the no-alignment case. The
          matrix exists -- it is built whenever FlexImaging data is
          present, because the opt-out path needs its inverse to carry
          the optical photo into micrometers -- but it was not applied to
          the raster, so nothing in the store is in optical pixels. The
          condition is
          :attr:`~thyra.converters.spatialdata.optical_image.OpticalImages.msi_in_pixel_space`
          rather than the matrix alone for exactly this reason.

        Multi-slice volumes additionally get ``z_spacing_um`` and
        ``z_spacing_source``. These are written **only** for volumes, so a
        2D store is byte-identical to what earlier versions produced and
        their absence is itself the signal that no z axis exists. That is
        also why ``convention_version`` does not move: the keys are purely
        additive, a consumer that does not do 3D is unaffected, and
        bumping the version would make every existing consumer log a
        "newer than I understand" warning on ordinary 2D datasets.

        Note ``z_spacing_um`` is an absolute micrometre distance even when
        ``unit="pixel"``, because the 3D route always scales z by it
        directly -- the optical affine only ever governs x and y.

        Returns:
            Dict suitable for storing under
            ``zarr.attrs["coordinate_systems"]``.
        """
        px, py = ctx.pixel_size_xy
        if ctx.optical.msi_in_pixel_space:
            unit = "pixel"
            pixel_size_um_x: Optional[float] = None
            pixel_size_um_y: Optional[float] = None
            # The element, not the file. This used to be
            # `optical.primary_filename` -- the .mis <ImageFile> stem,
            # lowercased -- which names nothing in the store: the element
            # is <dataset_id>_optical_<suffix>, so a consumer following
            # the documented meaning ("the canonical raster element that
            # defines pixel space") got a key that never resolves. It is
            # None when this store does not hold the alignment image
            # (optical images not included, or its pixels could not be
            # read): "global" is still that image's pixel grid, there is
            # just no element here that is it. The filename is in
            # `optical_images`.
            reference_element: Optional[str] = ctx.optical.alignment_element
        else:
            unit = "micrometer"
            pixel_size_um_x = float(px)
            pixel_size_um_y = float(py)
            reference_element = None

        global_cs: Dict[str, Any] = {
            "unit": unit,
            "pixel_size_um_x": pixel_size_um_x,
            "pixel_size_um_y": pixel_size_um_y,
            "reference_element": reference_element,
            "convention_version": COORDINATE_SYSTEMS_SCHEMA_VERSION,
            "produced_by": f"thyra/{thyra_version}",
        }

        # Additive keys; `convention_version` stays 1 for the same
        # reason as z_spacing_um below.  `raster_to_global_affine` is
        # the explicit 3x3 row-major affine from TIC raster indices to
        # "global" (the same mapping the TIC element's transform
        # expresses), so a consumer that reads only attrs still gets
        # the full placement.  `coordinate_offsets_px` preserves the
        # source's raw acquisition-index offsets, which 0-based
        # normalisation otherwise erases; `stage_offset_um` is their
        # physical equivalent, written only when "global" is in
        # micrometers so it cannot be misread in the optical-pixel
        # variant.
        if ctx.optical.msi_in_pixel_space:
            global_cs["raster_to_global_affine"] = [
                [float(v) for v in row] for row in ctx.optical.tic_to_image
            ]
        else:
            global_cs["raster_to_global_affine"] = [
                [px, 0.0, 0.0],
                [0.0, py, 0.0],
                [0.0, 0.0, 1.0],
            ]
        offsets = self._source_coordinate_offsets()
        if offsets is not None:
            global_cs["coordinate_offsets_px"] = [int(v) for v in offsets]
            if unit == "micrometer":
                global_cs["stage_offset_um"] = [
                    float(offsets[0]) * float(px),
                    float(offsets[1]) * float(py),
                ]

        if ctx.is_volume:
            global_cs["z_spacing_um"] = float(ctx.z_spacing_um)
            global_cs["z_spacing_source"] = ctx.z_spacing_source.value

        return {"global": global_cs}

    def _source_coordinate_offsets(self) -> Optional[Tuple[int, int, int]]:
        """The reader's raw coordinate offsets, if it reported any."""
        try:
            essential = self.reader.get_essential_metadata()
        except Exception as e:
            logger.debug("Could not read coordinate offsets: %s", str(e))
            return None
        offsets = getattr(essential, "coordinate_offsets", None)
        if offsets is None:
            return None
        x, y, z = offsets
        return (int(x), int(y), int(z))

    @staticmethod
    def _add_comprehensive_sections(
        attrs: Dict[str, Any], comprehensive_metadata_obj: Any
    ) -> None:
        """Add comprehensive metadata sections to attributes."""
        if comprehensive_metadata_obj.format_specific:
            attrs["format_specific_metadata"] = (
                comprehensive_metadata_obj.format_specific
            )

        if comprehensive_metadata_obj.acquisition_params:
            attrs["acquisition_parameters"] = (
                comprehensive_metadata_obj.acquisition_params
            )

        if comprehensive_metadata_obj.instrument_info:
            attrs["instrument_information"] = comprehensive_metadata_obj.instrument_info
