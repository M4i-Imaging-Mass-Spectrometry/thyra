"""The ``"acquisition"`` coordinate system: the raster where the source put it.

Readers re-zero the raster: grid index 0 is the first position the source
acquired, wherever on the slide that was. ``"global"`` and the per-dataset
system both start there, so two stores cut from one slide (two ``--region``
runs of one ``.d``, or two imzML files exported from parts of one
acquisition) land on top of each other.

``"acquisition"`` puts the offset back. It is ``"global"`` followed by a
translation to the source's own origin, in micrometres (D48):

* a source that numbers its raster positions: the index its own numbering
  gives grid position 0 (``coordinate_offsets``) times the pitch (Bruker,
  imzML, mzPeak, solariX, rapifleX, PHI);
* Waters, which re-zeroes on a fitted stage origin: that origin.

It is written only on stores whose ``"global"`` is in micrometres. A store
placed in a FlexImaging optical scan's pixels (D31) already shares the
scan's frame, and there the affine absorbs the offset. It is left out, too,
when the source reported no origin, because a zero there would claim one.
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional, Sequence, Tuple

from spatialdata.transformations import Identity
from spatialdata.transformations import Sequence as SequenceTransform
from spatialdata.transformations import Translation

#: The name of the coordinate system. Reserved: a ``dataset_id`` may not
#: take it, or the per-dataset system would silently become this one.
ACQUISITION_COORDINATE_SYSTEM = "acquisition"

#: The origin is a raster index offset times the pitch.
ORIGIN_FROM_RASTER_INDEX = "raster_index"
#: The origin is a fitted stage position, in micrometres (Waters).
ORIGIN_FROM_STAGE_POSITION = "stage_position"


@dataclass(frozen=True)
class AcquisitionOrigin:
    """Where grid index 0 sits in the source's own frame, in micrometres.

    Attributes:
        x_um: The x position of the centre of column 0.
        y_um: The y position of the centre of row 0.
        z_um: The z position of slice 0. Used only on a volume.
        source: :data:`ORIGIN_FROM_RASTER_INDEX` or
            :data:`ORIGIN_FROM_STAGE_POSITION`.
    """

    x_um: float
    y_um: float
    z_um: float
    source: str

    def along(self, axes: Sequence[str]) -> Tuple[float, ...]:
        """The origin's components on ``axes``, in that order."""
        values = {"x": self.x_um, "y": self.y_um, "z": self.z_um}
        return tuple(values[axis] for axis in axes)


def acquisition_origin(
    coordinate_offsets: Optional[Tuple[int, int, int]],
    stage_origin_um: Optional[Tuple[float, float]],
    pixel_size_xy: Tuple[float, float],
    z_spacing_um: Optional[float],
) -> Optional[AcquisitionOrigin]:
    """The origin a reader's metadata states, or ``None`` when it states none.

    A fitted stage origin wins over index offsets: a reader that has one
    reports it because its indices are counted from it.

    Args:
        coordinate_offsets: The index the source's numbering gives grid
            position 0, as ``EssentialMetadata.coordinate_offsets`` reports
            it.
        stage_origin_um: The fitted stage origin, as
            ``EssentialMetadata.stage_origin_um`` reports it.
        pixel_size_xy: The store's in-plane pitch as ``(x_um, y_um)``.
        z_spacing_um: The slice pitch; ``None`` leaves z at 0.
    """
    if stage_origin_um is not None:
        x, y = stage_origin_um
        return AcquisitionOrigin(
            float(x), float(y), 0.0, source=ORIGIN_FROM_STAGE_POSITION
        )
    if coordinate_offsets is None:
        return None
    px, py = pixel_size_xy
    ox, oy, oz = coordinate_offsets
    z = float(oz) * float(z_spacing_um) if z_spacing_um else 0.0
    return AcquisitionOrigin(
        float(ox) * float(px),
        float(oy) * float(py),
        z,
        source=ORIGIN_FROM_RASTER_INDEX,
    )


def element_transformations(
    dataset_id: str,
    transform: Any,
    origin: Optional[AcquisitionOrigin],
    axes: Tuple[str, ...] = ("x", "y"),
) -> Dict[str, Any]:
    """An element's transformations: its own two, plus ``"acquisition"``.

    ``transform`` places the element in micrometres; the acquisition
    system appends the origin's translation on ``axes``. A 2D element
    takes x and y only: a flat element's z translation is dropped by
    spatialdata anyway (see ``_create_pixel_shapes``).

    Args:
        dataset_id: Names the per-dataset system.
        transform: The element's transform to ``"global"``.
        origin: Where index 0 sits, or ``None`` to write no acquisition
            system (an unknown origin, or a store in optical pixels).
        axes: The axes the translation runs along.
    """
    transformations = {dataset_id: transform, "global": transform}
    if origin is None:
        return transformations
    shift = Translation(list(origin.along(axes)), axes=axes)
    transformations[ACQUISITION_COORDINATE_SYSTEM] = (
        shift
        if isinstance(transform, Identity)
        else SequenceTransform([transform, shift])
    )
    return transformations
