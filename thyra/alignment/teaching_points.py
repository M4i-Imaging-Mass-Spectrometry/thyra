# thyra/alignment/teaching_points.py
"""Teaching point alignment for FlexImaging optical-MSI registration.

This module handles the alignment between optical images and MSI data
using teaching point calibration data from FlexImaging .mis files.

Coordinate Systems:
- Image pixels: (x, y) in the optical reference image (origin at top-left)
- Stage coordinates (teaching): From teaching point calibration
- Stage coordinates (poslog): From acquisition position log
- MSI raster: (x, y) grid positions in the MSI dataset (0-based)

Key Challenge:
FlexImaging uses different coordinate frames for teaching (image calibration)
and acquisition (stage movement). These frames may have large offsets that
cannot be reliably computed without additional reference data.

The module provides:
1. Reliable image <-> teaching stage transformation from teaching points
2. Estimated MSI <-> image transformation (may need manual verification)
3. Methods to manually specify alignment offsets if automatic alignment fails
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .affine import AffineTransform

logger = logging.getLogger(__name__)


@dataclass
class TeachingPoint:
    """A single teaching point correspondence.

    Attributes:
        image_x: X coordinate in image pixels
        image_y: Y coordinate in image pixels
        stage_x: X coordinate in stage units (micrometers)
        stage_y: Y coordinate in stage units (micrometers)
    """

    image_x: int
    image_y: int
    stage_x: int
    stage_y: int

    @classmethod
    def from_dict(cls, data: Dict[str, Tuple[int, int]]) -> "TeachingPoint":
        """Create from parsed dictionary format.

        Args:
            data: Dictionary with 'image' and 'stage' tuples

        Returns:
            TeachingPoint instance
        """
        img = data["image"]
        stage = data["stage"]
        return cls(
            image_x=img[0],
            image_y=img[1],
            stage_x=stage[0],
            stage_y=stage[1],
        )


@dataclass
class RasterPosition:
    """A position from the poslog with raster and physical coordinates.

    Attributes:
        raster_x: Raster grid X coordinate
        raster_y: Raster grid Y coordinate
        phys_x: Physical stage X coordinate (micrometers)
        phys_y: Physical stage Y coordinate (micrometers)
    """

    raster_x: int
    raster_y: int
    phys_x: float
    phys_y: float


@dataclass
class AlignmentResult:
    """Result of teaching point alignment computation.

    Attributes:
        image_to_stage: Affine transform from image pixels to stage coords
        stage_to_image: Inverse transform (stage to image pixels)
        msi_to_image: Transform from MSI raster to image pixels
        image_to_msi: Transform from image pixels to MSI raster
        stage_offset: Estimated offset between teaching and poslog stages
        rmse: Root mean square error of teaching point fit
        warnings: List of alignment warnings
    """

    image_to_stage: AffineTransform
    stage_to_image: AffineTransform
    msi_to_image: Optional[AffineTransform] = None
    image_to_msi: Optional[AffineTransform] = None
    stage_offset: Optional[Tuple[float, float]] = None
    rmse: float = 0.0
    warnings: List[str] = field(default_factory=list)


@dataclass
class RegionMapping:
    """Mapping for a single acquisition region to image coordinates.

    The image coordinates use min/max bounds from the Area definition.
    The p1/p2 points in an Area are two corners of a bounding box - their
    order indicates scan direction but for coordinate mapping we use the
    actual min/max values to ensure monotonic mapping:
    - raster_min_x maps to image_min_x
    - raster_max_x maps to image_max_x
    - raster_min_y maps to image_min_y
    - raster_max_y maps to image_max_y

    Attributes:
        region_id: Region number (0, 1, 2, ...)
        name: Region name from Area definition
        raster_min_x: Minimum raster X in this region
        raster_max_x: Maximum raster X in this region
        raster_min_y: Minimum raster Y in this region
        raster_max_y: Maximum raster Y in this region
        image_min_x: Minimum image X for this region
        image_max_x: Maximum image X for this region
        image_min_y: Minimum image Y for this region
        image_max_y: Maximum image Y for this region
    """

    region_id: int
    name: str
    raster_min_x: int
    raster_max_x: int
    raster_min_y: int
    raster_max_y: int
    image_min_x: int
    image_max_x: int
    image_min_y: int
    image_max_y: int

    def _get_pixel_scale(self) -> float:
        """Get the consistent pixel scale for this region.

        Uses scale_x as the reference since it's typically more consistent
        across regions. This ensures square pixels with no gaps.

        Returns:
            Pixel scale in image pixels per raster step
        """
        raster_width = max(1, self.raster_max_x - self.raster_min_x)
        image_width = self.image_max_x - self.image_min_x
        return image_width / raster_width

    def raster_to_image(self, raster_x: int, raster_y: int) -> Tuple[float, float]:
        """Convert raster coordinates to image coordinates within this region.

        Uses consistent spacing (scale_x) for both X and Y to ensure square
        pixels with no gaps. Pixels may not perfectly fill the bounding box
        but will be correctly positioned relative to each other.

        Args:
            raster_x: Original raster X coordinate (not normalized)
            raster_y: Original raster Y coordinate (not normalized)

        Returns:
            Tuple of (image_x, image_y) in pixels
        """
        # Use consistent scale for both dimensions (square pixels, no gaps)
        scale = self._get_pixel_scale()

        # Position based on offset from first raster position
        image_x = self.image_min_x + (raster_x - self.raster_min_x) * scale
        image_y = self.image_min_y + (raster_y - self.raster_min_y) * scale

        return image_x, image_y

    def get_half_pixel_size(self) -> float:
        """Get the half-pixel size for this region in image coordinates.

        Returns a single value for square pixels. Uses the same scale as
        raster_to_image() to ensure pixel size matches spacing (no gaps).

        Returns:
            Half-pixel size in image pixels
        """
        return self._get_pixel_scale() / 2


@dataclass
class AreaAlignmentResult:
    """Result of area-based alignment computation.

    This is the preferred alignment method when Area definitions are available
    in the .mis file, as it provides direct region-to-image mapping.

    Attributes:
        region_mappings: List of RegionMapping objects for each region
        first_raster_x: First raster X offset from header
        first_raster_y: First raster Y offset from header
        pos_to_region: Mapping from (raster_x, raster_y) to region_id
    """

    region_mappings: List[RegionMapping]
    first_raster_x: int
    first_raster_y: int
    pos_to_region: Dict[Tuple[int, int], int] = field(default_factory=dict)

    def transform_point(
        self, norm_x: int, norm_y: int
    ) -> Optional[Tuple[float, float]]:
        """Transform normalized raster coordinates to image coordinates.

        Args:
            norm_x: Normalized (0-based) raster X coordinate
            norm_y: Normalized (0-based) raster Y coordinate

        Returns:
            Tuple of (image_x, image_y) or None if no mapping exists
        """
        # Convert normalized to original raster coords
        orig_x = norm_x + self.first_raster_x
        orig_y = norm_y + self.first_raster_y

        # Find which region this belongs to
        region_id = self.pos_to_region.get((orig_x, orig_y))
        if region_id is None:
            return None

        # Find the region mapping
        for mapping in self.region_mappings:
            if mapping.region_id == region_id:
                return mapping.raster_to_image(orig_x, orig_y)

        return None

    def get_half_pixel_size(self, norm_x: int, norm_y: int) -> Optional[float]:
        """Get the half-pixel size for a given position.

        Returns a single value for square pixels. The scale may vary
        between regions due to Area definitions.

        Args:
            norm_x: Normalized (0-based) raster X coordinate
            norm_y: Normalized (0-based) raster Y coordinate

        Returns:
            Half-pixel size in image pixels, or None if no mapping exists
        """
        # Convert normalized to original raster coords
        orig_x = norm_x + self.first_raster_x
        orig_y = norm_y + self.first_raster_y

        # Find which region this belongs to
        region_id = self.pos_to_region.get((orig_x, orig_y))
        if region_id is None:
            return None

        # Find the region mapping and get its half-pixel size
        for mapping in self.region_mappings:
            if mapping.region_id == region_id:
                return mapping.get_half_pixel_size()

        return None


class TeachingPointAlignment:
    """Computes alignment between optical images and MSI data.

    This class handles the coordinate system transformations needed to
    align optical images with MSI raster data using FlexImaging teaching
    points.

    The workflow:
    1. Parse teaching points from .mis metadata
    2. Compute image -> stage affine transformation
    3. Determine offset between teaching stage and poslog stage coords
    4. Compute final image -> MSI raster transformation

    Example:
        >>> aligner = TeachingPointAlignment()
        >>> result = aligner.compute_alignment(
        ...     teaching_points=reader.mis_metadata['teaching_points'],
        ...     poslog_positions=reader._positions,
        ...     raster_step=(20.0, 20.0),
        ... )
        >>> # Transform optical image coordinates to MSI raster
        >>> msi_coords = result.image_to_msi.transform_point(img_x, img_y)
    """

    def __init__(self):
        """Initialize the alignment calculator."""
        pass

    def compute_alignment(
        self,
        teaching_points: List[Dict[str, Tuple[int, int]]],
        poslog_positions: Optional[List[Dict[str, Any]]] = None,
        raster_step: Tuple[float, float] = (20.0, 20.0),
        raster_offset: Tuple[int, int] = (0, 0),
        flip_poslog_x: bool = False,
        flip_poslog_y: bool = False,
    ) -> AlignmentResult:
        """Compute alignment transformations from teaching points.

        Args:
            teaching_points: List of teaching point dictionaries with
                'image' and 'stage' keys containing (x, y) tuples
            poslog_positions: Optional list of position dictionaries from
                poslog parsing, used to estimate stage coordinate offset
            raster_step: (step_x, step_y) raster step size in micrometers
            raster_offset: (offset_x, offset_y) offset of first raster position
            flip_poslog_x: If True, negate poslog X coordinates (for inverted
                stage X-axis relative to teaching points)
            flip_poslog_y: If True, negate poslog Y coordinates (for inverted
                stage Y-axis relative to teaching points)

        Returns:
            AlignmentResult with computed transformations
        """
        warnings: List[str] = []

        # Parse and validate teaching points
        if len(teaching_points) < 3:
            n_pts = len(teaching_points)
            raise ValueError(f"At least 3 teaching points required, got {n_pts}")

        points = [TeachingPoint.from_dict(tp) for tp in teaching_points]
        logger.info(f"Processing {len(points)} teaching points")

        # Compute transformations and RMSE
        image_to_stage = self._compute_image_to_stage(points)
        rmse = self._compute_alignment_rmse(points, image_to_stage)

        if rmse > 10.0:  # More than 10 um error
            warnings.append(
                f"Teaching point fit has high error (RMSE={rmse:.2f} um). "
                "Alignment may be inaccurate."
            )

        stage_to_image = image_to_stage.inverse()

        # Compute MSI transforms if poslog available
        msi_result = self._compute_msi_alignment(
            points,
            poslog_positions,
            image_to_stage,
            stage_to_image,
            raster_step,
            flip_poslog_x,
            flip_poslog_y,
        )

        if msi_result["warning"]:
            warnings.append(msi_result["warning"])

        return AlignmentResult(
            image_to_stage=image_to_stage,
            stage_to_image=stage_to_image,
            msi_to_image=msi_result["msi_to_image"],
            image_to_msi=msi_result["image_to_msi"],
            stage_offset=msi_result["stage_offset"],
            rmse=rmse,
            warnings=warnings,
        )

    def _compute_alignment_rmse(
        self, points: List[TeachingPoint], transform: AffineTransform
    ) -> float:
        """Compute RMSE for alignment quality assessment."""
        src_pts = [(float(p.image_x), float(p.image_y)) for p in points]
        dst_pts = [(float(p.stage_x), float(p.stage_y)) for p in points]
        rmse, _ = transform.compute_residuals(src_pts, dst_pts)
        logger.info(f"Image->Stage RMSE: {rmse:.4f} um")
        return rmse

    def _apply_coordinate_flips(
        self,
        positions: List[Dict[str, Any]],
        flip_x: bool,
        flip_y: bool,
    ) -> List[Dict[str, Any]]:
        """Apply coordinate flips to poslog positions if needed."""
        if not flip_x and not flip_y:
            return positions

        flipped = []
        for pos in positions:
            flipped_pos = pos.copy()
            if flip_x:
                flipped_pos["phys_x"] = -pos["phys_x"]
            if flip_y:
                flipped_pos["phys_y"] = -pos["phys_y"]
            flipped.append(flipped_pos)
        return flipped

    def _compute_msi_alignment(
        self,
        points: List[TeachingPoint],
        poslog_positions: Optional[List[Dict[str, Any]]],
        image_to_stage: AffineTransform,
        stage_to_image: AffineTransform,
        raster_step: Tuple[float, float],
        flip_poslog_x: bool,
        flip_poslog_y: bool,
    ) -> Dict[str, Any]:
        """Compute MSI-to-image alignment from poslog positions."""
        result: Dict[str, Any] = {
            "stage_offset": None,
            "msi_to_image": None,
            "image_to_msi": None,
            "warning": None,
        }

        if not poslog_positions:
            return result

        flipped_positions = self._apply_coordinate_flips(
            poslog_positions, flip_poslog_x, flip_poslog_y
        )

        stage_offset = self._estimate_stage_offset(
            points, flipped_positions, raster_step
        )
        result["stage_offset"] = stage_offset

        if stage_offset is None:
            result["warning"] = (
                "Could not determine stage coordinate offset. "
                "MSI-to-image transform may require manual calibration."
            )
            return result

        logger.info(
            f"Estimated stage offset: ({stage_offset[0]:.1f}, "
            f"{stage_offset[1]:.1f}) um"
        )

        first_phys = (
            float(flipped_positions[0]["phys_x"]),
            float(flipped_positions[0]["phys_y"]),
        )

        msi_to_image, image_to_msi = self._compute_msi_transforms(
            image_to_stage,
            stage_to_image,
            stage_offset,
            raster_step,
            first_phys,
            flip_poslog_x,
        )
        result["msi_to_image"] = msi_to_image
        result["image_to_msi"] = image_to_msi

        return result

    def _compute_image_to_stage(self, points: List[TeachingPoint]) -> AffineTransform:
        """Compute affine transformation from image pixels to stage coords.

        Args:
            points: List of teaching points

        Returns:
            AffineTransform from image to stage coordinates
        """
        src_points = [(float(p.image_x), float(p.image_y)) for p in points]
        dst_points = [(float(p.stage_x), float(p.stage_y)) for p in points]

        return AffineTransform.from_points(src_points, dst_points)

    def _estimate_stage_offset(
        self,
        teaching_points: List[TeachingPoint],
        poslog_positions: List[Dict[str, Any]],
        raster_step: Tuple[float, float],
    ) -> Optional[Tuple[float, float]]:
        """Estimate offset between teaching stage and poslog stage coords.

        This computes the translation offset between the two coordinate
        systems by comparing their centers. The assumption is that the
        MSI acquisition region is approximately centered on the optical
        image region defined by the teaching points.

        Coordinate relationship:
        - teaching_stage = poslog_physical - offset
        - poslog_physical = teaching_stage + offset

        Args:
            teaching_points: Teaching point data
            poslog_positions: Position log entries
            raster_step: Raster step size (step_x, step_y) in um

        Returns:
            Estimated (offset_x, offset_y) or None if cannot determine
        """
        if not poslog_positions:
            return None

        # Extract poslog physical coordinates
        phys_x = np.array([pos["phys_x"] for pos in poslog_positions])
        phys_y = np.array([pos["phys_y"] for pos in poslog_positions])

        # Get teaching point stage coordinate center
        teaching_x = [p.stage_x for p in teaching_points]
        teaching_y = [p.stage_y for p in teaching_points]
        teaching_center_x = float(np.mean(teaching_x))
        teaching_center_y = float(np.mean(teaching_y))

        # Get poslog physical coordinate center
        poslog_center_x = float(np.mean(phys_x))
        poslog_center_y = float(np.mean(phys_y))

        logger.debug(
            f"Teaching stage center: "
            f"({teaching_center_x:.1f}, {teaching_center_y:.1f})"
        )
        logger.debug(
            f"Poslog physical center: "
            f"({poslog_center_x:.1f}, {poslog_center_y:.1f})"
        )

        # The offset is the difference between poslog physical coords
        # and teaching stage coords for the same physical location
        # offset = poslog_physical - teaching_stage
        offset_x = poslog_center_x - teaching_center_x
        offset_y = poslog_center_y - teaching_center_y

        logger.info(
            f"Stage coordinate offset (poslog - teaching): "
            f"({offset_x:.1f}, {offset_y:.1f}) um"
        )

        return (offset_x, offset_y)

    def _compute_msi_transforms(
        self,
        image_to_stage: AffineTransform,
        stage_to_image: AffineTransform,
        stage_offset: Tuple[float, float],
        raster_step: Tuple[float, float],
        first_phys: Tuple[float, float],
        flip_poslog_x: bool = False,
    ) -> Tuple[AffineTransform, AffineTransform]:
        """Compute transformations between MSI raster and image coordinates.

        The chain of transformations:
        MSI raster (0-based) -> poslog physical -> teaching stage -> image

        The poslog Y-axis is inverted: larger raster Y maps to smaller
        physical Y (step_y is effectively negative).

        Args:
            image_to_stage: Transform from image pixels to teaching stage
            stage_to_image: Inverse transform
            stage_offset: (offset_x, offset_y) where offset = poslog - teaching
            raster_step: (step_x, step_y) raster step in um (magnitudes)
            first_phys: (phys_x, phys_y) physical position at raster (0, 0)
            flip_poslog_x: If True, negate X scale (poslog X is inverted)

        Returns:
            Tuple of (msi_to_image, image_to_msi) transforms
        """
        step_x, step_y = raster_step
        offset_x, offset_y = stage_offset
        first_phys_x, first_phys_y = first_phys

        # Determine X scale sign based on flip
        # If flip_poslog_x is True, the poslog coordinates were negated,
        # so the raster step direction is also reversed
        x_scale = -step_x if flip_poslog_x else step_x

        # Normalized raster (0-based) to teaching stage coordinates:
        # 1. normalized_raster -> physical:
        #    phys_x = x_scale * norm_x + first_phys_x
        #    phys_y = -step_y * norm_y + first_phys_y  (Y is inverted!)
        # 2. physical -> teaching:
        #    teaching = physical - offset
        #
        # Combined: normalized raster -> teaching stage
        #    teaching_x = x_scale * norm_x + first_phys_x - offset_x
        #    teaching_y = -step_y * norm_y + first_phys_y - offset_y

        tx = first_phys_x - offset_x
        ty = first_phys_y - offset_y

        logger.debug(
            f"MSI transform: scale=({x_scale}, {-step_y}), "
            f"translation=({tx:.1f}, {ty:.1f})"
        )

        # Build MSI -> teaching stage transform
        # Note: step_y is negated because Y-axis is inverted in poslog
        msi_to_teaching = AffineTransform.from_scale_translate(
            scale_x=x_scale,
            scale_y=-step_y,  # Negative because poslog Y is inverted
            tx=tx,
            ty=ty,
        )

        # Chain: MSI -> teaching stage -> image
        msi_to_image = msi_to_teaching.compose(stage_to_image)

        # Inverse: image -> MSI
        image_to_msi = msi_to_image.inverse()

        return msi_to_image, image_to_msi

    def validate_alignment(
        self,
        result: AlignmentResult,
        image_shape: Tuple[int, int],
        msi_shape: Tuple[int, int],
    ) -> List[str]:
        """Validate alignment by checking if coordinates map sensibly.

        Args:
            result: AlignmentResult to validate
            image_shape: (height, width) of optical image
            msi_shape: (height, width) of MSI raster

        Returns:
            List of validation warnings (empty if OK)
        """
        warnings = []
        img_h, img_w = image_shape
        msi_h, msi_w = msi_shape

        if result.msi_to_image is None:
            warnings.append("MSI-to-image transform not available")
            return warnings

        # Check corners of MSI raster map to within image bounds
        corners = [
            (0, 0),
            (msi_w - 1, 0),
            (0, msi_h - 1),
            (msi_w - 1, msi_h - 1),
        ]

        for cx, cy in corners:
            ix, iy = result.msi_to_image.transform_point(cx, cy)

            # Allow some margin outside image bounds
            margin = 0.2  # 20% margin
            if ix < -img_w * margin or ix > img_w * (1 + margin):
                warnings.append(
                    f"MSI corner ({cx}, {cy}) maps to image X={ix:.0f}, "
                    f"outside valid range [0, {img_w}]"
                )
            if iy < -img_h * margin or iy > img_h * (1 + margin):
                warnings.append(
                    f"MSI corner ({cx}, {cy}) maps to image Y={iy:.0f}, "
                    f"outside valid range [0, {img_h}]"
                )

        return warnings

    def compute_area_alignment(
        self,
        areas: List[Dict[str, Any]],
        poslog_positions: List[Dict[str, Any]],
        first_raster_x: int,
        first_raster_y: int,
    ) -> AreaAlignmentResult:
        """Compute alignment using Area definitions from .mis file.

        This is the preferred alignment method when Area definitions are
        available, as they provide direct image pixel coordinates for each
        acquisition region without requiring coordinate system transformations.

        Args:
            areas: List of Area dictionaries with 'name', 'p1', 'p2' keys
                where p1 and p2 are (x, y) tuples of image pixel coordinates
            poslog_positions: List of position dictionaries from poslog parsing
            first_raster_x: First raster X offset from header
            first_raster_y: First raster Y offset from header

        Returns:
            AreaAlignmentResult with region mappings and coordinate transform
        """
        # Group positions by region
        regions: Dict[int, Dict[str, Any]] = {}
        for pos in poslog_positions:
            r = pos["region"]
            if r not in regions:
                regions[r] = {
                    "positions": [],
                    "raster_xs": [],
                    "raster_ys": [],
                }
            regions[r]["positions"].append(pos)
            regions[r]["raster_xs"].append(pos["raster_x"])
            regions[r]["raster_ys"].append(pos["raster_y"])

        # Compute bounds for each region
        for r in regions:
            rx = regions[r]["raster_xs"]
            ry = regions[r]["raster_ys"]
            regions[r]["min_x"] = min(rx)
            regions[r]["max_x"] = max(rx)
            regions[r]["min_y"] = min(ry)
            regions[r]["max_y"] = max(ry)

        # Build position to region lookup
        pos_to_region: Dict[Tuple[int, int], int] = {}
        for r, data in regions.items():
            for pos in data["positions"]:
                pos_to_region[(pos["raster_x"], pos["raster_y"])] = r

        # Match areas to regions by index (Area 01 = Region 0, etc.)
        region_mappings: List[RegionMapping] = []
        for i, area in enumerate(areas):
            if i not in regions:
                logger.warning(f"Area '{area['name']}' has no matching region {i}")
                continue

            # Area bounds in image pixels - use min/max of p1 and p2
            # p1 and p2 are two corners of the bounding box, their order
            # indicates scan direction but we use min/max for mapping
            p1_x, p1_y = area["p1"]
            p2_x, p2_y = area["p2"]

            mapping = RegionMapping(
                region_id=i,
                name=area["name"],
                raster_min_x=regions[i]["min_x"],
                raster_max_x=regions[i]["max_x"],
                raster_min_y=regions[i]["min_y"],
                raster_max_y=regions[i]["max_y"],
                image_min_x=min(p1_x, p2_x),
                image_max_x=max(p1_x, p2_x),
                image_min_y=min(p1_y, p2_y),
                image_max_y=max(p1_y, p2_y),
            )
            region_mappings.append(mapping)

            logger.info(
                f"Region {i} -> Area '{area['name']}': "
                f"raster ({mapping.raster_min_x}, {mapping.raster_min_y}) to "
                f"({mapping.raster_max_x}, {mapping.raster_max_y}), "
                f"image ({mapping.image_min_x}, {mapping.image_min_y}) to "
                f"({mapping.image_max_x}, {mapping.image_max_y})"
            )

        return AreaAlignmentResult(
            region_mappings=region_mappings,
            first_raster_x=first_raster_x,
            first_raster_y=first_raster_y,
            pos_to_region=pos_to_region,
        )
