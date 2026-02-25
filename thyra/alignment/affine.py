# thyra/alignment/affine.py
"""Affine transformation utilities for coordinate system alignment."""

import logging
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


@dataclass
class AffineTransform:
    """Represents a 2D affine transformation.

    The transformation is defined as:
        x' = a * x + b * y + tx
        y' = c * x + d * y + ty

    Or in matrix form:
        [x']   [a  b  tx] [x]
        [y'] = [c  d  ty] [y]
        [1 ]   [0  0  1 ] [1]

    Attributes:
        matrix: 3x3 transformation matrix
        scale_x: Scale factor in X direction
        scale_y: Scale factor in Y direction
        rotation: Rotation angle in degrees
        shear: Shear factor
        translation: (tx, ty) translation offset
    """

    matrix: NDArray[np.float64]

    @property
    def scale_x(self) -> float:
        """Get scale factor in X direction."""
        return float(np.sqrt(self.matrix[0, 0] ** 2 + self.matrix[1, 0] ** 2))

    @property
    def scale_y(self) -> float:
        """Get scale factor in Y direction."""
        return float(np.sqrt(self.matrix[0, 1] ** 2 + self.matrix[1, 1] ** 2))

    @property
    def rotation(self) -> float:
        """Get rotation angle in degrees."""
        angle = np.arctan2(self.matrix[1, 0], self.matrix[0, 0])
        return float(np.degrees(angle))

    @property
    def translation(self) -> Tuple[float, float]:
        """Get translation offset (tx, ty)."""
        return (float(self.matrix[0, 2]), float(self.matrix[1, 2]))

    @classmethod
    def from_points(
        cls,
        src_points: List[Tuple[float, float]],
        dst_points: List[Tuple[float, float]],
    ) -> "AffineTransform":
        """Compute affine transformation from point correspondences.

        Uses least squares to find the best-fit affine transformation
        that maps source points to destination points.

        Args:
            src_points: List of (x, y) source coordinates
            dst_points: List of (x, y) destination coordinates

        Returns:
            AffineTransform instance

        Raises:
            ValueError: If fewer than 3 point pairs provided
        """
        if len(src_points) < 3 or len(dst_points) < 3:
            raise ValueError("At least 3 point pairs required for affine")

        if len(src_points) != len(dst_points):
            raise ValueError("Source and destination point counts must match")

        n = len(src_points)

        # Build the system of equations: A * params = b
        # For each point: x' = a*x + b*y + tx, y' = c*x + d*y + ty
        # Rearranged: [x y 1 0 0 0] [a]   [x']
        #             [0 0 0 x y 1] [b] = [y']
        #                           [tx]
        #                           [c]
        #                           [d]
        #                           [ty]

        A = np.zeros((2 * n, 6))
        b = np.zeros(2 * n)

        for i, (src, dst) in enumerate(zip(src_points, dst_points)):
            sx, sy = src
            dx, dy = dst

            # Equation for x'
            A[2 * i, 0] = sx
            A[2 * i, 1] = sy
            A[2 * i, 2] = 1
            b[2 * i] = dx

            # Equation for y'
            A[2 * i + 1, 3] = sx
            A[2 * i + 1, 4] = sy
            A[2 * i + 1, 5] = 1
            b[2 * i + 1] = dy

        # Solve using least squares
        params, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)

        # Build 3x3 matrix
        matrix = np.array(
            [
                [params[0], params[1], params[2]],
                [params[3], params[4], params[5]],
                [0, 0, 1],
            ],
            dtype=np.float64,
        )

        transform = cls(matrix=matrix)

        # Log transformation properties
        logger.info(
            f"Computed affine transform: "
            f"scale=({transform.scale_x:.4f}, {transform.scale_y:.4f}), "
            f"rotation={transform.rotation:.2f} deg"
        )

        return transform

    @classmethod
    def identity(cls) -> "AffineTransform":
        """Create an identity transformation."""
        return cls(matrix=np.eye(3, dtype=np.float64))

    @classmethod
    def from_scale_translate(
        cls,
        scale_x: float,
        scale_y: float,
        tx: float = 0.0,
        ty: float = 0.0,
    ) -> "AffineTransform":
        """Create transformation from scale and translation parameters.

        Args:
            scale_x: Scale factor in X direction
            scale_y: Scale factor in Y direction
            tx: Translation in X
            ty: Translation in Y

        Returns:
            AffineTransform instance
        """
        matrix = np.array(
            [[scale_x, 0, tx], [0, scale_y, ty], [0, 0, 1]], dtype=np.float64
        )
        return cls(matrix=matrix)

    def transform_point(self, x: float, y: float) -> Tuple[float, float]:
        """Transform a single point.

        Args:
            x: X coordinate
            y: Y coordinate

        Returns:
            Transformed (x', y') coordinates
        """
        pt = np.array([x, y, 1])
        result = self.matrix @ pt
        return (float(result[0]), float(result[1]))

    def transform_points(self, points: NDArray[np.float64]) -> NDArray[np.float64]:
        """Transform multiple points.

        Args:
            points: Nx2 array of (x, y) coordinates

        Returns:
            Nx2 array of transformed coordinates
        """
        n = points.shape[0]
        # Add homogeneous coordinate
        homogeneous = np.hstack([points, np.ones((n, 1))])
        # Transform
        transformed = (self.matrix @ homogeneous.T).T
        return transformed[:, :2]

    def inverse(self) -> "AffineTransform":
        """Compute inverse transformation.

        Returns:
            Inverse AffineTransform
        """
        inv_matrix = np.linalg.inv(self.matrix).astype(np.float64)
        return AffineTransform(matrix=inv_matrix)

    def compose(self, other: "AffineTransform") -> "AffineTransform":
        """Compose with another transformation.

        The resulting transformation applies self first, then other.
        T_result = T_other * T_self

        Args:
            other: Transformation to apply after this one

        Returns:
            Composed AffineTransform
        """
        composed = other.matrix @ self.matrix
        return AffineTransform(matrix=composed)

    def to_spatialdata_matrix(self) -> NDArray[np.float64]:
        """Convert to SpatialData-compatible transformation matrix.

        SpatialData uses a 3x3 matrix in the format:
            [[scale_x, shear_x, offset_x],
             [shear_y, scale_y, offset_y],
             [0,       0,       1       ]]

        Returns:
            3x3 numpy array
        """
        return self.matrix.copy()

    def compute_residuals(
        self,
        src_points: List[Tuple[float, float]],
        dst_points: List[Tuple[float, float]],
    ) -> Tuple[float, NDArray[np.float64]]:
        """Compute transformation residuals for validation.

        Args:
            src_points: Source points
            dst_points: Expected destination points

        Returns:
            Tuple of (RMSE, array of individual residuals)
        """
        residuals = []
        for src, dst in zip(src_points, dst_points):
            transformed = self.transform_point(*src)
            residual = np.sqrt(
                (transformed[0] - dst[0]) ** 2 + (transformed[1] - dst[1]) ** 2
            )
            residuals.append(residual)

        residuals_array = np.array(residuals)
        rmse = float(np.sqrt(np.mean(residuals_array**2)))

        return rmse, residuals_array

    def __repr__(self) -> str:
        """String representation."""
        sx, sy = self.scale_x, self.scale_y
        rot = self.rotation
        trans = self.translation
        return (
            f"AffineTransform(scale=({sx:.4f}, {sy:.4f}), "
            f"rotation={rot:.2f} deg, translation={trans})"
        )
