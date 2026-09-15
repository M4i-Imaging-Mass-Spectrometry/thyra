"""Abstract base class for resampling strategies."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

import numpy as np
import numpy.typing as npt


@dataclass
class Spectrum:
    """Single mass spectrum with coordinates and metadata."""

    mz: npt.NDArray[np.floating[Any]]
    intensity: npt.NDArray[np.floating[Any]]
    coordinates: Tuple[int, int, int]  # (x, y, z) coordinates
    metadata: Optional[Dict[str, Any]] = None

    @property
    def is_centroid(self) -> bool:
        """Heuristic to detect centroid data."""
        if len(self.mz) < 100:
            return True
        # Check for zero-intensity gaps typical of centroid data
        zero_count = np.sum(self.intensity == 0)
        return bool(zero_count / len(self.intensity) > 0.5)


class ResamplingStrategy(ABC):
    """Abstract base class for resampling strategies."""

    @abstractmethod
    def resample(
        self, spectrum: Spectrum, target_axis: npt.NDArray[np.floating[Any]]
    ) -> Spectrum:
        """Resample spectrum to target mass axis.

        Args:
            spectrum: Input spectrum to resample
            target_axis: Target mass axis values

        Returns:
            Resampled spectrum with target_axis as mz values
        """
        pass
