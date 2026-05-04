# thyra/metadata/extractors/waters_extractor.py
"""Waters-specific metadata extractor for MSI data.

Extracts essential and comprehensive metadata from Waters .raw imaging
files using the MassLynx native library and pre-built imaging grid.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from ...core.base_extractor import MetadataExtractor
from ..types import ComprehensiveMetadata, EssentialMetadata

logger = logging.getLogger(__name__)


class WatersMetadataExtractor(MetadataExtractor):
    """Waters-specific metadata extractor.

    Uses the pre-built ImagingGrid for spatial dimensions and pixel sizes,
    and scans all MS spectra via the native library for mass range and
    peak count information.
    """

    def __init__(
        self,
        ml,  # MassLynxLib instance
        handle,  # opaque file handle
        data_path: Path,
        imaging_grid,  # ImagingGrid instance
        function_types: Dict[int, Any],
        ms_functions: List[int],
    ):
        """Initialize Waters metadata extractor.

        Args:
            ml: MassLynxLib instance for native library access.
            handle: Opaque file handle from ml.open_file().
            data_path: Path to the Waters .raw directory.
            imaging_grid: Pre-built ImagingGrid with spatial metadata.
            function_types: Map of function index to FunctionType.
            ms_functions: List of MS function indices.
        """
        super().__init__(ml)
        self._ml = ml
        self._handle = handle
        self._data_path = data_path
        self._imaging_grid = imaging_grid
        self._function_types = function_types
        self._ms_functions = ms_functions

    def _extract_essential_impl(self) -> EssentialMetadata:
        """Extract essential metadata.

        Uses imaging grid for dimensions/pixel_size,
        scans all MS function spectra for mass range and peak counts.
        """
        grid = self._imaging_grid
        dimensions = grid.dimensions  # (n_x, n_y, 1)

        # Coordinate bounds (0-based pixel indices)
        coordinate_bounds = (
            0.0,
            float(grid.pixel_count_x - 1),
            0.0,
            float(grid.pixel_count_y - 1),
        )

        # Pixel size in micrometers
        pixel_size: Optional[Tuple[float, float]] = None
        if grid.pixel_size_x > 0 and grid.pixel_size_y > 0:
            pixel_size = (grid.pixel_size_x, grid.pixel_size_y)

        # Scan all MS spectra for mass range, spectrum count, peak counts
        mass_range, n_spectra, total_peaks, peak_counts = self._scan_all_ms_spectra(
            dimensions
        )

        # Memory estimate: total_peaks * 2 values (mz + intensity) * 8 bytes
        estimated_memory_gb = (total_peaks * 2 * 8) / (1024**3)

        # Spectrum type -- determine from first MS function
        spectrum_type = self._detect_spectrum_type()

        return EssentialMetadata(
            dimensions=dimensions,
            coordinate_bounds=coordinate_bounds,
            mass_range=mass_range,
            pixel_size=pixel_size,
            n_spectra=n_spectra,
            total_peaks=total_peaks,
            estimated_memory_gb=estimated_memory_gb,
            source_path=str(self._data_path),
            spectrum_type=spectrum_type,
            peak_counts_per_pixel=peak_counts,
        )

    def _detect_spectrum_type(self) -> Optional[str]:
        """Detect whether data is centroid or profile."""
        if self._ms_functions:
            is_profile = self._ml.is_raw_spectrum_profile(
                self._handle, self._ms_functions[0]
            )
            if is_profile:
                return "profile spectrum"
            return "centroid spectrum"
        return None

    def _read_scan_mzs(self, func: int, scan: int) -> Optional[NDArray[np.floating]]:
        """Read m/z array for a single scan, returning None on failure."""
        scan_info = self._imaging_grid.scan_map.get((func, scan))
        if scan_info is None or not scan_info.has_position:
            return None

        coords = self._imaging_grid.get_coordinates(scan_info)
        if coords is None:
            return None

        try:
            mzs, _ = self._ml.read_spectrum(self._handle, func, scan)
        except Exception as e:
            logger.debug(f"Failed to read spectrum func={func} scan={scan}: {e}")
            return None

        return mzs if len(mzs) > 0 else None

    def _update_scan_stats(
        self,
        func: int,
        scan: int,
        mzs: NDArray[np.floating],
        stats: Dict[str, Any],
    ) -> None:
        """Update running statistics with data from one scan."""
        n_peaks = len(mzs)
        stats["min_mass"] = min(stats["min_mass"], float(mzs[0]))
        stats["max_mass"] = max(stats["max_mass"], float(mzs[-1]))
        stats["total_peaks"] += n_peaks
        stats["n_spectra"] += 1

        coords = self._imaging_grid.get_coordinates(
            self._imaging_grid.scan_map[(func, scan)]
        )
        x, y, z = coords
        n_x, n_y = stats["n_x"], stats["n_y"]
        pixel_idx = z * (n_x * n_y) + y * n_x + x
        if 0 <= pixel_idx < stats["n_pixels"]:
            stats["peak_counts"][pixel_idx] = n_peaks

    def _scan_all_ms_spectra(
        self,
        dimensions: Tuple[int, int, int],
    ) -> Tuple[Tuple[float, float], int, int, Optional[NDArray[np.int32]]]:
        """Single pass over all MS scans for mass range, counts, per-pixel peak counts.

        Returns:
            (mass_range, n_spectra, total_peaks, peak_counts_per_pixel)
        """
        n_x, n_y, n_z = dimensions
        n_pixels = n_x * n_y * n_z

        stats: Dict[str, Any] = {
            "min_mass": float("inf"),
            "max_mass": float("-inf"),
            "total_peaks": 0,
            "n_spectra": 0,
            "peak_counts": np.zeros(n_pixels, dtype=np.int32),
            "n_x": n_x,
            "n_y": n_y,
            "n_pixels": n_pixels,
        }

        total_scans = sum(
            self._ml.get_number_of_scans_in_function(self._handle, f)
            for f in self._ms_functions
        )

        with tqdm(
            total=total_scans,
            desc="Scanning Waters metadata",
            unit="scan",
        ) as pbar:
            for func in self._ms_functions:
                n_scans = self._ml.get_number_of_scans_in_function(self._handle, func)
                for scan in range(n_scans):
                    pbar.update(1)
                    mzs = self._read_scan_mzs(func, scan)
                    if mzs is not None:
                        self._update_scan_stats(func, scan, mzs, stats)

        if stats["min_mass"] == float("inf"):
            logger.warning("No valid spectra found in Waters data")
            stats["min_mass"], stats["max_mass"] = 0.0, 1000.0

        logger.info(
            f"Waters metadata scan complete: {stats['n_spectra']} spectra, "
            f"mass range {stats['min_mass']:.2f}-{stats['max_mass']:.2f}, "
            f"total peaks {stats['total_peaks']:,}"
        )

        return (
            (stats["min_mass"], stats["max_mass"]),
            stats["n_spectra"],
            stats["total_peaks"],
            stats["peak_counts"],
        )

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        """Extract comprehensive metadata including Waters-specific details."""
        essential = self.get_essential()

        return ComprehensiveMetadata(
            essential=essential,
            format_specific=self._extract_waters_specific(),
            acquisition_params=self._extract_acquisition_params(),
            instrument_info=self._extract_instrument_info(),
            raw_metadata=self._extract_raw_metadata(),
        )

    def _extract_waters_specific(self) -> Dict[str, Any]:
        """Extract Waters format-specific metadata."""
        return {
            "data_format": "waters_raw",
            "data_path": str(self._data_path),
            "is_imaging": True,
            "n_functions": self._ml.get_number_of_functions(self._handle),
            "ms_functions": self._ms_functions,
            "function_types": {
                str(f): ft.name if hasattr(ft, "name") else str(ft)
                for f, ft in self._function_types.items()
            },
            "pixel_count_x": self._imaging_grid.pixel_count_x,
            "pixel_count_y": self._imaging_grid.pixel_count_y,
            "lateral_width_um": self._imaging_grid.lateral_width,
            "lateral_height_um": self._imaging_grid.lateral_height,
        }

    def _extract_acquisition_params(self) -> Dict[str, Any]:
        """Extract acquisition parameters."""
        params: Dict[str, Any] = {}

        acq_date = self._ml.get_acquisition_date(self._handle)
        if acq_date:
            params["acquisition_date"] = acq_date

        # Profile/centroid info per MS function
        for func in self._ms_functions:
            is_profile = self._ml.is_raw_spectrum_profile(self._handle, func)
            params[f"function_{func}_is_profile"] = is_profile

        # Acquisition mass range per function
        for func in self._ms_functions:
            acq_range = self._ml.get_acquisition_range(self._handle, func)
            if acq_range:
                params[f"function_{func}_acq_range_start"] = acq_range[0]
                params[f"function_{func}_acq_range_end"] = acq_range[1]

        # Lock mass info
        params["is_lockmass_corrected"] = self._ml.is_lockmass_corrected(self._handle)
        lm_func = self._ml.get_lockmass_function(self._handle)
        params["lockmass_function"] = lm_func if lm_func >= 0 else None

        return params

    def _extract_instrument_info(self) -> Dict[str, Any]:
        """Extract instrument information."""
        return {
            "vendor": "Waters",
            "format": "MassLynx .raw",
        }

    def _extract_raw_metadata(self) -> Dict[str, Any]:
        """Extract raw metadata dictionary."""
        return {
            "data_path": str(self._data_path),
            "imaging_grid": {
                "pixel_count_x": self._imaging_grid.pixel_count_x,
                "pixel_count_y": self._imaging_grid.pixel_count_y,
                "pixel_size_x_um": self._imaging_grid.pixel_size_x,
                "pixel_size_y_um": self._imaging_grid.pixel_size_y,
                "lateral_width_um": self._imaging_grid.lateral_width,
                "lateral_height_um": self._imaging_grid.lateral_height,
            },
        }
