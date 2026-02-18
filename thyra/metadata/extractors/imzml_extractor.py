# thyra/metadata/extractors/imzml_extractor.py
import gc
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from pyimzml.ImzMLParser import ImzMLParser

from ...core.base_extractor import MetadataExtractor
from ...resampling.constants import ImzMLAccessions, SpectrumType
from ..types import ComprehensiveMetadata, EssentialMetadata

logger = logging.getLogger(__name__)


class ImzMLMetadataExtractor(MetadataExtractor):
    """ImzML-specific metadata extractor with optimized two-phase extraction."""

    def __init__(self, parser: ImzMLParser, imzml_path: Path):
        """Initialize ImzML metadata extractor.

        Args:
            parser: Initialized ImzML parser
            imzml_path: Path to the ImzML file
        """
        super().__init__(parser)
        self.parser = parser
        self.imzml_path = imzml_path

    def _extract_essential_impl(self) -> EssentialMetadata:
        """Extract essential metadata optimized for speed."""
        # Single coordinate scan for efficiency
        coords = np.array(self.parser.coordinates)

        if len(coords) == 0:
            raise ValueError("No coordinates found in ImzML file")

        dimensions = self._calculate_dimensions(coords)
        coordinate_bounds = self._calculate_bounds(coords)
        # Pass dimensions to collect per-pixel peak counts during scan
        mass_range, total_peaks, peak_counts = self._get_mass_range_complete(
            dimensions=dimensions
        )
        pixel_size = self._extract_pixel_size_fast()
        n_spectra = len(coords)
        estimated_memory = self._estimate_memory(n_spectra)

        # Check for centroid spectrum
        spectrum_type = self._detect_centroid_spectrum()

        return EssentialMetadata(
            dimensions=dimensions,
            coordinate_bounds=coordinate_bounds,
            mass_range=mass_range,
            pixel_size=pixel_size,
            n_spectra=n_spectra,
            total_peaks=total_peaks,
            estimated_memory_gb=estimated_memory,
            source_path=str(self.imzml_path),
            spectrum_type=spectrum_type,
            peak_counts_per_pixel=peak_counts,
        )

    def _extract_comprehensive_impl(self) -> ComprehensiveMetadata:
        """Extract comprehensive metadata with full XML parsing."""
        essential = self.get_essential()

        return ComprehensiveMetadata(
            essential=essential,
            format_specific=self._extract_imzml_specific(),
            acquisition_params=self._extract_acquisition_params(),
            instrument_info=self._extract_instrument_info(),
            raw_metadata=self._extract_raw_metadata(),
        )

    def _calculate_dimensions(self, coords: NDArray[np.int_]) -> Tuple[int, int, int]:
        """Calculate dataset dimensions from coordinates."""
        if len(coords) == 0:
            return (0, 0, 0)

        # Coordinates are 1-based in ImzML, convert to 0-based for calculation
        coords_0based = coords - 1

        max_coords = np.max(coords_0based, axis=0)
        return (
            int(max_coords[0]) + 1,
            int(max_coords[1]) + 1,
            int(max_coords[2]) + 1,
        )

    def _calculate_bounds(
        self, coords: NDArray[np.int_]
    ) -> Tuple[float, float, float, float]:
        """Calculate coordinate bounds (min_x, max_x, min_y, max_y)."""
        if len(coords) == 0:
            return (0.0, 0.0, 0.0, 0.0)

        # Convert to spatial coordinates (assuming 1-based indexing)
        x_coords = coords[:, 0].astype(float)
        y_coords = coords[:, 1].astype(float)

        return (
            float(np.min(x_coords)),
            float(np.max(x_coords)),
            float(np.min(y_coords)),
            float(np.max(y_coords)),
        )

    def _is_continuous_mode(self) -> bool:
        """Check if the ImzML file is in continuous mode.

        Continuous mode means all spectra share the same m/z axis.
        Detection checks file_description.param_by_name for "continuous" key.

        Returns:
            True if continuous mode, False otherwise (processed mode).
        """
        try:
            if not hasattr(self.parser, "metadata"):
                return False
            if self.parser.metadata is None:
                return False

            file_desc = getattr(self.parser.metadata, "file_description", None)
            if file_desc is None:
                return False

            param_by_name = getattr(file_desc, "param_by_name", None)
            if param_by_name is None or not isinstance(param_by_name, dict):
                return False

            return "continuous" in param_by_name
        except Exception:
            return False

    def _get_mass_range_complete(
        self,
        dimensions: Optional[Tuple[int, int, int]] = None,
    ) -> Tuple[Tuple[float, float], int, Optional[NDArray[np.int32]]]:
        """Complete mass range extraction.

        For continuous mode: reads only the first spectrum (all share same m/z axis).
        For processed mode: scans ALL spectra to find complete mass range.

        Also counts total peaks for COO matrix pre-allocation and
        optionally collects per-pixel peak counts for streaming conversion.

        Args:
            dimensions: Optional (n_x, n_y, n_z) grid dimensions.
                If provided, per-pixel peak counts will be collected.

        Returns:
            Tuple of ((min_mass, max_mass), total_peaks, peak_counts_per_pixel)
            peak_counts_per_pixel is None if dimensions not provided.
        """
        try:
            coords = self.parser.coordinates
            n_spectra = len(coords)

            # Check if continuous mode - all spectra share the same m/z axis
            # Detection: check file_description.param_by_name for "continuous" key
            is_continuous = self._is_continuous_mode()

            if is_continuous:
                return self._get_mass_range_continuous(n_spectra, dimensions)
            else:
                return self._get_mass_range_processed(coords, n_spectra, dimensions)

        except Exception as e:
            logger.error(f"Mass range extraction failed: {e}")
            return ((0.0, 1000.0), 0, None)

    def _get_mass_range_continuous(
        self,
        n_spectra: int,
        dimensions: Optional[Tuple[int, int, int]] = None,
    ) -> Tuple[Tuple[float, float], int, Optional[NDArray[np.int32]]]:
        """Get mass range for continuous mode - read only first spectrum.

        In continuous mode, all spectra share the same m/z axis, so we only
        need to read one spectrum to get the mass range and peak count.
        """
        logger.info(
            "Continuous mode detected - reading m/z axis from first spectrum only"
        )

        # Read first spectrum to get shared m/z axis
        mzs, _ = self.parser.getspectrum(0)
        n_peaks_per_spectrum = len(mzs)

        if n_peaks_per_spectrum == 0:
            logger.warning("First spectrum has no peaks")
            return ((0.0, 1000.0), 0, None)

        min_mass = float(np.min(mzs))
        max_mass = float(np.max(mzs))
        total_peaks = n_peaks_per_spectrum * n_spectra

        # For continuous mode, all pixels have the same peak count
        peak_counts = None
        if dimensions is not None:
            n_x, n_y, n_z = dimensions
            n_pixels = n_x * n_y * n_z
            peak_counts = np.full(n_pixels, n_peaks_per_spectrum, dtype=np.int32)
            logger.info(
                f"All {n_pixels:,} pixels have {n_peaks_per_spectrum:,} peaks (continuous mode)"
            )

        logger.info(f"Mass range: {min_mass:.2f} - {max_mass:.2f} m/z")
        logger.info(
            f"Total peaks: {total_peaks:,} ({n_peaks_per_spectrum:,} per spectrum)"
        )
        return ((min_mass, max_mass), total_peaks, peak_counts)

    def _get_mass_range_processed(
        self,
        coords: List,
        n_spectra: int,
        dimensions: Optional[Tuple[int, int, int]] = None,
    ) -> Tuple[Tuple[float, float], int, Optional[NDArray[np.int32]]]:
        """Get mass range for processed mode - scan all spectra.

        In processed mode, each spectrum can have different m/z values,
        so we must scan all spectra to find the complete mass range.
        """
        logger.info("Processed mode - scanning ALL spectra for complete mass range...")

        peak_counts = self._init_peak_counts_array(dimensions)

        min_mass, max_mass, total_peaks = self._scan_all_spectra(
            coords, n_spectra, dimensions, peak_counts
        )

        if min_mass == float("inf"):
            logger.warning("No valid spectra found")
            return ((0.0, 1000.0), 0, None)

        logger.info(f"Complete mass range: {min_mass:.2f} - {max_mass:.2f} m/z")
        logger.info(f"Total peaks: {total_peaks:,}")
        return ((min_mass, max_mass), total_peaks, peak_counts)

    def _init_peak_counts_array(
        self, dimensions: Optional[Tuple[int, int, int]]
    ) -> Optional[NDArray[np.int32]]:
        """Initialize per-pixel peak counts array if dimensions provided.

        Args:
            dimensions: Optional (n_x, n_y, n_z) grid dimensions.

        Returns:
            Array of zeros or None if dimensions not provided.
        """
        if dimensions is None:
            return None
        n_x, n_y, n_z = dimensions
        n_pixels = n_x * n_y * n_z
        logger.info(f"Collecting per-pixel peak counts ({n_pixels:,} pixels)")
        return np.zeros(n_pixels, dtype=np.int32)

    def _scan_all_spectra(
        self,
        coords: List,
        n_spectra: int,
        dimensions: Optional[Tuple[int, int, int]],
        peak_counts: Optional[NDArray[np.int32]],
    ) -> Tuple[float, float, int]:
        """Scan all spectra to find mass range and count peaks.

        Args:
            coords: List of spectrum coordinates.
            n_spectra: Total number of spectra.
            dimensions: Optional grid dimensions for pixel indexing.
            peak_counts: Optional array to store per-pixel peak counts.

        Returns:
            Tuple of (min_mass, max_mass, total_peaks).
        """
        from tqdm import tqdm

        min_mass = float("inf")
        max_mass = float("-inf")
        total_peaks = 0

        with tqdm(
            total=n_spectra,
            desc="Scanning mass range and counting peaks",
            unit="spectrum",
        ) as pbar:
            for idx in range(n_spectra):
                result = self._process_spectrum_for_range(
                    idx, coords, dimensions, peak_counts
                )
                if result is not None:
                    spec_min, spec_max, n_peaks = result
                    min_mass = min(min_mass, spec_min)
                    max_mass = max(max_mass, spec_max)
                    total_peaks += n_peaks

                if idx % 50000 == 0 and idx > 0:
                    gc.collect()

                pbar.update(1)

        return min_mass, max_mass, total_peaks

    def _process_spectrum_for_range(
        self,
        idx: int,
        coords: List,
        dimensions: Optional[Tuple[int, int, int]],
        peak_counts: Optional[NDArray[np.int32]],
    ) -> Optional[Tuple[float, float, int]]:
        """Process a single spectrum for mass range and peak count.

        Args:
            idx: Spectrum index.
            coords: List of spectrum coordinates.
            dimensions: Optional grid dimensions.
            peak_counts: Optional array for per-pixel counts.

        Returns:
            Tuple of (min_mz, max_mz, n_peaks) or None if failed.
        """
        try:
            mzs, intensities = self.parser.getspectrum(idx)
            n_peaks = len(mzs)

            # Store per-pixel count if tracking
            if peak_counts is not None and dimensions is not None:
                self._store_pixel_peak_count(
                    idx, coords, dimensions, peak_counts, n_peaks
                )

            # Release references
            del intensities

            if n_peaks > 0:
                result = (float(np.min(mzs)), float(np.max(mzs)), n_peaks)
                del mzs
                return result

            del mzs
            return None

        except Exception as e:
            logger.debug(f"Failed to read spectrum {idx}: {e}")
            return None

    def _store_pixel_peak_count(
        self,
        idx: int,
        coords: List,
        dimensions: Tuple[int, int, int],
        peak_counts: NDArray[np.int32],
        n_peaks: int,
    ) -> None:
        """Store peak count for a pixel.

        Args:
            idx: Spectrum index.
            coords: List of spectrum coordinates.
            dimensions: Grid dimensions (n_x, n_y, n_z).
            peak_counts: Array to store counts.
            n_peaks: Number of peaks in this spectrum.
        """
        # ImzML coordinates are 1-based
        x, y, z = coords[idx]
        x, y, z = x - 1, y - 1, max(z - 1, 0)
        n_x, n_y, n_z = dimensions
        pixel_idx = z * (n_x * n_y) + y * n_x + x
        if 0 <= pixel_idx < len(peak_counts):
            peak_counts[pixel_idx] = n_peaks

    def get_mass_range_for_resampling(self) -> Tuple[float, float]:
        """Get accurate mass range required for resampling.

        This performs a complete scan of all spectra to ensure no m/z
        values are missed when building the resampled axis.
        """
        mass_range, _, _ = self._get_mass_range_complete()
        return mass_range

    def _extract_pixel_size_fast(self) -> Optional[Tuple[float, float]]:
        """Fast pixel size extraction from imzmldict first."""
        if hasattr(self.parser, "imzmldict") and self.parser.imzmldict:
            # Check for pixel size parameters in the parsed dictionary
            x_size = self.parser.imzmldict.get("pixel size x")
            y_size = self.parser.imzmldict.get("pixel size y")

            if x_size is not None and y_size is not None:
                try:
                    return (float(x_size), float(y_size))
                except (ValueError, TypeError):
                    pass

        return None  # Defer to comprehensive extraction

    def _estimate_memory(self, n_spectra: int) -> float:
        """Estimate memory usage in GB."""
        # Rough estimate: assume average 1000 peaks per spectrum,
        # 8 bytes per float
        avg_peaks_per_spectrum = 1000
        bytes_per_value = 8  # float64
        estimated_bytes = (
            n_spectra * avg_peaks_per_spectrum * 2 * bytes_per_value
        )  # mz + intensity
        return estimated_bytes / (1024**3)  # Convert to GB

    def _extract_imzml_specific(self) -> Dict[str, Any]:
        """Extract ImzML format-specific metadata."""
        format_specific: Dict[str, Any] = {
            "imzml_version": "1.1.0",  # Default version
            "file_mode": (
                "continuous"
                if getattr(self.parser, "continuous", False)
                else "processed"
            ),
            "ibd_file": str(self.imzml_path.with_suffix(".ibd")),
            "uuid": None,
            "spectrum_count": len(self.parser.coordinates),
            "scan_settings": {},
        }

        # Extract UUID if available
        try:
            if hasattr(self.parser, "metadata") and hasattr(
                self.parser.metadata, "file_description"
            ):
                cv_params = getattr(
                    self.parser.metadata.file_description, "cv_params", []
                )
                if cv_params and len(cv_params) > 0:
                    format_specific["uuid"] = cv_params[0][2]
        except Exception as e:
            logger.debug(f"Could not extract UUID: {e}")

        return format_specific

    def _extract_acquisition_params(self) -> Dict[str, Any]:
        """Extract acquisition parameters from XML metadata."""
        params = {}

        # Extract pixel size with full XML parsing if not found in fast
        # extraction
        if not self.get_essential().has_pixel_size:
            pixel_size = self._extract_pixel_size_from_xml()
            if pixel_size:
                params["pixel_size_x_um"] = pixel_size[0]
                params["pixel_size_y_um"] = pixel_size[1]

        # Add other acquisition parameters from imzmldict
        if hasattr(self.parser, "imzmldict") and self.parser.imzmldict:
            acquisition_keys = [
                "scan direction",
                "scan pattern",
                "scan type",
                "laser power",
                "laser frequency",
                "laser spot size",
            ]
            for key in acquisition_keys:
                if key in self.parser.imzmldict:
                    params[key.replace(" ", "_")] = self.parser.imzmldict[key]

        return params

    def _extract_instrument_info(self) -> Dict[str, Any]:
        """Extract instrument information."""
        instrument = {}

        if hasattr(self.parser, "imzmldict") and self.parser.imzmldict:
            instrument_keys = [
                "instrument model",
                "instrument serial number",
                "software",
                "software version",
            ]
            for key in instrument_keys:
                if key in self.parser.imzmldict:
                    instrument[key.replace(" ", "_")] = self.parser.imzmldict[key]

        return instrument

    def _extract_raw_metadata(self) -> Dict[str, Any]:
        """Extract raw metadata from imzmldict and spectrum cvParams."""
        raw_metadata = {}

        if hasattr(self.parser, "imzmldict") and self.parser.imzmldict:
            raw_metadata = dict(self.parser.imzmldict)

        # Extract spectrum-level cvParams for centroid detection
        cv_params = self._extract_spectrum_cvparams()
        if cv_params:
            raw_metadata["cvParams"] = cv_params

        return raw_metadata

    def _extract_spectrum_cvparams(self) -> Optional[List[Dict[str, Any]]]:
        """Extract cvParams from first spectrum for centroid detection."""
        try:
            if not hasattr(self.parser, "metadata") or not self.parser.metadata:
                return None

            # Look for spectrum-level cvParams in the metadata
            if hasattr(self.parser.metadata, "file_description"):
                file_desc = self.parser.metadata.file_description
                if hasattr(file_desc, "param_by_name"):
                    params = file_desc.param_by_name
                    cv_params = []

                    # Check for centroid spectrum in file description
                    for name, value in params.items():
                        cv_params.append({"name": name, "value": value})

                    return cv_params

            return None
        except Exception as e:
            logger.debug(f"Could not extract spectrum cvParams: {e}")
            return None

    def _detect_centroid_spectrum(self) -> Optional[str]:
        """Detect spectrum type by looking for MS:1000127 (centroid) or MS:1000128 (profile)."""
        try:
            # Method 1: Check parser metadata first (no XML parsing needed)
            result = self._check_parser_metadata_for_centroid()
            if result:
                return result

            # Method 2: Stream-parse XML for spectrum type markers (memory efficient)
            result = self._check_xml_for_spectrum_type()
            if result:
                return result

            return None
        except Exception as e:
            logger.debug(f"Could not detect spectrum type: {e}")
            return None

    def _check_xml_for_spectrum_type(self) -> Optional[str]:
        """Check XML for spectrum type markers using streaming parser.

        Looks for PSI-MS controlled vocabulary accession codes:
        - MS:1000127 = centroid spectrum
        - MS:1000128 = profile spectrum

        Uses iterparse for memory-efficient streaming - stops as soon as a
        spectrum type marker is found, avoiding full file load for large datasets.
        """
        try:
            import xml.etree.ElementTree as ET  # nosec B405

            # Use iterparse for streaming - only parses until we find what we need
            # The spectrum type cvParam is typically in the fileDescription
            # section near the beginning of the file
            context = ET.iterparse(str(self.imzml_path), events=("end",))  # nosec B314

            elements_checked = 0
            max_elements = 10000  # Limit search to first 10k elements

            for event, elem in context:
                elements_checked += 1

                if elem.tag.endswith("cvParam"):
                    accession = elem.get("accession", "")
                    if accession == ImzMLAccessions.CENTROID_SPECTRUM:
                        logger.info(
                            f"Detected centroid spectrum from {ImzMLAccessions.CENTROID_SPECTRUM}"
                        )
                        del context
                        return SpectrumType.CENTROID
                    if accession == ImzMLAccessions.PROFILE_SPECTRUM:
                        logger.info(
                            f"Detected profile spectrum from {ImzMLAccessions.PROFILE_SPECTRUM}"
                        )
                        del context
                        return SpectrumType.PROFILE

                # Clear processed elements to save memory
                elem.clear()

                # Stop after checking enough elements - spectrum type info is at the start
                if elements_checked >= max_elements:
                    logger.debug(
                        f"Spectrum type marker not found in first {max_elements} elements"
                    )
                    break

            del context
        except Exception as e:
            logger.debug(f"XML streaming parse failed: {e}")
        return None

    def _get_xml_parser(self):
        """Get XML parser, preferring defusedxml for security."""
        try:
            # Use defusedxml for secure parsing
            import defusedxml.ElementTree as ET

            return ET
        except ImportError:
            # Fallback to standard library with warning
            import xml.etree.ElementTree as ET  # nosec B405

            logger.warning("defusedxml not available, using xml.etree.ElementTree")
            return ET

    def _check_parser_metadata_for_centroid(self) -> Optional[str]:
        """Check parser metadata for processed flag indicating centroid data."""
        if not (hasattr(self.parser, "metadata") and self.parser.metadata):
            return None

        if not hasattr(self.parser.metadata, "file_description"):
            return None

        file_desc = self.parser.metadata.file_description
        if not hasattr(file_desc, "param_by_name"):
            return None

        params = file_desc.param_by_name
        # If it's processed data, it's likely centroided
        if params.get("processed", False):
            logger.info("Assuming centroid spectrum for processed ImzML data")
            return SpectrumType.CENTROID

        return None

    def _extract_pixel_size_from_xml(self) -> Optional[Tuple[float, float]]:
        """Extract pixel size using full XML parsing as fallback."""
        try:
            if not hasattr(self.parser, "metadata") or not hasattr(
                self.parser.metadata, "root"
            ):
                return None

            root = self.parser.metadata.root

            # Define namespaces for XML parsing
            namespaces = {
                "mzml": "http://psi.hupo.org/ms/mzml",
                "ims": "http://www.maldi-msi.org/download/imzml/imagingMS.obo",
            }

            x_size = None
            y_size = None

            # Search for cvParam elements with the pixel size accessions
            for cvparam in root.findall(".//mzml:cvParam", namespaces):
                accession = cvparam.get("accession")
                if accession == "IMS:1000046":  # pixel size x
                    x_size = float(cvparam.get("value", 0))
                elif accession == "IMS:1000047":  # pixel size y
                    y_size = float(cvparam.get("value", 0))

            if x_size is not None and y_size is not None:
                logger.info(f"Detected pixel size from XML: x={x_size}μm, y={y_size}μm")
                return (x_size, y_size)

        except Exception as e:
            logger.warning(f"Failed to parse XML metadata for pixel size: {e}")

        return None
