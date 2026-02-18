"""Mass axis builder for efficient common mass axis construction.

This module provides optimized methods for building common mass axes
from large datasets with memory efficiency and progress tracking.
"""

import logging
from typing import Any, Callable, Dict, Iterator, List, Optional, Set

import numpy as np
from tqdm import tqdm

logger = logging.getLogger(__name__)


class MassAxisBuilder:
    """Efficient builder for common mass axes across multiple spectra.

    This class provides multiple strategies for building mass axes
    depending on dataset size and memory constraints.
    """

    def __init__(
        self,
        strategy: str = "auto",
        memory_limit_mb: float = 1024,
        progress_callback: Optional[Callable] = None,
    ):
        """Initialize the mass axis builder.

        Args:
            strategy: Building strategy ('auto', 'memory_efficient',
                'fast', 'streaming')
            memory_limit_mb: Memory limit for mass axis construction
            progress_callback: Optional callback for progress updates
        """
        self.strategy = strategy
        self.memory_limit_mb = memory_limit_mb
        self.progress_callback = progress_callback

        # Statistics
        self._stats = {
            "total_spectra_processed": 0,
            "total_peaks_processed": 0,
            "unique_masses_found": 0,
            "memory_usage_mb": 0.0,
        }

        logger.info(f"Initialized MassAxisBuilder with strategy: {strategy}")

    def build_from_spectra_iterator(
        self, spectra_iterator: Iterator, total_spectra: Optional[int] = None
    ) -> np.ndarray:
        """Build common mass axis from a spectra iterator.

        Args:
            spectra_iterator: Iterator yielding (coords, mzs, intensities)
                tuples
            total_spectra: Optional total number of spectra for progress
                tracking

        Returns:
            numpy array of unique m/z values in ascending order
        """
        if self.strategy == "auto":
            strategy = self._determine_optimal_strategy(total_spectra)
        else:
            strategy = self.strategy

        logger.info(f"Building mass axis using strategy: {strategy}")

        if strategy == "memory_efficient":
            return self._build_memory_efficient(spectra_iterator, total_spectra)
        elif strategy == "fast":
            return self._build_fast(spectra_iterator, total_spectra)
        elif strategy == "streaming":
            return self._build_streaming(spectra_iterator, total_spectra)
        else:
            raise ValueError(f"Unknown strategy: {strategy}")

    def _determine_optimal_strategy(self, total_spectra: Optional[int]) -> str:
        """Automatically determine the optimal strategy based on dataset size.

        Args:
            total_spectra: Total number of spectra

        Returns:
            Optimal strategy name
        """
        if total_spectra is None:
            return "streaming"  # Safe default for unknown size

        # Estimate memory requirements
        estimated_peaks_per_spectrum = 1500  # Conservative estimate
        estimated_total_peaks = total_spectra * estimated_peaks_per_spectrum
        estimated_memory_mb = (
            estimated_total_peaks * 8 / 1024 / 1024
        )  # 8 bytes per m/z value

        if estimated_memory_mb < self.memory_limit_mb / 2:
            return "fast"
        elif estimated_memory_mb < self.memory_limit_mb:
            return "memory_efficient"
        else:
            return "streaming"

    def _build_fast(
        self, spectra_iterator: Iterator, total_spectra: Optional[int]
    ) -> np.ndarray:
        """Fast mass axis building by collecting all m/z values in memory.

        This is the fastest method but uses the most memory.
        """
        logger.info("Building mass axis using fast strategy")

        all_mzs = []

        # Setup progress bar
        pbar = tqdm(
            total=total_spectra,
            desc="Building mass axis (fast)",
            unit="spectrum",
            disable=not total_spectra,
        )

        try:
            for coords, mzs, intensities in spectra_iterator:
                if mzs.size > 0:
                    all_mzs.append(mzs)
                    self._stats["total_peaks_processed"] += len(mzs)

                self._stats["total_spectra_processed"] += 1
                pbar.update(1)

                if self.progress_callback:
                    self.progress_callback(
                        self._stats["total_spectra_processed"], total_spectra
                    )

        finally:
            pbar.close()

        if not all_mzs:
            logger.warning("No m/z values found")
            return np.array([])

        # Combine and find unique values
        logger.info("Combining and sorting m/z values")
        combined_mzs = np.concatenate(all_mzs)
        unique_mzs = np.unique(combined_mzs)

        self._stats["unique_masses_found"] = len(unique_mzs)
        logger.info(
            f"Found {len(unique_mzs)} unique m/z values from " f"{len(all_mzs)} spectra"
        )

        return unique_mzs

    def _build_memory_efficient(
        self, spectra_iterator: Iterator, total_spectra: Optional[int]
    ) -> np.ndarray:
        """Memory-efficient mass axis building using chunked processing.

        This method processes spectra in chunks to limit memory usage.
        """
        logger.info("Building mass axis using memory-efficient strategy")

        chunk_size = 100  # Process 100 spectra at a time
        chunk_mzs = []
        all_unique_chunks = []

        pbar = tqdm(
            total=total_spectra,
            desc="Building mass axis (memory efficient)",
            unit="spectrum",
            disable=not total_spectra,
        )

        try:
            for coords, mzs, intensities in spectra_iterator:
                if mzs.size > 0:
                    chunk_mzs.append(mzs)
                    self._stats["total_peaks_processed"] += len(mzs)

                self._stats["total_spectra_processed"] += 1
                pbar.update(1)

                # Process chunk when it reaches target size
                if len(chunk_mzs) >= chunk_size:
                    chunk_unique = self._process_mz_chunk(chunk_mzs)
                    if len(chunk_unique) > 0:
                        all_unique_chunks.append(chunk_unique)
                    chunk_mzs = []

                if self.progress_callback:
                    self.progress_callback(
                        self._stats["total_spectra_processed"], total_spectra
                    )

            # Process final chunk
            if chunk_mzs:
                chunk_unique = self._process_mz_chunk(chunk_mzs)
                if len(chunk_unique) > 0:
                    all_unique_chunks.append(chunk_unique)

        finally:
            pbar.close()

        if not all_unique_chunks:
            logger.warning("No m/z values found")
            return np.array([])

        # Combine all chunk results
        logger.info(f"Combining {len(all_unique_chunks)} chunks")
        combined_unique = np.concatenate(all_unique_chunks)
        final_unique = np.unique(combined_unique)

        self._stats["unique_masses_found"] = len(final_unique)
        logger.info(f"Found {len(final_unique)} unique m/z values")

        return final_unique

    def _build_streaming(
        self, spectra_iterator: Iterator, total_spectra: Optional[int]
    ) -> np.ndarray:
        """Streaming mass axis building using set-based accumulation.

        This method uses a set to accumulate unique m/z values without
        storing all raw data in memory.
        """
        logger.info("Building mass axis using streaming strategy")

        unique_mzs_set: Set[float] = set()
        batch_size = 50
        batch_count = 0

        pbar = tqdm(
            total=total_spectra,
            desc="Building mass axis (streaming)",
            unit="spectrum",
            disable=not total_spectra,
        )

        try:
            for coords, mzs, intensities in spectra_iterator:
                if mzs.size > 0:
                    # Add m/z values to set (automatically handles uniqueness)
                    unique_mzs_set.update(mzs)
                    self._stats["total_peaks_processed"] += len(mzs)

                self._stats["total_spectra_processed"] += 1
                batch_count += 1
                pbar.update(1)

                # Periodically log progress for very large datasets
                if batch_count % batch_size == 0:
                    logger.debug(
                        f"Processed {self._stats['total_spectra_processed']} "
                        f"spectra, found {len(unique_mzs_set)} unique m/z "
                        f"values"
                    )

                if self.progress_callback:
                    self.progress_callback(
                        self._stats["total_spectra_processed"], total_spectra
                    )

        finally:
            pbar.close()

        if not unique_mzs_set:
            logger.warning("No m/z values found")
            return np.array([])

        # Convert set to sorted array
        logger.info("Converting set to sorted array")
        unique_mzs = np.array(sorted(unique_mzs_set))

        self._stats["unique_masses_found"] = len(unique_mzs)
        logger.info(f"Found {len(unique_mzs)} unique m/z values")

        return unique_mzs

    def _process_mz_chunk(self, mz_list: List[np.ndarray]) -> np.ndarray:
        """Process a chunk of m/z arrays to find unique values.

        Args:
            mz_list: List of m/z arrays

        Returns:
            Array of unique m/z values from the chunk
        """
        if not mz_list:
            return np.array([])

        try:
            combined = np.concatenate(mz_list)
            return np.unique(combined)
        except Exception as e:
            logger.warning(f"Error processing m/z chunk: {e}")
            return np.array([])

    def build_with_tolerance(
        self,
        spectra_iterator: Iterator,
        tolerance_ppm: float = 5.0,
        total_spectra: Optional[int] = None,
    ) -> np.ndarray:
        """Build mass axis with m/z tolerance for peak consolidation.

        Args:
            spectra_iterator: Iterator yielding (coords, mzs, intensities)
                tuples
            tolerance_ppm: Mass tolerance in ppm for peak consolidation
            total_spectra: Optional total number of spectra

        Returns:
            numpy array of consolidated m/z values
        """
        logger.info(f"Building mass axis with {tolerance_ppm} ppm tolerance")

        # First build regular mass axis
        unique_mzs = self.build_from_spectra_iterator(spectra_iterator, total_spectra)

        if len(unique_mzs) == 0:
            return unique_mzs

        # Apply tolerance-based consolidation
        consolidated_mzs = self._consolidate_with_tolerance(unique_mzs, tolerance_ppm)

        logger.info(
            f"Consolidated from {len(unique_mzs)} to {len(consolidated_mzs)} "
            f"m/z values"
        )
        return consolidated_mzs

    def _consolidate_with_tolerance(
        self, mzs: np.ndarray, tolerance_ppm: float
    ) -> np.ndarray:
        """Consolidate m/z values within specified tolerance.

        Args:
            mzs: Array of m/z values
            tolerance_ppm: Tolerance in ppm

        Returns:
            Array of consolidated m/z values
        """
        if len(mzs) <= 1:
            return mzs

        consolidated = []
        i = 0

        while i < len(mzs):
            current_mz = mzs[i]
            group_mzs = [current_mz]

            # Find all m/z values within tolerance
            j = i + 1
            while j < len(mzs):
                mz_diff_ppm = abs(mzs[j] - current_mz) / current_mz * 1e6
                if mz_diff_ppm <= tolerance_ppm:
                    group_mzs.append(mzs[j])
                    j += 1
                else:
                    break

            # Use average of the group as representative
            consolidated.append(np.mean(group_mzs))
            i = j

        return np.array(consolidated)

    def get_stats(self) -> Dict[str, Any]:
        """Get statistics about the mass axis building process."""
        return self._stats.copy()

    def reset_stats(self) -> None:
        """Reset all statistics."""
        self._stats = {
            "total_spectra_processed": 0,
            "total_peaks_processed": 0,
            "unique_masses_found": 0,
            "memory_usage_mb": 0.0,
        }
