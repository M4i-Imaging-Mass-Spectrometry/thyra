"""Batch processing utilities for efficient handling of large datasets.

This module provides utilities for processing large datasets in
manageable batches with progress tracking and memory management.
"""

import logging
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterator, List, Optional

import numpy as np
from tqdm import tqdm

logger = logging.getLogger(__name__)


@dataclass
class BatchInfo:
    """Information about a processing batch."""

    batch_id: int
    start_index: int
    end_index: int
    size: int
    estimated_memory_mb: float


class BatchProcessor:
    """Efficient batch processor for large spectrum datasets.

    This class provides intelligent batching strategies based on memory
    constraints and processing requirements.
    """

    def __init__(
        self,
        target_memory_mb: float = 512,
        min_batch_size: int = 10,
        max_batch_size: int = 1000,
        progress_callback: Optional[Callable] = None,
    ):
        """Initialize the batch processor.

        Args:
            target_memory_mb: Target memory usage per batch
            min_batch_size: Minimum batch size
            max_batch_size: Maximum batch size
            progress_callback: Optional callback for progress updates
        """
        self.target_memory_mb = target_memory_mb
        self.min_batch_size = min_batch_size
        self.max_batch_size = max_batch_size
        self.progress_callback = progress_callback

        # Statistics
        self._stats = {
            "total_batches": 0,
            "total_items_processed": 0,
            "average_batch_size": 0.0,
            "average_processing_time_ms": 0.0,
            "peak_memory_mb": 0.0,
        }

        logger.info(
            f"Initialized BatchProcessor (target: {target_memory_mb}MB, "
            f"batch size: {min_batch_size}-{max_batch_size})"
        )

    def calculate_optimal_batch_size(
        self, total_items: int, avg_item_size_bytes: float
    ) -> int:
        """Calculate optimal batch size based on memory constraints.

        Args:
            total_items: Total number of items to process
            avg_item_size_bytes: Average size per item in bytes

        Returns:
            Optimal batch size
        """
        # Calculate memory per item in MB
        item_size_mb = avg_item_size_bytes / 1024 / 1024

        # Calculate batch size to stay within target memory
        if item_size_mb > 0:
            batch_size = int(self.target_memory_mb / item_size_mb)
        else:
            batch_size = self.max_batch_size

        # Apply bounds
        batch_size = max(self.min_batch_size, min(batch_size, self.max_batch_size))
        batch_size = min(batch_size, total_items)

        logger.debug(
            f"Calculated optimal batch size: {batch_size} "
            f"(item size: {item_size_mb:.3f}MB)"
        )

        return batch_size

    def create_batches(
        self, total_items: int, batch_size: Optional[int] = None
    ) -> List[BatchInfo]:
        """Create batch information for processing.

        Args:
            total_items: Total number of items to process
            batch_size: Optional fixed batch size

        Returns:
            List of BatchInfo objects
        """
        if batch_size is None:
            # Use a reasonable default if no size provided
            batch_size = min(
                self.max_batch_size,
                max(self.min_batch_size, total_items // 10),
            )

        batches = []
        batch_id = 0

        for start_idx in range(0, total_items, batch_size):
            end_idx = min(start_idx + batch_size, total_items)
            actual_size = end_idx - start_idx

            batch_info = BatchInfo(
                batch_id=batch_id,
                start_index=start_idx,
                end_index=end_idx,
                size=actual_size,
                estimated_memory_mb=0.0,  # Will be updated during processing
            )

            batches.append(batch_info)
            batch_id += 1

        logger.info(
            f"Created {len(batches)} batches for {total_items} items "
            f"(avg size: {total_items / len(batches):.1f})"
        )

        return batches

    def process_spectrum_batches(
        self,
        spectrum_iterator: Iterator,
        total_spectra: int,
        processor_func: Callable,
        batch_size: Optional[int] = None,
    ) -> Iterator[Any]:
        """Process spectra in batches with progress tracking.

        Args:
            spectrum_iterator: Iterator yielding spectrum data
            total_spectra: Total number of spectra
            processor_func: Function to process each batch
            batch_size: Optional batch size

        Yields:
            Results from processor_func for each batch
        """
        if batch_size is None:
            # Estimate batch size based on typical spectrum size
            avg_spectrum_size = 1000 * 12  # 1000 peaks * 12 bytes per peak
            batch_size = self.calculate_optimal_batch_size(
                total_spectra, avg_spectrum_size
            )

        batches = self.create_batches(total_spectra, batch_size)

        # Setup progress tracking
        pbar = tqdm(
            total=total_spectra,
            desc="Processing batches",
            unit="spectrum",
            disable=getattr(self, "_quiet_mode", False),
        )

        try:
            current_batch = []
            current_batch_idx = 0
            spectrum_count = 0

            for spectrum_data in spectrum_iterator:
                current_batch.append(spectrum_data)
                spectrum_count += 1

                # Check if batch is complete
                if len(current_batch) >= batch_size or spectrum_count >= total_spectra:
                    # Process the batch
                    if current_batch_idx < len(batches):
                        batch_info = batches[current_batch_idx]
                        batch_info.size = len(current_batch)
                    else:
                        batch_info = BatchInfo(
                            batch_id=current_batch_idx,
                            start_index=spectrum_count - len(current_batch),
                            end_index=spectrum_count,
                            size=len(current_batch),
                            estimated_memory_mb=0.0,
                        )

                    # Process batch
                    result = processor_func(current_batch, batch_info)
                    yield result

                    # Update statistics
                    self._stats["total_batches"] += 1
                    self._stats["total_items_processed"] += len(current_batch)

                    # Update progress
                    pbar.update(len(current_batch))

                    # Progress callback
                    if self.progress_callback:
                        self.progress_callback(spectrum_count, total_spectra)

                    # Reset for next batch
                    current_batch = []
                    current_batch_idx += 1

                # Break if we've processed all spectra
                if spectrum_count >= total_spectra:
                    break

        finally:
            pbar.close()

        # Update final statistics
        if self._stats["total_batches"] > 0:
            self._stats["average_batch_size"] = (
                self._stats["total_items_processed"] / self._stats["total_batches"]
            )

    def _calculate_initial_batch_size(
        self, total_spectra: int, initial_batch_size: Optional[int]
    ) -> int:
        """Calculate initial batch size."""
        if initial_batch_size is None:
            return min(50, max(10, total_spectra // 100))
        return initial_batch_size

    def _should_process_batch(
        self,
        current_batch: list,
        batch_size: int,
        spectrum_count: int,
        total_spectra: int,
    ) -> bool:
        """Check if batch should be processed."""
        return len(current_batch) >= batch_size or spectrum_count >= total_spectra

    def _adjust_batch_size_adaptively(
        self, batch_times: list, current_batch_size: int
    ) -> int:
        """Adjust batch size based on performance."""
        if len(batch_times) >= 3:
            avg_time = np.mean(batch_times[-3:])

            if avg_time < 0.5:  # Fast processing, increase batch size
                return min(current_batch_size + 10, self.max_batch_size)
            elif avg_time > 2.0:  # Slow processing, decrease batch size
                return max(current_batch_size - 10, self.min_batch_size)

        return current_batch_size

    def adaptive_batch_processing(
        self,
        spectrum_iterator: Iterator,
        total_spectra: int,
        processor_func: Callable,
        initial_batch_size: Optional[int] = None,
    ) -> Iterator[Any]:
        """Adaptive batch processing that adjusts batch size based on performance.

        Args:
            spectrum_iterator: Iterator yielding spectrum data
            total_spectra: Total number of spectra
            processor_func: Function to process each batch
            initial_batch_size: Initial batch size

        Yields:
            Results from processor_func for each batch
        """
        import time

        batch_size = self._calculate_initial_batch_size(
            total_spectra, initial_batch_size
        )

        pbar = tqdm(
            total=total_spectra,
            desc="Adaptive processing",
            unit="spectrum",
            disable=getattr(self, "_quiet_mode", False),
        )

        try:
            current_batch = []
            spectrum_count = 0
            batch_times: List[float] = []

            for spectrum_data in spectrum_iterator:
                current_batch.append(spectrum_data)
                spectrum_count += 1

                if self._should_process_batch(
                    current_batch, batch_size, spectrum_count, total_spectra
                ):
                    start_time = time.time()

                    batch_info = BatchInfo(
                        batch_id=len(batch_times),
                        start_index=spectrum_count - len(current_batch),
                        end_index=spectrum_count,
                        size=len(current_batch),
                        estimated_memory_mb=0.0,
                    )

                    result = processor_func(current_batch, batch_info)
                    yield result

                    batch_time = time.time() - start_time
                    batch_times.append(batch_time)

                    batch_size = self._adjust_batch_size_adaptively(
                        batch_times, batch_size
                    )

                    pbar.update(len(current_batch))

                    if self.progress_callback:
                        self.progress_callback(spectrum_count, total_spectra)

                    current_batch = []

                if spectrum_count >= total_spectra:
                    break

        finally:
            pbar.close()

    def process_with_memory_monitoring(
        self,
        items: List[Any],
        processor_func: Callable,
        memory_limit_mb: Optional[float] = None,
    ) -> Iterator[Any]:
        """Process items with memory monitoring and adaptive batch sizing.

        Args:
            items: List of items to process
            processor_func: Function to process each batch
            memory_limit_mb: Optional memory limit override

        Yields:
            Results from processor_func for each batch
        """
        # Setup memory monitoring
        process, memory_limit_mb = self._setup_memory_monitoring(memory_limit_mb)

        # Initialize processing variables
        total_items = len(items)
        batch_size = self.min_batch_size
        processed = 0

        pbar = tqdm(
            total=total_items,
            desc="Memory-aware processing",
            unit="item",
            disable=getattr(self, "_quiet_mode", False),
        )

        try:
            yield from self._process_items_with_monitoring(
                items,
                processor_func,
                process,
                memory_limit_mb,
                batch_size,
                processed,
                total_items,
                pbar,
            )
        finally:
            pbar.close()

    def _setup_memory_monitoring(self, memory_limit_mb):
        """Setup memory monitoring with psutil if available."""
        try:
            import psutil

            PSUTIL_AVAILABLE = True
        except ImportError:
            PSUTIL_AVAILABLE = False
            psutil = None

        if memory_limit_mb is None:
            memory_limit_mb = self.target_memory_mb

        if PSUTIL_AVAILABLE:
            process = psutil.Process()
        else:
            process = None

        return process, memory_limit_mb

    def _process_items_with_monitoring(
        self,
        items,
        processor_func,
        process,
        memory_limit_mb,
        batch_size,
        processed,
        total_items,
        pbar,
    ):
        """Process items with memory monitoring and adaptive sizing."""
        while processed < total_items:
            # Check and adjust batch size based on memory
            memory_mb = self._get_current_memory(process)
            batch_size = self._adjust_batch_size(memory_mb, memory_limit_mb, batch_size)

            # Create and process batch
            batch, end_idx = self._create_batch(
                items, processed, batch_size, total_items
            )
            if not batch:
                break

            result = self._process_single_batch(
                batch, processor_func, processed, end_idx, memory_mb
            )
            yield result

            # Update progress
            processed += len(batch)
            pbar.update(len(batch))
            self._update_memory_stats(memory_mb)

    def _get_current_memory(self, process):
        """Get current memory usage."""
        if process:
            return process.memory_info().rss / 1024 / 1024
        return 100.0  # Fallback value

    def _adjust_batch_size(self, memory_mb, memory_limit_mb, batch_size):
        """Adjust batch size based on memory usage."""
        if memory_mb > memory_limit_mb * 0.8:  # 80% of limit
            batch_size = max(self.min_batch_size, batch_size // 2)
            logger.warning(
                f"High memory usage ({memory_mb:.1f}MB), "
                f"reducing batch size to {batch_size}"
            )
        elif memory_mb < memory_limit_mb * 0.4:  # 40% of limit
            batch_size = min(self.max_batch_size, batch_size * 2)
        return batch_size

    def _create_batch(self, items, processed, batch_size, total_items):
        """Create a batch of items to process."""
        end_idx = min(processed + batch_size, total_items)
        batch = items[processed:end_idx]
        return batch, end_idx

    def _process_single_batch(
        self, batch, processor_func, processed, end_idx, memory_mb
    ):
        """Process a single batch."""
        batch_info = BatchInfo(
            batch_id=processed // len(batch) if batch else 0,
            start_index=processed,
            end_index=end_idx,
            size=len(batch),
            estimated_memory_mb=memory_mb,
        )
        return processor_func(batch, batch_info)

    def _update_memory_stats(self, memory_mb):
        """Update peak memory statistics."""
        if memory_mb > self._stats["peak_memory_mb"]:
            self._stats["peak_memory_mb"] = memory_mb

    def get_stats(self) -> Dict[str, Any]:
        """Get batch processing statistics."""
        return self._stats.copy()

    def reset_stats(self) -> None:
        """Reset all statistics."""
        self._stats = {
            "total_batches": 0,
            "total_items_processed": 0,
            "average_batch_size": 0.0,
            "average_processing_time_ms": 0.0,
            "peak_memory_mb": 0.0,
        }
