import logging
from typing import Optional, Tuple, cast

import dask.array as da
import zarr

logger = logging.getLogger(__name__)


def optimize_zarr_chunks(
    zarr_path: str,
    array_path: str,
    output_path: Optional[str] = None,
    chunks: Optional[Tuple[int, ...]] = None,
    compressor=None,
) -> bool:
    """Optimize the chunking of a Zarr array.

    Parameters:
    -----------
    zarr_path : str
        Path to the Zarr store.
    array_path : str
        Path to the array within the Zarr store.
    output_path : Optional[str]
        Path for the output Zarr store. If None, overwrites the input.
    chunks : Optional[Tuple[int, ...]]
        Desired chunk size. If None, optimizes based on array shape.
    compressor : Optional
        Compression to use. If None, uses the input's compressor.

    Returns:
    --------
    bool: True if optimization was successful, False otherwise.
    """
    try:
        # Open the Zarr store
        zarr_store = zarr.open_group(zarr_path, mode="r")
        array = cast(zarr.Array, zarr_store[array_path])

        if chunks is None:
            # For MSI data, chunk along m/z dimension and larger spatial tiles
            shape = array.shape
            if len(shape) == 4:  # (c, z, y, x)
                chunks = (
                    min(10000, shape[0]),
                    1,
                    min(64, shape[2]),
                    min(64, shape[3]),
                )

        # Use original compressor if not provided
        if compressor is None:
            compressor = array.compressor

        # Create Dask array
        dask_array = da.from_array(array, chunks=array.chunks)

        # Rechunk
        rechunked = dask_array.rechunk(chunks)

        # Determine output path
        if output_path is None:
            output_path = zarr_path
            output_array_path = f"{array_path}_optimized"
        else:
            output_array_path = array_path

        # Create output store
        output_store = zarr.open_group(output_path, mode="a")

        # Write rechunked array
        da.to_zarr(
            rechunked,
            output_store,
            component=output_array_path,
            compressor=compressor,
            overwrite=True,
        )

        # If overwriting input, rename arrays
        if output_path == zarr_path:
            # Create temporary backup
            output_store[f"{array_path}_backup"] = output_store[array_path]

            # Delete original and rename optimized
            del output_store[array_path]
            output_store[array_path] = output_store[output_array_path]
            del output_store[output_array_path]

            # Delete backup
            del output_store[f"{array_path}_backup"]

        return True
    except Exception as e:
        logger.error(f"Error optimizing chunks: {e}")
        return False
