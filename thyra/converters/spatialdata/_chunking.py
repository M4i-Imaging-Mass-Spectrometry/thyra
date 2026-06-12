"""Image chunk/shard policy for Thyra-written SpatialData rasters.

Single source of truth for how Thyra chunks the raster images it writes to
SpatialData zarr -- and the designated seam for tile-aligned zarr v3 sharding
once scverse/spatialdata PR #1106 (``raster_write_kwargs``) lands.

Cross-repo contract: the viewer (Ousia) morphology tile server streams
``DEFAULT_TILE_SIZE`` (= 512) edge tiles. Once sharding lands, the zarr INNER
chunk edge must equal that so one served tile == one decompress. Thyra is a
separate package and CANNOT import Ousia, so the coupling is DUPLICATED BY
CONTRACT here -- keep ``INNER_TILE_EDGE`` numerically in sync with Ousia's
``DEFAULT_TILE_SIZE``; a silent divergence reintroduces cold-tile amplification
on the viewer. See the Ousia ADR docs/ADR-pyramid-tile-sharding.md.
"""

from typing import Tuple

# Outer chunk edge used TODAY (pre-sharding the chunk IS the outer block): a 4096
# block means a viewer 512-tile read decompresses at most one chunk per request.
IMAGE_CHUNK_EDGE = 4096

# Future zarr INNER chunk edge once #1106 lands == Ousia DEFAULT_TILE_SIZE (512).
# Duplicated by contract (Thyra cannot import Ousia).
INNER_TILE_EDGE = 512


def image_chunks(spatial_ndim: int, edge: int = IMAGE_CHUNK_EDGE) -> Tuple[int, ...]:
    """The zarr chunk tuple for a Thyra-written raster image (leading channel axis).

    ``spatial_ndim == 2`` (c, y, x) -> ``(1, edge, edge)``;
    ``spatial_ndim == 3`` (c, z, y, x) -> ``(1, 1, edge, edge)``.

    This is the ONE place the future inner=``INNER_TILE_EDGE`` + outer-shard policy
    will live when spatialdata PR #1106 ships ``raster_write_kwargs``. Today it
    returns the pre-sharding behaviour, byte-for-byte.
    """
    if spatial_ndim == 2:
        return (1, edge, edge)
    if spatial_ndim == 3:
        return (1, 1, edge, edge)
    raise ValueError(f"spatial_ndim must be 2 or 3, got {spatial_ndim}")
