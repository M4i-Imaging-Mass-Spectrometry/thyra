"""Foundation tests for the Thyra raster chunk/shard policy.

Locks the current chunk tuples (so the future zarr v3 sharding wiring, gated on
spatialdata PR #1106, cannot silently change output) and documents the
cross-repo INNER_TILE_EDGE contract with Ousia's serving tile size.
"""

import pytest

from thyra.converters.spatialdata._chunking import (
    IMAGE_CHUNK_EDGE,
    INNER_TILE_EDGE,
    image_chunks,
)


def test_image_chunks_2d_reproduces_todays_literal():
    # No-op-refactor lock: the helper must return exactly the (1, 4096, 4096)
    # literal it replaced at base_spatialdata_converter.py's optical parse.
    assert image_chunks(2) == (1, 4096, 4096)


def test_image_chunks_3d_prepends_a_z_chunk():
    assert image_chunks(3) == (1, 1, 4096, 4096)


def test_image_chunks_rejects_other_ndim():
    with pytest.raises(ValueError):
        image_chunks(4)


def test_inner_tile_edge_matches_ousia_serving_contract():
    # Duplicated-by-contract with Ousia DEFAULT_TILE_SIZE = 512 (Thyra cannot
    # import Ousia). If Ousia's serving tile changes, update this in lockstep.
    assert INNER_TILE_EDGE == 512
    # Outer chunk must tile exactly into inner tiles (no ragged remainder).
    assert IMAGE_CHUNK_EDGE % INNER_TILE_EDGE == 0
