"""A fingerprint of a store's mass axis (issue #444, design decision D42).

This module imports nothing from a reader or a converter (D25), and numpy
only inside the function, so ``import thyra.metadata.schema`` stays light.
"""

import hashlib
from typing import Any


def mass_axis_sha256(mz: Any) -> str:
    """Hash a mass axis so two stores can be compared without loading both.

    The hash is SHA-256 over the values as little-endian float64 bytes,
    whatever dtype the array is held or was written in and whatever the
    platform's byte order. The same values therefore always give the same
    digest: a float32 copy, a big-endian copy or a non-contiguous view of an
    axis hash like the axis itself. Values that differ by one ulp give a
    different digest, by design.

    Args:
        mz: The m/z values of a store's main MSI table (``var["mz"]``), any
            array-like of numbers.

    Returns:
        The lower-case hexadecimal digest, 64 characters.
    """
    import numpy as np

    return hashlib.sha256(np.ascontiguousarray(mz, dtype="<f8").tobytes()).hexdigest()
