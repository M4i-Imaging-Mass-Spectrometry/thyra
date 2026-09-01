# thyra/readers/mzpeak/__init__.py
"""Experimental mzPeak MSI reader.

mzPeak (HUPO-PSI) is a Parquet-in-ZIP container intended as the successor to
mzML/imzML at the raw/archival layer. Thyra treats it as an *input* only:
archives convert to SpatialData through the normal pipeline, and Thyra never
writes one.

The format is a moving draft (v0.9 at the time of writing), so the reader
codes defensively and announces itself as experimental in its log output.
"""

from .mzpeak_reader import MzPeakReader

__all__ = [
    "MzPeakReader",
]
