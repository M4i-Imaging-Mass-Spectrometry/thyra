"""Thresholds for resampling and instrument detection.

The controlled-vocabulary constants that used to live here --
:class:`ImzMLAccessions`, :class:`SpectrumType`, :class:`BinaryDataType`
and :func:`normalize_spectrum_type` -- moved to
:mod:`thyra.metadata.constants` under issue #381: they describe data, and
importing them should not mean importing a resampler. They are re-exported
below so ``from thyra.resampling.constants import SpectrumType`` and
``from thyra.resampling import SpectrumType`` keep resolving.
"""

from ..metadata.constants import (  # noqa: F401
    SPECTRUM_TYPE_ALIASES,
    BinaryDataType,
    ImzMLAccessions,
    SpectrumType,
    normalize_spectrum_type,
)


class Thresholds:
    """Threshold values for data classification."""

    # Peak density threshold for detecting profile data
    # Profile data typically has >5000 points per spectrum
    PROFILE_PEAK_DENSITY = 5000
