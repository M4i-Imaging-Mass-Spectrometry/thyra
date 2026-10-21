"""Which reader options exist, and which formats use each one.

One table both routes read: ``convert_msi`` refuses a ``reader_options`` key
that is not here and names a key the detected format does not use, and the
CLI derives its "is ignored" lines from the same entries (issues #260,
#519). A module of its own, with no imports, because the CLI reads it at
import time and ``thyra --version`` must not load the converter stack
(issue #381).
"""

from typing import Dict, Optional, Tuple

#: Every ``reader_options`` key ``convert_msi`` documents, and the detected
#: formats whose reader uses it; ``None`` means every format. ``convert_msi``
#: refuses a key that is not here and names, in one WARNING, a key the
#: detected format does not use; the CLI derives its "is ignored" table from
#: the same entries, so the two routes say the same thing (issue #519).
#: ``region`` is the exception to the warning: a reader that cannot select
#: a region refuses it (D44), because dropping it converts every region.
READER_OPTION_FORMATS: Dict[str, Optional[Tuple[str, ...]]] = {
    "intensity_threshold": None,
    "use_recalibrated_state": ("bruker",),
    "tdf_spectrum": ("bruker",),
    "region": None,
    "use_centroid": ("waters",),
    "spectrum_type": ("imzml",),
    "max_mass_axis_length": ("imzml", "mzpeak", "solarix", "waters"),
    "use_appended_calibration": ("phi",),
}

#: How each detected format is spelled in a message. ``detect_format``
#: returns a registry key, and "a imzml source" is not what anyone calls
#: the file they handed over.
FORMAT_NAMES = {
    "imzml": "imzML",
    "mzpeak": "mzPeak",
    "bruker": "Bruker timsTOF (.d)",
    "solarix": "Bruker solariX (.d)",
    "rapiflex": "Bruker Rapiflex",
    "waters": "Waters .raw",
    "phi": "PHI SmartSoft-TOF",
}


def format_name(input_format: str) -> str:
    """How to spell a detected format in a message."""
    return FORMAT_NAMES.get(input_format, input_format)


def reader_option_reaches(key: str, input_format: str) -> bool:
    """Whether the reader for ``input_format`` uses the option ``key``."""
    formats = READER_OPTION_FORMATS.get(key)
    return formats is None or input_format in formats
