# thyra/readers/phi/__init__.py
"""PHI SmartSoft-TOF MSI reader implementation.

Reads ToF-SIMS imaging data from PHI (Physical Electronics) nanoTOF
instruments, stored as a single ``.raw`` file holding an ASCII acquisition
header followed by a block-structured stream of individual ion events.

No vendor SDK or native library is required.
"""

from .event_stream import BlockIndex, DataSpan, iter_event_batches, scan_blocks
from .mass_axis import PhiMassAxis, build_mass_axis
from .phi_header import Calibrant, PhiHeader, parse_phi_header
from .phi_reader import PhiReader

__all__ = [
    "BlockIndex",
    "Calibrant",
    "DataSpan",
    "PhiHeader",
    "PhiMassAxis",
    "PhiReader",
    "build_mass_axis",
    "iter_event_batches",
    "parse_phi_header",
    "scan_blocks",
]
