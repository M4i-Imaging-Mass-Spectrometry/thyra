# thyra/readers/__init__.py
"""MSI data readers for various instrument formats.

This package provides reader implementations for different MSI data formats:
- ImzML: Open format for MSI data
- mzPeak: HUPO-PSI Parquet-in-ZIP archive (experimental, read-only)
- Bruker: timsTOF, Rapiflex and solariX data
- Waters: MassLynx .raw imaging data
- PHI: SmartSoft-TOF .raw ToF-SIMS data

Reader organization:
- bruker/: All Bruker formats (timsTOF, Rapiflex, solariX)
- imzml/: ImzML format reader
- mzpeak/: mzPeak archive reader (pyarrow, lazily imported)
- waters/: Waters .raw format reader (native DLL via ctypes)
- phi/: PHI SmartSoft-TOF .raw reader (pure Python, no vendor SDK)

Note that ``.raw`` is claimed by two vendors: Waters stores a directory,
PHI a single file. See :meth:`thyra.core.registry.MSIRegistry.detect_format`.
"""

# Import readers to trigger registration
from . import bruker, imzml, mzpeak, phi, waters  # noqa: F401
