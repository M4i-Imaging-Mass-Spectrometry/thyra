"""Bruker solariX (FT-ICR / MRMS) reader package.

Reads the processed ``peaks.sqlite`` store inside a solariX imaging ``.d``
directory. Pure Python -- no vendor SDK involved.
"""

from .solarix_reader import SolarixReader

__all__ = ["SolarixReader"]
