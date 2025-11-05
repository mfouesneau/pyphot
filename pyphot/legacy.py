"""Defines legacy aliases"""

from .phot import Filter as UnitFilter
from .libraries import (
    Library as UnitLibrary,
    HDF_Library as UnitHDFLibrary,
    Ascii_Library as UnitAsciiLibrary,
)

__all__ = [
    "UnitFilter",
    "UnitLibrary",
    "UnitHDFLibrary",
    "UnitAsciiLibrary",
]
