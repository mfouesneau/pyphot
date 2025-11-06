from . import typing
from . import legacy
from .version import __VERSION__

from .phot import Filter
from .libraries import Library, HDF_Library, Ascii_Library
from .licks import LickLibrary, LickIndex, reduce_resolution
from .vega import Vega
from .sun import Sun

__all__ = [
    "Ascii_Library",
    "Filter",
    "HDF_Library",
    "Library",
    "LickIndex",
    "LickLibrary",
    "Sun",
    "Vega",
    "__VERSION__",
    "legacy",
    "reduce_resolution",
    "typing",
]
