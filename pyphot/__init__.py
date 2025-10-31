from typing import Optional

from pyphot.future.unit_adapters.units_adapter import UnitsAdapter

from .licks import LickLibrary, LickIndex, reduce_resolution

from .sandbox import (
    UnitLibrary,
    UnitAscii_Library,
    UnitHDF_Library,
    UnitFilter,
    UnitLickIndex,
    UnitLickLibrary,
)

from .sandbox import get_library
from .ezunits import unit

from .vega import Vega
from .sun import Sun
from .version import __VERSION__

__all__ = [
    LickLibrary,
    LickIndex,
    reduce_resolution,
    UnitLibrary,
    UnitAscii_Library,
    UnitHDF_Library,
    UnitFilter,
    UnitLickIndex,
    UnitLickLibrary,
    get_library,
    unit,
    Vega,
    Sun,
    __VERSION__,
]
