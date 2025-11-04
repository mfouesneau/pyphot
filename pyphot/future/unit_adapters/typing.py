from __future__ import annotations
from typing import Union, TypeAlias


# Add each quantity type if available
try:
    from . import astropy_adapter
except ImportError:
    pass

try:
    from . import pint_adapter
except ImportError:
    pass

try:
    from . import ezunits_adapter
except ImportError:
    pass

from . import units_adapter


QuantityType = Union[
    "astropy_adapter.ap_Quantity",
    "pint_adapter.pint_Quantity",
    "ezunits_adapter.ez_Quantity",
]
#
# For UnitAdapterType, enumerate the actual adapter classes
UnitAdapterType = Union[
    "astropy_adapter.ApUnitsAdapter",
    "pint_adapter.PintUnitsAdapter",
    "ezunits_adapter.EzUnitsAdapter",
    "units_adapter.UnitsAdapter",
]

DimensionalityError = Union[
    "astropy_adapter.UnitConversionError",
    "pint_adapter.DimensionalityError",
    "ezunits_adapter.DimensionalityError",
    ValueError,
]

UndefinedUnitError = Union[
    "pint_adapter.UndefinedUnitError",
    "ezunits_adapter.DimensionalityError",
    ValueError,
]
