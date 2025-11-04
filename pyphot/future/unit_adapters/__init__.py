"""
Checks which backends are available and provides a global typing hint for the units adapters.
"""

from dataclasses import dataclass
from typing import Any
from collections import OrderedDict, namedtuple
from functools import reduce
import operator


__all__ = [
    "get_adapter",
    "backends",
    "UnitAdapterType",
    "find_default_units_backend",
    "DimensionalityError",
    "UndefinedUnitError",
    "QuantityType",
    "enforce_default_units",
]

from .units_adapter import enforce_default_units, UnitsAdapter

try:
    from .ezunits_adapter import EzUnitsAdapter
except ImportError:
    EzUnitsAdapter = None

try:
    from .astropy_adapter import ApUnitsAdapter, ap_Quantity
except ImportError:
    ApUnitsAdapter = None
    ap_Quantity = Any

try:
    from .pint_adapter import PintUnitsAdapter
except ImportError:
    PintUnitsAdapter = None


from .typing import (
    QuantityType,
    UnitAdapterType,
    DimensionalityError,
    UndefinedUnitError,
)

# order provides backend priority for defaults
backends = OrderedDict(
    [
        ("astropy", ApUnitsAdapter),
        ("pint", PintUnitsAdapter),
        ("ezunits", EzUnitsAdapter),
    ]
)

# Provide tools to setup the units adapters


def get_adapter(name) -> UnitAdapterType:
    """Get the units adapter for the given name.

    Parameters
    ----------
    name: str
        The name of the units adapter.

    Returns
    -------
    UnitAdapter
        The units adapter for the given name.

    Raises
    ------
    ValueError: If the given name is not a registered adapter.
    """
    try:
        backend = backends[name]
    except KeyError:
        raise RuntimeError(
            f"Unknown units adapter: {name}. Should be one of {list(backends.keys())}"
        )
    if backend is None:
        raise ImportError(
            f"{name} backend requires a package to be installed. run `pip install {name}`."
        )
    return backend()


def find_default_units_backend() -> UnitAdapterType:
    """Set the default units adapter."""
    for backend in backends.values():
        if backend is not None:
            return backend()
    raise RuntimeError("No units adapter found.")
