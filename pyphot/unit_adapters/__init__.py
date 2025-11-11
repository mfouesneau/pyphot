"""
This module provides a framework for using different unit packages such as astropy, pint, and others.
It provides a unified interface for working with units across different packages through adapters.

Registered Adapters
~~~~~~~~~~~~~~~~~~~

We provide the following adapters:

- :class:`~pyphot.unit_adapters.AstropyAdapter`: Adapter for the astropy unit system.
- :class:`~pyphot.unit_adapters.PintAdapter`: Adapter for the pint unit system.
- :class:`~pyphot.unit_adapters.EzUnitsAdapter`: Adapter for the legacy pyphot unit system based on pint (v0.1.0).

Note: The :class:`~pyphot.unit_adapters.EzUnitsAdapter` is provided as a fallback for older codes and in case other unit packages are missing.

API Overview
~~~~~~~~~~~~

All adapters must implement the :class:`~pyphot.unit_adapters.units_adapter.UnitsAdapter` interface, which defines the following methods:

- :meth:`~pyphot.unit_adapters.units_adapter.UnitsAdapter.U()`: queries the unit registry.
- :meth:`~pyphot.unit_adapters.units_adapter.UnitsAdapter.Q()`: make sure the returned object is a quantity, i.e. a value with unit (e.g. astropy differs between `U("kg")` and `Q("kg") = 1. * U("kg")`)

For source typing reasons, adapters also override a few methods:

- :meth:`~pyphot.unit_adapters.units_adapter.UnitsAdapter.has_unit`: checks if the object has a unit.
- :meth:`~pyphot.unit_adapters.units_adapter.UnitsAdapter.val_in_unit`: returns the value of the quantity in the specified unit or forces a default unit if not specified.

Decorator to impose units on arguments and keyword arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For convenience, adapters provide (by inheritance) a `decorate` decorator which can help defining functions with default units on arguments and keyword arguments.

For a more general decorator, use :func:`enforce_default_units` which will assumes the defined global adapter, if not a specific one provided.

.. note::
    Both decoration methods will update the docstring of the decorated function to include information about the imposed units.

Typing
~~~~~~

Hint typing information is provided at the module level: :class:`~pyphot.unit_adapters.UnitAdapterType`, :class:`~pyphot.unit_adapters.QuantityType`, which is a union of :class:`~pyphot.unit_adapters.units_adapter.UnitsAdapter` and varied :class:`~pyphot.unit_adapters.Quantity` definitions.
In addition, each adapter provides specific types in their `typing.Quantity` attribute.

Exceptions
~~~~~~~~~~

For convenience adapters also provide a list of potential exception types. At the module level
:class:`~pyphot.unit_adapters.UnitsAdapter.UndefinedUnitError`, and :class:`~pyphot.unit_adapters.UnitsAdapter.DimensionalityError` which applies to any backend units.

When needed :class:`~pyphot.unit_adapters.UnitsAdapter.typing.UndefinedUnitError` and :class:`~pyphot.unit_adapters.UnitsAdapter.typing.DimensionalityError` can be used to catch specific errors related to unit conversion or registration by the unit package.
"""

from typing import Any
from collections import OrderedDict


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

from .units_adapter import enforce_default_units

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
