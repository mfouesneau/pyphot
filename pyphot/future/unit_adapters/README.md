# Unit Adapters

This module provides a framework for using different unit packages such as astropy, pint, and others.
It provides a unified interface for working with units across different packages through adapters.

## Registered Adapters

We provide the following adapters:

- `AstropyAdapter`: Adapter for the astropy unit system.
- `PintAdapter`: Adapter for the pint unit system.
- `EzUnitsAdapter`: Adapter for the legacy pyphot unit system based on pint (v0.1.0).

Note: The `EzUnitsAdapter` is provided as a fallback for older codes and in case other unit packages are missing.

## API Overview

### General interface
All adapters must implement the `UnitsAdapter` interface, which defines the following methods:

- `UnitsAdapter.U()`: queries the unit registry.
- `UnitsAdapter.Q()`: make sure the returned object is a quantity, i.e. a value with unit (e.g. astropy differs between `U("kg")` and `Q("kg") = 1. * U("kg")`)

For source typing reasons, adapters also override a few methods:

- `UnitsAdapter.has_unit`: checks if the object has a unit.
- `UnitsAdapter.val_in_unit`: returns the value of the quantity in the specified unit or forces a default unit if not specified.

### Decorator to impose units on arguments and keyword arguments
For convenience, adapters provide (by inheritance) a `decorate` decorator which can help defining functions with default units on arguments and keyword arguments.

For a more general decorator, use `enforce_default_units` which will assumes the defined global adapter, if not a specific one provided.

Both decoration methods will update the docstring of the decorated function to include information about the imposed units.

### Typing
Hint typing information is provided at the module level: `UnitAdapterType`, `QuantityType`, which is a union of `UnitsAdapter` and varied `Quantity` definitions.
In addition, each adapter provides specific types in their `typing.Quantity` attribute.

### Exceptions

For convenience adapters also provide a list of potential exception types. At the module level
`UndefinedUnitError`, and `DimensionalityError` which applies to any backend units.

When needed `UnitsAdapter.typing.{UndefinedUnitError,DimensionalityError}` can be used to catch specific errors related to unit conversion or registration by the unit package.
