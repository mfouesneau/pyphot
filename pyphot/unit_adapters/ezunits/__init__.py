# -*- coding: utf-8 -*-
"""
pint
~~~~

Pint is Python module/package to define, operate and manipulate
**physical quantities**: the product of a numerical value and a
unit of measurement. It allows arithmetic operations between them
and conversions from and to different units.

:copyright: (c) 2012 by Hernan E. Grecco.
:license: BSD, see LICENSE for more details.
"""

from .pint import (
    UnitRegistry,
    DimensionalityError,
    UnitsContainer,
    UndefinedUnitError,
    logger,
    __version__,
)

# load a default registery.
## Example sage unit['m * s **-1']
unit = UnitRegistry()


def hasUnit(val):
    return hasattr(val, "units")


__all__ = [
    "UnitRegistry",
    "DimensionalityError",
    "UnitsContainer",
    "UndefinedUnitError",
    "unit",
    "hasUnit",
    "logger",
    "__version__",
]
