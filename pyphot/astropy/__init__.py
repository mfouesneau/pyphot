#!/usr/bin/env python
# -*- coding: utf-8 -*-
from astropy.units import Unit, def_unit, add_enabled_units

new_units = dict(
    flam='erg * s ** (-1) * AA ** (-1) * cm **(-2)',
    fnu='erg * s ** (-1) * Hz ** (-1) * cm **(-2)',
    photflam='photon * s ** (-1) * AA ** (-1) * cm **(-2)',
    photfnu='photon * s ** (-1) * Hz ** (-1) * cm **(-2)',
    angstroms='angstrom'
)

add_enabled_units([def_unit([k], Unit(v)) for k, v in new_units.items()])\
    .__enter__()

from .vega import Vega # noqa: E402
from .sun import Sun  # noqa: E402
from .sandbox import (UncertainFilter, UnitAscii_Library, UnitFilter,  # noqa: E402
                      UnitHDF_Library, UnitLibrary, UnitLickIndex,
                      UnitLickLibrary, get_library)

__all__ = [Vega, Sun, UncertainFilter, UnitAscii_Library, UnitFilter, UnitHDF_Library,
           UnitLibrary, UnitLickIndex, UnitLickLibrary, get_library]
 