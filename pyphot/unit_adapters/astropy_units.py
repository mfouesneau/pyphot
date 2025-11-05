"""
Declare missing photometric and spectral units for use with astropy.
"""

from astropy.units import Unit, def_unit, add_enabled_units
from astropy.units import Quantity
from astropy.units.core import Unit as U

__all__ = ["Unit", "Quantity", "U"]

new_units = dict(
    flam="erg * s ** (-1) * AA ** (-1) * cm **(-2)",
    fnu="erg * s ** (-1) * Hz ** (-1) * cm **(-2)",
    photflam="photon * s ** (-1) * AA ** (-1) * cm **(-2)",
    photfnu="photon * s ** (-1) * Hz ** (-1) * cm **(-2)",
    angstroms="angstrom",
    lsun="Lsun",
    ergs="erg",
)

add_enabled_units([def_unit([k], Unit(v)) for k, v in new_units.items()]).__enter__()
