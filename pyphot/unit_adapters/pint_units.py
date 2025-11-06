"""Remerging ezunits to the pint library"""

from pint import UnitRegistry
from pint import Quantity
from pint.errors import DimensionalityError, UndefinedUnitError

__all__ = ["unit", "Quantity", "UndefinedUnitError", "DimensionalityError"]

new_units = """
    #Astronomy related
    lsun = 3.828e26 * watt = Lsun
    msun = 1.98892e30 * kilogram = Msun
    rsun = 6.955e8 * meter = Rsun
    stephan_constant = 5.67037321e-8 * watt * meter ** -2 * kelvin ** -4
    magnitude = [luminosity] = mag
    flam = erg * s ** (-1) * AA ** (-1) * cm **(-2)
    fnu = erg * s ** (-1) * Hz ** (-1) * cm **(-2)
    Jy = 1e-26 * W * m ** (-2) * Hz ** (-1) = jansky
    photon = []
    photflam = photon * s ** (-1) * AA ** (-1) * cm **(-2)
    photfnu = photon * s ** (-1) * Hz ** (-1) * cm **(-2)


    #useful aliases
    micron = micrometer
    microns = micrometer
    AA = angstrom
""".splitlines()

unit = UnitRegistry()
unit.load_definitions(new_units)

# monkey patch pint to have a value and unit attribute
setattr(Quantity, "value", property(lambda self: self.magnitude))
setattr(Quantity, "unit", property(lambda self: self.units))
