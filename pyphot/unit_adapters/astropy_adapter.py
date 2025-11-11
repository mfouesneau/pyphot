"""
Adapter to use transparently with the internal ezunits (frozen pint version 0.1; pyphot legacy)
"""

from .units_adapter import UnitsAdapter, UnitTyping, raise_warning
from typing import Any, Union, Optional

from astropy.units import Quantity as ap_Quantity
from astropy.units.core import UnitConversionError
from .astropy_units import Unit as ap_Unit


class ApUnitsAdapter(UnitsAdapter):
    """Adapter for ezunits"""

    typing = UnitTyping(
        DimensionalityError=UnitConversionError,
        UndefinedUnitError=ValueError,
        Quantity=ap_Quantity,
    )

    @staticmethod
    def U(*args, **kwargs) -> ap_Unit:
        """Quantity from Unit Registry"""
        return ap_Unit(*args, **kwargs)

    @staticmethod
    def Q(*args, **kwargs) -> ap_Quantity:
        """Force quantity from Unit Registry"""
        return 1.0 * ap_Unit(*args, **kwargs)

    @staticmethod
    def has_unit(val: Union[Any, ap_Quantity]) -> bool:
        """Check if val is a value with unit information"""
        try:
            return isinstance(val, ap_Quantity)
        except TypeError:
            return hasattr(val, "units") or hasattr(val, "unit")

    @staticmethod
    def val_in_unit(
        varname: str,
        value: Union[Any, ap_Quantity],
        defaultunit: Optional[str] = None,
        warn: bool = True,
    ) -> ap_Quantity:
        """check units and convert to defaultunit or create the unit information

        Parameters
        ----------
        varname: str
            name of the variable
        value: object
            value of the variable, which may be unitless
        defaultunit: str
            default units is unitless
        warn: bool
            whether to warn if the variable does not have explicit units

        Returns
        -------
        quantity: ez_Quantity
            value with units

        Example
        -------
        >>> r = 0.5
        >>> print(val_in_unit('r', r, 'degree'))
        # UserWarning: Variable r does not have explicit units. Assuming `degree`
        <Quantity(0.5, 'degree')>
        >>> r = 0.5 * unit['degree']
        >>> print(val_in_unit('r', r, 'degree'))
        <Quantity(0.5, 'degree')>
        """
        if defaultunit is None:
            return value
        if not ApUnitsAdapter.has_unit(value):
            if warn:
                msg = "Variable {0:s} does not have explicit units. Assuming `{1:s}`\n"
                raise_warning(msg.format(varname, defaultunit), stacklevel=4)
            return value * ap_Unit(defaultunit)
        else:
            return value.to(defaultunit)
