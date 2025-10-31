"""
Adapter to use transparently with the internal ezunits (frozen pint version 0.1; pyphot legacy)
"""

from typing import Any, Union, NewType, Optional, cast

from pint import Quantity as pint_Quantity
from pint.errors import DimensionalityError, UndefinedUnitError

from .units_adapter import UnitsAdapter, UnitTyping, raise_warning
from .pint_units import unit as pint_Unit


class PintUnitsAdapter(UnitsAdapter):
    """Adapter for ezunits"""

    typing = UnitTyping(
        DimensionalityError=DimensionalityError,
        UndefinedUnitError=UndefinedUnitError,
        Quantity=pint_Quantity,
    )

    @staticmethod
    def U(*args, **kwargs) -> pint_Quantity:
        """Quantity from Unit Registry"""
        return cast(pint_Quantity, pint_Unit(*args, **kwargs))

    @staticmethod
    def Q(*args, **kwargs) -> pint_Quantity:
        """Quantity from Unit Registry"""
        return cast(pint_Quantity, pint_Unit(*args, **kwargs))

    @staticmethod
    def has_unit(val: Union[Any, pint_Quantity]) -> bool:
        """Check if val is a value with unit information"""
        try:
            return isinstance(val, pint_Quantity)
        except TypeError:
            return hasattr(val, "units") or hasattr(val, "unit")

    @staticmethod
    def val_in_unit(
        varname: str,
        value: Union[Any, pint_Quantity],
        defaultunit: Optional[str] = None,
        warn: bool = True,
    ) -> pint_Quantity:
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
        quantity: pint_Quantity
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
        if not PintUnitsAdapter.has_unit(value):
            if warn:
                msg = "Variable {0:s} does not have explicit units. Assuming `{1:s}`\n"
                raise_warning(msg.format(varname, defaultunit), stacklevel=4)
            return value * PintUnitsAdapter.U(defaultunit)
        else:
            return value.to(defaultunit)
