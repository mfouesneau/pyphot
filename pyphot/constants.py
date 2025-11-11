"""
A module containing fundamental constants used in photometry calculations.

This module provides access to fundamental constants used in photometry calculations while ensuring consistency with the units defined in the :mod:`config` module. (not all unit packages provide these constants)
"""

from dataclasses import dataclass
from typing import Any, Optional

from . import config

__all__ = ["constants"]


@dataclass
class ConstantDef:
    """Definition of a constant used in photometry calculations."""

    value: Any
    """The value of the constant."""
    unit: str
    """The unit of the constant."""
    long_name: Optional[str] = None
    """The long name of the constant."""
    source: Optional[str] = None
    """The source of the constant."""


@dataclass
class Constants:
    """Constants used in photometry calculations.

    This class provides access to fundamental constants used in photometry calculations while ensuring consistency with the units defined in the `config` module.
    (not all unit packages provide these constants)
    """

    _definitions = {
        "c": ConstantDef(
            299792458.0,
            "m * s**(-1)",
            "Speed of light",
            "CODATA 2018",
        ),
        "h": ConstantDef(
            6.62607015e-34,
            "J*s",
            "Planck constant",
            "CODATA 2018",
        ),
    }

    @classmethod
    def get(cls, name: str) -> Any:
        """the value of a constant from the definitions

        parameters
        ----------
        name : str
            The name of the constant.
        returns
        -------
        float
            The value of the constant.
        """
        def_ = cls._definitions[name]
        return def_.value * config.units.U(def_.unit)

    def __new__(cls):
        """Declare all constants as properties

        This ensures consistency with the unit backedn defined in the `config` module.
        """
        for name in cls._definitions:
            exec(f"cls.{name} = property(lambda x: Constants.get('{name}'))")
        obj = super().__new__(cls)
        return obj


constants = Constants()
"""Constants used in photometry calculations.

This class provides access to fundamental constants used in photometry calculations while ensuring consistency with the units defined in the `config` module.
"""
