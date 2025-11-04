import os
import inspect
from importlib import resources
from .unit_adapters import find_default_units_backend, get_adapter as get_units_adapter

__all__ = [
    "__ROOT__",
    "__default_passband_lib__",
    "__default_lick_lib__",
    "__vega_default_flavor__",
    "units",
    "set_units_backend",
    "set_vega_flavor",
]

# directories
__ROOT__ = "/".join(
    os.path.abspath(inspect.getfile(inspect.currentframe())).split("/")[:-1]  # type: ignore
)
# default library directory
# libsdir = os.path.abspath(os.path.join(__ROOT__, '../libs/'))
libsdir = resources.files("pyphot").joinpath("libs")

# default passband library
__default_passband_lib__ = libsdir.joinpath("new_filters.hd5")
__default_lick_lib__ = libsdir.joinpath("licks.dat")

# Set default vega flavor
__vega_default_flavor__ = "stis_003"

# provide a unit adapter by default that can be imported
units = find_default_units_backend()


def set_units_backend(backend: str):
    """Set the units backend to use throughout pyphot.

    Parameters
    ----------
    backend : str
        The name of the units backend to use.

    see :py:mod:`pyphot.unit_adapters`
    """
    global units
    units = get_units_adapter(backend)


def set_vega_flavor(flavor: str):
    """Set the vega flavor to use throughout pyphot.

    Parameters
    ----------
    flavor : str
        The name of the vega flavor to use.

    see :py:mod:`pyphot.vega`
    """
    global __vega_default_flavor__
    __vega_default_flavor__ = flavor
