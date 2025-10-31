import os
import inspect
from dataclasses import dataclass
from importlib import resources


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
