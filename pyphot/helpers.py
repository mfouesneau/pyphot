"""Internal helpers module of pyphot

This module provides internal shortcuts for pyphot code.
"""

import warnings
from functools import wraps
from math import pi

from .pbar import Pbar

# this is used to convert from bolometric luminosities to abs fluxes
# object to 10parsecs -- abs mag.
distc = 4.0 * pi * (3.0856775e19) ** 2


def progress_enumerate(it, *args, **kwargs):
    """Enumerate over a sequence with progression if requested

    Parameter
    ---------
    show_progress: bool
        set to show progress
    """
    progress = kwargs.pop("show_progress", False)
    if progress is True:
        yield from enumerate(Pbar(**kwargs).iterover(it), *args)
    else:
        yield from enumerate(it, *args)


def deprecated(message):
    """Deprecated warning decorator"""

    def decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            warnings.warn(message, DeprecationWarning, stacklevel=2)
            return function(*args, **kwargs)

        return wrapper

    return decorator
