"""
Pyphot - set of tools to compute synthetic photometry in a simple way, suitable for integration in larger projects.

This package is mostly organized around 2 main classes:

* :class:`~pyphot.phot.Filter` that handles filter definitions and manipulations.

* :class:`~pyphot.libraries.Library` that handles a collection of filters in a few formats. This library class is derived into
  - an ASCII file reader :class:`~pyphot.libraries.Ascii_Library`,
  - and a single HDF file reader :class:`~pyphot.libraries.HDF_Library`.

In addition,

* :class:`~pyphot.vega.Vega` provides an interface to a synthetic reference of Vega. See also :doc:`vega`.

* :class:`~pyphot.sun.Sun` provides an interface to a synthetic and empirical reference of the Sun spectrum.

* :class:`~pyphot.licks` provides an extension to computing lick indices (:class:`~pyphot.licks.LickIndex`, :class:`~pyphot.licks.LickLibrary`).

* :func:`pyphot.svo.get_pyphot_filter` provides an interface to the `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_ to download filters directly from the SVO database.

* :mod:`~pyphot.typing` provides a set of type hints for the pyphot package.

"""

from . import typing
from . import legacy

from .libraries import Library, HDF_Library, Ascii_Library, get_library
from .licks import LickLibrary, LickIndex, reduce_resolution
from .phot import Filter
from .sun import Sun
from .vega import Vega
from .version import __VERSION__

__all__ = [
    "Ascii_Library",
    "Filter",
    "HDF_Library",
    "Library",
    "LickIndex",
    "LickLibrary",
    "Sun",
    "Vega",
    "__VERSION__",
    "get_library",
    "legacy",
    "reduce_resolution",
    "typing",
]
