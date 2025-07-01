"""Handle vega spec/mags/fluxes manipulations

Works with both ascii and hd5 files for back-compatibility

Vega.wavelength and Vega.flux have now units!
"""

from __future__ import print_function

import warnings
from functools import wraps

import numpy

from .config import libsdir, __vega_default_flavor__
from .ezunits import unit
from .simpletable import SimpleTable

__all__ = ["Vega", "from_Vegamag_to_Flux", "from_Vegamag_to_Flux_SN_errors"]


_default_vega = {
    "mod_002": "{0}/alpha_lyr_mod_002.fits".format(libsdir),
    "mod_003": "{0}/alpha_lyr_mod_003.fits".format(libsdir),
    "mod_004": "{0}/alpha_lyr_mod_004.fits".format(libsdir),
    "stis_011": "{0}/alpha_lyr_stis_011.fits".format(libsdir),
    "stis_003": "{0}/alpha_lyr_stis_003.fits".format(libsdir),
    "legacy": "{0}/vega.hd5".format(libsdir),
}


class Vega(object):
    """
    Class that handles vega spectrum and references.  This class know where to
    find the Vega synthetic spectrum in order to compute fluxes
    and magnitudes in given filters

    Default Vega spectrum is Bohlin 2007, alpha_Lyr_mod_003.fits
    from the HST CDBS database, which is a synthetic spectrum of Vega

    Attributes
    ----------
    source: str
        filename of the vega library
    data: SimpleTable
        data table
    units: tuple
        detected units from file header
    flavor: str, (default theoretical)
        ["mod_002", "mod_003", "mod_004"] or ["stis_011"]

    An instance can be used as a context manager as:

    >>> filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_WFC3_F475W',\
                   'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
        with Vega() as v:
            vega_f, vega_mag, flamb = v.getSed(filters)
        print vega_f, vega_mag, flamb
    """

    def __init__(self, source=None, flavor=__vega_default_flavor__):
        """Constructor"""
        self.data = None
        self.units = None
        self._set_source_flavor(source, flavor)

    def _set_source_flavor(self, source=None, flavor=__vega_default_flavor__):
        """Set the source and flavor of the Vega spectrum"""
        if source is not None:
            self.source = source
        else:
            if flavor.lower() not in _default_vega:
                raise ValueError(
                    "Unknown Vega flavor: {0}. Available flavors {1}".format(
                        flavor, _default_vega.keys()
                    )
                )
            self.source = _default_vega[flavor]
        self.data = None
        self.units = None

    def _readfile(self):
        if self.data is not None:
            return
        fname = self.source
        # get extension
        ext = fname.split(".")[-1]
        if ext.lower() in ("hd5", "hdf", "hdf5"):
            self.data = SimpleTable(fname, "/spectrum", silent=True)
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.data = SimpleTable(fname, silent=True)

        try:
            self.data.header["WAVELENGTH_UNIT"] = self.data._units["WAVELENGTH"]
            self.data.header["FLUX_UNIT"] = self.data._units["FLUX"]
            self.units = self.data._units["WAVELENGTH"], self.data._units["FLUX"]
        except AttributeError:
            pass
        except KeyError:
            pass

        try:
            uw = self.data.header["WAVELENGTH_UNIT"].split("=")[0].rstrip()
            uf = self.data.header["FLUX_UNIT"].split("=")[0].rstrip()
            self.units = uw, uf
        except TypeError:
            uw = self.data.header["WAVELENGTH_UNIT"].split(b"=")[0].decode().rstrip()
            uf = self.data.header["FLUX_UNIT"].split(b"=")[0].decode().rstrip()
            self.units = uw, uf

    def __enter__(self):
        """Enter context"""
        self._readfile()
        return self

    def __exit__(self, *exc_info):
        """end context"""
        return False

    @property
    def wavelength(self):
        """wavelength (with units when found)"""
        self._readfile()
        try:
            return self.data.WAVELENGTH * unit[self.units[0].lower()]
        except Exception:
            return self.data.WAVELENGTH

    @property
    def flux(self):
        """flux(wavelength) values (with units when provided)"""
        self._readfile()
        try:
            return self.data.FLUX * unit[self.units[1].lower()]
        except Exception:
            return self.data.FLUX

    def getFlux(self, filters):
        """Return vega abs. fluxes in filters"""
        self._readfile()
        w = self.wavelength.to("AA").magnitude
        f = self.flux.magnitude
        r = numpy.array([k.getFlux(w, f) for k in filters])
        return r

    def getMag(self, filters):
        """Return vega abs. magnitudes in filters"""
        return -2.5 * numpy.log10(self.getFlux(filters))


def from_Vegamag_to_Flux(lamb, vega_mag):
    """function decorator that transforms vega magnitudes to fluxes (without vega reference)"""

    def deco(f):
        def vegamagtoFlux(mag, err, mask):
            f = numpy.power(10, -0.4 * (mag + vega_mag))
            e = f * (1.0 - numpy.power(10, -0.4 * err))
            return f, e, mask

        @wraps(f)
        def wrapper(*args, **kwargs):
            mag, err, mask = f(args[0], args[1], **kwargs)
            return vegamagtoFlux(mag, err, mask)

        return wrapper

    return deco


def from_Vegamag_to_Flux_SN_errors(lamb, vega_mag):
    """function decorator that transforms vega magnitudes to fluxes (without vega reference)"""

    def deco(f):
        def vegamagtoFlux(mag, errp, errm, mask):
            f = 10 ** (-0.4 * (mag + vega_mag))
            fp = 10 ** (-0.4 * (mag - errp + vega_mag))
            fm = 10 ** (-0.4 * (mag + errm + vega_mag))
            return f, fp - f, f - fm, mask

        @wraps(f)
        def wrapper(*args, **kwargs):
            mag, errp, errm, mask = f(args[0], args[1], **kwargs)
            return vegamagtoFlux(mag, errp, errm, mask)

        return wrapper

    return deco


def testUnit():
    """Unit test and example usage"""
    filters = [
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_WFC3_F475W",
        "HST_WFC3_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]
    with Vega() as v:
        vega_f, vega_mag, flamb = v.getSed(filters)
    print(vega_f, vega_mag, flamb)
