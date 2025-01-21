""" Handle vega spec/mags/fluxes manipulations

Works with both ascii and hd5 files for back-compatibility

Vega.wavelength and Vega.flux have now units!
"""
from __future__ import print_function
from functools import wraps
import numpy
from .config import libsdir
from ..simpletable import SimpleTable
from astropy.units import Unit


__all__ = ['Vega', 'from_Vegamag_to_Flux', 'from_Vegamag_to_Flux_SN_errors']

_default_vega = "{0}/vega.hd5".format(libsdir)


class Vega(object):
    """
    Class that handles vega spectrum and references.  This class know where to
    find the Vega synthetic spectrum (Bohlin 2007) in order to compute fluxes
    and magnitudes in given filters

    Attributes
    ----------
    source: str
        filename of the vega library
    data: SimpleTable
        data table
    units: tuple
        detected units from file header
    wavelength: array
        wavelength (with units when found)
    flux: array
        flux(wavelength) values (with units when provided)

    An instance can be used as a context manager as:

    >>> filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_WFC3_F475W',\
                   'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
        with Vega() as v:
            vega_f, vega_mag, flamb = v.getSed(filters)
        print vega_f, vega_mag, flamb
    """

    def __init__(self, source=_default_vega):
        """ Constructor """
        self.source = source
        self.data = None
        self.units = None

    def _readfile(self):
        if self.data is not None:
            return
        fname = self.source
        # get extension
        ext = fname.split('.')[-1]
        if ( ext.lower() in ('hd5', 'hdf', 'hdf5') ):
            self.data = SimpleTable(fname, '/spectrum', silent=True)
        else:
            self.data = SimpleTable(fname, silent=True)
        try:
            uw = self.data.header['WAVELENGTH_UNIT'].split('=')[0].rstrip()
            uf = self.data.header['FLUX_UNIT'].split('=')[0].rstrip()
            self.units = uw, uf
        except TypeError:
            uw = self.data.header['WAVELENGTH_UNIT'].split(b'=')[0].decode().rstrip()
            uf = self.data.header['FLUX_UNIT'].split(b'=')[0].decode().rstrip()
            self.units = uw, uf

    def __enter__(self):
        """ Enter context """
        self._readfile()
        return self

    def __exit__(self,  *exc_info):
        """ end context """
        return False

    @property
    def wavelength(self):
        """ wavelength (with units when found) """
        self._readfile()
        try:
            return self.data.WAVELENGTH * Unit(self.units[0].lower())
        except Exception:
            return self.data.WAVELENGTH

    @property
    def flux(self):
        """ flux(wavelength) values (with units when provided) """
        self._readfile()
        try:
            return self.data.FLUX * Unit(self.units[1].lower())
        except Exception:
            return self.data.FLUX

    def getFlux(self, filters):
        """ Return vega abs. fluxes in filters """
        self._readfile()
        w = self.wavelength.to('AA').magnitude
        f = self.flux.magnitude
        r = numpy.array([k.getFlux(w, f) for k in filters])
        return r

    def getMag(self, filters):
        """ Return vega abs. magnitudes in filters """
        return -2.5 * numpy.log10(self.getFlux(filters))


def from_Vegamag_to_Flux(lamb, vega_mag):
    """ function decorator that transforms vega magnitudes to fluxes (without vega reference) """
    def deco(f):
        def vegamagtoFlux(mag, err, mask):
            f = numpy.power(10, -0.4 * (mag + vega_mag))
            e = f * ( 1. - numpy.power(10, -0.4 * err) )
            return f, e, mask

        @wraps(f)
        def wrapper(*args, **kwargs):
            mag, err, mask = f(args[0], args[1], **kwargs)
            return vegamagtoFlux( mag, err, mask )

        return wrapper
    return deco


def from_Vegamag_to_Flux_SN_errors(lamb, vega_mag):
    """ function decorator that transforms vega magnitudes to fluxes (without vega reference) """
    def deco(f):
        def vegamagtoFlux(mag, errp, errm, mask):
            f = 10 ** (-0.4 * (mag + vega_mag))
            fp = 10 ** (-0.4 * (mag - errp + vega_mag))
            fm = 10 ** (-0.4 * (mag + errm + vega_mag))
            return f, fp - f, f - fm, mask

        @wraps(f)
        def wrapper(*args, **kwargs):
            mag, errp, errm, mask = f(args[0], args[1], **kwargs)
            return vegamagtoFlux( mag, errp, errm, mask )

        return wrapper
    return deco


def testUnit():
    """ Unit test and example usage """
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_WFC3_F475W',
               'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    with Vega() as v:
        vega_f, vega_mag, flamb = v.getSed(filters)
    print(vega_f, vega_mag, flamb)
