""" Handle the Sun Spectrum

"""
from __future__ import print_function
from functools import wraps
import numpy
from .config import libsdir
from ..simpletable import SimpleTable
from astropy.units import Unit
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits


__all__ = ['Sun']

_default_sun = {'observed': "{0}/sun_reference_stis_001.fits".format(libsdir),
                'theoretical': "{0}/sun_kurucz93.fits".format(libsdir)}
_default_distance = 1 * Unit('au')


class Sun(object):
    """
    Class that handles the Sun's spectrum and references.

    Observed solar spectrum comes from:
    ftp://ftp.stsci.edu/cdbs/current_calspec/sun_reference_stis_001.fits

    and theoretical spectrum comes from:
    ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/sun_kurucz93.fits

    The theoretical spectrum is scaled to match the observed spectrum from 1.5 -
    2.5 microns, and then it is used where the observed spectrum ends.
    The theoretical model of the Sun from Kurucz'93 atlas using the following
    parameters when the Sun is at 1 au.

     log_Z         T_eff        log_g           V_{Johnson}
     +0.0           5777        +4.44              -26.75

    Attributes
    ----------
    source: str
        filename of the sun library
    data: SimpleTable
        data table
    units: tuple
        detected units from file header
    wavelength: array
        wavelength (with units when found)
    flux: array
        flux(wavelength) values (with units when provided)
    distance: float
        distance to the observed Sun (default, 1 au)
    flavor: str, (default theoretical)
        either 'observed' using the stis reference,
        or  'theoretical' for the Kurucz model.
    """
    def __init__(self, source=None, distance=1 * Unit('au'), flavor='theoretical'):
        """ Constructor """
        self.data = None
        self.units = None
        self.distance = distance
        self.distance_conversion = 1.
        self._set_source_flavor(source, flavor)

    def _set_source_flavor(self, source, flavor):
        if not flavor.lower() in _default_sun:
            raise RuntimeError('Flavor must be either theoretical or observed')
        self.flavor = flavor
        if source is not None:
            self.source = source
        else:
            self.source = _default_sun[flavor]
            self.distance_conversion = float(((_default_distance / self.distance) ** 2))
            self.data = None
            self.units = None

    def _readfile(self):
        if self.data is not None:
            return
        fname = self.source
        # get extension
        self.data = SimpleTable(fname, silent=True)
        try:
            self.data.header['WAVELENGTH_UNIT'] = self.data._units['WAVELENGTH']
            self.data.header['FLUX_UNIT'] = self.data._units['FLUX']
            self.units = self.data._units['WAVELENGTH'], self.data._units['FLUX']
        except AttributeError:
            pass

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
        except Exception as e:
            print(e)
            return self.data.WAVELENGTH

    @property
    def flux(self):
        """ flux(wavelength) values (with units when provided) """
        self._readfile()
        try:
            return self.data.FLUX * Unit(self.units[1].lower()) * self.distance_conversion
        except Exception as e:
            print(e)
            return self.data.FLUX * self.distance_conversion
