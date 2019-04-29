"""
Sandbox of new developments

Use at your own risks
"""

from __future__ import print_function, division
import os
import numpy as np
from .phot import unit, _drop_units

from .vega import Vega
from .phot import Filter, set_method_default_units, Constants
from .phot import Library, Ascii_Library, HDF_Library, __default__

from .licks import reduce_resolution as _reduce_resolution
from .licks import LickIndex, LickLibrary


class UnitFilter(Filter):
    """ Evolution of Filter that makes sure the input spectra and output fluxes
    have units to avoid mis-interpretation.

    Note the usual (non SI) units of flux definitions:
        flam     = erg/s/cm**2/AA
        fnu      = erg/s/cm**2/Hz
        photflam = photon/s/cm**2/AA
        photnu   = photon/s/cm**2/Hz
    """
    def __init__(self, wavelength, transmit, name='', dtype="photon", unit=None):
        """Constructor"""
        try:   # get units from the inputs
            self._wavelength = wavelength.magnitude
            unit = str(wavelength.units)
        except AttributeError:
            self._wavelength = wavelength
        Filter.__init__(self, wavelength, transmit,
                        name=name, dtype=dtype, unit=unit)

    @classmethod
    def _validate_sflux(cls, slamb, sflux):
        """ clean data for inf in input """
        _sflux = _drop_units(sflux)
        _slamb = _drop_units(slamb)

        if True in np.isinf(sflux):
            indinf = np.where(np.isinf(_sflux))
            indfin = np.where(np.isfinite(_sflux))
            _sflux[indinf] = np.interp(_slamb[indinf], _slamb[indfin],
                                       _sflux[indfin], left=0, right=0)
        try:
            _unit = str(sflux.units)
            return _sflux * unit[_unit]
        except AttributeError:
            return _sflux

    @classmethod
    def _get_zero_like(cls, sflux, axis=-1):
        """ return a zero value corresponding to a flux calculation on sflux """
        _sflux = _drop_units(sflux)
        shape = _sflux.shape
        if axis < 0:
            axis = len(shape) + axis
        newshape = shape[:axis] + shape[axis + 1:]
        return np.zeros(newshape, _sflux.dtype)

    @property
    def lphot(self):
        """ Photon distribution based effective wavelength. Defined as

        lphot = int(lamb ** 2 * T * Vega dlamb) / int(lamb * T * Vega dlamb)

        which we calculate as

        lphot = get_flux(lamb * vega) / get_flux(vega)
        """
        if self.wavelength_unit is None:
            raise AttributeError('Needs wavelength units')

        with Vega() as v:
            wave = v.wavelength.magnitude
            # Cheating units to avoid making a new filter
            f_vega = self.get_flux(v.wavelength, v.flux, axis=-1)
            f_lamb_vega = self.get_flux(v.wavelength, wave * v.flux, axis=-1)
            f_lamb2_vega = self.get_flux(v.wavelength, wave ** 2 * v.flux, axis=-1)
        if 'photon' in self.dtype:
            lphot = (f_lamb_vega / f_vega)
        else:
            lphot = f_lamb2_vega / f_lamb_vega
        return lphot * unit[str(v.wavelength.units)]

    @set_method_default_units('AA', 'flam', output_unit='photon/s/cm**2/AA')
    def get_Nphotons(self, slamb, sflux, axis=-1):
        """getNphot the number of photons through the filter
        (Ntot / width  in the documentation)

        getflux() * leff / hc

        Parameters
        ----------
        slamb: ndarray(dtype=float, ndim=1)
            spectrum wavelength definition domain

        sflux: ndarray(dtype=float, ndim=1)
            associated flux in erg/s/cm2/AA

        Returns
        -------
        N: float
            Number of photons of the spectrum within the filter
        """
        passb = self.reinterp(slamb)
        wave = passb._wavelength
        dlambda = np.diff(wave)

        h = Constants.h.to('erg*s').magnitude     # h = 6.626075540e-27 erg * s
        c = Constants.c.to('AA/s').magnitude      # c = 2.99792458e18 cm / s
        vals = sflux.magnitude * wave * passb.transmit
        vals[~np.isfinite(vals)] = 0.
        Nphot = 0.5 * np.sum((vals[1:] + vals[:-1]) * dlambda) / (h * c)
        Nphot = Nphot * unit['photon/s/cm**2']
        return Nphot / passb.width   # photons / cm2 / s / A

    @set_method_default_units('AA', 'flam', output_unit='erg/s/cm**2/AA')
    def get_flux(self, slamb, sflux, axis=-1):
        """getFlux
        Integrate the flux within the filter and return the integrated energy
        If you consider applying the filter to many spectra, you might want to
        consider extractSEDs.

        Parameters
        ----------
        slamb: ndarray(dtype=float, ndim=1)
            spectrum wavelength definition domain

        sflux: ndarray(dtype=float, ndim=1)
            associated flux

        Returns
        -------
        flux: float
            Energy of the spectrum within the filter
        """
        passb = self.reinterp(slamb)
        ifT = passb.transmit
        _slamb = _drop_units(slamb)
        _sflux = _drop_units(passb._validate_sflux(slamb, sflux))

        _w_unit = str(slamb.units)
        _f_unit = str(sflux.units)

        # if the filter is null on that wavelength range flux is then 0
        # ind = ifT > 0.
        nonzero = np.where(ifT > 0)[0]
        if nonzero.size <= 0:
            return passb._get_zero_like(sflux)

        # avoid calculating many zeros
        nonzero_start = max(0, min(nonzero) - 5)
        nonzero_end = min(len(ifT), max(nonzero) + 5)
        ind = np.zeros(len(ifT), dtype=bool)
        ind[nonzero_start:nonzero_end] = True

        if True in ind:
            try:
                _sflux = _sflux[:, ind]
            except:
                _sflux = _sflux[ind]
            # limit integrals to where necessary
            if 'photon' in passb.dtype:
                a = np.trapz(_slamb[ind] * ifT[ind] * _sflux, _slamb[ind], axis=axis)
                b = np.trapz(_slamb[ind] * ifT[ind], _slamb[ind] )
                a = a * unit['*'.join((_w_unit, _f_unit, _w_unit))]
                b = b * unit['*'.join((_w_unit, _w_unit))]
            elif 'energy' in passb.dtype:
                a = np.trapz( ifT[ind] * _sflux, _slamb[ind], axis=axis )
                b = np.trapz( ifT[ind], _slamb[ind])
                a = a * unit['*'.join((_f_unit, _w_unit))]
                b = b * unit[_w_unit]
            if (np.isinf(a.magnitude).any() | np.isinf(b.magnitude).any()):
                print(self.name, "Warn for inf value")
            return a / b
        else:
            return passb._get_zero_like(sflux)

    def getFlux(self, slamb, sflux, axis=-1):
        """
        Integrate the flux within the filter and return the integrated energy
        If you consider applying the filter to many spectra, you might want to
        consider extractSEDs.

        Parameters
        ----------
        slamb: ndarray(dtype=float, ndim=1)
            spectrum wavelength definition domain

        sflux: ndarray(dtype=float, ndim=1)
            associated flux

        Returns
        -------
        flux: float
            Energy of the spectrum within the filter
        """
        return self.get_flux(slamb, sflux, axis=axis)

    @property
    def Vega_zero_flux(self):
        """ Vega flux zero point in erg/s/cm2/AA """
        with Vega() as v:
            f_vega = self.get_flux(v.wavelength, v.flux, axis=-1)
        return f_vega

    @property
    def Vega_zero_mag(self):
        """ vega magnitude zero point
        vegamag = -2.5 * log10(f_lamb) + 2.5 * log10(f_vega)
        vegamag = -2.5 * log10(f_lamb) - zpts
        """
        return -2.5 * np.log10(self.Vega_zero_flux.magnitude)


class UnitAscii_Library(Ascii_Library):

    def _load_filter(self, fname, interp=True, lamb=None, *args, **kwargs):
        """ Load a given filter from the library """
        try:
            fil = UnitFilter.from_ascii(fname, *args, **kwargs)
        except:
            content = self.content
            r = [k for k in content if fname in k]

            if len(r) <= 0:  # try all lower for filenames (ascii convention)
                r = [k for k in content if fname.lower() in k]

            if len(r) > 1:
                print("auto correction found multiple choices")
                print(r)
                raise ValueError('Refine name to one of {0}'.format(r))
            elif len(r) <= 0:
                raise ValueError('Cannot find filter {0}'.format(fname))
            else:
                fil = UnitFilter.from_ascii(r[0], *args, **kwargs)
        if (interp is True) and (lamb is not None):
            return fil.reinterp(lamb)
        else:
            return fil


class UnitHDF_Library(HDF_Library):

    def _load_filter(self, fname, interp=True, lamb=None):
        """ Load a given filter from the library

        Parameters
        ----------
        fname: str
            normalized names according to filtersLib

        interp: bool, optional
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        integrationFilter: bool, optional
            set True for specail integraion filter such as Qion or E_uv
            if set, lamb should be given

        Returns
        -------
        filter: Filter instance
            filter object
        """
        ftab = self.hdf
        if hasattr(fname, 'decode'):
            fnode    = ftab.get_node('/filters/' + fname.decode('utf8'))
        else:
            fnode    = ftab.get_node('/filters/' + fname)
        flamb    = fnode[:]['WAVELENGTH']
        transmit = fnode[:]['THROUGHPUT']
        dtype = 'photon'
        unit = None

        attrs = fnode.attrs
        if 'DETECTOR' in attrs:
            dtype = attrs['DETECTOR']
        if 'WAVELENGTH_UNIT' in attrs:
            unit = attrs['WAVELENGTH_UNIT']

        fil = UnitFilter(flamb, transmit, name=fnode.name, dtype=dtype, unit=unit)

        if interp & (lamb is not None):
            fil = fil.reinterp(lamb)
        return fil


class UnitLibray(Library):

    @classmethod
    def from_ascii(cls, filename, **kwargs):
        return UnitAscii_Library(filename, **kwargs)

    @classmethod
    def from_hd5(cls, filename, **kwargs):
        return UnitHDF_Library(filename, **kwargs)

def get_library(fname=__default__, **kwargs):
    """ Finds the appropriate class to load the library """
    if os.path.isfile(fname):
        return UnitHDF_Library(fname, **kwargs)
    else:
        return UnitAscii_Library(fname, **kwargs)


@set_method_default_units('AA', 'flam', output_unit='flam')
def reduce_resolution(wi, fi, fwhm0=0.55 * unit['AA'], sigma_floor=0.2 * unit['AA']):
    """ Adapt the resolution of the spectra to match the lick definitions

        Lick definitions have different resolution elements as function
        of wavelength. These definition are hard-coded in this function

    Parameters
    ---------
    wi: ndarray (n, )
        wavelength definition
    fi: ndarray (nspec, n) or (n, )
        spectra to convert
    fwhm0: float
        initial broadening in the spectra `fi`
    sigma_floor: float
        minimal dispersion to consider

    Returns
    -------
    flux_red: ndarray (nspec, n) or (n, )
        reduced spectra
    """
    flux_red = _reduce_resolution(wi.magnitude, fi.magnitude,
                                  fwhm0.to('AA').magnitude,
                                  sigma_floor.to('AA').magnitude)
    return flux_red * unit['flam']


class UnitLickIndex(LickIndex):
    """ Define a Lick Index similarily to a Filter object """

    @set_method_default_units('AA', 'flam')
    def get(self, wave, flux, **kwargs):
        """ compute spectral index after continuum subtraction

        Parameters
        ----------
        w: ndarray (nw, )
            array of wavelengths in AA
        flux: ndarray (N, nw)
            array of flux values for different spectra in the series
        degree: int (default 1)
            degree of the polynomial fit to the continuum
        nocheck: bool
            set to silently pass on spectral domain mismatch.
            otherwise raises an error when index is not covered

        Returns
        -------
        ew: ndarray (N,)
            equivalent width or magnitude array

        Raises
        ------
        ValueError: when the spectral coverage wave does not cover the index
        range
        """
        LickIndex.get(wave, flux.to('flam').magnitude, **kwargs)


class UnitLickLibrary(LickLibrary):
    """ Collection of Lick indices """

    def _load_filter(self, fname, **kwargs):
        """ Load a given filter from the library """
        with self as current_lib:
            return UnitLickIndex(fname, current_lib._content[fname])
