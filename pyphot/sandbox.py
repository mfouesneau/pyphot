"""
Sandbox of new developments

Use at your own risks
"""

from __future__ import print_function, division
import numpy as np
from .phot import unit, _drop_units

from .vega import Vega
from .phot import Filter, set_method_default_units, Constants


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
        c = Constants.c.to('cm/s').magnitude      # c = 2.99792458e18 cm / s
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
