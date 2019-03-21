"""
Sandbox of new developments
"""

from __future__ import print_function, division
import numpy as np
from .phot import Filter
from .phot import unit, _drop_units
from .simpletable import SimpleTable


class UncertainFilter(Filter):
    """ What could be a filter with uncertainties

    Attributes
    ----------
    wavelength: ndarray
        wavelength sequence defining the filter transmission curve

    mean_: Filter
        mean passband transmission

    samples_: sequence(Filter)
        samples from the uncertain passband transmission model

    name: string
        name of the passband

    dtype: str
        detector type, either "photon" or "energy" counter

    unit: str
        wavelength units
    """
    def __init__(self, wavelength, mean_transmit, samples,
                 name='', dtype='photon', unit=None):
        """ Constructor """
        self.mean_ = Filter(wavelength, mean_transmit,
                            name=name, dtype=dtype, unit=unit)
        self.samples_ = [Filter(wavelength, transmit_k,
                                name=name + '_{0:d}'.format(num),
                                dtype=dtype, unit=unit)
                         for (num, transmit_k) in enumerate(samples)]
        self.name = name
        self.dtype = self.mean_.dtype
        self.model_ = None

    @classmethod
    def from_gp_model(cls, model, xprime=None, n_samples=10, **kwargs):
        """ Generate a filter object from a sklearn GP model 

        Parameters
        ----------
        model: sklearn.gaussian_process.GaussianProcessRegressor
            model of the passband
        xprime: ndarray
            wavelength to express the model in addition to the training points
        n_samples: int
            number of samples to generate from the model.
        **kwawrgs: dict
            UncertainFilter keywords
        """
        if xprime is None:
            xpred = model.X_train_
        else:
            xpred = np.unique(np.hstack([_drop_units(xprime), model.X_train_.ravel()]))
            xpred = xpred.reshape(1, -1).T

        unit_ = kwargs.pop('unit', None)
        if unit_ is None:
            unit_ = str(getattr(xprime, 'units', None))

        mean_transmit, _ = model.predict(xpred, return_std=True)
        samples = model.sample_y(xpred, n_samples=n_samples)

        unc_filter = cls(xpred.ravel(),
                         mean_transmit,
                         samples.T, unit=unit_, **kwargs)

        unc_filter.model_ = model
        return unc_filter

    def info(self, show_zeropoints=True):
        """ display information about the current filter"""
        string = self.mean_.info(show_zeropoints)
        string = string.replace('Filter object information',
                                'Filter object mean information only')
        return string

    def set_dtype(self, dtype):
        """ Set the detector type (photon or energy)"""
        self.mean_.set_dtype(dtype)
        for filter_k in self.samples_:
            filter_k.set_dtype(dtype)
        self.dtype = self.mean_.dtype

    def set_wavelength_unit(self, unit):
        """ Set the wavelength units """
        self.mean_.set_wavelength_unit(unit)
        for filter_k in self.samples_:
            filter_k.set_wavelength_unit(unit)

    @property
    def wavelength(self):
        """ Unitwise wavelength definition """
        return self.mean_.wavelength

    @property
    def wavelength_unit(self):
        """ Unit wavelength definition """
        return self.mean_.wavelength_unit

    @property
    def _wavelength(self):
        """ Unitless wavelength definition """
        return self.mean_._wavelength

    @property
    def transmit(self):
        """ Transmission curves """
        return self._get_mean_and_samples_attribute('transmit')

    def _get_samples_attribute(self, attr, *args, **kwargs):
        """ Returns the attribute from all samples """
        try:
            vals = [getattr(fk, attr)(*args, **kwargs) for fk in self.samples_]
        except TypeError:
            vals = [getattr(fk, attr) for fk in self.samples_]
        try:
            unit_ = unit[str(vals[0].units)]
            return np.array([v.magnitude for v in vals]) * unit_
        except AttributeError:
            return np.array(vals)

    def _get_mean_attribute(self, attr, *args, **kwargs):
        """ Returns the attribute from the mean passband """
        attr = getattr(self.mean_, attr)
        try:
            return attr(*args, **kwargs)
        except TypeError:
            return attr

    def _get_mean_and_samples_attribute(self, attr, *args, **kwargs):
        """ Compute / extract mean and smapled filter attributes 

        Parameters
        ----------
        attr: str
            attribute to get (can be a callable attribute)
        args: sequence
            any argument of attr
        kwargs: dict
            any keywords for attr

        Returns
        -------
        mean_: object
            value from the mean passband
        samples_: sequence(object)
            values from each sampled passband
        """
        return (self._get_mean_attribute(attr, *args, **kwargs),
                self._get_samples_attribute(attr, *args, **kwargs))

    @property
    def lmax(self):
        """ Calculated as the last value with a transmission at least 1% of
        maximum transmission """
        return self._get_mean_and_samples_attribute('lmax')

    @property
    def lmin(self):
        """ Calculate das the first value with a transmission at least 1% of
        maximum transmission """
        return self._get_mean_and_samples_attribute('lmin')

    @property
    def width(self):
        """ Effective width
        Equivalent to the horizontal size of a rectangle with height equal
        to maximum transmission and with the same area that the one covered by
        the filter transmission curve.

        W = int(T dlamb) / max(T)
        """
        return self._get_mean_and_samples_attribute('width')

    @property
    def fwhm(self):
        """ the difference between the two wavelengths for which filter
        transmission is half maximum

        ..note::
            This calculation is not exact but rounded to the nearest passband
            data points
        """
        return self._get_mean_and_samples_attribute('fwhm')

    @property
    def lpivot(self):
        """ Unitwise wavelength definition """
        return self._get_mean_and_samples_attribute('lpivot')


    @property
    def cl(self):
        """ Unitwise wavelength definition """
        return self._get_mean_and_samples_attribute('cl')

    @property
    def leff(self):
        """ Unitwise Effective wavelength
        leff = int (lamb * T * Vega dlamb) / int(T * Vega dlamb)
        """
        return self._get_mean_and_samples_attribute('leff')

    @property
    def lphot(self):
        """ Photon distribution based effective wavelength. Defined as

        lphot = int(lamb ** 2 * T * Vega dlamb) / int(lamb * T * Vega dlamb)

        which we calculate as

        lphot = get_flux(lamb * vega) / get_flux(vega)
        """
        return self._get_mean_and_samples_attribute('lphot')

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
        mean, samples = self._get_mean_and_samples_attribute('get_Nphotons',
                                                             slamb, sflux,
                                                             axis=axis)
        return mean, samples

    @property
    def Vega_zero_photons(self):
        """ Vega number of photons per wavelength unit

        .. note::

            see `self.get_Nphotons`
        """
        return self._get_mean_and_samples_attribute('Vega_zero_photons')

    def getFlux(self, slamb, sflux, axis=-1):
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
        mean, samples = self._get_mean_and_samples_attribute('getFlux',
                                                             slamb, sflux,
                                                             axis=axis)
        return mean, samples

    def reinterp(self, lamb):
        """ reinterpolate filter onto a different wavelength definition """
        mean, samples = self._get_mean_and_samples_attribute('reinterp')
        mean_val = mean(lamb)
        samp_val = [sk(mean_val.wavelength) for sk in samples]
        samp_transmissions = [sk.transmit for sk in samp_val]

        return self.__class__(mean_val.wavelength, mean_val.transmit,
                              samp_transmissions, name=self.name,
                              dtype=mean_val.dtype, unit=mean_val.wavelength_unit)

    def apply_transmission(self, slamb, sflux):
        """
        Apply filter transmission to a spectrum (with reinterpolation of the filter)

        Parameters
        ----------
        slamb: ndarray
            spectrum wavelength definition domain

        sflux: ndarray
            associated flux

        Returns
        -------
        flux: float
            new spectrum values accounting for the filter
        """
        mean, samples = self._get_mean_and_samples_attribute('apply_transmission')
        mean_val = mean(slamb, sflux)
        samp_val = [sk(slamb, sflux) for sk in samples]
        return mean_val, samp_val

    @property
    def AB_zero_mag(self):
        """ AB magnitude zero point
        ABmag = -2.5 * log10(f_nu) - 48.60
              = -2.5 * log10(f_lamb) - 2.5 * log10(lpivot ** 2 / c) - 48.60
              = -2.5 * log10(f_lamb) - zpts
        """
        return self._get_mean_and_samples_attribute('AB_zero_mag')

    @property
    def AB_zero_flux(self):
        """ AB flux zero point in erg/s/cm2/AA """
        return self._get_mean_and_samples_attribute('AB_zero_flux')

    @property
    def AB_zero_Jy(self):
        """ AB flux zero point in Jansky (Jy) """
        return self._get_mean_and_samples_attribute('AB_zero_Jy')

    @property
    def Vega_zero_mag(self):
        """ Vega magnitude zero point
        Vegamag = -2.5 * log10(f_lamb) + 2.5 * log10(f_vega)
        Vegamag = -2.5 * log10(f_lamb) - zpts
        """
        return self._get_mean_and_samples_attribute('Vega_zero_mag')

    @property
    def Vega_zero_flux(self):
        """ Vega flux zero point in erg/s/cm2/AA """
        return self._get_mean_and_samples_attribute('Vega_zero_flux')

    @property
    def Vega_zero_Jy(self):
        """ Vega flux zero point in Jansky (Jy) """
        return self._get_mean_and_samples_attribute('Vega_zero_Jy')

    @property
    def ST_zero_mag(self):
        """ ST magnitude zero point
        STmag = -2.5 * log10(f_lamb) -21.1
        """
        return 21.1

    @property
    def ST_zero_flux(self):
        """ ST flux zero point in erg/s/cm2/AA """
        return 10 ** (-0.4 * self.ST_zero_mag) * unit['erg/s/cm ** 2/AA']

    @property
    def ST_zero_Jy(self):
        """ ST flux zero point in Jansky (Jy) """
        return self._get_mean_and_samples_attribute('ST_zero_Jy')

    def to_Table(self, **kwargs):
        """ Export filter to a SimpleTable object

        Parameters
        ----------
        fname: str
            filename

        Uses `SimpleTable` parameters
        """
        mean_transmit, transmit_ = self.transmit
        data_ = {'WAVELENGTH': self._wavelength,
                 'THROUGHPUT': mean_transmit}
        for num, filterk in enumerate(transmit_, 1):
            data_['THROUGHPUT_{0:d}'.format(num)] = filterk
        data = SimpleTable(data_)

        if self.wavelength_unit is not None:
            data.header['WAVELENGTH_UNIT'] = self.wavelength_unit
        data.header['DETECTOR'] = self.dtype
        data.header['COMPNAME'] = self.name
        data.header['NAME'] = self.name
        data.set_comment('THROUGHPUT', 'filter throughput definition')
        data.set_comment('WAVELENGTH', 'filter wavelength definition')
        for num in range(1, len(transmit_) + 1):
            data.set_comment('THROUGHPUT_{0:d}'.format(num), 'filter throughput sample')
        data.set_comment('WAVELENGTH', self.wavelength_unit or 'AA')
        return data

    @classmethod
    def from_ascii(cls, fname, dtype='csv', **kwargs):
        """ Load filter from ascii file """
        lamb = kwargs.pop('lamb', None)
        name = kwargs.pop('name', None)
        detector = kwargs.pop('detector', 'photon')
        unit_ = kwargs.pop('unit', None)

        if not isinstance(fname, SimpleTable):
            t = SimpleTable(fname, dtype=dtype, **kwargs)
        else:
            t = fname
        w = t['WAVELENGTH'].astype(float)
        r = t['THROUGHPUT'].astype(float)
        keys = [k for k in t.keys() if 'THROUGHPUT_' in k]

        # update properties from file header
        detector = t.header.get('DETECTOR', detector)
        unit_    = t.header.get('WAVELENGTH_UNIT', unit_)

        # try from the comments in the header first
        if name in (None, 'None', 'none', ''):
            name = [k.split()[1]
                    for k in t.header.get('COMMENT', '').split('\n')
                    if 'COMPNAME' in k]
            name = ''.join(name).replace('"', '').replace("'", '')
        # if that did not work try the table header directly
        if name in (None, 'None', 'none', ''):
            name = t.header['NAME']

        if len(keys) > 0:
            samp = np.array([t[key] for key in keys])
            _filter = cls(w, r, samp, name=name, dtype=detector, unit=unit_)
        else:
            _filter = Filter(w, r, name=name, dtype=detector, unit=unit_)

        # reinterpolate if requested
        if lamb is not None:
            _filter = _filter.reinterp(lamb)

        return _filter
