"""
Sandbox of new developments

Use at your own risks

Photometric package using Astropy Units
=======================================

Defines a Filter class and associated functions to extract photometry.

This also include functions to keep libraries up to date

.. note::

    integrations are done using :func:`trapezoid`
    Why not Simpsons? Simpsons principle is to take sequence of 3 points to
    make a quadratic interpolation. Which in the end, when filters have sharp
    edges, the error due to this "interpolation" are extremely large in
    comparison to the uncertainties induced by trapeze integration.
"""

from __future__ import division, print_function

import os
from functools import wraps

import numpy as np
import tables

try:
    from scipy.integrate import trapezoid
except ImportError:  # older scipy / numpy < 2.0
    from scipy.integrate import trapz as trapezoid  # type: ignore

from ..simpletable import SimpleTable
from .config import __vega_default_flavor__, libsdir
from .vega import Vega

# from .licks import LickIndex, LickLibrary

# directories
__default__ = libsdir.joinpath("new_filters.hd5")
__default_lick__ = libsdir.joinpath("licks.dat")

from astropy import constants
from astropy.units import Unit


class Constants(object):
    """A namespace for constants"""

    # Planck's constant in erg * sec
    h = constants.h.to("erg * s")
    # Speed of light in cm/s
    c = constants.c.to("AA/s")


def hasUnit(val):
    """Check is an object has units"""
    return hasattr(val, "unit") or hasattr(val, "units")


class set_method_default_units(object):
    """Decorator for classmethods that makes sure that
    the inputs of slamb, sflux are in given units

    expects the decorated method to be defined as
        >> def methodname(self, lamb, flux)
    """

    def __init__(self, wavelength_unit, flux_unit, output_unit=None):
        self.wavelength_unit = Unit(wavelength_unit)
        self.flux_unit = Unit(flux_unit)
        self.output_unit = output_unit

    @classmethod
    def force_units(cls, value, unit):
        if unit is None:
            return value
        try:
            return value.to(unit)
        except AttributeError:
            msg = "Warning: assuming {0:s} units to unitless object."
            print(msg.format(str(unit)))
            return value * unit

    def __call__(self, func):
        @wraps(func)
        def wrapper(filter_, slamb, sflux, *args, **kwargs):
            _slamb = set_method_default_units.force_units(slamb, self.wavelength_unit)
            _sflux = set_method_default_units.force_units(sflux, self.flux_unit)
            output = func(filter_, _slamb, _sflux, *args, **kwargs)
            return set_method_default_units.force_units(output, self.output_unit)

        return wrapper


def _drop_units(q):
    """Drop the unit definition silently"""
    try:
        return q.value
    except AttributeError:
        try:
            return q.value
        except AttributeError:
            return q


class UnitFilter(object):
    """Evolution of Filter that makes sure the input spectra and output fluxes
    have units to avoid mis-interpretation.

    Note the usual (non SI) units of flux definitions:
        flam     = erg/s/cm**2/AA
        fnu      = erg/s/cm**2/Hz
        photflam = photon/s/cm**2/AA
        photnu   = photon/s/cm**2/Hz

    Define a filter by its name, wavelength and transmission
    The type of detector (energy or photon counter) can be specified for
    adapting calculations. (default: photon)

    Attributes
    ----------
    name: str
        name of the filter

    cl: float
        central wavelength of the filter

    norm: float
        normalization factor of the filter

    lpivot: float
        pivot wavelength of the filter

    wavelength: ndarray
        wavelength sequence defining the filter transmission curve

    transmit: ndarray
        transmission curve of the filter

    dtype: str
        detector type, either "photon" or "energy" counter

    unit: str
        wavelength units
    """

    def __init__(
        self,
        wavelength,
        transmit,
        name="",
        dtype="photon",
        unit=None,
        vega=__vega_default_flavor__,
    ):
        """Constructor"""
        self.name = name
        self.set_dtype(dtype)
        try:  # get units from the inputs
            self._wavelength = wavelength.value
            unit = str(wavelength.unit)
        except AttributeError:
            self._wavelength = wavelength
        self.set_wavelength_unit(unit)
        self._vega_flavor = vega

        # make sure input data are ordered and cleaned of weird values.
        idx = np.argsort(self._wavelength)
        self._wavelength = self._wavelength[idx]
        self.transmit = np.clip(transmit[idx], 0.0, np.nanmax(transmit))

        self.norm = trapezoid(self.transmit, self._wavelength)
        self._lT = trapezoid(self._wavelength * self.transmit, self._wavelength)
        self._lpivot = self._calculate_lpivot()
        if self.norm > 0:
            self._cl = self._lT / self.norm
        else:
            self._cl = 0.0

    def _calculate_lpivot(self):
        if self.transmit.max() <= 0:
            return 0.0
        if "photon" in self.dtype:
            lpivot2 = self._lT / trapezoid(
                self.transmit / self._wavelength, self._wavelength
            )
        else:
            lpivot2 = self.norm / trapezoid(
                self.transmit / self._wavelength**2, self._wavelength
            )
        return np.sqrt(lpivot2)

    def set_wavelength_unit(self, unit):
        """Set the wavelength units"""
        try:  # get units from the inputs
            self.wavelength_unit = str(self._wavelength.unit)
        except AttributeError:
            self.wavelength_unit = unit

    def set_dtype(self, dtype):
        """Set the detector type (photon or energy)"""
        _d = dtype.lower()
        if "phot" in _d:
            self.dtype = "photon"
        elif "ener" in _d:
            self.dtype = "energy"
        else:
            raise ValueError("Unknown detector type {0}".format(dtype))

    def info(self, show_zeropoints=True):
        """display information about the current filter"""
        msg = """Filter object information:
    name:                 {s.name:s}
    detector type:        {s.dtype:s}
    wavelength units:     {s.wavelength_unit}
    central wavelength:   {s.cl:f}
    pivot wavelength:     {s.lpivot:f}
    effective wavelength: {s.leff:f}
    photon wavelength:    {s.lphot:f}
    minimum wavelength:   {s.lmin:f}
    maximum wavelength:   {s.lmax:f}
    norm:                 {s.norm:f}
    effective width:      {s.width:f}
    fullwidth half-max:   {s.fwhm:f}
    definition contains {s.transmit.size:d} points"""
        print(msg.format(s=self).replace("None", "unknown"))

        # zero points only if units
        if (self.wavelength_unit is None) or (not show_zeropoints):
            return

        print(
            """
    Zeropoints
        Vega: {s.Vega_zero_mag:f} mag,
              {s.Vega_zero_flux},
              {s.Vega_zero_Jy}
              {s.Vega_zero_photons}
          AB: {s.AB_zero_mag:f} mag,
              {s.AB_zero_flux},
              {s.AB_zero_Jy}
          ST: {s.ST_zero_mag:f} mag,
              {s.ST_zero_flux},
              {s.ST_zero_Jy}
        """.format(s=self)
        )

    def __repr__(self):
        return "Filter: {0:s}, {1:s}".format(self.name, object.__repr__(self))

    @property
    def wavelength(self):
        """Unitwise wavelength definition"""
        if self.wavelength_unit is not None:
            return self._wavelength * Unit(self.wavelength_unit)
        else:
            return self._wavelength

    @property
    def lmax(self):
        """Calculated as the last value with a transmission at least 1% of
        maximum transmission"""
        cond = (self.transmit / self.transmit.max()) > 1.0 / 100
        return max(self.wavelength[cond])

    @property
    def lmin(self):
        """Calculate das the first value with a transmission at least 1% of
        maximum transmission"""
        cond = (self.transmit / self.transmit.max()) > 1.0 / 100
        return min(self.wavelength[cond])

    @property
    def width(self):
        """Effective width
        Equivalent to the horizontal size of a rectangle with height equal
        to maximum transmission and with the same area that the one covered by
        the filter transmission curve.

        W = int(T dlamb) / max(T)
        """
        return (self.norm / max(self.transmit)) * Unit(self.wavelength_unit)

    @property
    def fwhm(self):
        """the difference between the two wavelengths for which filter
        transmission is half maximum

        ..note::
            This calculation is not exact but rounded to the nearest passband
            data points
        """
        vals = self.transmit / self.transmit.max() - 0.5
        zero_crossings = np.where(np.diff(np.sign(vals)))[0]
        lambs = self.wavelength[zero_crossings]
        return np.diff(lambs)[0]

    @property
    def lpivot(self):
        """Unitwise wavelength definition"""
        if self.wavelength_unit is not None:
            return self._lpivot * Unit(self.wavelength_unit)
        else:
            return self._lpivot

    @property
    def cl(self):
        """Unitwise wavelength definition"""
        if self.wavelength_unit is not None:
            return self._cl * Unit(self.wavelength_unit)
        else:
            return self._cl

    @property
    def leff(self):
        """Unitwise Effective wavelength
        leff = int (lamb * T * Vega dlamb) / int(T * Vega dlamb)
        """
        with Vega(flavor=self._vega_flavor) as v:
            s = self.reinterp(v.wavelength)
            w = s._wavelength
            if s.transmit.max() > 0:
                leff = trapezoid(w * s.transmit * v.flux.value, w, axis=-1)
                leff /= trapezoid(s.transmit * v.flux.value, w, axis=-1)
            else:
                leff = float("nan")
        if s.wavelength_unit is not None:
            leff = leff * Unit(s.wavelength_unit)
            if self.wavelength_unit is not None:
                return leff.to(self.wavelength_unit)
            return leff
        else:
            return leff

    @classmethod
    def _validate_sflux(cls, slamb, sflux):
        """clean data for inf in input"""
        _sflux = _drop_units(sflux)
        _slamb = _drop_units(slamb)

        if np.isinf(_sflux).any():
            indinf = np.where(np.isinf(_sflux))
            indfin = np.where(np.isfinite(_sflux))
            _sflux[indinf] = np.interp(
                _slamb[indinf], _slamb[indfin], _sflux[indfin], left=0, right=0
            )
        try:
            _unit = str(sflux.unit)
            return _sflux * Unit(_unit)
        except AttributeError:
            return _sflux

    @classmethod
    def _get_zero_like(cls, sflux, axis=-1):
        """return a zero value corresponding to a flux calculation on sflux"""
        # _sflux = _drop_units(sflux)
        # shape = _sflux.shape
        # if axis < 0:
        #     axis = len(shape) + axis
        # newshape = shape[:axis] + shape[axis + 1:]
        # return np.zeros(newshape, _sflux.dtype)
        return np.zeros_like(sflux).sum(axis=axis)

    @property
    def lphot(self):
        """Photon distribution based effective wavelength. Defined as

        lphot = int(lamb ** 2 * T * Vega dlamb) / int(lamb * T * Vega dlamb)

        which we calculate as

        lphot = get_flux(lamb * vega) / get_flux(vega)
        """
        if self.wavelength_unit is None:
            raise AttributeError("Needs wavelength units")

        with Vega(flavor=self._vega_flavor) as v:
            wave = v.wavelength.value
            # Cheating units to avoid making a new filter
            f_vega = self.get_flux(v.wavelength, v.flux, axis=-1)
            f_lamb_vega = self.get_flux(v.wavelength, wave * v.flux, axis=-1)
            f_lamb2_vega = self.get_flux(v.wavelength, wave**2 * v.flux, axis=-1)
        if "photon" in self.dtype:
            lphot = f_lamb_vega / f_vega
        else:
            lphot = f_lamb2_vega / f_lamb_vega
        return (lphot * Unit(str(v.wavelength.unit))).to(self.wavelength_unit)

    def _get_filter_in_units_of(self, slamb=None):
        w = self.wavelength
        if hasUnit(slamb) & hasUnit(w):
            return w.to(str(slamb.unit)).value
        else:
            print("Warning: assuming units are consistent")
            return self._wavelength

    @set_method_default_units("AA", "flam", output_unit="photon*s**-1*cm**-2*AA**-1")
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

        # h = 6.626075540e-27    # erg * s
        # c = 2.99792458e18      # AA / s
        h = Constants.h.to("erg * s").value
        c = Constants.c.to("AA/s").value
        vals = sflux.value * wave * passb.transmit
        vals[~np.isfinite(vals)] = 0.0
        Nphot = 0.5 * np.sum((vals[1:] + vals[:-1]) * dlambda) / (h * c)
        Nphot = Nphot * Unit("photon*s**-1*cm**-2")
        return Nphot / passb.width  # photons / cm2 / s / A

    @property
    def Vega_zero_photons(self):
        """Vega number of photons per wavelength unit

        .. note::

            see `self.get_Nphotons`
        """
        with Vega(flavor=self._vega_flavor) as v:
            return self.get_Nphotons(v.wavelength, v.flux)

    @set_method_default_units("AA", "flam", output_unit="erg*s**-1*cm**-2*AA**-1")
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

        _w_unit = str(slamb.unit)
        _f_unit = str(sflux.unit)

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
            except Exception:
                _sflux = _sflux[ind]
            # limit integrals to where necessary
            if "photon" in passb.dtype:
                a = trapezoid(_slamb[ind] * ifT[ind] * _sflux, _slamb[ind], axis=axis)
                b = trapezoid(_slamb[ind] * ifT[ind], _slamb[ind])
                a = a * Unit("*".join((_w_unit, _f_unit, _w_unit)))
                b = b * Unit("*".join((_w_unit, _w_unit)))
            elif "energy" in passb.dtype:
                a = trapezoid(ifT[ind] * _sflux, _slamb[ind], axis=axis)
                b = trapezoid(ifT[ind], _slamb[ind])
                a = a * Unit("*".join((_f_unit, _w_unit)))
                b = b * Unit(_w_unit)
            if np.isinf(a.value).any() | np.isinf(b.value).any():
                print(self.name, "Warn for inf value")
            return a / b
        else:
            return passb._get_zero_like(_sflux)

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

    def reinterp(self, lamb):
        """reinterpolate filter onto a different wavelength definition"""
        _wavelength = self._get_filter_in_units_of(lamb)
        _lamb = _drop_units(lamb)
        try:
            _unit = str(lamb.unit)
        except Exception:
            _unit = self.wavelength_unit
        ifT = np.interp(_lamb, _wavelength, self.transmit, left=0.0, right=0.0)
        return self.__class__(_lamb, ifT, name=self.name, dtype=self.dtype, unit=_unit)

    def __call__(self, slamb, sflux):
        return self.applyTo(slamb, sflux)

    def apply_transmission(self, slamb, sflux):
        """
        Apply filter transmission to a spectrum (with reinterpolation of the
        filter)

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
        _wavelength = self._get_filter_in_units_of(slamb)
        _lamb = _drop_units(slamb)
        ifT = np.interp(_lamb, _wavelength, self.transmit, left=0.0, right=0.0)
        return ifT * sflux

    def applyTo(self, slamb, sflux):
        """For compatibility but bad name"""
        return self.apply_transmission(slamb, sflux)

    @classmethod
    def from_ascii(cls, fname, dtype="csv", **kwargs):
        """Load filter from ascii file"""
        lamb = kwargs.pop("lamb", None)
        name = kwargs.pop("name", None)
        detector = kwargs.pop("detector", "photon")
        unit = kwargs.pop("unit", None)

        t = SimpleTable(fname, dtype=dtype, **kwargs)
        w = t["WAVELENGTH"].astype(float)
        r = t["THROUGHPUT"].astype(float)

        # update properties from file header
        detector = t.header.get("DETECTOR", detector)
        unit = t.header.get("WAVELENGTH_UNIT", unit)
        name = t.header.get("NAME", name)

        # try from the comments in the header first
        if name in (None, "None", "none", ""):
            name = [
                k.split()[1]
                for k in t.header.get("COMMENT", "").split("\n")
                if "COMPNAME" in k
            ]
            name = "".join(name).replace('"', "").replace("'", "")
        # if that did not work try the table header directly
        if name in (None, "None", "none", ""):
            name = t.header["NAME"]

        _filter = UnitFilter(w, r, name=name, dtype=detector, unit=unit)

        # reinterpolate if requested
        if lamb is not None:
            _filter = _filter.reinterp(lamb)

        return _filter

    def write_to(self, fname, **kwargs):
        """Export filter to a file

        Parameters
        ----------
        fname: str
            filename

        Uses `SimpleTable.write` parameters
        """
        data = self.to_Table()
        data.write(fname, **kwargs)

    def to_Table(self, **kwargs):
        """Export filter to a SimpleTable object

        Parameters
        ----------
        fname: str
            filename

        Uses `SimpleTable` parameters
        """
        data = SimpleTable(
            {"WAVELENGTH": self._wavelength, "THROUGHPUT": self.transmit}
        )

        if self.wavelength_unit is not None:
            data.header["WAVELENGTH_UNIT"] = self.wavelength_unit
        data.header["DETECTOR"] = self.dtype
        data.header["COMPNAME"] = str(self.name)
        data.header["NAME"] = str(self.name)
        data.set_comment("THROUGHPUT", "filter throughput definition")
        data.set_comment("WAVELENGTH", "filter wavelength definition")
        data.set_comment("WAVELENGTH", self.wavelength_unit or "AA")
        return data

    def to_dict(self):
        """Return a dictionary of the filter"""
        data = {"WAVELENGTH": self._wavelength, "THROUGHPUT": self.transmit}
        if self.wavelength_unit is not None:
            data["WAVELENGTH_UNIT"] = self.wavelength_unit
        data["DETECTOR"] = self.dtype
        data["NAME"] = self.name
        data["PIVOT"] = self._lpivot
        data["CENTRAL"] = self._cl
        data["EFFECTIVE"] = _drop_units(self.leff)
        data["NORM"] = self.norm
        return data

    @classmethod
    def make_integration_filter(cls, lmin, lmax, name="", dtype="photon", unit=None):
        """Generate an heavyside filter between lmin and lmax"""
        dyn = lmax - lmin
        try:
            unit = str(dyn.unit)
            dyn = _drop_units(dyn)
        except Exception:
            pass
        w = np.array([lmin - 0.01 * dyn, lmin, lmax, lmax + 0.01 * dyn])
        f = np.array([0.0, 1.0, 1.0, 0.0])
        return UnitFilter(w, f, name=name, dtype=dtype, unit=unit)

    @property
    def AB_zero_mag(self):
        """AB magnitude zero point

        ABmag = -2.5 * log10(f_nu) - 48.60
              = -2.5 * log10(f_lamb) - 2.5 * log10(lpivot ** 2 / c) - 48.60
              = -2.5 * log10(f_lamb) - zpts

        """
        if self.wavelength_unit is None:
            raise AttributeError("Needs wavelength units")

        C1 = Unit(self.wavelength_unit).to("AA") ** 2 / Constants.c.to("AA/s").value
        c1 = self._lpivot**2 * C1

        m = 2.5 * np.log10(_drop_units(c1)) + 48.6
        return m

    @property
    def AB_zero_flux(self):
        """AB flux zero point in erg/s/cm2/AA"""
        return 10 ** (-0.4 * self.AB_zero_mag) * Unit("erg*s**-1*cm**-2*AA**-1")

    @property
    def AB_zero_Jy(self):
        """AB flux zero point in Jansky (Jy)"""
        c = 1e-8 * Constants.c.to("m/s").value
        f = 1e5 / c * self.lpivot.to("AA").value ** 2 * self.AB_zero_flux.value
        return f * Unit("Jy")

    @property
    def Vega_zero_mag(self):
        """vega magnitude zero point
        vegamag = -2.5 * log10(f_lamb) + 2.5 * log10(f_vega)
        vegamag = -2.5 * log10(f_lamb) - zpts
        """
        flux = self.Vega_zero_flux.value
        if flux > 0:
            return -2.5 * np.log10(flux)
        else:
            return float("nan")

    @property
    def Vega_zero_flux(self):
        """Vega flux zero point in erg/s/cm2/AA"""
        with Vega(flavor=self._vega_flavor) as v:
            f_vega = self.get_flux(v.wavelength, v.flux, axis=-1)
        return f_vega

    @property
    def Vega_zero_Jy(self):
        """Vega flux zero point in Jansky (Jy)"""
        c = 1e-8 * Constants.c.to("m/s").value
        f = (
            1e5
            / c
            * (
                self.lpivot.to("AA").value ** 2
                * self.Vega_zero_flux.to("erg*s**-1*cm**-2*AA**-1").value
            )
        )
        return f * Unit("Jy")

    @property
    def ST_zero_mag(self):
        """ST magnitude zero point
        STmag = -2.5 * log10(f_lamb) -21.1
        """
        return 21.1

    @property
    def ST_zero_flux(self):
        """ST flux zero point in erg/s/cm2/AA"""
        return 10 ** (-0.4 * self.ST_zero_mag) * Unit("erg*s**-1*cm**-2*AA**-1")

    @property
    def ST_zero_Jy(self):
        """ST flux zero point in Jansky (Jy)"""
        c = 1e-8 * Constants.c.to("m/s").value
        f = 1e5 / c * self.lpivot.to("AA").value ** 2 * self.ST_zero_flux.value
        return f * Unit("Jy")


class UncertainFilter(UnitFilter):
    """What could be a filter with uncertainties

    Attributes
    ----------
    wavelength: ndarray
        wavelength sequence defining the filter transmission curve

    mean: Filter
        mean passband transmission

    samples: sequence(Filter)
        samples from the uncertain passband transmission model

    name: string
        name of the passband

    dtype: str
        detector type, either "photon" or "energy" counter

    unit: str
        wavelength units
    """

    def __init__(
        self, wavelength, mean_transmit, samples, name="", dtype="photon", unit=None
    ):
        """Constructor"""
        self.mean_ = UnitFilter(
            wavelength, mean_transmit, name=name, dtype=dtype, unit=unit
        )
        self.samples_ = [
            UnitFilter(
                wavelength,
                transmit_k,
                name=name + "_{0:d}".format(num),
                dtype=dtype,
                unit=unit,
            )
            for (num, transmit_k) in enumerate(samples)
        ]
        self.name = name
        self.dtype = self.mean_.dtype
        self.model_ = None

    @classmethod
    def from_gp_model(cls, model, xprime=None, n_samples=10, **kwargs):
        """Generate a filter object from a sklearn GP model

        Parameters
        ----------
        model: sklearn.gaussian_process.GaussianProcessRegressor
            model of the passband
        xprime: ndarray
            wavelength to express the model in addition to the training points
        n_samples: int
            number of samples to generate from the model.
        kwargs: dict
            UncertainFilter keywords
        """
        if xprime is None:
            xpred = model.X_train_
        else:
            xpred = np.unique(np.hstack([_drop_units(xprime), model.X_train_.ravel()]))
            xpred = xpred.reshape(1, -1).T

        unit_ = kwargs.pop("unit", None)
        if unit_ is None:
            unit_ = str(getattr(xprime, "units", None))

        mean_transmit, _ = model.predict(xpred, return_std=True)
        samples = model.sample_y(xpred, n_samples=n_samples)

        unc_filter = cls(xpred.ravel(), mean_transmit, samples.T, unit=unit_, **kwargs)

        unc_filter.model_ = model
        return unc_filter

    def info(self, show_zeropoints=True):
        """display information about the current filter"""
        string = self.mean_.info(show_zeropoints)
        string = string.replace(
            "Filter object information", "Filter object mean information only"
        )
        return string

    def set_dtype(self, dtype):
        """Set the detector type (photon or energy)"""
        self.mean_.set_dtype(dtype)
        for filter_k in self.samples_:
            filter_k.set_dtype(dtype)
        self.dtype = self.mean_.dtype

    def set_wavelength_unit(self, unit):
        """Set the wavelength units"""
        self.mean_.set_wavelength_unit(unit)
        for filter_k in self.samples_:
            filter_k.set_wavelength_unit(unit)

    @property
    def wavelength(self):
        """Unitwise wavelength definition"""
        return self.mean_.wavelength

    @property
    def wavelength_unit(self):
        """Unit wavelength definition"""
        return self.mean_.wavelength_unit

    @property
    def _wavelength(self):
        """Unitless wavelength definition"""
        return self.mean_._wavelength

    @property
    def transmit(self):
        """Transmission curves"""
        return self._get_mean_and_samples_attribute("transmit")

    def _get_samples_attribute(self, attr, *args, **kwargs):
        """Returns the attribute from all samples"""
        try:
            vals = [getattr(fk, attr)(*args, **kwargs) for fk in self.samples_]
        except TypeError:
            vals = [getattr(fk, attr) for fk in self.samples_]
        try:
            unit_ = Unit(str(vals[0].unit))
            return np.array([v.value for v in vals]) * unit_
        except AttributeError:
            return np.array(vals)

    def _get_mean_attribute(self, attr, *args, **kwargs):
        """Returns the attribute from the mean passband"""
        attr = getattr(self.mean_, attr)
        try:
            return attr(*args, **kwargs)
        except TypeError:
            return attr

    def _get_mean_and_samples_attribute(self, attr, *args, **kwargs):
        """Compute / extract mean and smapled filter attributes

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
        return (
            self._get_mean_attribute(attr, *args, **kwargs),
            self._get_samples_attribute(attr, *args, **kwargs),
        )

    @property
    def lmax(self):
        """Calculated as the last value with a transmission at least 1% of
        maximum transmission"""
        return self._get_mean_and_samples_attribute("lmax")

    @property
    def lmin(self):
        """Calculate das the first value with a transmission at least 1% of
        maximum transmission"""
        return self._get_mean_and_samples_attribute("lmin")

    @property
    def width(self):
        """Effective width
        Equivalent to the horizontal size of a rectangle with height equal
        to maximum transmission and with the same area that the one covered by
        the filter transmission curve.

        W = int(T dlamb) / max(T)
        """
        return self._get_mean_and_samples_attribute("width")

    @property
    def fwhm(self):
        """the difference between the two wavelengths for which filter
        transmission is half maximum

        ..note::
            This calculation is not exact but rounded to the nearest passband
            data points
        """
        return self._get_mean_and_samples_attribute("fwhm")

    @property
    def lpivot(self):
        """Unitwise wavelength definition"""
        return self._get_mean_and_samples_attribute("lpivot")

    @property
    def cl(self):
        """Unitwise wavelength definition"""
        return self._get_mean_and_samples_attribute("cl")

    @property
    def leff(self):
        """Unitwise Effective wavelength
        leff = int (lamb * T * Vega dlamb) / int(T * Vega dlamb)
        """
        return self._get_mean_and_samples_attribute("leff")

    @property
    def lphot(self):
        """Photon distribution based effective wavelength. Defined as

        lphot = int(lamb ** 2 * T * Vega dlamb) / int(lamb * T * Vega dlamb)

        which we calculate as

        lphot = get_flux(lamb * vega) / get_flux(vega)
        """
        return self._get_mean_and_samples_attribute("lphot")

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
        mean, samples = self._get_mean_and_samples_attribute(
            "get_Nphotons", slamb, sflux, axis=axis
        )
        return mean, samples

    @property
    def Vega_zero_photons(self):
        """Vega number of photons per wavelength unit

        .. note::

            see `self.get_Nphotons`
        """
        return self._get_mean_and_samples_attribute("Vega_zero_photons")

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
        mean, samples = self._get_mean_and_samples_attribute(
            "getFlux", slamb, sflux, axis=axis
        )
        return mean, samples

    def reinterp(self, lamb):
        """reinterpolate filter onto a different wavelength definition"""
        mean, samples = self._get_mean_and_samples_attribute("reinterp")
        mean_val = mean(lamb)
        samp_val = [sk(mean_val.wavelength) for sk in samples]
        samp_transmissions = [sk.transmit for sk in samp_val]

        return self.__class__(
            mean_val.wavelength,
            mean_val.transmit,
            samp_transmissions,
            name=self.name,
            dtype=mean_val.dtype,
            unit=mean_val.wavelength_unit,
        )

    def apply_transmission(self, slamb, sflux):
        """
        Apply filter transmission to a spectrum
        (with reinterpolation of the filter)

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
        mean, samples = self._get_mean_and_samples_attribute("apply_transmission")
        mean_val = mean(slamb, sflux)
        samp_val = [sk(slamb, sflux) for sk in samples]
        return mean_val, samp_val

    @property
    def AB_zero_mag(self):
        """AB magnitude zero point

        ABmag = -2.5 * log10(f_nu) - 48.60
              = -2.5 * log10(f_lamb) - 2.5 * log10(lpivot ** 2 / c) - 48.60
              = -2.5 * log10(f_lamb) - zpts

        """
        return self._get_mean_and_samples_attribute("AB_zero_mag")

    @property
    def AB_zero_flux(self):
        """AB flux zero point in erg/s/cm2/AA"""
        return self._get_mean_and_samples_attribute("AB_zero_flux")

    @property
    def AB_zero_Jy(self):
        """AB flux zero point in Jansky (Jy)"""
        return self._get_mean_and_samples_attribute("AB_zero_Jy")

    @property
    def Vega_zero_mag(self):
        """Vega magnitude zero point
        Vegamag = -2.5 * log10(f_lamb) + 2.5 * log10(f_vega)
        Vegamag = -2.5 * log10(f_lamb) - zpts
        """
        return self._get_mean_and_samples_attribute("Vega_zero_mag")

    @property
    def Vega_zero_flux(self):
        """Vega flux zero point in erg/s/cm2/AA"""
        return self._get_mean_and_samples_attribute("Vega_zero_flux")

    @property
    def Vega_zero_Jy(self):
        """Vega flux zero point in Jansky (Jy)"""
        return self._get_mean_and_samples_attribute("Vega_zero_Jy")

    @property
    def ST_zero_mag(self):
        """ST magnitude zero point
        STmag = -2.5 * log10(f_lamb) -21.1
        """
        return 21.1

    @property
    def ST_zero_flux(self):
        """ST flux zero point in erg/s/cm2/AA"""
        return 10 ** (-0.4 * self.ST_zero_mag) * Unit("erg*s-1*cm-2*AA-1")

    @property
    def ST_zero_Jy(self):
        """ST flux zero point in Jansky (Jy)"""
        return self._get_mean_and_samples_attribute("ST_zero_Jy")

    def to_Table(self, **kwargs):
        """Export filter to a SimpleTable object

        Parameters
        ----------
        fname: str
            filename

        Uses `SimpleTable` parameters
        """
        mean_transmit, transmit_ = self.transmit
        data_ = {"WAVELENGTH": self._wavelength, "THROUGHPUT": mean_transmit}
        for num, filterk in enumerate(transmit_, 1):
            data_["THROUGHPUT_{0:d}".format(num)] = filterk
        data = SimpleTable(data_)

        if self.wavelength_unit is not None:
            data.header["WAVELENGTH_UNIT"] = self.wavelength_unit
        data.header["DETECTOR"] = self.dtype
        data.header["COMPNAME"] = self.name
        data.header["NAME"] = self.name
        data.set_comment("THROUGHPUT", "filter throughput definition")
        data.set_comment("WAVELENGTH", "filter wavelength definition")
        for num in range(1, len(transmit_) + 1):
            data.set_comment("THROUGHPUT_{0:d}".format(num), "filter throughput sample")
        data.set_comment("WAVELENGTH", self.wavelength_unit or "AA")
        return data

    @classmethod
    def from_ascii(cls, fname, dtype="csv", **kwargs):
        """Load filter from ascii file"""
        lamb = kwargs.pop("lamb", None)
        name = kwargs.pop("name", None)
        detector = kwargs.pop("detector", "photon")
        unit_ = kwargs.pop("unit", None)

        if not isinstance(fname, SimpleTable):
            t = SimpleTable(fname, dtype=dtype, **kwargs)
        else:
            t = fname
        w = t["WAVELENGTH"].astype(float)
        r = t["THROUGHPUT"].astype(float)
        keys = [k for k in t.keys() if "THROUGHPUT_" in k]

        # update properties from file header
        detector = t.header.get("DETECTOR", detector)
        unit_ = t.header.get("WAVELENGTH_UNIT", unit_)

        # try from the comments in the header first
        if name in (None, "None", "none", ""):
            name = [
                k.split()[1]
                for k in t.header.get("COMMENT", "").split("\n")
                if "COMPNAME" in k
            ]
            name = "".join(name).replace('"', "").replace("'", "")
        # if that did not work try the table header directly
        if name in (None, "None", "none", ""):
            name = t.header["NAME"]

        if len(keys) > 0:
            samp = np.array([t[key] for key in keys])
            _filter = cls(w, r, samp, name=name, dtype=detector, unit=unit_)
        else:
            _filter = UnitFilter(w, r, name=name, dtype=detector, unit=unit_)

        # reinterpolate if requested
        if lamb is not None:
            _filter = _filter.reinterp(lamb)

        return _filter


class UnitLibrary(object):
    """Common grounds for filter libraries"""

    def __init__(self, source=__default__, *args, **kwargs):
        """Construct the library"""
        self.source = None

    def __repr__(self):
        msg = "Filter Library: {0}\n{1:s}"
        return msg.format(self.source, object.__repr__(self))

    def __enter__(self):
        """Enter context"""
        return self

    def __exit__(self, *exc_info):
        """end context"""
        return False

    def __len__(self):
        """Size of the library"""
        return len(self.content)

    def to_csv(self, directory="./", progress=True, **kwargs):
        """Export each filter into a csv file with its own name

        Parameters
        ----------
        directory: str
            directory to write into
        progress: bool
            show progress if set

        """
        from .helpers import progress_enumerate

        try:
            os.stat(directory)
        except Exception:
            os.mkdir(directory)
        with self as s:
            for _, k in progress_enumerate(
                s.content, desc="export", show_progress=progress
            ):
                f = s[k]
                if f.wavelength_unit is None:
                    f.wavelength_unit = "AA"
                f.write_to(
                    "{0:s}/{1:s}.csv".format(directory, f.name).lower(),
                    fmt="%.6f",
                    **kwargs,
                )

    def to_hdf(self, fname="filters.hd5", progress=True, **kwargs):
        """Export each filter into a csv file with its own name

        Parameters
        ----------
        directory: str
            directory to write into
        progress: bool
            show progress if set

        """
        from .helpers import progress_enumerate

        with self as s:
            for _, k in progress_enumerate(
                s.content, desc="export", show_progress=progress
            ):
                f = s[k]
                if f.wavelength_unit is None:
                    f.wavelength_unit = "AA"
                f.write_to(
                    "{0:s}".format(fname),
                    tablename="/filters/{0}".format(f.name),
                    createparents=True,
                    append=True,
                    silent=True,
                    **kwargs,
                )

    @classmethod
    def from_hd5(cls, filename, **kwargs):
        return UnitHDF_Library(filename, **kwargs)

    @classmethod
    def from_ascii(cls, filename, **kwargs):
        return UnitAscii_Library(filename, **kwargs)

    @property
    def content(self):
        """Get the content list"""
        return self.get_library_content()

    def __getitem__(self, name):
        """Make this object like a dictionary and load one or multiple filters"""
        with self as s:
            try:
                f = s._load_filter(name)
            except TypeError:
                f = [s._load_filter(k) for k in name]
        return f

    def _load_filter(self, *args, **kwargs):
        """Load a given filter from the library"""
        raise NotImplementedError

    def get_library_content(self):
        """get the content of the library"""
        raise NotImplementedError

    def load_all_filters(self, interp=True, lamb=None):
        """load all filters from the library"""
        raise NotImplementedError

    def add_filter(self, f):
        """add a filter to the library"""
        raise NotImplementedError

    def find(self, name, case_sensitive=True):
        r = []
        if case_sensitive:
            _n = name.lower()
            for k in self.get_library_content():
                if _n in k.lower():
                    r.append(k)
        else:
            for k in self.content:
                if name in k:
                    r.append(k)
        return r


class UnitAscii_Library(UnitLibrary):
    """Interface one or multiple directory or many files as a filter library

    >>> lib = Ascii_Library(['ground', 'hst', 'myfilter.csv'])
    """

    def __init__(self, source):
        self.source = source

    def _load_filter(self, fname, interp=True, lamb=None, *args, **kwargs):
        """Load a given filter from the library"""
        try:
            fil = UnitFilter.from_ascii(fname, *args, **kwargs)
        except Exception:
            content = self.content
            r = [k for k in content if fname in k]

            if len(r) <= 0:  # try all lower for filenames (ascii convention)
                r = [k for k in content if fname.lower() in k]

            if len(r) > 1:
                print("auto correction found multiple choices")
                print(r)
                raise ValueError("Refine name to one of {0}".format(r))
            elif len(r) <= 0:
                raise ValueError("Cannot find filter {0}".format(fname))
            else:
                fil = UnitFilter.from_ascii(r[0], *args, **kwargs)
        if (interp is True) and (lamb is not None):
            return fil.reinterp(lamb)
        else:
            return fil

    def get_library_content(self):
        """get the content of the library"""
        from glob import glob

        try:
            os.path.isdir(self.source)
            lst = glob(self.source + "/*")
        except TypeError:
            lst = [self.source]
        dircheck = True
        while dircheck is True:
            dircheck = False
            newlst = []
            for entry in lst:
                if os.path.isdir(entry):
                    newlst.extend(glob(entry + "/*"))
                    dircheck = True
                else:
                    newlst.append(entry)
            lst = newlst
        return lst

    def load_all_filters(self, interp=True, lamb=None):
        """load all filters from the library"""
        return [self._load_filter(k, interp=interp, lamb=lamb) for k in self.content]

    def load_filters(self, names, interp=True, lamb=None, filterLib=None):
        """load a limited set of filters

        Parameters
        ----------
        names: list[str]
            normalized names according to filtersLib

        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        filterLib: path
            path to the filter library hd5 file

        Returns
        -------
        filters: list[filter]
            list of filter objects
        """
        filters = [
            self._load_filter(fname, interp=interp, lamb=lamb) for fname in names
        ]
        return filters

    def add_filters(self, filter_object, fmt="%.6f", **kwargs):
        """Add a filter to the library permanently

        Parameters
        ----------
        filter_object: Filter object
            filter to add
        """
        if not isinstance(filter_object, UnitFilter):
            msg = "Argument of type Filter expected. Got type {0}"
            raise TypeError(msg.format(type(filter_object)))

        if filter_object.wavelength_unit is None:
            msg = "Filter wavelength must have units for storage."
            raise AttributeError(msg)
        fname = "{0:s}/{1:s}.csv".format(self.source, filter_object.name)
        filter_object.write_to(fname.lower(), fmt=fmt, **kwargs)


class UnitHDF_Library(UnitLibrary):
    """Storage based on HDF"""

    def __init__(self, source=__default__, mode="r"):
        self.source = source
        self.hdf = None
        self.mode = mode
        self._in_context = 0

    def __enter__(self):
        """Enter context"""
        if self.hdf is None:
            self.hdf = tables.open_file(self.source, self.mode)
        self._in_context += 1
        return self

    def __exit__(self, *exc_info):
        """end context"""
        if (self.hdf is not None) and (self._in_context < 2):
            self.hdf.close()
            self.hdf = None
        self._in_context -= 1
        return False

    def _load_filter(self, fname, interp=True, lamb=None):
        """Load a given filter from the library

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
        if hasattr(fname, "decode"):
            fnode = ftab.get_node("/filters/" + fname.decode("utf8"))
        else:
            fnode = ftab.get_node("/filters/" + fname)
        flamb = fnode[:]["WAVELENGTH"]
        transmit = fnode[:]["THROUGHPUT"]
        dtype = "photon"
        unit = None

        attrs = fnode.attrs
        if "DETECTOR" in attrs:
            dtype = attrs["DETECTOR"]
        if "WAVELENGTH_UNIT" in attrs:
            unit = attrs["WAVELENGTH_UNIT"]

        fil = UnitFilter(flamb, transmit, name=fnode.name, dtype=dtype, unit=unit)

        if interp & (lamb is not None):
            fil = fil.reinterp(lamb)
        return fil

    def get_library_content(self):
        """get the content of the library"""
        with self as s:
            try:
                filters = s.hdf.root.content.cols.TABLENAME[:]
            except Exception:
                filters = list(s.hdf.root.filters._v_children.keys())
        if hasattr(filters[0], "decode"):
            filters = [k.decode("utf8") for k in filters]
        return filters

    def load_all_filters(self, interp=True, lamb=None):
        """load all filters from the library

        Parameters
        ----------
        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        Returns
        -------
        filters: list[filter]
            list of filter objects
        """
        with self as s:
            filters = [
                s._load_filter(fname, interp=interp, lamb=lamb) for fname in s.content
            ]
        return filters

    def load_filters(self, names, interp=True, lamb=None, filterLib=None):
        """load a limited set of filters

        Parameters
        ----------
        names: list[str]
            normalized names according to filtersLib

        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        filterLib: path
            path to the filter library hd5 file

        Returns
        -------
        filters: list[filter]
            list of filter objects
        """
        with self as s:
            filters = [
                s._load_filter(fname, interp=interp, lamb=lamb) for fname in names
            ]
        return filters

    def add_filter(self, f, **kwargs):
        """Add a filter to the library permanently

        Parameters
        ----------
        f: Filter object
            filter to add
        """
        if not isinstance(f, UnitFilter):
            msg = "Argument of type Filter expected. Got type {0}"
            raise TypeError(msg.format(type(f)))

        if f.wavelength_unit is None:
            msg = "Filter wavelength must have units for storage."
            raise AttributeError(msg)

        append = kwargs.pop("append", True)

        f.write_to(
            "{0:s}".format(self.source),
            tablename="/filters/{0}".format(f.name),
            createparents=True,
            append=append,
            **kwargs,
        )


def get_library(fname=__default__, **kwargs):
    """Finds the appropriate class to load the library"""
    if os.path.isfile(fname):
        return UnitHDF_Library(fname, **kwargs)
    else:
        return UnitAscii_Library(fname, **kwargs)


def _reduce_resolution(wi, fi, fwhm0=0.55, sigma_floor=0.2):
    """Adapt the resolution of the spectra to match the lick definitions

        Lick definitions have different resolution elements as function
        of wavelength. These definition are hard-coded in this function

    Parameters
    ----------
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

    # all in AA
    w_lick_res = (4000.0, 4400.0, 4900.0, 5400.0, 6000.0)
    lick_res = (11.5, 9.2, 8.4, 8.4, 9.8)  # FWHM in AA

    w = np.asarray(wi)
    flux = np.atleast_2d(fi)

    # Linear interpolation of lick_res over w
    # numpy interp does constant instead of extrapolation
    # res = np.interp(w, w_lick_res, lick_res)

    # spline order: 1 linear, 2 quadratic, 3 cubic ...
    from scipy.interpolate import InterpolatedUnivariateSpline

    res = InterpolatedUnivariateSpline(w_lick_res, lick_res, k=1)(w)

    # Compute width from fwhm
    const = 2.0 * np.sqrt(2.0 * np.log(2))  # conversion fwhm --> sigma
    lick_sigma = np.sqrt((res**2 - fwhm0**2)) / const

    # Convolution by g=1/sqrt(2*pi*sigma^2) * exp(-r^2/(2*sigma^2))
    flux_red = np.zeros(flux.shape, dtype=flux.dtype)

    for i, sigma in enumerate(lick_sigma):
        maxsigma = 3.0 * sigma
        # sampling floor: min (0.2, sigma * 0.1)
        delta = min(sigma_floor, sigma * 0.1)
        delta_wj = np.arange(-maxsigma, +maxsigma, delta)
        wj = delta_wj + w[i]
        for k, fk in enumerate(flux):
            fluxj = np.interp(wj, w, fk, left=0.0, right=0.0)
            flux_red[k, i] = np.sum(
                fluxj * delta * np.exp(-0.5 * (delta_wj / sigma) ** 2)
            )

    flux_red /= lick_sigma * const

    return flux_red.reshape(np.shape(fi))


@set_method_default_units("AA", "flam", output_unit="flam")
def reduce_resolution(wi, fi, fwhm0=0.55 * Unit("AA"), sigma_floor=0.2 * Unit("AA")):
    """Adapt the resolution of the spectra to match the lick definitions

        Lick definitions have different resolution elements as function
        of wavelength. These definition are hard-coded in this function

    Parameters
    ----------
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
    flux_red = _reduce_resolution(
        wi.value, fi.value, fwhm0.to("AA").value, sigma_floor.to("AA").value
    )
    return flux_red * Unit("flam")


class UnitLickIndex(object):
    """Define a Lick Index similarily to a Filter object"""

    def __init__(self, name, lick, unit="AA"):
        """Constructor

        Parameters
        ----------
        name: str
            name of the index
        lick: dict
            expecting 'blue', 'red', 'band', and 'unit' definitions
            `blue` and `red` are used to continuum normalize the spectra
            `band` covers the index itself. `unit` gives the index measurement
            units, either magnitudes (mag) or equivalent width (ew)
        unit: str
            wavelength unit of the intervals
        """
        self.name = name
        self._lick = lick
        self.wavelength_unit = unit

    def to_dict(self):
        """return a dictionary of the current index"""
        d = {}
        d.update(**self._lick)
        return d

    def _get_wavelength_attrs_with_units(self, attrname, units="AA"):
        """return the unitwise definition corresponding to attrname"""
        attr = self._lick[attrname]
        if self.wavelength_unit is not None:
            if units is None:
                return attr * Unit(self.wavelength_unit)
            else:
                return (attr * Unit(self.wavelength_unit)).to(units)
        else:
            return attr

    @property
    def band(self):
        """Unitwise band definition"""
        return self._get_wavelength_attrs_with_units("band")

    @property
    def blue(self):
        """Unitwise band definition"""
        return self._get_wavelength_attrs_with_units("blue")

    @property
    def red(self):
        """Unitwise band definition"""
        return self._get_wavelength_attrs_with_units("red")

    @property
    def index_unit(self):
        return self._lick["unit"]

    def __repr__(self):
        txt = """LickIndex ({0}), {1}"""
        return txt.format(self.name, object.__repr__(self))

    def info(self):
        """display information about the current Index"""
        txt = """Lick Index {s.name}
    wavelength units:     {s.wavelength_unit}
    Index Band:           {s.band}
    Blue continuum band:  {s.blue}
    Red continuum band:   {s.red}
    Measurement unit:     {s.index_unit}""".format(s=self)
        print(txt)

    def __call__(self, *args, **kwargs):
        """compute spectral index after continuum subtraction

        Parameters
        ----------
        w: ndarray (nw, )
            array of wavelengths in AA
        flux: ndarray (N, nw)
            array of flux values for different spectra in the series
        degree: int (default 1)
            degree of the polynomial fit to the continuum

        Returns
        -------
        ew: ndarray (N,)
            equivalent width or magnitude array
        """
        return self.get(*args, **kwargs)

    def _get(self, wave, flux, **kwargs):
        """compute spectral index after continuum subtraction

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
        if hasUnit(wave):
            _w = wave.to("AA").value
        else:
            print("Warning: assuming units are in Angstroms")
            _w = _drop_units(wave)
        _f = _drop_units(flux)

        blue = self._get_wavelength_attrs_with_units("blue").value
        red = self._get_wavelength_attrs_with_units("red").value
        band = self._get_wavelength_attrs_with_units("band").value

        nocheck = kwargs.pop("nocheck", False)
        not_covered = (blue[0] < _w[0]) | (red[-1] > _w[-1])
        if not_covered:
            if not nocheck:
                raise ValueError("Spectrum does not cover this index.")
            else:
                return np.zeros(_f.shape[0]) * float("nan")
        else:
            return self._get_indice(_w, _f, blue, red, band, self.index_unit, **kwargs)

    @classmethod
    def _get_indice(cls, w, flux, blue, red, band=None, unit="ew", degree=1, **kwargs):
        """compute spectral index after continuum subtraction

        Parameters
        ----------
        w: ndarray (nw, )
            array of wavelengths in AA
        flux: ndarray (N, nw)
            array of flux values for different spectra in the series
        blue: tuple(2)
            selection for blue continuum estimate
        red: tuple(2)
            selection for red continuum estimate
        band: tuple(2), optional
            select region in this band only.
            default is band = (min(blue), max(red))
        unit: str
            `ew` or `mag` wether equivalent width or magnitude
        degree: int (default 1)
            degree of the polynomial fit to the continuum

        Returns
        -------
        ew: ndarray (N,)
            equivalent width array
        """
        wi, fi = cls.continuum_normalized_region_around_line(
            w, flux, blue, red, band=band, degree=degree
        )
        if unit in (0, "ew", "EW"):
            return trapezoid(1.0 - fi, wi, axis=-1)
        else:
            m = trapezoid(fi, wi, axis=-1)
            m = -2.5 * np.log10(m / np.ptp(wi))
            return m

    @classmethod
    def continuum_normalized_region_around_line(
        cls, wi, fi, blue, red, band=None, degree=1
    ):
        """
        cut out and normalize flux around a line

        Parameters
        ----------
        wi: ndarray (nw, )
            array of wavelengths in AA
        fi: ndarray (N, nw)
            array of flux values for different spectra in the series
        blue: tuple(2)
            selection for blue continuum estimate
        red: tuple(2)
            selection for red continuum estimate
        band: tuple(2), optional
            select region in this band only.
            default is band = (min(blue), max(red))
        degree: int
            degree of the polynomial fit to the continuum

        returns
        -------
        wnew: ndarray (nw1, )
            wavelength of the selection in AA
        f: ndarray (N, len(wnew))
            normalized flux in the selection region

        Example
        -------

        .. code-block:: python

            # indice of CaII
            # wavelength are always supposed in AA
            w, f = region_around_line(
                wavelength, flux, [3925, 3930],[3938, 3945]]
                )

        """
        w = np.asarray(wi)
        flux = np.atleast_2d(fi)
        # index is true in the region where we fit the polynomial
        indcont = ((w >= blue[0]) & (w <= blue[1])) | ((w >= red[0]) & (w <= red[1]))
        # index of the region we want to return
        if band is None:
            band = blue[0], red[1]

        indrange = (w > band[0]) & (w < band[1])
        wnew = w[indrange]
        wcont = w[indcont]

        # make a flux array of shape
        # (number of spectra, number of points in indrange)
        f = np.zeros((flux.shape[0], indrange.sum()))
        for i in range(flux.shape[0]):
            # fit polynomial of second order to the continuum region
            linecoeff = np.polyfit(wcont, flux[i, indcont], degree)
            # divide the flux by the polynomial and put the result in our new
            # flux array
            f[i, :] = flux[i, indrange] / np.polyval(linecoeff, wnew)
        return wnew, np.squeeze(f)

    @set_method_default_units("AA", "flam")
    def get(self, wave, flux, **kwargs):
        """compute spectral index after continuum subtraction

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
        return self._get(wave, flux.to("flam").value, **kwargs)


class UnitLickLibrary(object):
    """Collection of Lick indices"""

    def __init__(self, fname=__default_lick__, comment="#"):
        self.source = fname
        data, hdr = self._read_lick_list(fname, comment)
        self._content = data
        self._hdr = hdr

    @property
    def description(self):
        """any comment in the input file"""
        return self._hdr

    @classmethod
    def _read_lick_list(cls, fname=__default_lick__, comment="#"):
        """read the list of lick indices

        Parameters
        ----------
        fname: str
            file containing the indices' definitions
        comment: str
            character indicating comment in the file

        Returns
        -------
        data: dict
            dictionary of indices
            name: (band, blue, red, unit)
        """
        with open(fname, "r") as f:
            data = {}
            hdr = []
            for line in f:
                if line[0] != comment:
                    _line = line.split()
                    attr = dict(
                        band=(float(_line[1]), float(_line[2])),
                        blue=(float(_line[3]), float(_line[4])),
                        red=(float(_line[5]), float(_line[6])),
                        unit="mag" if int(_line[7]) > 0 else "ew",
                    )
                    name = _line[8]
                    data[name] = attr
                else:
                    hdr.append(line[1:-1])
        return data, hdr

    def __repr__(self):
        return "Lick Index Library: {0}\n{1:s}".format(
            self.source, object.__repr__(self)
        )

    def __enter__(self):
        """Enter context"""
        return self

    def __exit__(self, *exc_info):
        """end context"""
        return False

    def __len__(self):
        """Size of the library"""
        return len(self.content)

    def get_library_content(self):
        return list(self._content.keys())

    def __getitem__(self, name):
        """Make this object like a dictionary and load one or multiple filters"""
        with self as s:
            try:
                f = s._load_filter(name)
            except TypeError:
                f = [s._load_filter(k) for k in name]
        return f

    def _load_filter(self, fname, **kwargs):
        """Load a given filter from the library"""
        with self as current_lib:
            return UnitLickIndex(fname, current_lib._content[fname])

    @property
    def content(self):
        return self.get_library_content()

    def find(self, name, case_sensitive=True):
        r = []
        if not case_sensitive:
            _n = name.lower()
            for k in self.get_library_content():
                if _n in k.lower():
                    r.append(k)
        else:
            for k in self.content:
                if name in k:
                    r.append(k)
        return r
