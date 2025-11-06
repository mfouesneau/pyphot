"""
Photometric package
===================

This module defines the `Filter` class.

.. note::

    integrations are done using :func:`trapezoid`
    Why not Simpsons? Simpsons principle is to take sequence of 3 points to
    make a quadratic interpolation. Which in the end, when filters have sharp
    edges, the error due to this "interpolation" are extremely large in
    comparison to the uncertainties induced by trapeze integration.
"""

from typing import Literal, Union, Optional, Any, Dict, cast
import warnings
from io import IOBase
from os import PathLike

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.integrate import trapezoid

from . import config
from .constants import Constants
from .unit_adapters import QuantityType, enforce_default_units
from .vega import Vega
from .io.ascii import from_ascii, from_csv
from .io.hdf import from_hdf5
from .io import to_file

__all__ = ["Filter"]


def _drop_units(q: Any) -> Any:
    """Drop the unit definition silently"""
    try:
        return q.value
    except AttributeError:
        try:
            return q.value
        except AttributeError:
            return q


def _split_value_unit(q: Union[QuantityType, Any]) -> tuple[Any, Union[str, None]]:
    """Split a quantity into value and unit"""
    try:
        return q.value, str(q.unit)
    except AttributeError:
        return q, None


class Filter:
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

    vega: str
        Vega flavor to use for calculations, default is 'default'
        (see :class:`pyphot.vega.Vega` for details)
    """

    dtype: Literal["photon", "energy"]
    name: str
    norm: float
    transmit: npt.NDArray[np.floating]
    wavelength_unit: str

    _cl: float
    _lT: float
    _lpivot: float
    _vega_flavor: str
    _wavelength: npt.NDArray[np.floating]

    def __init__(
        self,
        wavelength: Union[npt.NDArray[np.floating], QuantityType],
        transmit: npt.NDArray[np.floating],
        *,
        name: str = "",
        dtype: Literal["photon", "energy"] = "photon",
        unit: Optional[str] = None,
        vega: Optional[str] = None,
    ):
        """Constructor"""
        self.name = name
        self.set_dtype(dtype)

        # set wavelength definition
        λ_, u_ = _split_value_unit(wavelength)
        self._wavelength = λ_
        self.set_wavelength_unit(u_ or unit)

        # make sure input data are ordered and cleaned of weird values.
        idx = np.argsort(self._wavelength)
        self._wavelength = self._wavelength[idx]
        self.transmit = np.clip(transmit[idx], 0.0, np.nanmax(transmit))

        # check for NaN values in input arrays
        if np.any(np.isnan(self._wavelength)):
            raise ValueError("Wavelength array contains NaN values")
        if np.any(np.isnan(self.transmit)):
            raise ValueError("Transmit array contains NaN values")

        # set vega flavor
        self.set_vega_flavor(vega, reset=False)

        self._reset_attributes()

    def set_dtype(self, dtype: Literal["photon", "energy"]):
        """Set the detector type (photon or energy)"""
        _d = dtype.lower()
        if "phot" in _d:
            self.dtype = "photon"
        elif "ener" in _d:
            self.dtype = "energy"
        else:
            raise ValueError("Unknown detector type {0}".format(dtype))

    def set_wavelength_unit(self, unit: Optional[str]):
        """Set the wavelength units"""
        if unit in (None, ""):
            raise ValueError(
                "Wavelength unit must be specified either in the wavelength array or as the `unit` keyword argument"
            )
        # check if unit is valid
        _ = config.units.U(unit)
        self.wavelength_unit = unit

    def set_vega_flavor(self, vega: Optional[str], reset: bool = True):
        """Set the Vega flavor to use for calculations."""
        if vega is None:
            vega = config.__vega_default_flavor__
        self._vega_flavor = vega
        if reset:
            self._reset_attributes()

    def info(self, show_zeropoints: bool = True):
        """print information about the current filter"""
        print(self._get_info(show_zeropoints))

    def reinterp(self, lamb: Union[npt.NDArray[np.floating], QuantityType]) -> "Filter":
        """reinterpolate filter onto a different wavelength definition

        Parameters
        ----------
        lamb : QuantityType | npt.NDArray[np.floating]
            Wavelength to reinterpolate onto
            If no unit is provided, the filter's wavelength unit is assumed.

        Returns
        -------
        Filter
            Filter reinterpolated onto the new wavelength definition
        """
        _lamb, _unit_lamb = _split_value_unit(lamb)
        _unit_lamb = _unit_lamb or self.wavelength_unit  # defaulting to filter's unit
        # make sure the interpolation target has units
        _wavelength = self._get_filter_wavelength_in_units_of(
            _lamb * config.units.U(_unit_lamb)
        )

        ifT = np.interp(_lamb, _wavelength, self.transmit, left=0.0, right=0.0)
        return self.__class__(
            _lamb,
            ifT,
            name=self.name,
            dtype=self.dtype,
            unit=_unit_lamb,
            vega=self._vega_flavor,
        )

    def apply_transmission(self, slamb: QuantityType, sflux: QuantityType):
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
        _wavelength = self._get_filter_wavelength_in_units_of(slamb)
        _lamb = _drop_units(slamb)
        ifT = np.interp(_lamb, _wavelength, self.transmit, left=0.0, right=0.0)
        return ifT * sflux

    @property
    def wavelength(self) -> QuantityType:
        """Unitwise wavelength definition"""
        if self.wavelength_unit is not None:
            return self._wavelength * config.units.U(self.wavelength_unit)
        else:
            raise ValueError("Wavelength unit is not defined")

    @property
    def lmax(self) -> QuantityType:
        """Calculated as the last value with a transmission at least 1% of
        maximum transmission"""
        cond = (self.transmit / self.transmit.max()) > 1.0 / 100
        w = cast(npt.NDArray, self.wavelength)
        return max(w[cond])

    @property
    def lmin(self) -> QuantityType:
        """Calculate das the first value with a transmission at least 1% of
        maximum transmission"""
        cond = (self.transmit / self.transmit.max()) > 1.0 / 100
        w = cast(npt.NDArray, self.wavelength)
        return min(w[cond])

    @property
    def width(self) -> QuantityType:
        r"""Effective width
        Equivalent to the horizontal size of a rectangle with height equal
        to maximum transmission and with the same area that the one covered by
        the filter transmission curve.

        .. math::
            W = \int_{\lambda_{min}}^{\lambda_{max}} T(\lambda) d\lambda / \max(T)

        """
        return (self.norm / max(self.transmit)) * config.units.U(self.wavelength_unit)

    @property
    def fwhm(self) -> QuantityType:
        """the difference between the two wavelengths for which filter
        transmission is half maximum

        ..note::
            This calculation is not exact but rounded to the nearest passband
            data points
        """
        vals = self.transmit / self.transmit.max() - 0.5
        zero_crossings = np.where(np.diff(np.sign(vals)))[0]
        lambs = cast(npt.NDArray, self.wavelength[zero_crossings])
        return np.diff(lambs)[0]

    @property
    def lpivot(self) -> QuantityType:
        r"""Unitwise wavelength definition

        the calculation depends on the type (`dtype`) of the filter.

        for `photon` filters:

        ..math::
            \lambda_{pivot}^2 = \frac{\int \lambda T(\lambda) d\lambda}{\int T(\lambda) d\lambda / \lambda}

        for `energy` filters:

        ..math::
            \lambda_{pivot}^2 = \frac{\int T(\lambda) d\lambda}{\int T(\lambda) d\lambda / \lambda^2}

        """
        if self.wavelength_unit is not None:
            return self._lpivot * config.units.U(self.wavelength_unit)
        else:
            raise ValueError("Wavelength unit is not defined")

    @property
    def cl(self) -> QuantityType:
        r"""Unitwise Central wavelength

        ..math::
            \lambda_{cl} = \frac{\int \lambda T(\lambda) d\lambda}{\int T(\lambda) d\lambda}
        """
        if self.wavelength_unit is not None:
            return self._cl * config.units.U(self.wavelength_unit)
        else:
            raise ValueError("Wavelength unit is not defined")

    @property
    def leff(self) -> QuantityType:
        """Unitwise Effective wavelength
        leff = int (lamb * T * Vega dlamb) / int(T * Vega dlamb)
        """
        U = config.units.U
        with Vega(flavor=self._vega_flavor) as v:
            s = self.reinterp(v.wavelength)
            w = s._wavelength
            if s.transmit.max() > 0:
                leff = trapezoid(w * s.transmit * v.flux.value, w, axis=-1)
                leff /= trapezoid(s.transmit * v.flux.value, w, axis=-1)
            else:
                leff = float("nan")
        if (s.wavelength_unit is not None) and (self.wavelength_unit is not None):
            leff = cast(QuantityType, leff * U(s.wavelength_unit))
            return leff.to(self.wavelength_unit)
        else:
            raise ValueError("Wavelength unit is not defined")

    @property
    def lphot(self) -> QuantityType:
        r"""Photon distribution based effective wavelength. Defined as

        .. math::

            lphot = \int \lambda^2 T \cdot Vega(\lambda) d\lambda / \int \lambda T \cdot Vega(\lambda) d\lambda

        which we calculate as

        lphot = get_flux(lamb * vega) / get_flux(vega)
        """
        if self.wavelength_unit is None:
            raise AttributeError("Needs wavelength units")

        with Vega(flavor=self._vega_flavor) as v:
            # check overlap with vega spectrum
            λv = v.wavelength.to("AA").value
            λf = self.wavelength.to("AA").value

            if (λv.max() < λf.min()) or (λv.min() > λf.max()):
                warnings.warn(
                    f"Vega spectrum ({v.source}) does not overlap with filter"
                )
                return np.nan * config.units.U(self.wavelength_unit)

            wave = v.wavelength.value
            # Cheating units to avoid making a new filter
            f_vega = self.get_flux(v.wavelength, v.flux, axis=-1)
            f_lamb_vega = self.get_flux(v.wavelength, wave * v.flux, axis=-1)
            f_lamb2_vega = self.get_flux(v.wavelength, wave**2 * v.flux, axis=-1)
        if "photon" in self.dtype:
            lphot = f_lamb_vega / f_vega
        else:
            lphot = f_lamb2_vega / f_lamb_vega
        return (lphot * config.units.U(str(v.wavelength.unit))).to(self.wavelength_unit)

    @property
    def AB_zero_mag(self) -> float:
        r"""AB magnitude zero point

        ABmag = -2.5 * log10(f_nu) - 48.60
              = -2.5 * log10(f_lamb) - 2.5 * log10(lpivot ** 2 / c) - 48.60
              = -2.5 * log10(f_lamb) - zpts

        """
        if self.wavelength_unit is None:
            raise AttributeError("Needs wavelength units")

        c_ = Constants.get("c").to("AA/s")
        C1 = config.units.Q(self.wavelength_unit).to("AA") ** 2 / c_
        c1 = _drop_units(self._lpivot**2 * C1)

        m = 2.5 * np.log10(c1) + 48.6
        return m

    @property
    def AB_zero_flux(self) -> QuantityType:
        """AB flux zero point in erg/s/cm2/AA"""
        unit_ = config.units.U("erg*s**-1*cm**-2*AA**-1")
        return 10 ** (-0.4 * self.AB_zero_mag) * unit_

    @property
    def AB_zero_Jy(self) -> QuantityType:
        r"""AB flux zero point in Jansky (Jy)

        .. math::

            {f_{Jy}} = \frac{10^5}{10^{-8}c} {\lambda_p^2} {f_\lambda}
        """
        c = 1e-8 * Constants.get("c").to("m/s").value
        f = 1e5 / c * self.lpivot.to("AA").value ** 2 * self.AB_zero_flux.value
        return f * config.units.U("Jy")

    @property
    def Vega_zero_photons(self) -> QuantityType:
        """Vega number of photons per wavelength unit

        .. note::

            see `self.get_Nphotons`
        """
        with Vega(flavor=self._vega_flavor) as v:
            return self.get_Nphotons(v.wavelength, v.flux)

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
        c = 1e-8 * Constants.get("c").to("m/s").value
        vega_zero_flux = self.Vega_zero_flux.to("erg*s**-1*cm**-2*AA**-1").value
        lpivot = self.lpivot.to("AA").value ** 2
        f = 1e5 / c * (lpivot * vega_zero_flux)
        return f * config.units.U("Jy")

    @property
    def ST_zero_mag(self) -> float:
        """ST magnitude zero point
        STmag = -2.5 * log10(f_lamb) -21.1
        """
        return 21.1

    @property
    def ST_zero_flux(self) -> QuantityType:
        """ST flux zero point in erg/s/cm2/AA"""
        unit_ = config.units.U("erg*s**-1*cm**-2*AA**-1")
        return 10 ** (-0.4 * self.ST_zero_mag) * unit_

    @property
    def ST_zero_Jy(self) -> QuantityType:
        """ST flux zero point in Jansky (Jy)"""
        c = 1e-8 * Constants.get("c").to("m/s").value
        f = 1e5 / c * self.lpivot.to("AA").value ** 2 * self.ST_zero_flux.value
        return f * config.units.U("Jy")

    @enforce_default_units(None, "AA", "flam", output="erg*s**-1*cm**-2*AA**-1")
    def get_flux(
        self,
        slamb: QuantityType,
        sflux: QuantityType,
        axis: int = -1,
    ) -> QuantityType:
        """Get integrated flux through the filter

        Parameters
        ----------
        slamb: Union[npt.NDArray[np.floating], QuantityType]
            spectrum wavelength definition domain

        sflux: Union[npt.NDArray[np.floating], QuantityType]
            associated flux or array of many fluxes

        axis: int, optional
            axis along which to integrate the flux

        Returns
        -------
        flux: float
            Energy of the spectrum within the filter
        """

        return self._get_flux(slamb, sflux, axis=axis)

    @enforce_default_units(None, "AA", "flam", output="photon*s**-1*cm**-2*AA**-1")
    def get_Nphotons(
        self,
        slamb: QuantityType,
        sflux: QuantityType,
        axis: int = -1,
    ) -> QuantityType:
        """Get integrated number of photons through the filter

        equivalent to self.get_flux(...) * leff / h / c

        Parameters
        ----------
        slamb: Union[npt.NDArray[np.floating], QuantityType]
            spectrum wavelength definition domain

        sflux: Union[npt.NDArray[np.floating], QuantityType]
            associated flux or array of many fluxes

        axis: int, optional
            axis along which to integrate the flux

        Returns
        -------
        Nphotons: QuantityType
            Integrated number of photons through the filter
            in photons / cm^2 / s / A
        """
        return self._get_n_photons(slamb, sflux, axis)

    @classmethod
    def make_integration_filter(
        cls,
        lmin: Union[float, QuantityType],
        lmax: Union[float, QuantityType],
        name: str = "",
        dtype: Literal["photon", "energy"] = "photon",
        unit: Optional[str] = None,
        **kwargs,
    ) -> "Filter":
        """Generate an heavyside filter between lmin and lmax

        Parameters
        ----------
        lmin: float or QuantityType
            Minimum wavelength of the filter
        lmax: float or QuantityType
            Maximum wavelength of the filter
        name: str, optional
            Name of the filter
        dtype: Literal["photon", "energy"], optional
            Type of the filter
        unit: Optional[str], optional
            Unit of the filter
        **kwargs
            Additional keyword arguments

        Returns
        -------
        filter: Filter
            Filter object
        """
        dyn, _units = _split_value_unit(lmax - lmin)
        λmin, _ = _split_value_unit(lmin)
        λmax, _ = _split_value_unit(lmax)
        w = np.array([λmin - 0.01 * dyn, λmin, λmax, λmax + 0.01 * dyn])
        if _units:
            w = w * config.units.U(_units)
        f = np.array([0.0, 1.0, 1.0, 0.0])
        return cls(w, f, name=name, dtype=dtype, **kwargs)

    def _reset_attributes(self):
        """reset calculated attributes"""

        self.norm = trapezoid(self.transmit, self._wavelength)
        self._lT = trapezoid(self._wavelength * self.transmit, self._wavelength)
        self._lpivot = self._calculate_lpivot()
        if self.norm > 0:
            self._cl = self._lT / self.norm
        else:
            self._cl = 0.0

    def _calculate_lpivot(self) -> float:
        r"""Calculate the pivot wavelength

        the calculation depends on the type (`dtype`) of the filter.

        for `photon` filters:

        ..math::
            \lambda_{pivot}^2 = \frac{\int \lambda T(\lambda) d\lambda}{\int T(\lambda) d\lambda / \lambda}

        for `energy` filters:

        ..math::
            \lambda_{pivot}^2 = \frac{\int T(\lambda) d\lambda}{\int T(\lambda) d\lambda / \lambda^2}
        """
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

    def _get_info(self, show_zeropoints: bool = True) -> str:
        """generate information about the current filter"""
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
        info = msg.format(s=self).replace("None", "unknown")

        # zero points only if units
        if (self.wavelength_unit is None) or (not show_zeropoints):
            return info

        info += """
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
        return info

    def _get_n_photons(
        self,
        slamb_in_AA: QuantityType,
        sflux_in_flam: QuantityType,
        axis: int = -1,
    ) -> QuantityType:
        """Get integrated number of photons through the filter

        Parameters
        ----------
        slamb: QuantityType
            spectrum wavelength definition domain assumed in AA

        sflux: QuantityType
            associated flux assumed in flam

        axis: int, optional
            axis along which to integrate the flux

        Returns
        -------
        Nphotons: QuantityType
            Integrated number of photons through the filter
        """
        passb = self.reinterp(slamb_in_AA)
        wave = passb._wavelength
        dlambda = np.diff(wave)

        h = Constants.get("h").to("erg * s").value
        c = Constants.get("c").to("AA/s").value

        vals = sflux_in_flam.value * wave * passb.transmit
        vals[~np.isfinite(vals)] = 0.0
        Nphot = 0.5 * np.sum((vals[1:] + vals[:-1]) * dlambda) / (h * c)
        Nphot = Nphot * config.units.U("photon*s**-1*cm**-2")
        return Nphot / passb.width

    def _get_flux(
        self,
        slamb_in_AA: QuantityType,
        sflux_in_flam: QuantityType,
        axis: int = -1,
    ) -> QuantityType:
        """Integrate the flux through the filter

        Parameters
        ----------
        slamb: QuantityType
            spectrum wavelength definition domain assumed in AA

        sflux: QuantityType
            associated flux assumed in flam

        axis: int, optional
            axis along which to integrate the flux

        Returns
        -------
        flux: QuantityType
            Energy of the spectrum within the filter in erg*s**-1*cm**-2*AA**-1
        """
        passb = self.reinterp(slamb_in_AA)
        ifT = passb.transmit

        _slamb = _drop_units(slamb_in_AA)
        _sflux = passb._validate_sflux(slamb_in_AA, sflux_in_flam)
        _sflux, _f_unit = _split_value_unit(_sflux)
        _f_unit = config.units.U(_f_unit)

        # if the filter is null on that wavelength range flux is then 0
        # ind = ifT > 0.
        nonzero = np.where(ifT > 0)[0]
        if nonzero.size <= 0:
            return passb._get_zero_like(sflux_in_flam) * _f_unit

        # avoid calculating many zeros
        nonzero_start = max(0, min(nonzero) - 5)
        nonzero_end = min(len(ifT), max(nonzero) + 5)
        ind = np.zeros(len(ifT), dtype=bool)
        ind[nonzero_start:nonzero_end] = True

        # check if the filter is null on that wavelength range
        # This case should not happen given the previous check but is here for safety
        if not np.any(ind):
            return passb._get_zero_like(_sflux) * _f_unit

        _sflux = np.atleast_2d(_sflux)[..., ind]
        # limit integrals to where necessary
        if "photon" in passb.dtype:
            a = trapezoid(_slamb[ind] * ifT[ind] * _sflux, _slamb[ind], axis=axis)
            b = trapezoid(_slamb[ind] * ifT[ind], _slamb[ind])
        elif "energy" in passb.dtype:
            a = trapezoid(ifT[ind] * _sflux, _slamb[ind], axis=axis)
            b = trapezoid(ifT[ind], _slamb[ind])
        else:
            raise ValueError(f"Unknown filter type {passb.dtype}")
        if np.isinf(a).any() | np.isinf(b).any():
            print(self.name, "Warn for inf value")
        return np.squeeze(a / b) * _f_unit

    @classmethod
    def _validate_sflux(
        cls,
        slamb: Union[npt.NDArray[np.floating], QuantityType],
        sflux: Union[npt.NDArray[np.floating], QuantityType],
    ) -> Union[npt.NDArray[np.floating], QuantityType]:
        """clean input spectrum for inf values by interpolation

        Parameters
        ----------
        slamb : array-like | Quantity[array]
            Wavelength array
        sflux : array-like | Quantity[array]
            Flux array

        Returns
        -------
        array-like | Quantity[array]
            Cleaned flux array with preserved units
        """
        _sflux, _sflux_unit = _split_value_unit(sflux)
        _slamb, _ = _split_value_unit(slamb)

        if np.isinf(_sflux).any():
            indinf = np.where(np.isinf(_sflux))
            indfin = np.where(np.isfinite(_sflux))
            _sflux[indinf] = np.interp(
                _slamb[indinf], _slamb[indfin], _sflux[indfin], left=0, right=0
            )
        if _sflux_unit is not None:
            return _sflux * config.units.U(_sflux_unit)
        else:
            return _sflux

    @classmethod
    def _get_zero_like(
        cls,
        sflux: Union[npt.NDArray[np.floating], QuantityType],
        axis: int = -1,
    ):
        """return a zero value corresponding to a flux calculation on sflux

        parameters
        ----------
        sflux : array-like | Quantity[array]
            Flux array
        axis : None | int | Tuple[int, ...]
            Axis or axes along which flux calculations are done. Default is -1.

        Returns
        -------
        array-like | Quantity[array]
            Zero value corresponding to a flux calculation on sflux
        """
        # weird way to define it but most robust with axis definitions
        return np.zeros_like(_drop_units(sflux)).sum(axis=axis)

    def _get_filter_wavelength_in_units_of(
        self, slamb: QuantityType
    ) -> npt.NDArray[np.floating]:
        """Return the wavelength in the units of slamb"""
        w = self.wavelength
        if config.units.has_unit(slamb) & config.units.has_unit(w):
            return w.to(str(slamb.unit)).value
        else:
            print("Warning: assuming units are consistent")
            return self._wavelength

    def __call__(self, slamb: QuantityType, sflux: QuantityType) -> QuantityType:
        """
        Apply filter transmission to a spectrum (with reinterpolation of the
        filter)

        Parameters
        ----------
        slamb: QuantityType
            spectrum wavelength definition domain

        sflux: QuantityType
            associated flux

        Returns
        -------
        flux: QuantityType
            new spectrum values accounting for the filter

        .. seealso:: :meth:`apply_transmission`
        """
        return self.apply_transmission(slamb, sflux)

    def __repr__(self):
        """string representation"""
        return "Filter: {0:s}, {1:s}".format(self.name, object.__repr__(self))

    # legacy aliases
    getFlux = get_flux
    applyTo = apply_transmission

    @classmethod
    def from_ascii(
        cls,
        fname: Union[str, PathLike, IOBase],
        *,
        dtype: str = "csv",
        **kwargs,
    ) -> "Filter":
        """Load a Filter from an ASCII file

        Parameters
        ----------
        fname: str
            path to the ASCII file

        dtype: str, optional
            type of the ASCII file (e.g. 'csv')

        kwargs: dict, optional
            additional keyword arguments to pass to the ASCII reader or to the
            Filter constructor

        Returns
        -------
        filter: Filter
            Filter object loaded from the ASCII file
        """
        # Filter constructor kwargs
        lamb = kwargs.pop("lamb", None)
        name = kwargs.pop("name", None)
        detector = kwargs.pop("detector", "photon")
        unit = kwargs.pop("unit", None)

        # parse file
        if dtype == "csv":
            df, hdr = from_csv(fname)
        elif dtype == "txt":
            df, hdr = from_ascii(fname)
        else:
            raise ValueError(f"Unknown dtype {dtype}")

        # parse the necessary data
        w = df["WAVELENGTH"].values.astype(float)
        r = df["THROUGHPUT"].values.astype(float)
        detector = hdr.header.get("DETECTOR", detector)
        unit = hdr.header.get("WAVELENGTH_UNIT", unit)
        name = hdr.header.get("NAME", name)

        # also check the header comments
        if name in (None, "None", "none", ""):
            name = [
                k.split()[1]
                for k in hdr.header.get("COMMENT", "").split("\n")
                if "COMPNAME" in k
            ]
            name = "".join(name).replace('"', "").replace("'", "")
        # if that did not work try COMPNAME in the table header directly
        if name in (None, "None", "none", ""):
            name = hdr.header.get("COMPNAME", name)

        # instanciate filter
        _filter = Filter(w, r, name=name, dtype=detector, unit=unit)
        # reinterpolate if requested
        if lamb is not None:
            _filter = _filter.reinterp(lamb)

        return _filter

    def to_Table(self, **kwargs) -> pd.DataFrame:
        """Export filter to a SimpleTable object

        Parameters
        ----------
        fname: str
            filename

        Uses `SimpleTable` parameters
        """
        data = pd.DataFrame(
            {
                "WAVELENGTH": self._wavelength,
                "THROUGHPUT": self.transmit,
            }
        )

        if self.wavelength_unit is not None:
            data.attrs["WAVELENGTH_UNIT"] = self.wavelength_unit
            data.attrs["units"] = {"WAVELENGTH": self.wavelength_unit}
        data.attrs["DETECTOR"] = self.dtype
        data.attrs["COMPNAME"] = str(self.name)
        data.attrs["NAME"] = str(self.name)
        data.attrs["COMMENT"] = {
            "THROUGHPUT": "filter throughput definition",
            "WAVELENGTH": "filter wavelength definition",
            "WAVELENGTH_UNIT": self.wavelength_unit or "AA",
        }
        return data

    def to_dict(self) -> dict:
        """Return a dictionary of the filter"""
        data: Dict[str, Any] = {
            "WAVELENGTH": self._wavelength,
            "THROUGHPUT": self.transmit,
        }
        if self.wavelength_unit is not None:
            data["WAVELENGTH_UNIT"] = self.wavelength_unit
        data["DETECTOR"] = self.dtype
        data["NAME"] = self.name
        data["PIVOT"] = self._lpivot
        data["CENTRAL"] = self._cl
        data["EFFECTIVE"] = _drop_units(self.leff)
        data["NORM"] = self.norm
        return data

    def write_to(self, fname: str, **kwargs):
        """Export filter to a file

        Parameters
        ----------
        fname: str
            filename

        Uses `SimpleTable.write` parameters
        """
        data = self.to_Table()
        to_file(data, fname, **kwargs)
