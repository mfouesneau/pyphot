"""Handle vega spec/mags/fluxes manipulations

Works with both ascii and hd5 files for back-compatibility

Vega.wavelength and Vega.flux have now units!
"""

from typing import Optional, Tuple, cast

import numpy.typing as npt
import pandas as pd

from . import config
from .unit_adapters import QuantityType
from . import io


__all__ = ["Vega"]

_default_vega = {
    "mod_002": f"{config.libsdir}/alpha_lyr_mod_002.fits",
    "mod_003": f"{config.libsdir}/alpha_lyr_mod_003.fits",
    "mod_004": f"{config.libsdir}/alpha_lyr_mod_004.fits",
    "stis_011": f"{config.libsdir}/alpha_lyr_stis_011.fits",
    "stis_003": f"{config.libsdir}/alpha_lyr_stis_003.fits",
    "legacy": f"{config.libsdir}/vega.hd5",
}


class Vega:
    """
    Class that handles vega spectrum and references. This class know where to
    find the Vega synthetic spectrum in order to compute fluxes and magnitudes
    in given filters

    An instance can be used as a context manager as:

    >>> filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_WFC3_F475W',\
                   'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
        with Vega(flavor='stis_011') as v:
            vega_f, vega_mag, flamb = v.getSed(filters)
        print(vega_f, vega_mag, flamb)

    .. seealso::
        See the documentation :doc:`vega` for more information about the
        different flavors of Vega spectra and associated references.

    """

    _data: Optional[pd.DataFrame] = None
    """Data table read from the source file"""
    units: Optional[Tuple[str, str]] = None
    """Units of the data (wavelength, flux)"""

    def __init__(self, *, source: Optional[str] = None, flavor: str = "legacy"):
        """Constructor

        Parameters
        ----------
        source : str, optional
            Source of the Vega spectrum. If not provided, the default flavor is used.
        flavor : str, optional
            Flavor of the Vega spectrum. If not provided, the default flavor is used.
        """
        self._data: Optional[pd.DataFrame] = None
        self.units: Optional[Tuple[str, str]] = None
        self._set_source_flavor(source=source, flavor=flavor)

    def _set_source_flavor(
        self,
        *,
        source: Optional[str] = None,
        flavor: Optional[str] = None,
    ):
        """Set the source and flavor of the Vega spectrum"""
        if source is not None:
            self.source = source
        else:
            if flavor is None:
                raise RuntimeError("Either `source` or `flavor` must be provided.")

            if flavor.lower() not in _default_vega:
                raise ValueError(
                    "Unknown Vega flavor: {}. Available flavors {}".format(
                        flavor, _default_vega.keys()
                    )
                )
            self.source = _default_vega[flavor]
        self._data = None
        self.units = None

    def _readfile(self, fname: Optional[str] = None) -> Tuple[pd.DataFrame, str, str]:
        """Read the data file and populate the data and units attributes"""
        if (self._data is not None) and (self.units is not None):
            return self._data, self.units[0], self.units[1]

        fname = fname or self.source

        # handle legacy files by extension
        ext = fname.split(".")[-1]
        if ext.lower() in ("hd5", "hdf", "hdf5"):
            df, _ = io.from_file(fname, tablename="/spectrum")
        else:
            df, _ = io.from_file(fname)
        self._data = df
        self.units = df.attrs["WAVELENGTH_UNIT"], df.attrs["FLUX_UNIT"]

        try:
            uw = self._data.attrs["WAVELENGTH_UNIT"].split("=")[0].rstrip()
            uf = self._data.attrs["FLUX_UNIT"].split("=")[0].rstrip()
            self.units = uw, uf
        except TypeError:
            uw = self._data.attrs["WAVELENGTH_UNIT"].split(b"=")[0].decode().rstrip()
            uf = self._data.attrs["FLUX_UNIT"].split(b"=")[0].decode().rstrip()
            self.units = uw, uf

        return self._data, self.units[0], self.units[1]

    @property
    def data(self) -> pd.DataFrame:
        if self._data is not None:
            return self._data
        else:
            data, _, _ = self._readfile()
            return data

    def __enter__(self):
        """Enter context"""
        self._readfile()
        return self

    def __exit__(self, *exc_info):
        """end context"""
        return False

    @property
    def wavelength(self) -> QuantityType:
        """wavelength (with units when found)"""
        data, λ_units, _ = self._readfile()
        λ = cast(npt.ArrayLike, data["WAVELENGTH"].to_numpy())
        return λ * config.units.U(λ_units.lower())

    @property
    def flux(self) -> QuantityType:
        """flux(wavelength) values (with units when provided)"""
        data, _, f_units = self._readfile()
        flux = cast(npt.ArrayLike, data["FLUX"].to_numpy())
        u_ = config.units.U(f_units.lower())
        return flux * u_
