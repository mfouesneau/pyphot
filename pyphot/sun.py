"""Handle the Sun Spectrum"""

from __future__ import print_function

import warnings
from typing import Optional, Literal, Tuple, cast
import numpy.typing as npt
import pandas as pd

from . import config
from .unit_adapters import QuantityType
from . import io

from .config import libsdir

__all__ = ["Sun"]

_default_sun = {
    "observed": "{0}/sun_reference_stis_001.fits".format(libsdir),
    "theoretical": "{0}/sun_kurucz93.fits".format(libsdir),
}


def get_library_default_distance() -> QuantityType:
    """Return the default distance associated with the library spectra of the Sun."""
    return 1.0 * config.units.U("au")


class Sun:
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

    +-----------+------------+----------+-------------+
    | log(Z)    | Teff       | log_g    | V_{Johnson} |
    +===========+============+==========+=============+
    | +0.0      | +5777      | +4.44    | -26.75      |
    +-----------+------------+----------+-------------+

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

    _data: Optional[pd.DataFrame]
    units: Optional[Tuple[str, str]]
    distance: QuantityType
    distance_conversion: float
    flavor: Optional[Literal["observed", "theoretical"]]

    def __init__(
        self,
        *,
        source: Optional[str] = None,
        distance: Optional[QuantityType] = None,
        flavor: Literal["observed", "theoretical"] = "theoretical",
    ):
        """Constructor"""
        self._data = None
        self.units = None
        if distance is None:
            distance = get_library_default_distance()
        self.distance = distance
        self.distance_conversion = 1.0
        self._set_source_flavor(source=source, flavor=flavor)

    def _set_source_flavor(
        self,
        *,
        source: Optional[str],
        flavor: Optional[Literal["observed", "theoretical"]],
    ):
        """Set the source and flavor of the Sun spectrum"""
        if source is not None:
            # if source provided
            self.source = source
            self.flavor = flavor
        else:
            if flavor is None:
                raise RuntimeError("Either `source` or `flavor` must be provided.")

            if flavor.lower() not in _default_sun:
                raise RuntimeError(
                    f"Flavor must be in {_default_sun.keys()}, got {flavor}."
                )
            _default_distance = get_library_default_distance()

            self.flavor = flavor
            self.source = _default_sun[flavor]
            self.distance_conversion = float(((_default_distance / self.distance) ** 2))
        self._data = None
        self.units = None

    def _readfile(self, fname: Optional[str] = None) -> Tuple[pd.DataFrame, str, str]:
        """Read the data file and populate the data and units attributes"""
        if (self._data is not None) and (self.units is not None):
            return self._data, self.units[0], self.units[1]

        fname = fname or self.source

        df, hdr = io.from_file(fname)
        self._data = df

        try:
            uw = self._data.attrs["WAVELENGTH_UNIT"].split("=")[0].rstrip()
            uf = self._data.attrs["FLUX_UNIT"].split("=")[0].rstrip()
            self.units = uw, uf
        except TypeError:
            uw = self._data.attrs["WAVELENGTH_UNIT"].split(b"=")[0].decode().rstrip()
            uf = self._data.attrs["FLUX_UNIT"].split(b"=")[0].decode().rstrip()
            self.units = uw, uf

            return self._data, self.units[0], self.units[1]

        return self._data, self.units[0], self.units[1]

    def __enter__(self):
        """Enter context"""
        self._readfile()
        return self

    def __exit__(self, *exc_info):
        """end context"""
        return False

    @property
    def data(self) -> pd.DataFrame:
        if self._data is not None:
            return self._data
        else:
            data, _, _ = self._readfile()
            return data

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
        flux = data["FLUX"].to_numpy() * self.distance_conversion
        return flux * config.units.U(f_units.lower())
