"""Lick indices calculations


This package provides function to compute spectral indices


A collection of many common indices is available in `licks.dat`

The Lick system of spectral line indices is one of the most commonly used
methods of determining ages and metallicities of unresolved (integrated light)
stellar populations.

The calibration of the Lick/ IDS system is complicated because the original
Lick spectra were not flux calibrated, so there are usually systematic effects
due to differences in continuum shape.  Proper calibration involves observing
many of the original Lick/IDS standard stars and deriving offsets to the
standard system.

references
~~~~~~~~~~

    Worthey G., Faber S. M., Gonzalez J. J., Burstein D., 1994, ApJS, 94, 687
    Worthey G., Ottaviani D. L., 1997, ApJS, 111, 377
    Puzia et al. 2002
    Zhang, Li & Han 2005, http://arxiv.org/abs/astro-ph/0508634v1


notes
~~~~~

    In Vazdekis et al. (2010), we propose a new Line Index System, hereafter
    LIS, with three new spectral resolutions at which to measure the Lick
    indices. Note that this new system should not be restricted to the Lick set
    of indices in a flux calibrated system. In fact, LIS can be used for any
    index in the literature (e.g., for the Rose (1984) indices), including
    newly defined indices (e.g., Cervantes & Vazdekis 2009).


    The LIS system is defined for 3 different spectral resolutions which are
    best suited for the following astrophysical cases:

    LIS-5.0AA: globular clusters
    LIS-8.4AA: low and intermediate-mass galaxies
    LIS-14.0AA: massive galaxies
    Conversions to transform the data from the Lick/IDS system to LIS can be
    found


    discussion of indices and information
    Johansson, Thomas & Maraston 2010
    http://wwwmpa.mpa-garching.mpg.de/~jonasj/milesff/milesff.pdf

"""

from dataclasses import dataclass
from typing import Union, Literal, Optional, Tuple, Any, Dict, List
import numpy as np
import numpy.typing as npt

from scipy.integrate import trapezoid
from scipy.interpolate import InterpolatedUnivariateSpline

from .unit_adapters import QuantityType, enforce_default_units
from . import config


def _reduce_resolution(
    wi_AA: npt.NDArray[np.floating],
    fi_flam: npt.NDArray[np.floating],
    fwhm0_AA: float = 0.55,
    sigma_floor_AA: float = 0.2,
) -> npt.NDArray[np.floating]:
    """Adapt the resolution of the spectra to match the lick definitions

    Lick definitions have different resolution elements as function
    of wavelength. These definition are hard-coded in this function

    Parameters
    ----------
    wi_AA: ndarray (n, )
        wavelength definition in Angstroms
    fi_flam: ndarray (nspec, n) or (n, )
        spectra to convert in flam
    fwhm0_AA: float
        initial broadening in the spectra `fi` in Angstroms
    sigma_floor_AA: float
        minimal dispersion to consider in Angstroms

    Returns
    -------
    flux_red: ndarray (nspec, n) or (n, )
        reduced spectra in flam
    """

    # all in AA
    w_lick_res = (4000.0, 4400.0, 4900.0, 5400.0, 6000.0)
    lick_res = (11.5, 9.2, 8.4, 8.4, 9.8)  # FWHM in AA

    w = np.asarray(wi_AA)
    flux = np.atleast_2d(fi_flam)

    # Linear interpolation of lick_res over w
    # numpy interp does constant instead of extrapolation
    # res = np.interp(w, w_lick_res, lick_res)

    # spline order: 1 linear, 2 quadratic, 3 cubic ...
    res = InterpolatedUnivariateSpline(w_lick_res, lick_res, k=1)(w)
    res = np.asarray(res)

    # Compute width from fwhm
    const = 2.0 * np.sqrt(2.0 * np.log(2))  # conversion fwhm --> sigma
    lick_sigma = np.sqrt((res**2 - fwhm0_AA**2)) / const

    # Convolution by g=1/sqrt(2*pi*sigma^2) * exp(-r^2/(2*sigma^2))
    flux_red = np.zeros(flux.shape, dtype=flux.dtype)

    for i, sigma in enumerate(lick_sigma):
        maxsigma = 3.0 * sigma
        # sampling floor
        delta = min(sigma_floor_AA, sigma * 0.1)
        delta_wj = np.arange(-maxsigma, +maxsigma, delta)
        wj = delta_wj + w[i]
        for k, fk in enumerate(flux):
            fluxj = np.interp(wj, w, fk, left=0.0, right=0.0)
            flux_red[k, i] = np.sum(
                fluxj * delta * np.exp(-0.5 * (delta_wj / sigma) ** 2)
            )

    flux_red /= lick_sigma * const

    return flux_red.reshape(np.shape(fi_flam))


def reduce_resolution(
    wi: Union[npt.NDArray[np.floating], QuantityType],
    fi: Union[npt.NDArray[np.floating], QuantityType],
    fwhm0: Union[float, QuantityType] = 0.55,
    sigma_floor: Union[float, QuantityType] = 0.2,
    warn: bool = True,
) -> QuantityType:
    """Adapt the resolution of the spectra to match the lick definitions

    Lick definitions have different resolution elements as function
    of wavelength. These definition are hard-coded in this function

    Parameters
    ----------
    wi: ndarray (n, ) with or without units
        wavelength definition
    fi: ndarray (nspec, n) or (n, ) with or without units
        spectra to convert
    fwhm0: float or QuantityType
        initial broadening in the spectra `fi`
    sigma_floor: float or QuantityType
        minimal dispersion to consider
    warn: bool
        Warn for unit assumptions

    Returns
    -------
    flux_red: ndarray (nspec, n) or (n, )
        reduced spectra
    """
    wi_AA = config.units.val_in_unit("wi", wi, "AA", warn=warn)
    fi_flam = config.units.val_in_unit("fi", fi, "flam", warn=warn)
    fwhm0_AA = config.units.val_in_unit("fwhm0", fwhm0, "AA", warn=warn)
    sigma_floor_AA = config.units.val_in_unit("sigma", sigma_floor, "AA", warn=warn)
    flux_red = _reduce_resolution(
        wi_AA.value,
        fi_flam.value,
        fwhm0_AA.value,
        sigma_floor_AA.value,
    )
    return flux_red * config.units.U("flam")


def _split_value_unit(q: Union[QuantityType, Any]) -> tuple[Any, Union[str, None]]:
    """Split a quantity into value and unit"""
    try:
        return q.value, str(q.unit)
    except AttributeError:
        return q, None


@dataclass
class LickDefinition:
    """Definition of a Lick Index

    expecting 'blue', 'red', 'band', and 'unit' definitions
    - `blue` and `red` are used to continuum normalize the spectra
    - `band` covers the index itself.
    - `index_unit` gives the index measurement type as either magnitude (mag) or equivalent width (ew)

    Consistency is checked such that all fields are set with consistent units (priority to band field)
    wavelength_unit is also associated
    """

    blue: QuantityType
    red: QuantityType
    band: QuantityType
    index_unit: Literal["mag", "ew"] = "mag"
    wavelength_unit: Optional[str] = None

    @staticmethod
    def _ensure_consistency(
        *,
        band: Union[QuantityType, Tuple[float, float]],
        blue: Union[QuantityType, Tuple[float, float]],
        red: Union[QuantityType, Tuple[float, float]],
        index_unit: Literal["mag", "ew"] = "mag",
        wavelength_unit: Optional[str] = None,
        _skip: bool = False,
    ) -> dict:
        """Check the definition is consistent

        All fields are set with consistent units (priority to band field)
        wavelength_unit is also associated

        """
        if _skip:
            return {}

        u = config.units

        if u.has_unit(band):
            _, band_unit = _split_value_unit(band)
        elif u.has_unit(blue):
            _, band_unit = _split_value_unit(blue)
        elif u.has_unit(red):
            _, band_unit = _split_value_unit(red)
        elif wavelength_unit:
            band_unit = wavelength_unit
        else:
            raise ValueError(
                "Could not determine the units of wavelength from the index definition"
                "blue/red/band need units or wavelength_unit must be provided"
            )

        if index_unit not in ("mag", "ew"):
            raise ValueError(
                f"Invalid index unit: {values['index_unit']}. Expected 'mag' or 'ew'"
            )
        return dict(
            band=u.val_in_unit("band", band, band_unit, warn=False),
            blue=u.val_in_unit("blue", blue, band_unit, warn=False),
            red=u.val_in_unit("red", red, band_unit, warn=False),
            wavelength_unit=band_unit,
            index_unit=index_unit,
        )

    def __init__(
        self,
        *,
        band: Union[QuantityType, Tuple[float, float]],
        blue: Union[QuantityType, Tuple[float, float]],
        red: Union[QuantityType, Tuple[float, float]],
        index_unit: Literal["mag", "ew"] = "mag",
        wavelength_unit: Optional[str] = None,
    ):
        """Build a LickDefinition from values

        Parameters
        ----------
        band : QuantityType or tuple of float
            Bandpass of the lick
        blue : QuantityType or tuple of float
            Blue limit of the lick
        red : QuantityType or tuple of float
            Red limit of the lick
        index_unit : Literal["mag", "ew"], optional
            Unit of the index, by default "mag"
        wavelength_unit : Optional[str], optional
            Unit of the wavelength, by default None

        Returns
        -------
        LickDefinition
            The LickDefinition object
        """
        sure = self._ensure_consistency(
            band=band,
            blue=blue,
            red=red,
            wavelength_unit=wavelength_unit,
            index_unit=index_unit,
        )
        self.__dict__.update(sure)


class LickIndex(object):
    """Define a Lick Index similarily to a Filter object"""

    _lick: LickDefinition
    name: str

    def __init__(
        self,
        name: str,
        lick: Union[dict, LickDefinition],
        unit="AA",
    ):
        """Constructor

        Parameters
        ----------
        name: str
            name of the index
        lick: dict or LickDefinition
            expecting 'blue', 'red', 'band', and 'unit' definitions
            `blue` and `red` are used to continuum normalize the spectra
            `band` covers the index itself. `unit` gives the index measurement
            units, either magnitudes (mag) or equivalent width (ew)
        unit: str
            wavelength unit of the intervals if not provided in the lick definition
        """
        self.name = name
        # lick integrity definition done by the LickDefinition class
        if isinstance(lick, dict):
            # make sure wavelength_unit is set
            lick.setdefault("wavelength_unit", unit)
            self._lick = LickDefinition(**lick)
        elif isinstance(lick, LickDefinition):
            # make sure wavelength_unit is set
            self._lick = lick
        else:
            raise TypeError("lick must be a dictionary or LickDefinition")

    @property
    def wavelength_unit(self) -> str:
        """Wavelength unit of the intervals"""
        if self._lick.wavelength_unit is None:
            raise ValueError("wavelength_unit is not defined")
        return self._lick.wavelength_unit

    def to_dict(self) -> dict:
        """return a dictionary of the current index"""
        d = {
            "blue": self._lick.blue,
            "red": self._lick.red,
            "band": self._lick.band,
            "unit": self._lick.index_unit,
            "wavelength_unit": self._lick.wavelength_unit,
        }

        return d

    def _get_wavelength_attrs_with_units(
        self,
        attrname: str,
        units: Optional[str] = None,
    ) -> QuantityType:
        """return the unitwise definition corresponding to attrname"""
        attr = getattr(self._lick, attrname)
        if units is not None:
            return attr.to(units)
        return attr

    @property
    def band(self) -> QuantityType:
        """Unitwise band definition"""
        return self._get_wavelength_attrs_with_units("band")

    @property
    def blue(self) -> QuantityType:
        """Unitwise band definition"""
        return self._get_wavelength_attrs_with_units("blue")

    @property
    def red(self) -> QuantityType:
        """Unitwise band definition"""
        return self._get_wavelength_attrs_with_units("red")

    @property
    def index_unit(self) -> Literal["mag", "ew"]:
        """Index definition type"""
        return self._lick.index_unit

    def __repr__(self) -> str:
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

    @classmethod
    def _get_indice(
        cls,
        w: npt.NDArray[np.floating],
        flux: npt.NDArray[np.floating],
        blue: tuple[float, float],
        red: tuple[float, float],
        band: tuple[float, float] | None = None,
        unit: str = "ew",
        degree: int = 1,
        **kwargs,
    ) -> npt.NDArray[np.floating]:
        """compute spectral index after continuum subtraction

        .. warning::
            This function is for internal use and not checking units.

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

    @enforce_default_units(None, "AA", "flam")
    def get(
        self, wave: QuantityType, flux: QuantityType, **kwargs
    ) -> npt.NDArray[np.floating]:
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
        ValueError: when the spectral coverage wave does not cover the index range
        """

        # let's make sure we work in AA and flam
        _w = config.units.val_in_unit("wave", wave, "AA").value
        _f = config.units.val_in_unit("flux", flux, "flam").value

        blue = self._get_wavelength_attrs_with_units("blue", "AA").value
        red = self._get_wavelength_attrs_with_units("red", "AA").value
        band = self._get_wavelength_attrs_with_units("band", "AA").value

        nocheck = kwargs.pop("nocheck", False)
        not_covered = (blue[0] < _w[0]) | (red[-1] > _w[-1])
        if not_covered:
            if not nocheck:
                raise ValueError("Spectrum does not cover this index.")
            else:
                return np.zeros(_f.shape[0]) * float("nan")
        else:
            return self._get_indice(_w, _f, blue, red, band, self.index_unit, **kwargs)


class LickLibrary:
    """Collection of Lick indices"""

    source: str
    _content: Dict[str, LickDefinition]
    _hdr: List[str]

    def __init__(self, fname: Optional[str] = None, comment: str = "#"):
        self.source: str = fname or str(config.__default_lick_lib__)
        data, hdr = self._read_lick_list(fname, comment)
        self._content = data
        self._hdr = hdr

    @property
    def description(self) -> List[str]:
        """any comment in the input file"""
        return self._hdr

    @classmethod
    def _read_lick_list(
        cls,
        fname: Optional[str] = None,
        comment="#",
    ) -> Tuple[Dict[str, LickDefinition], List[str]]:
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
        if fname is None:
            fname = str(config.__default_lick_lib__)
        with open(fname, "r") as f:
            data = {}
            hdr = []
            for line in f:
                if line[0] != comment:
                    _line = line.split()
                    attr = LickDefinition(
                        band=(float(_line[1]), float(_line[2])),
                        blue=(float(_line[3]), float(_line[4])),
                        red=(float(_line[5]), float(_line[6])),
                        index_unit="mag" if int(_line[7]) > 0 else "ew",
                        wavelength_unit="AA",
                    )
                    name = _line[8]
                    data[name] = attr
                else:
                    hdr.append(line[1:-1].strip())
        return data, hdr

    def __repr__(self) -> str:
        return "Lick Index Library: {0}\n{1:s}".format(
            self.source, object.__repr__(self)
        )

    def __enter__(self):
        """Enter context"""
        return self

    def __exit__(self, *exc_info):
        """end context"""
        return False

    def __len__(self) -> int:
        """Size of the library"""
        return len(self.content)

    def get_library_content(self) -> List[str]:
        return list(self._content.keys())

    def __getitem__(self, name) -> Union[LickIndex, List[LickIndex]]:
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
            return LickIndex(fname, current_lib._content[fname])

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
