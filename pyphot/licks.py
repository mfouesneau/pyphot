""" Lick indices calculations


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

.. todo::

    * fix units: all must be internally computed in AA, flux are not check for per AA

"""
from __future__ import print_function, division
import numpy as np

from .ezunits import unit, hasUnit
from .config import libsdir
try:
    from scipy.integrate import trapezoid
except ImportError:  # older scipy / numpy < 2.0
    from scipy.integrate import trapz as trapezoid

__default__ = libsdir.joinpath("licks.dat")


def _drop_units(q):
    """ Drop the unit definition silently """
    try:
        return q.magnitude
    except AttributeError:
        return q


def reduce_resolution(wi, fi, fwhm0=0.55, sigma_floor=0.2):
    """ Adapt the resolution of the spectra to match the lick definitions

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
    w_lick_res = (4000., 4400., 4900., 5400., 6000.)
    lick_res = (11.5, 9.2, 8.4, 8.4, 9.8)   # FWHM in AA

    w = np.asarray(wi)
    flux = np.atleast_2d(fi)

    # Linear interpolation of lick_res over w
    # numpy interp does constant instead of extrapolation
    # res = np.interp(w, w_lick_res, lick_res)

    # spline order: 1 linear, 2 quadratic, 3 cubic ...
    from scipy.interpolate import InterpolatedUnivariateSpline
    res = InterpolatedUnivariateSpline(w_lick_res, lick_res, k=1)(w)

    # Compute width from fwhm
    const = 2. * np.sqrt(2. * np.log(2))  # conversion fwhm --> sigma
    lick_sigma = np.sqrt((res ** 2 - fwhm0 ** 2)) / const

    # Convolution by g=1/sqrt(2*pi*sigma^2) * exp(-r^2/(2*sigma^2))
    flux_red = np.zeros(flux.shape, dtype=flux.dtype)

    for i, sigma in enumerate(lick_sigma):
        maxsigma = 3. * sigma
        # sampling floor: min (0.2, sigma * 0.1)
        delta = min(sigma_floor, sigma * 0.1)
        delta_wj = np.arange(-maxsigma, + maxsigma, delta)
        wj = delta_wj + w[i]
        for k, fk in enumerate(flux):
            fluxj = np.interp(wj, w, fk, left=0., right=0.)
            flux_red[k, i] = np.sum(fluxj * delta * np.exp(-0.5 * (delta_wj / sigma) ** 2))

    flux_red /= lick_sigma * const

    return flux_red.reshape(np.shape(fi))


class LickIndex(object):
    """ Define a Lick Index similarily to a Filter object """

    def __init__(self, name, lick, unit='AA'):
        """ Constructor

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
        """ return a dictionary of the current index """
        d = {}
        d.update(**self._lick)
        return d

    def _get_wavelength_attrs_with_units(self, attrname, units='AA'):
        """ return the unitwise definition corresponding to attrname """
        attr = self._lick[attrname]
        if self.wavelength_unit is not None:
            if units is None:
                return attr * unit[self.wavelength_unit]
            else:
                return (attr * unit[self.wavelength_unit]).to(units)
        else:
            return attr

    @property
    def band(self):
        """ Unitwise band definition """
        return self._get_wavelength_attrs_with_units('band')

    @property
    def blue(self):
        """ Unitwise band definition """
        return self._get_wavelength_attrs_with_units('blue')

    @property
    def red(self):
        """ Unitwise band definition """
        return self._get_wavelength_attrs_with_units('red')

    @property
    def index_unit(self):
        return self._lick['unit']

    def __repr__(self):
        return """LickIndex ({0}), {1}""".format(self.name, object.__repr__(self))

    def info(self):
        """ display information about the current Index"""
        txt = """Lick Index {s.name}
    wavelength units:     {s.wavelength_unit}
    Index Band:           {s.band}
    Blue continuum band:  {s.blue}
    Red continuum band:   {s.red}
    Measurement unit:     {s.index_unit}""".format(s=self)
        print(txt)

    def __call__(self, *args, **kwargs):
        """ compute spectral index after continuum subtraction

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
        if hasUnit(wave):
            _w = wave.to('AA').magnitude
        else:
            print("Warning: assuming units are in Angstroms")
            _w = _drop_units(wave)
        _f = _drop_units(flux)

        blue = self._get_wavelength_attrs_with_units('blue').magnitude
        red = self._get_wavelength_attrs_with_units('red').magnitude
        band = self._get_wavelength_attrs_with_units('band').magnitude
        
        nocheck = kwargs.pop('nocheck', False)
        not_covered = (blue[0] < _w[0]) | (red[-1] > _w[-1])
        if (not_covered):
            if (not nocheck):
                raise ValueError("Spectrum does not cover this index.")
            else:
                return np.zeros(_f.shape[0]) * float('nan') 
        else:
            return self._get_indice(_w, _f, blue, red, band, self.index_unit, **kwargs)

    @classmethod
    def _get_indice(cls, w, flux, blue, red, band=None, unit='ew', degree=1,
                    **kwargs):
        """ compute spectral index after continuum subtraction

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
        wi, fi = cls.continuum_normalized_region_around_line(w, flux, blue,
                                                             red, band=band,
                                                             degree=degree)
        if unit in (0, 'ew', 'EW'):
            return trapezoid(1. - fi, wi, axis=-1)
        else:
            m = trapezoid(fi, wi, axis=-1)
            m = -2.5 * np.log10(m / np.ptp(wi))
            return m

    @classmethod
    def continuum_normalized_region_around_line(cls, wi, fi, blue, red, band=None,
                                                degree=1):
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

        example
        ~~~~~~~

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
        indcont = (((w >= blue[0]) & (w <= blue[1])) |
                   ((w >= red[0]) & (w <= red[1]))
                   )
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
            # divide the flux by the polynomial and put the result in our new flux
            # array
            f[i, :] = flux[i, indrange] / np.polyval(linecoeff, wnew)
        return wnew, np.squeeze(f)


class LickLibrary(object):
    """ Collection of Lick indices """
    def __init__(self, fname=__default__, comment='#'):
        self.source = fname
        data, hdr = self._read_lick_list(fname, comment)
        self._content = data
        self._hdr = hdr

    @property
    def description(self):
        """ any comment in the input file """
        return self._hdr

    @classmethod
    def _read_lick_list(cls, fname=__default__, comment='#'):
        """ read the list of lick indices

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
        with open(fname, 'r') as f:
            data = {}
            hdr = []
            for line in f:
                if line[0] != comment:
                    elements = line.split()
                    attr = dict(
                        band=(float(elements[1]), float(elements[2])),
                        blue=(float(elements[3]), float(elements[4])),
                        red=(float(elements[5]), float(elements[6])),
                        unit='mag' if int(elements[7]) > 0 else 'ew',
                    )
                    name = elements[8]
                    data[name] = attr
                else:
                    hdr.append(line[1:-1])
        return data, hdr

    def __repr__(self):
        return "Lick Index Library: {0}\n{1:s}".format(self.source, object.__repr__(self))

    def __enter__(self):
        """ Enter context """
        return self

    def __exit__(self,  *exc_info):
        """ end context """
        return False

    def __len__(self):
        """ Size of the library """
        return len(self.content)

    def get_library_content(self):
        return list(self._content.keys())

    def __getitem__(self, name):
        """ Make this object like a dictionary and load one or multiple filters """
        with self as s:
            try:
                f = s._load_filter(name)
            except TypeError:
                f = [s._load_filter(k) for k in name]
        return f

    def _load_filter(self, fname, **kwargs):
        """ Load a given filter from the library """
        with self as s:
            return LickIndex(fname, s._content[fname])

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
