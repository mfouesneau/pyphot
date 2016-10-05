from __future__ import print_function, division
import numpy as np
from numpy import trapz

from .pbar import Pbar

# this is used to convert from bolometric luminosities to abs fluxes
# object to 10parsecs -- abs mag.
distc = 4. * np.pi * (3.0856775e19) ** 2


def progress_enumerate(it, *args, **kwargs):
    """ Enumerate over a sequence with progression if requested

    Parameter
    ---------
    show_progress: bool
        set to show progress
    """
    progress = kwargs.pop('show_progress', False)
    if progress is True:
        for a in enumerate(Pbar(**kwargs).iterover(it), *args):
            yield a
    else:
        for a in enumerate(it, *args):
            yield a


def extractPhotometry(lamb, spec, flist, absFlux=True, progress=True):
    """Extract seds from a one single spectrum

    Parameters
    ----------
    lamb: ndarray[float,ndim=1]
        wavelength of spec

    spec: ndarray[float, ndim=1]
        spectrum

    flist: list[filter]
        list of filter objects

    absflux: bool
        return SEDs in absolute fluxes if set

    progress: bool
        show progression if set

    Returns
    -------
    cls: ndarray[float, ndim=1]
        filters central wavelength

    seds: ndarray[float, ndim=1]
        integrated sed
    """
    cls  = np.empty( len(flist), dtype=float)
    seds = np.empty( len(flist), dtype=float)
    for e, k in progress_enumerate(flist, show_progress=progress, desc='Photometry'):
        xl  = k.transmit > 0.
        tmp = lamb[xl] * k.transmit[xl]
        s0  = spec[:, xl]
        # apply absolute flux conversion if requested
        if absFlux:
            s0 /= distc
        a = trapz( tmp[None, :] * s0, lamb[xl], axis=1 )
        seds[e] = a / k.lT   # divide by integral (lambda T dlambda)
        cls[e]  = k.cl

    return cls, seds


def extractSEDs(lamb, specs, flist, absFlux=True, progress=True):
    """ Extract seds from a grid

    Parameters
    ----------
    g0: ModelGrid instance
        initial spectral grid

    flist: sequence(filter)
        list of filter object instances

    absflux: bool
        return SEDs in absolute fluxes if set

    progress: bool
        show progression if set

    Returns
    -------
    cls: ndarray[float, ndim=1]
        filters central wavelength

    seds: ndarray[float, ndim=1]
        integrated sed

    grid: Table
        SED grid properties table from g0 (g0.grid)
    """
    seds = np.empty(( len(specs), len(flist) ), dtype=float)
    cls  = np.empty( len(flist), dtype=float)
    for e, k in progress_enumerate(flist, show_progress=progress, desc='Photometry'):
        xl  = k.transmit > 0.
        tmp = lamb[xl] * k.transmit[xl]
        s0  = specs[:, xl]
        # apply absolute flux conversion if requested
        if absFlux:
            s0 /= distc
        a = trapz( tmp[None, :] * s0, lamb[xl], axis=1 )
        seds[:, e] = a / k.lT
        cls[e] = k.cl

    return cls, seds


def STmag_to_flux( v ):
    """
    Convert an ST magnitude to erg/s/cm2/AA (Flambda)

    .. math::
        mag = -2.5 \log_{10}(F) - 21.10

        M0 = 21.10
        F0 = 3.6307805477010028 10^{-9} erg/s/cm2/AA

    Parameters
    ----------
    v: np.ndarray[float, ndim=N] or float
        array of magnitudes

    Returns
    -------
    flux: np.ndarray[float, ndim=N], or float
        array of fluxes
    """
    v0 = 21.1
    return 10. ** ( -0.4 * (v - v0) )


def STmag_from_flux( v ):
    """
    Convert to ST magnitude from erg/s/cm2/AA (Flambda)

    .. math::
        mag = -2.5 \log_{10}(F) - 21.10

        M0 = 21.10
        F0 = 3.6307805477010028 10^{-9} erg/s/cm2/AA

    Parameters
    ----------
    v: np.ndarray[float, ndim=N], or float
        array of fluxes

    Returns
    -------
    mag: np.ndarray[float, ndim=N], or float
        array of magnitudes
    """
    v0 = 21.1
    return -2.5 * np.log10( v ) - v0


def fluxToMag(flux):
    """ Return the magnitudes from flux values

    Parameters
    ----------
    flux: np.ndarray[float, ndim=N]
        array of fluxes

    Returns
    -------
    mag: np.ndarray[float, ndim=N]
        array of magnitudes
    """
    return -2.5 * np.log10(flux)


def fluxErrTomag(flux, fluxerr):
    """ Return the magnitudes and associated errors from fluxes and flux error
    values

    Parameters
    ----------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxerr: np.ndarray[float, ndim=1]
        array of flux errors

    Returns
    -------
    mag: np.ndarray[float, ndim=1]
        array of magnitudes

    err: np.ndarray[float, ndim=1]
        array of magnitude errors
    """
    mag = fluxToMag(flux)
    return mag, -2.5 * np.log10( 1. - fluxerr / flux )


def magToFlux(mag):
    """ Return the flux from magnitude values

    Parameters
    ----------
    mag: np.ndarray[float, ndim=N]
        array of magnitudes

    Returns
    -------
    flux:  np.ndarray[float, ndim=N]
        array of fluxes
    """
    return 10 ** (-0.4 * mag)


def magErrToFlux(mag, err):
    """ Return the flux and associated errors from magnitude and mag error values

    Parameters
    ----------
    mag: np.ndarray[float, ndim=1]
        array of magnitudes

    err: np.ndarray[float, ndim=1]
        array of magnitude errors

    Returns
    -------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxerr: np.ndarray[float, ndim=1]
        array of flux errors
    """
    flux = magToFlux(mag)
    return flux, flux * ( 1. - magToFlux(err) )
