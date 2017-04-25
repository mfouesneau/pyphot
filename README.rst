pyphot -- A tool for computing photometry from spectra
======================================================

This is a set of tools to compute synthetic photometry in a simple way, ideal to
integrate in larger projects.

full documentation at: http://mfouesneau.github.io/docs/pyphot/

The inputs are photonic or energetic response functions for the desired
photometric bands and stellar spectra. The modules are flexible to handle units 
in the wavelength definition through a simplified version of `pint` (link)

Filters are represented individually by a `Filter` object. Collections of
filters are handled with a `Library`. We provide an internal library that
contains a signitificant amount of common filters.

Each filter is minimally defined by a `wavelength` and `throughput`. Many
properties such as central of pivot wavelength are computed internally. When
units are provided for the wavelength, zero points in multiple units are also
accessible (AB, Vega magnitude, Jy, erg/s/cm2/AA). The default detector type is
assumed to be photonic, but energetic detectors are also handled for the
computations.

Installation
------------

Install with `pip`
`pip install git+git://github.com/mfouesneau/pyphot`

(`--user` if you want to install it in your user profile)

Manual installation

download the repository and run the setup
`python setup.py install`



Contributors
------------

Author:

Morgan Fouesneau

Direct contributions to the code base:

* Tim Morton
* Ariane Lancon (@lancon)
