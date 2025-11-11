pyphot -- A tool for computing photometry from spectra
======================================================

.. image:: https://img.shields.io/pypi/v/pyphot.svg
    :target: https://pypi.org/project/pyphot/

.. image:: https://zenodo.org/badge/70060728.svg
   :target: https://zenodo.org/badge/latestdoi/70060728

.. image:: https://static.pepy.tech/badge/pyphot
   :target: https://pepy.tech/project/pyphot

.. image:: https://static.pepy.tech/badge/pyphot/month
   :target: https://pepy.tech/project/pyphot

.. image:: https://img.shields.io/badge/python-3.9,_3.10,_3.11,_3.12,_3.13-blue.svg

This is a set of tools to compute synthetic photometry in a simple way, ideal to
integrate in larger projects.

full documentation at: http://mfouesneau.github.io/pyphot/

The inputs are photonic or energetic response functions for the desired
photometric bands and stellar spectra. The modules are flexible to handle units
in the wavelength definition through a simplified version of `pint` (link) and `astropy.units`

Filters are represented individually by a `Filter` object. Collections of
filters are handled with a `Library`. We provide an internal library that
contains a signitificant amount of common filters.

Each filter is minimally defined by a `wavelength` and `throughput`. Many
properties such as central of pivot wavelength are computed internally. When
units are provided for the wavelength, zero points in multiple units are also
accessible (AB, Vega magnitude, Jy, erg/s/cm2/AA). The default detector type is
assumed to be photonic, but energetic detectors are also handled for the
computations.

.. image:: https://mybinder.org/badge.svg
  :target: https://mybinder.org/v2/gh/mfouesneau/pyphot/master?filepath=examples%2FQuickStart.ipynb

.. image:: https://img.shields.io/badge/render%20on-nbviewer-orange.svg
  :target: https://nbviewer.jupyter.org/github/mfouesneau/pyphot/tree/master/examples/



.. note::
    `cphot <https://github.com/mfouesneau/cphot>`_ is a spin-off project that provides a C++ version of pyphot

What's new?
-----------
* [November 11, 2025] major refactoring, remove code duplication, added type hints (v2.0.0) -- See `release notes <https://mfouesneau.github.io/pyphot/whats_new.html>`_
* [July 2, 2025] Added Vega flavors to the Vega class, allowing users to select different Vega reference spectra for photometric calculations. (v1.7.0)
* [January 21, 2025] Moved installation to use `pyproject.toml`; Updated documentation; Updates for deprecated calls. (v1.6.0)
* [January 21, 2025] Support for Python 3.13 (v1.5.0);
* [June 26, 2024] Dropping support for Python <= 3.8 (due to HDF5 modules). Minor updates for Scipy 1.14.0 and Numpy 2.0
* [November 22, 2021] new filters, SVO interface, automated tests and documentation.
* [November 6, 2019] astropy version available in beta (`from pyphot import astropy as pyphot`)
* [April 29, 2019] sandbox contains fully unit aware passbands and lick indices libraries
* [April 15, 2019] merged UncertainFilter to main, sandbox contains passbands accounting for spectrum units
* [March 4, 2019] added flux calculations in photon/s/cm2
* [March 4, 2019] added many properties per filter (lphot, lmin, lmax)
* [June 12, 2018] adding Sun reference spectra (see `:class:Sun`)
* [Apr. 26, 2018] includes Gaia nominal, DR2 and revised DR2 passbands

Installation
------------
Before installation, make sure you have HDF5 version 1.8.4 or above (this is required for pytables, see details at: https://github.com/PyTables/PyTables). We will remove this dependency in a future release.

For OSX:

.. code::

  brew install hdf5

For Debian-based distributions:

.. code::

  sudo apt-get install libhdf5-serial-dev



* Installation from PyPI

.. code::

  pip install pyphot

  pip install git+https://github.com/mfouesneau/pyphot   # if you want the unreleased version

* Manual installation

Download the repository and run the setup

.. code::

  pip install .



Contributors
------------

Author:

Morgan Fouesneau

Direct contributions to the code base:

* Tim Morton (@timothydmorton)
* Ariane Lancon (@lancon)
