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

.. image:: https://img.shields.io/badge/python-3.9,_3.10,_3.11,_3.12,_3.13,_3.14-blue.svg

This is a set of tools for computing synthetic photometry in a simple way, ideal for integration into larger projects.

Full documentation at: http://mfouesneau.github.io/pyphot/

The inputs are photonic or energetic response functions for the desired photometric bands and stellar spectra. The modules are flexible to handle units in the wavelength definition through a simplified version of `pint <https://pint.readthedocs.io/en/stable/>`_  and `astropy.units <https://docs.astropy.org/en/stable/units/index.html>`_

Filters are represented individually by a `Filter` object. Collections of filters are handled with a `Library`. We provide an internal library that contains a significant number of common filters.

Each filter is minimally defined by a `wavelength` and `throughput`. Many properties, such as the central or pivot wavelength, are computed internally. When units are provided for the wavelength, zero points in multiple units are also accessible (AB, Vega magnitude, Jy, erg/s/cm2/AA). The default detector type is assumed to be photonic, but energetic detectors are also handled for the computations.

.. image:: https://mybinder.org/badge.svg
  :target: https://mybinder.org/v2/gh/mfouesneau/pyphot/master?filepath=examples%2FQuickStart.ipynb

.. image:: https://img.shields.io/badge/render%20on-nbviewer-orange.svg
  :target: https://nbviewer.jupyter.org/github/mfouesneau/pyphot/tree/master/examples/

What's new?
-----------
See `release notes <https://mfouesneau.github.io/pyphot/whats_new.html>`_

Installation
------------
Before installation, make sure you have HDF5 version 1.8.4 or later (required for pytables; see details at https://github.com/PyTables/PyTables). We will remove this dependency in a future release.

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

  git clone https://github.com/mfouesneau/pyphot
  cd pyphot
  pip install .



Contributors
------------

Author:

Morgan Fouesneau

Direct contributions to the code base:

* Tim Morton (@timothydmorton)
* Ariane Lancon (@lancon)

Related projects
----------------

- `cphot <https://github.com/mfouesneau/cphot>`_ is a spin-off project that provides a C++ version of pyphot
