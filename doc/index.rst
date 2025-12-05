.. pyphot documentation master file, created by
   sphinx-quickstart on Wed Oct  5 11:25:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PYPHOT -- A tool for computing photometry from spectra
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

This is a set of tools to compute synthetic photometry in a simple way, suitable
for integration in larger projects.

.. important::
   New major release (v2.0.0) which introduces new features and **breaking changes**.

   See :doc:`What's new <whats_new>`.

.. info::
   Support for python 3.14 is currently blocked by `pytables` not supporting it yet.


The inputs are response functions for the desired photometric passbands and stellar spectra.

Filters are represented individually by a `Filter` object. Collections of filters are handled with a `Library`. We provide an internal library that contains a significant number of common filters. Note that we use the word "filter" in a generic sense: in many cases, it corresponds in fact to a total throughput, which includes the actual filter transmission, the transmission of the optics, the wavelength-dependent efficiency of the detector, and (for ground-based photometric systems) the transmission of the atmosphere.

Each filter is minimally defined by a `wavelength`, a `throughput`, and a detector type, which may be either `energy` or `photon` (default). Many properties, such as central or pivot wavelengths, are computed internally.

Units are provided for the wavelengths of the filters, zero points in multiple units are also accessible (AB, Vega magnitude, Jy, erg/s/cm2/AA). The default detector is assumed to be a photon-counting device, but energy-sensitive detectors are also handled for the computations.

.. note::

    Pyphot handles units from `Astropy Units <https://docs.astropy.org/en/stable/units/index.html>`_ (as default), `Pint <https://pint.readthedocs.io/en/stable/>`_, and `EzUnits` (:mod:`~pyphot.unit_adapters.ezunits`).
    You can switch between them using the :func:`pyphot.config.set_units_backend` function.

    .. code-block:: python

        import pyphot.config
        pyphot.config.set_units_backend('astropy')
        pyphot.config.set_units_backend('pint')
        pyphot.config.set_units_backend('ezunits')



Package main content
~~~~~~~~~~~~~~~~~~~~

This package is mostly organized around 2 main classes:

* :class:`~pyphot.phot.Filter` that handles filter definitions and manipulations.

* :class:`~pyphot.libraries.Library` that handles a collection of filters in a few formats.


Additionally, and for convenience

* the library class is derived into an ascii file reader
  :class:`~pyphot.libraries.Ascii_Library`, and a single hdf file reader
  :class:`~pyphot.libraries.HDF_Library`.

* :class:`~pyphot.vega.Vega` provides an interface to a synthetic reference of Vega.

* :class:`~pyphot.sun.Sun` provides an interface to a synthetic and empirical reference of the Sun spectrum.

* :class:`~pyphot.licks` provides an extension to computing lick indices (:class:`~pyphot.licks.LickIndex`, :class:`~pyphot.licks.LickLibrary`).

* :func:`pyphot.svo.get_pyphot_filter` provides an interface to the `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_ to download filters directly from the SVO database.

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 1

   Home <index>
   What's New <whats_new>
   QuickStart
   examples
   photometry
   libcontent
   licks
   vega
   sun
   contributing
   modules


Installation
~~~~~~~~~~~~
Pyphot is available on `PyPI <https://pypi.org/project/pyphot/>`_ and can be installed using any common package manager.

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
.. image:: https://img.shields.io/pypi/v/pyphot.svg
    :target: https://pypi.org/project/pyphot/

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
~~~~~~~~~~~~

Author: Morgan Fouesneau (`@mfouesneau <https://github.com/mfouesneau>`_)

Direct contributions to the code base:

* Ariane Lançon (`@lancon <https://github.com/lancon>`_)
* Tim Morton (`@timothydmorton <https://github.com/timothydmorton>`_)

If you want to contribute, please look at our :doc:`contributing guide <contributing>`.

Related projects
~~~~~~~~~~~~~~~~

* `cphot <https://github.com/mfouesneau/cphot>`_ is a spin-off project that provides a C++ version of pyphot


Citation
~~~~~~~~

If you use this software in your work, please cite it using the following:

.. tabs::

        .. tab:: bibtex

                .. code-block::

                        @software{fouesneau2025pyphot,
                        author       = {Fouesneau, Morgan},
                        title        = {pyphot},
                        month        = jan,
                        year         = 2025,
                        publisher    = {Zenodo},
                        version      = {pyphot\_v2.0.0},
                        doi          = {10.5281/zenodo.14712174},
                        url          = {https://doi.org/10.5281/zenodo.14712174},
                        }

        .. tab:: citation cff

                `Learn more about CITATION files.  <https://citation-file-format.github.io/>`_

                .. literalinclude:: ../CITATION.cff


Indices and tables
~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
