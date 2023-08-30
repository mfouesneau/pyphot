.. pyphot documentation master file, created by
   sphinx-quickstart on Wed Oct  5 11:25:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PYPHOT -- A tool for computing photometry from spectra
======================================================

This is a set of tools to compute synthetic photometry in a simple way, ideal to
integrate in larger projects.

The inputs are photonic or energetic response functions for the desired
photometric bands and stellar spectra. The modules are flexible to handle units 
in the wavelength definition through a simplified version of `pint` (link)

Filters are represented individually by a `Filter` object. Collections of
filters are handled with a `Library`. We provide an internal library that
contains a signitificant amount of common filters.

Each filter is minimally defined by a `wavelength` and `throughput`. Many
properties such as central of pivot wavelength are computed internally. 

**When units** are provided for the wavelength of the filters, zero points in
multiple units are also accessible (AB, Vega magnitude, Jy, erg/s/cm2/AA). The
default detector type is assumed to be photonic, but energetic detectors are
also handled for the computations.

.. tip::

        All provided filters have defined units and detector type
        that are used transparently throughout the package.

.. image:: https://mybinder.org/badge.svg 
  :target: https://mybinder.org/v2/gh/mfouesneau/pyphot/master?filepath=examples%2FQuickStart.ipynb
 
.. image:: https://img.shields.io/badge/render%20on-nbviewer-orange.svg
  :target: https://nbviewer.jupyter.org/github/mfouesneau/pyphot/tree/master/examples/
  

Package main content
~~~~~~~~~~~~~~~~~~~~

This package is mostly organized around 2 main classes:

* :class:`pyphot.phot.Filter` that handles filter definitions and manipulations.
 
* :class:`pyphot.phot.Library` that handles a collection of filters in a few formats. 

Both classes are able to manipulate units through a lightweight version of
**Pint** :class:`pyphot.ezunits` (...link....)

Additionally, and for convenience 

* the library class is derived into an ascii file reader
  :class:`pyphot.phot.Ascii_Library`, and a single hdf file reader
  :class:`pyphot.phot.HDF_Library`. 

* :class:`pyphot.vega.Vega` provides an interface to a synthetic reference of
  Vega.

* :class:`pyphot.sun.Sun` provides an interface to a synthetic and empirical reference of
  the Sun spectrum.

* :class:`pyphot.licks` provides an extension to computing lick indices.

Contents
---------

.. toctree::
   :maxdepth: 2

   Home <index>
   more_examples
   photometry
   libcontent
   licks
   modules


Installation
~~~~~~~~~~~~

* Using pip: Use the `--user` option if you don't have permissions to install libraries

.. code-block:: none

        pip install git+https://github.com/mfouesneau/pyphot

* Manually:

.. code-block:: none

        git clone https://github.com/mfouesneau/pyphot
        cd pyphot
        python setup.py intall



Quick Start
~~~~~~~~~~~

.. code-block:: python

        import pyphot
        # get the internal default library of passbands filters
        lib = pyphot.get_library()
        print("Library contains: ", len(lib), " filters")
        # find all filter names that relates to IRAC 
        # and print some info
        f = lib.find('irac')
        for name in f:
            lib[name].info(show_zeropoints=True)

.. code-block:: none

        Library contains:  196  filters
        Filter object information:
            name:                 SPITZER_IRAC_45
            detector type:        photon
            wavelength units:     AA
            central wavelength:   45110.141614 angstrom
            pivot wavelength:     45020.219955 angstrom
            effective wavelength: 44425.747085 angstrom
            norm:                 4664.680820
            definition contains 417 points

            Zeropoints
                Vega: 28.933084 mag,
                      2.671569250882836e-12 erg / angstrom * centimeter ** 2 * second,
                      175.8794962126167 Jy
                  AB: 25.674986 mag,
                      5.370385702161592e-11 erg / angstrom * centimeter ** 2 * second,
                      3535.5277855945205 Jy
                  ST: 21.100000 mag,
                      3.6307805477010028e-09 erg / angstrom * centimeter ** 2 * second,
                      239027.9995089771 Jy 

        [...]

Suppose one has a calibrated spectrum and wants to compute the vega magnitude
throug the HST WFC3 F110W passband, 

.. code-block:: python

        # convert to magnitudes
        import numpy as np
        f = lib['hst_wfc3_f110w']
        # compute the integrated flux through the filter f
        # note that it work on many spectra at once
        fluxes = f.get_flux(lamb, spectra, axis=1)
        # convert to vega magnitudes
        mags = -2.5 * np.log10(fluxes) - f.Vega_zero_mag
        # or similarly
        mags = -2.5 * np.log10(fluxes / f.Vega_zero_flux)


If one wants to use a given transmission curve as filter, defined by `lamb_T` and
`T`, one would use the :class:`pyphot.phot.Filter` directly as 

.. code-block:: python

        # convert to magnitudes
        from pyphot import Filter
        # if lamb_T has units the Filter object will use those.
        f = Filter(lamb_T, T, name='my_filter', dtype='photon', unit='Angstrom')
        # compute the integrated flux through the filter f
        fluxes = f.get_flux(lamb, spectra, axis=1)
        ...

       
Internal Vega reference
~~~~~~~~~~~~~~~~~~~~~~~

As mentioned in the above, sometimes a spectrum of reference of Vega is
necessary. 

We use the synthetic spectrum provided by Bohlin 2007, a common reference
througout many photometric suites.

The interface to the Vega template is given through the :class:`pyphot.vega.Vega` class.


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

