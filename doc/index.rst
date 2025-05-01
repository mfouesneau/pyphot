.. pyphot documentation master file, created by
   sphinx-quickstart on Wed Oct  5 11:25:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PYPHOT -- A tool for computing photometry from spectra
======================================================

This is a set of tools to compute synthetic photometry in a simple way, suitable for
integration in larger projects.

The inputs are response functions for the desired photometric passbands, and stellar spectra. The modules are flexible to handle units in the wavelength definition through 
a simplified version of `pint` (link)

Filters are represented individually by a `Filter` object. Collections of
filters are handled with a `Library`. We provide an internal library that
contains a signitificant amount of common filters. Note that we use the word "filter" 
in a generic sense: in many cases it corresponds in fact to a total throughput, that includes
the actual filter transmission, the transmission of the optics, the wavelength-dependent efficiency
of the detector, and (for ground-based photometric systems) the transmission of the 
atmosphere.

Each filter is minimally defined by a `wavelength`, a `throughput`, and a detector
type that may be either `energy` or `photon` (default). Many properties such as central 
or pivot wavelengths are computed internally. 

Units are provided for the wavelengths of the filters, zero points in
multiple units are also accessible (AB, Vega magnitude, Jy, erg/s/cm2/AA). The
default detector is assumed to be a photon-counting device, but energy-sensitive
detectors are also handled for the computations.

.. note::

        All provided filters have defined units and detector type
        that are used transparently throughout the package.

.. image:: https://mybinder.org/badge.svg 
  :target: https://mybinder.org/v2/gh/mfouesneau/pyphot/master?filepath=examples%2FQuickStart.ipynb
 
.. image:: https://img.shields.io/badge/render%20on-nbviewer-orange.svg
  :target: https://nbviewer.jupyter.org/github/mfouesneau/pyphot/tree/master/examples/
  

Package main content
~~~~~~~~~~~~~~~~~~~~

This package is mainly organized around two backend to handle quantities with units: a lightweight version of
`Pint`_ (:class:`pyphot.ezunits`) and the `Astropy Units <https://docs.astropy.org/en/stable/units/index.html>`_. 
As these do no provide the same functionalities, the package is designed to be able to switch between the two by importing either from :class:`pyphot.sandbox` or :class:`pyphot.astropy`.

.. _Pint: https://pint.readthedocs.io/en/stable/

.. tip::

        You can transparently switch between pint and astropy units with aliasing the import line:

        .. code:: python

                from pyphot import astropy as pyphot
                # or 
                from pyphot import sandbox as pyphot


This package is mostly organized around 2 main classes (from either root module above):

* :class:`~pyphot.sandbox.UnitFilter` that handles filter definitions and manipulations.
 
* :class:`~pyphot.sandbox.UnitLibrary` that handles a collection of filters in a few formats. 


Additionally, and for convenience 

* the library class is derived into an ascii file reader
  :class:`~pyphot.phot.Ascii_Library`, and a single hdf file reader
  :class:`~pyphot.phot.HDF_Library`. 

* :class:`~pyphot.svo` provides an interface to the `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_ to download filters directly from the SVO database.

* :class:`~pyphot.vega.Vega` provides an interface to a synthetic reference of Vega.

* :class:`~pyphot.sun.Sun` provides an interface to a synthetic and empirical reference of the Sun spectrum.

* :class:`~pyphot.licks` provides an extension to computing lick indices.

Contents
---------

.. toctree::
   :maxdepth: 1

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

Additional examples can be found on the :doc:`more_examples` page.

If one wishes to find out details about one of the available transmission curves:

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

        Filter object information:
            name:                 SPITZER_IRAC_45
            detector type:        photon
            wavelength units:     AA
            central wavelength:   45110.141614 angstrom
            pivot wavelength:     45020.219955 angstrom
            effective wavelength: 44425.747085 angstrom
            photon wavelength:    44603.204646 angstrom
            minimum wavelength:   39250.000000 angstrom
            maximum wavelength:   50640.000000 angstrom
            norm:                 4664.680820
            effective width:      8714.143135 angstrom
            fullwidth half-max:   10110.000000
            definition contains 417 points

            Zeropoints
                Vega: 28.933083 mag,
                      2.6715713304827696e-12 erg / angstrom * centimeter ** 2 * second,
                      180.61811118349118 Jy
                      3.246733893355585 photon / angstrom * centimeter ** 2 * second
                  AB: 25.674986 mag,
                      5.370385702161592e-11 erg / angstrom * centimeter ** 2 * second,
                      3630.780547701007 Jy
                  ST: 21.100000 mag,
                      3.6307805477010028e-09 erg / angstrom * centimeter ** 2 * second,
                      245467.79536259372 Jy

        [...]

Suppose one has a calibrated spectrum and wants to compute the vega magnitude
through the HST/WFC3 F110W passband: 

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

As mentioned in the above, sometimes a spectrum of reference of Vega is necessary. 

We use the synthetic spectrum provided by `Bohlin 2007 <https://ui.adsabs.harvard.edu/abs/2007ASPC..364..315B/abstract>`_, a common reference
througout many photometric suites.

The interface to the Vega template is given through the :class:`pyphot.vega.Vega` class.

Internal Sun reference
~~~~~~~~~~~~~~~~~~~~~~~

We also provide observed and theoretical references for the solar spectrum following 
`Colina, Bohlin, & Castelli 1996 <https://ui.adsabs.harvard.edu/abs/1996AJ....112..307C/abstract>`_.

The observed solar spectrum comes from CALSPEC 
`sun_reference_stis_001.fits <ftp://ftp.stsci.edu/cdbs/current_calspec/sun_reference_stis_001.fits>`_ 
which provides the ultraviolet to near-infrared absolute flux distribution of
the Sun covering the 0.12-2.5 μm wavelength range. The solar reference spectrum
combines absolute flux measurements from satellites and from the ground with a
model spectrum for the near-infrared. 

The theoretical spectrum comes from the Kurucz'93 atlas:
`<sun_kurucz93.fits ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/sun_kurucz93.fits>`_
The theoretical spectrum is scaled to match the observed spectrum from 1.5 - 2.5 microns, and then it is used where the observed spectrum ends.
The theoretical model of the Sun uses the following parameters when the Sun is at 1 au:


.. list-table:: Solar Parameters
        :header-rows: 1

        * - log_Z
          - T_eff
          - log_g
          - V_{Johnson}
        * - +0.0
          - 5777
          - +4.44
          - -26.75

The interface to the Sun templates is given through the :class:`pyphot.sun.Sun` class.

Contributors
~~~~~~~~~~~~

Author: Morgan Fouesneau (`@mfouesneau <https://github.com/mfouesneau>`_)

Direct contributions to the code base:

* Ariane Lancon (`@lancon <https://github.com/lancon>`_)
* Tim Morton (`@timothydmorton <https://github.com/timothydmorton>`_)

Related projects
~~~~~~~~~~~~~~~~

* `cphot <https://github.com/mfouesneau/cphot>`_ is a spin-off project that provides a C++ version of pyphot


Citation
~~~~~~~~

If you use this software in your work, please cite it using the following:

.. tabs::

        .. tab:: bibtex

                .. code-block::

                        @software{zenodopyphot,
                        author       = {Fouesneau, Morgan},
                        title        = {pyphot},
                        month        = jan,
                        year         = 2025,
                        publisher    = {Zenodo},
                        version      = {pyphot\_v1.6.0},
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

