.. pyphot documentation master file, created by
   sphinx-quickstart on Wed Oct  5 11:25:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PYPHOT -- A tool for computing photometry from spectra
======================================================

This is a set of tools to compute synthetic photometry in a simple way, suitable
for integration in larger projects.

The inputs are response functions for the desired photometric passbands, and
stellar spectra.

Filters are represented individually by a `Filter` object. Collections of
filters are handled with a `Library`. We provide an internal library that
contains a signitificant amount of common filters. Note that we use the word
"filter" in a generic sense: in many cases it corresponds in fact to a total
throughput, that includes the actual filter transmission, the transmission of
the optics, the wavelength-dependent efficiency of the detector, and (for
ground-based photometric systems) the transmission of the atmosphere.

Each filter is minimally defined by a `wavelength`, a `throughput`, and a
detector type that may be either `energy` or `photon` (default). Many properties
such as central or pivot wavelengths are computed internally.

Units are provided for the wavelengths of the filters, zero points in multiple
units are also accessible (AB, Vega magnitude, Jy, erg/s/cm2/AA). The default
detector is assumed to be a photon-counting device, but energy-sensitive
detectors are also handled for the computations.

.. note::

    **New in v2.0.0**

    Pyphot handles units from `Astropy Units <https://docs.astropy.org/en/stable/units/index.html>`_ (as default), `Pint <https://pint.readthedocs.io/en/stable/>`_, and `EzUnits`_ (:mod:`~pyphot.unit_adapters.ezunits`).
    You can switch between them using the :func:`pyphot.config.set_units_backend` function.

    .. code-block:: python

        import pyphot.config
        pyphot.config.set_units_backend('astropy')
        pyphot.config.set_units_backend('pint')
        pyphot.config.set_units_backend('ezunits')

.. image:: https://mybinder.org/badge.svg
  :target: https://mybinder.org/v2/gh/mfouesneau/pyphot/master?filepath=examples%2FQuickStart.ipynb

.. image:: https://img.shields.io/badge/render%20on-nbviewer-orange.svg
  :target: https://nbviewer.jupyter.org/github/mfouesneau/pyphot/tree/master/examples/


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

Contents
---------

.. toctree::
   :maxdepth: 1

   Home <index>
   more_examples
   photometry
   libcontent
   licks
   vega
   modules


Installation
~~~~~~~~~~~~
Pyphot is available on `PyPI <https://pypi.org/project/pyphot/>`_ and can be installed using any common package manager.



* Using pip: Use the `--user` option if you don't have permissions to install libraries

.. image:: https://img.shields.io/pypi/v/pyphot.svg
    :target: https://pypi.org/project/pyphot/

.. code-block:: none

        pip install git+https://github.com/mfouesneau/pyphot

* Manually from source on GitHub:

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

Suppose one has a calibrated spectrum and wants to compute the corresponding vega magnitude
through the HST/WFC3 F110W passband:

.. code-block:: python

        import numpy as np
        f = lib['hst_wfc3_f110w']
        # compute the integrated flux through the filter f
        # (note that 'spectra' can be a 2D array holding many spectra)
        fluxes = f.get_flux(lamb, spectra, axis=1)
        # convert fluxes to vega magnitudes
        mags = -2.5 * np.log10(fluxes) - f.Vega_zero_mag
        # or similarly
        mags = -2.5 * np.log10(fluxes / f.Vega_zero_flux)


To use one's own transmission curve, defined by `lamb_T` and
`T`, instead of one of the pre-defined ones, one would use
the :class:`pyphot.phot.Filter` directly:

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

By default, we use the synthetic spectrum `alpha_lyr_stis_003` provided by `Bohlin 2007 <https://ui.adsabs.harvard.edu/abs/2007ASPC..364..315B/abstract>`_, a common reference
througout many photometric suites.

The interface to the Vega template is given through the :class:`pyphot.vega.Vega` class.

Since v.1.7.0, additional flavors and description of the internal Vega reference can be found on the :doc:`vega` page.

To use a specific Vega flavor for the photometric calculations in Pyphot, you can set the `vega` keyword parameter  when creating a passband or use the `set_vega_flavor` method to update it.
For example, to use the `alpha_lyr_stis_011` flavor when creating a passband filter, you can do the following:

.. code-block:: python

    from pyphot.astropy import UnitFilter, Unit as u

    # Create a passband using the Vega flavor
    pb = UnitFilter(
        [4000, 5000, 6000] * u.AA,
        [0.1, 0.8, 0.1],
        name="Example Passband",
        dtype="photon"
        vega="stis_011"  # Specify the Vega flavor,
    )

Or alternatively, you can set/reset the Vega flavor after:

.. code-block:: python

    from pyphot import svo

    pb = svo.get_pyphot_filter("GALEX/GALEX.FUV")
    # Set the Vega flavor to use
    pb.set_vega_flavor("stis_011")


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

* Ariane Lançon (`@lancon <https://github.com/lancon>`_)
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
