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
   modules


Installation
~~~~~~~~~~~~
Pyphot is available on `PyPI <https://pypi.org/project/pyphot/>`_ and can be installed using any common package manager.

* Using pip: (Use the `--user` option if you don't have permissions to install libraries)

.. image:: https://img.shields.io/pypi/v/pyphot.svg
    :target: https://pypi.org/project/pyphot/

.. code-block:: none

        pip install git+https://github.com/mfouesneau/pyphot

* Manually from latest source on GitHub:

.. code-block:: none

        git clone https://github.com/mfouesneau/pyphot
        cd pyphot
        python setup.py intall


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

How to contribute
~~~~~~~~~~~~~~~~~

This project is open source and new contributions and contributors are very welcome!

Please open a new issue or new pull request for bugs, feedback, or new features you would like to see. If there is an issue you would like to work on, please leave a comment, and we will be happy to assist.

Please have a look at our `code of conduct <https://github.com/mfouesneau/pyphot/blob/main/CODE_OF_CONDUCT.md>`_.


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
