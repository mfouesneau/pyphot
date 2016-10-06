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

**When units** are provided for the wavelength, zero points in multiple units are also
accessible (AB, Vega magnitude, Jy, erg/s/cm2/AA). The default detector type is
assumed to be photonic, but energetic detectors are also handled for the
computations.

.. tip::

        All provided filters have defined units and detector type
        that are used transparently throughout the package.


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

Contents
---------

.. toctree::
   :maxdepth: 1

   libcontent
   modules




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
        

Details on predicting photometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is sometimes not obvious that there are important differences between
photometric systems. But even less known the difference between detector types
which requires also special care.

In this section, we review the important details for computing the luminosity
and the magnitude of a star through a photometric passband filter.

If we consider a filter throughput defined in wavelength by the dimensionless
function :math:`T(\lambda)`, 
this function tells you what fraction of the arriving photons at wavelength
:math:`\lambda` actually get through the instrument.  Therefore, the total number of
photons, per unit time per unit area, received in this filter is

.. math::

        \begin{equation}
        N_{tot} = \frac{1}{hc} \int_\lambda f_\lambda\,\lambda\,T(\lambda)\,d\lambda,
        \end{equation}

where :math:`f_\lambda` is the wavelength dependent flux density of an object
given in energy per unit time per unit area per unit wavelength.

Consequently, interpreting :math:`\lambda T(\lambda)` as a distribution leads to
the statistical mean of the flux density, :math:`\overline{f_\lambda}` 

.. math::

        \begin{equation}
        \overline{f_\lambda}(T) = \frac{\int_\lambda \lambda f_\lambda T(\lambda) d\lambda}{\int_\lambda \lambda T(\lambda) d\lambda}.
        \end{equation}

Note that this is not the mean flux density because of the :math:`\lambda` factor in
the integrals, it is the mean photon rate density in this filter commonly
expressed in stellar physics literature as :math:`erg/s/cm^2/\unicode{x212B}` or 
:math:`W/m^2/\unicode{x212B}`.

.. code-block:: python

        # computing the flux of a spectrum
        flux = lib['hst_wfc3_f110w'].get_flux(lamb, spec)
        # lamb may have units, otherwise assuming consistent definitions.

        # computing the flux of many spectra
        fluxes = lib['hst_wfc3_f110w'].get_flux(lamb, spectra, axis=1)


Finally, at least for instruments using CCD or CCD-like cameras, i.e., counting
photons, we obtain the usual definition of magnitude 

.. math::

        \begin{equation}
        mag_\lambda(T) = -2.5\,\log_{10}\left(\overline{f_\lambda}\right) - ZP\left(\overline{f_\lambda}\right),
        \end{equation}

where :math:`ZP(\overline{f_\lambda})` gives the pass-band reference value
(zeropoint) for a given photometric/magnitude system.

However, the zeropoints themselves depend on the adopted photometric system used
to report the measurements. They may vary fundamentally from one to another.
Below we briefly describe the main systems used in large surveys.



Vega magnitude system
---------------------

This system is defined such that the star Alpha Lyr (Vega) has magnitude 0 in
any pass-band filter. In other words, the zeropoints are set by the magnitude of
vega, :math:`-2.5 \log_{10} \overline{f_\lambda}(Vega)`, or

.. math:: 

        \begin{equation}
        mag_{Vega}(T) = -2.5\,\log_{10}\left(\overline{f_\lambda} / \overline{f_\lambda}(Vega)\right).
        \end{equation}

.. code-block:: python

        # convert to magnitudes
        import numpy as np
        f = lib['hst_wfc3_f110w']
        fluxes = f.get_flux(lamb, spectra, axis=1)
        mags = -2.5 * np.log10(fluxes) - f.Vega_zero_mag
        # or similarly
        mags = -2.5 * np.log10(fluxes / f.Vega_zero_flux)

       

Johnson system
--------------

The Johson system is defined such that the star Alpha Lyr (Vega) has :math:`V=0.03`
and all colors equal to zero. It is very similar to the Vega magnitude system,
but using mean flux definition (instead of photon counts), i.e., **energy
counter** detectors

.. math::

        \begin{equation}
        \widetilde{f_\lambda}(T) = \frac{\int_\lambda f_\lambda T(\lambda) d\lambda}{\int_\lambda T(\lambda) d\lambda},
        \label{eq:Johnsonmag}
        \end{equation}

(the true definition of mean flux throughout a given transmission filter.)

.. note::

        Table A2 of Bessell et al. 1998 gives zero points for the UBVRIJHKL(+Kp and L’) filters in the Counsins-Glass-Johnson system.

If one defines the **effective wavelength** :math:`\lambda_{\rm eff}` as the
photon weighted mean wavelength:

.. math::

        \lambda_{\rm eff} = \frac{\int \lambda f_\lambda T(\lambda) d\lambda}{\int f_\lambda T(\lambda) d\lambda},

.. code-block:: python

        # the effective wavelength for vega is given by
        lib['ground_johnson_u'].leff


then the difference between the Johnson and Vega systems within the same filter is given by

.. math:: 

        \begin{equation}
        \widetilde{mag}_\lambda - \overline{mag}_\lambda = 0.03 - 2.5 \log_{10} \frac{\lambda_{\rm eff}(Vega)}{\lambda_{\rm eff}(star)},
        \end{equation}

where we explicit which equation was used to compute magnitudes.



.. code-block:: python

        # The switch between the energy and the photon count equation is done
        # through the `Filter.set_dtype` method, and becomes transparent for any
        # use. So if you define you own filter either use the constructor or the
        # method

        # define a constant filter in energy count from 100 to 110 AA
        f = Filter(np.arange(100, 110), np.ones(10), \
                        dtype='energy', unit='AA')
        # manually set the detector type
        f.set_dtype('photon')



AB magnitude system
-------------------

This system is defined such that, when monochromatic flux :math:`f_\nu` is measured in
:math:`erg\,s^{-1}\,cm^{-2} Hz^{-1}`,

.. math::

        mag_{AB}(T) = -2.5\, \log_{10}(\overline{f_\nu}) - 48.60

where the value of the constant is selected to define :math:`m_{AB}=V` for a
flat-spectrum source. In this system, an object with constant flux per unit
frequency interval has zero color.

Koornneef et al. gives the respective definition of :math:`\overline{f_\nu}(T)`:

.. math::

        \begin{equation}
        \overline{f_\nu}(T) = \frac{\int_\nu f_\nu T(\nu) d\nu / \nu}{\int_\nu T(\nu) d\nu / \nu}
         = \frac{\int_\lambda f_\nu T(\lambda) d\lambda / \lambda}{\int_\lambda T(\lambda) d\lambda / \lambda}
        \end{equation}

To go back to wavelength units, we have :math:`d\nu = (c/\lambda^2) d\lambda`.

If one defines the **pivot wavelength** :math:`\lambda_p` to convert between
:math:`\overline{f_\nu}` and :math:`\overline{f_\lambda}` as

.. math::

        \begin{equation}
        \overline{f_\nu} = \frac{\lambda_p^2}{c} \overline{f_\lambda},
        \end{equation}

one can easily show that

.. math::

        \begin{equation}
        \lambda_p^2 = \frac{\int_\lambda T(\lambda)\,\lambda\,d\lambda}{\int_\lambda T(\lambda)\,d\lambda /\lambda}.
        \end{equation}

Therefore for filters with AB magnitudes, one can compute 

.. math::

        \begin{equation}
        mag_{AB}(T) = -2.5\, \log_{10}(\overline{f_\lambda}) - 2.5\log_{10}\left(\lambda_p^2/c\right) - 48.6,
        \end{equation}

where care must be taken to use the speed of light :math:`c` and :math:`\lambda_p` in matching units.


.. code-block:: python

        # convert to magnitudes
        import numpy as np
        f = lib['hst_wfc3_f110w']
        fluxes = f.get_flux(lamb, spectra, axis=1)
        mags = -2.5 * np.log10(fluxes) - f.AB_zero_mag
        # or similarly
        mags = -2.5 * np.log10(fluxes / f.AB_zero_flux)



ST magnitude system
-------------------

This system is defined such as a source with flat :math:`f_\lambda` will have
the same magnitude in every filter. 

Koornneef et al. (1986; same as above) defines 

.. math::

        \begin{equation}
        mag_{ST}(T) = -2.5\, \log_{10}(\overline{f_\lambda}) - 21.1,
        \end{equation}


.. code-block:: python

        # convert to magnitudes
        import numpy as np
        f = lib['hst_wfc3_f110w']
        fluxes = f.get_flux(lamb, spectra, axis=1)
        mags = -2.5 * np.log10(fluxes) - f.ST_zero_mag
        # or similarly
        mags = -2.5 * np.log10(fluxes / f.ST_zero_flux)


Jansky definition
-----------------

References
----------

* Bessel, M. S. 1990, PASP, 91, 589;

* Bessel, M. S. 1983, PASP, 95, 480;

* Bessel, M. S. 1990, PASP, 102, 1181;

* Hayes, D. S., \& Latham, D. W. 1975, ApJ, 197, 593;

* Johnson, H. L. \& Morgan, W. W. 1953, ApJ, 117, 313

* Oke, J.B. 1974, ApJS, 27, 21;

* Koornneef, Bohlin, Buser, Horne, Turnshek : Synthetic photometry and the calibration of HST.





Internal Vega reference
~~~~~~~~~~~~~~~~~~~~~~~

As mentioned in the above, sometimes a spectrum of reference of Vega is
necessary. 

We use the synthetic spectrum provided by Bohlin 2007, a common reference
througout many photometric suites.

The interface to the Vega template is given through the :class:`pyphot.vega.Vega` class.





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

