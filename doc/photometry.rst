Details on predicting photometry
================================

It is sometimes not obvious that there are important differences between photometric systems. But even less known is the difference between detector count types (energy or photons) which requires also special care.

In this section, we review the important details for computing the luminosity and the magnitude of a star through a photometric passband. We do not discuss calibration, which in principle is covered by instrument documentations.

Synthetic photometry uses transmission curves that typically include a variety of wavelength-dependent components (a filter transmission, the response of the optical elements of the telescope and instrument, the detector sensitivity, sometimes a telluric absorption component, ...). Transmission curves are usually named after the associated filter, and it is up to the user to verify that the other components are also included. In this documentation, we use the terms filter throughput, transmission curve, response function, or simply filter, interchangeably, to refer to a transmission curve, generally leaving the verification of the nature of these transmission curves to the user. 

Most modern detectors count photons, rather than cumulating the energies of these photons.  If we consider a filter throughput (a.k.a, transmission curve, or response function) defined in wavelength by the dimensionless function :math:`T(\lambda)`, this function tells you what fraction of the arriving photons at wavelength :math:`\lambda` actually get detected (on average).  Therefore, the total number of photons, per unit time per unit collecting area, expected to be detected in this filter is

.. math::

        \begin{equation}
        N_{tot} = \frac{1}{hc} \int_\lambda f_\lambda\,\lambda\,T(\lambda)\,d\lambda,
        \end{equation}

where :math:`f_\lambda` is the wavelength dependent flux density of an object given in energy per unit time, per unit collecting area, and per unit wavelength.

Consequently, interpreting :math:`\lambda T(\lambda)` as a distribution leads to the statistical mean of the flux density, :math:`\overline{f_\lambda}` 

.. math::

        \begin{equation}
        \overline{f_\lambda}(T) = \frac{\int_\lambda \lambda f_\lambda T(\lambda) d\lambda}{\int_\lambda \lambda T(\lambda) d\lambda}.
        \end{equation}

Note that although this is formally a weighted mean of a flux density, with weights proportional to :math:`\lambda T(\lambda)` (the denominator ensuring that the sum of the weights is 1), it actually measures the mean photon rate density in this filter. The result is commonly expressed in the same units as :math:`f_\lambda`, i.e. :math:`erg/s/cm^2/\unicode{x212B}` or :math:`W/m^2/\unicode{x212B}`.

.. code-block:: python

        # computing the mean flux of a spectrum 
        # single_spectrum: np.array of shape (n_lambda,)
        flux = lib['hst_wfc3_f110w'].get_flux(λ, single_spectrum)
        # λ may have units, otherwise assuming consistent definitions.

        # computing the mean flux of many spectra
        # spectra: np.array of shape (n_spectra, n_lambda)
        fluxes = lib['hst_wfc3_f110w'].get_flux(λ, spectra, axis=1)


Finally, at least for instruments using CCD or CCD-like cameras, i.e., counting photons, we obtain the usual definition of magnitude 

.. math::

        \begin{equation}
        mag_\lambda(T) = -2.5\,\log_{10}\left(\overline{f_\lambda}\right) - ZP\left(\overline{f_\lambda}\right),
        \end{equation}

where :math:`ZP(\overline{f_\lambda})` gives the passband reference value (zeropoint) for a given photometric/magnitude system.

However, the zeropoints themselves depend on the photometric system adopted to report the measurements. They may vary fundamentally from one to another.  Below we briefly describe the main systems used in large surveys.


Vega magnitude system
~~~~~~~~~~~~~~~~~~~~~

This system is defined such that the star Alpha Lyr (Vega) has magnitude 0 in any pass-band filter. In other words, the zeropoints are set by the magnitude of vega, :math:`-2.5 \log_{10} \overline{f_\lambda}(Vega)`, or

.. math:: 

        \begin{equation}
        mag_{Vega}(T) = -2.5\,\log_{10}\left(\overline{f_\lambda} / \overline{f_\lambda}(Vega)\right).
        \end{equation}

.. code-block:: python
        # convert to magnitudes
        import numpy as np
        f = lib['hst_wfc3_f110w']
        fluxes = f.get_flux(lamb, spectra, axis=1)
        # careful taking the log of non-dimensionless Quantity! -> use .value
        mags = -2.5 * np.log10(fluxes.value) - f.Vega_zero_mag
        # or similarly
        mags = -2.5 * np.log10(fluxes / f.Vega_zero_flux)


We use the synthetic spectrum `alpha_stis_003` provided by `Bohlin 2007 <https://ui.adsabs.harvard.edu/abs/2007ASPC..364..315B/abstract>`_, a common reference througout many photometric suites. Additional flavors and description of the internal Vega reference can be found on the :doc:`vega` page.
       

Johnson system
~~~~~~~~~~~~~~

The Johnson system is defined such that the star Alpha Lyr (Vega) has :math:`V=0.03` mag and all colors equal to zero. It is very similar to the Vega magnitude system, but using mean flux definition (instead of photon counts), as appropriate for historical **energy counter** detectors

.. math::

        \begin{equation}
        \widetilde{f_\lambda}(T) = \frac{\int_\lambda f_\lambda T(\lambda) d\lambda}{\int_\lambda T(\lambda) d\lambda},
        \label{eq:Johnsonmag}
        \end{equation}

(this is the mean flux weighted simply by the normalized throughout.)

.. note::

        Table A2 of `Bessell et al. (1998) <https://ui.adsabs.harvard.edu/abs/1998A%26A...333..231B>`_ gives zero points for the UBVRIJHKL(+Kp and L') filters in the Counsins-Glass-Johnson system.

If one defines the **effective wavelength** :math:`\lambda_{\rm eff}` as the photon weighted mean wavelength:

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

        from pyphot import Filter
        from pyphot.config import units
        import numpy as np
        # define a constant filter in energy count from 100 to 110 AA
        f = Filter(np.arange(100, 110), np.ones(10), dtype='energy', unit='AA')
        # or equivalently
        f = Filter(np.arange(100, 110) * units.U("AA"), np.ones(10), dtype='energy')
        # manually set the detector type
        f.set_dtype('photon')



AB magnitude system
~~~~~~~~~~~~~~~~~~~

This system is defined such that, when monochromatic flux :math:`f_\nu` is measured in :math:`erg\,s^{-1}\,cm^{-2} Hz^{-1}`,

.. math::

        mag_{AB}(T) = -2.5\, \log_{10}(\overline{f_\nu}) - 48.60

where the value of the constant is selected to define :math:`m_{AB}=V` for a flat-spectrum source. In this system, an object with constant flux per unit frequency interval has zero color.

`Koornneef et al. (1986) <https://ui.adsabs.harvard.edu/abs/1986HiA.....7..833K>`_ gives the respective definition of :math:`\overline{f_\nu}(T)`:

.. math::

        \begin{equation}
        \overline{f_\nu}(T) = \frac{\int_\nu f_\nu T(\nu) d\nu / \nu}{\int_\nu T(\nu) d\nu / \nu}
         = \frac{\int_\lambda f_\nu T(\lambda) d\lambda / \lambda}{\int_\lambda T(\lambda) d\lambda / \lambda}
        \end{equation}

To go back to wavelength units, we have :math:`d\nu = (c/\lambda^2) d\lambda`.

If one defines the **pivot wavelength** :math:`\lambda_p` to convert between :math:`\overline{f_\nu}` and :math:`\overline{f_\lambda}` as

.. math::

        \begin{equation}
        \overline{f_\nu} = \frac{\lambda_p^2}{c} \overline{f_\lambda},
        \end{equation}

one can show that

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
~~~~~~~~~~~~~~~~~~~

This system is defined such as a source with flat :math:`f_\lambda` will have
the same magnitude in every filter. 

`Koornneef et al. (1986; same as above) <https://ui.adsabs.harvard.edu/abs/1986HiA.....7..833K>`_ defines 

.. math::

        \begin{equation}
        mag_{ST}(T) = -2.5\, \log_{10}(\overline{f_\lambda}) - 21.1,
        \end{equation}


.. code-block:: python

        # convert to magnitudes
        import numpy as np
        f = lib['hst_wfc3_f110w']
        fluxes = f.get_flux(lamb, spectra, axis=1)
        # careful taking the log of non-dimensionless Quantity! -> use .value
        mags = -2.5 * np.log10(fluxes.value) - f.ST_zero_mag
        # or similarly
        mags = -2.5 * np.log10(fluxes / f.ST_zero_flux)


Jansky definition
~~~~~~~~~~~~~~~~~

The jansky (symbol Jy) is a non-SI unit of spectral flux density, it is equivalent to :math:`10^{−26} W.m^{-2}.Hz^{-1}` or :math:`10^{-23} erg/s/cm^2/Hz`.

.. math::

        \begin{equation}
        {f_{Jy}} = \frac{10^5}{10^{-8}c} {\lambda_p^2} {f_\lambda},
        \end{equation}

where :math:`c` is the speed of light in :math:`m/s`,  :math:`\lambda_p` is the pivot wavelength in :math:`Å`, and :math:`{f_\lambda}` the flux (Vega, AB, or ST) in flam (:math:`erg.s^{-1}.cm^{-2}.Å^{-1}`).

.. code-block:: python

        f = lib['hst_wfc3_f110w']
        print(f.AB_zero_Jy, f.Vega_zero_Jy, f.ST_zero_Jy)


References
~~~~~~~~~~

* Bessell, M. S. 1983, PASP, 95, 480, "VRI photometry : an addendum." `1983PASP...95..480B <https://ui.adsabs.harvard.edu/abs/1983PASP...95..480B>`_;

* Bessell, M. S. 1990, PASP, 102, 1181, "UBVRI passbands" `1990PASP..102.1181B <https://ui.adsabs.harvard.edu/abs/1990PASP..102.1181B>`_;

* Bessell, M. S., Castelli, F., \& Plez, B. 1998, A&A, 333, 231, "Model atmospheres broad-band colors, bolometric corrections and temperature calibrations for O - M stars." `1998A&A...333..231B <https://ui.adsabs.harvard.edu/abs/1998A%26A...333..231B/abstract>`_;

* Hayes, D. S., \& Latham, D. W. 1975, ApJ, 197, 593, "A rediscussion of the atmospheric extinction and the absolute spectral-energy distribution of Vega." `1975ApJ...197..593H <https://ui.adsabs.harvard.edu/abs/1975ApJ...197..593H>`_;

* Johnson, H. L. \& Morgan, W. W. 1953, ApJ, 117, 313, "Fundamental stellar photometry for standards of spectral type on the Revised System of the Yerkes Spectral Atlas." `1953ApJ...117..313J <https://ui.adsabs.harvard.edu/abs/1953ApJ...117..313J>`_;

* Koornneef, Bohlin, Buser, Horne, Turnshek : Synthetic photometry and the calibration of HST. `1986HiA.....7..833K <https://ui.adsabs.harvard.edu/abs/1986HiA.....7..833K>`_

* Oke, J.B. 1974, ApJS, 27, 21, "Absolute Spectral Energy Distributions for White Dwarfs" `1974ApJS...27...21O <https://ui.adsabs.harvard.edu/abs/1974ApJS...27...21O>`_;



