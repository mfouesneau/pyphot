More example usage
==================


Notebooks
~~~~~~~~~

check notebooks in the `examples directory`_ in the repository.

.. _examples directory: https://github.com/mfouesneau/pyphot/tree/master/examples

You can also run these online

.. image:: https://mybinder.org/badge.svg 
  :target: https://mybinder.org/v2/gh/mfouesneau/pyphot/master?filepath=examples%2FQuickStart.ipynb
 
.. image:: https://img.shields.io/badge/render%20on-nbviewer-orange.svg
  :target: https://nbviewer.jupyter.org/github/mfouesneau/pyphot/tree/master/examples/
  

Basic example
~~~~~~~~~~~~~

This example reads a spectrum from a fits file and calculate some photometry
using units to the calculations.

.. code:: python 

        from astropy.io import fits
        from pyphot import unit
        import pyphot

        wavelength, flux = fits.getdata("foo.fits").T
        wavelength  = wavelength * unit['AA']
        flux = flux * unit['erg/s/cm**2/AA']
        lib = pyphot.get_library()
        f = lib['HST_WFC3_F110W']
        mag = -2.5 * np.log10(f.get_flux(wavelength, flux)) - f.Vega_zero_mag

Unit conversions
~~~~~~~~~~~~~~~~

Pyphot contains a simple unit package based on `pint`_. (One day maybe `astropy`)
It allows handling wavelength and flux units very easily.

.. code:: python

        from pyphot import Vega 

        vega = Vega()
        wavelength_angstrom = vega.wavelength
        wavelength_nm = vega.wavelength.to('nm')
        flux_flam = vega.flux
        flux_other = flux_flam.to('W / (m**2 * angstrom)')


.. _pint: https://pint.readthedocs.io/en/0.9/


Number of photons through a given passband
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following examples gives the number of photons that a passband would collect
for the Sun at 10 pc.

.. code:: python

        from pyphot import Sun, unit
        with Sun(distance=10 * unit['pc']) as v:
            print(vband.get_Nphotons(v.wavelength, v.flux))

        11.985932608649655 photon / angstrom


Defining your own passband
~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's suppose you're interested in defining bandpasses (here tophat functions)
in Pyphot to determine the flux/magnitude in regions. This may be useful to determine the S/N ratio or integrated flux over a given wavelength range.

.. code:: python

        from pyphot import (unit, Filter)

        wave = [4499, 4500, 4700, 4701] * unit['AA']
        transmit = [0., 1., 1., 0.]
        tophat = Filter(wave, transmit, name='tophat', dtype='photon', unit='Angstrom')
        flux_tophat = tophat.get_flux(wavelength, flux)

.. tip::

        * The unit in the filter declaration is redundant with the unit on wave
          in this example.
        * the transmission is unitless. Its scaling factor will
          not affect the flux/magnitude calculations as these are implicitly
          normalizing the passband  (see the equations) but will affect the
          number of photon :math:`photons/s/cm^2`


 
Dust attenuated SED modeling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The associated notebook can be found `online`_. It requires a bit more than just
`Pyphot`, but it demonstrates how to integrate `Pyphot` into larger projects.

.. image:: https://mybinder.org/badge.svg 
  :target: https://mybinder.org/v2/gh/mfouesneau/GaiaSprint2018/master

.. _online: https://github.com/mfouesneau/GaiaSprint2018/blob/master/dust_attenuated_seds.ipynb


Flux density units
~~~~~~~~~~~~~~~~~~

* What are the expected units on the spectra?

This main version of `Pyphot` assumes that the spectral density is in the units
of the wavelength of the spectrum. Also it expects that the flux is provided in
`flam`, i.e., :math:`erg/s/cm^2/\unicode{x212B}`.

Conversions can be done using `pyphot.unit`, however it may be tedious
sometimes.

If you are willing to live at the bleeding edge, `pyphot.sandbox` contains the
next generation of `pyphot` which handles the flux density conversions. With a
single import line, you can switch to the new version.

.. code:: python

        from pyphot import sandbox as pyphot


Check the examples associated with this new version online
in the `examples directory`_ in the repository.
