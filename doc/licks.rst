Extention to Lick indices
=========================

We also include functions to compute lick indices and provide a series of
commonly use ones.

Lick index
~~~~~~~~~~

The Lick system of spectral line indices is one of the most commonly used
methods of determining ages and metallicities of unresolved (integrated light)
stellar populations.

The calibration of the Lick / IDS system is complicated because the original
Lick spectra were not flux calibrated, so there are usually systematic effects
due to differences in continuum shape.  Proper calibration involves observing
many of the original Lick/IDS standard stars and deriving offsets to the
standard system.

In `Vazdekis et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010MNRAS.404.1639V>`_, they propose a new Line Index System, hereafter
LIS, with three new spectral resolutions at which to measure the Lick
indices. Note that this new system should not be restricted to the Lick set
of indices in a flux calibrated system. In fact, LIS can be used for any
index in the literature (e.g., for the Rose (1984) indices), including
newly defined indices (e.g., Cervantes & Vazdekis 2009).


The LIS system is defined for 3 different spectral resolutions which are
best suited for the following astrophysical cases:

* LIS-5.0AA: globular clusters
* LIS-8.4AA: low and intermediate-mass galaxies
* LIS-14.0AA: massive galaxies

Conversions to transform the data from the Lick/IDS system to LIS can be found
in Johansson, Thomas & Maraston (2010), which provides a discussion of indices
and the information content of them.


Quick start example
~~~~~~~~~~~~~~~~~~~

The lick extension is very similar to the broadband usage.


.. code-block:: python

        # convert to magnitudes
        import numpy as np
        from pyphot import LickLibrary
        # using the internal collection of indices
        lib = LickLibrary()
        f = lib['CN_1']
        # work on many spectra at once
        index = f.get(lamb, spectra, axis=1)



Calculations
~~~~~~~~~~~~

Suppose one has a spectrum :math:`f_\lambda` defined over the wavelength
:math:`\lambda`.  First, we must adapt the resolution of the spectrum to match
one of the LIS range.  Indeed Lick definitions have different resolution
elements as function of wavelength:

+---------------------------------------------+------+------+------+-----+-------+
| :math:`\lambda` in :math:`\unicode{x212B}`  | 4000 | 4400 | 4900 | 5400|  6000 |
+---------------------------------------------+------+------+------+-----+-------+
| resolution (FWHM in :math:`\unicode{x212B}`)| 11.5 | 9.2  | 8.4  | 8.4 |  9.8  |
+---------------------------------------------+------+------+------+-----+-------+


.. code-block:: python

        new_f = licks.reduce_resolution(w, f, sigma0=0.55, sigma_floor=0.2)


Indices are defined on continuum normalized spectra. Therefore all indices come
with 3 intervals: a `band` that gives the index range but also a `blue` and a
`red` interval on each side which are used to fit a polynomial function as local
continuum.

see: :func:`licks.LickIndex.continuum_normalized_region_around_line`.


Finally any index is calculated by integrated the continuum normalized flux.
Some indices are given in `magnitudes`, and some in `equivalent width` units.


References
~~~~~~~~~~

* Johansson, Thomas & Maraston, 2010, MNRAS, 406, 165, "Empirical calibrations of optical absorption-line indices based on the stellar library MILES" `2010MNRAS.406..165J <https://ui.adsabs.harvard.edu/abs/2010MNRAS.406..165J>`_

* Puzia et al. 2002, A&A, 395, 45, "Integrated spectroscopy of bulge globular clusters and fields. I. The data base and comparison of individual Lick indices in clusters and bulge" `2002A&A...395...45P <https://ui.adsabs.harvard.edu/abs/2002A&A...395...45P>`_ 

* Vazdekis et al., 2010, MNRAS, 404, 1639, "Evolutionary stellar population synthesis with MILES - I. The base models and a new line index system" `2010MNRAS.404.1639V <https://ui.adsabs.harvard.edu/abs/2010MNRAS.404.1639V>`_

* Worthey G., Faber S. M., Gonzalez J. J., Burstein D., 1994, ApJS, 94, 687, "Old Stellar Populations. V. Absorption Feature Indices for the Complete Lick/IDS Sample of Stars" `1994ApJS...94..687W <https://ui.adsabs.harvard.edu/abs/1994ApJS...94..687W>`_

* Worthey G., Ottaviani D. L., 1997, ApJS, 111, 377, "Hγ and Hδ Absorption Features in Stars and Stellar Populations" `1997ApJS..111..377W <https://ui.adsabs.harvard.edu/abs/1997ApJS..111..377W>`_

* Zhang, Li & Han 2005 MNRAS, 364, 503, "Evolutionary population synthesis for binary stellar population at high spectral resolution: integrated spectral energy distributions and absorption-feature indices", `2005MNRAS.364..503Z <https://ui.adsabs.harvard.edu/abs/2005MNRAS.364..503Z>`_


This page shows the content of the provided library with respective properties
of the passband filters. The code to generate the table is also provided below.


Library content
~~~~~~~~~~~~~~~

.. literalinclude:: licks_table.py
   :language: python

.. csv-table:: Current internal library
        :file: licks_table.csv
