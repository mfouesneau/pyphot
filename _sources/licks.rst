Extention to Lick indices
=========================

We also include functions to compute lick indices and provide a series of commonly use ones.

Lick index
~~~~~~~~~~

The Lick system of spectral line indices is one of the most commonly used methods of determining ages and metallicities of unresolved (integrated light) stellar populations.

The calibration of the Lick / IDS system is complicated because the original Lick spectra were not flux calibrated, so there are usually systematic effects due to differences in continuum shape.  Proper calibration involves observing many of the original Lick/IDS standard stars and deriving offsets to the standard system.

In `Vazdekis et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010MNRAS.404.1639V>`_, they propose a new Line Index System, hereafter
LIS, with three new spectral resolutions at which to measure the Lick indices. Note that this new system should not be restricted to the Lick set of indices in a flux calibrated system. In fact, LIS can be used for any index in the literature (e.g., for the `Rose (1984) indices <https://ui.adsabs.harvard.edu/abs/1985AJ.....90.1927R>`_), including newly defined indices (e.g., `Cervantes & Vazdekis 2009 <https://ui.adsabs.harvard.edu/abs/2009MNRAS.392..691C>`_).


The LIS system is defined for 3 different spectral resolutions which are best suited for the following astrophysical cases:

* LIS-5.0AA: globular clusters
* LIS-8.4AA: low and intermediate-mass galaxies
* LIS-14.0AA: massive galaxies

Conversions to transform the data from the Lick/IDS system to LIS can be found in `Johansson, Thomas & Maraston (2010) <https://ui.adsabs.harvard.edu/abs/2010MNRAS.404.1639V>`_, which provides a discussion of indices
and the information content of them.

We provide a brief summary of the compiled Lick definitions provided with pyphot.

.. list-table:: Indices compiled in pyphot and brief information
    :widths: 15 38
    :header-rows: 1

    * - Index
      - Information content
    * - CN_1
      - CN band; nitrogen (and carbon) abundance; metallicity; giant-star contribution [1][2].
    * - CN_2
      - CN band; nitrogen/carbon; metallicity; giant-star contribution [1][2].
    * - Ca1_LB13
      - Ca II triplet (λ≈8498); metallicity; surface gravity; giant fraction/IMF; telluric/sky sensitive.
    * - Ca2_LB13
      - Ca II triplet (λ≈8542); metallicity; gravity; giant fraction/IMF.
    * - Ca3_LB13
      - Ca II triplet (λ≈8662); metallicity; gravity; giant fraction/IMF.
    * - Ca4227
      - Calcium line; [Ca/Fe]; metallicity; impacted by CN [1][2].
    * - Ca4455
      - Weak Ca/Fe blend; behaves largely like Fe; metallicity [1][2].
    * - CaH
      - CaH molecular band; cool-dwarf gravity/IMF sensitivity; temperature.
    * - CaHK_LB13
      - Ca H+K region; age, metallicity, [Ca/Fe]; continuum-shape sensitive.
    * - CaH_1
      - CaH band; dwarf/IMF sensitivity; temperature.
    * - CaH_2
      - CaH band; dwarf/IMF sensitivity; temperature.
    * - Fe4383
      - Strong Fe blend; iron abundance; metallicity [1][5].
    * - Fe4531
      - Fe+Ti blend; metallicity; some α sensitivity [1][2].
    * - Fe4668
      - Often denoted C4668; strong C and Fe feature; very Z-sensitive; carbon abundance [1][5].
    * - Fe5015
      - Fe-dominated (with Ti); metallicity [1][5].
    * - Fe5270
      - Fe line; iron abundance; used with Fe5335 to form [1][5].
    * - Fe5335
      - Fe line; iron abundance; pairs with Fe5270 [1][5].
    * - Fe5406
      - Fe feature; iron abundance [1][5].
    * - Fe5709
      - Weak Fe feature; metallicity [1].
    * - Fe5782
      - Weak Fe feature; metallicity [1].
    * - G4300
      - CH G-band; carbon abundance; metallicity; some age sensitivity [1].
    * - H_beta
      - Balmer Hβ; primary age indicator; mild Z dependence [1][5].
    * - Hbeta0
      - Optimized Hβo; age indicator with reduced metallicity sensitivity.
    * - HbetaEW
      - Hβ equivalent width; age; beware emission fill-in.
    * - HdeltaEW
      - Hδ equivalent width; sensitive to A-star light and recent SF.
    * - Hdelta_A
      - Broad HδA Balmer index; strong age sensitivity [3][5].
    * - Hdelta_F
      - Narrow HδF; age indicator with reduced Z cross-talk vs HδA [3][5].
    * - HgammaEW
      - Hγ equivalent width; age; sensitive to A-star light.
    * - Hgamma_A
      - Broad HγA; age-sensitive Balmer index [3][5].
    * - Hgamma_F
      - Narrow HγF; age with reduced Z sensitivity [3][5].
    * - Mg4780
      - Mg/Ti blend; metallicity; gravity/IMF sensitivity in some works.
    * - Mg_1
      - Molecular MgH + continuum break; metallicity and [Mg/Fe] (α) [1][5].
    * - Mg_2
      - As Mg1; widely used metallicity and [Mg/Fe] indicator [1][5].
    * - Mg_b
      - Mg triplet; [Mg/Fe] (α) and metallicity; some gravity sensitivity [1][5].
    * - NaI
      - Na I 8190 doublet; dwarf/IMF and [Na/Fe]; sensitive to tellurics/sky [6].
    * - NaI_F13
      - Na I 8190 definition (Ferreras+13 style); IMF/[Na/Fe] sensitive [6].
    * - NaI_LB13
      - Na I 8190 (La Barbera+13 style); IMF/[Na/Fe] sensitive [6].
    * - NaI_V12
      - Na I 8190 (Vazdekis+12 style); IMF/[Na/Fe] sensitive [6].
    * - Na_D
      - Na D doublet; [Na/Fe]; strong ISM/dust absorption component [1][6].
    * - OIIEW
      - [O II] λλ3726,3729 nebular emission; SF/AGN activity tracer.
    * - TiO2SDSS_LB13
      - TiO2 (SDSS variant); cool-dwarf/IMF sensitivity; metallicity, temperature [6].
    * - TiO3
      - TiO band; cool-dwarf/IMF sensitivity; temperature/metallicity.
    * - TiOCaH
      - Composite TiO+CaH gravity-sensitive band; IMF/cool-star content.
    * - TiO_1
      - Lick TiO1; cool-star fraction, IMF sensitivity; metallicity [1][5].
    * - TiO_2
      - Lick TiO2; cool-star/IMF sensitivity; metallicity [1][5].
    * - TiO_3
      - Red TiO band; cool-dwarf/IMF sensitivity; temperature.
    * - TiO_4
      - Red TiO band; cool-dwarf/IMF sensitivity; temperature.
    * - aTiO
      - TiO-based composite; cool-dwarf/IMF sensitivity; α-abundance/continuum effects possible.
    * - bTiO
      - Blue TiO/Mg region; metallicity and gravity/IMF sensitivity.

.. note::
    This table is compiled from various sources.

    1. Classic Lick/IDS definitions and usage: `Worthey et al. 1994 <https://ui.adsabs.harvard.edu/abs/1994ApJ...425..659W>`_
    2. Sensitivities to abundance changes: `Korn et al. 2005 <https://ui.adsabs.harvard.edu/abs/2005A%26A...438..685K>`_
    3. Balmer Hδ/Hγ A/F index definitions: `Worthey & Ottaviani 1997 <https://ui.adsabs.harvard.edu/abs/1997ApJS..111..377W>`_
    4. Context that `Vazdekis et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010MNRAS.404.1639V>`_ mentions Lick/IDS but does not tabulate them.
    5. Flux-calibrated Lick models and sensitivities: `Thomas, Maraston & Johansson 2011 <https://ui.adsabs.harvard.edu/abs/2011MNRAS.412.2183T>`_
    6. IMF/gravity-sensitive Na I, TiO variants in ETGs: `Spiniello et al. 2014 <https://ui.adsabs.harvard.edu/abs/2014MNRAS.438.1483S>`_


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

see: :func:`pyphot.licks.LickIndex.continuum_normalized_region_around_line`.


Finally any index is calculated by integrated the continuum normalized flux.
Some indices are given in `magnitudes`, and some in `equivalent width` units.


References
~~~~~~~~~~
* Cervantes & Vazdekis, 2009, MNRAS, 392, 691, "An optimized Hbeta index for disentangling stellar population ages" `2009MNRAS.392..691C <https://ui.adsabs.harvard.edu/abs/2009MNRAS.392..691C>`_

* Johansson, Thomas & Maraston, 2010, MNRAS, 406, 165, "Empirical calibrations of optical absorption-line indices based on the stellar library MILES" `2010MNRAS.406..165J <https://ui.adsabs.harvard.edu/abs/2010MNRAS.406..165J>`_

* Korn, Maraston & Thomas, 2005, A&A, 438, 685, "The sensitivity of Lick indices to abundance variations" `2005A&A...438..685K <https://ui.adsabs.harvard.edu/abs/2005A%26A...438..685K>`_


* Puzia et al. 2002, A&A, 395, 45, "Integrated spectroscopy of bulge globular clusters and fields. I. The data base and comparison of individual Lick indices in clusters and bulge" `2002A&A...395...45P <https://ui.adsabs.harvard.edu/abs/2002A&A...395...45P>`_

* Rose, 1985, AJ, 90, 1927, "Constraints on stellar populations in elliptical galaxies." `1985AJ.....90.1927R <https://ui.adsabs.harvard.edu/abs/1985AJ.....90.1927R>`_

* Spiniello et al 2014, MNRAS, 438, 1483, "The stellar IMF in early-type galaxies from a non-degenerate set of optical line indices" `2014MNRAS.438.1483S <https://ui.adsabs.harvard.edu/abs/2014MNRAS.438.1483S>`_

* Thomas, Maraston & Johansson, 2011, MNRAS, 412, 2183, "Flux-calibrated stellar population models of Lick absorption-line indices with variable element abundance ratios" `2011MNRAS.412.2183T <https://ui.adsabs.harvard.edu/abs/2011MNRAS.412.2183T>`_

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
