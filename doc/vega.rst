Details on the internal Vega reference spectra
==============================================

Vega (aka :math:`\alpha` Lyrae, HD 172167) serves as the fundamental calibration
standard for stellar photometry due to its exceptional brightness, and favorable
observational characteristics.  It is the fifth brightest star in the sky, after
Sirius, Canopus, Alpha Centauri and Arcturus and the second brightest star in
the northern celestial hemisphere.
As a nearby A0V star located approximately 25 light-years from Earth, Vega was
historically chosen as the primary standard because of its position near the
north celestial pole, making it easily accessible to northern hemisphere
observatories year-round (`Johnson & Morgan 1953
<https://ui.adsabs.harvard.edu/abs/1953ApJ...117..313J>`_). 

Vega is a Delta Scuti variable star (dwarf Cepheids), one whose variations in
luminosity result from both radial and non-radial pulsations of its surface. 
Some parts of its surface contract while others simultaneously expand
(non-radial pulsations), and the star also contracts and expands by changing its
radius to maintain its spherical shape (radial pulsations).

Delta Scuti variables are commonly used as standard candles to establish
distances because of their relatively flat spectral energy distribution across
optical wavelengths and "stability" over decades make it an ideal
reference for establishing magnitude zero-points across multiple filter systems
(`Oke & Gunn 1983 <https://ui.adsabs.harvard.edu/abs/1983ApJ...266..713O>`_).

However, the fact that Vega is actually a variable star with small but
measurable brightness variations has led the astronomical community to adopt
more stable references. The International Astronomical Union (IAU) and major
observatories now commonly use synthetic standards based on theoretical stellar
atmosphere models or carefully selected ensembles of stable stars. These
synthetic standards provide consistent, reproducible reference points that are
not subject to the intrinsic variability observed in individual stars like Vega,
ensuring long-term stability in photometric calibrations across different
observatories and epochs.

Modern space-based observations have refined Vega's spectral energy distribution
to unprecedented precision, with the CALSPEC database providing the definitive
reference spectrum used by HST and other major observatories (`Bohlin et al.
2014 <https://ui.adsabs.harvard.edu/abs/2014PASP..126..711B>`_). Despite the
discovery of Vega's infrared excess due to a circumstellar debris disk and minor
photometric variability, it remains the cornerstone of photometric systems, with
the AB magnitude system providing an alternative that maintains compatibility
while addressing some of Vega's limitations (`Fukugita et al. 1996
<https://ui.adsabs.harvard.edu/abs/1996AJ....111.1748F>`_).

Vega Flavors in Pyphot
----------------------

Since version 1.7.0, Pyphot includes a set of Vega flavors one can use transparently as photometric standards. We summarize the available flavors below:

* `alpha_lyr_004` from `Bohlin, Colina, & Finley (1995) <https://ui.adsabs.harvard.edu/abs/1995AJ....110.1316B>`_. It corrsponds to pure hydrogen white dwarf models are the absolute flux standards.
* `alpha_lyr_005` from Colina, Bohlin & Castelli (1996), Instrument Science Report (`OSG-CAL-96-01 <https://www.stsci.edu/instruments/observatory/PDF/scs8.rev.pdf>`_). It corresponds to the combination of ultraviolet covered by average IUE spectrum; Optical up to 1.05 microns covered by Hayes (1985) average spectrum (IAU Symp 111, p 225); near-infrared covered by ATLAS12 based model computed by Dr. Castelli, rebinned to 25A, and normalized to Hayes (1985) Johnson V flux. V filter as in Buser & Kurucz 1979, AA 70, 555. (All wavelengths are vacuum.)
* `alpha_stis_003` (legacy) from `Bohlin (2007) <https://ui.adsabs.harvard.edu/abs/2007ASPC..364..315B>`_ a special model from Kurucz 9400K Vega spectrum T/g=9400/3.9 [M/H]=-0.5 fit against STIS data with microturbulence of 0. km/s, and where V = 0.023 mag.
* `alpha_mod_001` from `Bohlin (2014) <https://ui.adsabs.harvard.edu/abs/2014AJ....147..127B>`_ Special Model from Kurucz 9400K Vega spectrum T/g=9400/3.9 [M/H]=-0.5 at R=500 which reconciles visible and IR absolute flux, and where Vega Flux(5556A)=3.44e-9
* `alpha_mod_002` is identical to `alpha_mod_001` but normalized to STIS flux at 5545-5570A, corresponding to a scaled flux by $0.994242$.
* `alpha_mod_003` is a different temperature, logg using Kurucz 9550K Vega spectrum T/g=9550/3.95 [M/H]=-0.5 at R=500 and normalized to Flux(5556A) = 3.44e-9 (over 5545-5570A).
* `alpha_mod_004` from `Bohlin, Hubeny, Rauch (2020) <https://ui.adsabs.harvard.edu/abs/2020AJ....160...21B>`_ is similar to `alpha_mod_003` but with a Flux(5556A)=3.47e-9. In the details, it also includes additional lines.
* `alpha_stis_011` is a special model from Bohlin which is a composite flux of a special Kurucz 9550K model from 900-1152A (Kurucz 2003), IUE data from 1152-1675A, STIS CCD fluxes from 1675-10200A (Bohlin & Gilliland 2004a), and the 9550K model longward of 10200A. It differs significantly in the uv-optical range by accounting for the new TMAP AND TLUSTY WD NLTE models (`Bohlin, Hubeny, Rauch (2020) <https://ui.adsabs.harvard.edu/abs/2020AJ....160...21B>`_) but corresponds to `alpha_mod_004` above 1 micron.

Impact of Vega Flavors on zeropoints provided by Pyphot:

* Changes in temperature will introduce a wavelength dependent shift in the vega zero points.
* Changes in logg, metallicity, turbulence will introduce non trivial variations in the vega zero points.

.. note:: 

    By default, Pyphot uses the `alpha_stis_003` flavor as the Vega standard, which may be updated to `alpha_stis_011` in the future.


TODO: 
=====

* Need to set a vega flavor keyword in the passband class (calculations call `with Vega() as v: ...`)
* Add text/example on how to change the flavor of Vega used in the calculations.
* make a plots that shows major differences and their scale: ~percent and below for minor updates, up to 10% in the UV for  `alpha_mod_004` vs `alpha_stis_011`.


References
~~~~~~~~~~

.. custom ADS format: %l, %Y, %q, %V, %p, "%T", `%R <%u>`_

* Bohlin, R. C., Gordon, K. D., & Tremblay, P.-E., 2014, PASP, 126, 711, "Techniques and Review of Absolute Flux Calibration from the Ultraviolet to the Mid-Infrared", `2014PASP..126..711B <https://ui.adsabs.harvard.edu/abs/2014PASP..126..711B>`_
* Colina, Bohlin & Castelli 1996, Instrument Science Report, "Absolute Flux Calibrated Spectrum of Vega" `OSG-CAL-96-01 <https://www.stsci.edu/instruments/observatory/PDF/scs8.rev.pdf>`_
* Fukugita, M., Ichikawa, T., Gunn, J. E., Doi, M., Shimasaku, K., & Schneider, D. P., 1996, AJ, 111, 1748, "The Sloan Digital Sky Survey Photometric System", `1996AJ....111.1748F <https://ui.adsabs.harvard.edu/abs/1996AJ....111.1748F>`_
* Johnson, H. L. \& Morgan, W. W. 1953, ApJ, 117, 313, "Fundamental stellar photometry for standards of spectral type on the Revised System of the Yerkes Spectral Atlas." `1953ApJ...117..313J <https://ui.adsabs.harvard.edu/abs/1953ApJ...117..313J>`_; 
* Oke, J. B. and Gunn, J. E., 1983, ApJ, 266, 713, "Secondary standard stars for absolute spectrophotometry.”  `1983ApJ...266..713O <https://ui.adsabs.harvard.edu/abs/1983ApJ...266..713O>`_;
* Bohlin, R. C., Colina, L., & Finley, D. S., 1995, AJ, 110, 1316, "White Dwarf Standard Stars: G191-B2B, GD 71, GD 153, HZ 43", `1995AJ....110.1316B <https://ui.adsabs.harvard.edu/abs/1995AJ....110.1316B>`_
* Bohlin, R. C., 2007, ASPC, 364, 315, "HST Stellar Standards with 1% Accuracy in Absolute Flux", `2007ASPC..364..315B <https://ui.adsabs.harvard.edu/abs/2007ASPC..364..315B>`_
* Bohlin, R. C., 2014, AJ, 147, 127, "Hubble Space Telescope CALSPEC Flux Standards: Sirius (and Vega)", `2014AJ....147..127B <https://ui.adsabs.harvard.edu/abs/2014AJ....147..127B>`_
* Bohlin, R. C., Hubeny, I., & Rauch, T., 2020, AJ, 160, 21, "New Grids of Pure-hydrogen White Dwarf NLTE Model Atmospheres and the HST/STIS Flux Calibration", `2020AJ....160...21B <https://ui.adsabs.harvard.edu/abs/2020AJ....160...21B>`_


