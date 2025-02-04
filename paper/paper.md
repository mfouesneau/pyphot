---
title: 'pyphot: A tool for computing photometry from spectra'
tags:
  - Python
  - astronomy
  - photometry
authors:
  - name: Morgan Fouesneau
    orcid: 0000-0001-9256-5516
    affiliation: 1 
    corresponding: true # (This is how to denote the corresponding author)
  - name: Ariane Lançon
    affiliation: 2
affiliations:
 - name: Max Planck Institute for Astronomy, Königstuhl 17, 69117 Heidelberg, Germany
   index: 1
 - name: Université de Strasbourg, CNRS, Observatoire astronomique de Strasbourg, UMR 7550, F-67000 Strasbourg, France
   index: 2
date: 4 February 2025
bibliography: paper.bib

---

# Summary

Pyphot is a Python package designed for performing synthetic photometry in astronomy. It provides a comprehensive set of tools for calculating fluxes and magnitudes of astronomical objects given their spectral energy distributions (SEDs) and filter transmission curves. Key features include a library of commonly used filter systems (e.g., Johnson-Cousins, SDSS), the ability to define custom filter systems, tools for handling SEDs, and integration with other astronomical Python packages. Pyphot also extends to calculations of spectral line indices from the Line Index System, which is commonly used for determining ages, mass distribution properties, and metallicities of unresolved (integrated light) stellar populations. Pyphot is targeted towards astronomers and astrophysicists working with photometric data.

# Statement of need

Synthetic photometry is a crucial task in astronomy, used to estimate the flux an object would have if observed through a given transmission curve. Prior to `pyphot`, a standardized, readily available tool for performing synthetic photometry in Python was lacking. Researchers often relied on ad-hoc scripts or software with limited flexibility and interoperability, hindering reproducibility and collaboration. Pyphot addresses this need by providing a user-friendly, well-documented, and versatile package for synthetic photometry. It simplifies the process of calculating magnitudes and fluxes, enables the use of standard and custom filter systems, and promotes reproducible research practices within the astronomical community. By offering a centralized resource for synthetic photometry, `pyphot` facilitates collaboration and ensures consistent results across different studies.

Photometry involves measuring the flux or intensity of light emitted by astronomical objects. It is a fundamental technique used in astronomy to study any object that emits light, such as stars, galaxies, and quasars. Photometry provides information about the brightness, color, and spectral energy distribution of these objects, which is essential for understanding their physical properties and evolutionary history.

The type of detector and the units of flux can vary depending on the wavelength range being considered. One of the challenges in computing photometry is ensuring that one performs the calculations in the correct units and handling conversions between different units (e.g., Vega to AB magnitude or vice-versa). The online documentation[^3] of `pyphot` provides the mathematical details and references

To address this challenge, the `pyphot` tool provides a comprehensive set of functions and utilities for computing photometry from spectra. It includes built-in support for unit handling and conversions based on either the implementations of physical units in the `Astropy` package [@astropy] (`astropy.units`) or in the `Pint` package[^2]. It also provides a database of pre-defined filter transmission curves, lick index definitions [@worthey1994], and interfaces with the [SVO filter profile service](http://svo2.cab.inta-csic.es/theory/fps/) [@svofps2020] to enable a broader range of filters.

In addition to calculating photometry, `pyphot` provides a set of tools to define and calculate lick indices. The Lick system of spectral line indices is one of the most commonly used methods for estimating ages and metallicities of unresolved (integrated light) stellar populations. The calibration of the Lick system is complicated because the original Lick spectra were not flux calibrated, so there are usually systematic effects due to differences in continuum shape. Proper calibration involves observing many original lick standard stars and deriving offsets to the standard system.

@Vazdekis2010 proposed a new Line Index System, hereafter LIS, with three new spectral resolutions at which to measure the Lick indices. Note that this new system should not be restricted to the Lick set of indices in a flux-calibrated system. In fact, LIS can be used for any index in the literature (e.g., for the @Rose1985 indices), including newly defined indices (e.g., @Cervantes2009). An additional challenge comes from the definition of the LIS for different spectral resolutions, initially designed for studies of globular clusters, low and intermediate-mass galaxies, and massive galaxies. `pyphot` deals with adapting the spectral resolution internally to match the definitions and provides a pre-defined compiled list of indices.

Several Python packages offer tools for synthetic photometry in astronomy. `synphot` [@synphot] is a powerful library for simulating photometric data and spectra, incorporating custom filters and spectra. `pysynphot` [@pysynphot], initially designed as an object-oriented successor to IRAF's `STSDAS.SYNPHOT`, is now out of support. Both `synphot` and `pysynphot` are closely tied to the Space Telescope Science Institute (STScI) ecosystem and define their own data format. A key difference between these and `pyphot` lies in their scope and integration with the broader astronomical Python ecosystem. While `synphot` and `pysynphot` are feature-rich for HST-specific analyses, `pyphot` offers a more streamlined and flexible approach, readily integrating with packages like Astropy and the SVO-profile service and focusing on ease of use for general synthetic photometry tasks. `pyphot` also prioritizes simplicity and integration into larger projects, while `synphot` offers more advanced features for spectral manipulation and analysis.

By using `pyphot`, researchers and students can perform photometric computations with confidence, knowing that the units are handled correctly and conversions are applied accurately. This ensures the reliability and reproducibility of scientific results in the field of astronomy. `pyphot` has supported 29 scientific publications[^1] since 2019.

[^1]: Nasa ADS [search for "pyphot"](https://ui.adsabs.harvard.edu/search/fq=%7B!type%3Daqp%20v%3D%24fq_database%7D&fq_database=(database%3Aastronomy%20OR%20database%3Aphysics)&q=ack%3A%22pyphot%22%20or%20pyphot&sort=date%20desc%2C%20bibcode%20desc&p_=0)

[^2]: https://pint.readthedocs.io/en/stable/

[^3]: [`pyphot` mathematical details on predicting photometry](https://mfouesneau.github.io/pyphot/photometry.html)

# Acknowledgements

We acknowledge contributions from Tim Morton, Cliff Johnson, and Karl Gordon during the genesis of this project.  Turning this project into an independent package was in part motivated by the 2016 NYC Gaia Sprint, hosted by the Center for Computational Astrophysics at the Simons Foundation in New York City. 

# References