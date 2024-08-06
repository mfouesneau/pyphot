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
 - name: Observatoire astronomique de Strasbourg, Université de Strasbourg, CNRS, UMR 7550, 67000 Strasbourg, France
   index: 2
date: 5 August 2024
bibliography: paper.bib

---

# Summary

Photometry involves measuring the flux or intensity of light emitted by
astronomical objects. It is a fundamental technique used in astronomy to study
any object that emits light, such as stars, galaxies, and quasars. Photometry
provides information about the brightness, color, and spectral energy
distribution of these objects, which is essential for understanding their
physical properties and evolutionary history.

The type of detector and the units of flux can vary depending on the wavelength
range being considered.  Converting between these units requires knowledge of
the spectral response of the instrument used for observation. One of the
challenges in computing photometry is ensuring that one performs the
calculations in the correct units and handling conversions between different
units (e.g. Vega to AB magnitude or vice-versa). 

To address this challenge, the `pyphot` tool provides a comprehensive set of
functions and utilities for computing photometry from spectra. It includes
built-in support for unit handling and conversions based on either the
implementations of physical units in the `Astropy` package [@astropy]
(`astropy.units`) or in the `Pint` package[^1]. It also provides a database of
pre-defined filter transmission curves, lick index definitions
[@worthey1994] and interfaces with the [SVO filter profile
service](http://svo2.cab.inta-csic.es/theory/fps/) [@svofps2020] to
enable a broader range of filters.

By using `pyphot`, researchers and students can perform photometric computations
with confidence, knowing that the units are handled correctly and conversions
are applied accurately. This ensures the reliability and reproducibility of
scientific results in the field of astronomy. The source code for Gala has been archived to Zenodo with
the linked DOI: @zenodopyphot

`pyphot` has supported 25 of scientific publications[^2] since 2019.

[^1]: https://pint.readthedocs.io/en/stable/

[^2]: Nasa ADS [search for "pyphot"](https://ui.adsabs.harvard.edu/search/fq=%7B!type%3Daqp%20v%3D%24fq_database%7D&fq_database=(database%3Aastronomy%20OR%20database%3Aphysics)&q=ack%3A%22pyphot%22%20or%20pyphot&sort=date%20desc%2C%20bibcode%20desc&p_=0)


# Acknowledgements

We acknowledge contributions from Tim Morton, Cliff Johnson, and Karl Gordon
during the genesis of this project.  This project was in part motivated by the
2016 NYC Gaia Sprint, hosted by the Center for Computational Astrophysics at the
Simons Foundation in New York City.

# References