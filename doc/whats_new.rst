Version History
===============

This page lists the recent changes in versions of PyPhot.

v2.0.0
------
[November 11, 2025]

Breaking changes
~~~~~~~~~~~~~~~~
- major refactoring of the codebase: remove code duplication due to the split astropy vs. pint
    - :code:`import pyphot.astropy as pyphot` replaced by :code:`import pyphot` (default backend is 'astropy')
    - set units with :code:`pyphot.config.set_units_backend('astropy')` or :code:`pyphot.config.set_units_backend('pint')` (default is 'astropy')
- adapters to interface transparently with varied unit frameworks (e.g., astropy.units, pint)
  - :func:`pyphot.config.set_units_backend` to switch between different unit frameworks
  - default priority is `astropy`, `pint`, `ezunits[legacy]`

Non-breaking changes
~~~~~~~~~~~~~~~~~~~~
- added type hints (incl. :mod:`pyphot.typing`)
   - hint typing checked with pyright, not mypy.
- a more robust and general unit ensuring decorator :func:`pyphot.units_adapter.enforce_default_units` which also updates the docstring accordingly
- vega flavor can now be switched globally
    - use :func:`pyphot.config.set_vega_flavor` (see :doc:`vega` documentation)
    - remains available at the filter definition (since v.1.7.0)
- internal simpletable removed in favor of pandas
    - :mod:`pyphot.io` for ascii, fits, hdf5, and votable adapted from https://github.com/mfouesneau/simpletable
- removed astropy dependency for SVO interface (still currently needed for fits)
- more comprehensive unit tests (coverage > 70%)
- updated documentation, :doc:`QuickStart`, and :doc:`examples`


v1.7.0
------
[July 2, 2025]
- Added Vega flavors to the Vega class, allowing users to select different Vega reference spectra for photometric calculations.

v1.6.0
------
[January 21, 2025]
- Moved installation to use `pyproject.toml`
- Updated documentation
- Updates for deprecated calls.

v1.5.0
------
[January 21, 2025]
- add support for Python 3.13

v1.4.6
------
[June 26, 2024]
- Dropping support for Python <= 3.8 (due to HDF5 modules).
- Minor updates for Scipy 1.14.0 and Numpy 2.0

For versions prior to v1.4.6, please refer to the `release notes`_.

.. _release notes: https://github.com/mfouesneau/pyphot/releases
