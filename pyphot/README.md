# Future developments of pyphot

## New features

- add hint typing in the code (not tested with mypy; pyright ok)
- adapters to interface transparently with varied unit frameworks (e.g., astropy.units, pint)
  - `future.config.set_units_backend` to switch between different unit frameworks
  - default priority is `astropy`, `pint`, `ezunits[legacy]`
- a more robust and general unit ensuring decorator `future.units_adapter.enforce_default_units` which also updates the docstring accordingly
- vega flavor can now be switched globally `future.config.set_vega_flavor`
  - still available locally (since v.1.7.0)
- added specific unittests - coverage > 70%
- internal simpletable removed in favor of pandas
  - python.io for ascii, fits, hdf5adapted from https://github.com/mfouesneau/simpletable



```{warning} configuration update only affects newly created objects

Any configuration changes will only affect newly created objects.

If the unit backend is changed, you may have a mix of units definitions that are not compatible with each other.

If you change the vega flavor, you may have a mix of Vega definitions, which will not raise errors.

```

## Incoming

- io.hdf does not allow export yet.

- update documentation with new features/refactoring

- remove dependency on astropy for the svo module
  - astropy.io.votable is used for parsing VOTable files from the SVO service
