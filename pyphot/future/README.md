# Future developments of pyphot

## New features

- add hint typing in the code (not tested with mypy; pyright ok)
- adapters to interface transparently with varied unit frameworks (e.g., astropy.units, pint)
  - `future.config.set_units_backend` to switch between different unit frameworks
  - default priority is `astropy`, `pint`, `ezunits[legacy]`
- a more robust and general unit ensuring decorator `future.units_adapter.enforce_default_units` which also updates the docstring accordingly
- vega flavor can now be switched globally `future.config.set_vega_flavor`
- added specific unittests
  - test_units_adapters, test_units_adapters_decorator
  - test_future_filter, test_future_library
  - test_future_vega_sun
  - test_future_licks

```{warning} configuration update only affects newly created objects

Any configuration changes will only affect newly created objects.

If the unit backend is changed, you may have a mix of units definitions that are not compatible with each other.

If you change the vega flavor, you may have a mix of Vega definitions, which will not raise errors.

```

## Incoming

- missing test on `Ascii_Library`

- replace internal simpletable by https://github.com/mfouesneau/simpletable
  - better support of varied formats
  - simpler API based on Pandas
  - DictDataFrame interface probably best here (avoid pandas)

- remove dependency on astropy for the svo module
  - astropy.io.votable is used for parsing VOTable files from the SVO service
