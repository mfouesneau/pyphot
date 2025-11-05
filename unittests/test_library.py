"""Test the photometric library features"""

import pytest
import pathlib
from typing import cast

from pyphot.libraries import HDF_Library, Ascii_Library, get_library
from pyphot import config
from pyphot.phot import _drop_units, Filter

lib = get_library()

reference_info = {
    "GROUND_JOHNSON_V": {
        "dtype": "energy",
        "wavelength_unit": "AA",
        "cl": 5504.666021,
        "lpivot": 5469.853351,
        "leff": 5438.689749,
        "lphot": 5460.082104,
        "lmin": 4766.670000,
        "lmax": 6899.950000,
        "norm": 870.650149,
        "width": 870.651020,
        "fwhm": 827.370000,
        "Vega_zero_mag": 21.099102,
        "AB_zero_mag": 21.097827,
        "ST_zero_mag": 21.100000,
    }
}


def test_hdf_library():
    """Test HDF library"""
    try:
        import tables
    except ImportError:
        pytest.skip("tables not installed")
    lib = get_library()
    assert isinstance(lib, HDF_Library)
    assert lib.source is not None
    assert config.__default_passband_lib__ == pathlib.Path(lib.source)
    with tables.File(lib.source, "r") as f:
        assert len(lib) == len(f.root.filters._v_children)


@pytest.mark.parametrize("name", lib.content)
def test_integrity_in_internal_library(name):
    """Make sure all internal filters can be loaded"""
    # get the internal default library of passbands filters
    fk = cast(Filter, lib[name])

    fmt = "{0:s},{1:s},{2:s},{3:.3f},{4:.3f},{5:.3f},{6:.5f},{7:.5g},{8:.5g},{9:.5f},{10:.5g},{11:.5g},{12:.5f},{13:.5g},{14:.5g}\n"
    # check main attributes exist
    rec = (
        fk.name,
        fk.dtype,
        fk.wavelength_unit,
        fk.cl.value,
        fk.lpivot.value,
        fk.leff.value,
        fk.Vega_zero_mag,
        fk.Vega_zero_flux.value,
        fk.Vega_zero_Jy.value,
        fk.AB_zero_mag,
        fk.AB_zero_flux.value,
        fk.AB_zero_Jy.value,
        fk.ST_zero_mag,
        fk.ST_zero_flux.value,
        fk.ST_zero_Jy.value,
    )
    print(fmt.format(*rec))  # checks value types


def test_ascii_library():
    """Test ascii library using internal files"""
    where = pathlib.Path(str(config.__default_passband_lib__)).parent
    where = str(where / "ascii_sources")
    lib = Ascii_Library(where)
    assert isinstance(lib, Ascii_Library)
    assert lib.source is not None
    assert len(lib) > 0
    which = lib.find("johnson_v")
    assert "ground_johnson_v" in which[0]
    pb = lib.load_filters(which)[0]

    refs = reference_info["GROUND_JOHNSON_V"]
    for key, value in refs.items():
        other = _drop_units(getattr(pb, key))
        if isinstance(value, str):
            assert value == other
        else:
            assert abs(value - other) < 1e-6
