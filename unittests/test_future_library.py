"""Test the photometric library features"""

from re import I
import pytest
import pathlib
from typing import cast

from pyphot.future.libraries import HDF_Library, Ascii_Library, get_library
from pyphot import config
from pyphot.future.phot import Filter

lib = get_library()


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
