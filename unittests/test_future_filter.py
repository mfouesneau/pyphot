"""Test Filter definitions"""

import pytest
from typing import cast
import numpy as np

from pyphot.future.libraries import get_library
from pyphot.future import config
from pyphot.future.phot import _drop_units, Filter
from pyphot.future.vega import Vega

config.set_vega_flavor("stis_003")  # make sure we compare apples to apples
# assume internal default library for this
LIBRARY = get_library()


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


@pytest.mark.parametrize("name", list(reference_info.keys()))
def test_filter_attributes(name):
    """Test filter attributes"""
    pb = cast(Filter, LIBRARY[name])
    pb.info()  # checking for any issues in the calculations
    refs = reference_info[name]
    assert pb.name == name
    for key, value in refs.items():
        other = _drop_units(getattr(pb, key))
        if isinstance(value, str):
            assert value == other
        else:
            assert abs(value - other) < 1e-6


@pytest.mark.parametrize("pbname", ["GROUND_JOHNSON_V", "HST_WFC3_F110W"])
def test_convert_mags(pbname):
    # get the internal default library of passbands filters
    lib = get_library()
    vega = Vega()
    f = cast(Filter, lib[pbname])
    # compute the integrated flux through the filter f
    # note that it work on many spectra at once
    fluxes = f.get_flux(vega.wavelength, vega.flux, axis=-1)
    # Note that fluxes is now with units of erg/s/cm2/AA
    # pyphot gives Vega in flam and can convert between flux density units.

    # convert to vega magnitudes
    mags = -2.5 * np.log10(fluxes.value) - f.Vega_zero_mag
    # print("Vega magnitude of Vega in {0:s} is : {1:f} mag".format(f.name, mags))
    assert abs(mags) < 1e-5, "Vega magnitude is not within expected range"
    flux = fluxes / f.Vega_zero_flux
    # print("Vega flux of Vega in {0:s} is : {1:f}".format(f.name, flux))
    assert abs(flux - 1.0) < 1e-5, "Vega flux is not within expected range"


def test_vega_reference_non_overlapping():
    """Test the code is robust to the lack of Vega reference

    This happens in the far IR (e.g., Herschel SPIRE), where calibration has to
    replace the Vega reference with a synthetic atmosphere.

    """
    pbname = "HERSCHEL_SPIRE_PLW"
    lib = get_library()
    f = cast(Filter, lib[pbname])

    import math

    assert math.isnan(f.leff.value)
    assert f.Vega_zero_flux == 0.0 * config.units.U("flam")


def generate_ascii_filter_test_data():
    from textwrap import dedent

    data = """
    # WAVELENGTH_UNIT	nm
    # NAME           	TESTDATA
    # DETECTOR       	photon
    #
    ## WAVELENGTH	None	nm
    ## THROUGHPUT	None	filter throughput definition
    #
    WAVELENGTH,THROUGHPUT
    300.000000,0.000000
    309.000000,0.000000
    399.000000,0.000000
    400.000000,1.000000
    500.000000,1.000000
    501.000000,0.000000
    600.000000,0.000000
    """
    return dedent(data[1:])


def test_filter_from_ascii_data():
    from io import StringIO

    datastr = StringIO(generate_ascii_filter_test_data())
    dtype = "csv"
    kwargs = {}

    r = Filter.from_ascii(datastr, dtype=dtype, **kwargs)
    assert r is not None
    assert r.name == "TESTDATA"
    assert r.dtype == "photon"
    assert r.wavelength_unit == "nm"
    assert abs(r.cl.value - 450.0) < 1e-6
    assert abs(r.lmin.value - 400.0) < 1e-6
    assert abs(r.lmax.value - 500.0) < 1e-6
