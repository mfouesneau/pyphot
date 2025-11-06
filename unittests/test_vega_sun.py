"""testing the vega/sun interfaces"""

import pytest

from typing import cast
import numpy as np

from pyphot.libraries import get_library
from pyphot import vega, sun, config
from pyphot.phot import Filter
from pyphot.unit_adapters import backends


@pytest.mark.parametrize("flavor", list(vega._default_vega.keys()))
def test_instanciate_vega(flavor):
    """Checks that all registered flavors can be instantiated"""
    _ = vega.Vega(flavor=flavor)


@pytest.mark.parametrize("flavor", list(sun._default_sun.keys()))
def test_instanciate_sun(flavor):
    """Checks that all registered flavors can be instantiated"""
    _ = sun.Sun(flavor=flavor)


@pytest.mark.parametrize("AdapterName", list(backends.keys()))
def test_sun_magnitudes(AdapterName):
    if AdapterName is not None:
        try:
            config.set_units_backend(AdapterName)
        except ValueError:
            pytest.skip(f"Skipping test for {AdapterName} adapter")
    sun_obs = sun.Sun(flavor="observed")
    sun_th = sun.Sun()  # default is theoric spectrum
    sun_th_10pc = sun.Sun(distance=10 * config.units.U("pc"))

    lib = get_library()

    f = cast(Filter, lib["GROUND_JOHNSON_V"])
    expectations = -26.76, -26.76, +4.81

    for name, sun_flavor, expect_mag in zip(
        ("observed", "theoretical", "th. 10pc"),
        (sun_obs, sun_th, sun_th_10pc),
        expectations,
    ):
        flux = f.get_flux(sun_flavor.wavelength, sun_flavor.flux)
        vegamag = f.Vega_zero_mag
        sun_vega_mag = -2.5 * np.log10(flux.value) - vegamag
        assert abs(sun_vega_mag - expect_mag) < 1e-2
