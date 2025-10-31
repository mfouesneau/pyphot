"""testing the vega/sun interfaces"""

import pytest

from typing import cast

from pyphot.future.libraries import get_library
from pyphot.future import vega, sun, config
from pyphot.future.phot import Filter


@pytest.mark.parametrize("flavor", list(vega._default_vega.keys()))
def test_instanciate_vega(flavor):
    """Checks that all registered flavors can be instantiated"""
    _ = vega.Vega(flavor=flavor)


@pytest.mark.parametrize("flavor", list(sun._default_sun.keys()))
def test_instanciate_sun(flavor):
    """Checks that all registered flavors can be instantiated"""
    _ = sun.Sun(flavor=flavor)


def test_sun_magnitudes(self):
    import numpy as np

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
        print(
            "{0:12s} {1:0.5e} {2:+3.4f}".format(
                name, flux.magnitude, -2.5 * np.log10(flux.magnitude) - vegamag
            )
        )
        sun_vega_mag = -2.5 * np.log10(flux.magnitude) - vegamag

        assert abs(sun_vega_mag - expect_mag) < 1e-2
