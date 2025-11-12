"""testing the vega/sun interfaces"""

import pytest
import pandas as pd
from unittest.mock import patch

from typing import cast
import numpy as np

from pyphot.libraries import get_library
from pyphot import sun, config
from pyphot.phot import Filter
from pyphot.unit_adapters import backends


@pytest.mark.parametrize("flavor", list(sun._default_sun.keys()))
def test_instanciate_sun(flavor):
    """Checks that all registered flavors can be instantiated"""
    s = sun.Sun(flavor=flavor)
    assert s.data is not None
    assert s._data is not None  # should not be after the previous call
    assert s.data is not None


@pytest.mark.parametrize("AdapterName", list(backends.keys()))
def test_sun_magnitudes(AdapterName):
    if AdapterName is not None:
        try:
            config.set_units_backend(AdapterName)
        except ImportError:
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


def test_get_library_default_distance():
    """Test the get_library_default_distance function"""
    distance = sun.get_library_default_distance()

    # Should return 1 AU
    assert distance.value == 1.0
    # Should have appropriate units (AU)
    assert str(distance.unit).lower() in ["au", "astronomical_unit"]


def test_sun_constructor_defaults():
    """Test Sun constructor with default parameters"""
    sun_obj = sun.Sun()

    assert sun_obj.flavor == "theoretical"
    assert sun_obj.source == sun._default_sun["theoretical"]
    assert sun_obj._data is None
    assert sun_obj.units is None
    assert sun_obj.distance_conversion == 1.0


def test_sun_constructor_with_flavor():
    """Test Sun constructor with different flavors"""
    # Test theoretical flavor
    sun_theo = sun.Sun(flavor="theoretical")
    assert sun_theo.flavor == "theoretical"
    assert sun_theo.source == sun._default_sun["theoretical"]

    # Test observed flavor
    sun_obs = sun.Sun(flavor="observed")
    assert sun_obs.flavor == "observed"
    assert sun_obs.source == sun._default_sun["observed"]


def test_sun_constructor_with_custom_source():
    """Test Sun constructor with custom source file"""
    custom_source = "/path/to/custom_sun.fits"
    sun_obj = sun.Sun(source=custom_source, flavor="observed")

    assert sun_obj.source == custom_source
    assert sun_obj.flavor == "observed"


def test_sun_constructor_with_distance():
    """Test Sun constructor with custom distance"""
    custom_distance = 10.0 * config.units.U("pc")
    sun_obj = sun.Sun(distance=custom_distance)

    assert sun_obj.distance == custom_distance
    # Distance conversion should be calculated
    assert sun_obj.distance_conversion != 1.0


def test_sun_constructor_invalid_flavor():
    """Test Sun constructor with invalid flavor"""
    with pytest.raises(RuntimeError, match="Flavor must be in"):
        sun.Sun(flavor="invalid_flavor")  # type: ignore / on purpose


def test_sun_constructor_no_source_no_flavor():
    """Test Sun constructor with neither source nor flavor"""
    with pytest.raises(
        RuntimeError, match="Either `source` or `flavor` must be provided"
    ):
        sun.Sun(source=None, flavor=None)  # type: ignore / on purpose


def test_set_source_flavor_with_source():
    """Test _set_source_flavor method with custom source"""
    sun_obj = sun.Sun()
    custom_source = "/custom/path.fits"

    sun_obj._set_source_flavor(source=custom_source, flavor="observed")

    assert sun_obj.source == custom_source
    assert sun_obj.flavor == "observed"
    assert sun_obj._data is None
    assert sun_obj.units is None


def test_set_source_flavor_with_flavor():
    """Test _set_source_flavor method with flavor only"""
    sun_obj = sun.Sun(flavor="observed")  # Start with observed

    # Change to theoretical
    sun_obj._set_source_flavor(source=None, flavor="theoretical")

    assert sun_obj.flavor == "theoretical"
    assert sun_obj.source == sun._default_sun["theoretical"]


def test_set_source_flavor_invalid_flavor():
    """Test _set_source_flavor method with invalid flavor"""
    sun_obj = sun.Sun()

    with pytest.raises(RuntimeError, match="Flavor must be in"):
        sun_obj._set_source_flavor(source=None, flavor="invalid")  # type: ignore / on purpose


def test_set_source_flavor_no_parameters():
    """Test _set_source_flavor method with no parameters"""
    sun_obj = sun.Sun()

    with pytest.raises(
        RuntimeError, match="Either `source` or `flavor` must be provided"
    ):
        sun_obj._set_source_flavor(source=None, flavor=None)


@patch("pyphot.sun.io.from_file")
def test_readfile_method(mock_from_file):
    """Test _readfile method"""
    # Mock data
    mock_data = pd.DataFrame(
        {"WAVELENGTH": [4000, 5000, 6000], "FLUX": [1.0, 1.5, 1.2]}
    )
    mock_data.attrs = {
        "WAVELENGTH_UNIT": "ANGSTROMS",
        "FLUX_UNIT": "FLAM",
    }
    mock_hdr = {}
    mock_from_file.return_value = (mock_data, mock_hdr)

    sun_obj = sun.Sun()
    data, w_unit, f_unit = sun_obj._readfile()

    assert isinstance(data, pd.DataFrame)
    assert w_unit == "ANGSTROMS"
    assert f_unit == "FLAM"
    assert sun_obj._data is not None
    assert sun_obj.units == ("ANGSTROMS", "FLAM")


@patch("pyphot.sun.io.from_file")
def test_readfile_method_bytes_units(mock_from_file):
    """Test _readfile method with bytes units (fallback case)"""
    # Mock data with bytes units
    mock_data = pd.DataFrame(
        {"WAVELENGTH": [4000, 5000, 6000], "FLUX": [1.0, 1.5, 1.2]}
    )
    mock_data.attrs = {
        "WAVELENGTH_UNIT": b"ANGSTROMS=something",
        "FLUX_UNIT": b"FLAM=something",
    }
    mock_hdr = {}
    mock_from_file.return_value = (mock_data, mock_hdr)

    sun_obj = sun.Sun()
    data, w_unit, f_unit = sun_obj._readfile()

    assert w_unit == "ANGSTROMS"
    assert f_unit == "FLAM"


@patch("pyphot.sun.io.from_file")
def test_readfile_cached(mock_from_file):
    """Test that _readfile uses cached data on subsequent calls"""
    # Mock data
    mock_data = pd.DataFrame(
        {"WAVELENGTH": [4000, 5000, 6000], "FLUX": [1.0, 1.5, 1.2]}
    )
    mock_data.attrs = {
        "WAVELENGTH_UNIT": "ANGSTROMS",
        "FLUX_UNIT": "FLAM",
    }
    mock_hdr = {}
    mock_from_file.return_value = (mock_data, mock_hdr)

    sun_obj = sun.Sun()

    # First call
    data1, w_unit1, f_unit1 = sun_obj._readfile()
    # Second call
    data2, w_unit2, f_unit2 = sun_obj._readfile()

    # Should only call the file reader once
    assert mock_from_file.call_count == 1
    # Results should be the same
    assert data1 is data2
    assert w_unit1 == w_unit2
    assert f_unit1 == f_unit2


def test_context_manager():
    """Test Sun as context manager"""
    with patch("pyphot.sun.io.from_file") as mock_from_file:
        # Mock data
        mock_data = pd.DataFrame(
            {"WAVELENGTH": [4000, 5000, 6000], "FLUX": [1.0, 1.5, 1.2]}
        )
        mock_data.attrs = {
            "WAVELENGTH_UNIT": "ANGSTROMS",
            "FLUX_UNIT": "FLAM",
        }
        mock_hdr = {}
        mock_from_file.return_value = (mock_data, mock_hdr)

        sun_obj = sun.Sun()

        with sun_obj as s:
            assert s is sun_obj
            assert s._data is not None

        # Should have read the file during context entry
        assert mock_from_file.call_count == 1


@patch("pyphot.sun.io.from_file")
def test_data_property(mock_from_file):
    """Test data property"""
    # Mock data
    mock_data = pd.DataFrame(
        {"WAVELENGTH": [4000, 5000, 6000], "FLUX": [1.0, 1.5, 1.2]}
    )
    mock_data.attrs = {
        "WAVELENGTH_UNIT": "ANGSTROMS",
        "FLUX_UNIT": "FLAM",
    }
    mock_hdr = {}
    mock_from_file.return_value = (mock_data, mock_hdr)

    sun_obj = sun.Sun()

    # Access data property
    data = sun_obj.data

    assert isinstance(data, pd.DataFrame)
    assert len(data) == 3
    assert "WAVELENGTH" in data.columns
    assert "FLUX" in data.columns


@patch("pyphot.sun.io.from_file")
def test_wavelength_property(mock_from_file):
    """Test wavelength property"""
    # Mock data
    mock_data = pd.DataFrame(
        {"WAVELENGTH": [4000, 5000, 6000], "FLUX": [1.0, 1.5, 1.2]}
    )
    mock_data.attrs = {
        "WAVELENGTH_UNIT": "ANGSTROMS",
        "FLUX_UNIT": "FLAM",
    }
    mock_hdr = {}
    mock_from_file.return_value = (mock_data, mock_hdr)

    sun_obj = sun.Sun()

    # Access wavelength property
    wavelength = sun_obj.wavelength

    # Should have units
    assert hasattr(wavelength, "unit")
    assert hasattr(wavelength, "value")
    # Should have the right values
    np.testing.assert_array_equal(wavelength.value, [4000, 5000, 6000])


@patch("pyphot.sun.io.from_file")
def test_flux_property(mock_from_file):
    """Test flux property"""
    # Mock data
    mock_data = pd.DataFrame(
        {"WAVELENGTH": [4000, 5000, 6000], "FLUX": [1.0, 1.5, 1.2]}
    )
    mock_data.attrs = {
        "WAVELENGTH_UNIT": "ANGSTROMS",
        "FLUX_UNIT": "FLAM",
    }
    mock_hdr = {}
    mock_from_file.return_value = (mock_data, mock_hdr)

    sun_obj = sun.Sun()

    # Access flux property
    flux = sun_obj.flux

    # Should have units
    assert hasattr(flux, "unit")
    assert hasattr(flux, "value")
    # Should have the right values (multiplied by distance_conversion)
    expected_flux = np.array([1.0, 1.5, 1.2]) * sun_obj.distance_conversion
    np.testing.assert_array_equal(flux.value, expected_flux)


@patch("pyphot.sun.io.from_file")
def test_flux_property_with_distance_conversion(mock_from_file):
    """Test flux property with distance conversion"""
    # Mock data
    mock_data = pd.DataFrame(
        {"WAVELENGTH": [4000, 5000, 6000], "FLUX": [2.0, 3.0, 2.4]}
    )
    mock_data.attrs = {
        "WAVELENGTH_UNIT": "ANGSTROMS",
        "FLUX_UNIT": "FLAM",
    }
    mock_hdr = {}
    mock_from_file.return_value = (mock_data, mock_hdr)

    # Create Sun with custom distance (should trigger distance conversion)
    custom_distance = 10.0 * config.units.U("pc")
    sun_obj = sun.Sun(distance=custom_distance)

    # Access flux property
    flux = sun_obj.flux

    # Should apply distance conversion
    expected_flux = np.array([2.0, 3.0, 2.4]) * sun_obj.distance_conversion
    np.testing.assert_array_equal(flux.value, expected_flux)
    # Distance conversion should not be 1.0 for custom distance
    assert sun_obj.distance_conversion != 1.0


def test_distance_conversion_calculation():
    """Test distance conversion calculation"""
    # Test with 10 parsecs distance
    distance_10pc = 10.0 * config.units.U("pc")
    sun_obj = sun.Sun(distance=distance_10pc)

    # Distance conversion should be (1 AU / 10 pc)^2
    default_distance = sun.get_library_default_distance()
    expected_conversion = float((default_distance / distance_10pc).to("") ** 2)

    assert abs(sun_obj.distance_conversion - expected_conversion) < 1e-10


def test_sun_default_sources():
    """Test that default sources dictionary is properly defined"""
    assert "observed" in sun._default_sun
    assert "theoretical" in sun._default_sun

    # Check that paths contain expected keywords
    assert "sun_reference_stis" in sun._default_sun["observed"]
    assert "sun_kurucz93" in sun._default_sun["theoretical"]


def test_sun_all_exports():
    """Test that __all__ contains expected exports"""
    assert sun.__all__ == ["Sun"]
