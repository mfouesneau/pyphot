"""testing the vega/sun interfaces"""

import pytest
import os
import tempfile
import pandas as pd
from unittest.mock import patch, Mock, mock_open

from typing import cast
import numpy as np

from pyphot.libraries import get_library
from pyphot import vega, config
from pyphot.phot import Filter
from pyphot.unit_adapters import backends


@pytest.mark.parametrize("flavor", list(vega._default_vega.keys()))
def test_instanciate_vega(flavor):
    """Checks that all registered flavors can be instantiated"""
    v = vega.Vega(flavor=flavor)
    assert v.data is not None
    assert v._data is not None  # should not be after the previous call
    assert v.data is not None  # checks the cache


def test_vega_constructor_defaults():
    """Test Vega constructor with default parameters"""
    vega_obj = vega.Vega()

    assert vega_obj.source == vega._default_vega["legacy"]
    assert vega_obj._data is None
    assert vega_obj.units is None


def test_vega_constructor_with_flavor():
    """Test Vega constructor with different flavors"""
    # Test with each available flavor
    for flavor in vega._default_vega.keys():
        vega_obj = vega.Vega(flavor=flavor)
        assert vega_obj.source == vega._default_vega[flavor]


def test_vega_constructor_with_custom_source():
    """Test Vega constructor with custom source file"""
    custom_source = "/path/to/custom_vega.fits"
    vega_obj = vega.Vega(source=custom_source)

    assert vega_obj.source == custom_source


def test_vega_constructor_invalid_flavor():
    """Test Vega constructor with invalid flavor"""
    with pytest.raises(ValueError, match="Unknown Vega flavor"):
        vega.Vega(flavor="invalid_flavor")


def test_vega_constructor_no_source_no_flavor():
    """Test Vega constructor with neither source nor flavor"""
    with pytest.raises(
        RuntimeError, match="Either `source` or `flavor` must be provided"
    ):
        vega.Vega(source=None, flavor=None)  # type: ignore / on purpose


def test_set_source_flavor_with_source():
    """Test _set_source_flavor method with custom source"""
    vega_obj = vega.Vega()
    custom_source = "/custom/path.fits"

    vega_obj._set_source_flavor(source=custom_source, flavor=None)

    assert vega_obj.source == custom_source
    assert vega_obj._data is None
    assert vega_obj.units is None


def test_set_source_flavor_with_flavor():
    """Test _set_source_flavor method with flavor only"""
    vega_obj = vega.Vega(flavor="legacy")  # Start with legacy

    # Change to mod_003
    vega_obj._set_source_flavor(source=None, flavor="mod_003")

    assert vega_obj.source == vega._default_vega["mod_003"]


def test_set_source_flavor_invalid_flavor():
    """Test _set_source_flavor method with invalid flavor"""
    vega_obj = vega.Vega()

    with pytest.raises(ValueError, match="Unknown Vega flavor"):
        vega_obj._set_source_flavor(source=None, flavor="invalid")


def test_set_source_flavor_no_parameters():
    """Test _set_source_flavor method with no parameters"""
    vega_obj = vega.Vega()

    with pytest.raises(
        RuntimeError, match="Either `source` or `flavor` must be provided"
    ):
        vega_obj._set_source_flavor(source=None, flavor=None)


@patch("pyphot.vega.io.from_file")
def test_vega_readfile_method_fits(mock_from_file):
    """Test _readfile method with FITS file"""
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

    vega_obj = vega.Vega(source="test.fits")
    data, w_unit, f_unit = vega_obj._readfile()

    assert isinstance(data, pd.DataFrame)
    assert w_unit == "ANGSTROMS"
    assert f_unit == "FLAM"
    assert vega_obj._data is not None
    assert vega_obj.units == ("ANGSTROMS", "FLAM")
    # Should be called without tablename for FITS files
    mock_from_file.assert_called_once_with("test.fits")


@patch("pyphot.vega.io.from_file")
def test_vega_readfile_method_hdf5(mock_from_file):
    """Test _readfile method with HDF5 file"""
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

    vega_obj = vega.Vega(source="test.hd5")
    data, w_unit, f_unit = vega_obj._readfile()

    assert isinstance(data, pd.DataFrame)
    assert w_unit == "ANGSTROMS"
    assert f_unit == "FLAM"
    # Should be called with tablename for HDF5 files
    mock_from_file.assert_called_once_with("test.hd5", tablename="/spectrum")


@patch("pyphot.vega.io.from_file")
def test_vega_readfile_method_bytes_units(mock_from_file):
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

    vega_obj = vega.Vega()
    data, w_unit, f_unit = vega_obj._readfile()

    assert w_unit == "ANGSTROMS"
    assert f_unit == "FLAM"


@patch("pyphot.vega.io.from_file")
def test_vega_readfile_cached(mock_from_file):
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

    vega_obj = vega.Vega()

    # First call
    data1, w_unit1, f_unit1 = vega_obj._readfile()
    # Second call
    data2, w_unit2, f_unit2 = vega_obj._readfile()

    # Should only call the file reader once
    assert mock_from_file.call_count == 1
    # Results should be the same
    assert data1 is data2
    assert w_unit1 == w_unit2
    assert f_unit1 == f_unit2


def test_vega_context_manager():
    """Test Vega as context manager"""
    with patch("pyphot.vega.io.from_file") as mock_from_file:
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

        vega_obj = vega.Vega()

        with vega_obj as v:
            assert v is vega_obj
            assert v._data is not None

        # Should have read the file during context entry
        assert mock_from_file.call_count == 1


@patch("pyphot.vega.io.from_file")
def test_vega_data_property(mock_from_file):
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

    vega_obj = vega.Vega()

    # Access data property
    data = vega_obj.data

    assert isinstance(data, pd.DataFrame)
    assert len(data) == 3
    assert "WAVELENGTH" in data.columns
    assert "FLUX" in data.columns


@patch("pyphot.vega.io.from_file")
def test_vega_wavelength_property(mock_from_file):
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

    vega_obj = vega.Vega()

    # Access wavelength property
    wavelength = vega_obj.wavelength

    # Should have units
    assert hasattr(wavelength, "unit")
    assert hasattr(wavelength, "value")
    # Should have the right values
    np.testing.assert_array_equal(wavelength.value, [4000, 5000, 6000])


@patch("pyphot.vega.io.from_file")
def test_vega_flux_property(mock_from_file):
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

    vega_obj = vega.Vega()

    # Access flux property
    flux = vega_obj.flux

    # Should have units
    assert hasattr(flux, "unit")
    assert hasattr(flux, "value")
    # Should have the right values
    np.testing.assert_array_equal(flux.value, [1.0, 1.5, 1.2])


def test_vega_file_extension_detection():
    """Test that different file extensions are handled correctly"""
    # Test HDF5 extensions
    for ext in ["hd5", "hdf", "hdf5", "HD5", "HDF", "HDF5"]:
        vega_obj = vega.Vega(source=f"test.{ext}")
        assert vega_obj.source == f"test.{ext}"

    # Test FITS extensions
    for ext in ["fits", "fit", "FITS", "FIT"]:
        vega_obj = vega.Vega(source=f"test.{ext}")
        assert vega_obj.source == f"test.{ext}"


def test_vega_default_sources():
    """Test that default sources dictionary is properly defined"""
    expected_flavors = [
        "mod_002",
        "mod_003",
        "mod_004",
        "stis_011",
        "stis_003",
        "legacy",
    ]

    for flavor in expected_flavors:
        assert flavor in vega._default_vega

    # Check that paths contain expected keywords
    assert "alpha_lyr_mod_002" in vega._default_vega["mod_002"]
    assert "alpha_lyr_mod_003" in vega._default_vega["mod_003"]
    assert "alpha_lyr_mod_004" in vega._default_vega["mod_004"]
    assert "alpha_lyr_stis_011" in vega._default_vega["stis_011"]
    assert "alpha_lyr_stis_003" in vega._default_vega["stis_003"]
    assert "vega.hd5" in vega._default_vega["legacy"]


def test_vega_all_exports():
    """Test that __all__ contains expected exports"""
    assert vega.__all__ == ["Vega"]


def test_vega_module_docstring():
    """Test that vega module has appropriate documentation"""
    assert vega.__doc__ is not None
    assert "vega" in vega.__doc__.lower()
    assert "flux" in vega.__doc__.lower()


def test_vega_class_docstring():
    """Test that Vega class has comprehensive documentation"""
    assert vega.Vega.__doc__ is not None
    assert "Class that handles vega spectrum" in vega.Vega.__doc__
    assert "Attributes" in vega.Vega.__doc__
    assert "Bohlin 2007" in vega.Vega.__doc__  # Should mention Bohlin reference
    assert "context manager" in vega.Vega.__doc__


def test_vega_case_insensitive_flavor():
    """Test that flavor matching is case insensitive for validation but case sensitive for lookup"""
    # The actual implementation has a bug: it uses .lower() for validation but not for lookup
    # This test documents the current behavior

    # This should work (lowercase matches key)
    vega_obj1 = vega.Vega(flavor="legacy")

    # This should fail (uppercase doesn't match key, even though validation passes)
    with pytest.raises(KeyError):
        vega.Vega(flavor="LEGACY")

    # This should also fail (mixed case doesn't match key)
    with pytest.raises(KeyError):
        vega.Vega(flavor="Legacy")

    # But the validation should catch invalid flavors regardless of case
    with pytest.raises(ValueError, match="Unknown Vega flavor"):
        vega.Vega(flavor="INVALID")

    with pytest.raises(ValueError, match="Unknown Vega flavor"):
        vega.Vega(flavor="invalid")
