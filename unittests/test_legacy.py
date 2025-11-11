"""Test legacy module functionality"""

import pytest
from typing import Any

import pyphot.legacy as legacy
from pyphot.phot import Filter
from pyphot.libraries import Library, HDF_Library, Ascii_Library

mappings = {
    "UnitFilter": Filter,
    "UnitLibrary": Library,
    "UnitHDFLibrary": HDF_Library,
    "UnitAsciiLibrary": Ascii_Library,
}


@pytest.mark.parametrize("legacy_name", mappings.keys())
def test_legacy_imports_available(legacy_name: str):
    """Test that all legacy imports are available in the module"""
    assert hasattr(legacy, legacy_name)


@pytest.mark.parametrize("legacy_name, target", mappings.items())
def test_alias_points_to_reference(legacy_name: str, target: Any):
    """Test that aliases point to the correct references"""
    source = getattr(legacy, legacy_name)
    assert source is target
    assert source.__name__ == target.__name__
    assert source.__module__ == target.__module__
    assert source.__doc__ == target.__doc__
    assert source.__bases__ == target.__bases__
    assert source.__mro__ == target.__mro__


def test_all_exports_defined():
    """Test that __all__ contains all expected exports"""
    expected_exports = [
        "UnitFilter",
        "UnitLibrary",
        "UnitHDFLibrary",
        "UnitAsciiLibrary",
    ]
    assert legacy.__all__ == expected_exports


def test_all_exports_accessible():
    """Test that all items in __all__ are accessible in the module"""
    for export in legacy.__all__:
        assert hasattr(legacy, export)
        # Make sure the attribute is not None
        assert getattr(legacy, export) is not None


def test_from_legacy_import_star():
    """Test that 'from legacy import *' works correctly"""
    # This simulates what happens with 'from pyphot.legacy import *'
    namespace = {}
    exec("from pyphot.legacy import *", namespace)

    # Check that all expected symbols are imported
    for symbol in legacy.__all__:
        assert symbol in namespace
        assert namespace[symbol] is getattr(legacy, symbol)


def test_no_additional_attributes():
    """Test that legacy module doesn't define unexpected attributes"""
    expected_attrs = set(
        legacy.__all__
        + [
            "__doc__",
            "__file__",
            "__loader__",
            "__name__",
            "__package__",
            "__spec__",
            "__all__",
            "__builtins__",
        ]
    )
    actual_attrs = set(dir(legacy))

    # Allow for some standard module attributes that might be present
    unexpected = actual_attrs - expected_attrs
    # Filter out any private attributes that are commonly added by Python
    unexpected = {attr for attr in unexpected if not attr.startswith("_")}

    assert len(unexpected) == 0, f"Unexpected attributes found: {unexpected}"
