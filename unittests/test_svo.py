"""Test the future SVO filter interface"""

import pytest

import numpy as np

from pyphot import config
from pyphot import svo
from pyphot.phot import Filter
from pyphot.unit_adapters import backends

pb_refs = [
    "2MASS/2MASS.J",
    "2MASS/2MASS.H",
    "2MASS/2MASS.Ks",
    "HST/ACS_WFC.F475W",
    "HST/ACS_WFC.F814W",
]


@pytest.mark.parametrize("name", pb_refs)
@pytest.mark.parametrize("unit_backend", list(backends.keys()))
def test_svo_filter_interface(unit_backend, name):
    """Test the future SVO filter interface"""
    if unit_backend is not None:
        try:
            config.set_units_backend(unit_backend)
        except ImportError:
            pytest.skip(f"Skipping test for {unit_backend} adapter")

    pb = svo.get_pyphot_filter(name)
    pb.info()
    assert isinstance(pb, Filter)
    assert pb.name == name.replace("/", "_")


def compare_filters(filter1: Filter, filter2: Filter):
    """Compare two filters by their transmission curves.

    Parameters
    ----------
    filter1 : Filter
        First filter to compare.
    filter2 : Filter
        Second filter to compare.

    Returns
    -------
    float
        The difference between the two filters.
    """
    attrs = set(filter1.__dict__.keys()) | set(filter2.__dict__.keys())
    for attr in attrs:
        val1, val2 = getattr(filter1, attr), getattr(filter2, attr)
        if "unit" in attr:
            if val1 != val2:
                # this test is independent of the unit backend
                assert config.units.U(val1) == config.units.U(val2), (
                    f"Unit mismatch for attribute {attr}: {val1} != {val2}"
                )
        elif isinstance(val1, str) and isinstance(val2, str):
            assert val1 == val2, (
                f"Value mismatch for attribute {attr}: {val1} != {val2}"
            )
        else:
            assert np.allclose(val1, val2), (
                f"Value mismatch for attribute {attr}: {val1} != {val2}"
            )


@pytest.mark.parametrize("name", pb_refs)
@pytest.mark.parametrize("unit_backend", list(backends.keys()))
def test_orig_intern_consistency(name, unit_backend):
    """Test the consistency of the original astropy and internal filter retrieval.
    The switch to internal allows for more flexibility in handling filter units.
    """
    if unit_backend is not None:
        try:
            config.set_units_backend(unit_backend)
        except ImportError:
            pytest.skip(f"Skipping test for {unit_backend} adapter")

    pb1 = svo.get_pyphot_filter(name)
    pb2 = svo._get_pyphot_filter_astropy(name)
    compare_filters(pb1, pb2)
