"""Test the future SVO filter interface"""

from pyphot.svo import get_pyphot_filter
from pyphot.phot import Filter
from pyphot import config


def test_svo_filter_interface(unit_backend=None):
    """Test the future SVO filter interface"""
    lst = "2MASS/2MASS.J 2MASS/2MASS.H 2MASS/2MASS.Ks HST/ACS_WFC.F475W HST/ACS_WFC.F814W".split()
    if unit_backend is not None:
        config.set_units_backend(unit_backend)

    for name in lst:
        pb = get_pyphot_filter(name)
        pb.info()
        assert isinstance(pb, Filter)
        assert pb.name == name.replace("/", "_")
