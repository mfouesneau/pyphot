"""Testing units adapters"""

import pytest

from typing import List, Tuple, cast
from pyphot.future.unit_adapters import get_adapter, backends, UnitAdapterType

# make sure we test all available backends
test_backends = [name for name, adapter in backends.items() if adapter is not None]


@pytest.mark.parametrize("Adapter", list(backends.values()))
def test_adapter_decorator(Adapter: UnitAdapterType):
    """Test the ezunits related decorator"""

    if Adapter is None:
        pytest.skip(f"Skipping test for {Adapter.__name__} adapter")

    U = Adapter.U
    decorator = Adapter.decorate

    def test_func(x, y, scale=1.0):
        return (x + y) * scale

    assert test_func(1, 2) == 3

    # add units on the input arguments
    fn1 = decorator("kg", "kg")(test_func)
    assert fn1(1, 2) == 3 * U("kg")
    assert fn1(1000 * U("g"), 2) == 3 * U("kg")

    # add units on the keywords
    fn2 = decorator("kg", "kg", scale="m * s**(-1)")(test_func)
    assert fn2(1, 2, scale=1.0) == 3 * U("kg * m * s**(-1)")
    assert fn2(1, 2) == 3 * U("kg * m * s**(-1)")
    assert fn2(1 * U("kg"), 2) == 3 * U("kg * m * s**(-1)")

    # test invalid unit dimensionality
    try:
        fn2(1 * U("mm"), 2)
    except Adapter.typing.DimensionalityError:
        pass

    # test invalid units
    try:
        _ = 1 * U("non_existing_unit")
    except Adapter.typing.UndefinedUnitError:
        pass


@pytest.mark.parametrize("Adapter", list(backends.values()))
def test_pyphot_units(Adapter: UnitAdapterType):
    """Checking correct definitions of needed units"""
    if Adapter is None:
        pytest.skip(f"Skipping test for {Adapter.__name__} adapter")

    U = Adapter.U
    # Astronomy related
    values: List[Tuple[str, Adapter.typing.Quantity]] = [
        ("lsun", 3.828e26 * U("watt")),
        ("erg", 1e-7 * U("joule")),
        ("flam", 1.0 * U("erg * s ** (-1) * AA ** (-1) * cm **(-2)")),
        ("fnu", 1.0 * U("erg * s ** (-1) * Hz ** (-1) * cm **(-2)")),
        ("photflam", 1.0 * U("photon * s ** (-1) * AA ** (-1) * cm **(-2)")),
        ("photfnu", 1.0 * U("photon * s ** (-1) * Hz ** (-1) * cm **(-2)")),
        ("angstroms", 1e-10 * U("m")),
    ]

    for what, ref in values:
        # forcing cast as StructuredQuantity does not define value
        ratio = cast(Adapter.typing.Quantity, (U(what)) / ref.to(what))
        assert (
            abs(ratio.value - 1) < 1e-10
        ), f"Error: '{what}' does not match reference {ref} ({ratio.value})"
