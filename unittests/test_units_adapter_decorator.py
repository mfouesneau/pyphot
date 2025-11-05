"""Testing units decorator"""

import pytest

from pyphot import config
from pyphot.unit_adapters import enforce_default_units, backends


@pytest.mark.parametrize("backend", list(backends.keys()))
def test_decorator(backend):
    """Test the ezunits related decorator"""

    # set default backedn
    config.set_units_backend(backend)
    # get the default backend decorator
    decorator = enforce_default_units

    # get the specfic backend for comparison
    adapter = config.get_units_adapter(backend)
    U = adapter.U

    # define a test function with some docstring
    def test_func(x, y, scale=1.0):
        """
        Test function

        Parameters
        ----------
        x : float
            First argument
        y : float
            Second argument
        scale : float, optional
            Scaling factor, by default 1.0

        Returns
        -------
        float
            Sum of x and y scaled by scale
        """
        return (x + y) * scale

    # define a test function with some docstring
    def test_func2(x, y, scale=1.0):
        """
        Test function

        Parameters
        ----------
        x : float
            First argument
        y : float
            Second argument
        scale : float, optional
            Scaling factor, by default 1.0

        Returns
        -------
        float
            Sum of x and y scaled by scale
        """
        return (x + y) * scale

    # test the naked function
    assert test_func(1, 2) == 3

    # add units on the input arguments
    fn1 = decorator("kg", "kg")(test_func)
    assert fn1(1, 2) == 3 * U("kg")
    assert fn1(1000 * U("g"), 2) == 3 * U("kg")

    # add units on the keywords
    fn2 = decorator("kg", "kg", scale="m * s**(-2)")(test_func)
    assert fn2(1, 2, scale=2.0) == 6 * U("kg * m * s**(-2)")
    assert fn2(1, 2) == 3 * U("kg * m * s**(-2)")
    assert fn2(1 * U("kg"), 2) == 3 * U("kg * m * s**(-2)")

    # test invalid unit dimensionality
    try:
        fn2(1 * U("mm"), 2)
    except adapter.typing.DimensionalityError:
        pass

    # test invalid units
    try:
        _ = 1 * U("non_existing_unit")
    except adapter.typing.UndefinedUnitError:
        pass
