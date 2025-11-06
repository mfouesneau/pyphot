"""Test helpers module functionality"""

import warnings
from unittest.mock import Mock, patch
from math import pi

from pyphot.helpers import progress_enumerate, deprecated, distc


class TestProgressEnumerate:
    """Test the progress_enumerate function"""

    def test_progress_enumerate_without_progress(self):
        """Test progress_enumerate with show_progress=False"""
        items = [1, 2, 3, 4, 5]
        result = list(progress_enumerate(items, show_progress=False))
        expected = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]
        assert result == expected

    def test_progress_enumerate_default_no_progress(self):
        """Test progress_enumerate with default parameters (no progress)"""
        items = ["a", "b", "c"]
        result = list(progress_enumerate(items))
        expected = [(0, "a"), (1, "b"), (2, "c")]
        assert result == expected

    def test_progress_enumerate_with_start_parameter(self):
        """Test progress_enumerate with start parameter"""
        items = [10, 20, 30]
        result = list(progress_enumerate(items, 5, show_progress=False))
        expected = [(5, 10), (6, 20), (7, 30)]
        assert result == expected

    @patch("pyphot.helpers.Pbar")
    def test_progress_enumerate_with_progress(self, mock_pbar_class):
        """Test progress_enumerate with show_progress=True"""
        # Setup mock
        mock_pbar_instance = Mock()
        mock_pbar_class.return_value = mock_pbar_instance
        mock_pbar_instance.iterover.return_value = [1, 2, 3]

        items = [1, 2, 3]
        result = list(progress_enumerate(items, show_progress=True))
        expected = [(0, 1), (1, 2), (2, 3)]

        assert result == expected
        mock_pbar_class.assert_called_once()
        mock_pbar_instance.iterover.assert_called_once_with(items)

    @patch("pyphot.helpers.Pbar")
    def test_progress_enumerate_with_progress_and_kwargs(self, mock_pbar_class):
        """Test progress_enumerate with show_progress=True and additional kwargs"""
        # Setup mock
        mock_pbar_instance = Mock()
        mock_pbar_class.return_value = mock_pbar_instance
        mock_pbar_instance.iterover.return_value = ["x", "y"]

        items = ["x", "y"]
        result = list(progress_enumerate(items, show_progress=True, desc="Testing"))
        expected = [(0, "x"), (1, "y")]

        assert result == expected
        mock_pbar_class.assert_called_once_with(desc="Testing")
        mock_pbar_instance.iterover.assert_called_once_with(items)

    def test_progress_enumerate_empty_sequence(self):
        """Test progress_enumerate with empty sequence"""
        items = []
        result = list(progress_enumerate(items, show_progress=False))
        assert result == []

    def test_progress_enumerate_single_item(self):
        """Test progress_enumerate with single item"""
        items = [42]
        result = list(progress_enumerate(items, show_progress=False))
        expected = [(0, 42)]
        assert result == expected


class TestDeprecatedDecorator:
    """Test the deprecated decorator"""

    def test_deprecated_decorator_warns(self):
        """Test that deprecated decorator issues a DeprecationWarning"""
        test_message = "This function is deprecated"

        @deprecated(test_message)
        def test_function():
            return "test_result"

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = test_function()

            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert str(w[0].message) == test_message
            assert result == "test_result"

    def test_deprecated_decorator_preserves_function_behavior(self):
        """Test that deprecated decorator preserves original function behavior"""

        @deprecated("Test deprecation")
        def add_numbers(a, b):
            return a + b

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = add_numbers(3, 4)
            assert result == 7

    def test_deprecated_decorator_with_args_and_kwargs(self):
        """Test deprecated decorator with function that takes args and kwargs"""

        @deprecated("Function with args and kwargs is deprecated")
        def complex_function(a, b, c=None, d=10):
            return a + b + (c or 0) + d

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = complex_function(1, 2, c=3, d=4)
            assert result == 10

    def test_deprecated_decorator_preserves_function_metadata(self):
        """Test that deprecated decorator preserves function metadata using wraps"""
        original_docstring = "Original function docstring"
        original_name = "original_function"

        @deprecated("Test deprecation")
        def original_function():
            """Original function docstring"""
            pass

        # The @wraps decorator should preserve the original function's metadata
        assert original_function.__name__ == original_name
        assert original_function.__doc__ == original_docstring

    def test_deprecated_decorator_stacklevel(self):
        """Test that deprecated decorator uses correct stacklevel"""
        test_message = "Stacklevel test"

        @deprecated(test_message)
        def test_function():
            return True

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            test_function()

            # Verify that the warning points to the correct location (stacklevel=2)
            assert len(w) == 1
            # The filename should point to this test file, not the helpers.py file
            assert "test_helpers.py" in w[0].filename

    def test_deprecated_multiple_calls(self):
        """Test that deprecated decorator warns on each call"""

        @deprecated("Multiple calls test")
        def test_function():
            return "result"

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            # Call the function multiple times
            test_function()
            test_function()
            test_function()

            # Should have one warning per call
            assert len(w) == 3
            for warning in w:
                assert issubclass(warning.category, DeprecationWarning)
                assert str(warning.message) == "Multiple calls test"


def test_distc_value():
    """Test that distc constant has the expected value"""
    # distc = 4.0 * pi * (3.0856775e19) ** 2
    expected_value = 4.0 * pi * (3.0856775e19) ** 2
    assert distc == expected_value
