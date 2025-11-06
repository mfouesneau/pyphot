"""Test pbar module functionality"""

import sys
import io
from unittest.mock import Mock, patch
from array import array

from pyphot.pbar import Pbar


def test_pbar_init_defaults():
    """Test Pbar initialization with default parameters"""
    pbar = Pbar()

    assert pbar.time is True
    assert pbar.eta is True
    assert pbar.rate is True
    assert pbar.desc == ""
    assert pbar.units == "iters"
    assert pbar.file == sys.stdout
    assert pbar._last_print_len == 0
    assert pbar.keep is True
    assert pbar.mininterval == 0.5
    assert pbar.miniters == 1
    assert pbar._auto_width is True
    assert pbar.length == 10
    assert pbar._maxval is None


def test_pbar_init_custom_parameters():
    """Test Pbar initialization with custom parameters"""
    custom_file = io.StringIO()
    pbar = Pbar(
        maxval=100,
        desc="Testing",
        time=False,
        eta=False,
        rate=False,
        length=20,
        file=custom_file,
        keep=False,
        mininterval=1.0,
        miniters=5,
        units="items",
    )

    assert pbar._maxval == 100
    assert pbar.desc == "Testing"
    assert pbar.time is False
    assert pbar.eta is False
    assert pbar.rate is False
    assert pbar.length == 20
    assert pbar.file == custom_file
    assert pbar.keep is False
    assert pbar.mininterval == 1.0
    assert pbar.miniters == 5
    assert pbar.units == "items"
    assert pbar._auto_width is False


def test_pbar_init_backward_compatibility():
    """Test Pbar initialization with backward compatibility txt parameter"""
    pbar = Pbar(txt="Legacy description")
    assert pbar.desc == "Legacy description"


def test_format_interval():
    """Test the format_interval static method"""
    # Test seconds only
    assert Pbar.format_interval(45) == "00:45"

    # Test minutes and seconds
    assert Pbar.format_interval(125) == "02:05"

    # Test hours, minutes, and seconds
    assert Pbar.format_interval(3665) == "01:01:05"

    # Test days, hours, minutes, and seconds
    assert Pbar.format_interval(90061) == "1d 01:01:01"

    # Test zero
    assert Pbar.format_interval(0) == "00:00"

    # Test large values
    assert Pbar.format_interval(172861) == "2d 01:01"


def test_build_str_meter_no_total():
    """Test build_str_meter with no total (indeterminate progress)"""
    pbar = Pbar(time=True, eta=True, rate=True)

    # Use 0 as total to trigger the "not total" behavior
    result = pbar.build_str_meter(50, 0, 10.5)

    assert "50" in result
    assert "00:10" in result  # elapsed time
    assert "4.76 iters/sec" in result  # rate
    assert "|" not in result  # no progress bar for indeterminate


def test_build_str_meter_with_total():
    """Test build_str_meter with known total"""
    pbar = Pbar(time=True, eta=True, rate=True, length=10)

    result = pbar.build_str_meter(25, 100, 10.0)

    assert "|" in result
    assert "25/100" in result
    assert "25%" in result
    assert "##--" in result  # progress bar (25% of 10 chars = 2.5 ≈ 2 chars)
    assert "00:10" in result  # elapsed time
    assert "2.50 iters/sec" in result  # rate
    assert "00:30" in result  # eta (remaining 75 items at 2.5/sec = 30 sec)


def test_build_str_meter_complete():
    """Test build_str_meter at 100% completion"""
    pbar = Pbar(length=5)

    result = pbar.build_str_meter(100, 100, 20.0)

    assert "100/100" in result
    assert "100%" in result
    assert "#####" in result  # full progress bar
    assert "00:20" in result


def test_build_str_meter_overtime():
    """Test build_str_meter when n > total"""
    pbar = Pbar(length=5)

    # When n > total, the code sets total to None internally
    # This actually causes the behavior to switch to indeterminate mode
    result = pbar.build_str_meter(150, 100, 20.0)

    assert "150" in result
    assert "00:20" in result
    assert "7.50 iters/sec" in result


def test_build_str_meter_options():
    """Test build_str_meter with different option combinations"""
    # Test with time only
    pbar = Pbar(time=True, eta=False, rate=False)
    result = pbar.build_str_meter(50, 100, 10.0)
    assert "time:" in result
    assert "eta:" not in result
    assert "iters/sec" not in result

    # Test with eta only
    pbar = Pbar(time=False, eta=True, rate=False)
    result = pbar.build_str_meter(50, 100, 10.0)
    assert "time:" not in result
    assert "eta:" in result
    assert "iters/sec" not in result

    # Test with rate only
    pbar = Pbar(time=False, eta=False, rate=True)
    result = pbar.build_str_meter(50, 100, 10.0)
    assert "time:" not in result
    assert "eta:" not in result
    assert "iters/sec" in result

    # Test with no options
    pbar = Pbar(time=False, eta=False, rate=False)
    result = pbar.build_str_meter(50, 100, 10.0)
    assert "[" not in result  # no info section


def test_build_str_meter_zero_elapsed():
    """Test build_str_meter with zero elapsed time"""
    pbar = Pbar(rate=True)

    result = pbar.build_str_meter(10, 100, 0)

    assert "?" in result  # rate should be ? when elapsed is 0


def test_build_str_meter_zero_progress():
    """Test build_str_meter with zero progress"""
    pbar = Pbar(eta=True)

    result = pbar.build_str_meter(0, 100, 5.0)

    assert "0/100" in result
    assert "eta: ?" in result  # eta should be ? when no progress made


def test_print_status():
    """Test print_status method"""
    output = io.StringIO()
    pbar = Pbar(file=output)

    # First status
    pbar.print_status("Testing 123")
    assert output.getvalue() == "\rTesting 123"

    # Reset for next test
    output.truncate(0)
    output.seek(0)

    # Longer status
    pbar.print_status("Much longer testing message")
    result = output.getvalue()
    assert result.startswith("\rMuch longer testing message")

    # Reset for next test
    output.truncate(0)
    output.seek(0)

    # Shorter status (should pad with spaces)
    pbar.print_status("Short")
    result = output.getvalue()
    assert result.startswith("\rShort")
    assert " " in result  # should contain padding spaces


@patch("pyphot.pbar._time")
def test_iterover_basic(mock_time):
    """Test iterover method with basic functionality"""
    # Provide more mock time values to handle all calls
    mock_time.time.side_effect = [0.0, 0.6, 1.2, 1.8, 2.4, 3.0]

    output = io.StringIO()
    pbar = Pbar(
        file=output, mininterval=0.1, miniters=1
    )  # Lower mininterval for testing

    items = [1, 2, 3]
    result = list(pbar.iterover(items))

    assert result == [1, 2, 3]
    output_str = output.getvalue()
    assert "3/3" in output_str
    assert "100%" in output_str


@patch("pyphot.pbar._time")
def test_iterover_with_total(mock_time):
    """Test iterover method with explicit total"""
    mock_time.time.side_effect = [0.0, 0.6, 1.2, 1.8]

    output = io.StringIO()
    pbar = Pbar(file=output, mininterval=0.1)  # Lower mininterval

    items = [1, 2]
    result = list(pbar.iterover(items, total=5))

    assert result == [1, 2]
    output_str = output.getvalue()
    assert "2/5" in output_str


@patch("pyphot.pbar._time")
def test_iterover_no_keep(mock_time):
    """Test iterover method with keep=False"""
    mock_time.time.side_effect = [0.0, 0.6, 1.2, 1.8]

    output = io.StringIO()
    pbar = Pbar(file=output, keep=False, mininterval=0.1)

    items = [1, 2]
    list(pbar.iterover(items))

    output_str = output.getvalue()
    assert "\r" in output_str  # should contain carriage return


def test_iterover_with_description():
    """Test iterover method with description"""
    output = io.StringIO()
    pbar = Pbar(file=output, desc="Processing")

    items = [1]
    list(pbar.iterover(items))

    output_str = output.getvalue()
    assert "Processing:" in output_str


def test_iterover_no_length():
    """Test iterover with object that has no length"""
    output = io.StringIO()
    pbar = Pbar(file=output, maxval=3)

    def generator():
        yield 1
        yield 2
        yield 3

    result = list(pbar.iterover(generator()))
    assert result == [1, 2, 3]


@patch("pyphot.pbar._time")
def test_iterover_mininterval_miniters(mock_time):
    """Test iterover respects mininterval and miniters"""
    mock_time.time.side_effect = [0.0, 0.1, 0.2, 0.3, 1.0]  # Fast iterations

    output = io.StringIO()
    pbar = Pbar(file=output, mininterval=0.5, miniters=2)

    items = [1, 2, 3, 4]
    list(pbar.iterover(items))

    # Should not update on every iteration due to mininterval/miniters


def test_context_manager():
    """Test Pbar as context manager"""
    pbar = Pbar()

    with pbar as p:
        assert p is pbar

    # Should exit without error


@patch("pyphot.pbar._time")
def test_update_method(mock_time):
    """Test update method for manual progress updates"""
    mock_time.time.return_value = 5.0

    output = io.StringIO()
    pbar = Pbar(file=output)
    pbar._start_t = 0.0  # Set start time

    pbar.update(50, desc="Updated", total=100)

    assert pbar.desc == "Updated"
    output_str = output.getvalue()
    assert "Updated:" in output_str
    assert "50/100" in output_str


@patch("pyphot.pbar._time")
def test_decorator_functionality(mock_time):
    """Test decorator method"""
    mock_time.time.side_effect = [0.0, 1.0, 2.0, 3.0]

    output = io.StringIO()
    pbar = Pbar(file=output, maxval=10)  # Set maxval to avoid None total

    @pbar.decorator
    def test_function(x):
        return x * 2

    # Call the decorated function
    result1 = test_function(5)
    result2 = test_function(10)

    assert result1 == 10
    assert result2 == 20
    assert pbar.desc == "test_function"

    # Check that output contains progress information
    output_str = output.getvalue()
    assert len(output_str) > 0  # Should have printed something


@patch("pyphot.pbar.signal")
@patch("pyphot.pbar.ioctl")
def test_buffer_width_success(mock_ioctl, mock_signal):
    """Test _buffer_width method when terminal info is available"""
    # Mock successful ioctl call
    mock_ioctl.return_value = None

    pbar = Pbar()
    pbar.handle_resize = Mock()
    pbar.handle_resize.return_value = None

    # Mock the handle_resize to set term_width
    def side_effect(signum, frame):
        pbar.term_width = 80

    pbar.handle_resize.side_effect = side_effect

    width = pbar._buffer_width()

    assert pbar._auto_width is True
    mock_signal.signal.assert_called_once()


@patch("pyphot.pbar.signal")
def test_buffer_width_exception(mock_signal):
    """Test _buffer_width method when terminal info is not available"""
    # Mock signal.signal to raise an exception
    mock_signal.signal.side_effect = Exception("No terminal")

    pbar = Pbar()
    width = pbar._buffer_width()

    assert pbar.term_width == 79
    assert pbar._auto_width is False


@patch("pyphot.pbar.ioctl")
def test_handle_resize(mock_ioctl):
    """Test handle_resize method"""
    # Mock ioctl to return width and height info
    mock_ioctl.return_value = array("h", [25, 100, 0, 0]).tobytes()

    output = io.StringIO()
    pbar = Pbar(file=output)

    pbar.handle_resize(None, None)

    assert pbar.term_width == 100


def test_all_exports():
    """Test that __all__ contains expected exports"""
    from pyphot.pbar import __all__

    assert __all__ == ["Pbar"]


def test_pbar_with_custom_units():
    """Test Pbar with custom units"""
    pbar = Pbar(units="files", rate=True)

    result = pbar.build_str_meter(10, 50, 5.0)

    assert "files/sec" in result


@patch("pyphot.pbar._time")
def test_auto_width_feature(mock_time):
    """Test auto width adjustment feature"""
    mock_time.time.return_value = 1.0

    pbar = Pbar(length=10, _auto_width=True)
    pbar.term_width = 100  # Mock terminal width
    pbar._buffer_width = Mock(return_value=100)

    result = pbar.build_str_meter(25, 100, 10.0)

    # Should adjust bar length based on terminal width
    assert "|" in result
    assert "#" in result


def test_module_docstring():
    """Test that pbar module has appropriate documentation"""
    import pyphot.pbar as pbar_module

    assert pbar_module.__doc__ is not None
    assert "progressbar" in pbar_module.__doc__.lower()
    assert "example" in pbar_module.__doc__.lower()


def test_pbar_class_docstring():
    """Test that Pbar class has comprehensive documentation"""
    assert Pbar.__doc__ is not None
    assert "progress string" in Pbar.__doc__
    assert "Attributes" in Pbar.__doc__
    assert "time:" in Pbar.__doc__
    assert "eta:" in Pbar.__doc__
    assert "rate:" in Pbar.__doc__
