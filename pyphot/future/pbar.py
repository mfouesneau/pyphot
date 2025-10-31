"""
Simple progressbar
==================

This package implement a unique progress bar class that can be used to decorate
an iterator, a function or even standalone.

The format of the meter is flexible and can display along with the progress
meter, the running time, an eta, and the rate of the iterations.

An example is::
    description    [----------] k/n  10% [time: 00:00:00, eta: 00:00:00, 2.7 iters/sec]
"""

import time as _time
import sys
import signal
from array import array
from fcntl import ioctl
import termios


__all__ = ["Pbar"]


class Pbar(object):
    """
    make a progress string  in a shape of::

        [----------] k/n  10% [time: 00:00:00, eta: 00:00:00, 2.7 iters/sec]


    Attributes
    ----------

    time: bool, optional (default: True)
        if set, add the runtime information

    eta: bool, optional (default: True)
        if set, add an estimated time to completion

    rate: bool, optional (default: True)
        if set, add the rate information

    length: int, optional (default: None)
        number of characters showing the progress meter itself
        if None, the meter will adapt to the buffer width

        TODO: make it variable with the buffer length

    keep: bool, optional (default: True)
        If not set, deletes its traces from screen after completion

    file: buffer
        the buffer to write into

    mininterval: float (default: 0.5)
        minimum time in seconds between two updates of the meter

    miniters: int, optional (default: 1)
        minimum iteration number between two updates of the meter

    units: str, optional (default: 'iters')
        unit of the iteration
    """

    def __init__(
        self,
        maxval=None,
        desc=None,
        time=True,
        eta=True,
        rate=True,
        length=None,
        file=None,
        keep=True,
        mininterval=0.5,
        miniters=1,
        units="iters",
        **kwargs,
    ):
        self.time = time
        self.eta = eta
        self.rate = rate
        self.desc = desc or ""
        self.units = units
        self.file = file or sys.stdout
        self._last_print_len = 0
        self.keep = keep
        self.mininterval = mininterval
        self.miniters = miniters
        self._auto_width = True
        self.length = 10
        if length is not None:
            self.length = length
            self._auto_width = False
        # backward compatibility
        self._start_t = _time.time()
        self._maxval = maxval
        if "txt" in kwargs:
            self.desc = kwargs["txt"]

    def _buffer_width(self):
        """returns the width of the buffer when available"""
        try:
            self.handle_resize(None, None)
            signal.signal(signal.SIGWINCH, self.handle_resize)
            self._auto_width = True
        except Exception:
            self.term_width = 79
            self._auto_width = False

        return self.term_width

    def handle_resize(self, signum, frame):
        ioinfo = ioctl(self.file, termios.TIOCGWINSZ, "\0" * 8)  # type: ignore
        h, w = array("h", ioinfo)[:2]
        self.term_width = w

    @staticmethod
    def format_interval(t):
        """make a human readable time interval decomposed into days, hours,
        minutes and seconds

        Parameters
        ----------
        t: int
            interval in seconds

        Returns
        -------
        txt: str
            string representing the interval
            (format:  <days>d <hrs>:<min>:<sec>)
        """
        mins, s = divmod(int(t), 60)
        h, m = divmod(mins, 60)
        d, h = divmod(h, 24)

        txt = "{m:02d}:{s:02d}"
        if h:
            txt = "{h:02d}:" + txt
        if d:
            txt = "{d:d}d " + txt
        return txt.format(d=d, h=h, m=m, s=s)

    def build_str_meter(self, n, total, elapsed):
        """
        make a progress string  in a shape of::

            [----------] k/n  10% [time: 00:00:00, eta: 00:00:00, 2.7 iters/sec]

        Parameters
        ----------
        n: int
            number of finished iterations

        total: int
            total number of iterations, or None

        elapsed: int
            number of seconds passed since start

        Returns
        -------
        txt: str
            string representing the meter
        """
        if n > total:
            total = None

        vals = {"n": n}
        vals["elapsed"] = self.format_interval(elapsed)
        vals["rate"] = "{0:5.2f}".format((n / elapsed)) if elapsed else "?"
        vals["units"] = self.units

        if not total:
            txt = "{n:d}"
        else:
            txt = "|{bar:s}| {n:d}/{total:d} {percent:s}"

        if self.time or self.eta or self.rate:
            txt += " ["
            info = []
            if self.time:
                info.append("time: {elapsed:s}")
            if self.eta and total:
                info.append("eta: {left:s}")
            if self.rate:
                info.append("{rate:s} {units:s}/sec")
            txt += ", ".join(info) + "]"

        if not total:
            return txt.format(**vals)

        frac = float(n) / total
        bar_length = int(frac * self.length)
        vals["bar"] = "#" * bar_length + "-" * (self.length - bar_length)
        vals["percent"] = "{0:3.0%}".format(frac)
        vals["left"] = self.format_interval(elapsed / n * (total - n)) if n else "?"
        vals["total"] = total

        if self._auto_width:
            full_length = self._buffer_width()
            current_length = len(txt.format(**vals))
            new_length = full_length - current_length + self.length - 1 - len(self.desc)
            frac = float(n) / total
            bar_length = int(frac * new_length)
            vals["bar"] = "#" * bar_length + "-" * (new_length - bar_length)

        return txt.format(**vals)

    def print_status(self, s):
        """print a status s  on the last file line and clean the rest of the line

        Parameters
        ----------
        s: str
            message to write
        """
        self.file.write("\r" + s + " " * max(self._last_print_len - len(s), 0))
        self.file.flush()
        self._last_print_len = len(s)

    def iterover(self, iterable, total=None):
        """
        Get an iterable object, and return an iterator which acts exactly like the
        iterable, but prints a progress meter and updates it every time a value is
        requested.

        Parameters
        ----------
        iterable: generator or iterable object
            object to iter over.

        total: int, optional
            the number of iterations is assumed to be the length of the
            iterator.  But sometimes the iterable has no associated length or
            its length is not the actual number of future iterations. In this
            case, total can be set to define the number of iterations.

        Returns
        -------
        gen: generator
            pass the values from the initial iterator
        """
        if total is None:
            try:
                total = len(iterable)
            except TypeError:
                total = self._maxval

        prefix = "{0:s}:".format(self.desc) if self.desc else ""

        self.print_status(prefix + self.build_str_meter(0, total, 0))
        last_print_n = 0

        start_t = last_print_t = _time.time()

        n = 0
        for obj in iterable:
            yield obj
            n += 1
            if n - last_print_n >= self.miniters:
                cur_t = _time.time()
                if cur_t - last_print_t >= self.mininterval:
                    self.print_status(
                        prefix + self.build_str_meter(n, total, cur_t - start_t)
                    )
                    last_print_n = n
                    last_print_t = cur_t

        if not self.keep:
            self.print_status("")
            sys.stdout.write("\r")
        else:
            if last_print_n < n:
                cur_t = _time.time()
                self.print_status(
                    prefix + self.build_str_meter(n, total, cur_t - start_t)
                )
            self.file.write("\n")

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        return False

    def update(self, n, desc=None, total=None):
        """Kept for backward compatibility and the decorator feature

        Parameters
        ----------
        n: int
            force iteration number n

        desc: str
            update description string

        total: int
            update the total number of iterations
        """
        if total is None:
            total = self._maxval
        if desc is not None:
            self.desc = desc
        prefix = "{0:s}:".format(self.desc) if self.desc else ""
        cur_t = _time.time()
        self.print_status(
            prefix + self.build_str_meter(n, total, cur_t - self._start_t)
        )

    def decorator(self, func):
        """Provide a function decorator allowing for counting calls and rates"""
        self._deco_iter = 0
        self.desc = func.__name__

        def deco(*args, **kwargs):
            # start the time at the first call
            if self._deco_iter == 0:
                self._start_t = _time.time()
            self.update(self._deco_iter)
            self._deco_iter += 1
            return func(*args, **kwargs)

        return deco
