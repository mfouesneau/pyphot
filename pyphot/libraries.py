"""Library defintions as Collections of Filters"""

import os
import pathlib
from typing import List, Literal, Optional, Sequence, Union

import numpy as np
import numpy.typing as npt
import tables

from . import config
from .helpers import progress_enumerate
from .phot import Filter, QuantityType


class Library:
    """Common grounds for filter libraries"""

    source: Optional[str]
    """Source of the library"""

    def __init__(
        self,
        source: Optional[str] = None,
        *args,
        **kwargs,
    ):
        """Construct the library"""
        self.source = source or str(config.__default_passband_lib__)

    def __repr__(self) -> str:
        msg = "Filter Library: {0}\n{1:s}"
        return msg.format(self.source, object.__repr__(self))

    def __enter__(self):
        """Enter context"""
        return self

    def __exit__(self, *exc_info):
        """end context"""
        return False

    def __len__(self) -> int:
        """Size of the library"""
        return len(self.content)

    def to_csv(self, directory="./", progress=True, **kwargs):
        """Export each filter into a csv file with its own name

        Parameters
        ----------
        directory: str
            directory to write into
        progress: bool
            show progress if set

        """
        os.makedirs(directory, exist_ok=True)

        def export_filter(f: Filter):
            if f.wavelength_unit in (None, ""):
                f.wavelength_unit = "AA"
            f.write_to(
                f"{directory:s}/{f.name:s}.csv".lower(),
                fmt="%.6f",
                **kwargs,
            )

        with self as s:
            for _, k in progress_enumerate(
                s.content, desc="export", show_progress=progress
            ):
                f = s[k]
                # s[k] can return a list
                if isinstance(f, Filter):
                    export_filter(f)
                else:
                    for fk in f:
                        export_filter(fk)

    def to_hdf(self, fname="filters.hd5", progress=True, **kwargs):
        """Export each filter into a csv file with its own name

        Parameters
        ----------
        directory: str
            directory to write into
        progress: bool
            show progress if set

        """

        def export_filter(f):
            if f.wavelength_unit in (None, ""):
                f.wavelength_unit = "AA"
            f.write_to(
                f"{fname:s}",
                tablename=f"/filters/{f.name}",
                createparents=True,
                append=True,
                silent=True,
                **kwargs,
            )

        with self as s:
            for _, k in progress_enumerate(
                s.content, desc="export", show_progress=progress
            ):
                f = s[k]
                if isinstance(f, Filter):
                    export_filter(f)
                else:
                    for fk in f:
                        export_filter(fk)

    @classmethod
    def from_hd5(cls, filename, **kwargs) -> "HDF_Library":
        """Read in an HDF5 library"""
        return HDF_Library(filename, **kwargs)

    @classmethod
    def from_ascii(cls, filename, **kwargs) -> "Ascii_Library":
        """Read in an ASCII library"""
        return Ascii_Library(filename, **kwargs)

    @property
    def content(self) -> List[str]:
        """Get the content list"""
        return self.get_library_content()

    def __getitem__(
        self, name: Union[str, Sequence[str]]
    ) -> Union[Filter, List[Filter]]:
        """Make this object like a dictionary and load one or multiple filters"""
        with self as s:
            try:
                f = s._load_filter(name)
            except TypeError:
                f = [s._load_filter(k) for k in name]
        return f

    def _load_filter(self, *args, **kwargs) -> Filter:
        """Load a given filter from the library"""
        raise NotImplementedError

    def get_library_content(self) -> List[str]:
        """get the content of the library"""
        raise NotImplementedError

    def load_all_filters(
        self,
        *,
        interp: bool = True,
        lamb: Optional[Union[npt.NDArray[np.floating], QuantityType]] = None,
    ) -> List[Filter]:
        """load all filters from the library

        Parameters
        ----------
        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        Returns
        -------
        filters: list[filter]
            list of filter objects
        """
        raise NotImplementedError

    def add_filter(self, f: Filter):
        """add a filter to the library"""
        raise NotImplementedError

    def find(self, name: str, case_sensitive=True) -> List[str]:
        """Search for a filter in the library"""
        r = []
        if case_sensitive:
            _n = name.lower()
            for k in self.get_library_content():
                if _n in k.lower():
                    r.append(k)
        else:
            for k in self.content:
                if name in k:
                    r.append(k)
        return r


class Ascii_Library(Library):
    """Interface one or multiple directory or many files as a filter :class:`Library`

    >>> lib = Ascii_Library(["ground", "hst", "myfilter.csv"])
    """

    def _load_filter(
        self,
        fname: str,
        interp: bool = True,
        lamb: Union[None, npt.NDArray[np.floating], QuantityType] = None,
        *args,
        **kwargs,
    ):
        """Load a given filter from the library

        Parameters
        ----------
        fname : str
            Name of the filter to load.
        interp : bool, optional
            Whether to interpolate the filter to a given wavelength grid.
        lamb : array_like, optional
            Wavelength grid to interpolate the filter to.
        args : tuple, optional
            Additional arguments to pass to the filter constructor.
        kwargs : dict, optional
            Additional keyword arguments to pass to `Filter.from_ascii`.

        Returns
        -------
        Filter
            The loaded filter.
        """
        try:  # attempt to load filter from ascii file
            fil = Filter.from_ascii(fname, *args, **kwargs)
        except Exception:
            content = self.content
            r = [k for k in content if fname in k]
            # if nothing matched, try all lowercase for names
            if len(r) <= 0:
                r = [k for k in content if fname.lower() in k]

            if len(r) > 1:
                raise ValueError(
                    "Auto correction found multiple choices."
                    "Refine name to one of {}".format(r)
                )
            elif len(r) <= 0:
                raise ValueError(f"Cannot find filter {fname}")
            else:
                fil = Filter.from_ascii(r[0], *args, **kwargs)
        if (interp is True) and (lamb is not None):
            return fil.reinterp(lamb)
        else:
            return fil

    def get_library_content(self) -> List[str]:
        """get the content of the library"""
        from glob import glob

        if self.source is None:
            raise ValueError("Library source not set")

        # Assume source is either a directory or a pattern or a single file
        try:
            os.path.isdir(self.source)
            lst = glob(self.source + "/*")
        except TypeError:
            lst = [self.source]

        # expand directories
        dircheck = True
        while dircheck is True:
            dircheck = False
            newlst = []
            for entry in lst:
                if os.path.isdir(entry):
                    newlst.extend(glob(entry + "/*"))
                    dircheck = True
                else:
                    newlst.append(entry)
            lst = newlst

        return lst

    def load_all_filters(
        self,
        *,
        interp: bool = True,
        lamb: Optional[Union[npt.NDArray[np.floating], QuantityType]] = None,
    ) -> List[Filter]:
        """load all filters from the library

        Parameters
        ----------
        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        Returns
        -------
        filters: list[filter]
            list of filter objects
        """
        return self.load_filters(self.content, interp=interp, lamb=lamb)

    def load_filters(
        self,
        names: List[str],
        *,
        interp: bool = True,
        lamb: Optional[Union[npt.NDArray, QuantityType]] = None,
    ) -> List[Filter]:
        """load a limited set of filters

        Parameters
        ----------
        names: list[str]
            normalized names according to filtersLib

        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        filterLib: path
            path to the filter library hd5 file

        Returns
        -------
        filters: list[filter]
            list of filter objects
        """
        filters = [
            self._load_filter(fname, interp=interp, lamb=lamb)
            for fname in names
        ]
        return filters

    def add_filters(
        self,
        filter_object: Filter,
        fmt="%.6f",
        **kwargs,
    ):
        """Add a filter to the library permanently

        Parameters
        ----------
        filter_object: Filter object
            filter to add
        """
        if not isinstance(filter_object, Filter):
            msg = "Argument of type Filter expected. Got type {0}"
            raise TypeError(msg.format(type(filter_object)))

        if filter_object.wavelength_unit is None:
            msg = "Filter wavelength must have units for storage."
            raise AttributeError(msg)
        fname = f"{self.source:s}/{filter_object.name:s}.csv"
        filter_object.write_to(fname.lower(), fmt=fmt, **kwargs)


class HDF_Library(Library):
    """:class:`Library` for storage based on HDF files"""

    hdf: Optional[tables.File]
    """Source file stream of the library"""
    mode: "Literal['r', 'w', 'a', 'r+']" = "r"
    """Mode of the library (file). It can be one of the following:

            * *'r'*: Read-only; no data can be modified.
            * *'w'*: Write; a new file is created (an existing file
              with the same name would be deleted).
            * *'a'*: Append; an existing file is opened for reading
              and writing, and if the file does not exist it is created.
            * *'r+'*: It is similar to 'a', but the file must already
              exist.
    """
    _in_context: int
    """Number of times the library is in context (potentially nested)"""

    def __init__(
        self,
        source: Optional[str] = None,
        mode: "Literal['r', 'w', 'a', 'r+']" = "r",
    ):
        super().__init__(source)
        self.hdf = None
        self.mode = mode
        self._in_context = 0

    def __enter__(self):
        """Enter context"""
        if self.source is None:
            raise ValueError("Source must be provided")

        if self.hdf is None:
            self.hdf = tables.open_file(self.source, self.mode)
        self._in_context += 1
        return self

    def __exit__(self, *exc_info):
        """end context"""
        if (self.hdf is not None) and (self._in_context < 2):
            self.hdf.close()
            self.hdf = None
        self._in_context -= 1
        return False

    def _load_filter(
        self,
        fname: str,
        interp: bool = True,
        lamb: Union[None, npt.NDArray[np.floating], QuantityType] = None,
    ) -> Filter:
        """Load a given filter from the library

        Parameters
        ----------
        fname: str
            normalized names according to filtersLib

        interp: bool, optional
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        integrationFilter: bool, optional
            set True for specail integraion filter such as Qion or E_uv
            if set, lamb should be given

        Returns
        -------
        filter: Filter instance
            filter object
        """
        with self as s:
            ftab = s.hdf
            if ftab is None:
                raise ValueError("Library not initialized")

            if hasattr(fname, "decode"):
                fnode = ftab.get_node("/filters/" + fname.decode("utf8"))  # type: ignore
            else:
                fnode = ftab.get_node("/filters/" + fname)
            flamb = fnode[:]["WAVELENGTH"]
            transmit = fnode[:]["THROUGHPUT"]
            dtype = "photon"
            unit = None
            attrs = fnode.attrs
            if "DETECTOR" in attrs:
                dtype = attrs["DETECTOR"]
            if "WAVELENGTH_UNIT" in attrs:
                unit = attrs["WAVELENGTH_UNIT"]

            fil = Filter(
                flamb,
                transmit,
                name=fnode.name,
                dtype=dtype,
                unit=unit,
            )

        if (lamb is not None) and interp:
            fil = fil.reinterp(lamb)
        return fil

    def get_library_content(self) -> List[str]:
        """get the content of the library"""
        with self as s:
            if s.hdf is None:
                raise ValueError("Library not initialized")
            try:
                filters = s.hdf.root.content.cols.TABLENAME[:]
            except Exception:
                filters = list(s.hdf.root.filters._v_children.keys())
        if hasattr(filters[0], "decode"):
            filters = [k.decode("utf8") for k in filters]
        return filters

    def load_all_filters(
        self,
        *,
        interp: bool = True,
        lamb: Optional[Union[npt.NDArray[np.floating], QuantityType]] = None,
    ) -> List[Filter]:
        """load all filters from the library

        Parameters
        ----------
        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        Returns
        -------
        filters: list[filter]
            list of filter objects
        """
        return self.load_filters(self.content, interp=interp, lamb=lamb)

    def load_filters(
        self,
        names: List[str],
        *,
        interp: bool = True,
        lamb: Optional[Union[npt.NDArray, QuantityType]] = None,
    ) -> List[Filter]:
        """load a limited set of filters

        Parameters
        ----------
        names: list[str]
            normalized names according to filtersLib

        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        filterLib: path
            path to the filter library hd5 file

        Returns
        -------
        filters: list[filter]
            list of filter objects
        """
        with self as s:
            filters = [
                s._load_filter(fname, interp=interp, lamb=lamb)
                for fname in names
            ]
        return filters

    def add_filter(self, f: Filter, **kwargs):
        """Add a filter to the library permanently

        Parameters
        ----------
        f: Filter object
            filter to add
        """
        if not isinstance(f, Filter):
            msg = "Argument of type Filter expected. Got type {0}"
            raise TypeError(msg.format(type(f)))

        if f.wavelength_unit is None:
            msg = "Filter wavelength must have units for storage."
            raise AttributeError(msg)

        append = kwargs.pop("append", True)

        f.write_to(
            f"{self.source:s}",
            tablename=f"/filters/{f.name}",
            createparents=True,
            append=append,
            **kwargs,
        )


def get_library(fname: Optional[str] = None, **kwargs):
    """Finds the appropriate class to load the library"""
    fname = fname or str(config.__default_passband_lib__)
    library_path = pathlib.Path(fname)
    if (library_path.suffix in (".hdf5", ".hd5", ".h5")) and (
        library_path.is_file()
    ):
        return HDF_Library(fname, **kwargs)
    else:
        return Ascii_Library(fname, **kwargs)
