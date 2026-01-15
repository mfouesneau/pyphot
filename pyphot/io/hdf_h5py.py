"""Read and write HDF5 files with pytables preserving metadata (tables, https://www.pytables.org/)

.. important::
    This module relies on `pytables <https://www.pytables.org/>`_

"""

from os import PathLike
from os.path import basename
from typing import Literal, Optional, Union

import h5py
import numpy as np
import pandas as pd

from .header import HeaderInfo


def _decode_string_ifneeded(s: str) -> str:
    """Silently decode a string if it is bytes"""
    if isinstance(s, bytes):
        return s.decode("utf-8")
    return s


def from_hdf5(
    filename: str,
    tablename: Optional[str] = None,
    *,
    silent: bool = True,
    **kwargs,
) -> tuple[pd.DataFrame, HeaderInfo]:
    """Generate the corresponding ascii Header that contains all necessary info

    Parameters
    ----------
    filename: str
        file to read from

    tablename: str
        node containing the table

    silent: bool
        skip verbose messages

    Returns
    -------
    hdr: str
        string that will be be written at the beginning of the file
    """
    with h5py.File(filename, **kwargs) as source:
        tablename = tablename or "/"
        if not tablename.startswith("/"):
            tablename = "/" + tablename

        node = source[tablename]
        if not node:
            raise ValueError(
                f"Table '{tablename}' not found in file '{filename}'"
            )
        if not isinstance(node, h5py.Dataset):
            raise TypeError(
                f"Node '{tablename}' is not a dataset (found {type(node)})"
            )

        attrs = node.attrs
        if not silent:
            print(f"\tLoading table: {tablename}")

        header = {}
        aliases = {}

        # read header
        exclude = ["NROWS", "VERSION", "CLASS", "EXTNAME", "TITLE"]
        for k, v in attrs.items():
            if k not in exclude:
                if not k.startswith("FIELD") and not k.startswith("ALIAS"):
                    header[k] = _decode_string_ifneeded(v)
                elif k.startswith("ALIAS"):
                    c0, c1 = _decode_string_ifneeded(v).split("=")
                    aliases[c0] = c1

        title = attrs.get("TITLE", "")
        if title not in ["", "None", "Noname", None]:
            header["NAME"] = _decode_string_ifneeded(title)
        else:
            header["NAME"] = f"{filename:s}/{node.name:s}"

        # read column meta
        units = {}
        desc = {}

        colnames = node.dtype.names
        for k, colname in enumerate(colnames):
            _u = attrs.get(f"FIELD_{k:d}_UNIT", None)
            _u = attrs.get(f"{colname:s}_UNIT", _u)
            _d = attrs.get(f"FIELD_{k:d}_DESC", None)
            _d = attrs.get(f"{colname:s}_DESC", _d)
            if _u is not None:
                units[colname] = _decode_string_ifneeded(_u)
            if _d is not None:
                desc[colname] = _decode_string_ifneeded(_d)

        data = node[:]

        hdr = HeaderInfo(
            header=header,
            alias=aliases,
            units=units,
            comments=desc,
        )

        return pd.DataFrame(data), hdr
    raise ValueError(
        "Something went wrong without much information from pytables"
    )


def to_hdf5(
    df: pd.DataFrame,
    filename: Union[str, h5py.File, PathLike],
    *,
    tablename: Optional[str] = None,
    header_info: Optional[HeaderInfo] = None,
    mode: "Literal['r', 'w', 'a', 'r+', 'w-', 'x']" = "w",
    append: bool = False,
    **kwargs,
) -> None:
    """
    Write a pandas DataFrame to an HDF5 file.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to write.
    filename : str or tables.File or PathLike
        The filename or open HDF5 file to write to.
    tablename : str, optional
        The name of the table to write to.
    header_info : HeaderInfo, optional
        The header information to write. Default is to use from df.attrs
    mode : {'r', 'w', 'a', 'r+'}, default 'w'
        The mode to open the file in.
    append : bool, default False
        Whether to append data to an existing file.
    **kwargs
        Additional keyword arguments to pass to h5py.File.

    Raises
    ------
    Exception
        If the HDF backend does not implement stream.
    tables.FileModeError
        If the file is already opened in a different mode.
    ValueError
        If something went wrong without much information from pytables.
    """
    if hasattr(filename, "read"):
        raise Exception("HDF backend does not implement stream")

    # ensure mode is valid for appending data
    mode = "a" if append is True else mode

    # open output file, or if provided tables.File, check it's in the correct mode
    if isinstance(filename, h5py.File):
        if (filename.mode != mode) & (mode != "r"):
            raise RuntimeError(
                f"The file {basename(filename.filename)} is already opened in a different mode (mode {filename.mode})"
            )
        hd5 = filename
    else:
        hd5 = h5py.File(str(filename), mode=mode)

    if header_info is None:
        # attempt to get it from the dataframe attributes
        header_info = HeaderInfo(
            header={
                k: v
                for k, v in df.attrs.items()
                if k not in ["aliases", "units", "comments"]
            },
            alias=df.attrs.get("aliases", {}),
            units=df.attrs.get("units", {}),
            comments=df.attrs.get("comments", {}),
        )

    # check table name and path
    tablename_path = (
        tablename
        or df.name
        or header_info.header.get("name", None)
        or header_info.header.get("NAME", None)
    )
    if tablename_path in ("", None, "Noname", "None"):
        tablename_path = "/data"
    if not tablename_path.startswith("/"):
        tablename_path = "/" + tablename_path

    if append:
        try:
            tab = hd5[tablename_path]
            if not isinstance(tab, h5py.Dataset):
                raise TypeError(f"Node is not a table (got {type(tab)})")
            if tab.dtype is None:
                raise ValueError("Table dtype description is missing")

            dtypes = tab.dtype  # pyright: ignore / it is there
            data = df.to_records(index=False, column_dtypes=dtypes)
            data_length = data.shape[0]

            # try resize the table first
            # this requires the table to be resizable
            try:
                tab.resize(tab.shape[0] + data_length, axis=0)
                # write the data to the table
                tab[-data_length:] = data
            except TypeError:
                # Only chunked datasets can be resized
                # We need to make a new data array with all data and replace the node
                attrs = {k: v for k, v in tab.attrs.items()}
                new_data = np.hstack([tab[:], data])
                del hd5[tablename_path]
                tab = hd5.create_dataset(
                    tablename_path, data=new_data, dtype=dtypes
                )
                tab.attrs.update(attrs)
            tab.flush()
        except KeyError:
            print(
                f"Warning: Table {tablename_path} does not exist. A new table will be created."
            )
            append = False

    if not append:
        # we need the path and the name separately
        w_ = tablename_path.split("/")
        where = "/".join(w_[:-1])
        tablename = str(w_[-1])
        if tablename in (None, ""):
            raise ValueError(
                "Table name cannot be empty. Did you leave a trailing slash?"
            )
        if where in ("", None):
            where = "/"

        data = df.to_records(index=False)
        tab = hd5.create_dataset(f"{where}/{tablename}", data=data, **kwargs)

        # update hdr attrs with header_info.header
        for k, v in header_info.header.items():
            tab.attrs[k] = v
        if "TITLE" not in header_info.header:
            tab.attrs["TITLE"] = tablename

        # add column descriptions and units
        for e, colname in enumerate(df.columns):
            _u = header_info.units.get(colname, None)
            _d = header_info.comments.get(colname, None)
            if _u is not None:
                tab.attrs["FIELD_{0:d}_UNIT"] = _u
            if _d is not None:
                tab.attrs["FIELD_{0:d}_DESC"] = _d

        # add aliases
        for i, (k, v) in enumerate(header_info.alias.items()):
            tab.attrs[f"ALIAS{i:d}"] = f"{k:s}={v:s}"

        tab.flush()

    if not isinstance(filename, h5py.File):
        hd5.flush()
        hd5.close()
