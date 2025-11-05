from typing import Optional
import numpy.typing as npt
import pandas as pd
from .header import HeaderInfo

import tables


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
):
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
    with tables.open_file(filename, **kwargs) as source:
        tablename = tablename or "/"
        if not tablename.startswith("/"):
            tablename = "/" + tablename

        node = source.get_node(tablename)
        attrs = node._v_attrs
        if not silent:
            print("\tLoading table: {0}".format(tablename))

        header = {}
        aliases = {}

        # read header
        exclude = ["NROWS", "VERSION", "CLASS", "EXTNAME", "TITLE"]
        for k in attrs._v_attrnames:
            if k not in exclude:
                if not k.startswith("FIELD") and not k.startswith("ALIAS"):
                    header[k] = _decode_string_ifneeded(attrs[k])
                elif k.startswith("ALIAS"):
                    c0, c1 = _decode_string_ifneeded(attrs[k]).split("=")
                    aliases[c0] = c1

        if attrs["TITLE"] not in ["", "None", "Noname", None]:
            header["NAME"] = attrs["TITLE"]
        else:
            header["NAME"] = f"{filename:s}/{node._v_name:s}"

        # read column meta
        units = {}
        desc = {}

        colnames = node.colnames  # type: ignore / it does exist
        print(attrs)
        for k, colname in enumerate(colnames):
            _u = getattr(attrs, f"FIELD_{k:d}_UNIT", None)
            _u = getattr(attrs, f"{colname:s}_UNIT", _u)
            _d = getattr(attrs, f"FIELD_{k:d}_DESC", None)
            _d = getattr(attrs, f"{colname:s}_DESC", _d)
            if _u is not None:
                units[colname] = _decode_string_ifneeded(_u)
            if _d is not None:
                desc[colname] = _decode_string_ifneeded(_d)

        data = node[:]  # type: ignore / it is defined.

        hdr = HeaderInfo(
            header=header,
            alias=aliases,
            units=units,
            comments=desc,
        )

        return pd.DataFrame(data), hdr
    raise ValueError("Something went wrong without much information from pytables")
