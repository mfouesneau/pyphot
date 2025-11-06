"""Adapted from SimpleTable (v2.0; https://github.com/mfouesneau/simpletable) with minimal depedencies"""

from typing import Optional
import pandas as pd
import warnings

from . import fits
from . import ascii
from . import hdf
from . import votable

from .header import HeaderInfo

__all__ = ["fits", "ascii", "hdf", "votable", "HeaderInfo", "from_file"]


def from_file(
    fname: str,
    *,
    format: Optional[str] = None,
    **kwargs,
) -> tuple[pd.DataFrame, HeaderInfo]:
    """Read a file into a DataFrame and a Header

    Parameters
    ----------
    fname : str
        File name to read.
    format : str, optional
        File format to read. If not provided, the format is inferred from the file extension.
    **kwargs
        Additional keyword arguments to pass to the reader.

    Returns
    -------
    df : pd.DataFrame
        DataFrame containing the data.
    hdr : HeaderInfo
        Header information.

    """
    # handle legacy files by extension
    if format:
        ext = format.lower()
    else:
        ext = fname.split(".")[-1].lower()
    if ext in ("hd5", "hdf", "hdf5"):
        df, hdr = hdf.from_hdf5(fname, **kwargs)
    elif ext in ("csv"):
        df, hdr = ascii.from_csv(fname, **kwargs)
    elif ext in ("txt", "dat"):
        df, hdr = ascii.from_ascii(fname, **kwargs)
    elif ext in ("vot", "votable", "xml"):
        df, hdr = votable.from_votable(fname, **kwargs)
    elif ext in ("fits",):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df, hdr = fits.from_fits(fname, **kwargs)
    else:
        raise ValueError(f"Unsupported file format: {ext}")
    df = pd.DataFrame(df)
    df.attrs.update(hdr.header)
    df.attrs["units"] = hdr.units
    df.attrs["WAVELENGTH_UNIT"] = hdr.units["WAVELENGTH"]
    df.attrs["FLUX_UNIT"] = hdr.units["FLUX"].split("=")[0].rstrip()

    return df, hdr


def to_file(
    df: pd.DataFrame,
    fname: str,
    *,
    format: Optional[str] = None,
    header_info: Optional[HeaderInfo] = None,
    **kwargs,
):
    """Write a file from a DataFrame and a Header

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write.
    fname : str
        File name to read.
    header_info: Optional[HeaderInfo] = None
        Header information to write.
    format : str, optional
        File format to read. If not provided, the format is inferred from the file extension.
    **kwargs
        Additional keyword arguments to pass to the reader.

    """
    # handle legacy files by extension
    if format:
        ext = format.lower()
    else:
        ext = fname.split(".")[-1].lower()
    if ext in ("hd5", "hdf", "hdf5"):
        hdf.to_hdf5(df, fname, header_info=header_info, **kwargs)
    elif ext in ("csv"):
        ascii.to_csv(df, fname, **kwargs)
    elif ext in ("txt", "dat"):
        ascii.to_ascii(df, fname, **kwargs)
    elif ext in ("fits",):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fits.to_fits(df, fname, header_info=header_info, **kwargs)
    else:
        raise ValueError(f"Unsupported file format: {ext}")
