from typing import Any, Tuple, Union, Optional
from astropy.io import fits
from astropy.io.fits import Header
from astropy.io.fits import TableHDU, BinTableHDU
import numpy.typing as npt
import pandas as pd
from .header import HeaderInfo


try:
    from astropy.io import fits

    def fix_endian_issue(arr: Union[npt.NDArray, Any]) -> npt.NDArray:
        """Fix endian issue in array which happens often when reading FITS files"""
        return arr.astype(arr.dtype.newbyteorder())

    def fits_read_header(hdr: Header) -> HeaderInfo:
        """
        Convert pyfits header into dictionary with relevant values

        Parameters
        ----------

        hdr: pyftis.Header
            fits unit

        Returns
        -------
        headerinfo: HeaderInfo
            extracted information from header
        """
        header = {}
        alias = {}
        units = {}
        comments = {}

        # generic cards
        genTerms = [
            "XTENSION",
            "BITPIX",
            "NAXIS",
            "NAXIS1",
            "NAXIS2",
            "PCOUNT",
            "GCOUNT",
            "TFIELDS",
            "EXTNAME",
        ]
        fieldTerms = ["TTYPE", "TFORM", "TUNIT", "ALIAS"]

        for card in hdr.cards["TTYPE*"]:
            name = card.value
            comments[name] = card.comment
            u = hdr.get(card.keyword.replace("TYPE", "UNIT"), None)
            if u is not None:
                units[name] = u

        # for k, val, _ in hdr.ascard['ALIAS*']:
        for card in hdr.cards["ALIAS*"]:
            k = card.keyword
            val = card.value
            al, orig = val.split("=")
            alias[al] = orig

        # other specific keywords: COMMENT, HISTORY
        header_comments = []
        header_history = []
        for k, v in hdr.items():
            if (k not in genTerms) and (k[:5] not in fieldTerms):
                if k == "COMMENT":
                    header_comments.append(v)
                elif k == "HISTORY":
                    header_history.append(v)
                else:
                    header[k] = v

        # COMMENT, HISTORY polish
        if len(header_comments) > 0:
            header["COMMENT"] = "\n".join(header_comments)
        if len(header_history) > 0:
            header["HISTORY"] = "\n".join(header_history)

        if "EXTNAME" in hdr:
            header["NAME"] = hdr["EXTNAME"]

        return HeaderInfo(header, alias, units, comments)

    def from_fits(
        filename: str, extension_number: int = 1
    ) -> Tuple[npt.NDArray, HeaderInfo]:
        """Load a DataFrame from a FITS file.

        Parameters
        ----------
        filename : str
            The path to the FITS file.
        extension_number : int, optional
            The extension number to load, by default 1.

        Returns
        -------
        Tuple[npt.NDArray, HeaderInfo]
            The loaded DataFrame and its header information.
        """
        with fits.open(filename) as hdu:
            extension: Union[TableHDU, BinTableHDU] = hdu[extension_number]  # pyright: ignore[reportAssignmentType]
            hdr_info = fits_read_header(extension.header)
            data = fix_endian_issue(extension.data)
        return data, hdr_info

except ImportError:
    from_fits = None  # pyright: ignore[reportAssignmentType]


def fits_generate_header(df: pd.DataFrame) -> fits.Header:
    """Generate the corresponding fits Header that contains all necessary info

    Parameters
    ----------
    df: pd.DataFrame
        DataFrame or HeaderInfo instance

    Returns
    -------
    hdr: fits.Header
        header instance
    """
    # extract necessary pieces
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pd.DataFrame or HeaderInfo instance")

    others = ["aliases", "units", "comments"]
    info = HeaderInfo(
        header={k: v for k, v in df.attrs.items() if k not in others},
        alias=df.attrs.get("aliases", {}),
        units=df.attrs.get("units", {}),
        comments=df.attrs.get("comments", {}),
    )

    # Set column cards
    columns = list(df.columns)
    cards = []

    for e, colname in enumerate(columns):
        # TTYPE<n=1..N> = <colname> // comment or ""
        cards.append((f"TTYPE{e + 1:d}", colname, info.comments.get(colname, "")))
        # TUNIT<n=1..N> = <colunit> // unit of <colname>
        colunit = info.units.get(colname, "")
        if colunit not in ("", "None", None):
            cards.append((f"TUNIT{e + 1:d}", colunit, f"unit of {colname:s}"))

    # add aliases ALIAS<n=1..N> = <from>=<to> // alias of <from>
    for e, (to_, from_) in enumerate(info.alias.items()):
        cards.append(
            ("ALIAS{0:d}".format(e + 1), f"{to_:s}={from_:s}", f"alias of {from_:s}")
        )

    # add NAME if any
    table_name = info.header.get("NAME", None)
    if table_name not in ["", "None", None, "No Name"]:
        cards.append(("EXTNAME", table_name, ""))

    # generate header
    hdr = fits.Header(cards)

    # add other header information, COMMENT, and HISTORY cards
    for key, value in info.header.items():
        if value in ("", "None", None):
            continue
        if key not in ("NAME", "COMMENT", "HISTORY"):
            hdr.update({key: value})
        elif key == "COMMENT":
            txt = value.split("\n")
            for j in txt:
                hdr.add_comment(j)
        elif key == "HISTORY":
            txt = value.split("\n")
            for j in txt:
                hdr.add_history(j)
    return hdr


def fits_generate_hdu(df: pd.DataFrame, index: bool = True) -> BinTableHDU:
    """Generate a FITS BinTableHDU from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to convert.
    index : bool, optional
        Whether to include the index in the table, by default True.

    Returns
    -------
    BinTableHDU
        The generated HDU.
    """
    table_name = df.attrs.get("NAME", None)
    records = df.to_records(index=index)
    header = fits_generate_header(df)
    hdu = fits.BinTableHDU(records, header=header, name=table_name)
    # astropy bug: missing some keywords (e.g. TUNIT*) from above.
    # forcing full propagation of header information
    hdu.header.update(header)
    return hdu


def to_fits(
    df: pd.DataFrame,
    filename: str,
    extension_number: int = 1,
    header_info: Optional[HeaderInfo] = None,
    output_verify: str = "exception",
    checksum: bool = False,
    index: bool = True,
    overwrite: bool = False,
    append: bool = False,
    **kwargs,
) -> None:
    """Save a DataFrame to a FITS file.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to save.
        header information taken from data.attrs or header_info if provided
    filename : str
        The path to the FITS file.
    extension_number : int, optional
        The extension number to save, by default 1.
    header_info : Optional[HeaderInfo], optional
        Header information to save with the FITS file
        by default None and taken from data.attrs
        override data.attrs if provided
    output_verify : str
        Output verification option.  Must be one of ``"fix"``, ``"silentfix"``,
        ``"ignore"``, ``"warn"``, or ``"exception"``.  May also be any
        combination of ``"fix"`` or ``"silentfix"`` with ``"+ignore"``,
        ``+warn``, or ``+exception" (e.g. ``"fix+warn"``).  See :ref:`verify`
        for more info.
    checksum : bool, optional
        If `True`, adds both ``DATASUM`` and ``CHECKSUM`` cards to the
        headers of all HDU's written to the file
    index : bool, optional
        If `True`, includes the index in the FITS file. Default is `True`.
    append : bool, optional
        If `True`, appends the DataFrame to the FITS file. Default is `False`.
    overwrite : bool, optional
        If `True`, overwrites the DataFrame in the FITS file. Default is `False`.
    **kwargs: dict
        Additional keyword arguments to pass to the FITS writer.
    """
    hdu = fits_generate_hdu(df, index=index)

    name, closed, noexist_or_empty = fits.convenience._stat_filename_or_fileobj(
        filename
    )

    if not append or noexist_or_empty:
        hdu.writeto(
            filename,
            checksum=checksum,
            output_verify=output_verify,
            overwrite=overwrite,
        )
    else:
        if not closed:
            fits_file = fits.convenience.fitsopen(filename, mode="append")
            fits_file.append(hdu)
            # Set a flag in the HDU so that only this HDU gets a checksum when
            # writing the file.
            hdu._output_checksum = checksum
            # Keep status of the file as it was
            fits_file.close(output_verify=output_verify, closed=closed)
        else:
            fits_file = fits.convenience._File(filename, mode="append")
            hdu._output_checksum = checksum
            hdu._writeto(fits_file)
            fits_file.close()
