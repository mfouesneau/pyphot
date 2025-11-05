"""Export dataframe to ASCII format while preserving attrs"""

import pandas as pd
from typing import Callable, Hashable, Sequence, List, Tuple, cast, Union
from pandas.io.common import get_handle
from pandas._typing import (
    CompressionOptions,
    BaseBuffer,
    FilePath,
    IndexLabel,
    StorageOptions,
    OpenFileErrors,
)
from io import IOBase
from os import PathLike

from typing import Any, Dict, Hashable
from dataclasses import dataclass


@dataclass
class HeaderInfo:
    """Extracted information from FITS header

    Attributes
    ----------
        header: dict
            header dictionary

        alias: dict
            aliases

        units: dict
            units

        comments: dict
            comments/description of keywords
    """

    header: Dict[Hashable, Any]
    alias: Dict[Hashable, str]
    units: Dict[Hashable, str]
    comments: Dict[Hashable, str]


def ascii_read_header(
    fname: str | FilePath | IOBase,
    *,
    commentchar: str = "#",
    delimiter: str = ",",
    commented_header: bool = True,
    **kwargs,
) -> Tuple[int, HeaderInfo, List[str]]:
    """
    Read ASCII/CSV header

    Parameters
    ----------
    fname: str, FilePath, BaseBuffer
        File, filename, or generator to read.
        Note that generators should return byte strings for Python >=3.

    comments: str, optional
        The character used to indicate the start of a comment;
        default: '#'. ("" is equivalent to None)

    delimiter: str, optional
        The string used to separate values.  By default, this is any
        whitespace.

    commented_header: bool, optional
        if set, the last line of the header is expected to be the column titles (with comment character)
        otherwise, the first line of the data will be the column titles

    Returns
    -------
    nlines: int
        number of lines from the header

    info: HeaderInfo
        header information (header, alias, units, comments)

    names: List[str]
        sequence or str, first data line after header, expected to be the
        column names.
    """

    # define some internal functions
    def parseStrNone(v: str) -> str | None:
        """robust parse"""
        _v = v.strip().split()
        if len(_v) == 0:
            return None
        else:
            _v = " ".join(_v)
            if (_v.lower()) == "none" or (_v.lower() == "null"):
                return None
            return _v

    def parseColInfo(line: str) -> tuple[str, str | None, str | None]:
        """parse column info"""
        line = line.replace(commentchar, "").strip()
        tokens = line.split("\t")
        colname = tokens[0].strip()
        if len(tokens) > 1:
            colunit = parseStrNone(tokens[1].strip())
        else:
            colunit = None
        if len(tokens) > 2:
            colcomm = tokens[2].strip()
        else:
            colcomm = None
        return colname, colunit, colcomm

    if hasattr(fname, "read"):
        stream: IOBase = cast(IOBase, fname)
    else:
        stream: IOBase = open(fname, "r")  # pyright: ignore

    if commentchar is None:
        commentchar = "#"

    # initialize storage
    alias = {}
    units = {}
    desc = {}
    header = {}
    comment = []
    history = []

    done = False
    line = ""
    oldline = ""  # contains the last processed line
    nlines = 0  # contains the total number of lines processes in the header
    names = []  # contains the names of the columns

    while not done:
        line = str(stream.readline().rstrip())  # getting rid of '\n'
        nlines += 1

        if line.startswith(f"{commentchar}{commentchar}"):  # column info
            colname, colunit, colcomm = parseColInfo(line)
            if colunit is not None:
                units[colname] = colunit
            if colcomm is not None:
                desc[colname] = colcomm
        elif line.startswith(commentchar):  # normal header part or alias
            # header is expected as "# key \t value"
            line = line.replace(commentchar, "").strip()
            tokens = line.split("\t")
            if not line:  # skip empty lines
                continue
            if line and (len(tokens) == 1):  # assume no key as comment
                comment.append(line)
            else:
                key = tokens[0].strip()
                value = " ".join(tokens[1:]).strip()  # remove trailing spaces
                # COMMENT or HISTORY needs to be appended
                if key in ("COMMENT",):
                    comment.append(f"{value:s}")
                elif key in ("HISTORY",):
                    history.append(f"{value:s}")
                elif "alias" in key.lower():
                    # take care of aliases
                    al, orig = value.split("=")
                    alias[al] = orig
                else:
                    header[key] = value
        else:
            done = True
            if commented_header and (oldline is not None):
                names = oldline.split(delimiter)
                # remove the last line from the header part
                nlines -= 1
                if comment[-1] == oldline:
                    del comment[-1]
            else:
                names = line.split(delimiter)
        oldline = line.replace(commentchar, "").strip()

    header["COMMENT"] = "\n".join(comment)
    header["HISTORY"] = "\n".join(history)

    if not hasattr(fname, "read"):
        stream.close()
    else:
        # if stream, rewind by the length of the last line + \n
        if commented_header:
            stream.seek(stream.tell() - len(line) - 1)
        nlines = 0  # make sure the value is set to the current position

    info = HeaderInfo(
        header=header,
        alias=alias,
        units=units,
        comments=desc,
    )

    return nlines, info, names


def ascii_generate_header(
    df: pd.DataFrame,
    comments: str | None = "#",
    delimiter: str | None = " ",
    commented_header: bool = True,
) -> str:
    """Generate the corresponding ascii Header that contains all necessary info

    Parameters
    ----------

    df: pd.DataFrame
        table to export

    comments: str
        string to prepend header lines

    delimiter: str, optional
        The string used to separate values.  By default, this is any
        whitespace.

    commented_header: bool, optional
        if set, the last line of the header is expected to be the column titles

    Returns
    -------
    hdr: str
        string that will be be written at the beginning of the file
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

    if comments is None:
        comments = ""

    if delimiter is None:
        delimiter = " "

    hdr = []
    columns = list(df.keys())

    # table header keys
    length = max(len(str(k)) for k in info.header.keys())
    fmt = "{{0:s}} {{1:{0:d}s}}\t{{2:s}}".format(length)
    for key, value in info.header.items():
        for vk in str(value).split("\n"):
            if len(vk) > 0:
                hdr.append(fmt.format(comments, str(key).upper(), vk.strip()))

    # column metadata
    hdr.append(comments)  # add empty line
    length = max(len(str(k)) for k in columns)
    fmt = "{{0:s}}{{0:s}} {{1:{0:d}s}}\t{{2:s}}\t{{3:s}}".format(length)
    for colname in columns:
        unit = info.units.get(colname, "None")
        desc = info.comments.get(colname, "None")
        hdr.append(fmt.format(comments, colname, unit, desc))

    # aliases
    if info.alias:
        hdr.append(comments)  # add empty line
        for to_, from_ in info.alias.items():
            hdr.append(f"{comments:s} alias\t{to_:s}={from_:s}")

    # column names
    hdr.append(comments)  # add empty line
    if commented_header:
        hdr.append("{0:s} {1:s}".format(comments, delimiter.join(columns)))
    else:
        hdr.append("{0:s}".format(delimiter.join(columns)))

    return "\n".join(hdr) + "\n"


def to_csv(
    self: pd.DataFrame,
    filepath_or_buffer: FilePath | BaseBuffer,
    *,
    sep: "str" = ",",
    commentchar: str = "#",
    na_rep: "str" = "",
    float_format: "str | Callable | None" = None,
    columns: "Sequence[Hashable] | None" = None,
    header: bool | List[str] = True,
    index: bool = True,
    index_label: IndexLabel | None = None,
    mode: str = "w",
    encoding: str | None = None,
    compression: CompressionOptions = "infer",
    quoting: int | None = None,
    quotechar: str = '"',
    lineterminator: str | None = None,
    chunksize: int | None = None,
    date_format: str | None = None,
    doublequote: bool = True,
    escapechar: str | None = None,
    decimal: str = ".",
    errors: OpenFileErrors = "strict",
    storage_options: StorageOptions | None = None,
) -> str | None:
    r"""
    Write object to a comma-separated values (csv) file while preserving attrs

    Fallsback to `pd.DataFrame.to_csv` if no attrs content

    Parameters
    ----------
    path_or_buf : str, path object, file-like object, or None, default None
        String, path object (implementing os.PathLike[str]), or file-like
        object implementing a write() function. If None, the result is
        returned as a string. If a non-binary file object is passed, it should
        be opened with `newline=''`, disabling universal newlines. If a binary
        file object is passed, `mode` might need to contain a `'b'`.
    sep : str, default ','
        String of length 1. Field delimiter for the output file.
    commentchar : str, default '#'
        Character starting a comment line for the output file.
    na_rep : str, default ''
        Missing data representation.
    float_format : str, Callable, default None
        Format string for floating point numbers. If a Callable is given, it takes
        precedence over other numeric formatting parameters, like decimal.
    columns : sequence, optional
        Columns to write.
    header : bool or list of str, default True
        Write out the column names. If a list of strings is given it is
        assumed to be aliases for the column names.
    index : bool, default True
        Write row names (index).
    index_label : str or sequence, or False, default None
        Column label for index column(s) if desired. If None is given, and
        `header` and `index` are True, then the index names are used. A
        sequence should be given if the object uses MultiIndex. If
        False do not print fields for index names. Use index_label=False
        for easier importing in R.
    mode : {{'w', 'x', 'a'}}, default 'w'
        Forwarded to either `open(mode=)` or `fsspec.open(mode=)` to control
        the file opening. Typical values include:

        - 'w', truncate the file first.
        - 'x', exclusive creation, failing if the file already exists.
        - 'a', append to the end of file if it exists.

    encoding : str, optional
        A string representing the encoding to use in the output file,
        defaults to 'utf-8'. `encoding` is not supported if `path_or_buf`
        is a non-binary file object.
        compression : str or dict, default 'infer'
            For on-the-fly compression of the output data. If 'infer' and 'path_or_buf' is
            path-like, then detect compression from the following extensions: '.gz',
            '.bz2', '.zip', '.xz', '.zst', '.tar', '.tar.gz', '.tar.xz' or '.tar.bz2'
            (otherwise no compression).
            Set to ``None`` for no compression.
            Can also be a dict with key ``'method'`` set
            to one of {``'zip'``, ``'gzip'``, ``'bz2'``, ``'zstd'``, ``'xz'``, ``'tar'``} and
            other key-value pairs are forwarded to
            ``zipfile.ZipFile``, ``gzip.GzipFile``,
            ``bz2.BZ2File``, ``zstandard.ZstdCompressor``, ``lzma.LZMAFile`` or
            ``tarfile.TarFile``, respectively.
            As an example, the following could be passed for faster compression and to create
            a reproducible gzip archive:
            ``compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}``.

            May be a dict with key 'method' as compression mode
            and other entries as additional compression options if
            compression mode is 'zip'.

            Passing compression options as keys in dict is
            supported for compression modes 'gzip', 'bz2', 'zstd', and 'zip'.
    quoting : optional constant from csv module
        Defaults to csv.QUOTE_MINIMAL. If you have set a `float_format`
        then floats are converted to strings and thus csv.QUOTE_NONNUMERIC
        will treat them as non-numeric.
    quotechar : str, default '\"'
        String of length 1. Character used to quote fields.
    lineterminator : str, optional
        The newline character or character sequence to use in the output
        file. Defaults to `os.linesep`, which depends on the OS in which
        this method is called ('\\n' for linux, '\\r\\n' for Windows, i.e.).
    chunksize : int or None
        Rows to write at a time.
    date_format : str, default None
        Format string for datetime objects.
    doublequote : bool, default True
        Control quoting of `quotechar` inside a field.
    escapechar : str, default None
        String of length 1. Character used to escape `sep` and `quotechar`
        when appropriate.
    decimal : str, default '.'
        Character recognized as decimal separator. E.g. use ',' for
        European data.
    errors : str, default 'strict'
        Specifies how encoding and decoding errors are to be handled.
        See the errors argument for :func:`open` for a full list
        of options.
    storage_options : dict, optional
        Extra options that make sense for a particular storage connection, e.g.
        host, port, username, password, etc. For HTTP(S) URLs the key-value pairs
        are forwarded to ``urllib.request.Request`` as header options. For other
        URLs (e.g. starting with "s3://", and "gcs://") the key-value pairs are
        forwarded to ``fsspec.open``. Please see ``fsspec`` and ``urllib`` for more
        details, and for more examples on storage options refer `here
        <https://pandas.pydata.org/docs/user_guide/io.html?
        highlight=storage_options#reading-writing-remote-files>`_.

    Returns
    -------
    None or str
        If path_or_buf is None, returns the resulting csv format as a
        string. Otherwise returns None.
    """
    kwargs = {
        "sep": sep,
        "na_rep": na_rep,
        "float_format": float_format,
        "columns": columns,
        "header": header,
        "index": index,
        "index_label": index_label,
        "mode": mode,
        "encoding": encoding,
        "compression": compression,
        "quoting": quoting,
        "quotechar": quotechar,
        "lineterminator": lineterminator,
        "chunksize": chunksize,
        "date_format": date_format,
        "doublequote": doublequote,
        "escapechar": escapechar,
        "decimal": decimal,
        "errors": errors,
        "storage_options": storage_options,
    }

    if not self.attrs:
        return pd.DataFrame.to_csv(self, filepath_or_buffer, **kwargs)  # pyright: ignore[reportCallIssue, reportArgumentType]

    with get_handle(
        filepath_or_buffer,
        mode,
        encoding=encoding,
        errors=errors,
        compression=compression,
        storage_options=storage_options,
    ) as handles:
        # write header here
        if header:
            handles.handle.write(
                ascii_generate_header(
                    self,
                    comments=commentchar,
                    delimiter=sep,
                    commented_header=(sep != ","),
                )
            )
        kwargs["header"] = False
        return pd.DataFrame.to_csv(self, handles.handle, **kwargs)


def to_ascii(
    self: pd.DataFrame,
    filepath_or_buffer: FilePath | BaseBuffer,
    *,
    sep: "str" = " ",
    commentchar: str = "#",
    **kwargs,
) -> str | None:
    r"""
    Write object to an ASCII values file while preserving attrs

    Equivalent to `to_csv` with default `sep` set to a space.

    See also
    --------
    to_csv: Write object to a CSV file while preserving attrs
    """
    return to_csv(self, filepath_or_buffer, sep=sep, commentchar=commentchar, **kwargs)


def from_csv(
    filepath_or_buffer: Union[str, IOBase, PathLike],
    *,
    commented_header: bool = False,
    **kwargs,
):
    r"""
    Read a CSV file into a DataFrame while preserving header information

    Equivalent to `pd.read_csv` with preserved header information.

    Also supports optionally iterating or breaking of the file
    into chunks.

    Additional help can be found in the online docs for
    `IO Tools <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html>`_.

    Parameters
    ----------
    filepath_or_buffer : str, path object or file-like object
        Any valid string path is acceptable.
    commented_header: bool, default False
        Whether the column definition header line starts with a comment character.
    commentchar: str, default '#'
        Character to treat as a comment character.
    sep : str, default ','
        Character or regex pattern to treat as the delimiter. If ``sep=None``, the
        C engine cannot automatically detect the separator, but the Python
        parsing engine can, meaning the latter will be used and automatically
        detect the separator from only the first valid row of the file by
        Python's builtin sniffer tool, ``csv.Sniffer``.
        In addition, separators longer than 1 character and different from
        ``'\s+'`` will be interpreted as regular expressions and will also force
        the use of the Python parsing engine. Note that regex delimiters are prone
        to ignoring quoted data. Regex example: ``'\r\t'``.

    Returns
    -------
    DataFrame: pd.DataFrame
        The parsed data as a pd.DataFrame.
    header : HeaderInfo
        The header information extracted from the file.

    See also
    --------
    pd.read_csv: Read a CSV file into a DataFrame
    """
    kwargs.setdefault("delimiter", ",")
    kwargs.setdefault("comment", "#")

    nlines, hdr, names = ascii_read_header(
        filepath_or_buffer,
        commented_header=commented_header,
        delimiter=kwargs["delimiter"],
        commentchar=kwargs["comment"],
    )
    kwargs.setdefault("names", names)
    kwargs.setdefault("skiprows", nlines)
    kwargs.setdefault("header", None)
    df = pd.read_csv(filepath_or_buffer, **kwargs)  # type: ignore / safe

    return df, hdr


def from_ascii(
    filepath_or_buffer: Union[str, IOBase, PathLike],
    *,
    commented_header: bool = False,
    **kwargs,
):
    """Read an ASCII file into a DataFrame.

    from_csv with delimiter set to " " by default

    Parameters
    ----------
    filepath_or_buffer : str, path object or file-like object
        Any valid string path is acceptable. The string could be a URL. Valid
        URL schemes include http, ftp, s3, and file. For file URLs, a host is
        expected. A local file could be:
        ``file://localhost/path/to/table.csv``.
    commented_header : bool, default False
        Whether the header is commented or not.
    **kwargs : dict
        Additional keyword arguments passed to ``pd.read_csv``.

    Returns
    -------
    DataFrame: pd.DataFrame
        The parsed data as a pd.DataFrame.
    header : HeaderInfo
        The header information extracted from the file.

    See also
    --------
    from_csv: Read a CSV file into a DataFrame
    """
    kwargs.setdefault("delimiter", " ")
    return from_csv(filepath_or_buffer, commented_header=commented_header, **kwargs)
