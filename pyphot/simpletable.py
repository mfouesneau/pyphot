""" This file implements a Table class
    that is designed to be the basis of any format

Requirements
------------

* FIT format:
    * astropy:
        provides a replacement to pyfits
        pyfits can still be used instead but astropy is now the default

* HDF5 format:
    * pytables

RuntimeError will be raised when writing to a format associated with missing
package.


.. code-block::python

    >>> t = SimpleTable('path/mytable.csv')
    # get a subset of columns only
    >>> s = t.get('M_* logTe logLo U B V I J K')
    # set some aliases
    >>> t.set_alias('logT', 'logTe')
    >>> t.set_alias('logL', 'logLLo')
    # make a query on one or multiple column
    >>> q = s.selectWhere('logT logL', '(J > 2) & (10 ** logT > 5000)')
    # q is also a table object
    >>> q.Plotter.plot('logT', 'logL', ',')
    # makes a simple plot (see :module:`plotter`)
    >>> s.write('newtable.fits')
    # export the initial subtable to a new file
"""
from __future__ import (absolute_import, division, print_function)

__version__ = '3.0'
__all__ = ['AstroHelpers', 'AstroTable', 'SimpleTable', 'stats']

import sys
import math
from copy import deepcopy
import re
import itertools
from functools import wraps, partial
import numpy as np
from numpy import deg2rad, rad2deg, sin, cos, sqrt, arcsin, arctan2
from numpy.lib import recfunctions
import types

try:
    from astropy.io import fits as pyfits
except ImportError:
    try:
        import pyfits       
    except ImportError:
        pyfits = None

try:
    import tables
except ImportError:
    tables = None

try:
    import pandas as _pd
except ImportError:
    _pd = None

try:
    from astropy.table import Table as _astropytable
except ImportError:
    _astropytable = None

try:
    from .plotter import Plotter
except ImportError:
    Plotter = None

# ==============================================================================
# Python 3 compatibility behavior
# ==============================================================================
# remap some python 2 built-ins on to py3k behavior or equivalent
# Most of them become generators
import operator

PY3 = sys.version_info[0] > 2

if PY3:
    iteritems = operator.methodcaller('items')
    itervalues = operator.methodcaller('values')
    basestring = (str, bytes)
else:
    range = xrange
    from itertools import izip as zip
    iteritems = operator.methodcaller('iteritems')
    itervalues = operator.methodcaller('itervalues')
    basestring = (str, unicode)


# ==============================================================================
# Specials -- special functions
# ==============================================================================

def pretty_size_print(num_bytes):
    """
    Output number of bytes in a human readable format

    keywords
    --------
    num_bytes: int
        number of bytes to convert

    returns
    -------
    output: str
        string representation of the size with appropriate unit scale
    """
    if num_bytes is None:
        return

    KiB = 1024
    MiB = KiB * KiB
    GiB = KiB * MiB
    TiB = KiB * GiB
    PiB = KiB * TiB
    EiB = KiB * PiB
    ZiB = KiB * EiB
    YiB = KiB * ZiB

    if num_bytes > YiB:
        output = '%.3g YB' % (num_bytes / YiB)
    elif num_bytes > ZiB:
        output = '%.3g ZB' % (num_bytes / ZiB)
    elif num_bytes > EiB:
        output = '%.3g EB' % (num_bytes / EiB)
    elif num_bytes > PiB:
        output = '%.3g PB' % (num_bytes / PiB)
    elif num_bytes > TiB:
        output = '%.3g TB' % (num_bytes / TiB)
    elif num_bytes > GiB:
        output = '%.3g GB' % (num_bytes / GiB)
    elif num_bytes > MiB:
        output = '%.3g MB' % (num_bytes / MiB)
    elif num_bytes > KiB:
        output = '%.3g KB' % (num_bytes / KiB)
    else:
        output = '%.3g Bytes' % (num_bytes)

    return output


def _fits_read_header(hdr):
    """
    Convert pyfits header into dictionary with relevant values

    Parameters
    ----------

    hdr: pyftis.Header
        fits unit

    Returns
    -------
    header: dict
        header dictionary

    alias: dict
        aliases

    units: dict
        units

    comments: dict
        comments/description of keywords
    """
    header = {}
    alias = {}
    units = {}
    comments = {}

    # generic cards
    genTerms = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1',
                'NAXIS2', 'PCOUNT', 'GCOUNT', 'TFIELDS',
                'EXTNAME']
    fieldTerms = ['TTYPE', 'TFORM', 'TUNIT', 'ALIAS']

    # read col comments
    # for k, name, comment in hdr.ascard['TTYPE*']:
    try:
        for card in hdr.cards['TTYPE*']:
            name = card.value
            comments[name] = card.comment
            u = hdr.get(card.keyword.replace('TYPE', 'UNIT'), None)
            if u is not None:
                units[name] = u

        # for k, val, _ in hdr.ascard['ALIAS*']:
        for card in hdr.cards['ALIAS*']:
            k = card.keyword
            val = card.value
            al, orig = val.split('=')
            alias[al] = orig
    except:   #pyfits stsci
        for card in hdr.ascard['TTYPE*']:
            name = card.value
            comments[name] = card.comment
            u = hdr.get(card.key.replace('TYPE', 'UNIT'), None)
            if u is not None:
                units[name] = u

        # for k, val, _ in hdr.ascard['ALIAS*']:
        for card in hdr.ascard['ALIAS*']:
            k = card.key
            val = card.value
            al, orig = val.split('=')
            alias[al] = orig

    # other specific keywords: COMMENT, HISTORY
    header_comments = []
    header_history = []
    for k, v in hdr.items():
        if (k not in genTerms) and (k[:5] not in fieldTerms):
            if (k == 'COMMENT'):
                header_comments.append(v)
            elif (k == 'HISTORY'):
                header_history.append(v)
            else:
                header[k] = v

    # COMMENT, HISTORY polish
    if len(header_comments) > 0:
        header['COMMENT'] = '\n'.join(header_comments)
    if len(header_history) > 0:
        header['HISTORY'] = '\n'.join(header_history)

    if 'EXTNAME' in hdr:
        header['NAME'] = hdr['EXTNAME']

    return header, alias, units, comments


def _fits_generate_header(tab):
    """ Generate the corresponding fits Header that contains all necessary info

    Parameters
    ----------

    tab: SimpleTable instance
        table

    Returns
    -------
    hdr: pyfits.Header
        header instance
    """
    # get column cards

    cards = []

    # names units and comments
    for e, k in enumerate(tab.keys()):
        cards.append(('TTYPE{0:d}'.format(e + 1), k, tab._desc.get(k, '')))
        u = tab._units.get(k, '')
        if u not in ['', 'None', None]:
            cards.append(('TUNIT{0:d}'.format(e + 1), tab._units.get(k, ''),
                          'unit of {0:s}'.format(k)))

    # add aliases
    for e, v in enumerate(tab._aliases.items()):
        cards.append( ('ALIAS{0:d}'.format(e + 1), '='.join(v), '') )

    if tab.header['NAME'] not in ['', 'None', None, 'No Name']:
        cards.append(('EXTNAME', tab.header['NAME'], ''))

    hdr = pyfits.Header(cards)

    for k, v in tab.header.items():
        if (v not in ['', 'None', None]) & (k != 'NAME'):
            if (k != 'COMMENT') & (k != 'HISTORY'):
                hdr.update(k, v)
            else:
                txt = v.split('\n')
                for j in txt:
                    if k == 'COMMENT':
                        hdr.add_comment(j)
                    elif k == 'HISTORY':
                        hdr.add_history(j)
    return hdr


def _fits_writeto(filename, data, header=None, output_verify='exception',
                  clobber=False, checksum=False):
    """
    Create a new FITS file using the supplied data/header.
    Patched version of pyfits to correctly include provided header

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to write to.  If opened, must be opened in a writeable binary
        mode such as 'wb' or 'ab+'.

    data : array, record array, or groups data object
        data to write to the new file

    header : `Header` object, optional
        the header associated with ``data``. If `None`, a header
        of the appropriate type is created for the supplied data. This
        argument is optional.

    output_verify : str
        Output verification option.  Must be one of ``"fix"``, ``"silentfix"``,
        ``"ignore"``, ``"warn"``, or ``"exception"``.  May also be any
        combination of ``"fix"`` or ``"silentfix"`` with ``"+ignore"``,
        ``+warn``, or ``+exception" (e.g. ``"fix+warn"``).  See :ref:`verify`
        for more info.

    clobber : bool, optional
        If `True`, and if filename already exists, it will overwrite
        the file.  Default is `False`.

    checksum : bool, optional
        If `True`, adds both ``DATASUM`` and ``CHECKSUM`` cards to the
        headers of all HDU's written to the file
    """

    hdu = pyfits.convenience._makehdu(data, header)
    hdu.header.update(header.cards)
    if hdu.is_image and not isinstance(hdu, pyfits.PrimaryHDU):
        hdu = pyfits.PrimaryHDU(data, header=header)
    hdu.writeto(filename, clobber=clobber, output_verify=output_verify,
                checksum=checksum)


def _fits_append(filename, data, header=None, checksum=False, verify=True,
                 **kwargs):
    """
    Append the header/data to FITS file if filename exists, create if not.

    If only ``data`` is supplied, a minimal header is created.
    Patched version of pyfits to correctly include provided header

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to write to.  If opened, must be opened for update (rb+) unless it
        is a new file, then it must be opened for append (ab+).  A file or
        `~gzip.GzipFile` object opened for update will be closed after return.

    data : array, table, or group data object
        the new data used for appending

    header : `Header` object, optional
        The header associated with ``data``.  If `None`, an appropriate header
        will be created for the data object supplied.

    checksum : bool, optional
        When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards to the header
        of the HDU when written to the file.

    verify : bool, optional
        When `True`, the existing FITS file will be read in to verify it for
        correctness before appending.  When `False`, content is simply appended
        to the end of the file.  Setting ``verify`` to `False` can be much
        faster.

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.
    """

    name, closed, noexist_or_empty = pyfits.convenience._stat_filename_or_fileobj(filename)

    if noexist_or_empty:
        #
        # The input file or file like object either doesn't exits or is
        # empty.  Use the writeto convenience function to write the
        # output to the empty object.
        #
        _fits_writeto(filename, data, header, checksum=checksum, **kwargs)
    else:
        hdu = pyfits.convenience._makehdu(data, header)
        hdu.header.update(header.cards)

        if isinstance(hdu, pyfits.PrimaryHDU):
            hdu = pyfits.ImageHDU(data, header)

        if verify or not closed:
            f = pyfits.convenience.fitsopen(filename, mode='append')
            f.append(hdu)

            # Set a flag in the HDU so that only this HDU gets a checksum when
            # writing the file.
            hdu._output_checksum = checksum
            f.close(closed=closed)
        else:
            f = pyfits.convenience._File(filename, mode='append')
            hdu._output_checksum = checksum
            hdu._writeto(f)
            f.close()


def _ascii_read_header(fname, comments='#', delimiter=None, commentedHeader=True,
                       *args, **kwargs):
    """
    Read ASCII/CSV header

    Parameters
    ----------
    fname: str or stream
        File, filename, or generator to read.
        Note that generators should return byte strings for Python 3k.

    comments: str, optional
        The character used to indicate the start of a comment;
        default: '#'.

    delimiter: str, optional
        The string used to separate values.  By default, this is any
        whitespace.

    commentedHeader: bool, optional
        if set, the last line of the header is expected to be the column titles

    Returns
    -------
    nlines: int
        number of lines from the header

    header: dict
        header dictionary

    alias: dict
        aliases

    units: dict
        units

    comments: dict
        comments/description of keywords

    names: sequence
        sequence or str, first data line after header, expected to be the column
        names.
    """
    if hasattr(fname, 'read'):
        stream = fname
    else:
        stream = open(fname, 'r')

    header = {}
    alias = {}
    units = {}
    desc = {}

    def parseStrNone(v):
        """ robust parse """
        _v = v.split()
        if (len(_v) == 0):
            return None
        else:
            _v = ' '.join(_v)
            if (_v.lower()) == 'none' or (_v.lower() == 'null'):
                return None
            else:
                return _v

    done = False
    oldline = None
    lasthdr = None
    nlines = 0
    header.setdefault('COMMENT', '')
    header.setdefault('HISTORY', '')
    while done is False:
        line = stream.readline()[:-1]  # getting rid of '\n'
        nlines += 1
        if (line[0] == comments):  # header part
            if (len(line) > 2):
                if line[1] == comments:  # column meta data
                    # column meta is expected to start with ##
                    k = line[2:].split('\t')
                    colname = k[0].strip()
                    colunit = None
                    colcomm = None
                    if len(k) > 1:
                        colunit = parseStrNone(k[1])
                    if len(k) > 2:
                        colcomm = parseStrNone(k[2])

                    if colunit is not None:
                        units[colname] = colunit
                    if colcomm is not None:
                        desc[colname] = colcomm
                else:
                    # header is expected as "# key \t value"
                    k = line[1:].split('\t')
                    if len(k) > 1:
                        key = k[0].strip()  # remove trainling spaces
                        val = ' '.join(k[1:]).strip()

                        if key in ('', None, 'None', 'NONE', 'COMMENT'):
                            header['COMMENT'] = header['COMMENT'] + '\n' + val
                        if key in ('HISTORY', ):
                            header['HISTORY'] = header['HISTORY'] + '\n' + val
                        elif 'alias' in key.lower():
                            # take care of aliases
                            al, orig = val.split('=')
                            alias[al] = orig
                        else:
                            header[key] = val
                        lasthdr = key
                    else:
                        header['COMMENT'] = header['COMMENT'] + '\n' + line[1:]
        else:
            done = True
            if commentedHeader and (oldline is not None):
                names = oldline.split(delimiter)
                nlines -= 1
                if lasthdr == names[0]:
                    header.pop(lasthdr)
            else:
                names = line.split(delimiter)
        oldline = line[1:]

    if not hasattr(fname, 'read'):
        stream.close()
    else:
        stream.seek(stream.tell() - len(line))
        nlines = 0  # make sure the value is set to the current position

    return nlines, header, units, desc, alias, names


def _hdf5_write_data(filename, data, tablename=None, mode='w', append=False,
                     header={}, units={}, comments={}, aliases={}, **kwargs):
    """ Write table into HDF format

    Parameters
    ----------
    filename : file path, or tables.File instance
        File to write to.  If opened, must be opened and writable (mode='w' or 'a')

    data: recarray
        data to write to the new file

    tablename: str
        path of the node including table's name

    mode: str
        in ('w', 'a') mode to open the file

    append: bool
        if set, tends to append data to an existing table

    header: dict
        table header

    units: dict
        dictionary of units

    alias: dict
        aliases

    comments: dict
        comments/description of keywords

    .. note::
        other keywords are forwarded to :func:`tables.open_file`
    """

    if hasattr(filename, 'read'):
        raise Exception("HDF backend does not implement stream")

    if append is True:
        mode = 'a'
    silent = kwargs.pop('silent', False)

    if isinstance(filename, tables.File):
        if (filename.mode != mode) & (mode != 'r'):
            raise tables.FileModeError('The file is already opened in a different mode')
        hd5 = filename
    else:
        hd5 = tables.open_file(filename, mode=mode)

    # check table name and path
    tablename = tablename or header.get('NAME', None)
    if tablename in ('', None, 'Noname', 'None'):
        tablename = '/data'

    w = tablename.split('/')
    where = '/'.join(w[:-1])
    name = w[-1]
    if where in ('', None):
        where = '/'
    if where[0] != '/':
        where = '/' + where

    if append:
        try:
            t = hd5.get_node(where + name)
            t.append(data.astype(t.description._v_dtype))
            t.flush()
        except tables.NoSuchNodeError:
            if not silent:
                print(("Warning: Table {0} does not exists.  \n A new table will be created").format(where + '/' + name))
            append = False

    if not append:
        # t = hd5.createTable(where, name, data, **kwargs)
        t = hd5.create_table(where, name, data, **kwargs)

        # update header
        for k, v in header.items():
            if (k == 'FILTERS') & (float(t.attrs['VERSION']) >= 2.0):
                t.attrs[k.lower()] = v
            else:
                t.attrs[k] = v
        if 'TITLE' not in header:
            t.attrs['TITLE'] = name

        # add column descriptions and units
        for e, colname in enumerate(data.dtype.names):
            _u = units.get(colname, None)
            _d = comments.get(colname, None)
            if _u is not None:
                t.attrs['FIELD_{0:d}_UNIT'] = _u
            if _d is not None:
                t.attrs['FIELD_{0:d}_DESC'] = _d

        # add aliases
        for i, (k, v) in enumerate(aliases.items()):
            t.attrs['ALIAS{0:d}'.format(i)] = '{0:s}={1:s}'.format(k, v)

        t.flush()

    if not isinstance(filename, tables.File):
        hd5.flush()
        hd5.close()


def _hdf5_read_data(filename, tablename=None, silent=False, *args, **kwargs):
    """ Generate the corresponding ascii Header that contains all necessary info

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
    source = tables.open_file(filename, *args, **kwargs)

    if tablename is None:
        node = source.listNodes('/')[0]
        tablename = node.name
    else:
        if tablename[0] != '/':
            node = source.get_node('/' + tablename)
        else:
            node = source.get_node(tablename)
    if not silent:
        print("\tLoading table: {0}".format(tablename))

    hdr = {}
    aliases = {}

    # read header
    exclude = ['NROWS', 'VERSION', 'CLASS', 'EXTNAME', 'TITLE']
    for k in node.attrs._v_attrnames:
        if (k not in exclude):
            if (k[:5] != 'FIELD') & (k[:5] != 'ALIAS'):
                hdr[k] = node.attrs[k]
            elif k[:5] == 'ALIAS':
                c0, c1 = node.attrs[k].split('=')
                aliases[c0] = c1

    empty_name = ['', 'None', 'Noname', None]
    if node.attrs['TITLE'] not in empty_name:
        hdr['NAME'] = node.attrs['TITLE']
    else:
        hdr['NAME'] = '{0:s}/{1:s}'.format(filename, node.name)

    # read column meta
    units = {}
    desc = {}

    for (k, colname) in enumerate(node.colnames):
        _u = getattr(node.attrs, 'FIELD_{0:d}_UNIT'.format(k), None)
        _d = getattr(node.attrs, 'FIELD_{0:d}_DESC'.format(k), None)
        if _u is not None:
            units[colname] = _u
        if _d is not None:
            desc[colname] = _d

    data = node[:]

    source.close()

    return hdr, aliases, units, desc, data


def _ascii_generate_header(tab, comments='#', delimiter=' ',
                           commentedHeader=True):
    """ Generate the corresponding ascii Header that contains all necessary info

    Parameters
    ----------

    tab: SimpleTable instance
        table

    comments: str
        string to prepend header lines

    delimiter: str, optional
        The string used to separate values.  By default, this is any
        whitespace.

    commentedHeader: bool, optional
        if set, the last line of the header is expected to be the column titles

    Returns
    -------
    hdr: str
        string that will be be written at the beginning of the file
    """
    hdr = []

    if comments is None:
        comments = ''

    # table header
    length = max(map(len, tab.header.keys()))
    fmt = '{{0:s}} {{1:{0:d}s}}\t{{2:s}}'.format(length)
    for k, v in tab.header.items():
        for vk in v.split('\n'):
            if len(vk) > 0:
                hdr.append(fmt.format(comments, k.upper(), vk.strip()))

    # column metadata
    hdr.append(comments)  # add empty line
    length = max(map(len, tab.keys()))
    fmt = '{{0:s}}{{0:s}} {{1:{0:d}s}}\t{{2:s}}\t{{3:s}}'.format(length)
    for colname in tab.keys():
        unit = tab._units.get(colname, 'None')
        desc = tab._desc.get(colname, 'None')
        hdr.append(fmt.format(comments, colname, unit, desc))

    # aliases
    if len(tab._aliases) > 0:
        hdr.append(comments)  # add empty line
        for k, v in tab._aliases.items():
            hdr.append('{0:s} alias\t{1:s}={2:s}'.format(comments, k, v))

    # column names
    hdr.append(comments)
    if commentedHeader:
        hdr.append('{0:s} {1:s}'.format(comments, delimiter.join(tab.keys())))
    else:
        hdr.append('{0:s}'.format(delimiter.join(tab.keys())))

    return '\n'.join(hdr)


def _latex_writeto(filename, tab, comments='%'):
    """ Write the data into a latex table format

    Parameters
    ----------
    filename: str
        file or unit to write into

    tab: SimpleTable instance
        table

    comments: str
        string to prepend header lines

    delimiter: str, optional
        The string used to separate values.  By default, this is any
        whitespace.

    commentedHeader: bool, optional
        if set, the last line of the header is expected to be the column titles
    """
    txt = "\\begin{table}\n\\begin{center}\n"

    # add caption
    tabname = tab.header.get('NAME', None)
    if tabname not in ['', None, 'None']:
        txt += "\\caption{{{0:s}}}\n".format(tabname)

    # tabular
    txt += '\\begin{{tabular}}{{{0:s}}}\n'.format('c' * tab.ncols)
    txt += tab.pprint(delim=' & ', fields='MAG*', headerChar='', endline='\\\\\n', all=True, ret=True)
    txt += '\\end{tabular}\n'

    # end table
    txt += "\\end{center}\n"

    # add notes if any
    if len(tab._desc) > 0:
        txt += '\% notes \n\\begin{scriptsize}\n'
        for e, (k, v) in enumerate(tab._desc.items()):
            if v not in (None, 'None', 'none', ''):
                txt += '{0:d} {1:s}: {2:s} \\\\\n'.format(e, k, v)
        txt += '\\end{scriptsize}\n'
    txt += "\\end{table}\n"
    if hasattr(filename, 'write'):
        filename.write(txt)
    else:
        with open(filename, 'w') as unit:
            unit.write(txt)


def _convert_dict_to_structured_ndarray(data):
    """convert_dict_to_structured_ndarray

    Parameters
    ----------

    data: dictionary like object
        data structure which provides iteritems and itervalues

    returns
    -------
    tab: structured ndarray
        structured numpy array
    """
    newdtype = []
    try:
        for key, dk in iteritems(data):
            _dk = np.asarray(dk)
            dtype = _dk.dtype
            # unknown type is converted to text
            if dtype.type == np.object_:
                if len(data) == 0:
                    longest = 0
                else:
                    longest = len(max(_dk, key=len))
                    _dk = _dk.astype('|%iS' % longest)
            if _dk.ndim > 1:
                newdtype.append((str(key), _dk.dtype, (_dk.shape[1],)))
            else:
                newdtype.append((str(key), _dk.dtype))
        tab = np.rec.fromarrays(itervalues(data), dtype=newdtype)
    except AttributeError:  # not a dict
        # hope it's a tuple ((key, value),) pairs.
        from itertools import tee
        d1, d2 = tee(data)
        for key, dk in d1:
            _dk = np.asarray(dk)
            dtype = _dk.dtype
            # unknown type is converted to text
            if dtype.type == np.object_:
                if len(data) == 0:
                    longest = 0
                else:
                    longest = len(max(_dk, key=len))
                    _dk = _dk.astype('|%iS' % longest)
            if _dk.ndim > 1:
                newdtype.append((str(key), _dk.dtype, (_dk.shape[1],)))
            else:
                newdtype.append((str(key), _dk.dtype))
        tab = np.rec.fromarrays((dk for (_, dk) in d2), dtype=newdtype)

    return tab


def __indent__(rows, header=None, units=None, headerChar='-',
               delim=' | ', endline='\n', **kwargs):
    """Indents a table by column.

    Parameters
    ----------
    rows: sequences of rows
        one sequence per row.

    header: sequence of str
        row consists of the columns' names

    units: sequence of str
        Sequence of units

    headerChar: char
        Character to be used for the row separator line

    delim: char
        The column delimiter.

    returns
    -------
    txt: str
        string represation of rows
    """
    length_data = list(map(max, zip(*[list(map(len, k)) for k in rows])))
    length = length_data[:]

    if (header is not None):
        length_header = list(map(len, header))
        length = list(map(max, zip(length_data, length_header)))

    if (units is not None):
        length_units = list(map(len, units))
        length = list(map(max, zip(length_data, length_units)))

    if headerChar not in (None, '', ' '):
        rowSeparator = headerChar * (sum(length) + len(delim) * (len(length) - 1)) + endline
    else:
        rowSeparator = ''

    # make the format
    fmt = ['{{{0:d}:{1:d}s}}'.format(k, l) for (k, l) in enumerate(length)]
    fmt = delim.join(fmt) + endline
    # write the string
    txt = rowSeparator
    if header is not None:
        txt += fmt.format(*header)  # + endline
        txt += rowSeparator
    if units is not None:
        txt += fmt.format(*units)  # + endline
        txt += rowSeparator
    for r in rows:
        txt += fmt.format(*r)  # + endline
    txt += rowSeparator
    return txt


def pprint_rec_entry(data, num=0, keys=None):
        """ print one line with key and values properly to be readable

        Parameters
        ----------
        data: recarray
            data to extract entry from

        num: int, slice
            indice selection

        keys: sequence or str
            if str, can be a regular expression
            if sequence, the sequence of keys to print
        """
        if (keys is None) or (keys == '*'):
            _keys = data.dtype.names
        elif type(keys) in basestring:
            _keys = [k for k in data.dtype.names if (re.match(keys, k) is not None)]
        else:
            _keys = keys

        length = max(map(len, _keys))
        fmt = '{{0:{0:d}s}}: {{1}}'.format(length)
        data = data[num]

        for k in _keys:
            print(fmt.format(k, data[k]))


def pprint_rec_array(data, idx=None, fields=None, ret=False, all=False,
                     headerChar='-', delim=' | ', endline='\n' ):
        """ Pretty print the table content
            you can select the table parts to display using idx to
            select the rows and fields to only display some columns
            (ret is only for insternal use)

        Parameters
        ----------
        data: array
            array to show

        idx: sequence, slide
            sub selection to print

        fields: str, sequence
            if str can be a regular expression, and/or list of fields separated
            by spaces or commas

        ret: bool
            if set return the string representation instead of printing the result

        all: bool
            if set, force to show all rows

        headerChar: char
            Character to be used for the row separator line

        delim: char
            The column delimiter.
        """
        if (fields is None) or (fields == '*'):
            _keys = data.dtype.names
        elif type(fields) in basestring:
            if ',' in fields:
                _fields = fields.split(',')
            elif ' ' in fields:
                _fields = fields.split()
            else:
                _fields = [fields]
            lbls = data.dtype.names
            _keys = []
            for _fk in _fields:
                _keys += [k for k in lbls if (re.match(_fk, k) is not None)]
        else:
            lbls = data.dtype.names
            _keys = []
            for _fk in _fields:
                _keys += [k for k in lbls if (re.match(_fk, k) is not None)]

        nfields = len(_keys)
        nrows = len(data)
        fields = list(_keys)

        if idx is None:
            if (nrows < 10) or (all is True):
                rows = [ [ str(data[k][rk]) for k in _keys ] for rk in range(nrows)]
            else:
                _idx = range(6)
                rows = [ [ str(data[k][rk]) for k in _keys ] for rk in range(5) ]
                if nfields > 1:
                    rows += [ ['...' for k in range(nfields) ] ]
                else:
                    rows += [ ['...' for k in range(nfields) ] ]
                rows += [ [ str(data[k][rk]) for k in fields ] for rk in range(-5, 0)]
        elif isinstance(idx, slice):
            _idx = range(idx.start, idx.stop, idx.step or 1)
            rows = [ [ str(data[k][rk]) for k in fields ] for rk in _idx]
        else:
            rows = [ [ str(data[k][rk]) for k in fields ] for rk in idx]

        out = __indent__(rows, header=_keys, units=None, delim=delim,
                         headerChar=headerChar, endline=endline)
        if ret is True:
            return out
        else:
            print(out)


def elementwise(func):
    """
    Quick and dirty elementwise function decorator it provides a quick way
    to apply a function either on one element or a sequence of elements
    """
    @wraps(func)
    def wrapper(it, **kwargs):
        if hasattr(it, '__iter__') & (type(it) not in basestring):
            _f = partial(func, **kwargs)
            return map(_f, it)
        else:
            return func(it, **kwargs)
    return wrapper


class AstroHelpers(object):
    """ Helpers related to astronomy data """

    @staticmethod
    @elementwise
    def hms2deg(_str, delim=':'):
        """ Convert hex coordinates into degrees

        Parameters
        ----------
        str: string or sequence
            string to convert

        delimiter: str
            character delimiting the fields

        Returns
        -------
        deg: float
            angle in degrees
        """
        if _str[0] == '-':
            neg = -1
            _str = _str[1:]
        else:
            neg = 1
        _str = _str.split(delim)
        return neg * ((((float(_str[-1]) / 60. +
                         float(_str[1])) / 60. +
                        float(_str[0])) / 24. * 360.))

    @staticmethod
    @elementwise
    def deg2dms(val, delim=':'):
        """ Convert degrees into hex coordinates

        Parameters
        ----------
        deg: float
            angle in degrees

        delimiter: str
            character delimiting the fields

        Returns
        -------
        str: string or sequence
            string to convert
        """
        if val < 0:
            sign = -1
        else:
            sign = 1
        d = int( sign * val )
        m = int( (sign * val - d) * 60. )
        s = (( sign * val - d) * 60.  - m) * 60.
        return '{0}{1}{2}{3}{4}'.format( sign * d, delim, m, delim, s)

    @staticmethod
    @elementwise
    def deg2hms(val, delim=':'):
        """ Convert degrees into hex coordinates

        Parameters
        ----------
        deg: float
            angle in degrees

        delimiter: str
            character delimiting the fields

        Returns
        -------
        str: string or sequence
            string to convert
        """
        if val < 0:
            sign = -1
        else:
            sign = 1
        h = int( sign * val / 45. * 3.)   # * 24 / 360
        m = int( (sign * val / 45. * 3. - h) * 60. )
        s = (( sign * val / 45. * 3. - h) * 60.  - m) * 60.
        return '{0}{1}{2}{3}{4}'.format( sign * h, delim, m, delim, s)

    @staticmethod
    @elementwise
    def dms2deg(_str, delim=':'):
        """ Convert hex coordinates into degrees

        Parameters
        ----------
        str: string or sequence
            string to convert

        delimiter: str
            character delimiting the fields

        Returns
        -------
        deg: float
            angle in degrees
        """
        if _str[0] == '-':
            neg = -1
            _str = _str[1:]
        else:
            neg = 1
        _str = _str.split(delim)
        return (neg * ((float(_str[-1]) / 60. + float(_str[1])) / 60. + float(_str[0])))

    @staticmethod
    @elementwise
    def euler(ai_in, bi_in, select, b1950=False, dtype='f8'):
        """
        Transform between Galactic, celestial, and ecliptic coordinates.
        Celestial coordinates (RA, Dec) should be given in equinox J2000
        unless the b1950 is True.

        +-------+--------------+------------+----------+----------+-----------+
        |select | From         | To         |   select |   From   |  To       |
        +-------+--------------+------------+----------+----------+-----------+
        |1      |RA-Dec (2000) | Galactic   |     4    | Ecliptic |  RA-Dec   |
        +-------+--------------+------------+----------+----------+-----------+
        |2      |Galactic      | RA-DEC     |     5    | Ecliptic |  Galactic |
        +-------+--------------+------------+----------+----------+-----------+
        |3      |RA-Dec        | Ecliptic   |     6    | Galactic |  Ecliptic |
        +-------+--------------+------------+----------+----------+-----------+

        Parameters
        ----------

        long_in: float, or sequence
            Input Longitude in DEGREES, scalar or vector.

        lat_in: float, or sequence
            Latitude in DEGREES

        select: int
            Integer from 1 to 6 specifying type of coordinate transformation.

        b1950: bool
            set equinox set to 1950


        Returns
        -------
        long_out: float, seq
            Output Longitude in DEGREES

        lat_out: float, seq
            Output Latitude in DEGREES


        .. note::

            Written W. Landsman,  February 1987
            Adapted from Fortran by Daryl Yentis NRL
            Converted to IDL V5.0   W. Landsman   September 1997
            Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
            Add option to specify SELECT as a keyword W. Landsman March 2003
            Converted from IDL to numerical Python: Erin Sheldon, NYU, 2008-07-02
        """

        # Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
        ai = np.array(ai_in, ndmin=1, copy=True, dtype=dtype)
        bi = np.array(bi_in, ndmin=1, copy=True, dtype=dtype)

        PI = math.pi
        # HALFPI = PI / 2.0
        D2R = PI / 180.0
        R2D = 1.0 / D2R

        twopi   = 2.0 * PI
        fourpi  = 4.0 * PI

        #   J2000 coordinate conversions are based on the following constants
        #   (see the Hipparcos explanatory supplement).
        #  eps = 23.4392911111d           Obliquity of the ecliptic
        #  alphaG = 192.85948d            Right Ascension of Galactic North Pole
        #  deltaG = 27.12825d             Declination of Galactic North Pole
        #  lomega = 32.93192d             Galactic longitude of celestial equator
        #  alphaE = 180.02322d            Ecliptic longitude of Galactic North Pole
        #  deltaE = 29.811438523d         Ecliptic latitude of Galactic North Pole
        #  Eomega  = 6.3839743d           Galactic longitude of ecliptic equator
        # Parameters for all the different conversions
        if b1950:
            # equinox = '(B1950)'
            psi    = np.array([ 0.57595865315, 4.9261918136,
                                0.00000000000, 0.0000000000,
                                0.11129056012, 4.7005372834], dtype=dtype)
            stheta = np.array([ 0.88781538514, -0.88781538514,
                                0.39788119938, -0.39788119938,
                                0.86766174755, -0.86766174755], dtype=dtype)
            ctheta = np.array([ 0.46019978478, 0.46019978478,
                                0.91743694670, 0.91743694670,
                                0.49715499774, 0.49715499774], dtype=dtype)
            phi    = np.array([ 4.9261918136,  0.57595865315,
                                0.0000000000, 0.00000000000,
                                4.7005372834, 0.11129056012], dtype=dtype)
        else:
            # equinox = '(J2000)'
            psi    = np.array([ 0.57477043300, 4.9368292465,
                                0.00000000000, 0.0000000000,
                                0.11142137093, 4.71279419371], dtype=dtype)
            stheta = np.array([ 0.88998808748, -0.88998808748,
                                0.39777715593, -0.39777715593,
                                0.86766622025, -0.86766622025], dtype=dtype)
            ctheta = np.array([ 0.45598377618, 0.45598377618,
                                0.91748206207, 0.91748206207,
                                0.49714719172, 0.49714719172], dtype=dtype)
            phi    = np.array([ 4.9368292465,  0.57477043300,
                                0.0000000000, 0.00000000000,
                                4.71279419371, 0.11142137093], dtype=dtype)

        # zero offset
        i  = select - 1
        a  = ai * D2R - phi[i]

        b = bi * D2R
        sb = sin(b)
        cb = cos(b)
        cbsa = cb * sin(a)
        b  = -stheta[i] * cbsa + ctheta[i] * sb
        w, = np.where(b > 1.0)
        if w.size > 0:
            b[w] = 1.0
        bo = arcsin(b) * R2D
        a  = arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a) )
        ao = ( (a + psi[i] + fourpi) % twopi) * R2D
        return ao, bo

    @staticmethod
    def sphdist(ra1, dec1, ra2, dec2):
        """measures the spherical distance between 2 points

        Parameters
        ----------
        ra1: float or sequence
            first right ascensions in degrees

        dec1: float or sequence
            first declination in degrees
        ra2: float or sequence
            second right ascensions in degrees
        dec2: float or sequence
            first declination in degrees

        Returns
        -------
        Outputs: float or sequence
            returns a distance in degrees
        """
        dec1_r = deg2rad(dec1)
        dec2_r = deg2rad(dec2)
        return 2. * rad2deg(arcsin(sqrt((sin((dec1_r - dec2_r) / 2)) ** 2 +
                                        cos(dec1_r) * cos(dec2_r) * (
                                            sin((deg2rad(ra1 - ra2)) / 2)) **
                                        2)))

    @staticmethod
    def conesearch(ra0, dec0, ra, dec, r, outtype=0):
        """ Perform a cone search on a table

        Parameters
        ----------
        ra0: ndarray[ndim=1, dtype=float]
            column name to use as RA source in degrees

        dec0: ndarray[ndim=1, dtype=float]
            column name to use as DEC source in degrees

        ra: float
            ra to look for (in degree)

        dec: float
            ra to look for (in degree)

        r: float
            distance in degrees

        outtype: int
            type of outputs
                0 -- minimal, indices of matching coordinates
                1 -- indices and distances of matching coordinates
                2 -- full, boolean filter and distances

        Returns
        -------
        t: tuple
            if outtype is 0:
                only return indices from ra0, dec0
            elif outtype is 1:
                return indices from ra0, dec0 and distances
            elif outtype is 2:
                return conditional vector and distance to all ra0, dec0
        """
        @elementwise
        def getDist( pk ):
            """ get spherical distance between 2 points """
            return AstroHelpers.sphdist(pk[0], pk[1], ra, dec)

        dist = np.array(list(getDist(zip(ra0, dec0))))
        v = (dist <= r)

        if outtype == 0:
            return np.ravel(np.where(v))
        elif outtype == 1:
            return np.ravel(np.where(v)), dist[v]
        else:
            return v, dist


# ==============================================================================
# SimpleTable -- provides table manipulations with limited storage formats
# ==============================================================================
class SimpleTable(object):
    """ Table class that is designed to be the basis of any format wrapping
    around numpy recarrays

    Attributes
    ----------

    fname: str or object
        if str, the file to read from. This may be limited to the format
        currently handled automatically. If the format is not correctly handled,
        you can try by providing an object.__

        if object with a structure like dict, ndarray, or recarray-like
            the data will be encapsulated into a Table

    caseless: bool
        if set, column names will be caseless during operations

    aliases: dict
        set of column aliases (can be defined later :func:`set_alias`)

    units: dict
        set of column units (can be defined later :func:`set_unit`)

    desc: dict
        set of column description or comments (can be defined later :func:`set_comment`)

    header: dict
        key, value pair corresponding to the attributes of the table
    """

    def __init__(self, fname, *args, **kwargs):

        dtype = kwargs.pop('dtype', None)
        dtype = kwargs.pop('format', dtype)
        self.caseless = kwargs.get('caseless', False)
        self._aliases = kwargs.get('aliases', {})
        self._units = kwargs.get('units', {})
        self._desc = kwargs.get('desc', {})

        if (isinstance(fname, (dict, tuple, list, types.GeneratorType))) or (dtype in [dict, 'dict']):
            try:
                self.header = fname.pop('header', {})
            except (AttributeError, TypeError):
                self.header = kwargs.pop('header', {})
            self.data = _convert_dict_to_structured_ndarray(fname)
        elif (type(fname) in basestring) or (dtype is not None):
            if (type(fname) in basestring):
                extension = fname.split('.')[-1]
            else:
                extension = None
            if (extension == 'csv') or dtype == 'csv':
                kwargs.setdefault('delimiter', ',')
                commentedHeader = kwargs.pop('commentedHeader', False)
                n, header, units, comments, aliases, names = _ascii_read_header(fname, commentedHeader=commentedHeader, **kwargs)
                if 'names' in kwargs:
                    n -= 1
                kwargs.setdefault('names', names)
                if _pd is not None:   # pandas is faster
                    kwargs.setdefault('comment', '#')
                    kwargs.setdefault('skiprows', n)
                    self.data = _pd.read_csv(fname, *args, **kwargs).to_records()
                else:
                    kwargs.setdefault('skip_header', n)
                    kwargs.setdefault('comments', '#')
                    self.data = np.recfromcsv(fname, *args, **kwargs)
                self.header = header
                self._units.update(**units)
                self._desc.update(**comments)
                self._aliases.update(**aliases)
                kwargs.setdefault('names', True)
            elif (extension in ('tsv', 'dat', 'txt')) or (dtype in ('tsv', 'dat', 'txt')):
                commentedHeader = kwargs.pop('commentedHeader', True)
                n, header, units, comments, aliases, names = _ascii_read_header(fname, commentedHeader=commentedHeader, **kwargs)
                kwargs.setdefault('names', names)
                if _pd is not None:   # pandas is faster
                    kwargs.setdefault('delimiter', '\s+')
                    kwargs.setdefault('comment', '#')
                    self.data = _pd.read_csv(fname, *args, **kwargs).to_records()
                else:
                    kwargs.setdefault('delimiter', None)
                    kwargs.setdefault('comments', '#')
                    kwargs.setdefault('skip_header', n)
                    self.data = np.recfromtxt(fname, *args, **kwargs)
                self.header = header
                self._units.update(**units)
                self._desc.update(**comments)
                self._aliases.update(**aliases)
            elif (extension == 'fits') or dtype == 'fits':
                if pyfits is None:
                    raise RuntimeError('Cannot read this format, Astropy or pyfits not found')
                if ('extname' not in kwargs) and ('ext' not in kwargs) and (len(args) == 0):
                    args = (1, )
                self.data = np.array(pyfits.getdata(fname, *args, **kwargs))
                header, aliases, units, comments = _fits_read_header(pyfits.getheader(fname, *args, **kwargs))
                self.header = header
                self._desc.update(**comments)
                self._units.update(**units)
                self._aliases.update(**aliases)
            elif (extension in ('hdf5', 'hd5', 'hdf')) or (dtype in ('hdf5', 'hd5', 'hdf')):
                if tables is None:
                    raise RuntimeError('Cannot read this format, pytables not found')
                hdr, aliases, units, desc, data = _hdf5_read_data(fname, *args, **kwargs)
                self.data = data
                self.header = hdr
                self._units.update(**units)
                self._desc.update(**desc)
                self._aliases.update(**aliases)
            elif (extension in ('vot', 'votable')) or (dtype in ('vot', 'votable')):
                # Votable case
                if _astropytable is None:
                    raise RuntimeError('Cannot read this votable format, astropy not found')
                data = _astropytable.read(fname, format='votable', *args, **kwargs)
                units = [(k, data[k].unit.name) for k in data.keys()]
                desc = [(k, data[k].description) for k in data.keys()]
                self.data = data.as_array()
                self.header = {}
                self._units.update(units)
                self._desc.update(desc)
            else:
                raise Exception('Format {0:s} not handled'.format(extension))
        elif type(fname) == np.ndarray:
            self.data = fname
            self.header = {}
        elif type(fname) == pyfits.FITS_rec:
            self.data = np.array(fname)
            self.header = {}
        elif isinstance(fname, SimpleTable):
            cp = kwargs.pop('copy', True)
            if cp:
                self.data = deepcopy(fname.data)
                self.header = deepcopy(fname.header)
                self._aliases = deepcopy(fname._aliases)
                self._units = deepcopy(fname._units)
                self._desc = deepcopy(fname._desc)
            else:
                self.data = fname.data
                self.header = fname.header
                self._aliases = fname._aliases
                self._units = fname._units
                self._desc = fname._desc
        elif hasattr(fname, 'dtype'):
            self.data = np.array(fname)
            self.header = {}
        else:
            raise Exception('Type {0!s:s} not handled'.format(type(fname)))
        if 'NAME' not in self.header:
            if type(fname) not in basestring:
                self.header['NAME'] = 'No Name'
            else:
                self.header['NAME'] = fname

    def pprint_entry(self, num, keys=None):
        """ print one line with key and values properly to be readable

        Parameters
        ----------
        num: int, slice
            indice selection

        keys: sequence or str
            if str, can be a regular expression
            if sequence, the sequence of keys to print
        """
        if (keys is None) or (keys == '*'):
            _keys = self.keys()
        elif type(keys) in basestring:
            _keys = [k for k in (self.keys() + tuple(self._aliases.keys()))
                     if (re.match(keys, k) is not None)]
        else:
            _keys = keys

        length = max(map(len, _keys))
        fmt = '{{0:{0:d}s}}: {{1}}'.format(length)
        data = self[num]

        for k in _keys:
            print(fmt.format(k, data[self.resolve_alias(k)]))

    def pprint(self, idx=None, fields=None, ret=False, all=False,
               full_match=False, headerChar='-', delim=' | ', endline='\n',
               **kwargs):
        """ Pretty print the table content
            you can select the table parts to display using idx to
            select the rows and fields to only display some columns
            (ret is only for insternal use)

        Parameters
        ----------

        idx: sequence, slide
            sub selection to print

        fields: str, sequence
            if str can be a regular expression, and/or list of fields separated
            by spaces or commas

        ret: bool
            if set return the string representation instead of printing the result

        all: bool
            if set, force to show all rows

        headerChar: char
            Character to be used for the row separator line

        delim: char
            The column delimiter.
        """
        if full_match is True:
            fn = re.fullmatch
        else:
            fn = re.match

        if (fields is None) or (fields == '*'):
            _keys = self.keys()
        elif type(fields) in basestring:
            if ',' in fields:
                _fields = fields.split(',')
            elif ' ' in fields:
                _fields = fields.split()
            else:
                _fields = [fields]
            lbls = self.keys() + tuple(self._aliases.keys())
            _keys = []
            for _fk in _fields:
                _keys += [k for k in lbls if (fn(_fk, k) is not None)]
        else:
            lbls = self.keys() + tuple(self._aliases.keys())
            _keys = []
            for _fk in _fields:
                _keys += [k for k in lbls if (fn(_fk, k) is not None)]

        nfields = len(_keys)

        fields = list(map( self.resolve_alias, _keys ))

        if idx is None:
            if (self.nrows < 10) or all:
                rows = [ [ str(self[k][rk]) for k in _keys ] for rk in range(self.nrows)]
            else:
                _idx = range(6)
                rows = [ [ str(self[k][rk]) for k in _keys ] for rk in range(5) ]
                if nfields > 1:
                    rows += [ ['...' for k in range(nfields) ] ]
                else:
                    rows += [ ['...' for k in range(nfields) ] ]
                rows += [ [ str(self[k][rk]) for k in fields ] for rk in range(-5, 0)]
        elif isinstance(idx, slice):
            _idx = range(idx.start, idx.stop, idx.step or 1)
            rows = [ [ str(self[k][rk]) for k in fields ] for rk in _idx]
        else:
            rows = [ [ str(self[k][rk]) for k in fields ] for rk in idx]

        if len(self._units) == 0:
            units = None
        else:
            units = [ '(' + str( self._units.get(k, None) or '') + ')' for k in fields ]

        out = __indent__(rows, header=_keys, units=units, delim=delim,
                         headerChar=headerChar, endline=endline)
        if ret is True:
            return out
        else:
            print(out)

    def write(self, fname, **kwargs):
        """ write table into file

        Parameters
        ----------
        fname: str
            filename to export the table into

        .. note::
            additional keywords are forwarded to the corresponding libraries
            :func:`pyfits.writeto` or :func:`pyfits.append`
            :func:`np.savetxt`
        """
        extension = kwargs.pop('extension', None)
        if extension is None:
            extension = fname.split('.')[-1]
        if (extension == 'csv'):
            comments = kwargs.pop('comments', '#')
            delimiter = kwargs.pop('delimiter', ',')
            commentedHeader = kwargs.pop('commentedHeader', False)
            hdr = _ascii_generate_header(self, comments=comments, delimiter=delimiter,
                                         commentedHeader=commentedHeader)
            header = kwargs.pop('header', hdr)
            np.savetxt(fname, self.data, delimiter=delimiter, header=header,
                       comments='', **kwargs)
        elif (extension in ['txt', 'dat']):
            comments = kwargs.pop('comments', '#')
            delimiter = kwargs.pop('delimiter', ' ')
            commentedHeader = kwargs.pop('commentedHeader', True)
            hdr = _ascii_generate_header(self, comments=comments, delimiter=delimiter,
                                         commentedHeader=commentedHeader)
            header = kwargs.pop('header', hdr)
            np.savetxt(fname, self.data, delimiter=delimiter, header=header,
                       comments='', **kwargs)
        elif (extension == 'fits'):
            hdr0 = kwargs.pop('header', None)
            append = kwargs.pop('append', False)
            hdr = _fits_generate_header(self)
            if hdr0 is not None:
                hdr.update(**hdr0)
            if append:
                _fits_append(fname, self.data, hdr, **kwargs)
            else:
                # patched version to correctly include the header
                _fits_writeto(fname, self.data, hdr, **kwargs)
        elif (extension in ('hdf', 'hdf5', 'hd5')):
            _hdf5_write_data(fname, self.data, header=self.header,
                             units=self._units, comments=self._desc,
                             aliases=self._aliases, **kwargs)
        else:
            raise Exception('Format {0:s} not handled'.format(extension))

    def to_records(self, **kwargs):
        """ Construct a numpy record array from this dataframe """
        return self.data

    def to_pandas(self, **kwargs):
        """ Construct a pandas dataframe

        Parameters
        ----------
        data : ndarray 
            (structured dtype), list of tuples, dict, or DataFrame
        keys: sequence, optional
            ordered subset of columns to export
        index : string, list of fields, array-like
            Field of array to use as the index, alternately a specific set of
            input labels to use
        exclude : sequence, default None
            Columns or fields to exclude
        columns : sequence, default None
            Column names to use. If the passed data do not have names
            associated with them, this argument provides names for the
            columns. Otherwise this argument indicates the order of the columns
            in the result (any names not found in the data will become all-NA
            columns)
        coerce_float : boolean, default False
            Attempt to convert values to non-string, non-numeric objects (like
            decimal.Decimal) to floating point, useful for SQL result sets

        Returns
        -------
        df : DataFrame
        """
        try:
            from pandas import DataFrame
            keys = kwargs.pop('keys', None)
            return DataFrame.from_dict(self.to_dict(keys=keys), **kwargs)
        except ImportError as error:
            print("Pandas import error")
            raise error

    def to_dict(self, keys=None, contiguous=False):
        """ Construct a dictionary from this dataframe with contiguous arrays

        Parameters
        ----------
        keys: sequence, optional
            ordered subset of columns to export

        contiguous: boolean
            make sure each value is a contiguous numpy array object
            (C-aligned)

        Returns
        -------
        data: dict
            converted data
        """
        if keys is None:
            keys = self.keys()
        if contiguous:
            return {k: np.ascontiguousarray(self[k]) for k in keys}
        return {k: self[k] for k in keys}

    def to_xarray(self, **kwargs):
        """ Construct an xarray dataset

        Each column will be converted into an independent variable in the
        Dataset. If the dataframe's index is a MultiIndex, it will be expanded
        into a tensor product of one-dimensional indices (filling in missing
        values with NaN). This method will produce a Dataset very similar to
        that on which the 'to_dataframe' method was called, except with
        possibly redundant dimensions (since all dataset variables will have
        the same dimensionality).
        """
        try:
            from xarray import Dataset
            return Dataset.from_dataframe(self.to_pandas(**kwargs))
        except ImportError as error:
            print("xray import error")
            raise error

    def to_vaex(self, **kwargs):
        """
        Create an in memory Vaex dataset

        Parameters
        ----------
        name: str
            unique for the dataset
        keys: sequence, optional
            ordered subset of columns to export

        Returns
        -------
        df: vaex.DataSetArrays
            vaex dataset
        """
        try:
            import vaex
            return vaex.from_arrays(**self.to_dict(contiguous=True, **kwargs))
        except ImportError as error:
            print("Vaex import error")
            raise error

    def to_dask(self, **kwargs):
        """ Construct a Dask DataFrame

        This splits an in-memory Pandas dataframe into several parts and constructs
        a dask.dataframe from those parts on which Dask.dataframe can operate in
        parallel.

        Note that, despite parallelism, Dask.dataframe may not always be faster
        than Pandas.  We recommend that you stay with Pandas for as long as
        possible before switching to Dask.dataframe.

        Parameters
        ----------
        keys: sequence, optional
            ordered subset of columns to export
        npartitions : int, optional
            The number of partitions of the index to create. Note that depending on
            the size and index of the dataframe, the output may have fewer
            partitions than requested.
        chunksize : int, optional
            The size of the partitions of the index.
        sort: bool
            Sort input first to obtain cleanly divided partitions or don't sort and
            don't get cleanly divided partitions
        name: string, optional
            An optional keyname for the dataframe.  Defaults to hashing the input

        Returns
        -------
        dask.DataFrame or dask.Series
            A dask DataFrame/Series partitioned along the index
        """
        try:
            from dask import dataframe
            keys = kwargs.pop('keys', None)
            return dataframe.from_pandas(self.to_pandas(keys=keys), **kwargs)
        except ImportError as error:
            print("Dask import error")
            raise error

    def to_astropy_table(self, **kwargs):
        """
        A class to represent tables of heterogeneous data.

        `astropy.table.Table` provides a class for heterogeneous tabular data,
        making use of a `numpy` structured array internally to store the data
        values.  A key enhancement provided by the `Table` class is the ability
        to easily modify the structure of the table by adding or removing
        columns, or adding new rows of data.  In addition table and column
        metadata are fully supported.

        Parameters
        ----------
        masked : bool, optional
            Specify whether the table is masked.
        names : list, optional
            Specify column names
        dtype : list, optional
            Specify column data types
        meta : dict, optional
            Metadata associated with the table.
        copy : bool, optional
            Copy the input data (default=True).
        rows : numpy ndarray, list of lists, optional
            Row-oriented data for table instead of ``data`` argument
        copy_indices : bool, optional
            Copy any indices in the input data (default=True)
        **kwargs : dict, optional
            Additional keyword args when converting table-like object

        Returns
        -------
        df: astropy.table.Table
            dataframe
        """
        try:
            from astropy.table import Table
            keys = kwargs.pop('keys', None)
            return Table(self.to_records(keys=keys), **kwargs)
        except ImportError as e:
            print("Astropy import error")
            raise e

    def _repr_html_(self):
        return self.to_pandas().head()._repr_html_()

    def set_alias(self, alias, colname):
        """
        Define an alias to a column

        Parameters
        ----------
        alias: str
            The new alias of the column

        colname: str
            The column being aliased
        """
        if (colname not in self.keys()):
            raise KeyError("Column {0:s} does not exist".format(colname))
        self._aliases[alias] = colname

    def reverse_alias(self, colname):
        """
        Return aliases of a given column.

        Given a colname, return a sequence of aliases associated to this column
        Aliases are defined by using .define_alias()
        """
        _colname = self.resolve_alias(colname)
        if (_colname not in self.keys()):
            raise KeyError("Column {0:s} does not exist".format(colname))

        return tuple([ k for (k, v) in self._aliases.iteritems() if (v == _colname) ])

    def resolve_alias(self, colname):
        """
        Return the name of an aliased column.

        Given an alias, return the column name it aliases. This
        function is a no-op if the alias is a column name itself.

        Aliases are defined by using .define_alias()
        """
        # User aliases
        if hasattr(colname, '__iter__') & (type(colname) not in basestring):
            return [ self.resolve_alias(k) for k in colname ]
        else:
            if self.caseless is True:
                maps = dict( [ (k.lower(), v) for k, v in self._aliases.items() ] )
                maps.update( (k.lower(), k) for k in self.keys() )
                return maps.get(colname.lower(), colname)
            else:
                return self._aliases.get(colname, colname)

    def set_unit(self, colname, unit):
        """ Set the unit of a column referenced by its name

        Parameters
        ----------
        colname: str
            column name or registered alias

        unit: str
            unit description
        """
        if isinstance(unit, basestring) and isinstance(colname, basestring):
            self._units[self.resolve_alias(colname)] = str(unit)
        else:
            for k, v in zip(colname, unit):
                self._units[self.resolve_alias(k)] = str(v)

    def set_comment(self, colname, comment):
        """ Set the comment of a column referenced by its name

        Parameters
        ----------
        colname: str
            column name or registered alias

        comment: str
            column description
        """
        if isinstance(comment, basestring) and isinstance(colname, basestring):
            self._desc[self.resolve_alias(colname)] = str(comment)
        else:
            for k, v in zip(colname, comment):
                self._desc[self.resolve_alias(k)] = str(v)

    def keys(self, regexp=None, full_match=False):
        """
        Return the data column names or a subset of it

        Parameters
        ----------
        regexp: str
            pattern to filter the keys with

        full_match: bool
            if set, use :func:`re.fullmatch` instead of :func:`re.match`

        Try to apply the pattern at the start of the string, returning
        a match object, or None if no match was found.

        returns
        -------
        seq: sequence
            sequence of keys
        """
        if (regexp is None) or (regexp == '*'):
            return self.colnames
        elif type(regexp) in basestring:
            if full_match is True:
                fn = re.fullmatch
            else:
                fn = re.match

            if regexp.count(',') > 0:
                _re = regexp.split(',')
            elif regexp.count(' ') > 0:
                _re = regexp.split()
            else:
                _re = [regexp]

            lbls = self.colnames + tuple(self._aliases.keys())
            _keys = []
            for _rk in _re:
                _keys += [k for k in lbls if (fn(_rk, k) is not None)]

            return _keys
        elif hasattr(regexp, '__iter__'):
            _keys = []
            for k in regexp:
                _keys += self.keys(k)
            return _keys
        else:
            raise ValueError('Unexpected type {0} for regexp'.format(type(regexp)))

    @property
    def name(self):
        """ name of the table given by the Header['NAME'] attribute """
        return self.header.get('NAME', None)

    @property
    def colnames(self):
        """ Sequence of column names """
        return self.data.dtype.names

    @property
    def ncols(self):
        """ number of columns """
        return len(self.colnames)

    @property
    def nrows(self):
        """ number of lines """
        return len(self.data)

    @property
    def nbytes(self):
        """ number of bytes of the object """
        n = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k)
                for k in self.__dict__.values())
        return n

    def __len__(self):
        """ number of lines """
        return self.nrows

    @property
    def shape(self):
        """ shape of the data """
        return self.data.shape

    @property
    def dtype(self):
        """ dtype of the data """
        return self.data.dtype

    @property
    def Plotter(self):
        """ Plotter instance related to this dataset.
        Requires plotter add-on to work """
        if Plotter is None:
            raise AttributeError('the add-on was not found, this property is not available')
        else:
            return Plotter(self, label=self.name)

    def __getitem__(self, v):
        return np.asarray(self.data.__getitem__(self.resolve_alias(v)))

    def take(self, indices, axis=None, out=None, mode='raise'):
        """
        Take elements from an array along an axis.

        This function does the same thing as "fancy" indexing (indexing arrays
        using arrays); however, it can be easier to use if you need elements
        along a given axis.

        Parameters
        ----------
        indices : array_like
            The indices of the values to extract.
            Also allow scalars for indices.

        axis : int, optional
            The axis over which to select values. By default, the flattened
            input array is used.

        out : ndarray, optional
            If provided, the result will be placed in this array. It should
            be of the appropriate shape and dtype.

        mode : {'raise', 'wrap', 'clip'}, optional
            Specifies how out-of-bounds indices will behave.

            * 'raise' -- raise an error (default)
            * 'wrap' -- wrap around
            * 'clip' -- clip to the range

            'clip' mode means that all indices that are too large are replaced
            by the index that addresses the last element along that axis. Note
            that this disables indexing with negative numbers.

        Returns
        -------
        subarray : ndarray
            The returned array has the same type as `a`.
        """
        return self.data.take(indices, axis, out, mode)

    def compress(self, condition, axis=None, out=None):
        """
        Return selected slices of an array along given axis.

        When working along a given axis, a slice along that axis is returned in
        `output` for each index where `condition` evaluates to True. When
        working on a 1-D array, `compress` is equivalent to `extract`.

        Parameters
        ----------
        condition : 1-D array of bools
            Array that selects which entries to return. If len(condition)
            is less than the size of `a` along the given axis, then output is
            truncated to the length of the condition array.

        axis : int, optional
            Axis along which to take slices. If None (default), work on the
            flattened array.

        out : ndarray, optional
            Output array.  Its type is preserved and it must be of the right
            shape to hold the output.

        Returns
        -------
        compressed_array : ndarray
            A copy of `a` without the slices along axis for which `condition`
            is false.
        """
        return self.data.compress(condition, axis, out)

    def get(self, v, full_match=False):
        """ returns a table from columns given as v

        this function is equivalent to :func:`__getitem__` but preserve the
        Table format and associated properties (units, description, header)

        Parameters
        ----------
        v: str
            pattern to filter the keys with

        full_match: bool
            if set, use :func:`re.fullmatch` instead of :func:`re.match`

        """
        new_keys = self.keys(v)
        t = self.__class__(self[new_keys])
        t.header.update(**self.header)
        t._aliases.update((k, v) for (k, v) in self._aliases.items() if v in new_keys)
        t._units.update((k, v) for (k, v) in self._units.items() if v in new_keys)
        t._desc.update((k, v) for (k, v) in self._desc.items() if v in new_keys)
        return t

    def __setitem__(self, k, v):
        if k in self:
            return self.data.__setitem__(self.resolve_alias(k), v)
        else:
            object.__setitem__(self, k, v)

    def __getattr__(self, k):
        try:
            return self.data.__getitem__(self.resolve_alias(k))
        except:
            return object.__getattribute__(self, k)

    def __iter__(self):
        newtab = self.select('*', [0])
        for d in self.data:
            newtab.data[0] = d
            yield newtab
        # return self.data.__iter__()

    def iterkeys(self):
        """ Iterator over the columns of the table """
        for k in self.colnames:
            yield k

    def itervalues(self):
        """ Iterator over the lines of the table """
        for l in self.data:
            yield l

    def items(self):
        """ Iterator on the (key, value) pairs """
        for k in self.colnames:
            yield k, self[k]

    def info(self):
        """ prints information on the table """
        s = "\nTable: {name:s}\n       nrows={s.nrows:d}, ncols={s.ncols:d}, mem={size:s}"
        s = s.format(name=self.header.get('NAME', 'Noname'), s=self,
                     size=pretty_size_print(self.nbytes))

        s += '\n\nHeader:\n'
        vals = list(self.header.items())
        length = max(map(len, self.header.keys()))
        fmt = '\t{{0:{0:d}s}} {{1}}\n'.format(length)
        for k, v in vals:
            s += fmt.format(k, v)

        vals = [(k, self._units.get(k, ''), self._desc.get(k, ''))
                for k in self.colnames]
        lengths = [(len(k), len(self._units.get(k, '')), len(self._desc.get(k, '')))
                   for k in self.colnames]
        lengths = list(map(max, (zip(*lengths))))

        s += '\nColumns:\n'

        fmt = '\t{{0:{0:d}s}} {{1:{1:d}s}} {{2:{2:d}s}}\n'.format(*(k + 1 for k in lengths))
        for k, u, c in vals:
            s += fmt.format(k, u, c)

        print(s)

        if len(self._aliases) > 0:
            print("\nTable contains alias(es):")
            for k, v in self._aliases.items():
                print('\t{0:s} --> {1:s}'.format(k, v))

    def __repr__(self):
        s = object.__repr__(self)
        s += "\nTable: {name:s}\n       nrows={s.nrows:d}, ncols={s.ncols:d}, mem={size:s}"
        return s.format(name=self.header.get('NAME', 'Noname'), s=self,
                        size=pretty_size_print(self.nbytes))

    def __getslice__(self, i, j):
        return self.data.__getslice__(i, j)

    def __contains__(self, k):
        if hasattr(k, 'decode'):
            _k = k.decode('utf8')
        else:
            _k = k
        return (_k in self.colnames) or (_k in self._aliases)

    def __array__(self):
        return self.data

    def __call__(self, *args, **kwargs):
        if (len(args) > 0) or (len(kwargs) > 0):
            return self.evalexpr(*args, **kwargs)
        else:
            return self.info()

    def sort(self, keys, copy=False):
        """
        Sort the table inplace according to one or more keys. This operates on
        the existing table (and does not return a new table).

        Parameters
        ----------

        keys: str or seq(str)
            The key(s) to order by

        copy: bool
            if set returns a sorted copy instead of working inplace
        """
        if not hasattr(keys, '__iter__'):
            keys = [keys]

        if copy is False:
            self.data.sort(order=keys)
        else:
            t = self.__class__(self, copy=True)
            t.sort(keys, copy=False)
            return t

    def match(self, r2, key):
        """ Returns the indices at which the tables match
        matching uses 2 columns that are compared in values

        Parameters
        ----------
        r2:  Table
            second table to use

        key: str
            fields used for comparison.

        Returns
        -------
        indexes: tuple
            tuple of both indices list where the two columns match.
        """
        return np.where( np.equal.outer( self[key], r2[key] ) )

    def stack(self, r, *args, **kwargs):
        """
        Superposes arrays fields by fields inplace

        t.stack(t1, t2, t3, default=None, inplace=True)

        Parameters
        ----------
        r: Table
        """
        if not hasattr(r, 'data'):
            raise AttributeError('r should be a Table object')
        defaults = kwargs.get('defaults', None)
        inplace = kwargs.get('inplace', False)

        data = [self.data, r.data] + [k.data for k in args]
        sdata = recfunctions.stack_arrays(data, defaults, usemask=False,
                                          asrecarray=True)

        if inplace:
            self.data = sdata
        else:
            t = self.__class__(self)
            t.data = sdata
            return t

    def join_by(self, r2, key, jointype='inner', r1postfix='1', r2postfix='2',
                defaults=None, asrecarray=False, asTable=True):
        """
        Join arrays `r1` and `r2` on key `key`.

        The key should be either a string or a sequence of string corresponding
        to the fields used to join the array.
        An exception is raised if the `key` field cannot be found in the two input
        arrays.
        Neither `r1` nor `r2` should have any duplicates along `key`: the presence
        of duplicates will make the output quite unreliable. Note that duplicates
        are not looked for by the algorithm.

        Parameters
        ----------
        key: str or seq(str)
            corresponding to the fields used for comparison.

        r2: Table
            Table to join with

        jointype: str in {'inner', 'outer', 'leftouter'}
            * 'inner'     : returns the elements common to both r1 and r2.
            * 'outer'     : returns the common elements as well as the elements of r1 not in r2 and the elements of not in r2.
            * 'leftouter' : returns the common elements and the elements of r1 not in r2.

        r1postfix: str
            String appended to the names of the fields of r1 that are present in r2

        r2postfix:  str
            String appended to the names of the fields of r2 that are present in r1

        defaults:   dict
            Dictionary mapping field names to the corresponding default values.

        Returns
        -------
        tab: Table
            joined table

        .. note::

            * The output is sorted along the key.

            * A temporary array is formed by dropping the fields not in the key
              for the two arrays and concatenating the result. This array is
              then sorted, and the common entries selected. The output is
              constructed by filling the fields with the selected entries.
              Matching is not preserved if there are some duplicates...
        """
        arr = recfunctions.join_by(key, self.data, r2.data, jointype=jointype,
                                   r1postfix=r1postfix, r2postfix=r2postfix,
                                   defaults=defaults, usemask=False,
                                   asrecarray=True)

        return SimpleTable(arr)

    @property
    def empty_row(self):
        """ Return an empty row array respecting the table format """
        return np.rec.recarray(shape=(1,), dtype=self.data.dtype)

    def add_column(self, name, data, dtype=None, unit=None, description=None):
        """
        Add one or multiple columns to the table

        Parameters
        ----------
        name: str or sequence(str)
           The name(s) of the column(s) to add

        data: ndarray, or sequence of ndarray
            The column data, or sequence of columns

        dtype: dtype
            numpy dtype for the data to add

        unit: str
            The unit of the values in the column

        description: str
            A description of the content of the column
        """

        _data = np.array(data, dtype=dtype)
        dtype = _data.dtype

        # unknown type is converted to text
        if dtype.type == np.object_:
            if len(data) == 0:
                longest = 0
            else:
                longest = len(max(data, key=len))
                _data = np.asarray(data, dtype='|%iS' % longest)

        dtype = _data.dtype

        if len(self.data.dtype) > 0:
            # existing data in the table
            if type(name) in basestring:
                # _name = name.encode('utf8')
                _name = str(name)
            else:
                # _name = [k.encode('utf8') for k in name]
                _name = [str(k) for k in name]

            self.data = recfunctions.append_fields(self.data, _name, _data,
                                                   dtypes=dtype, usemask=False,
                                                   asrecarray=True)

        else:
            if _data.ndim > 1:
                newdtype = (str(name), _data.dtype, (_data.shape[1],))
            else:
                newdtype = (str(name), _data.dtype)
            self.data = np.array(_data, dtype=[newdtype])

        if unit is not None:
            self.set_unit(name, unit)

        if description is not None:
            self.set_comment(name, description)

    def append_row(self, iterable):
        """
        Append one row in this table.

        see also: :func:`stack`

        Parameters
        ----------
        iterable: iterable
            line to add
        """
        if (len(iterable) != self.ncols):
            raise AttributeError('Expecting as many items as columns')
        r = self.empty_row
        for k, v in enumerate(iterable):
            r[0][k] = v
        self.stack(r)

    def remove_columns(self, names):
        """
        Remove several columns from the table

        Parameters
        ----------
        names: sequence
            A list containing the names of the columns to remove
        """
        self.pop_columns(names)

    def pop_columns(self, names):
        """
        Pop several columns from the table

        Parameters
        ----------

        names: sequence
            A list containing the names of the columns to remove

        Returns
        -------

        values: tuple
            list of columns
        """

        if not hasattr(names, '__iter__') or type(names) in basestring:
            names = [names]

        p = [self[k] for k in names]

        _names = set([ self.resolve_alias(k) for k in names ])
        self.data = recfunctions.drop_fields(self.data, _names)
        for k in names:
            self._aliases.pop(k, None)
            self._units.pop(k, None)
            self._desc.pop(k, None)

        return p

    def find_duplicate(self, index_only=False, values_only=False):
        """Find duplication in the table entries, return a list of duplicated
        elements Only works at this time is 2 lines are *the same entry* not if
        2 lines have *the same values*
        """
        dup = []
        idd = []
        for i in range(len(self.data)):
            if (self.data[i] in self.data[i + 1:]):
                if (self.data[i] not in dup):
                    dup.append(self.data[i])
                    idd.append(i)
        if index_only:
            return idd
        elif values_only:
            return dup
        else:
            return zip(idd, dup)

    def evalexpr(self, expr, exprvars=None, dtype=float):
        """ evaluate expression based on the data and external variables
            all np function can be used (log, exp, pi...)

        Parameters
        ----------
        expr: str
            expression to evaluate on the table
            includes mathematical operations and attribute names

        exprvars: dictionary, optional
            A dictionary that replaces the local operands in current frame.

        dtype: dtype definition
            dtype of the output array

        Returns
        -------
        out : NumPy array
            array of the result
        """
        _globals = {}
        for k in ( list(self.colnames) + list(self._aliases.keys()) ):
            if k in expr:
                _globals[k] = self[k]

        if exprvars is not None:
            if (not (hasattr(exprvars, 'keys') & hasattr(exprvars, '__getitem__' ))):
                raise AttributeError("Expecting a dictionary-like as condvars")
            for k, v in ( exprvars.items() ):
                _globals[k] = v

        # evaluate expression, to obtain the final filter
        r    = np.empty( self.nrows, dtype=dtype)
        r[:] = eval(expr, _globals, np.__dict__)

        return r

    def where(self, condition, condvars=None, *args, **kwargs):
        """ Read table data fulfilling the given `condition`.
        Only the rows fulfilling the `condition` are included in the result.

        Parameters
        ----------
        condition: str
            expression to evaluate on the table
            includes mathematical operations and attribute names

        condvars: dictionary, optional
            A dictionary that replaces the local operands in current frame.

        Returns
        -------
        out: ndarray/ tuple of ndarrays
        result equivalent to :func:`np.where`

        """
        ind = np.where(self.evalexpr(condition, condvars, dtype=bool ), *args, **kwargs)
        return ind

    def select(self, fields, indices=None, **kwargs):
        """
        Select only a few fields in the table

        Parameters
        ----------
        fields: str or sequence
            fields to keep in the resulting table

        indices: sequence or slice
            extract only on these indices

        returns
        -------
        tab: SimpleTable instance
            resulting table
        """
        _fields = self.keys(fields)

        if fields == '*':
            if indices is None:
                return self
            else:
                tab = self.__class__(self[indices])
                for k in self.__dict__.keys():
                    if k not in ('data', ):
                        setattr(tab, k, deepcopy(self.__dict__[k]))
                return tab
        else:
            d = {}
            for k in _fields:
                _k = self.resolve_alias(k)
                if indices is not None:
                    d[k] = self[_k][indices]
                else:
                    d[k] = self[_k]
            d['header'] = deepcopy(self.header)
            tab = self.__class__(d)
            for k in self.__dict__.keys():
                if k not in ('data', ):
                    setattr(tab, k, deepcopy(self.__dict__[k]))
            return tab

    def selectWhere(self, fields, condition, condvars=None, **kwargs):
        """ Read table data fulfilling the given `condition`.
            Only the rows fulfilling the `condition` are included in the result.

        Parameters
        ----------
        fields: str or sequence
            fields to keep in the resulting table

        condition: str
            expression to evaluate on the table
            includes mathematical operations and attribute names

        condvars: dictionary, optional
            A dictionary that replaces the local operands in current frame.

        Returns
        -------
        tab: SimpleTable instance
            resulting table
        """
        if condition in [True, 'True', None]:
            ind = None
        else:
            ind = self.where(condition, condvars, **kwargs)[0]

        tab = self.select(fields, indices=ind)

        return tab

    def groupby(self, *key):
        """
        Create an iterator which returns (key, sub-table) grouped by each value
        of key(value)

        Parameters
        ----------
        key: str
            expression or pattern to filter the keys with

        Returns
        -------
        key: str or sequence
            group key

        tab: SimpleTable instance
           sub-table of the group
           header, aliases and column metadata are preserved (linked to the
           master table).
        """
        _key = self.keys(key)
        getter = operator.itemgetter(*_key)

        for k, grp in itertools.groupby(self.data, getter):
            t = self.__class__(np.dstack(grp))
            t.header = self.header
            t._aliases = self._aliases
            t._units = self._units
            t._desc = self._desc
            yield (k, t)

    def stats(self, fn=None, fields=None, fill=None):
        """ Make statistics on columns of a table

        Parameters
        ----------
        fn: callable or sequence of callables
            functions to apply to each column
            default: (np.mean, np.std, np.nanmin, np.nanmax)

        fields: str or sequence
            any key or key expression to subselect columns
            default is all columns

        fill: value
            value when not applicable
            default np.nan

        returns
        -------
        tab: Table instance
            collection of statistics, one column per function in fn and 1 ligne
            per column in the table
        """
        from collections import OrderedDict

        if fn is None:
            fn = (stats.mean, stats.std,
                stats.min, stats.max,
                stats.has_nan)

        d = OrderedDict()
        d.setdefault('FIELD', [])
        for k in fn:
            d.setdefault(k.__name__, [])

        if fields is None:
            fields = self.colnames
        else:
            fields = self.keys(fields)

        if fill is None:
            fill = np.nan

        for k in fields:
            d['FIELD'].append(k)
            for fnk in fn:
                try:
                    val = fnk(self[k])
                except:
                    val = fill
                d[fnk.__name__].append(val)

        return self.__class__(d, dtype=dict)

    # method aliases
    remove_column = remove_columns

    # deprecated methods
    addCol = add_column
    addLine = append_row
    setComment = set_comment
    setUnit = set_unit
    delCol = remove_columns


class AstroTable(SimpleTable):
    """
    Derived from the Table, this class add implementations of common astro
    tools especially conesearch
    """
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self._ra_name, self._dec_name = self.__autoRADEC__()
        if (len(args) > 0):
            if isinstance(args[0], AstroTable):
                self._ra_name = args[0]._ra_name
                self._dec_name = args[0]._dec_name
        self._ra_name = kwargs.get('ra_name', self._ra_name)
        self._dec_name = kwargs.get('dec_name', self._dec_name)

    def __autoRADEC__(self):
        """ Tries to identify the columns containing RA and DEC coordinates """
        if 'ra' in self:
            ra_name = 'ra'
        elif 'RA' in self:
            ra_name = 'RA'
        else:
            ra_name = None
        if 'dec' in self:
            dec_name = 'dec'
        elif 'DEC' in self:
            dec_name = 'DEC'
        else:
            dec_name = None
        return ra_name, dec_name

    def set_RA(self, val):
        """ Set the column that defines RA coordinates """
        assert(val in self), 'column name {} not found in the table'.format(val)
        self._ra_name = val

    def set_DEC(self, val):
        """ Set the column that defines DEC coordinates """
        assert(val in self), 'column name {} not found in the table'.format(val)
        self._dec_name = val

    def get_RA(self, degree=True):
        """ Returns RA, converted from hexa/sexa into degrees """
        if self._ra_name is None:
            return None
        if (not degree) or (self.dtype[self._ra_name].kind != 'S'):
            return self[self._ra_name]
        else:
            if (len(str(self[0][self._ra_name]).split(':')) == 3):
                return np.asarray(AstroHelpers.hms2deg(self[self._ra_name],
                                                       delim=':'))
            elif (len(str(self[0][self._ra_name]).split(' ')) == 3):
                return np.asarray(AstroHelpers.hms2deg(self[self._ra_name],
                                                       delim=' '))
            else:
                raise Exception('RA Format not understood')

    def get_DEC(self, degree=True):
        """ Returns RA, converted from hexa/sexa into degrees """
        if self._dec_name is None:
            return None
        if (not degree) or (self.dtype[self._dec_name].kind != 'S'):
            return self[self._dec_name]
        else:
            if (len(str(self[0][self._dec_name]).split(':')) == 3):
                return np.asarray(AstroHelpers.dms2deg(self[self._dec_name],
                                                       delim=':'))
            elif (len(str(self[0][self._dec_name]).split(' ')) == 3):
                return np.asarray(AstroHelpers.dms2deg(self[self._dec_name],
                                                       delim=' '))
            else:
                raise Exception('RA Format not understood')

    def info(self):
        s = "\nTable: {name:s}\n       nrows={s.nrows:d}, ncols={s.ncols:d}, mem={size:s}"
        s = s.format(name=self.header.get('NAME', 'Noname'), s=self,
                     size=pretty_size_print(self.nbytes))

        s += '\n\nHeader:\n'
        vals = list(self.header.items())
        length = max(map(len, self.header.keys()))
        fmt = '\t{{0:{0:d}s}} {{1}}\n'.format(length)
        for k, v in vals:
            s += fmt.format(k, v)

        vals = [(k, self._units.get(k, ''), self._desc.get(k, ''))
                for k in self.colnames]
        lengths = [(len(k), len(self._units.get(k, '')), len(self._desc.get(k, '')))
                   for k in self.colnames]
        lengths = list(map(max, (zip(*lengths))))

        if (self._ra_name is not None) & (self._dec_name is not None):
            s += "\nPosition coordinate columns: {0}, {1}\n".format(self._ra_name,
                                                                    self._dec_name)

        s += '\nColumns:\n'

        fmt = '\t{{0:{0:d}s}} {{1:{1:d}s}} {{2:{2:d}s}}\n'.format(*(k + 1 for k in lengths))
        for k, u, c in vals:
            s += fmt.format(k, u, c)

        print(s)

        if len(self._aliases) > 0:
            print("\nTable contains alias(es):")
            for k, v in self._aliases.items():
                print('\t{0:s} --> {1:s}'.format(k, v))

    def coneSearch(self, ra, dec, r, outtype=0):
        """ Perform a cone search on a table

        Parameters
        ----------
        ra0: ndarray[ndim=1, dtype=float]
            column name to use as RA source in degrees

        dec0: ndarray[ndim=1, dtype=float]
            column name to use as DEC source in degrees

        ra: float
            ra to look for (in degree)

        dec: float
            ra to look for (in degree)

        r: float
            distance in degrees

        outtype: int
            type of outputs
                0 -- minimal, indices of matching coordinates
                1 -- indices and distances of matching coordinates
                2 -- full, boolean filter and distances

        Returns
        -------
        t: tuple
            if outtype is 0:
                only return indices from ra0, dec0
            elif outtype is 1:
                return indices from ra0, dec0 and distances
            elif outtype is 2:
                return conditional vector and distance to all ra0, dec0
        """
        if (self._ra_name is None) or (self._dec_name is None):
            raise AttributeError('Coordinate columns not set.')

        ra0  = self.get_RA()
        dec0 = self.get_DEC()
        return AstroHelpers.conesearch(ra0, dec0, ra, dec, r, outtype=outtype)

    def zoneSearch(self, ramin, ramax, decmin, decmax, outtype=0):
        """ Perform a zone search on a table, i.e., a rectangular selection

        Parameters
        ----------
        ramin: float
            minimal value of RA

        ramax: float
            maximal value of RA

        decmin: float
            minimal value of DEC

        decmax: float
            maximal value of DEC

        outtype: int
            type of outputs
                0 or 1 -- minimal, indices of matching coordinates
                2 -- full, boolean filter and distances

        Returns
        -------
        r: sequence
            indices or conditional sequence of matching values
        """

        assert( (self._ra_name is not None) & (self._dec_name is not None) ), 'Coordinate columns not set.'

        ra0  = self.get_RA()
        dec0 = self.get_DEC()
        ind = (ra0 >= ramin) & (ra0 <= ramax) & (dec0 >= decmin) & (dec0 <= decmax)
        if outtype <= 2:
            return ind
        else:
            return np.where(ind)

    def where(self, condition=None, condvars=None, cone=None, zone=None, **kwargs):
        """ Read table data fulfilling the given `condition`.
        Only the rows fulfilling the `condition` are included in the result.

        Parameters
        ----------
        condition: str
            expression to evaluate on the table
            includes mathematical operations and attribute names

        condvars: dictionary, optional
            A dictionary that replaces the local operands in current frame.

        Returns
        -------
        out: ndarray/ tuple of ndarrays
        result equivalent to :func:`np.where`
        """
        if cone is not None:
            if len(cone) != 3:
                raise ValueError('Expecting cone keywords as a triplet (ra, dec, r)')
        if zone is not None:
            if len(zone) != 4:
                raise ValueError('Expecting zone keywords as a tuple of 4 elements (ramin, ramax, decmin, decmax)')

        if condition is not None:
            ind = super(self.__class__, self).where(condition, **kwargs)
            if ind is None:
                if (cone is None) & (zone is None):
                    return None
        else:
            ind = True

        blobs = []
        if (cone is not None) and (zone is not None):  # cone + zone
            ra, dec, r = cone
            ind, d = self.coneSearch(ra, dec, r, outtype=2)
            ind = ind & self.zoneSearch(zone[0], zone[1], zone[2], zone[3], outtype=2)
            d = d[ind]
            blobs.append(d)
        elif (cone is not None):
            ra, dec, r = cone
            _ind, d = self.coneSearch(ra, dec, r, outtype=2)
            ind = ind & _ind.astype(bool)
            blobs.append(d[ind])
        elif (zone is not None):
            _ind = self.zoneSearch(zone[0], zone[1], zone[2], zone[3], outtype=1)
            ind = ind & _ind

        ind = np.where(ind)[0]

        return ind, blobs

    def selectWhere(self, fields, condition=None, condvars=None, cone=None, zone=None, **kwargs):
        """ Read table data fulfilling the given `condition`.
            Only the rows fulfilling the `condition` are included in the result.
            conesearch is also possible through the keyword cone formatted as (ra, dec, r)
            zonesearch is also possible through the keyword zone formatted as (ramin, ramax, decmin, decmax)

            Combination of multiple selections is also available.
        """
        ind, blobs = self.where(condition, condvars, cone, zone, **kwargs)
        tab = self.select(fields, indices=ind)

        if cone is not None:
            tab.add_column('separation', np.squeeze(blobs), unit='degree')

        if self._ra_name in tab:
            tab.set_RA(self._ra_name)

        if self._dec_name in tab:
            tab.set_DEC(self._dec_name)

        return tab


class stats(object):
    @classmethod
    def has_nan(s, v):
        return (True in np.isnan(v))

    @classmethod
    def mean(s, v):
        return np.nanmean(v)

    @classmethod
    def max(s, v):
        return np.nanmax(v)

    @classmethod
    def min(s, v):
        return np.nanmin(v)

    @classmethod
    def std(s, v):
        return np.nanstd(v)

    @classmethod
    def var(s, v):
        return np.var(v)

    @classmethod
    def p16(s, v):
        try:
            return np.nanpercentile(v, 16)
        except AttributeError:
            return np.percentile(v, 16)

    @classmethod
    def p84(s, v):
        try:
            return np.nanpercentile(v, 84)
        except AttributeError:
            return np.percentile(v, 84)

    @classmethod
    def p50(s, v):
        try:
            return np.nanmedian(v)
        except AttributeError:
            return np.percentile(v, 50)


'''
# =============================================================================
# Adding some plotting functions
# =============================================================================

try:
    import pylab as plt

    def plot_function(tab, fn, *args, **kwargs):
        """ Generate a plotting method of tab from a given function

        Parameters
        ----------
        tab: SimpleTable instance
            table instance

        fn: str or callable
            if str, will try a function in matplotlib
            if callable, calls the function directly

        xname: str
            expecting a column name from the table

        yname: str, optional
            if provided, another column to use for the plot

        onlywhere: sequence or str, optional
            if provided, selects only data with this condition
            the condition can be a ndarray slice or a string.
            When a string is given, the evaluation calls :func:`SimpleTable.where`

        ax: matplotlib.Axes instance
            if provided make sure it uses the axis to do the plots if a mpl
            function is used.

        Returns
        -------
        r: object
            anything returned by the called function
        """
        if not hasattr(fn, '__call__'):
            ax = kwargs.pop('ax', None)
            if ax is None:
                ax = plt.gca()
            _fn = getattr(ax, fn, None)
            if _fn is None:
                raise AttributeError('function neither callable or found in matplotlib')
        else:
            _fn = fn

        onlywhere = kwargs.pop('onlywhere', None)
        if type(onlywhere) in basestring:
            select = tab.where(onlywhere)
        else:
            select = onlywhere

        _args = ()
        for a in args:
            if (hasattr(a, '__iter__')):
                try:
                    b = tab[a]
                    if select is not None:
                        b = b.compress(select)
                    if (len(b.dtype) > 1):
                        b = list((b[k] for k in b.dtype.names))
                    _args += (b, )
                except Exception as e:
                    print(e)
                    _args += (a, )
            else:
                _args += (a, )

        return _fn(*_args, **kwargs)

    def attached_function(fn, doc=None, errorlevel=0):
        """ eclare a function as a method to the class table"""

        def _fn(self, *args, **kwargs):
            try:
                return plot_function(self, fn, *args, **kwargs)
            except Exception as e:
                if errorlevel < 1:
                    pass
                else:
                    raise e

        if doc is not None:
            _fn.__doc__ = doc

        return _fn

    SimpleTable.plot_function = plot_function
    SimpleTable.plot = attached_function('plot', plt.plot.__doc__)
    SimpleTable.hist = attached_function('hist', plt.hist.__doc__)
    SimpleTable.hist2d = attached_function('hist2d', plt.hist2d.__doc__)
    SimpleTable.hexbin = attached_function('hexbin', plt.hexbin.__doc__)
    SimpleTable.scatter = attached_function('scatter', plt.scatter.__doc__)

    # newer version of matplotlib
    if hasattr(plt, 'violinplot'):
        SimpleTable.violinplot = attached_function('violinplot', plt.violinplot.__doc__)
    if hasattr(plt, 'boxplot'):
        SimpleTable.boxplot = attached_function('boxplot', plt.boxplot.__doc__)

except Exception as e:
    print(e)
'''
