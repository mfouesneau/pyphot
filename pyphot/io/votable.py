"""
VOTable parser for astronomical tabular data.

VOTable is the standard XML format for astronomical tabular data.
This module implements a custom `VOTableParser` that uses XML parsing and not other dependencies.
`from_votable` provides the standard interface of io operations (pandas.DataFrame, HeaderInfo)

"""

from io import IOBase
from os import PathLike
from typing import Dict, List, Any, Union, Tuple
import numpy as np
import numpy.typing as npt
import pandas as pd
import xml.etree.ElementTree as ET
import requests

from .header import HeaderInfo


class VOTableParser:
    """
    A custom VOTable parser using XML parsing.

    VOTable is the standard XML format for astronomical tabular data.
    This example shows how to parse the structure and extract data.
    """

    def __init__(
        self,
        source: Union[str, bytes, PathLike, IOBase],
        is_url: bool = False,
    ):
        """
        Initialize VOTable parser

        Parameters
        ----------
        source : str or bytes
            Either a file path, URL, or XML string/bytes
        is_url : bool
            If True, treat source as URL to fetch
        """
        if is_url and isinstance(source, (str, bytes)):
            response = requests.get(source)
            response.raise_for_status()
            xml_content = response.content
        elif isinstance(source, str) and source.startswith("<?xml"):
            xml_content = source
        elif isinstance(source, bytes) and source.startswith(b"<?xml"):
            xml_content = source.decode("utf-8")
        elif isinstance(source, IOBase):
            xml_content = source.read()
        else:
            # Assume it's a file path
            with open(source, encoding="utf-8") as f:
                xml_content = f.read()

        self.root = ET.fromstring(xml_content)
        self.namespace = self._get_namespace()
        self.tables = []
        self._parse_tables()

    def _get_namespace(self) -> Dict[str, str]:
        """Extract namespace from VOTable root element"""
        # VOTable typically uses namespace like http://www.ivoa.net/xml/VOTable/v1.3
        namespace = {}
        if self.root.tag.startswith("{"):
            uri = self.root.tag.split("}")[0][1:]
            namespace["vot"] = uri
        return namespace

    def _find_element(self, parent, tag: str):
        """Find element accounting for namespace"""
        if self.namespace:
            return parent.find(f"{{{self.namespace['vot']}}}{tag}")
        return parent.find(tag)

    def _findall_elements(self, parent, tag: str):
        """Find all elements accounting for namespace"""
        if self.namespace:
            return parent.findall(f".//{{{self.namespace['vot']}}}{tag}")
        return parent.findall(f".//{tag}")

    def _parse_tables(self):
        """Parse all TABLE elements in the VOTable"""
        table_elements = self._findall_elements(self.root, "TABLE")

        for table_elem in table_elements:
            table_data = self._parse_single_table(table_elem)
            self.tables.append(table_data)

    def _parse_single_table(self, table_elem) -> Dict[str, Any]:
        """Parse a single TABLE element"""
        table_info = {
            "name": table_elem.get("name", "unnamed"),
            "description": "",
            "fields": [],
            "params": [],
            "data": [],
        }

        # Get table description
        desc_elem = self._find_element(table_elem, "DESCRIPTION")
        if desc_elem is not None:
            table_info["description"] = desc_elem.text or ""

        # Parse FIELD elements (column definitions)
        field_elements = self._findall_elements(table_elem, "FIELD")
        for field_elem in field_elements:
            field_info = self._parse_field(field_elem)
            table_info["fields"].append(field_info)

        # Parse PARAM elements (parameters/metadata)
        param_elements = self._findall_elements(table_elem, "PARAM")
        for param_elem in param_elements:
            param_info = self._parse_param(param_elem)
            table_info["params"].append(param_info)

        # Parse DATA element
        data_elem = self._find_element(table_elem, "DATA")
        if data_elem is not None:
            table_info["data"] = self._parse_data(data_elem, table_info["fields"])

        return table_info

    def _parse_field(self, field_elem) -> Dict[str, Any]:
        """Parse FIELD element (column definition)"""
        field_info = {
            "name": field_elem.get("name", ""),
            "id": field_elem.get("ID", ""),
            "ucd": field_elem.get("ucd", ""),
            "datatype": field_elem.get("datatype", ""),
            "unit": field_elem.get("unit", ""),
            "arraysize": field_elem.get("arraysize", ""),
            "description": "",
        }

        # Get field description
        desc_elem = self._find_element(field_elem, "DESCRIPTION")
        if desc_elem is not None:
            field_info["description"] = desc_elem.text or ""

        return field_info

    def _parse_param(self, param_elem) -> Dict[str, Any]:
        """Parse PARAM element (parameter/metadata)"""
        param_info = {
            "name": param_elem.get("name", ""),
            "id": param_elem.get("ID", ""),
            "value": param_elem.get("value", ""),
            "datatype": param_elem.get("datatype", ""),
            "unit": param_elem.get("unit", ""),
            "ucd": param_elem.get("ucd", ""),
            "description": "",
        }

        desc_elem = self._find_element(param_elem, "DESCRIPTION")
        if desc_elem is not None:
            param_info["description"] = desc_elem.text or ""

        return param_info

    def _parse_data(self, data_elem, fields: List[Dict]) -> List[List[Any]]:
        """Parse DATA element - supports TABLEDATA format"""
        tabledata_elem = self._find_element(data_elem, "TABLEDATA")
        if tabledata_elem is None:
            return []

        rows = []
        tr_elements = self._findall_elements(tabledata_elem, "TR")

        for tr_elem in tr_elements:
            td_elements = self._findall_elements(tr_elem, "TD")
            row_data = []

            for i, td_elem in enumerate(td_elements):
                value = td_elem.text or ""

                # Convert based on field datatype
                if i < len(fields):
                    converted_value = self._convert_value(value, fields[i]["datatype"])
                    row_data.append(converted_value)
                else:
                    row_data.append(value)

            rows.append(row_data)

        return rows

    def _convert_value(self, value: str, datatype: str) -> Any:
        """Convert string value based on VOTable datatype"""
        if not value or value.strip() == "":
            return None

        try:
            if datatype in ["float", "double"]:
                return float(value)
            elif datatype in ["int", "short", "long"]:
                return int(value)
            elif datatype == "boolean":
                return value.lower() in ["true", "1", "yes"]
            else:
                return value
        except ValueError:
            return value

    def get_table_as_dict(self, table_index: int = 0) -> Dict[str, List]:
        """
        Convert table data to dictionary with column names as keys

        Parameters
        ----------
        table_index : int
            Index of table to convert (default: 0 for first table)

        Returns
        -------
        dict
            Dictionary with column names as keys and data as lists
        """
        if table_index >= len(self.tables):
            raise IndexError(f"Table index {table_index} out of range")

        table = self.tables[table_index]
        result = {}

        for i, field in enumerate(table["fields"]):
            column_name = field["name"] or field["id"] or f"col_{i}"
            column_data = []

            for row in table["data"]:
                if i < len(row):
                    column_data.append(row[i])
                else:
                    column_data.append(None)

            result[column_name] = column_data

        return result

    def print_table_info(self, table_index: int = 0):
        """Print information about a table"""
        if table_index >= len(self.tables):
            print(f"Table index {table_index} out of range")
            return

        table = self.tables[table_index]
        print(f"Table: {table['name']}")
        print(f"Description: {table['description']}")
        print(f"Number of fields: {len(table['fields'])}")
        print(f"Number of rows: {len(table['data'])}")

        print("\nFields:")
        for field in table["fields"]:
            print(f"  {field['name']} ({field['datatype']}) - {field['description']}")

        print("\nParameters:")
        for param in table["params"]:
            print(f"  {param['name']}: {param['value']} ({param['unit']})")


def _converter(val: Any, subtype: str) -> npt.NDArray:
    """Convert string arrays to appropriate subtype arrays"""
    data = np.array(val, dtype=object)
    # Replace all the None with an appropriate fill value
    mask = data == None  # noqa: E711
    kind = np.dtype(subtype).kind
    data[mask] = {"U": "", "S": b""}.get(kind, 0)
    return np.ma.array(data.astype(subtype), mask=mask)


def from_votable(
    fname: Union[str, bytes, IOBase, PathLike],
    *,
    table_index: int = 0,
    is_url: bool = False,
) -> Tuple[pd.DataFrame, HeaderInfo]:
    """Read a VOTable file and return a pandas DataFrame and header information.

    Parameters
    ----------
    fname : str, bytes, IOBase, PathLike
        The filename or file-like object to read.
    table_index : int, optional
        The index of the table to read, by default 0.
    is_url : bool, optional
        Whether the file is a URL, by default False.

    Returns
    -------
    Tuple[pd.DataFrame, HeaderInfo]
        A tuple containing the pandas DataFrame and header information.
    """
    parser = VOTableParser(fname, is_url=is_url)
    table = parser.tables[table_index]
    fields = table["fields"]
    columns = (field["name"] for field in fields)
    data_ = parser.get_table_as_dict(table_index)

    # Check subtypes if any and generate the converter function.
    converters = {}
    for field in fields:
        name = field["name"]
        dtype = field["datatype"]
        if dtype == "boolean":
            dtype = "bool"
        converters[name] = lambda x: _converter(x, dtype)

    data = {name: converters[name](data_[name]) for name in columns}
    df = pd.DataFrame.from_dict(data)
    df.name = table["name"]

    hdr = {}
    units = {}
    comments = {}

    for par in table["params"]:
        hdr[par["name"]] = {key: value for key, value in par.items() if key != "name"}
    for field in fields:
        hdr[field["name"]] = {
            key: value for key, value in field.items() if key != "name"
        }

    for field in fields:
        if field["unit"] not in (None, ""):
            units[field["name"]] = field["unit"]
        if field["description"] not in (None, ""):
            comments[field["name"]] = field["description"]
    header = HeaderInfo(
        header=hdr,
        alias={},
        units=units,
        comments=comments,
    )
    return df, header


def to_votable(*args, **kwargs):
    raise NotImplementedError("to_votable is not implemented yet.")
