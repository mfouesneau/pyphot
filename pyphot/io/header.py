"""Defines the HeaderInfo class that contains the metadata of a file."""

from typing import Any, Dict, Hashable
from dataclasses import dataclass


@dataclass
class HeaderInfo:
    """Extracted information from FITS header"""

    header: Dict[Hashable, Any]
    """Header dictionary containing any metadata from a file input"""
    alias: Dict[Hashable, str]
    """Alias dictionary which contains potential mappings of data columns to aliases"""
    units: Dict[Hashable, str]
    """Units dictionary containing potential mappings of data columns to units"""
    comments: Dict[Hashable, str]
    """Comments/description dictionary containing potential mappings of data columns to comments"""
