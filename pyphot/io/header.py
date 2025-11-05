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
