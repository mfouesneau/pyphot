"""Link to the SVO filter profile service

See also their website http://svo2.cab.inta-csic.es/theory/fps/

.. important::
    If your research benefits from the use of the SVO Filter Profile Service, include the following acknowledgement in your publication:

    This research has made use of the SVO Filter Profile Service(http://svo2.cab.inta-csic.es/theory/fps/) supported from the Spanish MINECO through grant AYA2017-84089 and described in `Rodrigo et al, (2012) <https://ui.adsabs.harvard.edu/abs/2012ivoa.rept.1015R/abstract>`_ and `Rodrigo et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020sea..confE.182R/abstract>`_

    with references:

    * The SVO Filter Profile Service. Rodrigo, C., Solano, E., Bayo, A., 2012; `2012ivoa.rept.1015R <https://ui.adsabs.harvard.edu/abs/2012ivoa.rept.1015R/abstract>`_
    * The SVO Filter Profile Service. Rodrigo, C., Solano, E., 2020; `2020sea..confE.182R <https://ui.adsabs.harvard.edu/abs/2020sea..confE.182R/abstract>`_

Example
-------

>>> lst = "2MASS/2MASS.J 2MASS/2MASS.H 2MASS/2MASS.Ks HST/ACS_WFC.F475W HST/ACS_WFC.F814W".split()
    objects = [get_pyphot_filter(k) for k in lst]

.. note::
    This module uses :mod:`pyphot.io.votable` to parse the SVO filter profile service response. (i.e. it does not depend on Astropy)
"""

from io import BytesIO
from typing import List, Literal, cast

import requests

from . import config
from .phot import Filter
from .unit_adapters import QuantityType
from .io.votable import from_votable

QUERY_URL: str = "http://svo2.cab.inta-csic.es/theory/fps/fps.php"
DETECTOR_TYPE: List[Literal["energy", "photon"]] = [
    "energy",
    "photon",
]  # svo returns 0, 1

__all__ = ["get_pyphot_filter"]


def _get_pyphot_filter_astropy(identifier: str) -> Filter:
    """Query the SVO filter profile service and return the filter object
    This function uses the astropy.io.votable module to parse the response from the SVO filter profile service.
    This does not play well with other unit backends than astropy.units

    Parameters
    ----------
    identifier : str
        SVO identifier of the filter profile
        e.g., 2MASS/2MASS.Ks HST/ACS_WFC.F475W
        The identifier is the first column on the webpage of the facilities.

    Returns
    -------
    filter : Filter
        Filter object
    """
    from astropy.io import votable

    query = {"ID": identifier}
    response = requests.get(QUERY_URL, params=query)
    response.raise_for_status()
    table = votable.parse_single_table(BytesIO(response.content))
    params = {p.name: p.value for p in table.params}
    tab = table.to_table()
    return Filter(
        tab["Wavelength"].to("nm"),
        tab["Transmission"],
        name=params["filterID"].replace("/", "_"),
        dtype=DETECTOR_TYPE[int(params["DetectorType"])],  # type: ignore
    )


def get_pyphot_filter(identifier: str) -> Filter:
    """Query the SVO filter profile service and return the filter object

    Parameters
    ----------
    identifier : str
        SVO identifier of the filter profile
        e.g., 2MASS/2MASS.Ks HST/ACS_WFC.F475W
        The identifier is the first column on the webpage of the facilities.

    Returns
    -------
    filter : Filter
        Filter object
    """
    query = {"ID": identifier}
    response = requests.get(QUERY_URL, params=query)
    response.raise_for_status()
    df, header = from_votable(BytesIO(response.content))

    wavelength = cast(
        QuantityType,
        df["Wavelength"].to_numpy() * config.units.U(header.units["Wavelength"]),
    )
    transmission = df["Transmission"].to_numpy()
    name = header.header["filterID"]["value"].replace("/", "_")
    dtype = DETECTOR_TYPE[int(header.header["DetectorType"]["value"])]

    return Filter(wavelength.to("nm"), transmission, name=name, dtype=dtype)
