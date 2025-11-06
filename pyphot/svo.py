"""Link to the SVO filter profile service

http://svo2.cab.inta-csic.es/theory/fps/

If your research benefits from the use of the SVO Filter Profile Service, include the following acknowledgement in your publication:

> This research has made use of the SVO Filter Profile Service
> (http://svo2.cab.inta-csic.es/theory/fps/) supported from the Spanish MINECO
> through grant AYA2017-84089.

and please include the following references in your publication:

* The SVO Filter Profile Service. Rodrigo, C., Solano, E., Bayo, A., 2012; https://ui.adsabs.harvard.edu/abs/2012ivoa.rept.1015R/abstract
* The SVO Filter Profile Service. Rodrigo, C., Solano, E., 2020; https://ui.adsabs.harvard.edu/abs/2020sea..confE.182R/abstract

Example
-------

>>> lst = "2MASS/2MASS.J 2MASS/2MASS.H 2MASS/2MASS.Ks HST/ACS_WFC.F475W HST/ACS_WFC.F814W".split()
    objects = [get_pyphot_filter(k) for k in lst]
"""

from io import BytesIO
from typing import List, Literal

import requests
from astropy.io import votable

from .phot import Filter

QUERY_URL: str = "http://svo2.cab.inta-csic.es/theory/fps/fps.php"
DETECTOR_TYPE: List[Literal["energy", "photon"]] = [
    "energy",
    "photon",
]  # svo returns 0, 1


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
    table = votable.parse_single_table(BytesIO(response.content))
    params = {p.name: p.value for p in table.params}
    tab = table.to_table()
    return Filter(
        tab["Wavelength"].to("nm"),
        tab["Transmission"],
        name=params["filterID"].replace("/", "_"),
        dtype=DETECTOR_TYPE[int(params["DetectorType"])],  # type: ignore
    )
