from typing import cast

import pyphot
from pyphot import Filter

# define header and table format (as csv)
hdr = (
    "name",
    "detector type",
    "wavelength units",
    "central wavelength",
    "pivot wavelength",
    "effective wavelength",
    "Vega mag",
    "Vega flux",
    "Vega Jy",
    "AB mag",
    "AB flux",
    "AB Jy",
    "ST mag",
    "ST flux",
    "ST Jy",
)
fmt = "{0:s},{1:s},{2:s},{3:.3f},{4:.3f},{5:.3f},{6:.5f},{7:.5g},{8:.5g},{9:.5f},{10:.5g},{11:.5g},{12:.5f},{13:.5g},{14:.5g}\n"

lib = pyphot.get_library()

with open("table.csv", "w") as output:
    output.write(",".join(hdr) + "\n")

    for k in sorted(lib.content):
        fk = cast(Filter, lib[k])
        rec = (
            fk.name,
            fk.dtype,
            fk.wavelength_unit,
            fk.cl.value,
            fk.lpivot.value,
            fk.leff.value,
            fk.Vega_zero_mag,
            fk.Vega_zero_flux.value,
            fk.Vega_zero_Jy.value,
            fk.AB_zero_mag,
            fk.AB_zero_flux.value,
            fk.AB_zero_Jy.value,
            fk.ST_zero_mag,
            fk.ST_zero_flux.value,
            fk.ST_zero_Jy.value,
        )
        output.write(fmt.format(*rec))
