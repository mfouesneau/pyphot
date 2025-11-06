import pyphot

# define header and table format (as csv)
hdr = (
    "name",
    "wavelength units",
    "index units",
    "min",
    "maxmin blue",
    "max blue",
    "min red",
    "max red",
)
fmt = "{0:s},{1:s},{2:s},{3:.3f},{4:.3f},{5:.3f},{6:.5f},{7:.3f},{8:.3f}\n"

lib = pyphot.LickLibrary()

with open("licks_table.csv", "w") as output:
    output.write(",".join(hdr) + "\n")

    for k in sorted(lib.content):
        fk = lib[k]
        # wavelength have units
        band = fk.band.magnitude
        blue = fk.blue.magnitude
        red = fk.red.magnitude
        rec = (
            fk.name,
            fk.wavelength_unit,
            fk.index_unit,
            band[0],
            band[1],
            blue[0],
            blue[1],
            red[0],
            red[1],
        )
        output.write(fmt.format(*rec))
