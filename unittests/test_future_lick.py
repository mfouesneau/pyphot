"""Testing the LickIndex and LickLibrary features"""

import pytest

from typing import cast, Union

from pyphot.future import licks
from pyphot.future import config

references = {
    "Mg_1": dict(
        blue=(4895.125, 4957.625),
        band=(5069.125, 5134.125),
        red=(5301.125, 5366.125),
        index_unit="mag",
        wavelength_unit="AA",
    ),
    "Hbeta0": dict(
        wavelength_unit="AA",
        band=[4839.275, 4877.097],
        blue=[4821.175, 4838.404],
        red=[4897.445, 4915.845],
        index_unit="ew",
    ),
}

lib = licks.LickLibrary()


def check_equality_ref_values(
    lick: Union[licks.LickDefinition, licks.LickIndex], info: dict
):
    for what in ("blue", "band", "red"):
        value = getattr(lick, what)
        assert cast(licks.QuantityType, value).unit == config.units.U(
            info["wavelength_unit"]
        )
        assert abs(value[0].value - info[what][0]) < 1e-6
        assert abs(value[1].value - info[what][1]) < 1e-6


@pytest.mark.parametrize("info", references.values())
def test_lickdefinition(info):
    """Make sure the definition is robust"""
    # check all ranges are set with wavelength_unit units
    # without altering the values
    lref = licks.LickDefinition(**info)  # type: ignore
    check_equality_ref_values(lref, info)

    # check consistencywhen attributes are with units
    l1 = licks.LickDefinition(
        blue=lref.blue.to("nm"),
        band=lref.band.to("AA"),  # band sets the consistency unit
        red=lref.red.to("cm"),
        index_unit=lref.index_unit,
    )
    for what in ("blue", "band", "red"):
        value = getattr(lref, what)
        assert cast(licks.QuantityType, value).unit == lref.band.unit
        assert abs(value[0].value - info[what][0]) < 1e-6
        assert abs(value[1].value - info[what][1]) < 1e-6
    assert l1.wavelength_unit == lref.band.unit


def test_library_against_references():
    lib = licks.LickLibrary()
    assert len(lib) > 0

    cns = lib.find("CN")
    assert cns == ["CN_1", "CN_2"]

    for name, info in references.items():
        index = cast(licks.LickIndex, lib[name])
        check_equality_ref_values(index, info)


@pytest.mark.parametrize("name", lib.content)
def test_all_indices(name: str):
    # convert to magnitudes
    from pyphot.future.vega import Vega

    vega = Vega()
    # using the internal collection of indices
    fk = cast(licks.LickIndex, lib["CN_1"])
    # work on many spectra at once
    index = fk.get(vega.wavelength, vega.flux, axis=-1)

    # define header and table format (as csv)
    """
    hdr = (
        "name",
        "wavelength units",
        "index units",
        "min",
        "maxmin blue",
        "max blue",
        "min red",
        "max red",
        "value(vega)",
    )
    """
    fmt = "{0:s},{1:s},{2:s},{3:.3f},{4:.3f},{5:.3f},{6:.5f},{7:.3f},{8:.3f},{9:.3f}\n"
    rec = (
        fk.name,
        fk.wavelength_unit,
        fk.index_unit,
        fk.band.value[0],
        fk.band.value[1],
        fk.blue.value[0],
        fk.blue.value[1],
        fk.red.value[0],
        fk.red.value[1],
        index,
    )
    print(fmt.format(*rec))
