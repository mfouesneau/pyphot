import unittest
import sys
sys.path.append('../')

class TestQuick(unittest.TestCase):
    def test_dependencies(self):
        from pyphot import sandbox as pyphot

    def test_library_irac(self):
        from pyphot import sandbox as pyphot
        # get the internal default library of passbands filters
        lib = pyphot.get_library()
        print("Library contains: ", len(lib), " filters")
        # find all filter names that relates to IRAC
        # and print some info
        f = lib.find('irac')
        for name in f:
            lib[name].info(show_zeropoints=True)

    def test_vega(self):
        from pyphot.vega import Vega
        vega = Vega()

    def test_convert_mags(self):
        import numpy as np
        from pyphot.vega import Vega
        from pyphot import sandbox as pyphot
        # get the internal default library of passbands filters
        lib = pyphot.get_library()
        vega = Vega()
        f = lib['HST_WFC3_F110W']
        # compute the integrated flux through the filter f
        # note that it work on many spectra at once
        fluxes = f.get_flux(vega.wavelength, vega.flux, axis=-1)

        # Note that fluxes is now with units of erg/s/cm2/AA
        # pyphot gives Vega in flam and can convert between flux density units.
        print(fluxes, vega.wavelength, vega.flux)

        # convert to vega magnitudes
        mags = -2.5 * np.log10(fluxes.magnitude) - f.Vega_zero_mag
        print("Vega magnitude of Vega in {0:s} is : {1:f} mag".format(f.name, mags))
        mags = -2.5 * np.log10(fluxes.magnitude) - f.AB_zero_mag
        print("AB magnitude of Vega in {0:s} is : {1:f} mag".format(f.name, mags))
        mags = -2.5 * np.log10(fluxes.magnitude) - f.ST_zero_mag
        print("ST magnitude of Vega in {0:s} is : {1:f} mag".format(f.name, mags))

    def test_load_all_filters(self):
        from pyphot import sandbox as pyphot
        # get the internal default library of passbands filters
        lib = pyphot.get_library()
        lib.load_all_filters()

    def test_all_filters_in_library(self):
        from pyphot import sandbox as pyphot

        # define header and table format (as csv)
        hdr = ("name", "detector type", "wavelength units",
            "central wavelength", "pivot wavelength", "effective wavelength",
            "Vega mag", "Vega flux", "Vega Jy",
            "AB mag", "AB flux", "AB Jy",
            "ST mag", "ST flux", "ST Jy")
        fmt = "{0:s},{1:s},{2:s},{3:.3f},{4:.3f},{5:.3f},{6:.5f},{7:.5g},{8:.5g},{9:.5f},{10:.5g},{11:.5g},{12:.5f},{13:.5g},{14:.5g}\n"

        l = pyphot.get_library()

        with open('table.csv', 'w') as output:
            output.write(','.join(hdr) + '\n')

            for k in sorted(l.content):
                fk = l[k]
                rec = (fk.name, fk.dtype, fk.wavelength_unit,
                    fk.cl.magnitude, fk.lpivot.magnitude, fk.leff.magnitude,
                    fk.Vega_zero_mag, fk.Vega_zero_flux.magnitude, fk.Vega_zero_Jy.magnitude,
                    fk.AB_zero_mag, fk.AB_zero_flux.magnitude, fk.AB_zero_Jy.magnitude,
                    fk.ST_zero_mag, fk.ST_zero_flux.magnitude, fk.ST_zero_Jy.magnitude)
                output.write(fmt.format(*rec))

    def test_all_indices_in_library(self):
        # convert to magnitudes
        from pyphot.sandbox import UnitLickLibrary as LickLibrary
        from pyphot.vega import Vega

        vega = Vega()
        # using the internal collection of indices
        lib = LickLibrary()
        f = lib['CN_1']
        # work on many spectra at once
        index = f.get(vega.wavelength, vega.flux, axis=-1)
        print("The index of Vega in {0:s} is {1:f} {2:s}".format(f.name, index, f.index_unit))

        # define header and table format (as csv)
        hdr = ("name", "wavelength units", "index units", "min", "max" "min blue", "max blue", "min red", "max red")
        fmt = "{0:s},{1:s},{2:s},{3:.3f},{4:.3f},{5:.3f},{6:.5f},{7:.3f},{8:.3f}\n"

        with open('licks_table.csv', 'w') as output:
            output.write(','.join(hdr) + '\n')

            for k in sorted(lib.content):
                fk = lib[k]
                # wavelength have units
                band = fk.band.magnitude
                blue = fk.blue.magnitude
                red = fk.red.magnitude
                rec = (fk.name, fk.wavelength_unit, fk.index_unit, band[0], band[1],
                    blue[0], blue[1], red[0], red[1])
                output.write(fmt.format(*rec))

    def test_Sun(self):
        import numpy as np
        from pyphot import sandbox as pyphot
        from pyphot import Sun
        from pyphot import unit
        sun_obs = Sun(flavor='observed')
        sun_th = Sun()   # default is theoric spectrum
        sun_th_10pc = Sun(distance=10 * unit['pc'])
        lib = pyphot.get_library()
        f = lib['GROUND_JOHNSON_V']
        expectations = -26.76, -26.76, +4.81
        for name, sun, expect_mag in zip(('observed', 'theoretical', 'th. 10pc'), (sun_obs,sun_th, sun_th_10pc), expectations):
            flux = f.get_flux(sun.wavelength, sun.flux)
            vegamag = f.Vega_zero_mag
            print('{0:12s} {1:0.5e} {2:+3.4f}'.format(name, flux.magnitude, -2.5 * np.log10(flux.magnitude) - vegamag))
            self.assertAlmostEqual(-2.5 * np.log10(flux.magnitude) - vegamag, expect_mag, places=2)

        filter_names = ['GROUND_JOHNSON_B', 'GROUND_JOHNSON_V', 'GROUND_BESSELL_J', 'GROUND_BESSELL_K']
        filter_names +=  lib.find('GaiaDR2')

        filters = lib.load_filters(filter_names, lamb=sun_th.wavelength)
        mags = {}
        for name, fn in zip(filter_names, filters):
            flux = fn.get_flux(sun_th.wavelength, sun_th.flux)
            vegamag = fn.Vega_zero_mag
            mag = -2.5 * np.log10(flux.magnitude) - vegamag
            mags[name] = mag
            print('{0:>25s} {1:+3.4f} mag'.format(name, mag))
        colors = (('GROUND_JOHNSON_B', 'GROUND_JOHNSON_V'),
                  ('GROUND_JOHNSON_V', 'GROUND_BESSELL_K'),
                  ('GROUND_BESSELL_J', 'GROUND_BESSELL_K'),
                  ('GaiaDR2_BP', 'GaiaDR2_RP'),
                  ('GaiaDR2_BP', 'GaiaDR2_G'),
                  ('GaiaDR2_G', 'GaiaDR2_RP'),
                  ('GaiaDR2v2_BP', 'GaiaDR2v2_RP'),
                  ('GaiaDR2v2_BP', 'GaiaDR2v2_G'),
                  ('GaiaDR2v2_G', 'GaiaDR2v2_RP'),
                  ('GaiaDR2_weiler_BPbright', 'GaiaDR2_weiler_RP'),
                  ('GaiaDR2_weiler_BPfaint', 'GaiaDR2_weiler_RP'),
                  ('GaiaDR2_weiler_BPbright', 'GaiaDR2_weiler_G'),
                  ('GaiaDR2_weiler_BPfaint', 'GaiaDR2_weiler_G'),
                  ('GaiaDR2_weiler_G', 'GaiaDR2_weiler_RP'))

        color_values = {}

        for color in colors:
            color_values[color] = mags[color[0]] - mags[color[1]]
            print('{0:>25s} - {1:<25s} = {2:3.4f} mag'.format(color[0], color[1], mags[color[0]] - mags[color[1]]))


class TestAstropyQuick(unittest.TestCase):
    def test_dependencies(self):
        from pyphot import astropy as pyphot

    def test_library_irac(self):
        from pyphot import astropy as pyphot
        # get the internal default library of passbands filters
        lib = pyphot.get_library()
        print("Library contains: ", len(lib), " filters")
        # find all filter names that relates to IRAC
        # and print some info
        f = lib.find('irac')
        for name in f:
            lib[name].info(show_zeropoints=True)

    def test_vega(self):
        from pyphot.astropy import Vega
        vega = Vega()

    def test_convert_mags(self):
        import numpy as np
        from pyphot.astropy import Vega
        from pyphot import astropy as pyphot
        # get the internal default library of passbands filters
        lib = pyphot.get_library()
        vega = Vega()
        f = lib['HST_WFC3_F110W']
        # compute the integrated flux through the filter f
        # note that it work on many spectra at once
        fluxes = f.get_flux(vega.wavelength, vega.flux, axis=-1)

        # Note that fluxes is now with units of erg/s/cm2/AA
        # pyphot gives Vega in flam and can convert between flux density units.
        print(fluxes, vega.wavelength, vega.flux)

        # convert to vega magnitudes
        mags = -2.5 * np.log10(fluxes.value) - f.Vega_zero_mag
        print("Vega value of Vega in {0:s} is : {1:f} mag".format(f.name, mags))
        mags = -2.5 * np.log10(fluxes.value) - f.AB_zero_mag
        print("AB value of Vega in {0:s} is : {1:f} mag".format(f.name, mags))
        mags = -2.5 * np.log10(fluxes.value) - f.ST_zero_mag
        print("ST value of Vega in {0:s} is : {1:f} mag".format(f.name, mags))

    def test_load_all_filters(self):
        from pyphot import astropy as pyphot
        # get the internal default library of passbands filters
        lib = pyphot.get_library()
        lib.load_all_filters()

    def test_all_filters_in_library(self):
        from pyphot import astropy as pyphot

        # define header and table format (as csv)
        hdr = ("name", "detector type", "wavelength units",
            "central wavelength", "pivot wavelength", "effective wavelength",
            "Vega mag", "Vega flux", "Vega Jy",
            "AB mag", "AB flux", "AB Jy",
            "ST mag", "ST flux", "ST Jy")
        fmt = "{0:s},{1:s},{2:s},{3:.3f},{4:.3f},{5:.3f},{6:.5f},{7:.5g},{8:.5g},{9:.5f},{10:.5g},{11:.5g},{12:.5f},{13:.5g},{14:.5g}\n"

        l = pyphot.get_library()

        with open('table.csv', 'w') as output:
            output.write(','.join(hdr) + '\n')

            for k in sorted(l.content):
                fk = l[k]
                rec = (fk.name, fk.dtype, fk.wavelength_unit,
                    fk.cl.value, fk.lpivot.value, fk.leff.value,
                    fk.Vega_zero_mag, fk.Vega_zero_flux.value, fk.Vega_zero_Jy.value,
                    fk.AB_zero_mag, fk.AB_zero_flux.value, fk.AB_zero_Jy.value,
                    fk.ST_zero_mag, fk.ST_zero_flux.value, fk.ST_zero_Jy.value)
                output.write(fmt.format(*rec))

    def test_all_indices_in_library(self):
        # convert to magnitudes
        from pyphot.astropy import UnitLickLibrary as LickLibrary
        from pyphot.astropy import Vega

        vega = Vega()
        # using the internal collection of indices
        lib = LickLibrary()
        f = lib['CN_1']
        # work on many spectra at once
        index = f.get(vega.wavelength, vega.flux, axis=-1)
        print("The index of Vega in {0:s} is {1:f} {2:s}".format(f.name, index, f.index_unit))

        # define header and table format (as csv)
        hdr = ("name", "wavelength units", "index units", "min", "max" "min blue", "max blue", "min red", "max red")
        fmt = "{0:s},{1:s},{2:s},{3:.3f},{4:.3f},{5:.3f},{6:.5f},{7:.3f},{8:.3f}\n"

        with open('licks_table.csv', 'w') as output:
            output.write(','.join(hdr) + '\n')

            for k in sorted(lib.content):
                fk = lib[k]
                # wavelength have units
                band = fk.band.value
                blue = fk.blue.value
                red = fk.red.value
                rec = (fk.name, fk.wavelength_unit, fk.index_unit, band[0], band[1],
                    blue[0], blue[1], red[0], red[1])
                output.write(fmt.format(*rec))

    def test_Sun(self):
        import numpy as np
        from pyphot import astropy as pyphot
        from pyphot.astropy import Sun
        from pyphot.astropy import Unit
        sun_obs = Sun(flavor='observed')
        sun_th = Sun()   # default is theoric spectrum
        sun_th_10pc = Sun(distance=10 * Unit('pc'))
        lib = pyphot.get_library()
        f = lib['GROUND_JOHNSON_V']
        expectations = -26.76, -26.76, +4.81
        for name, sun, expect_mag in zip(('observed', 'theoretical', 'th. 10pc'), (sun_obs,sun_th, sun_th_10pc), expectations):
            flux = f.get_flux(sun.wavelength, sun.flux)
            vegamag = f.Vega_zero_mag
            print('{0:12s} {1:0.5e} {2:+3.4f}'.format(name, flux.value, -2.5 * np.log10(flux.value) - vegamag))
            self.assertAlmostEqual(-2.5 * np.log10(flux.value) - vegamag, expect_mag, places=2)

        filter_names = ['GROUND_JOHNSON_B', 'GROUND_JOHNSON_V', 'GROUND_BESSELL_J', 'GROUND_BESSELL_K']
        filter_names +=  lib.find('GaiaDR2')

        filters = lib.load_filters(filter_names, lamb=sun_th.wavelength)
        mags = {}
        for name, fn in zip(filter_names, filters):
            flux = fn.get_flux(sun_th.wavelength, sun_th.flux)
            vegamag = fn.Vega_zero_mag
            mag = -2.5 * np.log10(flux.value) - vegamag
            mags[name] = mag
            print('{0:>25s} {1:+3.4f} mag'.format(name, mag))
        colors = (('GROUND_JOHNSON_B', 'GROUND_JOHNSON_V'),
                  ('GROUND_JOHNSON_V', 'GROUND_BESSELL_K'),
                  ('GROUND_BESSELL_J', 'GROUND_BESSELL_K'),
                  ('GaiaDR2_BP', 'GaiaDR2_RP'),
                  ('GaiaDR2_BP', 'GaiaDR2_G'),
                  ('GaiaDR2_G', 'GaiaDR2_RP'),
                  ('GaiaDR2v2_BP', 'GaiaDR2v2_RP'),
                  ('GaiaDR2v2_BP', 'GaiaDR2v2_G'),
                  ('GaiaDR2v2_G', 'GaiaDR2v2_RP'),
                  ('GaiaDR2_weiler_BPbright', 'GaiaDR2_weiler_RP'),
                  ('GaiaDR2_weiler_BPfaint', 'GaiaDR2_weiler_RP'),
                  ('GaiaDR2_weiler_BPbright', 'GaiaDR2_weiler_G'),
                  ('GaiaDR2_weiler_BPfaint', 'GaiaDR2_weiler_G'),
                  ('GaiaDR2_weiler_G', 'GaiaDR2_weiler_RP'))

        color_values = {}

        for color in colors:
            color_values[color] = mags[color[0]] - mags[color[1]]
            print('{0:>25s} - {1:<25s} = {2:3.4f} mag'.format(color[0], color[1], mags[color[0]] - mags[color[1]]))


if __name__ == '__main__':
    unittest.main()