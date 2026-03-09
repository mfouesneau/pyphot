Internal Sun reference spectra
==============================

Although the Sun is not a necessary reference for photometric calculations,  we provide observed and theoretical references for the solar spectrum following `Colina, Bohlin, & Castelli 1996 <https://ui.adsabs.harvard.edu/abs/1996AJ....112..307C/abstract>`_.

The observed solar spectrum comes from CALSPEC `sun_reference_stis_001.fits <ftp://ftp.stsci.edu/cdbs/current_calspec/sun_reference_stis_001.fits>`_ which provides the ultraviolet to near-infrared absolute flux distribution of the Sun covering the 0.12-2.5 μm wavelength range. The solar reference spectrum combines absolute flux measurements from satellites and from the ground with a model spectrum for the near-infrared.

The theoretical spectrum comes from the Kurucz'93 atlas: `<sun_kurucz93.fits ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/sun_kurucz93.fits>`_ The theoretical spectrum is scaled to match the observed spectrum from 1.5 - 2.5 microns, and then it is used where the observed spectrum ends. The theoretical model of the Sun uses the following parameters when the Sun is at 1 au:

.. list-table:: Solar Parameters
        :header-rows: 1

        * - log_Z
          - T_eff
          - log_g
          - V_{Johnson}
        * - +0.0
          - 5777
          - +4.44
          - -26.75

The interface to the Sun templates is given through the :class:`pyphot.sun.Sun` class.
