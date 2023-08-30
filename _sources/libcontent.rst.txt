
Provided Filter library
=======================

This page shows the content of the provided library with respective properties
of the passband filters. The code to generate the table is also provided below.


.. literalinclude:: pyphot_table.py
   :language: python

Table description
-----------------

* `name`:  the identification name of the filter in the library.
* `detector type`: energy or photon counter.
* `wavelength units`:  filter defined with these units and all wavelength
  properties: `central wavelength`, `pivot wavelength`, and `effective wavelength`.
* `<X> mag`: magnitude in Vega, AB or ST system (w.r.t. the detector type)
* `<X> flux`: flux in :math:`erg/s/cm^2/\unicode{x212B}` in the `X` system
* `<X> Jy`: flux in :math:`Jy` (Jansky) in the `X` system


.. csv-table:: Current internal library
        :file: table.csv

