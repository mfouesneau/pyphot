
Provided Filter library
=======================

This page shows the content of the provided library with respective properties
of the passband filters. The code to generate the table is also provided below.

.. important::

  The internal library is not exhaustive and is meant to be a starting point for
  users to build their own library or contribute to this one.

  :class:`~pyphot.svo` provides an interface to the `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_ which provides a more exhaustive list of passbands.

  .. code-block:: python

        from pyphot import svo 
        lst = ("2MASS/2MASS.J", "2MASS/2MASS.H")
        pbands = [svo.get_pyphot_filter(k) for k in lst]

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

