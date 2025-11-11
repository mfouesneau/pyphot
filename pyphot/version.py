"""Expose the version of pyphot the package metadata

.. code-block:: python

    >>> import pyphot
    >>> pyphot.__VERSION__
    '2.0.0'

"""

import importlib.metadata

__VERSION__ = importlib.metadata.version("pyphot")
