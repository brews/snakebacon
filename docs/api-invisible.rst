.. Generate API reference pages, but don't display these in tables.
.. This extra page is a work around for sphinx not having any support for
.. hiding an autosummary table. Props to xarray for help with this.

.. currentmodule:: snakebacon

.. autosummary::
   :toctree: generated/


   CalibCurve
   DatedProxyRecord
   ChronRecord
   ProxyRecord

   read_14c
   read_chron
   read_proxy

