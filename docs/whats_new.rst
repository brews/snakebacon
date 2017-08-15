.. currentmodule:: snakebacon

What's New
==========

.. _whats-new.0.0.2:

v0.0.2
-------------------

Enhancements
~~~~~~~~~~~~
- DateRecord class has been renamed ChronRecord to avoid confusion with DatedProxyRecord (`Issue #3 <https://github.com/brews/snakebacon/issues/3>`_). The read_chron() is also now read_chron().

- AgeDepthModels now can plot prior distributions.

- Can now query AgeDepthModels to get prior distribution for chronology dates, sediment rates, and sediment rate memory (`Issue #5 <https://github.com/brews/snakebacon/issues/5>`_).

- Now can use snakebacon.mcmcbackends.bacon.fetch_calibcurve() and the @registercurve() decorator to define and fetch calibration curves by name (`Issue #6 <https://github.com/brews/snakebacon/issues/6>`_).

- McmcSetup and McmcResults now throw bacon its tasks through mcmcbackends.Bacon. This should make snakebacon more bacon-agnostic. Have a look at McmcSetup's mcmcbackend arg.

- suggest_accumulation_rate() is no longer a method, but a humble function.

Bug fixes
~~~~~~~~~
- Fixed common AssertionError in AgeDepthModel.agedepth caused by float precision bug (`Issue #4 <https://github.com/brews/snakebacon/issues/4>`_).

- Fixed problem with ProxyRecords and DatedProxyRecords not making copies of Pandas DataFrame attributes. (`Issue #2 <https://github.com/brews/snakebacon/issues/2>`_).

.. _whats-new.0.0.1:

v0.0.1
-------------------

This is the initial, alpha release of snakebacon.
