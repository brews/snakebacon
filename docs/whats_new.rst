.. currentmodule:: snakebacon

What's New
==========

.. _whats-new.0.0.2:

v0.0.2
-------------------

Enhancements
~~~~~~~~~~~~
- Use snakebacon.mcmcbackends.bacon.fetch_calibcurve() and the @registercurve() decorator to define and fetch calibration curves but name (`Issue #6 <https://github.com/brews/snakebacon/issues/6>`_).

- McmcSetup and McmcResults now throw bacon its tasks through mcmcbackends.Bacon. This should make snakebacon more bacon-agnostic. Have a look at McmcSetup's mcmcbackend arg.


Bug fixes
~~~~~~~~~
- Fixed common AssertionError in AgeDepthModel.agedepth caused by float precision bug (`Issue #4 <https://github.com/brews/snakebacon/issues/4>`_).


.. _whats-new.0.0.1:

v0.0.1
-------------------

This is the initial, alpha release of snakebacon.
