"""
Cython wrapper for the C/C++ from Bacon age-modelling software (http://chrono.qub.ac.uk/blaauw/bacon.html) written by 
Maarten Blaauw (maarten.blaauw@qub.ac.uk) and  Andres Christen (jac@cimat.mx).
"""

from .baconwrap import run_baconmcmc
from .calibcurves import fetch_calibcurve
from .utils import calibrate_dates