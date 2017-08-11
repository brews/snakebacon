import os

from .Curves import here as curvespath
from snakebacon.records import read_14c


available_curves = dict()


def fetch_calibcurve(curvename):
    """Get CalibCurve from name string"""
    f = available_curves[curvename]
    return f()


def registercurve(curvename):
    """Decorator to register functions returning CalibCurves"""
    def decor(func):
        available_curves[curvename] = func
        return func
    return decor


@registercurve('IntCal13')
def fetch_intcal13():
    return read_14c(os.path.join(curvespath, 'intcal13.14C'))


@registercurve('Marine13')
def fetch_marine13():
    return read_14c(os.path.join(curvespath, 'marine13.14C'))


@registercurve('SHCal13')
def fetch_shcal13():
    return read_14c(os.path.join(curvespath, 'shcal13.14C'))
