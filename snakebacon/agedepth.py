import logging as logging
from .mcmc import McmcResults


log = logging.getLogger(__name__)


class AgeDepthModel:

    def __init__(self):
        pass

    def date(self):
        pass

    def plot(self):
        pass


def agedepth(d, x, deltac, theta0, c0):
    """Get true age for a depth

    Parameters
    ----------
    d : float
        Sediment depth (in cm).
    x : 1 or 2darray
        i-length array of sedimentation rates (yr/cm). Can also be (i, j) array where i is along sediment core segments
        and j is iterations or realizations of the core.
    deltac : float
        Change in depth for a uniform depth segments (cm).
    theta0 : float or 1darray
        Age abscissa (in yrs).  If array, dimension should be iterations or realizations of the sediment
        core.
    c0 : Uniform depth segment abscissa (in cm).

    Returns
    -------
    Numeric giving true age at given depth.
    """
    # TODO(brews): Funciton needs to be tested. Carefully.
    # TODO(brews): Function cannot handle hiatus
    # See lines 77 - 100 of hist2.cpp
    assert d >= c0
    out = theta0.copy()
    i = int(np.floor((d - c0) / deltac))
    for j in range(i):
        out += x[j] * deltac
    ci = c0 + i * deltac
    assert ci <= d
    try:
        next_x = x[i]
    except IndexError:
        # Extrapolating
        next_x = x[i - 1]
    out += next_x * (d - ci)
    return out
