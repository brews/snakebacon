import logging as logging
import numpy as np
import matplotlib.pylab as plt

from .mcmc import McmcResults


log = logging.getLogger(__name__)


class AgeDepthModel:

    def __init__(self, *args, **kwargs):
        self.mcmcfit = McmcResults(*args, **kwargs)
        try:
            self.mcmcfit.burnin(kwargs['burnin'])
        except KeyError:
            self.mcmcfit.burnin(200)
        self.thick = (kwargs['depth_max'] - kwargs['depth_min']) / kwargs['k']
        self.depth = np.arange(kwargs['depth_min'], kwargs['depth_max'] + 1)
        self.cal_age = np.array([self.agedepth(d=dx) for dx in self.depth])
        # ds = np.array([agedepth(d=dx, x=test['x'], deltac=5, x0=test['theta'], c0=1.5) for dx in d])

    def date(self, proxy):
        """Date a proxy record"""
        pass

    def plot(self):
        """Age-depth plot"""
        plt.hist2d(np.repeat(self.depth, self.cal_age.shape[1]), self.cal_age.flatten(), (len(self.depth), 50), cmin=1)
        plt.step(self.depth, np.median(self.cal_age, axis=1), where='mid', color='red')
        lower_env = np.percentile(self.cal_age, q=2.5, axis=1)
        upper_env = np.percentile(self.cal_age, q=97.5, axis=1)
        plt.step(self.depth, upper_env, where='mid', color='red', linestyle=':')
        plt.step(self.depth, lower_env, where='mid', color='red', linestyle=':')
        plt.ylabel('Age (cal yr BP)')
        plt.xlabel('Depth (cm)')
        plt.grid()


    def agedepth(self, d):
        """Get calendar age for a depth

        Parameters
        ----------
        d : float
            Sediment depth (in cm).

        Returns
        -------
        Numeric giving true age at given depth.
        """
        # TODO(brews): Function cannot handle hiatus
        # See lines 77 - 100 of hist2.cpp
        x = self.mcmcfit.sediment_rate
        theta0 = self.mcmcfit.headage  # Age abscissa (in yrs).  If array, dimension should be iterations or realizations of the sediment
        deltac = self.thick
        c0 = min(self.depth) # Uniform depth segment abscissa (in cm).
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
