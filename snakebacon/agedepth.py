import logging as logging
from copy import deepcopy
import numpy as np
import matplotlib.pylab as plt

from .mcmc import McmcResults
from .records import DatedProxyRecord


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
        self.age_ensemble = np.array([self.agedepth(d=dx) for dx in self.depth])
        self.age_median = np.median(self.age_ensemble, axis=1)
        self.conf_interv = {2.5:np.percentile(self.age_ensemble, q=2.5, axis=1),
                            97.5:np.percentile(self.age_ensemble, q=97.5, axis=1)}

    def date(self, proxy, how='median', n=500):
        """Date a proxy record

        Parameters
        ----------
        proxy : ProxyRecord
        how : str
            How to perform the dating. 'median' returns the average of the MCMC ensemble. 'ensemble' returns a 'n'
            randomly selected members of the MCMC ensemble. Default is 'median'.
        n : int
            If 'how' is 'ensemble', the function will randomly select 'n' MCMC ensemble members, with replacement.

        Returns
        -------
        A DatedProxyRecord, copied from 'proxy'.
        """
        assert how in ['median', 'ensemble']
        ens_members = self.mcmcfit.n_members()
        if how == 'ensemble':
            select_idx = np.random.choice(range(ens_members), size=n, replace=True)
        out = []
        for d in proxy.data.depth.values:
            age = self.agedepth(d)
            if how == 'median':
                age = np.median(age)
            elif how == 'ensemble':
                age = age[select_idx]
            out.append(age)
        return DatedProxyRecord(proxy.data.copy(), out)

    def plot(self, agebins=50):
        """Age-depth plot"""
        plt.hist2d(np.repeat(self.depth, self.age_ensemble.shape[1]), self.age_ensemble.flatten(), (len(self.depth), agebins), cmin=1)
        plt.step(self.depth, self.age_median, where='mid', color='red')
        plt.step(self.depth, self.conf_interv[97.5], where='mid', color='red', linestyle=':')
        plt.step(self.depth, self.conf_interv[2.5], where='mid', color='red', linestyle=':')
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
