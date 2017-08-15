import logging as logging

import matplotlib.pylab as plt
import numpy as np

from .mcmc import McmcSetup
from .records import DatedProxyRecord


log = logging.getLogger(__name__)


class AgeDepthModel:
    def __init__(self, coredates, *, mcmc_kwargs, hold=False, burnin=200):
        self.burnin = int(burnin)
        self.mcmcsetup = McmcSetup(coredates, **mcmc_kwargs)
        self._mcmcfit = None
        self._thick = None
        self._depth = None
        self._age_ensemble = None
        if not hold:
            self.fit()

    @property
    def mcmcfit(self):
        if self._mcmcfit is not None:
            return self._mcmcfit
        else:
            raise NeedFitError('AgeDepthModel instance needs to be fit() first')

    @property
    def thick(self):
        if self._thick is not None:
            return self._thick
        else:
            raise NeedFitError('AgeDepthModel instance needs to be fit() first')

    @property
    def depth(self):
        if self._depth is not None:
            return self._depth
        else:
            raise NeedFitError('AgeDepthModel instance needs to be fit() first')

    @property
    def age_ensemble(self):
        if self._age_ensemble is not None:
            return self._age_ensemble
        else:
            raise NeedFitError('AgeDepthModel instance needs to be fit() first')

    def __repr__(self):
        return '%s(coredates=%r, mcmc_kwargs=%r, burnin=%r)' % (type(self).__name__, self.mcmcsetup.coredates, self.mcmcsetup.mcmc_kwargs, self.burnin)

    def age_median(self):
        return np.median(self.age_ensemble, axis=1)

    def age_percentile(self, p):
        return np.percentile(self.age_ensemble, q=p, axis=1)

    def fit(self):
        """Fit MCMC AgeDepthModel"""
        self._mcmcfit = self.mcmcsetup.run()
        self._mcmcfit.burnin(self.burnin)
        dmin = min(self._mcmcfit.depth_segments)
        dmax = max(self._mcmcfit.depth_segments)
        self._thick = (dmax - dmin) / len(self.mcmcfit.depth_segments)
        self._depth = np.arange(dmin, dmax + 1)
        self._age_ensemble = np.array([self.agedepth(d=dx) for dx in self.depth])

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
        DatedProxyRecord
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

    def plot(self, agebins=50, p=(2.5, 97.5), ax=None):
        """Age-depth plot"""
        if ax is None:
            ax = plt.gca()
        ax.hist2d(np.repeat(self.depth, self.age_ensemble.shape[1]), self.age_ensemble.flatten(),
                   (len(self.depth), agebins), cmin=1)
        ax.step(self.depth, self.age_median(), where='mid', color='red')
        ax.step(self.depth, self.age_percentile(p[0]), where='mid', color='red', linestyle=':')
        ax.step(self.depth, self.age_percentile(p[1]), where='mid', color='red', linestyle=':')
        ax.set_ylabel('Age (cal yr BP)')
        ax.set_xlabel('Depth (cm)')
        ax.grid()
        return ax

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
        c0 = min(self.depth)  # Uniform depth segment abscissa (in cm).
        assert d > c0 or np.isclose(c0, d, atol = 1e-4)
        out = theta0.astype(float)
        i = int(np.floor((d - c0) / deltac))
        for j in range(i):
            out += x[j] * deltac
        ci = c0 + i * deltac
        assert ci < d or np.isclose(ci, d, atol = 1e-4)
        try:
            next_x = x[i]
        except IndexError:
            # Extrapolating
            next_x = x[i - 1]
        out += next_x * (d - ci)
        return out

    def prior_dates(self):
        self.mcmcsetup.prior_dates()

    def prior_sediment_rate(self):
        self.mcmcsetup.prior_sediment_rate()

    def prior_sediment_memory(self):
        self.mcmcsetup.prior_sediment_memory()


class NeedFitError(Exception):
    pass
