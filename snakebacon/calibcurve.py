import logging
import numpy as np
import scipy.stats as stats
import pandas as pd

log = logging.getLogger(__name__)


class CalibCurve:

    def __init__(self, calbp, c14age, error, delta14c, sigma):
        self.calbp = calbp
        self.c14age = c14age
        self.error = error
        self.delta14c = delta14c
        self.sigma = sigma

    def d_cal(self, rcmean, w2, t_a=3, t_b=4, cutoff=0.001, normal_distr=False):
        """ Line 943 of Bacon.R
        cc : calib_curve (3-col format)
        rcmean : det['age'][i] - d_R
        w2 : dat['error'][i]^2 + d_STD**2
        """
        if normal_distr:
            # TODO(brews): Test this. Line 946 of Bacon.R.
            std = np.sqrt(self.error ** 2 + w2)
            dens = stats.norm(loc=rcmean, scale=std).pdf(self.c14age)
        else:
            # TODO(brews): Test this. Line 947 of Bacon.R.
            dens = (t_b + ((rcmean - self.c14age) ** 2) / (2 * (self.error ** 2 + w2))) ** (-1 * (t_a + 0.5))
        cal = np.array([self.calbp.copy(), dens]).T
        cal[:, 1] = cal[:, 1] / cal[:, 1].sum()
        # "ensure that also very precise dates get a range of probabilities"
        cutoff_mask = cal[:, 1] > cutoff
        if cutoff_mask.sum() > 5:
            out = cal[cutoff_mask, :]
        else:
            calx = np.linspace(cal[:, 0].min(), cal[:, 0].max(), num=50)
            caly = np.interp(calx, cal[:, 0], cal[:, 1])
            out = np.array([calx, caly / caly.sum()]).T
        return out


def read_14c(fl):
    """Create CalibCurve instance from Bacon curve file
    """
    indata = pd.read_csv(fl, index_col=None, skiprows=11, header=None,
                         names=['calbp', 'c14age', 'error', 'delta14c', 'sigma'])
    outcurve = CalibCurve(calbp=indata['calbp'].values,
                          c14age=indata['c14age'].values,
                          error=indata['error'].values,
                          delta14c=indata['delta14c'].values,
                          sigma='sigma')
    return outcurve
