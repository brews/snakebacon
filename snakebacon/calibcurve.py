import logging
import numpy as np
import scipy.stats as stats
import pandas as pd

log = logging.getLogger(__name__)


class CalibCurve:
    """A calibration curve
    """

    def __init__(self, calbp, c14age, error, delta14c=None, sigma=None):
        """Create a calibration curve instance

        Parameters
        ----------
        calbp : ndarray
        c14age : ndarray
        error : ndarray
        delta14c : ndarray
        sigma : ndarray
        
        """
        self.calbp = calbp
        self.c14age = c14age
        self.error = error
        if delta14c is None:
            delta14c = np.zeros(calbp.shape)
        self.delta14c = delta14c  # d_R
        if sigma is None:
            sigma = np.zeros(calbp.shape)
        self.sigma = sigma # d_R variance?


    def c14age_from_age(self, theta):
        """Interpolate C14 mean age and variance from true age

        Parameters
        ----------
        theta : float
            A true age value in calendar years BP.

        Returns
        -------
        A tuple, (mu, std), giving the mean and standard deviation of corresponding C14 age, interpolated from the curve.

        Dertived from bacon @ cal.h lines 156 - 193.

        """
        idx = np.searchsorted(self.calbp, theta)
        mu = self.c14age[idx - 1] + (theta - self.calbp[idx - 1]) * (self.c14age[idx] - self.c14age[idx - 1]) / (
        self.c14age[idx] - self.calbp[idx - 1])
        sig = self.error[idx - 1] + (theta - self.calbp[idx - 1]) * (self.error[idx] - self.error[idx - 1]) / (
        self.c14age[idx] - self.calbp[idx - 1])
        return (mu, sig)


    def d_cal(self, rcmean, w2, cutoff=0.001, normal_distr=False, t_a=3, t_b=4):
        """Get calendar date probabilities
        
        Parameters
        ----------
        rcmean : scalar
            Reservoir-adjusted age.
        w2 : scalar
            r'$w^2_j(\theta)$' from pg 461 & 463 of Blaauw and Christen 2011.
        cutoff : scalar, optional
            Unknown.
        normal_distr : Bool, optional
            Use normal distribution for date errors. If False, then use Student's t-distribution.
        t_a : scalar, optional
            Student's t-distribution parameter, a. t_a - 1 must equal t_b.
        t_b : scalar, optional
            Student's t-distribution parameter, b. t_a - 1 must equal t_b.

        
        #Line 943 of Bacon.R
        #cc : calib_curve (3-col format)
        #rcmean : det['age'][i] - d_R
        #w2 : dat['error'][i]^2 + d_STD**2
        """
        assert t_a - 1 == t_b
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
