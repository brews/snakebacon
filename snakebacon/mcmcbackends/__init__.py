import numpy as np
import scipy.stats as stats

from .bacon import run_baconmcmc, fetch_calibcurve, calibrate_dates


class Bacon:
    def runmcmc(*args, **kwargs):
        return run_baconmcmc(*args, **kwargs)

    def prior_dates(*args, **kwargs):
        """Get the prior distribution of calibrated radiocarbon dates"""
        try:
            chron = args[0]
        except IndexError:
            chron = kwargs['coredates']

        d_r = np.array(kwargs['d_r'])
        d_std = np.array(kwargs['d_std'])
        t_a = np.array(kwargs['t_a'])
        t_b = np.array(kwargs['t_b'])

        try:
            normal_distr = kwargs['normal_distr']
        except KeyError:
            normal_distr = None

        cc_int = kwargs['cc']

        ccdict = {0: 'ConstCal', 1: 'IntCal3', 2: 'Marine13',
                  3: 'SHCal13', 4: 'ConstCal'}
        # There is a better way to do this.
        if 'cc1' in kwargs:
            ccdict[1] = str(kwargs['cc1'])
        if 'cc2' in kwargs:
            ccdict[2] = str(kwargs['cc2'])
        if 'cc3' in kwargs:
            ccdict[3] = str(kwargs['cc3'])
        if 'cc4' in kwargs:
            ccdict[4] = str(kwargs['cc4'])

        cc = []
        for i in cc_int:
            i = int(i)
            cc.append(fetch_calibcurve(ccdict[i]))

        d, p = calibrate_dates(chron, calib_curve=cc, d_r=d_r, d_std=d_std,
                               t_a=t_a, t_b=t_b, normal_distr=normal_distr)
        return d, p


    def prior_sediment_rate(*args, **kwargs):
        """Get the prior density of sediment rates

        Returns
        -------
        y : ndarray
            Array giving the density.
        x : ndarray
            Array of sediment accumulation values (yr/cm) over which the density was evaluated.
        """
        # PlotAccPrior @ Bacon.R ln 113 -> ln 1097-1115
        # alpha = acc_shape, beta = acc_shape / acc_mean
        # TODO(brews): Check that these stats are correctly translated to scipy.stats distribs.
        acc_mean = kwargs['acc_mean']
        acc_shape = kwargs['acc_shape']
        x = np.linspace(0, 6 * np.max(acc_mean), 100)
        y = stats.gamma.pdf(x, a=acc_shape,
                            scale=1 / (acc_shape/acc_mean))
        return y, x


    def prior_sediment_memory(*args, **kwargs):
        """Get the prior density of sediment memory

        Returns
        -------
        y : ndarray
            Array giving the density.
        x : ndarray
            Array of Memory (ratio) values over which the density was evaluated.
        """
        # "plot the prior for the memory (= accumulation rate varibility between neighbouring depths)"
        # PlotMemPrior @ Bacon.R ln 114 -> ln 1119 - 1141
        # w_a = mem_strength * mem_mean, w_b = mem_strength * (1 - mem_mean)
        # TODO(brews): Check that these stats are correctly translated to scipy.stats distribs.
        mem_shape = kwargs['mem_strength']  # aka. `mem_shape`
        mem_mean = kwargs['mem_mean']
        x = np.linspace(0, 1, 100)
        y = stats.beta.pdf(x, a=mem_shape * mem_mean,
                           b=mem_shape * (1 - mem_mean))
        return y, x

