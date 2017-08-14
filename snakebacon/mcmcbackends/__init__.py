import numpy as np

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
        """Get the prior distribution of sediment rates"""
        pass

    def prior_sediment_memory(*args, **kwargs):
        """Get the prior distribution of sediment memory"""
        pass
