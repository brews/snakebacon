import logging
import numpy as np
import scipy.stats as stats
import pandas as pd


log = logging.getLogger(__name__)


class Core:

    def __init__(self, age, error, depth, labid, depth_units='meters'):
        # TODO(brews): Add support for `pint` unit handling. Note that Bacon uses cm for depth.
        self.labid = labid
        self.age = age
        self.error = error
        self.depth = depth

    def suggest_accumulation_rate(self):
        """From core age-depth data, suggest accum. rate
        Follow's Bacon's method @ Bacon.R ln 30 - 44
        """
        # Suggested round vals.
        sugg = np.tile([1, 2, 5], (4, 1)) * np.reshape(np.repeat([0.1, 1.0, 10, 100], 3), (4, 3))
        # Get ballpark accumulation rates, uncalibrated dates.
        ballpacc = stats.linregress(x=self.depth, y=self.age * 1.1).slope
        ballpacc = np.abs(sugg - ballpacc)
        sugg = sugg.flat[ballpacc.argmin()]  # Suggest rounded acc.rate with lowest abs diff.
        return sugg

    def suggest_thick(self, reswarn, thick=5, d_min=None, d_max=None):
        """ Bacon.R lines #76 - #87
        """
        # TODO(brews): Missing py equivalent to R's `pretty()`, then can finish. See https://stackoverflow.com/questions/43075617/python-function-equivalent-to-rs-pretty
        # @ Bacon.R line #72.
        # "assign depths, possibly suggest alternative value for thick"
        # Bacon.R now checks to see if we want suggested values @ line #76
        # thick, d, k = suggest_thick(k, d, reswarn, thick, d_min, d_max)
        if d_min is None:
            d_min = self.depth.min()
        if d_max is None:
            d_max = self.depth.max()
        d = np.arange(np.floor(d_min), np.ceil(d_max), thick)
        k = len(d)
        # ans = 'n'  # brews: What is this for?
        # if len(reswarn) == 2:
        #     if k < np.min(reswarn):
        #         # Stopped on line #80 in Bacon.R. Need to write 'pretty()' in python.
        #         sugg = np.min( thick * (k/np.min(reswarn)) )
        #     elif k > np.max(reswarn):
        #         pass
        # thick = sugg
        # d_out = np.arange(np.floor(d_min), np.ceil(d_max), thick)
        # k_out = len(d_out)
        # return (thick_out, d_out, k_out)

    def calibrate_dates(self, calib_curve, d_r, d_std, t_a=3, t_b=4, cutoff=0.001, normal_distr=False):
        """Python version of .bacon.calib() on line 908 in Bacon.R
        """
        # .bacon.calib - line 908

        # rcmean = 4128; w2 = 4225; t_a=3; t_b=4
        # test = d_cal(cc = calib_curve.rename(columns = {0:'a', 1:'b', 2:'c'}), rcmean = 4128, w2 = 4225, t_a=t_a, t_b=t_b, cutoff=cutoff, normal = normal)

        # Line 959 of Bacon.R
        # calib = list(dets.iloc[:, 3])
        # Now Bacon goes and checks the ncol in the dets See line #960 in Bacon.R

        # Line #973
        calib_probs = []
        # I think we can do the below without a loop.
        for i in range(len(self.depth)):
            # TODO(brews): Rename columns, or have the columns passed in with names.
            age_realizations = calib_curve.d_cal(rcmean=self.age[i] - d_r, w2=self.error[i] ** 2 + d_std ** 2,
                                                 t_a=t_a, t_b=t_b, cutoff=cutoff, normal_distr=normal_distr)
            calib_probs.append(age_realizations)
        return self.depth, calib_probs


def read_corefile(fl):
    """Create proxy instance from Bacon proxy file
    """
    indata = pd.read_table(fl, sep=r'\,\s*', index_col=None, engine='python')
    outcore = Core(age=indata['age'].values,
                   error=indata['error'].values,
                   depth=indata['depth'].values,
                   labid=indata['labID'].values,
                   depth_units='meters')
    return outcore
