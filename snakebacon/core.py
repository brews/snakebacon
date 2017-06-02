import logging
import numpy as np
import scipy.stats as stats
import pandas as pd

import matplotlib.pylab as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


log = logging.getLogger(__name__)


class Core:
    """A sediment core
    """

    def __init__(self, age, error, depth, labid, depth_units='meters'):
        """Create a sediment Core instance
        
        Parameters
        ----------
        age : ndarray
        error : ndarray
        depth : ndarray
        labid : ndarray
        depth_units : string, optional
        
        """
        # TODO(brews): Add support for `pint` unit handling. Note that Bacon uses cm for depth.
        self.labid = labid
        self.age = age
        self.error = error  # Note this is called "std" in output .bacon file.
        self.depth = depth

    def suggest_accumulation_rate(self):
        """From core age-depth data, suggest accumulation rate
        
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
                
        Parameters
        ----------
        reswarn : Unknown
        thick : int, optional
            Sediment segment thickness.
        d_min : float, optional
            Minimum depth.
        d_max : float, optional
            Maximum depth.

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

    def calibrate_dates(self, calib_curve, d_r, d_std, cutoff=0.001, normal_distr=False, t_a=3, t_b=4):
        """Get probability of calendar dates for each depth segment in core
        
        Parameters
        ----------
        calib_curve : Curve
            Radiocarbon calibration curve.
        d_r : scalar or ndarray
            Carbon reservoir offset.
        d_std : scalar or ndarray
            Carbon reservoir offset error standard deviation.
        cutoff : scalar, optional
            Unknown.
        normal_distr : Bool, optional
            Use normal distribution for date errors. If False, then use Student's t-distribution.
        t_a : scalar, optional
            Student's t-distribution parameter, a. t_a - 1 must equal t_b.
        t_b : scalar, optional
            Student's t-distribution parameter, b. t_a - 1 must equal t_b.
            
        Returns
        -------
        out : (ndarray, list)
            out[0] is ndarray of sediment segment depths. out[1] is a list of ndarray of probabilities with one ndarray 
            per core segment. len(out[0]) == len(out[1])

        Python version of .bacon.calib() on line 908 in Bacon.R
        """
        # .bacon.calib - line 908

        # rcmean = 4128; w2 = 4225; t_a=3; t_b=4
        # test = d_cal(cc = calib_curve.rename(columns = {0:'a', 1:'b', 2:'c'}), rcmean = 4128, w2 = 4225, t_a=t_a,
        # t_b=t_b, cutoff=cutoff, normal = normal)

        # Line 959 of Bacon.R
        # calib = list(dets.iloc[:, 3])
        # Now Bacon goes and checks the ncol in the dets See line #960 in Bacon.R

        # Line #973
        assert t_b - 1 == t_a
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

def plot_acc_prior(mem_shape, mem_mean, thick):
    """Plot accumulation rate varibility between neighbouring depths ("memory") prior
    
    # PlotAccPrior @ Bacon.R ln 113 -> ln 1097-1115
    """
    x = np.linspace(0, 1, 100)
    y = stats.beta.pdf(x, a=mem_shape * mem_mean, b=mem_shape * (1 - mem_mean))
    plt.ylabel('Density')
    plt.xlabel('Memory (ratio)')
    plt.annotate('mem_strength: {0}\nmem_mean: {1}\nK: {2}'.format(mem_shape, mem_mean, thick), xy=(0.9, 0.9), xycoords='axes fraction',
                 horizontalalignment='right', verticalalignment='top')
    plt.plot(x, y)

def calib_plot(x, normalize=False):
    """Plot calibration curve
    
    I believe this plots the output from a Core instance's calibrate_dates() output.
    
        # calib.plot(info,
    #            date.res = date.res, rotate.axes = rotate.axes, width = width, cutoff = cutoff, rev.d = rev.d, rev.yr = rev.yr, normalise.dists = normalise.dists, C14.col = C14.col, C14.border = C14.border, cal.col = cal.col, cal.border = cal.border)
    # @ Bacon.R ln 1009 - 1052
    ### produce plots of the calibrated distributions    
    """
    width = 30
    fig, ax = plt.subplots()
    ax.set_ylim(4000, 7000)
    ax.set_xlim(0, 100)
    pat = []
    # d = x[0][i]
    for i, d in enumerate(x[0]):
        probs = x[1][i]
        z = np.array([probs[:, 0], width * probs[:, 1] / np.sum(probs[:, 1]) ])  # Normalize
        z = z[:, z[0].argsort(kind='mergesort')] # np.interp requires `xp` arg to be sorted.
        zy = np.linspace(np.min(z[0]), np.max(z[0]), num=200)
        zp = np.interp(x=zy, xp=z[0], fp=z[1])
        pol = np.vstack([np.concatenate([d + zp, d - zp[::-1]]),np.concatenate([zy, zy[::-1]])])
        pat.append(Polygon(pol.T, alpha = 0.75))
    p = PatchCollection(pat)
    ax.add_collection(p)
    ax.set_ylabel('cal yr BP')
    ax.set_xlabel('Depth')

def plot_acc_prior(acc_shape, acc_mean):
    """Plot prior distribution accumulation rate
    # PlotMemPrior @ Bacon.R ln 114 -> ln 1119 - 1141
    ## plot the prior for the memory (= accumulation rate varibility between neighbouring depths)
    """
    x = np.linspace(0, 3 * max([acc_mean]), 100)
    y = stats.gamma.pdf(x, a=acc_shape, scale =1 / (acc_shape / acc_mean))
    plt.ylabel('Density')
    plt.xlabel('Acc. rate (yr/cm)')
    plt.annotate('acc_shape: {0}\nacc_mean: {1}'.format(acc_shape, acc_mean), xy=(0.9, 0.9), xycoords='axes fraction',
                 horizontalalignment='right', verticalalignment='top')
    plt.plot(x, y)
