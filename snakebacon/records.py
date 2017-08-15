import logging

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon


log = logging.getLogger(__name__)


def read_14c(fl):
    """Create CalibCurve instance from Bacon curve file
    """
    indata = pd.read_csv(fl, index_col=None, skiprows=11, header=None,
                         names=['calbp', 'c14age', 'error', 'delta14c', 'sigma'])
    outcurve = CalibCurve(calbp=indata['calbp'],
                          c14age=indata['c14age'],
                          error=indata['error'],
                          delta14c=indata['delta14c'],
                          sigma=indata['sigma'])
    return outcurve


def read_dates(fl):
    """Create proxy instance from Bacon proxy file
    """
    indata = pd.read_table(fl, sep=r'\s*\,\s*', index_col=None, engine='python')
    outcore = DateRecord(age=indata['age'],
                         error=indata['error'],
                         depth=indata['depth'],
                         labid=indata['labID'])
    return outcore


def read_proxy(fl):
    """Read a file to create a proxy record instance
    """
    outcore = ProxyRecord(data=pd.read_table(fl, sep=r'\s*\,\s*', index_col=None, engine='python'))
    return outcore


class ProxyRecord:
    def __init__(self, data):
        """Create a proxy record instance

        Parameters
        ----------
        data : DataFrame
            Pandas dataframe containing columns with proxy sample measurements. Must also have 'depth' column.
        """
        assert 'depth' in data.columns.values
        assert 'age' not in data.columns.values
        self.data = pd.DataFrame(data).copy()

    def __repr__(self):
        return '%s(data=%r)' % (type(self).__name__, self.data)


class DatedProxyRecord(ProxyRecord):
    def __init__(self, data, age):
        """Create a dated proxy record instance

        Parameters
        ----------
        data : DataFrame
            Pandas dataframe containing columns with proxy sample measurements. Must also have 'depth' column.
        age : iterable
            Iterable containing calendar year, or a list of years (cal yr BP) for corresponding to each sample depth in
            data.depth.
        """
        super().__init__(data)
        assert len(data.depth) == len(age)
        self.age = np.array(age)

    def __repr__(self):
        return '%s(data=%r, age=%r)' % (type(self).__name__, self.data, self.age)

    def n_members(self):
        """Get number of MCMC ensemble members in calendar age estimates"""
        try:
            n = len(self.age[0])
        except TypeError:
            n = 1
        return n

    def to_pandas(self):
        """Convert record to pandas.DataFrame"""
        agedepthdf = pd.DataFrame(self.age, index=self.data.depth)
        agedepthdf.columns = list(range(self.n_members()))
        out = (agedepthdf.join(self.data.set_index('depth'))
               .reset_index()
               .melt(id_vars=self.data.columns.values, var_name='mciter', value_name='age'))
        out['mciter'] = pd.to_numeric(out.loc[:, 'mciter'])
        if self.n_members() == 1:
            out = out.drop('mciter', axis=1)
        return out


class DateRecord:
    def __init__(self, obj=None, **kwargs):
        """Create a sediment core date instance

        Parameters
        ----------
        obj : obj, optional
            Object with iterable attributes 'labid', 'age', 'error', and 'depth'. Assumes that depth is in increasing order
            order. Cannot use **kwargs if passing obj.
        **kwargs : optional
            Must include objects with iterables for 'labid', 'age', 'error', and 'depth'. Assumes depth is in in
            increasing order. Only parsed if obj is None.

        Returns
        -------
        A DateRecord instance.
        """
        if obj is not None:
            self.labid = np.array(obj.labid)
            self.age = np.array(obj.age)
            self.error = np.array(obj.error)  # Note this is called "std" in output .bacon file.
            self.depth = np.array(obj.depth)
        else:
            self.labid = np.array(kwargs['labid'])
            self.age = np.array(kwargs['age'])
            self.error = np.array(kwargs['error'])
            self.depth = np.array(kwargs['depth'])

    def __repr__(self):
        # return '%s(%r)' % (type(self).__name__, self)
        return '%s(age=%r, error=%r, depth=%r, labid=%r)' % (type(self).__name__, self.age, self.error, self.depth, self.labid)


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
        self.calbp = np.array(calbp)
        self.c14age = np.array(c14age)
        self.error = np.array(error)
        if delta14c is None:
            delta14c = np.zeros(calbp.shape)
        self.delta14c = np.array(delta14c)  # d_R
        if sigma is None:
            sigma = np.zeros(calbp.shape)
        self.sigma = np.array(sigma)  # d_R variance?

    def __repr__(self):
        return '%s(calbp=%r, c14age=%r, error=%r, delta14c=%r, sigma=%r)' % (type(self).__name__, self.calbp, self.c14age, self.error, self.delta14c, self.sigma)


def plot_accmem_prior(mem_shape, mem_mean, thick):
    """Plot accumulation rate varibility between neighbouring depths ("memory") prior

    # PlotAccPrior @ Bacon.R ln 113 -> ln 1097-1115
    """
    x = np.linspace(0, 1, 100)
    y = stats.beta.pdf(x, a=mem_shape * mem_mean, b=mem_shape * (1 - mem_mean))
    plt.ylabel('Density')
    plt.xlabel('Memory (ratio)')
    plt.annotate('mem_strength: {0}\nmem_mean: {1}\nK: {2}'.format(mem_shape, mem_mean, thick), xy=(0.9, 0.9),
                 xycoords='axes fraction',
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
        z = np.array([probs[:, 0], width * probs[:, 1] / np.sum(probs[:, 1])])  # Normalize
        z = z[:, z[0].argsort(kind='mergesort')]  # np.interp requires `xp` arg to be sorted.
        zy = np.linspace(np.min(z[0]), np.max(z[0]), num=200)
        zp = np.interp(x=zy, xp=z[0], fp=z[1])
        pol = np.vstack([np.concatenate([d + zp, d - zp[::-1]]), np.concatenate([zy, zy[::-1]])])
        pat.append(Polygon(pol.T, alpha=0.75))
    p = PatchCollection(pat)
    ax.add_collection(p)
    ax.set_ylabel('cal yr BP')
    ax.set_xlabel('Depth')


def plot_accrate_prior(acc_shape, acc_mean):
    """Plot prior distribution accumulation rate
    # PlotMemPrior @ Bacon.R ln 114 -> ln 1119 - 1141
    ## plot the prior for the memory (= accumulation rate varibility between neighbouring depths)
    """
    x = np.linspace(0, 3 * max([acc_mean]), 100)
    y = stats.gamma.pdf(x, a=acc_shape, scale=1 / (acc_shape / acc_mean))
    plt.ylabel('Density')
    plt.xlabel('Acc. rate (yr/cm)')
    plt.annotate('acc_shape: {0}\nacc_mean: {1}'.format(acc_shape, acc_mean), xy=(0.9, 0.9), xycoords='axes fraction',
                 horizontalalignment='right', verticalalignment='top')
    plt.plot(x, y)
