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


def read_chron(fl):
    """Create ChronRecord instance from Bacon file
    """
    indata = pd.read_table(fl, sep=r'\s*\,\s*', index_col=None, engine='python')
    outcore = ChronRecord(age=indata['age'],
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


class ChronRecord:
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
        A ChronRecord instance.
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
