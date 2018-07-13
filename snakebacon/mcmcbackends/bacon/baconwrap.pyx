import os
import datetime
import shutil
import contextlib
import tempfile
import numpy as np
import pandas as pd
from bacon cimport notmain
from libc.stdlib cimport malloc, free
from .Curves import here as curvespath


@contextlib.contextmanager
def try_chdir(path):
    """Use in `with as`, returning to CWD afterwards"""
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


def run_baconmcmc(ssize=2000, **kwargs):
    """Run bacon MCMC, given parameters.

    Parameters
    ----------
        ssize : int
        ???  # TODO(brews): I have no idea what this actually does.
        **kwargs :
        Bacon MCMC run parameters passed to `write_baconin()`.

    Returns
    -------
    Output from bacon MCMC.
    """
    cwd = os.getcwd()
    infile_str = 'intobacon.txt'
    outfile_str = 'outofbacon.bacon'
    curves_dir = os.path.basename(curvespath)
    with tempfile.TemporaryDirectory() as tmpdir:
        write_baconin(os.path.join(tmpdir, infile_str), **kwargs)
        shutil.copytree(curvespath, os.path.join(tmpdir, curves_dir))
        with try_chdir(tmpdir):
            _baconmain(infile_str, outfile_str, ssize)
            out = read_baconout(outfile_str)
    return out


def run_baconmcmcfiles(inpath, outpath, ssize=2000):
    """Run bacon MCMC using bacon file at inpath and write results to outpath"""
    cwd = os.getcwd()
    inpath_fl = os.path.basename(inpath)
    outpath_fl = os.path.basename(outpath)
    curves_dir = os.path.basename(curvespath)
    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(curvespath, os.path.join(tmpdir, curves_dir))
        shutil.copy2(inpath, tmpdir)
        with try_chdir(tmpdir):
            _baconmain(inpath_fl, outpath_fl, ssize)
        shutil.copy2(os.path.join(tmpdir, outpath_fl), outpath)
    print('Done.')


def read_baconout(path):
    """Read output file from bacon MCMC

    Parameters
    ----------
    path : str
        Path of file output from bacon MCMC.

    Returns
    -------
    Dictionary with four members. 'theta' (or 'theta0') array (i) of calendar age of sedimentation core head for i MCMC
    iterations retained. 'x', a 2d array (j, i) of sediment accumulation rates for each fixed-length segment (j) down
    the sediment core of each MCMC iteration member (i). 'w', array (i) of memory or coherence of accumulation rates
    along sediment core for each MCMC iteration member. 'objective' (i.e. 'Us'), an array (i) of objective or energy
    function values used in the twalk MCMC.
    """
    d = pd.read_table(path, delim_whitespace=True, header=None)
    # TODO(brews): Function cannot handle hiatus
    out = {'theta': d.iloc[:, 0].values,  # `theta0` or often just `theta`, array (i) of Age of sedimentation core head.
           'x': d.iloc[:, 1:-2].values.T,  # `x`, 2d array (i, j) of sediment accumulation rates for each segment (i) down the sediment core of each MCMC iteration (j).
           'w': d.iloc[:, -2].values,  # `w`, array (i) of memory or coherence of accumulation rates along sediment core.
           'objective': d.iloc[:, -1].values}  # `Us`, array (i) of objective or energy function used in the twalk MCMC.
    return out


def _baconin_str(*, core_labid, core_age, core_error, core_depth, depth_min, depth_max, cc, d_r, d_std, k, th01, th02,
                 mem_strength, mem_mean, acc_shape, acc_mean, t_a=None, t_b=None, cc1='IntCal13', cc2='Marine13',
                 cc3='SHCal13', cc4='ConstCal', minyr=-1000, maxyr=1e6, normal=False, postbomb=0):
    """Make list of strings to write to .bacon file

    Parameters
    ----------
    core_labid : n-length iterable of strs
        Laboratory ID for each core sample.
    core_age : n-length iterable of floats
        Carbon-14 age of each core sample.
    core_error : n-length iterable of floats
        Carbon-14 age error (standard deviation) for each core sample.  # TODO(brews): Clarify what this means.
    core_depth : n-length iterable of floats
        Depth (cm) at which each core sample was taken.
    depth_min : float
        Minimum core depth (cm) to be considered. If outside the range of core_depth, will be extrapolated.
    depth_max : float
        Maximum core depth (cm) to be considered. If outside the range of core_depth, will be extrapolated.
    cc : list of ints
        Int indicating which calibration curve to use for each core sample (i.e. 'cc1', 'cc2', ... 'cc4'). If list
        contains single value, this value is used for all core samples.
    cc1 : str, optional
        String indicating calibration curve option. Must be one of 'IntCal13', 'Marine13', 'SHCal13', 'ConstCal'.
        Default is 'IntCal13'.
    cc2 : str, optional
        String indicating calibration curve option. Must be one of 'IntCal13', 'Marine13', 'SHCal13', 'ConstCal'.
        Default is 'Marine13'
    cc3 : str, optional
        String indicating calibration curve option. Must be one of 'IntCal13', 'Marine13', 'SHCal13', 'ConstCal'.
        Default is 'SHCal13'.
    cc4 : str, optional
        String indicating calibration curve option. Must be one of 'IntCal13', 'Marine13', 'SHCal13', 'ConstCal'.
        Default is 'ConstCal'.
    d_r : n-length iterable of floats
        Delta carbon reservoir (delta R) values for each core sample. If single value is given, it will be used for all
        core samples.
    d_std : n-length iterable of floats
        Delta carbon reservoir (delta R) standard deviation for each core sample. If single value is given, it will be
        used for all core samples.
    t_a : n-length iterable of floats
        Parameter 'a' used in t-walk MCMC. If single value is given, it will be used for all core samples. Must be one
        less than 't_b'. Default is 3.
    t_b : n-length iterable of floats
        Parameter 'b' used in t-walk MCMC. If single value is given, it will be used for all core samples. Must be one
        greater than 't_a'. Default is 4.
    k : int
        Number of fixed-length segments to divide core into for MCMC.
    minyr : int, optional
        Lowest year considered in MCMC.  # TODO(brews): See line 368 of Bacon.R. Note how this is calculated.
    maxyr : int, optional
        Highest year considered in MCMC.  # TODO(brews): See line 370 of Bacon.R. Note how this is calculated.
    th01 : int
        Initial ballpark guess for top-most fixed segment age. # TODO(brews): Check on this @ ln 474 of Bacon.R
    th02 : int
        Another Initial ballpark guess for fixed segment age. # TODO(sbm): Check on this @ ln 474 of Bacon.R
    acc_mean : float
        Mean sediment accumulation rate (cm/yr).
    acc_shape : float
        Sediment accumulation rate distribution shape.
    mem_strength: float
        Mean sediment accumulation rate strength from one fixed-length segment to the next.  # TODO(brews): Find more details on this.
    mem_mean: float
        Mean sediment accumulation rate memory from one fixed-length segment to the next. Must be between 1 and 0.
    normal: bool, optional
        Whether to use a normal for radiocarbon age. Default is 'False', using Student's t based distribution.
    postbomb: int, optional
        Integer indicating which post-bomb curve to use for radiocarbon calibration. `0` is no curve. `1` is
        'postbomb_NH1.14C'. `2` is 'postbomb_NH2.14C'. `3` is 'postbomb_NH3.14C'. `4` is 'postbomb_SH1-2.14C'.
        `5` is 'postbomb_SH3.14C'.

    Returns
    -------
    List of strings for each line of a text file to be read into the bacon MCMC.

    See .write.Bacon.file @ Bacon.R ln 385-480
    """
    # import snakebacon as snek
    # c = snek.read_corefile('/home/sbm/Desktop/dtda sandbox/bacon/LinBacon_2.2/Cores/MSB2K/MSB2K.csv')
    # Test with:
    # write_baconin('testout.txt', core_labid = c.labid, core_age = c.age, core_error = c.error, core_depth = c.depth, depth_min = 1.5, depth_max = 99.5, cc=[1], cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4 = 'ConstCal', d_r = [0], d_std = [0], t_a=[3], t_b=[4], k= 20, minyr=-1000, maxyr = 1e6, th01=4147, th02=4145, acc_mean = 20, acc_shape = 1.5, mem_strength = 4, mem_mean = 0.7)
    # TODO(brews): Function cannot handle hiatus.
    core_labid = np.array(core_labid, dtype='<U16')
    core_depth = np.array(core_depth, dtype='float64')
    if len(cc) == 1:
        cc = np.repeat(cc, len(core_depth))
    if len(d_r) == 1:
        d_r = np.repeat(d_r, len(core_depth))
    if len(d_std) == 1:
        d_std = np.repeat(d_std, len(core_depth))
    if depth_min < min(core_depth):
        core_labid = np.insert(np.array(core_labid), 0, 'NA')
        core_age = np.insert(np.array(core_age), 0, [np.min(core_age)])
        core_error = np.insert(np.array(core_error), 0, np.max([1e4, np.max(100 * core_error)]))
        core_depth = np.insert(np.array(core_depth), 0, [depth_min])
        cc = np.insert(np.array(cc), 0, [0])
        d_r = np.insert(np.array(d_r), 0, [0])
        d_std = np.insert(np.array(d_std), 0, [0])
    if depth_max > max(core_depth):
        core_labid = np.append(np.array(core_labid), 'NA')
        core_age = np.append(np.array(core_age), [np.max(core_age)])
        core_error = np.append(np.array(core_error), np.max([1e4, np.max(100 * core_error)]))
        core_depth = np.append(np.array(core_depth), [depth_max])
        cc = np.append(np.array(cc), [0])
        d_r = np.append(np.array(d_r), [0])
        d_std = np.append(np.array(d_std), [0])

    # Non-14C values do not have C-reservoir adjustments.
    nonc14mask = cc == 0
    d_r[nonc14mask] = 0
    d_std[nonc14mask] = 0

    if t_a is None:
        t_a = [3]
    if t_b is None:
        t_b = [4]
    if len(t_a) == 1:
        t_a = np.repeat(t_a, len(core_depth))
    if len(t_b) == 1:
        t_b = np.repeat(t_b, len(core_depth))

    outlines = list()
    outlines.append('## Ran on {0}\n\n'.format(datetime.datetime.today().strftime('%c')))
    outlines.append('Cal 0 : ConstCal;\n')
    outlines.append('Cal 1 : {0}, {1};\n'.format(cc1, postbomb))
    outlines.append('Cal 2 : {0};\n'.format(cc2))  # TODO(brews): This is very weird @ ln 400:406.
    outlines.append('Cal 3 : {0}, {1};\n'.format(cc3, postbomb))
    # TODO(brews): Find out if this last bit is a way to specify a calibration curve directory. ln 409
    # outlines.append('Cal 4 : {0}, {1};'.format(cc3, postbomb))
    # Something with os.path.sep on windows, see Bacon.R 406:409
    outlines.append('\n##   id.   yr    std   depth  d.R  d.STD     t.a   t.b   cc\n')
    str_template = 'Det {count} : {id} , {age}, {error}, {depth}, {r}, {std}, {a}, {b}, {cc};\n'
    for i, ln_depth in enumerate(core_depth):
        # ln_cc = # Bacon.R @ ln 448.
        outlines.append(str_template.format(count=i, id=core_labid[i],
                                            age=core_age[i], error=core_error[i], depth=ln_depth, r=d_r[i],
                                            std=d_std[i], a=t_a[i], b=t_b[i], cc=cc[i]))
    # TODO(brews): if for hiatus @ Bacon.R ln 451:469
    wrapup_header = '\n##\t\t K   MinYr   MaxYr   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax\n'
    outlines.append(wrapup_header)
    dist_opt = 'FixT'
    if normal:
        dist_opt = 'FixNor'
    wrapup_body = 'Bacon 0: {dist_opt}, {k}, {minyr}, {maxyr}, {th0}, {th0p}, {w_a}, {w_b}, {alpha}, {beta}, {dmin}, {dmax};\n'
    outlines.append(wrapup_body.format(dist_opt=dist_opt, k=k, minyr=minyr, maxyr=maxyr,
                                       th0=th01, th0p=th02,
                                       w_a=mem_strength * mem_mean, w_b=mem_strength * (1 - mem_mean),
                                       alpha=acc_shape, beta=acc_shape / acc_mean,
                                       dmin=depth_min, dmax=depth_max))
    return outlines


def write_baconin(path, **kwargs):
    """Write .bacon file to be read into bacon MCMC
    """
    outlines = _baconin_str(**kwargs)
    with open(path, 'w') as fl:
        fl.writelines(outlines)


def _baconmain(str infile, str outfile, int ssize):
    """Run bacon MCMC on input file, and put output into outfile
    
    The underlying C/C++ from Bacon assumes there is a directory, 'Curve' in runtime CWD that holds special format 
    calibration curve. 
    
    Parameters
    ----------
        infile : str
            Path of existing bacon-format file to be input into MCMC.
        outfile : str
            Path of file where bacon MCMC results are dumped.
        ssize : int
            Sample size of input data.
            
    Returns
    -------
    int relating to burn-in and sub-sample thinning parameters. See bacon.cpp lines 41-44 and 156.
    """
    cdef char **outgoing_argv
    cdef bytes binfile
    cdef bytes boutfile
    cdef bytes bssize
    cdef char* cinfile
    cdef char* coutfile
    cdef char* cssize
    cdef int argcount = 4
    outgoing_argv = <char**> malloc(argcount * sizeof(char*))
    try:
        binfile = infile.encode('utf8')
        cinfile = binfile

        boutfile = outfile.encode('utf8')
        coutfile = boutfile

        bssize = str(ssize).encode('utf8')
        cssize = bssize
        
        outgoing_argv[1] = cinfile
        outgoing_argv[2] = coutfile
        outgoing_argv[3] = cssize
        return notmain(argcount, outgoing_argv)
    finally:
        free(outgoing_argv)


# _baconmain('MSB2K_20.bacon', 'out.bacon', 2000)