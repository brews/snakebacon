import os
import sys
import datetime
import shutil
import contextlib
import tempfile
import numpy as np
import pandas as pd
from libc.stdlib cimport malloc, free
import inspect


if not hasattr(sys.modules[__name__], '__file__'):
    __file__ = inspect.getfile(inspect.currentframe())


HERE = os.path.abspath(os.path.dirname(__file__))
CURVEDIR_PATH = os.path.join(HERE, 'Curves')


@contextlib.contextmanager
def try_chdir(path):
    """Use in `with as`, returning to CWD afterwards"""
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


def run_bacon(inpath, outpath, ssize):
    """TODO: Something of a test function for now."""
    cwd = os.getcwd()
    inpath_fl = os.path.basename(inpath)
    outpath_fl = os.path.basename(outpath)
    curves_dir = os.path.basename(CURVEDIR_PATH)
    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copytree(CURVEDIR_PATH, os.path.join(tmpdir, curves_dir))
        shutil.copy2(inpath, tmpdir)
        with try_chdir(tmpdir):
            _baconmain(inpath_fl, outpath_fl, ssize)
        shutil.copy2(os.path.join(tmpdir, outpath_fl), outpath)
    print('done')


def read_baconout(path):
    """Read output from _baconmain"""
    d = pd.read_table(path, delim_whitespace=True, header=None)
    # TODO(brews): Function cannot handle hiatus
    # TODO(brews): Not sure about the outgoing structure here. Might transpose depending on which is easier in later analysis.
    out = {'theta': d.iloc[:, 0].values,  # `theta0` or often just `theta`, array (i) of Age of sedimentation core head.
           'x': d.iloc[:, 1:-2].values.T,  # `x`, 2d array (i, j) of sediment accumulation rates for each segment (i) down the sediment core of each MCMC iteration (j).
           'w': d.iloc[:, -2].values,  # `w`, array (i) of memory or coherence of accumulation rates along sediment core.
           'objective': d.iloc[:, -1].values}  # `Us`, array (i) of objective or energy function used in the twalk MCMC.
    return out


def _baconin_str(*, core_labid, core_age, core_error, core_depth, depth_min, depth_max, cc, cc1, cc2, cc3, cc4, d_r,
                 d_std, t_a, t_b, k, th01, th02, mem_strength, mem_mean, acc_shape, acc_mean, minyr=-1000, maxyr=1e6,
                  normal=False, postbomb=0):
    """Get string to write to .bacon file
    See .write.Bacon.file @ Bacon.R ln 385-480
    """
    # import snakebacon as snek
    # c = snek.read_corefile('/home/sbm/Desktop/dtda sandbox/bacon/LinBacon_2.2/Cores/MSB2K/MSB2K.csv')
    # Test with:
    # write_baconin('testout.txt', core_labid = c.labid, core_age = c.age, core_error = c.error, core_depth = c.depth, depth_min = 1.5, depth_max = 99.5, cc=[1], cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4 = 'ConstCal', d_r = [0], d_std = [0], t_a=[3], t_b=[4], k= 20, minyr=-1000, maxyr = 1e6, th01=4147, th02=4145, acc_mean = 20, acc_shape = 1.5, mem_strength = 4, mem_mean = 0.7)
    # TODO(brews): Function cannot handle hiatus.
    if depth_min < min(core_depth):
        extrap = [np.nan, max(core_age), np.max([1e4, np.max(100 * core_error)]), depth_max, 0]
        dets = np.array([np.array([core_labid, core_age, core_error, core_depth]).T, extrap])
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
    # TODO(brews): Cannot do multiple callibration curves via dets[,5]. See Bacon.R @ ln 424:449
    if len(d_r) == 1:
        d_r = np.repeat(d_r, len(core_labid))
    if len(d_std) == 1:
        d_std = np.repeat(d_std, len(core_labid))
    if len(t_a) == 1:
        t_a = np.repeat(t_a, len(core_labid))
    if len(t_b) == 1:
        t_b = np.repeat(t_b, len(core_labid))
    if len(cc) == 1:
        cc = np.repeat(cc, len(core_labid))
    for i, ln_labid in enumerate(core_labid):
        # ln_cc = # Bacon.R @ ln 448.
        outlines.append(str_template.format(count=i, id=ln_labid,
                                            age=core_age[i], error=core_error[i], depth=core_depth[i], r=d_r[i],
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
    """Write .bacon file to be read into baconmain
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
    cdef extern from "bacon.cpp":
        int notmain(int argc, char *argv[])
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