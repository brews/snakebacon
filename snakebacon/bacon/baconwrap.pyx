import pandas as pd
from libc.stdlib cimport malloc, free


def read_baconout(path):
    """Read output from _baconmain"""
    d = pd.read_table(path, delim_whitespace=True, header=None)
    # TODO(brews): Not sure about the outgoing structure here. Might transpose depending on which is easier in later analysis.
    out = {'theta': d.iloc[:, 0].values,  # `theta0` or often just `theta`, array (i) of Age of sedimentation core head.
           'x': d.iloc[:, 1:-2].values,  # `x`, 2d array (i, j) of sediment accumulation rates for each segment (j) down the sediment core of each MCMC iteration (i).
           'w': d.iloc[:, -2].values,  # `w`, array (i) of memory or coherence of accumulation rates along sediment core.
           'objective': d.iloc[:, -1].values}  # `Us`, array (i) of objective or energy function used in the twalk MCMC.
    return out


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