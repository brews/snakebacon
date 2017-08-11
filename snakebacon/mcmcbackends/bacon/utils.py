import numpy as np


def calibrate_dates(chron, calib_curve, d_r, d_std, cutoff=0.001, normal_distr=False, t_a=3, t_b=4):
    """Get probability of calendar dates for each depth segment in core

    Parameters
    ----------
    chron : DatedProxy-like
    calib_curve : CalibCurve or list of CalibCurves
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
    depth : ndarray
        Depth of Fixed-length sediment segment.
    probs : list of ndarrays
        Probabilities with one ndarray per fixed-length core segment.

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
    # TODO(brews): Check that `normal_dist` is used and documented correctly in docstring above.
    assert t_b - 1 == t_a
    n = len(chron)
    t_a = np.array(t_a)
    t_b = np.array(t_b)
    d_r = np.array(d_r)
    d_std = np.array(d_std)
    if len(t_a) == 1:
        t_a = np.repeat(t_a, n)
    if len(t_b) == 1:
        t_b = np.repeat(t_b, n)
    if len(d_r) == 1:
        d_r = np.repeat(d_r, n)
    if len(d_std) == 1:
        d_std = np.repeat(d_std, n)

    calib_probs = []
    # TODO(brews): I think we can do the below without a loop.
    rcmean = chron.age - d_r
    w2 = chron.error ** 2 + d_std ** 2
    try:
        for i in range(n):
            age_realizations = calib_curve[i].d_cal(rcmean=rcmean[i], w2=w2[i],
                                                    t_a=t_a[i], t_b=t_b[i],
                                                    cutoff=cutoff,
                                                    normal_distr=normal_distr)
            calib_probs.append(age_realizations)
    except TypeError:
        calib_probs = calib_curve.d_cal(rcmean=rcmean, w2=w2, t_a=t_a, t_b=t_b,
                                        cutoff=cutoff,
                                        normal_distr=normal_distr)
    return np.array(chron.depth), calib_probs
