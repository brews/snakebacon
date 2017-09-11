import numpy as np
import scipy.stats as stats


def d_cal(calibcurve, rcmean, w2, cutoff=0.0001, normal_distr=False, t_a=3, t_b=4):
    """Get calendar date probabilities

    Parameters
    ----------
    calibcurve : CalibCurve
        Calibration curve.
    rcmean : scalar
        Reservoir-adjusted age.
    w2 : scalar
        r'$w^2_j(\theta)$' from pg 461 & 463 of Blaauw and Christen 2011.
    cutoff : scalar, optional
        Unknown.
    normal_distr : Bool, optional
        Use normal distribution for date errors. If False, then use Student's t-distribution.
    t_a : scalar, optional
        Student's t-distribution parameter, a. t_b - 1 must equal t_b.
    t_b : scalar, optional
        Student's t-distribution parameter, b. t_b - 1 must equal t_b.


    #Line 943 of Bacon.R
    #cc : calib_curve (3-col format)
    #rcmean : det['age'][i] - d_R
    #w2 : dat['error'][i]^2 + d_STD**2
    """
    assert t_b - 1 == t_a
    if normal_distr:
        # TODO(brews): Test this. Line 946 of Bacon.R.
        std = np.sqrt(calibcurve.error ** 2 + w2)
        dens = stats.norm(loc=rcmean, scale=std).pdf(calibcurve.c14age)
    else:
        # TODO(brews): Test this. Line 947 of Bacon.R.
        dens = (t_b + ((rcmean - calibcurve.c14age) ** 2) / (2 * (calibcurve.error ** 2 + w2))) ** (-1 * (t_a + 0.5))
    cal = np.array([calibcurve.calbp.copy(), dens]).T
    cal[:, 1] = cal[:, 1] / cal[:, 1].sum()
    # "ensure that also very precise dates get a range of probabilities"
    cutoff_mask = cal[:, 1] > cutoff
    if cutoff_mask.sum() > 5:
        out = cal[cutoff_mask, :]
    else:
        calx = np.linspace(cal[:, 0].min(), cal[:, 0].max(), num=50)
        caly = np.interp(calx, cal[:, 0], cal[:, 1])
        out = np.array([calx, caly / caly.sum()]).T
    return out


def calibrate_dates(chron, calib_curve, d_r, d_std, cutoff=0.0001, normal_distr=False, t_a=[3], t_b=[4]):
    """Get density of calendar dates for chron date segment in core

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
    t_a : scalar or ndarray, optional
        Student's t-distribution parameter, a. t_a - 1 must equal t_b.
    t_b : scalar or ndarray, optional
        Student's t-distribution parameter, b. t_a - 1 must equal t_b.

    Returns
    -------
    depth : ndarray
        Depth of dated sediment sample.
    probs : list of 2d arrays
        Density of calendar age for each dated sediment sample. For each
        sediment sample, the 2d array has two columns, the first is the
        calendar age. The second column is the density for that calendar age.

    """
    # Python version of .bacon.calib() on line 908 in Bacon.R

    # .bacon.calib - line 908

    # rcmean = 4128; w2 = 4225; t_a=3; t_b=4
    # test = d_cal(cc = calib_curve.rename(columns = {0:'a', 1:'b', 2:'c'}), rcmean = 4128, w2 = 4225, t_a=t_a,
    # t_b=t_b, cutoff=cutoff, normal = normal)

    # Line 959 of Bacon.R
    # calib = list(dets.iloc[:, 3])
    # Now Bacon goes and checks the ncol in the dets See line #960 in Bacon.R

    # Line #973
    # TODO(brews): Check that `normal_dist` is used and documented correctly in docstring above.
    # TODO(brews): Check whether we call returned values densities, freqs or what options we should have.
    n = len(chron.depth)
    calib_curve = np.array(calib_curve)
    t_a = np.array(t_a)
    t_b = np.array(t_b)
    assert t_b - 1 == t_a
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
    if len(calib_curve) == 1:
        calib_curve = np.repeat(calib_curve, n)

    calib_probs = []
    rcmean = chron.age - d_r
    w2 = chron.error ** 2 + d_std ** 2
    for i in range(n):
        age_realizations = d_cal(calib_curve[i], rcmean=rcmean[i], w2=w2[i],
                                 t_a=t_a[i], t_b=t_b[i],
                                 cutoff=cutoff, normal_distr=normal_distr)
        calib_probs.append(age_realizations)
    return np.array(chron.depth), calib_probs
