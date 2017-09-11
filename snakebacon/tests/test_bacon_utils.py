import unittest
import os

import numpy as np

from snakebacon.mcmcbackends.bacon import fetch_calibcurve
from snakebacon.mcmcbackends.bacon.utils import d_cal, calibrate_dates
from snakebacon import read_chron
import snakebacon.records as curve


here = os.path.abspath(os.path.dirname(__file__))


class TestBaconUtilsDcal(unittest.TestCase):
    def setUp(self):
        np.random.seed(123)
        self.minitestcurve = curve.CalibCurve(calbp=np.arange(1, 4),
                                              c14age=np.arange(3, 0, -1),
                                              error=np.arange(21, 24),
                                              delta14c=None, sigma=None)
        self.smalltestcurve = curve.CalibCurve(calbp=np.arange(1, 7),
                                               c14age=np.arange(6, 0, -1),
                                               error=np.arange(21, 27),
                                               delta14c=None, sigma=None)

    def test_d_cal_normal(self):
        goalmini = np.array([[2.959184, 0.01904568], [3.000000, 0.01900584]])
        goalsmall = np.array([[4, 0.1624059], [5, 0.1553523], [6, 0.1487062]])
        test_minivictim = d_cal(self.minitestcurve, rcmean=5, w2=2,
                                normal_distr=True)
        test_smallvictim = d_cal(self.smalltestcurve, rcmean=5, w2=2,
                                 normal_distr=True)
        np.testing.assert_allclose(test_minivictim[-2:], goalmini, atol=1e-7)
        np.testing.assert_allclose(test_smallvictim[-3:], goalsmall, atol=1e-7)

    def test_d_cal_notnormal(self):
        goalmini = np.array([[2.959184, 0.01990764], [3.000000, 0.01990352]])
        goalsmall = np.array([[4, 0.1667562], [5, 0.1662155], [6, 0.1655462]])
        test_minivictim = d_cal(self.minitestcurve, rcmean=5, w2=2,
                                cutoff=0.1, normal_distr=False)
        test_smallvictim = d_cal(self.smalltestcurve, rcmean=5, w2=2,
                                 cutoff=0.1, normal_distr=False)
        np.testing.assert_allclose(test_minivictim[-2:], goalmini, atol=1e-7)
        np.testing.assert_allclose(test_smallvictim[-3:], goalsmall, atol=1e-7)



class TestBaconUtils(unittest.TestCase):

    def test_calibrate_dates(self):
        pgoal_age_mean = 6726.8
        pgoal_age_std = 192.94564322687978

        pgoal_density_mean = 0.011642671907680679
        pgoal_density_std = 0.010893804880746285

        dgoal = np.array([77.5, 79.5, 99.5])
        dgoal_n = 40
        pgoal_n = dgoal_n

        chron = read_chron(os.path.join(here, 'MSB2K.csv'))

        d_r = [0]
        d_std = [0]
        t_a = [3]
        t_b = [4]

        cc = [fetch_calibcurve('IntCal13')]

        d_target, p_target = calibrate_dates(chron, calib_curve=cc,
                                             d_r=d_r, d_std=d_std,
                                             t_a=t_a, t_b=t_b)
        np.testing.assert_allclose(d_target[-3:], dgoal)
        np.testing.assert_equal(len(d_target), dgoal_n)

        np.testing.assert_allclose(np.mean(p_target[-1][:, 0]),
                                   pgoal_age_mean, atol=7)
        np.testing.assert_allclose(np.std(p_target[-1][:, 0]),
                                   pgoal_age_std, atol=2)
        np.testing.assert_equal(len(p_target), pgoal_n)


        np.testing.assert_allclose(np.mean(p_target[-1][:, 1]),
                                   pgoal_density_mean, atol=1e-2)
        np.testing.assert_allclose(np.std(p_target[-1][:, 1]),
                                   pgoal_density_std, atol=1e-2)


if __name__ == '__main__':
    unittest.main()
