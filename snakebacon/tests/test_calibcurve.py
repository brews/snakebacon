import unittest

import numpy as np

import snakebacon.records as curve


class TestCalibCurveMethods(unittest.TestCase):
    def setUp(self):
        np.random.seed(123)
        self.minitestcurve = curve.CalibCurve(calbp=np.arange(1, 4), c14age=np.arange(3, 0, -1),
                                              error=np.arange(21, 24), delta14c=None, sigma=None)
        self.smalltestcurve = curve.CalibCurve(calbp=np.arange(1, 7), c14age=np.arange(6, 0, -1),
                                               error=np.arange(21, 27), delta14c=None, sigma=None)

    def test_d_cal_normal(self):
        goalmini = np.array([[2.959184, 0.01904568], [3.000000, 0.01900584]])
        goalsmall = np.array([[4, 0.1624059], [5, 0.1553523], [6, 0.1487062]])
        test_minivictim = self.minitestcurve.d_cal(rcmean=5, w2=2, normal_distr=True)
        test_smallvictim = self.smalltestcurve.d_cal(rcmean=5, w2=2, normal_distr=True)
        np.testing.assert_allclose(test_minivictim[-2:], goalmini, atol=1e-7)
        np.testing.assert_allclose(test_smallvictim[-3:], goalsmall, atol=1e-7)

    def test_d_cal_notnormal(self):
        goalmini = np.array([[2.959184, 0.01990764], [3.000000, 0.01990352]])
        goalsmall = np.array([[4, 0.1667562], [5, 0.1662155], [6, 0.1655462]])
        test_minivictim = self.minitestcurve.d_cal(rcmean=5, w2=2, cutoff=0.1, normal_distr=False)
        test_smallvictim = self.smalltestcurve.d_cal(rcmean=5, w2=2, cutoff=0.1, normal_distr=False)
        np.testing.assert_allclose(test_minivictim[-2:], goalmini, atol=1e-7)
        np.testing.assert_allclose(test_smallvictim[-3:], goalsmall, atol=1e-7)


class TestRead14c(unittest.TestCase):
    @unittest.skip('Test not written')
    def test_read_14c(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
