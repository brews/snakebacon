import unittest
import os

import numpy as np

from snakebacon.mcmcbackends import Bacon
from snakebacon import read_dates


here = os.path.abspath(os.path.dirname(__file__))


class TestBaconUtils(unittest.TestCase):

    def test_calibrate_dates(self):
        pgoal_age_mean = 6702.5
        pgoal_age_std = 121.23496470353207

        pgoal_density_mean = 0.011642671907680679
        pgoal_density_std = 0.010893804880746285

        dgoal = np.array([77.5, 79.5, 99.5])
        dgoal_n = 40
        pgoal_n = dgoal_n

        chron = read_dates(os.path.join(here, 'MSB2K.csv'))
        mcmc_kwargs = dict(depth_min=1.5, depth_max=99.5, cc=[1],
                           cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4='ConstCal',
                           d_r=[0], d_std=[0], t_a=[3], t_b=[4], k=20,
                           minyr=-1000, maxyr=1e6, th01=4147, th02=4145,
                           acc_mean=20, acc_shape=1.5, mem_strength=4, mem_mean=0.7)
        d_target, p_target = Bacon.prior_dates(chron, **mcmc_kwargs)
        np.testing.assert_allclose(d_target[-3:], dgoal)
        np.testing.assert_equal(len(d_target), dgoal_n)

        np.testing.assert_allclose(np.mean(p_target[-1][:, 0]), pgoal_age_mean, atol=7)
        np.testing.assert_allclose(np.std(p_target[-1][:, 0]), pgoal_age_std, atol=2)
        np.testing.assert_equal(len(p_target), pgoal_n)


        np.testing.assert_allclose(np.mean(p_target[-1][:, 1]), pgoal_density_mean, atol=1e-2)
        np.testing.assert_allclose(np.std(p_target[-1][:, 1]), pgoal_density_std, atol=1e-2)


if __name__ == '__main__':
    unittest.main()
