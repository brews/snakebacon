import unittest
import os

import numpy as np

from snakebacon.mcmcbackends import Bacon
from snakebacon import read_chron


here = os.path.abspath(os.path.dirname(__file__))


class TestBaconMethods(unittest.TestCase):

    def setUp(self):
        self.mcmc_kws = dict(depth_min=1.5, depth_max=99.5, cc=[1],
                                cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4='ConstCal',
                                d_r=[0], d_std=[0], t_a=[3], t_b=[4], k=20,
                                minyr=-1000, maxyr=1e6, th01=4147, th02=4145,
                                acc_mean=20, acc_shape=1.5, mem_strength=4, mem_mean=0.7)

    def test_prior_dates(self):
        pgoal_age_mean = 6726.8
        pgoal_age_std = 192.94564322687978

        pgoal_density_mean = 0.011642671907680679
        pgoal_density_std = 0.010893804880746285

        dgoal = np.array([77.5, 79.5, 99.5])
        dgoal_n = 40
        pgoal_n = dgoal_n

        chron = read_chron(os.path.join(here, 'MSB2K.csv'))
        d_target, p_target = Bacon.prior_dates(chron, **self.mcmc_kws)
        np.testing.assert_allclose(d_target[-3:], dgoal)
        np.testing.assert_equal(len(d_target), dgoal_n)

        np.testing.assert_allclose(np.mean(p_target[-1][:, 0]), pgoal_age_mean, atol=7)
        np.testing.assert_allclose(np.std(p_target[-1][:, 0]), pgoal_age_std, atol=2)
        np.testing.assert_equal(len(p_target), pgoal_n)

        np.testing.assert_allclose(np.mean(p_target[-1][:, 1]), pgoal_density_mean, atol=1e-2)
        np.testing.assert_allclose(np.std(p_target[-1][:, 1]), pgoal_density_std, atol=1e-2)

    def test_prior_sediment_rate(self):
        np.random.seed(123)

        goal_mean = 0.008194080517375957
        goal_std = 0.01172658825323754
        goal_n = 100

        victim, x = Bacon.prior_sediment_rate(**self.mcmc_kws)

        np.testing.assert_equal(len(victim), goal_n)
        np.testing.assert_allclose(victim.mean(), goal_mean, atol=1e-3)
        np.testing.assert_allclose(victim.std(), goal_std, atol=1e-3)

    def test_prior_sediment_memory(self):
        np.random.seed(123)

        goal_mean = 0.98457848264590286
        goal_std = 0.71613816177236256
        goal_n = 100

        victim, x = Bacon.prior_sediment_memory(**self.mcmc_kws)

        np.testing.assert_equal(len(victim), goal_n)
        np.testing.assert_allclose(victim.mean(), goal_mean, atol=1e-3)
        np.testing.assert_allclose(victim.std(), goal_std, atol=1e-3)


if __name__ == '__main__':
    unittest.main()
