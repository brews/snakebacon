import unittest
from copy import deepcopy
from os import path

import numpy as np

from snakebacon import read_chron
from snakebacon.mcmc import McmcResults, McmcSetup


here = path.abspath(path.dirname(__file__))


mcmc_kws = dict(depth_min=1.5, depth_max=99.5, cc=[1],
                   cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4='ConstCal',
                   d_r=[0], d_std=[0], t_a=[3], t_b=[4], k=20,
                   minyr=-1000, maxyr=1e6, th01=4147, th02=4145,
                   acc_mean=20, acc_shape=1.5, mem_strength=4, mem_mean=0.7)
fullrun_setup = McmcSetup(read_chron(path.join(here, 'MSB2K.csv')), **mcmc_kws)
fullrun_victim = McmcResults(fullrun_setup)


class TestMcmcResults(unittest.TestCase):
    def setUp(self):
        self.testdummy = deepcopy(fullrun_victim)

    def test_init(self):
        niter_goal = 3432
        nsegs_goal = 20
        objective_mean_goal = 218.26212744522144
        theta_mean_goal = 4549.6406701631704
        w_mean_goal = 0.06277018029195805
        x0_mean_goal = 14.487713009032634
        xneg1_mean_goal = 15.65397709810606
        # Fuzzy to deal with vars across platforms.
        np.testing.assert_allclose(len(self.testdummy.headage), niter_goal, atol=50)
        np.testing.assert_allclose(len(self.testdummy.sediment_memory), niter_goal, atol=50)
        np.testing.assert_allclose(len(self.testdummy.sediment_rate[0]), niter_goal, atol=50)
        np.testing.assert_allclose(len(self.testdummy.sediment_rate[-1]), niter_goal, atol=50)
        self.assertEqual(nsegs_goal, len(self.testdummy.sediment_rate))
        np.testing.assert_allclose(self.testdummy.objective.mean(), objective_mean_goal, atol=1e-1)
        np.testing.assert_allclose(self.testdummy.headage.mean(), theta_mean_goal, atol=15)
        np.testing.assert_allclose(self.testdummy.sediment_memory.mean(), w_mean_goal, atol=1e-2)
        np.testing.assert_allclose(self.testdummy.sediment_rate[0].mean(), x0_mean_goal, atol=2)
        np.testing.assert_allclose(self.testdummy.sediment_rate[-1].mean(), xneg1_mean_goal, atol=2)

    def test_n_members(self):
        n_goal = 3432
        n = self.testdummy.n_members()
        # Fuzzy to deal with vars across platforms.
        np.testing.assert_allclose(n, n_goal, atol=50)

    def test_burnin(self):
        n_goal = 3432 - 200
        self.testdummy.burnin(n=200)
        # Fuzzy to deal with vars across platforms.
        np.testing.assert_allclose(len(self.testdummy.headage), n_goal, atol=50)
        np.testing.assert_allclose(len(self.testdummy.sediment_memory), n_goal, atol=50)
        np.testing.assert_allclose(len(self.testdummy.sediment_rate[0]), n_goal, atol=50)


if __name__ == '__main__':
    unittest.main()
