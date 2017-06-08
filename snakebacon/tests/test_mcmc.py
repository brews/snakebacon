import unittest
from copy import deepcopy
from os import path

import numpy as np

from snakebacon import read_dates
from snakebacon.mcmc import McmcResults

here = path.abspath(path.dirname(__file__))

fullrun_victim = McmcResults(read_dates(path.join(here, 'MSB2K.csv')),
                             depth_min=1.5, depth_max=99.5, cc=[1],
                             cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4='ConstCal',
                             d_r=[0], d_std=[0], t_a=[3], t_b=[4], k=20,
                             minyr=-1000, maxyr=1e6, th01=4147, th02=4145,
                             acc_mean=20, acc_shape=1.5, mem_strength=4, mem_mean=0.7)


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
        self.assertEqual(niter_goal, len(self.testdummy.headage))
        self.assertEqual(niter_goal, len(self.testdummy.sediment_memory))
        self.assertEqual(niter_goal, len(self.testdummy.sediment_rate[0]))
        self.assertEqual(niter_goal, len(self.testdummy.sediment_rate[-1]))
        self.assertEqual(nsegs_goal, len(self.testdummy.sediment_rate))
        np.testing.assert_allclose(self.testdummy.objective.mean(), objective_mean_goal, atol=1e-2)
        np.testing.assert_allclose(self.testdummy.headage.mean(), theta_mean_goal, atol=1)
        np.testing.assert_allclose(self.testdummy.sediment_memory.mean(), w_mean_goal, atol=1e-2)
        np.testing.assert_allclose(self.testdummy.sediment_rate[0].mean(), x0_mean_goal, atol=1e-2)
        np.testing.assert_allclose(self.testdummy.sediment_rate[-1].mean(), xneg1_mean_goal, atol=1e-2)

    def test_n_members(self):
        n_goal = 3432
        n = self.testdummy.n_members()
        self.assertEqual(n_goal, n)

    def test_burnin(self):
        n_goal = 3432 - 200
        self.testdummy.burnin(n=200)
        self.assertEqual(n_goal, len(self.testdummy.headage))
        self.assertEqual(n_goal, len(self.testdummy.sediment_memory))
        self.assertEqual(n_goal, len(self.testdummy.sediment_rate[0]))


if __name__ == '__main__':
    unittest.main()
