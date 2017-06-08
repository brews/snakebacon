import unittest
from copy import deepcopy
from os import path

import numpy as np
import pandas as pd

from snakebacon import read_dates, ProxyRecord
from snakebacon.agedepth import AgeDepthModel

here = path.abspath(path.dirname(__file__))

fullrun_agemodel = AgeDepthModel(read_dates(path.join(here, 'MSB2K.csv')),
                                 depth_min=1.5, depth_max=99.5, cc=[1],
                                 cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4='ConstCal',
                                 d_r=[0], d_std=[0], t_a=[3], t_b=[4], k=20,
                                 minyr=-1000, maxyr=1e6, th01=4147, th02=4145,
                                 acc_mean=20, acc_shape=1.5, mem_strength=4, mem_mean=0.7)


class TestAgeDepth(unittest.TestCase):
    def setUp(self):
        self.testdummy = deepcopy(fullrun_agemodel)

    def test_init(self):
        goal_thick = 4.9
        goal_depthrange = (1.5, 99.5)
        goal_agemedian_edges = (4552.5750000000007, 6714.1464476499987)
        goal_agemedian_025 = (4430.3599999999997, 6582.4278499575012)
        goal_agemedian_975 = (4680.6210000000001, 6868.9317877499998)
        goal_ageensemble_shape = (99, 3232)

        self.assertTrue(self.testdummy.mcmcfit is not None)
        self.assertEqual(goal_thick, self.testdummy.thick)
        self.assertTupleEqual(goal_depthrange, (min(self.testdummy.depth), max(self.testdummy.depth)))

        np.testing.assert_allclose(self.testdummy.age_median[0], goal_agemedian_edges[0], atol=1e-2)
        np.testing.assert_allclose(self.testdummy.age_median[-1], goal_agemedian_edges[-1], atol=1e-2)

        np.testing.assert_allclose(self.testdummy.conf_interv[2.5][0], goal_agemedian_025[0], atol=1e-2)
        np.testing.assert_allclose(self.testdummy.conf_interv[2.5][-1], goal_agemedian_025[-1], atol=1e-2)

        np.testing.assert_allclose(self.testdummy.conf_interv[97.5][0], goal_agemedian_975[0], atol=1e-2)
        np.testing.assert_allclose(self.testdummy.conf_interv[97.5][-1], goal_agemedian_975[-1], atol=1e-2)

        self.assertTupleEqual(goal_ageensemble_shape, (len(self.testdummy.age_ensemble),
                                                       len(self.testdummy.age_ensemble[0])))

    def test_date(self):
        np.random.seed(123)
        goal_median_idx0 = 4552.5750000000007
        goal_median_nmember = 3
        goal_ens_2_idx0 = 4503.01
        goal_ens_2_nmember = 2
        goal_ens_20_nmember = 20
        testproxy = ProxyRecord(pd.DataFrame({'depth': np.arange(1.5, 4.5), 'a': np.arange(20, 23)}))
        victim_median = self.testdummy.date(testproxy, how='median')
        victim_ens_2 = self.testdummy.date(testproxy, how='ensemble', n=2)
        victim_ens_20 = self.testdummy.date(testproxy, how='ensemble', n=20)

        np.testing.assert_allclose(victim_median.age[0], goal_median_idx0, atol=1)
        self.assertEqual(goal_median_nmember, len(victim_median.age))

        np.testing.assert_allclose(victim_ens_2.age[0][0], goal_ens_2_idx0, atol=1)
        self.assertEqual(goal_ens_2_nmember, len(victim_ens_2.age[0]))

        self.assertEqual(goal_ens_20_nmember, len(victim_ens_20.age[0]))

    def test_agedepth(self):
        goal_len = 3232
        goal_mean = 4567.1221008870671
        goal_var = 3982.7408290300741
        victim = self.testdummy.agedepth(2.5)
        self.assertEqual(goal_len, len(victim))
        np.testing.assert_allclose(victim.mean(), goal_mean, atol=1e-2)
        np.testing.assert_allclose(victim.var(), goal_var, atol=1e-2)


if __name__ == '__main__':
    unittest.main()
