import unittest
from os import path

import numpy as np

import snakebacon as snek
from snakebacon.mcmcbackends.bacon import baconwrap

here = path.abspath(path.dirname(__file__))


class TestBaconwrap(unittest.TestCase):
    def test__baconin_str(self):
        # Not very clever test, but okay first pass.
        victim = baconwrap._baconin_str(core_labid=np.array(['a', 'b', 'c']), core_age=np.array([1, 2, 3]),
                                        core_error=np.array([1, 2, 1]), core_depth=np.array([2, 3, 4]), depth_min=1.5,
                                        depth_max=99.5, cc=[1], cc1='IntCal13', cc2='Marine13', cc3='SHCal13',
                                        cc4='ConstCal', d_r=[0], d_std=[0], t_a=[3], t_b=[4], k=20,
                                        minyr=-1000, maxyr=1e6, th01=4147, th02=4145, acc_mean=20, acc_shape=1.5,
                                        mem_strength=4, mem_mean=0.7)
        goal = ['## Ran on Tue 09 May 2017 02:57:05 PM \n\n',
                'Cal 0 : ConstCal;\n',
                'Cal 1 : IntCal13, 0;\n',
                'Cal 2 : Marine13;\n',
                'Cal 3 : SHCal13, 0;\n',
                '\n##   id.   yr    std   depth  d.R  d.STD     t.a   t.b   cc\n',
                'Det 0 : NA , 1, 10000.0, 1.5, 0, 0, 3, 4, 0;\n',
                'Det 1 : a , 1, 1.0, 2.0, 0, 0, 3, 4, 1;\n',
                'Det 2 : b , 2, 2.0, 3.0, 0, 0, 3, 4, 1;\n',
                'Det 3 : c , 3, 1.0, 4.0, 0, 0, 3, 4, 1;\n',
                'Det 4 : NA , 3, 1000000.0, 99.5, 0, 0, 3, 4, 0;\n',
                '\n##\t\t K   MinYr   MaxYr   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax\n',
                'Bacon 0: FixT, 20, -1000, 1000000.0, 4147, 4145, 2.8, 1.2000000000000002, 1.5, 0.075, 1.5, 99.5;\n']
        self.assertCountEqual(victim[1:], goal[1:])  # Skip first line because datetime won't match.

    def test_run_baconmcmc(self):
        testcore_path = path.join(here, 'MSB2K.csv')
        c = snek.read_chron(testcore_path)
        fullrun_victim = baconwrap.run_baconmcmc(core_labid=c.labid, core_age=c.age, core_error=c.error,
                                                 core_depth=c.depth, depth_min=1.5, depth_max=99.5, cc=[1],
                                                 cc1='IntCal13', cc2='Marine13', cc3='SHCal13', cc4='ConstCal',
                                                 d_r=[0], d_std=[0], t_a=[3], t_b=[4], k=20,
                                                 minyr=-1000, maxyr=1e6, th01=4147, th02=4145,
                                                 acc_mean=20, acc_shape=1.5, mem_strength=4, mem_mean=0.7)

        niter_goal = 3432
        nsegs_goal = 20
        objective_mean_goal = 218.26212744522144
        theta_mean_goal = 4549.6406701631704
        w_mean_goal = 0.06277018029195805
        x0_mean_goal = 14.487713009032634
        xneg1_mean_goal = 15.65397709810606
        # Fuzzy to deal with vars across platforms.
        np.testing.assert_allclose(len(fullrun_victim['objective']), niter_goal, atol=50)
        np.testing.assert_allclose(len(fullrun_victim['w']), niter_goal, atol=50)
        np.testing.assert_allclose(len(fullrun_victim['x'][0]), niter_goal, atol=50)
        np.testing.assert_allclose(len(fullrun_victim['x'][-1]), niter_goal, atol=50)
        self.assertEqual(nsegs_goal, len(fullrun_victim['x']))
        np.testing.assert_allclose(fullrun_victim['objective'].mean(), objective_mean_goal, atol=1e-1)
        np.testing.assert_allclose(fullrun_victim['theta'].mean(), theta_mean_goal, atol=15)
        np.testing.assert_allclose(fullrun_victim['w'].mean(), w_mean_goal, atol=1e-2)
        np.testing.assert_allclose(fullrun_victim['x'][0].mean(), x0_mean_goal, atol=2)
        np.testing.assert_allclose(fullrun_victim['x'][-1].mean(), xneg1_mean_goal, atol=2)


# class TestRead14c(unittest.TestCase):
#
#     @unittest.skip('Test not written')
#     def test_read_14c(self):
#         self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
