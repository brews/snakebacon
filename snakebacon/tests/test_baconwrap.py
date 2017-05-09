import unittest
import numpy as np
from snakebacon.bacon import baconwrap

# TODO(brews): Build cython as setup via distutils.core.run_setup
# Be sure to run python setup.py build_ext --inplace

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
                'Det 0 : a , 1, 1, 2, 0, 0, 3, 4, 1;\n',
                'Det 1 : b , 2, 2, 3, 0, 0, 3, 4, 1;\n',
                'Det 2 : c , 3, 1, 4, 0, 0, 3, 4, 1;\n',
                '\n##\t\t K   MinYr   MaxYr   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax\n',
                'Bacon 0: FixT, 20, -1000, 1000000.0, 4147, 4145, 2.8, 1.2000000000000002, 1.5, 0.075, 1.5, 99.5;\n']
        self.assertCountEqual(victim[1:], goal[1:])  # Skip first line because datetime won't match.


# class TestRead14c(unittest.TestCase):
#
#     @unittest.skip('Test not written')
#     def test_read_14c(self):
#         self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()