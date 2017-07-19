import unittest

import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal

from snakebacon import DatedProxyRecord


class TestDatedProxyRecord(unittest.TestCase):
    def setUp(self):
        self.testdummy_median = DatedProxyRecord(pd.DataFrame({'depth': [1, 2], 'a': [2, 3]}), age=[10, 11])
        self.testdummy_ensemble = DatedProxyRecord(pd.DataFrame({'depth': [1, 2], 'a': [2, 3]}),
                                                   age=[[10, 10, 10], [11, 11, 11]])

    def test_n_members(self):
        goal_median = 1
        goal_ensemble = 3
        self.assertEqual(goal_median, self.testdummy_median.n_members())
        self.assertEqual(goal_ensemble, self.testdummy_ensemble.n_members())

    def test_to_pandas(self):
        goal_median = pd.DataFrame({'a': [2, 3], 'age': [10, 11], 'depth': [1, 2]})
        goal_ensemble = pd.DataFrame({'a': [2, 3] * 3, 'depth': [1, 2] * 3,
                                      'mciter': np.repeat(np.arange(3), 2),
                                      'age': [10, 11] * 3, })
        assert_frame_equal(goal_median.sort_index(axis=1), self.testdummy_median.to_pandas().sort_index(axis=1))
        assert_frame_equal(goal_ensemble.sort_index(axis=1), self.testdummy_ensemble.to_pandas().sort_index(axis=1))


if __name__ == '__main__':
    unittest.main()
