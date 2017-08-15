import unittest

import numpy as np

from snakebacon import suggest_accumulation_rate
from snakebacon.records import ChronRecord


class TestUtils(unittest.TestCase):

    def test_suggest_accumulation_rate(self):
        goal = 20
        testcore = ChronRecord(age=np.array([20, 50]), error=None, depth=np.array([1.5, 3]), labid=None)
        test_victim = suggest_accumulation_rate(testcore)
        self.assertEqual(test_victim, goal)


if __name__ == '__main__':
    unittest.main()
