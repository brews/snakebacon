import unittest
import numpy as np
from snakebacon.records import DateRecord

class TestDateRecordMethods(unittest.TestCase):

    def test_suggest_accumulation_rate(self):
        goal = 20
        testcore = DateRecord(age=np.array([20, 50]), error=None, depth=np.array([1.5, 3]), labid=None)
        test_victim = testcore.suggest_accumulation_rate()
        self.assertEqual(test_victim, goal)

    @unittest.skip('Test not written')
    def test_suggest_thick(self):
        pass

    @unittest.skip('Test not written')
    @unittest.skip('')
    def test_calibrate_dates(self):
        pass


if __name__ == '__main__':
    unittest.main()
