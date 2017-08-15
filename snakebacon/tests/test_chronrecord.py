import unittest

import numpy as np

from snakebacon.records import ChronRecord


class TestChronRecordMethods(unittest.TestCase):

    @unittest.skip('Test not written. Need methods to test equality.')
    def test__repr__(self):
        testcore = ChronRecord(age=np.array([20, 50]), error=None, depth=np.array([1.5, 3]), labid=None)
        self.assertTrue(testcore == eval(repr(testcore)))


if __name__ == '__main__':
    unittest.main()
