""" Test for EBCIC

This is part of https://github.com/KazKobara/ebcic/ .
"""

import doctest
import unittest
import ebcic


# unittest
class TestFunc(unittest.TestCase):
    def test_func(self):
        expected = 0
        actual = ebcic.test_all()
        self.assertEqual(expected, actual)


# doctest
def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite('ebcic'))
    return tests


if __name__ == '__main__':
    unittest.main()
