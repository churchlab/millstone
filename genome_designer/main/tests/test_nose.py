"""Tests to see how nose (our testing helper) works.
"""

import time

from django.test import TestCase


class TestIsParallel(TestCase):
    """Testing to see whether nose runs tests in parallel.
    """

    def test_one(self):
        print 'test_one 1'
        time.sleep(2)

        print 'test_one 2'
        time.sleep(3)

        print 'test_one 3'
        time.sleep(4)

    def test_two(self):
        print 'test_two 1'
        time.sleep(1)
        print 'test_two 2'
        time.sleep(1)
        print 'test_two 3'
        time.sleep(1)
