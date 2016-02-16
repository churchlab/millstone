"""
Simple test while setting up travis-ci.
"""

from django.test import TestCase


class TestTravis(TestCase):

    def test_travis(self):
        self.assertTrue(True)
