"""Tests for optmage_util.py.
"""

from django.test import TestCase

from utils.optmage_util import print_mage_oligos


class TestOptmageUtil(TestCase):

    def test_print_mage_oligos(self):
        print_mage_oligos(None)
