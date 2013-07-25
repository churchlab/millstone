"""
Tests for the main app.

We use this module to discover all other tests of the form 'test*.py'
located within the same directory containing this module.
"""

import os
import unittest

PWD = os.path.dirname(os.path.realpath(__file__ ))

def suite():
    return unittest.TestLoader().discover(PWD, pattern='test*.py')
