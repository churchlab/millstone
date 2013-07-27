"""
Tests for models.py.
"""

import os

from django.test import TestCase

from main.models import clean_filesystem_location
import settings

FAKE_MEDIA_ROOT = 'fake/media/root'

class TestModels(TestCase):

    def setUp(self):
        """Override.
        """
        self.original_media_root = settings.MEDIA_ROOT
        settings.MEDIA_ROOT = FAKE_MEDIA_ROOT

    def tearDown(self):
        """Override.
        """
        settings.MEDIA_ROOT = self.original_media_root

    def test_clean_filesystem_location(self):
        FAKE_ABS_ROOT = '/root/of/all/evil'
        EXPECTED_CLEAN_URL = 'projects/blah'
        dirty_full_url = os.path.join(FAKE_ABS_ROOT, FAKE_MEDIA_ROOT,
                EXPECTED_CLEAN_URL)
        self.assertEqual(EXPECTED_CLEAN_URL,
                clean_filesystem_location(dirty_full_url))
