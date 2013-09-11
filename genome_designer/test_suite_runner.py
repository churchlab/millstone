"""
Extends the nose runner to create a temporary filesystem during tests that is
subsequently torn down.
"""

import os
import shutil

from django_nose import NoseTestSuiteRunner

import settings

TEST_FILESYSTEM_DIR = 'temp_test_data'


class TempFilesystemTestSuiteRunner(NoseTestSuiteRunner):
    """Subclasses the nose test runner in order to use a temp filesystem.
    """

    def setup_databases(self):
        self.setup_temp_filesystem()
        return super(TempFilesystemTestSuiteRunner, self).setup_databases()

    def teardown_databases(self, *args, **kwargs):
        self.teardown_temp_filesystem()
        return super(TempFilesystemTestSuiteRunner, self).teardown_databases(
                *args, **kwargs)

    def setup_temp_filesystem(self):
        settings.MEDIA_ROOT = TEST_FILESYSTEM_DIR
        if os.path.exists(settings.MEDIA_ROOT):
            shutil.rmtree(settings.MEDIA_ROOT)
        os.mkdir(settings.MEDIA_ROOT)

    def teardown_temp_filesystem(self):
        assert settings.MEDIA_ROOT == TEST_FILESYSTEM_DIR
        if os.path.exists(settings.MEDIA_ROOT):
            shutil.rmtree(settings.MEDIA_ROOT)
