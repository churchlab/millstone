"""
Extends the nose runner to create a temporary filesystem during tests that is
subsequently torn down.
"""

import os
import shutil

from django.conf import settings
from django_nose import NoseTestSuiteRunner
from djcelery_testworker import run_celery_test_worker


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
        settings.MEDIA_ROOT = settings.TEST_FILESYSTEM_DIR
        if os.path.exists(settings.MEDIA_ROOT):
            shutil.rmtree(settings.MEDIA_ROOT)
        os.mkdir(settings.MEDIA_ROOT)

    def teardown_temp_filesystem(self):
        assert settings.MEDIA_ROOT == settings.TEST_FILESYSTEM_DIR
        if os.path.exists(settings.MEDIA_ROOT):
            shutil.rmtree(settings.MEDIA_ROOT)


class CustomTestSuiteRunner(TempFilesystemTestSuiteRunner):
    """Our custom TestSuiteRunner.
    """

    def setup_test_environment(self, **kwargs):
        setup_test_environment_common()

        # Don't use Celery for tests. If we ever want to test with Celery,
        # we'll need to create a different TestSuiteRunner.
        settings.CELERY_EAGER_PROPAGATES_EXCEPTIONS = True
        settings.CELERY_ALWAYS_EAGER = True
        settings.BROKER_BACKEND = 'memory'

        return super(CustomTestSuiteRunner, self).setup_test_environment()


class IntegrationTestSuiteRunner(TempFilesystemTestSuiteRunner):
    """TestSuiteRunner for integration tests.

    Main feature: Starts a celery worker connected to the test database.
    """

    def setup_test_environment(self, **kwargs):
        setup_test_environment_common()

        self.celeryd = run_celery_test_worker()

        return super(IntegrationTestSuiteRunner, self).setup_test_environment()

    def teardown_test_environment(self):
        self.celeryd.kill()
        return super(IntegrationTestSuiteRunner, self).teardown_test_environment()


def setup_test_environment_common():
    """Common setup procedures.
    """
    # Catch leftover .pyc from, say, changing git branches.
    assert_no_orphaned_pyc_files('.')

    # Make sure the startup function works.
    # As of implementation, this function adds a custom function to
    # our Postgres database.
    from django.db.models.signals import post_syncdb
    import main
    post_syncdb.connect(handle_post_syncdb_startup, sender=main.models)


def handle_post_syncdb_startup(sender, **kwargs):
    from main import startup
    startup.run()


def assert_no_orphaned_pyc_files(start_dir):
    """Walk from the start directory through python pkg dirs, looking for
    .pyc files with no corresponding .py file.

    This is to avoid the bug where we delete a module but the .pyc file stays
    and allows tests to pass on the local machine.

    Raises AssertionError if an orphaned .pyc file is found.
    """
    orphaned_files = []
    for (dirpath, dirnames, filenames) in os.walk(start_dir, topdown=True):
        # NOTE: Python docs say it is okay to modify dirnames in-place in order
        # to prune the directories that we are walking.
        if not '__init__.py' in filenames:
            # Don't recurse any further.
            while len(dirnames):
                dirnames.pop()
            # No need to check this dir.
            continue


        # Build a set of names ending in .py.
        py_set = set()
        for filename in filenames:
            if os.path.splitext(filename)[1] == '.py':
                py_set.add(os.path.splitext(filename)[0])

        # Check that all .pyc files have a corresponding .py file.
        for filename in filenames:
            if os.path.splitext(filename)[1] == '.pyc':
                if not os.path.splitext(filename)[0] in py_set:
                    full_orphan_path = os.path.join(dirpath, filename)
                    orphaned_files.append(full_orphan_path)

    if len(orphaned_files):
        raise AssertionError, (
                "The following files are orphaned .pyc files:\n\n%s\n\n" %
                        '\n'.join(orphaned_files) +
                "Perhaps you moved the file and meant to delete them?")
