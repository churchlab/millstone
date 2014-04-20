#!/usr/bin/env python

"""Script to run integration tests.

django_nose pukes with the --testrunner flag, so one way to use our
IntegrationTestRunner with django_nose is to override the TEST_RUNNER
value in settings.
"""

import os
import sys

if __name__ == "__main__":
    # Update the path so we can import settings.
    sys.path.append(os.path.join(os.path.dirname(
            os.path.realpath(__file__)), '../'))
    sys.path.append('../')

    # Set settings module so we can do Django stuff.
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")

    from django.conf import settings
    from django.core.management import execute_from_command_line

    # Override with our custom test runner.
    settings.TEST_RUNNER = 'test_suite_runner.IntegrationTestSuiteRunner'

    # Add 'test' to argv so from here on the command looks like
    #     manage.py test ...
    sys.argv = sys.argv[:1] + ['test', '--settings=tests.test_settings'] + sys.argv[1:]

    execute_from_command_line(sys.argv)
