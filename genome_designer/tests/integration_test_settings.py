"""Settings for tests.
"""

from tests.test_settings import *

TEST_RUNNER = 'test_suite_runner.IntegrationTestSuiteRunner'

ENABLED_VARIANT_CALLERS = [
    'freebayes',
    #'pindel',
    #'delly',
    'lumpy'
]
