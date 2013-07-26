"""
Tests celery tasks.
"""

from django.test import TestCase
from django.test.utils import override_settings

from main.tasks import add

class TestCelery(TestCase):
    """Tests async tasks.
    """

    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS = True,
            CELERY_ALWAYS_EAGER = True, BROKER_BACKEND = 'memory')
    def test_add(self):
        """Basic test using the placeholder add() tasks.
        """
        async_result = add.delay(2, 2)
        async_result.wait()
        self.assertEqual(4, async_result.result)
