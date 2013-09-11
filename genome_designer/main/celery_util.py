"""
Methods for interfacing with the Celery task queue management library.
"""

from errno import errorcode

from celery.task.control import inspect
from django.conf import settings


CELERY_ERROR_KEY = 'ERROR'

def get_celery_worker_status():
    """Checks whether celery is running and reports the error if not.

    Source: http://stackoverflow.com/questions/8506914/detect-whether-celery-is-available-running
    """
    if settings.BROKER_BACKEND == 'memory':
        # We are testing with in-memory celery. Celery is effectively running.
        return {}

    try:
        insp = inspect()
        d = insp.stats()
        if not d:
            d = { CELERY_ERROR_KEY: 'No running Celery workers were found.' }
    except IOError as e:
        msg = "Error connecting to the backend: " + str(e)
        if len(e.args) > 0 and errorcode.get(e.args[0]) == 'ECONNREFUSED':
            msg += ' Check that the RabbitMQ server is running.'
        d = { CELERY_ERROR_KEY: msg }
    except ImportError as e:
        d = { CELERY_ERROR_KEY: str(e)}
    return d
