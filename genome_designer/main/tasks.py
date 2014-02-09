"""
Analysis tasks handled by django-celery.

The current usage of celery is a bit hacky in that rather than decorating
individual tasks, we provide a generic_task() method that accepts a function
name and arguments, and then calls it.

NOTE: In order to use a function with generic_task, it is necessary to import
it here.
"""

from celery import task

@task()
def add(x, y):
    """Example task used for testing that celery is setup.
    """
    return x + y
