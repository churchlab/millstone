"""
Analysis tasks handled by django-celery.

The current usage of celery is a bit hacky in that rather than decorating
individual tasks, we provide a generic_task() method that accepts a function
name and arguments, and then calls it.

NOTE: In order to use a function with generic_task, it is necessary to import
it here.
"""

from celery import task

from scripts.alignment_pipeline import *
from scripts.snp_callers import *
from s3 import project_s3_persisted

@task()
def generic_task(fn_name, project, args_list):
    """Generic placeholder that allows running any task without having to
    explicitly register it.

    NOTE: The parent modules for any methods called here must be imported here.
    """
    with project_s3_persisted(project):
        exec('fn = ' + fn_name)
        fn(*args_list)

@task()
def add(x, y):
    """Example task used for testing that celery is setup.
    """
    return x + y
