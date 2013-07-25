from celery import task

@task()
def add(x, y):
    """Example task.
    """
    return x + y
