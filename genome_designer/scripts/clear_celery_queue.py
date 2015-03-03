#!/usr/bin/env python

"""
Script that clears celery queue, discarding all tasks.
"""

from celery.task.control import discard_all

discard_all()
