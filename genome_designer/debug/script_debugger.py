"""
Convenience module for debugging scripts.
"""
import os
import sys

# Make this script runnable from command line.
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from main.celery_util import assert_celery_running
from main.models import *
from variants import materialized_view_manager

def debug_materialized_view_manager(ag):
    ref = ReferenceGenome.objects.get(alignmentgroup__uid=ag)
    mvm = materialized_view_manager.MeltedVariantMaterializedViewManager(ref)
    mvm.create()


def debug_celery():
    assert_celery_running()


def main():
    debug_celery()


if __name__ == '__main__':
    debug_materialized_view_manager(sys.argv[1])
