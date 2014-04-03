"""
Convenience module for debugging scripts.
"""
import os
import sys

# Make this script runnable from command line.
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from main.models import *
from variants import materialized_view_manager

def main():
    ref = ReferenceGenome.objects.get(alignmentgroup__uid='c359d2cf')
    mvm = materialized_view_manager.MeltedVariantMaterializedViewManager(ref)
    mvm.create()

if __name__ == '__main__':
    main()
