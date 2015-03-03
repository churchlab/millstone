"""
Script to get all histo files.
"""

import os
import shutil
import sys

# Setup Django environment.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'


from main.model_utils import get_dataset_with_type
from main.models import *

OUTPUT_DIR = 'histo_files'


def main():
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    for sa in ExperimentSampleToAlignment.objects.all():
        histo_dataset = get_dataset_with_type(sa,
                Dataset.TYPE.LUMPY_INSERT_METRICS_HISTOGRAM)
        histo_dataset_full_path = histo_dataset.get_absolute_location()

        # Update Dataset name.
        histo_dataset_name = (
                os.path.split(os.path.split(histo_dataset_full_path)[0])[1] + '.txt')

        # Copy.
        new_full_path = os.path.join(OUTPUT_DIR, histo_dataset_name)
        shutil.copyfile(histo_dataset_full_path, new_full_path)


if __name__ == '__main__':
    main()
