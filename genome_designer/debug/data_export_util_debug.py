"""
Debugging data_export_util.py.
"""

import csv
import os
import sys

# Setup Django environment.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from data_export_util import export_melted_variant_view
from main.models import ReferenceGenome

OUTPUT_CSV = 'fix_recoli_melted_variant_data.csv'

def main():
    ref_genome = ReferenceGenome.objects.get(uid='29517b77')
    with open(OUTPUT_CSV, 'w') as csvfile_fh:
        variant_id_list = []
        export_melted_variant_view(ref_genome, variant_id_list, csvfile_fh)

if __name__ == '__main__':
    main()
