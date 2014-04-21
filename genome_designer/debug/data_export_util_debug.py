"""
Debugging data_export_util.py.
"""

import csv
import os

# Since this script is intended to be used from the terminal, setup the
# environment first so that django and model imports work.
PWD = os.path.dirname(os.path.realpath(__file__ ))
scripts_path = os.path.join(PWD, '../scripts')
os.sys.path.append(scripts_path)
from utils import setup_django_env
setup_django_env()

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
