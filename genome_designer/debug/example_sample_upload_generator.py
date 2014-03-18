"""
Example script to generate a sample input.
"""

from collections import defaultdict
import csv
import os
import re

from well_id_generator import WellIdGenerator

DATA_DIR = '/home/ubuntu/recoli/fastq'
OUTPUT_FILE = '2014_03_17_recoli_test.tsv'

FIELD_NAMES = [
    'Sample_Name',
    'Plate_or_Group',
    'Well',
    'Read_1_Path',
    'Read_2_Path'
]

NAME_REGEX = '(recoli_lib4_[\w]+)_([1,2]{1}).*'


def _generate_sample_name(filename):
    match = re.match(NAME_REGEX, filename)
    return match.group(1), match.group(2)


def main():
    # First grab all the data.
    sample_name_to_data_dict = defaultdict(lambda: {})
    for fastq_filename in os.listdir(DATA_DIR):
        name, read = _generate_sample_name(fastq_filename)
        full_path = os.path.join(DATA_DIR, fastq_filename)
        sample_name_to_data_dict[name][read] = full_path

    well_gen = WellIdGenerator()

    # Write the output.
    with open(OUTPUT_FILE, 'w') as tsv_fh:
        writer = csv.DictWriter(tsv_fh, FIELD_NAMES, delimiter='\t')
        writer.writeheader()

        sorted_keys = sorted(sample_name_to_data_dict.keys())
        for key in sorted_keys:
            data = sample_name_to_data_dict[key]
            writer.writerow({
                'Sample_Name': key,
                'Plate_or_Group': '1',
                'Well': well_gen.next(),
                'Read_1_Path': data['1'],
                'Read_2_Path': data['2']
            })


if __name__ == '__main__':
    main()
