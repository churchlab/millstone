"""
Generate the sample input file that can be passed to Millstone.
"""

import csv
import os
import re


MILLSTONE_TARGETS_FILE = 'millstone_targets_file.tsv'

PWD = os.path.dirname(os.path.realpath(__file__))

def iterate_samples(writer):
    barcode_set = set()
    file_count = 0
    for maybe_sample_lib_dir in os.listdir('.'):
        if re.match('Sample_LIB', maybe_sample_lib_dir):
            root_dir = os.path.join(PWD, maybe_sample_lib_dir)
            sample_contents = os.listdir(root_dir)
            read_1_path = None
            read_2_path = None
            for f in sample_contents:
                if re.search('L001_R1.fastq.bz2', f) or re.search('L001_R1.fastq$', f):
                    read_1_name_for_find_barcode = f
                    read_1_path = os.path.join(root_dir, f)
                elif re.search('L001_R2.fastq.bz2', f) or re.search('L001_R2.fastq$', f):
                    read_2_path = os.path.join(root_dir, f)
            assert read_1_path, "Failed for %s" % root_dir
            assert read_2_path, "Failed for %s" % root_dir

            # Use barcode for sample name.
            barcode = re.search('[ACGT]{6}', read_1_name_for_find_barcode).group()
            if barcode in barcode_set:
                print barcode
            barcode_set.add(barcode)
            file_count += 1

            writer.writerow({
                'Sample_Name': barcode,
                'Plate_or_Group': 1,
                'Well': 0,
                'Read_1_Path': read_1_path,
                'Read_2_Path': read_2_path,
            })
    print len(barcode_set)
    print file_count

def main():
    FIELD_NAMES = ['Sample_Name', 'Plate_or_Group',	'Well', 'Read_1_Path', 'Read_2_Path']

    with open(MILLSTONE_TARGETS_FILE, 'w') as output_fh:
        writer = csv.DictWriter(output_fh, FIELD_NAMES, delimiter='\t')
        writer.writeheader()
        iterate_samples(writer)


if __name__ == '__main__':
    main()
