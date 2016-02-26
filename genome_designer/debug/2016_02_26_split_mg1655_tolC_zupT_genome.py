import os

from Bio import SeqIO


TEST_DATA_DIR = '../test_data/full_vcf_test_set'

YGIB_OUTPUT = os.path.join(TEST_DATA_DIR, 'mg1655_partial_ygiB_to_zupT.gb')

TOLC_OUTPUT = os.path.join(TEST_DATA_DIR, 'mg1655_partial_tolC.gb')

full_seq_record = SeqIO.read(
        os.path.join(TEST_DATA_DIR, 'mg1655_tolC_through_zupT.gb'),
        'genbank')


# Get the cut position.
cut_position = None
for f in full_seq_record.features:
    if f.type == 'CDS' and f.qualifiers['gene'][0] == 'ygiB':
        cut_position = f.location.start - 50
        break
assert cut_position

mg1655_partial_tolC_seq_record = full_seq_record[:cut_position]
with open(TOLC_OUTPUT, 'w') as fh:
    SeqIO.write(mg1655_partial_tolC_seq_record, fh, 'genbank')

mg1655_partial_ygiB_seq_record = full_seq_record[cut_position:]
with open(YGIB_OUTPUT, 'w') as fh:
    SeqIO.write(mg1655_partial_ygiB_seq_record, fh, 'genbank')
