"""
Re-running de novo assembly, this time including reads that map to mobile elements.
"""

import os
import sys

# Setup Django environment.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

from Bio import SeqIO

from experimental.de_novo_assembly import run_velvet
from main.models import *



def identify_intervals(ag):
    # First identify intervals that map to mobile elements.
    genbank_filepath = get_dataset_with_type(ag.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()

    # Extract the proper genome record.
    genome_record = None
    with open(genbank_filepath) as input_fh:
        genome_record_list = SeqIO.parse(input_fh, 'genbank')
        for rec in genome_record_list:
            if rec.name == 'CP006698':
                genome_record = rec
    assert genome_record

    # Pick out the intervals we want:
    #    * mobile elements
    #    * lon gene
    intervals = []
    found_lon = False
    for f in genome_record.features:
        if f.type == 'mobile_element':
            intervals.append((f.location.start, f.location.end))
        if (f.type == 'gene' and 'gene' in f.qualifiers and
                f.qualifiers['gene'][0] == 'lon'):
            found_lon = True
            intervals.append((f.location.start, f.location.end))
    assert found_lon
    assert 47 == len(intervals)

    # Add buffer to each interval in case reads start before or after.
    buffer_size = 150
    def _add_buffer(i):
        return (
                max(i[0] - buffer_size, 0),
                min(i[1] + buffer_size, len(genome_record))
        )
    intervals = [_add_buffer(i) for i in intervals]

    return intervals


def main():
    ag = AlignmentGroup.objects.get(uid='edc74a3d')

    intervals = identify_intervals(ag)

    for idx, sa in enumerate(ag.experimentsampletoalignment_set.all()):
        print idx + 1, 'of', ag.experimentsampletoalignment_set.count()
        run_velvet(sa, force_include_reads_in_intervals=intervals)


if __name__ == '__main__':
    main()
