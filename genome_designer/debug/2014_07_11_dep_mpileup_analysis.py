"""
ucript to manually bin alignments against position. The idea is to manually
look for mutations that freebayes might miss because of the low probability
of the events that we are dealing with.
"""

import os
import sys

# Setup Django environment.
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'


from main.models import *
from utils import generate_safe_filename_prefix_from_label
from utils.bam_utils import filter_bam_file_by_row
from utils.samtools_utils import run_mpileup


def get_pEVOL_bam_path(bam_path):
    return os.path.splitext(bam_path)[0] + '.pEVOL_only.bam'


def filter_pEVOL_reads():
    p = Project.objects.get(uid='9a5fc6ec')
    sample_alignments = ExperimentSampleToAlignment.objects.filter(
            alignment_group__reference_genome__project=p)

    def is_pEVOL(line):
        parts = line.split('\t')
        rnext_col = parts[2]
        return rnext_col == 'pEVOL-bipA'

    for sa in sample_alignments:
        bam_dataset = get_dataset_with_type(sa, Dataset.TYPE.BWA_ALIGN)
        bam_path = bam_dataset.get_absolute_location()
        pEVOL_bam_path = get_pEVOL_bam_path(bam_path)
        filter_bam_file_by_row(bam_path, is_pEVOL, pEVOL_bam_path)


def generate_mpileups(idx_range=None):
    output_dir = '/dep_data/2014_07_11_pileup_analysis'
    p = Project.objects.get(uid='9a5fc6ec')
    sample_alignments = ExperimentSampleToAlignment.objects.filter(
            alignment_group__reference_genome__project=p)
    for idx, sa in enumerate(sample_alignments):
        if idx_range and not idx in idx_range:
            continue
        print 'Running %d of %d: %s' % (
                idx + 1, len(sample_alignments), sa.experiment_sample.label)
        bam_dataset = get_dataset_with_type(sa, Dataset.TYPE.BWA_ALIGN)
        bam_path = bam_dataset.get_absolute_location()
        pEVOL_bam_path = get_pEVOL_bam_path(bam_path)

        ref_genome_fasta_location = get_dataset_with_type(
                sa.alignment_group.reference_genome,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

        output_filename = generate_safe_filename_prefix_from_label(
                sa.experiment_sample.label + '_' + sa.uid) + '.mpileup'
        output_path = os.path.join(output_dir, output_filename)

        run_mpileup(pEVOL_bam_path, ref_genome_fasta_location, output_path,
                coverage_only=False)


def main():
    # filter_pEVOL_reads()
    generate_mpileups([6,7])


if __name__ == '__main__':
    main()
