import os

import numpy as np
import pysam

from genome_finish.millstone_de_novo_fns import get_altalign_reads
from main.models import Dataset
from utils.bam_utils import index_bam
from utils.import_util import add_dataset_to_entity


def make_altalign_dataset(sample_alignment):

    sample_alignment_bam = sample_alignment.dataset_set.get(
                type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()
    alignment_file_prefix = os.path.join(
        sample_alignment.get_model_data_dir(),
        'bwa_align')
    altalign_bam = '.'.join([
                alignment_file_prefix,
                'altalign',
                'bam'
                ])
    get_altalign_reads(sample_alignment_bam, altalign_bam)

    return add_dataset_to_entity(
                sample_alignment,
                Dataset.TYPE.BWA_ALTALIGN,
                Dataset.TYPE.BWA_ALTALIGN,
                filesystem_location=altalign_bam)


def per_base_cov_stats_opt(sample_alignment):

    sample_alignment_bam = sample_alignment.dataset_set.get(
            type=Dataset.TYPE.BWA_ALIGN).get_absolute_location()
    bamfile = pysam.AlignmentFile(sample_alignment_bam, 'rb')

    chrom_list = bamfile.references
    chrom_lens = bamfile.lengths
    c_starts = [0] * len(chrom_list)
    c_ends = chrom_lens

    chrom_to_cov_list = {}

    for chrom, c_start, c_end in zip(chrom_list, c_starts, c_ends):

        depth_arr = np.zeros(c_end)
        for pileup_col in bamfile.pileup(chrom,
                start=c_start, end=c_end, truncate=True):

            # number of segments aligned to this position
            depth = pileup_col.nsegments
            depth_arr[pileup_col.reference_pos] = depth

        chrom_to_cov_list[chrom] = {
                'depths': depth_arr
        }


    # Do altaligns
    altalign_dataset_query = sample_alignment.dataset_set.filter(
            type=Dataset.TYPE.BWA_ALTALIGN)

    if altalign_dataset_query.count():
        assert altalign_dataset_query.count() == 1
        altalign_dataset = altalign_dataset_query[0]
    else:
        altalign_dataset = make_altalign_dataset(sample_alignment)

    altalign_bam_path = altalign_dataset.get_absolute_location()
    index_bam(altalign_bam_path)

    bamfile = pysam.AlignmentFile(
            altalign_dataset.get_absolute_location(), 'rb')

    chrom_list = bamfile.references
    chrom_lens = bamfile.lengths
    c_starts = [0] * len(chrom_list)
    c_ends = chrom_lens


    for chrom, c_start, c_end in zip(chrom_list, c_starts, c_ends):
        depth_arr = np.zeros(c_end)

        for pileup_col in bamfile.pileup(chrom,
                start=c_start, end=c_end, truncate=True):

            # number of segments aligned to this position
            depth = pileup_col.nsegments
            depth_arr[pileup_col.reference_pos] = depth

        chrom_to_cov_list[chrom]['altaligns'] = depth_arr

    return chrom_to_cov_list
