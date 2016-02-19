import os

from Bio import SeqIO
from celery import task
import numpy as np
import pysam

from genome_finish.graph_contig_placement import get_fasta
from genome_finish.millstone_de_novo_fns import get_altalign_reads
from main.models import Dataset
from utils.bam_utils import index_bam
from utils.data_export_util import export_var_dict_list_as_vcf
from utils.import_util import add_dataset_to_entity

CALLER_STRING = 'COV_DETECT_DELETION'


@task
def cov_detect_deletion_make_vcf(sample_alignment):
    """Uses coverage data to call large deletions and
    creates a VCF_COV_DETECT_DELETIONS dataset for the sample alignment

    Args:
        sample_alignment: ExperimentSampleToAlignment instance
    """

    print "Generating coverage data\n"
    chrom_regions = get_deleted_regions(sample_alignment)
    var_dict_list = make_var_dict_list(
            chrom_regions,
            get_fasta(sample_alignment.alignment_group.reference_genome))

    if var_dict_list:

        dataset_query = sample_alignment.dataset_set.filter(
                type=Dataset.TYPE.VCF_COV_DETECT_DELETIONS)

        if dataset_query:
            dataset_query.delete()

        vcf_path = os.path.join(
            sample_alignment.get_model_data_dir(),
            'cov_detect_deletion.vcf')

        # Write variant dicts to vcf
        export_var_dict_list_as_vcf(var_dict_list, vcf_path,
                sample_alignment, caller_string=CALLER_STRING)

        # Make dataset for contigs vcf
        new_dataset = add_dataset_to_entity(
                sample_alignment,
                Dataset.TYPE.VCF_COV_DETECT_DELETIONS,
                Dataset.TYPE.VCF_COV_DETECT_DELETIONS,
                vcf_path)

        new_dataset.save()


def make_var_dict_list(chrom_regions, ref_fasta):

    with open(ref_fasta, 'r') as fh:
        record_dict = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))

    var_dict_list = []
    for chrom, regions in chrom_regions.items():
        if not regions:
            continue

        record = record_dict[chrom]
        for region in regions:
            var_dict = {
                'chromosome': chrom,
                'pos': region[0],
                'ref_seq': str(record.seq[region[0]:region[1]]),
                'alt_seq': '',
            }
            var_dict_list.append(var_dict)

    return var_dict_list


def get_deleted_regions(sample_alignment, cov_cutoff=5):

    chrom_to_cov_list = per_base_cov_stats_opt(sample_alignment)

    chrom_regions = {}
    for chrom, cov_dict in chrom_to_cov_list.items():
        depths = cov_dict['depths']
        altaligns = cov_dict['altaligns']
        unique = depths - altaligns

        low_cov_regions = get_low_cov_regions(
                depths, cov_cutoff)

        smoothed_regions = smoothed_deletions(
                unique, low_cov_regions)

        chrom_regions[chrom] = smoothed_regions

    return chrom_regions


def get_low_cov_regions(cov_arr, min_cov):
    split_indices = np.where(np.diff(cov_arr < min_cov) != 0)[0] + 1
    split_indices = list(split_indices)

    if cov_arr[0] < min_cov:
        split_indices.insert(0, 0)
    if cov_arr[-1] < min_cov:
        split_indices.append(len(cov_arr))

    assert len(split_indices) % 2 == 0

    return [(split_indices[i], split_indices[i+1]) for i in
               range(0, len(split_indices), 2)]


def smoothed_deletions(unique_arr, low_cov_regions):

    SMOOTHING_COV_CUTOFF = 3
    EXP_COV_DECAY_HALF_LIFE = 500
    LARGE_DEL_MAX_SMOOTH_DIST = 1000
    LARGE_DEL_MIN_DEL_LEN = 2000

    cov_mean = np.mean(unique_arr)

    def is_smoothable(curr_region, next_region):

        dist = next_region[0] - curr_region[1]
        curr_len = curr_region[1] - curr_region[0]
        next_len = next_region[1] - next_region[0]
        del_cov = np.mean(unique_arr[curr_region[1]:next_region[0]])

        # Cast to int from numpy int64
        dist = int(dist)

        # Allow smoothing Between large deleted regions
        if (dist < LARGE_DEL_MAX_SMOOTH_DIST and
                min(curr_len, next_len) > LARGE_DEL_MIN_DEL_LEN):
            return True

        return del_cov < max(
                cov_mean * 2 ** (-dist / EXP_COV_DECAY_HALF_LIFE),
                SMOOTHING_COV_CUTOFF)

    regions = low_cov_regions[:]
    is_same = False
    i = 0
    while(len(regions) > 1):
        curr_region = regions[i]
        next_region = regions[i+1]

        if is_smoothable(curr_region, next_region):
            del regions[i]
            del regions[i]
            concat_region = (curr_region[0], next_region[1])
            regions.insert(i, concat_region)
            is_same = False

        if i >= len(regions) - 2:
            if is_same:
                return regions
            is_same = True
            i = 0
        else:
            i += 1

    return regions


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
        print "Writing bam dataset of alternative alignment reads\n"
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
