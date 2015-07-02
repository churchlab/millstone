"""Functions for running lumpy and processing output.
"""

import glob
import os
import subprocess

from django.conf import settings
import vcf

from main.model_utils import get_dataset_with_type
from main.models import Dataset
from pipeline.variant_calling.common import add_vcf_dataset
from pipeline.variant_calling.common import get_common_tool_params
from pipeline.variant_calling.common import TOOL_LUMPY
from pipeline.read_alignment import get_discordant_read_pairs
from pipeline.read_alignment import get_split_reads
from pipeline.variant_calling.common import process_vcf_dataset
from utils import uppercase_underscore


def run_lumpy(
        fasta_ref, sample_alignments, vcf_output_dir, vcf_output_filename,
        alignment_type, **kwargs):
    """Runs lumpy.
    """
    # NOTE: Only supporting single sample alignment for now. Previously we
    # tried to use lumpy for multiple sample alignments but the machine would
    # run out of memory so we are going to limit functionality to single
    # alignment only for now.
    assert len(sample_alignments) == 1

    # Get relevant files. Note this is written to handle more than 1 sample
    # although right now we are not running lumpy on more than one sample at
    # a time as enforced by the assert above.
    bam_file_list = []
    bam_disc_file_list = []
    bam_sr_file_list = []
    for sa in sample_alignments:
        bam_dataset = get_dataset_with_type(sa, Dataset.TYPE.BWA_ALIGN)
        bam_file_list.append(bam_dataset.get_absolute_location())

        # Get or create discordant reads.
        bam_disc_dataset = get_discordant_read_pairs(sa)
        bam_disc_file_list.append(bam_disc_dataset.get_absolute_location())

        # Get or create split reads.
        bam_sr_dataset = get_split_reads(sa)
        bam_sr_file_list.append(bam_sr_dataset.get_absolute_location())

    # Run Lumpy Express.
    subprocess.check_call([
        settings.LUMPY_EXPRESS_BINARY,
        '-B', ','.join(bam_file_list),
        '-S', ','.join(bam_sr_file_list),
        '-D', ','.join(bam_disc_file_list),
        '-o', vcf_output_filename,
        '-P' # get probability distributions, required for merge
    ])

    return True  # success


def filter_lumpy_vcf(original_vcf_path, new_vcf_path):
    """Filters lumpy vcf to get rid of noisy values.

    Args:
        original_vcf_path: Full path to starting vcf.
        new_vcf_path: Path where new vcf will be written.

    NOTE: Now that we are using lumpyexpress, we might not need this anymore.
    """
    with open(original_vcf_path) as orig_vcf_fh:
        with open(new_vcf_path, 'w') as new_vcf_fh:
            vcf_reader = vcf.Reader(orig_vcf_fh)
            vcf_writer = vcf.Writer(new_vcf_fh, vcf_reader)
            for record in vcf_reader:
                # If record fails any filter, continue to next record without
                # writing.
                if int(record.INFO['DP']) < 10:
                    continue
                vcf_writer.write_record(record)


def merge_lumpy_vcf(alignment_group):
    """Merge lumpy outputs run on individual samples.

    If no lumpy, returns None.
    """
    common_params = get_common_tool_params(alignment_group)
    partial_vcf_output_dir = os.path.join(
            common_params['output_dir'], TOOL_LUMPY)

    # Glob all the partial (region-specific) vcf files.
    # Assert that there is at least one.
    vcf_output_filename_prefix = os.path.join(partial_vcf_output_dir,
            uppercase_underscore(common_params['alignment_type']) +
            '.partial.*.vcf')
    partial_vcf_files = glob.glob(vcf_output_filename_prefix)
    if not len(partial_vcf_files):
        return None

    # l_sort combines the vcfs.
    l_sort_output_vcf_filepath = os.path.join(
            partial_vcf_output_dir, 'lumpy_l_sort_output.vcf')
    l_sort_cmd_list = [settings.LUMPY_L_SORT_BINARY] + partial_vcf_files
    with open(l_sort_output_vcf_filepath, 'w') as l_sort_output_fh:
        subprocess.check_call(l_sort_cmd_list, stdout=l_sort_output_fh)

    # l_merge merges lines representing the same variant.
    merged_vcf_filepath = os.path.join(
            partial_vcf_output_dir,
            uppercase_underscore(common_params['alignment_type']) + '.vcf')
    l_merge_cmd_list = [
        settings.LUMPY_L_MERGE_BINARY,
        '-i', l_sort_output_vcf_filepath]
    with open(merged_vcf_filepath, 'w') as l_merge_output_fh:
        subprocess.check_call(l_merge_cmd_list, stdout=l_merge_output_fh)

    # Create Dataset pointing to merged vcf file.
    vcf_dataset_type = Dataset.TYPE.VCF_LUMPY
    vcf_dataset = add_vcf_dataset(
            alignment_group, vcf_dataset_type, merged_vcf_filepath)

    # Parse VCF to add variants to database.
    process_vcf_dataset(alignment_group, vcf_dataset_type)

    # # Remove the partial vcfs.
    # for filename in partial_vcf_files:
    #     os.remove(filename)

    return vcf_dataset
