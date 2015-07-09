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
    print 'RUNNING LUMPY...'

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

    lumpy_cmd = [
        settings.LUMPY_EXPRESS_BINARY,
        '-B', ','.join(bam_file_list),
        '-S', ','.join(bam_sr_file_list),
        '-D', ','.join(bam_disc_file_list),
        '-o', vcf_output_filename,
        '-P' # get probability distributions, required for merge
    ]

    print ' '.join(lumpy_cmd)

    # Run Lumpy Express.
    lumpy_error_output = vcf_output_filename + '.error'
    with open(lumpy_error_output, 'w') as error_output_fh:
        subprocess.check_call(lumpy_cmd, stderr=error_output_fh)

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
    l_merge_output_path = os.path.join(
            partial_vcf_output_dir, 'lumpy_l_merge_output.vcf')
    l_merge_cmd_list = [
        settings.LUMPY_L_MERGE_BINARY,
        '-i', l_sort_output_vcf_filepath]
    with open(l_merge_output_path, 'w') as l_merge_output_fh:
        subprocess.check_call(l_merge_cmd_list, stdout=l_merge_output_fh)

    # Post-processing following l-merge.
    merged_vcf_filepath = os.path.join(
            partial_vcf_output_dir,
            uppercase_underscore(common_params['alignment_type']) + '.vcf')
    process_vcf_post_l_merge(l_merge_output_path, merged_vcf_filepath)

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


def process_vcf_post_l_merge(l_merge_output_vcf_path, processed_vcf_path):
    """Processes vcf following l_merge.

    The output of l_merge doesn't have a column per sample with GT information,
    which is the format that vcf_parser expects. Instead, l_merge places the
    information into the INFO string. So we need to parse this and output
    the properly formatted vcf file.
    """
    with open(l_merge_output_vcf_path) as l_merge_output_fh:
        with open(processed_vcf_path, 'w') as processed_vcf_fh:
            vcf_reader = vcf.Reader(l_merge_output_fh)

            # Fix info strings.
            _update_info_string_number(vcf_reader, 'SVTYPE', -1)
            _update_info_string_number(vcf_reader, 'SVLEN', -1)

            # Fix format header.
            orig = vcf_reader.formats['SU']
            vcf_reader.formats['DP'] = vcf.parser._Format(
                    'DP', orig.num, orig.type, orig.desc)
            del vcf_reader.formats['SU']

            # Make column headers match what's expected by vcf_parser.
            # l_merge output is missing FORMAT column header, and columns
            # for each sample.
            if not 'FORMAT' in vcf_reader._column_headers:
                vcf_reader._column_headers.append('FORMAT')
            vcf_reader.samples = [
                    x['ID'] for x in vcf_reader.metadata['SAMPLE']]

            # Writer object using Reader as template.
            vcf_writer = vcf.Writer(processed_vcf_fh, vcf_reader)

            # Format each record with correct setting.
            for record in vcf_reader:
                # Per-sample values.
                record.FORMAT = 'GT:DP'

                # vcf.model._Call requires data as a hashable type so follow
                # vcf internal code pattern of making a tuple.
                calldata_tuple_type = vcf.model.make_calldata_tuple(
                        record.FORMAT.split(':'))

                samples_with_sv = [
                        x.split(':')[0] for x in record.INFO['SNAME']]

                if 'SULIST' in record.INFO:
                    dp_list = [x.split(':')[0] for x in record.INFO['SULIST']]
                else:
                    dp_list = record.INFO['SU']

                # Parse the record
                record_samples = []
                for sample_id in vcf_reader.samples:
                    try:
                        sample_idx = samples_with_sv.index(sample_id)

                        sample_data = calldata_tuple_type(
                                GT='1/1',
                                DP=dp_list[sample_idx])
                    except ValueError:
                        sample_data = calldata_tuple_type(GT='./.', DP=0)
                    record_samples.append(
                            vcf.model._Call(record, sample_id, sample_data))
                record.samples = record_samples

                vcf_writer.write_record(record)


def _update_info_string_number(vcf_reader, key, new_number):
    """Sets new number on info strings.
    """
    orig = vcf_reader.infos[key]
    vcf_reader.infos[key] = vcf.parser._Info(
            orig.id, -1, orig.type, orig.desc)
