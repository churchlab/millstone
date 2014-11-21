"""Common functions for variant calling.
"""

import os

import vcf

from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.model_utils import get_dataset_with_type


def common_postprocess_vcf(vcf_reader):
    # Do postprocessing in common to Pindel and Delly VCFs.
    modified_header_lines = []

    # These properties should be part of VA_DATA, although SV tools will
    #   output at most one property in each row and set Number=1 in the
    #   VCF header line. The easiest way to store these properties in
    #   VA_DATA is just to postprocess the VCF here and change these
    #   header lines to all say Number=A.
    va_properties = ['SVTYPE', 'SVLEN']

    def modify_header(header_line):
        if any([prop in header_line for prop in va_properties]):
            header_line = header_line.replace('Number=1', 'Number=A')
            modified_header_lines.append(header_line)
        return header_line
    vcf_reader._header_lines = map(modify_header, vcf_reader._header_lines)

    # Also add a field for METHOD.
    method_header_line = '##INFO=<ID=METHOD,Number=1,Type=String,' + \
        'Description="Type of approach used to detect SV">\n'
    modified_header_lines.append(method_header_line)
    vcf_reader._header_lines.append(method_header_line)

    # Now update the header lines in vcf_reader.infos map as well.
    parser = vcf.parser._vcf_metadata_parser()
    for header_line in modified_header_lines:
        key, val = parser.read_info(header_line)
        vcf_reader.infos[key] = val



# Returns a dictionary of common parameters required for all variant callers
# (freebayes, pindel, delly, lumpy).
def get_common_tool_params(alignment_group):
    alignment_type = Dataset.TYPE.BWA_ALIGN
    return {
        'alignment_group': alignment_group,
        'alignment_type': alignment_type,
        'fasta_ref': _get_fasta_ref(alignment_group),
        'output_dir': _create_output_dir(alignment_group),
        'sample_alignments': _find_valid_sample_alignments(
                alignment_group, alignment_type),
    }


def _get_fasta_ref(alignment_group):
    # Grab the reference genome fasta for the alignment.
    return get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()


def _create_output_dir(alignment_group):
    # Prepare a directory to put the output files.
    # We'll put them in
    #     /projects/<project_uid>/alignment_groups/vcf/<variant tool>/
    #     <alignment_type>.vcf
    # We'll save these for now, maybe it's not necessary later.
    vcf_dir = os.path.join(alignment_group.get_model_data_dir(), 'vcf')
    ensure_exists_0775_dir(vcf_dir)
    return vcf_dir


def _find_valid_sample_alignments(alignment_group, alignment_type):
    """ Returns a list sample alignment objects for an alignment,
        skipping those that failed. """
    sample_alignment_list = (
            alignment_group.experimentsampletoalignment_set.all())

    # Filter out mis-aligned files.
    # TODO: Should we show in the UI that some alignments failed and are
    # being skipped?
    def _is_successful_alignment(sample_alignment):
        bam_dataset = get_dataset_with_type(sample_alignment, alignment_type)
        return bam_dataset.status == Dataset.STATUS.READY
    sample_alignment_list = [sample_alignment for sample_alignment in
            sample_alignment_list if _is_successful_alignment(sample_alignment)]

    if len(sample_alignment_list) == 0:
        raise Exception('No successful alignments, Freebayes cannot proceed.')

    bam_files = [
            get_dataset_with_type(sa, alignment_type).get_absolute_location()
            for sa in sample_alignment_list]

    # Keep only valid bam_files
    valid_bam_files = []
    for bam_file in bam_files:
        if bam_file is None:
            continue
        if not os.stat(bam_file).st_size > 0:
            continue
        valid_bam_files.append(bam_file)
    assert len(valid_bam_files) == len(sample_alignment_list), (
            "Expected %d bam files, but found %d" % (
                    len(sample_alignment_list), len(bam_files)))
    return sample_alignment_list
