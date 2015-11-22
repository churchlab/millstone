"""Wrapper for running pindel.
"""

import os
import shutil
import subprocess

from django.conf import settings
import vcf

from main.model_utils import get_dataset_with_type
from pipeline.read_alignment import get_insert_size_mean_and_stdev
from pipeline.variant_calling.common import common_postprocess_vcf


def run_pindel(fasta_ref, sample_alignments, vcf_output_dir,
        vcf_output_filename, alignment_type, **kwargs):
    """Run pindel to find SVs."""
    if not os.path.isdir('%s/pindel' % settings.TOOLS_DIR):
        raise Exception('Pindel is not installed. Aborting.')

    bam_files = [
            get_dataset_with_type(sa, alignment_type).get_absolute_location()
            for sa in sample_alignments]

    samples = [sa.experiment_sample for sa in sample_alignments]
    insert_sizes = [get_insert_size_mean_and_stdev(sa) for sa in
            sample_alignments]

    assert len(bam_files) == len(insert_sizes)

    # Create pindel config file
    pindel_config = os.path.join(vcf_output_dir, 'pindel_config.txt')
    at_least_one_config_line_written = False
    with open(pindel_config, 'w') as fh:
        for bam_file, sample, insert_size in zip(
                bam_files, samples, insert_sizes):

            # Skip bad alignments.
            mean, stdev = insert_size
            if mean == -1:
                continue
            fh.write('%s %s %s\n' % (bam_file, mean, sample.uid))
            at_least_one_config_line_written = True

    if not at_least_one_config_line_written:
        raise Exception
        return False # failure

    # Build the full pindel command.
    pindel_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    subprocess.check_call(['%s/pindel/pindel' % settings.TOOLS_DIR,
        '-f', fasta_ref,
        '-i', pindel_config,
        '-c', 'ALL',
        '-o', pindel_root
    ])

    # convert all different structural variant types to vcf
    subprocess.check_call(['%s/pindel/pindel2vcf' % settings.TOOLS_DIR,
        '-P', pindel_root,
        '-r', fasta_ref,
        '-R', 'name',
        '-d', 'date',
        '-mc', '1',  # just need one read to show 1/1 in vcf
    ])

    postprocess_pindel_vcf(vcf_output_filename)

    return True # success


def postprocess_pindel_vcf(vcf_file):
    """Process vcfs output by Pindel, so that the information is
    customized to whatever is needed in Millstone, and the format is
    the same as that of Freebayes.

    Args:
        vcf_file: This is typically the output of pindel2vcf. This file is
            over-written by this function.
    """
    # First create a backup file.
    shutil.copyfile(vcf_file, vcf_file + '.prepostprocess')

    temp_vcf_filename = vcf_file + '.tmp'

    with open(vcf_file) as input_fh:
        # Use the original file as the template.
        vcf_reader = vcf.Reader(input_fh)

        # Update some headers stuff.
        common_postprocess_vcf(vcf_reader)

        with open(temp_vcf_filename, 'w') as output_fh:
            vcf_writer = vcf.Writer(output_fh, vcf_reader)
            for record in vcf_reader:
                assert 'SVLEN' in record.__dict__['INFO']

                # pindel uses negative SVLEN for deletions; make them positive
                # always have one entry
                svlen = abs(record.__dict__['INFO']['SVLEN'][0])
                record.__dict__['INFO']['SVLEN'] = [svlen]

                # Filter out small variants. These should be handled by Freebayes.
                if svlen < 10:
                    continue

                # Update METHOD field.
                record.__dict__['INFO']['METHOD'] = 'PINDEL'

                # If REF is super long, truncate to first char.
                if len(record.REF) > 10:
                    record.REF = record.REF[0]

                # Set ALT to SVTYPE.
                record.ALT = ['<{svtype}>'.format(svtype=record.__dict__['INFO']['SVTYPE'][0])]

                vcf_writer.write_record(record)

    # Replace original filepath with modified temp file.
    shutil.move(temp_vcf_filename, vcf_file)
