"""Wrapper for delly.
"""

import os
import subprocess

from django.conf import settings
import vcf

from main.model_utils import get_dataset_with_type
from pipeline.variant_calling.common import common_postprocess_vcf


def run_delly(fasta_ref, sample_alignments, vcf_output_dir,
        vcf_output_filename, alignment_type, **kwargs):
    """Run delly to find SVs."""
    if not os.path.isdir('%s/delly' % settings.TOOLS_DIR):
        raise Exception('Delly is not installed. Aborting.')

    delly_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    transformations = ['DEL', 'DUP', 'INV']
    vcf_outputs = map(lambda transformation:
            '%s_%s.vcf' % (delly_root, transformation), transformations)

    # Rename bam files, because Delly uses the name of the file as sample uid.
    # Use cp instead of mv, because other sv callers will be reading from the
    #   original bam file.
    bam_files = [
            get_dataset_with_type(sa, alignment_type).get_absolute_location()
            for sa in sample_alignments]
    samples = [sa.experiment_sample for sa in sample_alignments]

    new_bam_files = []
    for bam_file, sample in zip(bam_files, samples):
        new_bam_file = os.path.join(
                os.path.dirname(bam_file), sample.uid + '.bam')
        subprocess.check_call(['cp', bam_file, new_bam_file])
        subprocess.check_call(['cp', bam_file + '.bai', new_bam_file + '.bai'])
        new_bam_files.append(new_bam_file)

    # run delly for each type of transformation
    for transformation, vcf_output in zip(transformations, vcf_outputs):

        # not checked_call, because delly errors if it doesn't find any SVs
        subprocess.call(['%s/delly/delly' % settings.TOOLS_DIR,
            '-t', transformation,
            '-o', vcf_output,
            '-g', fasta_ref] + new_bam_files)

    # combine the separate vcfs for each transformation
    vcf_outputs = filter(lambda file: os.path.exists(file), vcf_outputs)
    if vcf_outputs:
        temp_vcf = os.path.join(vcf_output_dir, 'temp_vcf')
        with open(temp_vcf, 'w') as fh:
            subprocess.check_call(['vcf-concat'] + vcf_outputs, stdout=fh)
        with open(vcf_output_filename, 'w') as fh:
            subprocess.check_call(['vcf-sort', temp_vcf], stdout=fh)
        subprocess.check_call(['rm', temp_vcf])
    else:
        # hack: create empty vcf
        subprocess.check_call(['touch', delly_root])
        subprocess.check_call(['%s/pindel/pindel2vcf' % settings.TOOLS_DIR,
            '-p', delly_root,  # TODO does this work?
            '-r', fasta_ref,
            '-R', 'name',
            '-d', 'date'])

    # Delete temporary renamed bam files
    for bam_file in new_bam_files:
        subprocess.check_call(['rm', bam_file])
        subprocess.check_call(['rm', bam_file + '.bai'])

    postprocess_delly_vcf(vcf_output_filename)

    return True # success


def postprocess_delly_vcf(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file))

    common_postprocess_vcf(vcf_reader)

    vcf_writer = vcf.Writer(open(vcf_file + '.tmp', 'a'), vcf_reader)
    for record in vcf_reader:
        record.__dict__['INFO']['METHOD'] = 'DELLY'
        vcf_writer.write_record(record)

    subprocess.check_call(['mv', vcf_file + '.tmp', vcf_file])
