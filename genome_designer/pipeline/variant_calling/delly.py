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
    assert os.path.exists(settings.DELLY_BIN), (
            'Delly is not installed. Aborting.')

    delly_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    transformations = ['DEL', 'DUP', 'INV']
    vcf_outputs = map(lambda transformation:
            '%s_%s.vcf' % (delly_root, transformation), transformations)

    # Create symlinks to bam files which use uid because Delly uses the name of
    # the file as sample uid in the output report.
    new_bam_files = []
    bam_files = [
            get_dataset_with_type(sa, alignment_type).get_absolute_location()
            for sa in sample_alignments]
    samples = [sa.experiment_sample for sa in sample_alignments]
    for bam_file, sample in zip(bam_files, samples):
        new_bam_file = os.path.join(
                os.path.dirname(bam_file), sample.uid + '.bam')
        _clean_symlink(bam_file, new_bam_file)
        _clean_symlink(bam_file + '.bai', new_bam_file + '.bai')
        new_bam_files.append(new_bam_file)

    # run delly for each type of transformation
    for transformation, vcf_output in zip(transformations, vcf_outputs):

        # not checked_call, because delly errors if it doesn't find any SVs
        subprocess.call([
            settings.DELLY_BIN,
            '-t', transformation,
            '-o', vcf_output,
            '-g', fasta_ref] + new_bam_files)

    # combine the separate vcfs for each transformation
    vcf_outputs = [f for f in vcf_outputs if os.path.exists(f)]
    if vcf_outputs:
        temp_vcf = os.path.join(vcf_output_dir, 'temp_vcf')
        os.putenv('PERL5LIB', os.path.join(settings.VCFTOOLS_DIR, 'perl'))
        with open(temp_vcf, 'w') as fh:
            subprocess.check_call([settings.VCF_CONCAT_BINARY] + vcf_outputs,
                    stdout=fh)
        with open(vcf_output_filename, 'w') as fh:
            subprocess.check_call([settings.VCF_SORT_BINARY, temp_vcf],
                    stdout=fh)
        os.remove(temp_vcf)
    else:
        # hack: create empty vcf
        subprocess.check_call(['touch', delly_root])
        subprocess.check_call(['%s/pindel/pindel2vcf' % settings.TOOLS_DIR,
            '-p', delly_root,  # TODO does this work?
            '-r', fasta_ref,
            '-R', 'name',
            '-d', 'date'])

    # Delete temporary bam file symlinks.
    for f in new_bam_files:
        os.remove(f)
        os.remove(f + '.bai')

    postprocess_delly_vcf(vcf_output_filename)

    return True # success


def _clean_symlink(src, dest):
    """Creates symlink, deleting dest if it already exists.
    """
    if os.path.exists(dest):
        os.remove(dest)
    os.symlink(src, dest)


def postprocess_delly_vcf(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file))

    common_postprocess_vcf(vcf_reader)

    vcf_writer = vcf.Writer(open(vcf_file + '.tmp', 'a'), vcf_reader)
    for record in vcf_reader:
        record.__dict__['INFO']['METHOD'] = 'DELLY'
        vcf_writer.write_record(record)

    subprocess.check_call(['mv', vcf_file + '.tmp', vcf_file])
