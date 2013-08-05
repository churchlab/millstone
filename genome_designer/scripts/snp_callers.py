"""
Functions for calling SNPs.
"""

import os
import subprocess

import vcf

from main.models import clean_filesystem_location
from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.models import get_dataset_with_type
from scripts.util import fn_runner
from scripts.vcf_parser import parse_alignment_group_vcf
from settings import DEBUG_CONCURRENT
from settings import PWD, TOOLS_DIR

# Load the tools dir from settings.py
TOOLS_DIR = os.path.join(PWD,TOOLS_DIR)

# For now, we always use this dataset type for storing the vcf.
VCF_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES


def run_snp_calling_pipeline(alignment_group, concurrent=DEBUG_CONCURRENT):
    """Calls SNPs for all of the alignments in the alignment_group.
    """
    args = [alignment_group]
    fn_runner(run_snp_calling_pipeline_internal, args, concurrent)


def run_snp_calling_pipeline_internal(alignment_group):
    """Internal method to provide async interface.
    """
    run_freebayes(alignment_group, Dataset.TYPE.BWA_ALIGN)

    # Parse the resulting vcf.
    parse_alignment_group_vcf(alignment_group, VCF_DATASET_TYPE)


def run_freebayes(alignment_group, alignment_type):
    """Run freebayes using the bam alignment files keyed by the alignment_type
    for all Genomes of the passed in ReferenceGenome.

    NOTE: If a Genome doesn't have a bam alignment file with this
    alignment_type, then it won't be used.
    """
    # Grab the reference genome fasta for the alignment.
    fasta_ref = get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # We'll store it as a Dataset on the Genome,
    # implicitly validating the alignment_type.
    vcf_dataset_type = VCF_DATASET_TYPE

    # Prepare a directory to put the output files.
    # We'll put them in /projects/<project_uid>/alignment_groups/vcf/freebayes/
    #     <alignment_type>.vcf
    # We'll save these for now, maybe it's not necessary later.
    vcf_dir = os.path.join(alignment_group.get_model_data_dir(), 'vcf')
    ensure_exists_0775_dir(vcf_dir)
    freebayes_vcf_dir = os.path.join(vcf_dir, 'freebayes')
    ensure_exists_0775_dir(freebayes_vcf_dir)
    vcf_output_filename = os.path.join(
            freebayes_vcf_dir, alignment_type + '.vcf')

    # Get handles for each of the bam files.
    sample_alignment_list = (
            alignment_group.experimentsampletoalignment_set.all())
    bam_files = map(
            lambda sample_alignment: get_dataset_with_type(
                    sample_alignment, alignment_type).get_absolute_location(),
            sample_alignment_list
    )

    # Keep only valid bam_files
    bam_files = filter(lambda bam_file: bam_file is not None, bam_files)
    assert len(bam_files) == len(sample_alignment_list), (
            "Expected %d bam files, but found %d" % (
                    len(sample_alignment_list), len(bam_files)))

    # Build up the bam part of the freebayes binary call.
    bam_part = []
    for bam_file in bam_files:
        bam_part.append('--bam')
        bam_part.append(bam_file)

    # Build the full command and execute it for all bam files at once.
    full_command = (['%s/freebayes/freebayes' %  TOOLS_DIR] + bam_part + [
        '--fasta-reference', fasta_ref,
        '--pvar', '0.001',
        '--ploidy', '2',
        '--min-alternate-fraction', '.3',
        '--hwe-priors-off',
        '--binomial-obs-priors-off',
        '--use-mapping-quality',
        '--min-base-quality', '25',
        '--min-mapping-quality', '30'
    ])

    with open(vcf_output_filename, 'w') as fh:
        subprocess.call(full_command, stdout=fh)

    # If a Dataset already exists, delete it, might have been a bad run.
    existing_set = Dataset.objects.filter(
            type=vcf_dataset_type,
            label=vcf_dataset_type,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    if len(existing_set) > 0:
        existing_set[0].delete()

    dataset = Dataset.objects.create(
            type=vcf_dataset_type,
            label=vcf_dataset_type,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    alignment_group.dataset_set.add(dataset)
