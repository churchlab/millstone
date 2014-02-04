"""
Functions for calling SNPs.
"""

import os
import subprocess

from celery import task

from main.models import clean_filesystem_location
from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.models import get_dataset_with_type
from main.s3 import project_files_needed
from variant_effects import run_snpeff
from scripts.vcf_parser import parse_alignment_group_vcf
from settings import TOOLS_DIR

# TODO: These VCF types should be set somewhere else. snpeff_util and vcf_parser
# also use them, but where should they go? settings.py seems logical, but it
# cannot import from models.py... -dbg

# Dataset type to use for snp calling.
VCF_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES
# Dataset type to use for snp annotation.
VCF_ANNOTATED_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES_SNPEFF


@task
@project_files_needed
def call_snvs(alignment_group):
    """Calls SNVs for all of the alignments in the alignment_group.
    """

    run_freebayes(alignment_group, Dataset.TYPE.BWA_ALIGN)

    # For now, automatically run snpeff if a genbank annotation is available.
    # If no annotation, then skip it, and pass the unannotated vcf type.
    if alignment_group.reference_genome.is_annotated():
        run_snpeff(alignment_group, Dataset.TYPE.BWA_ALIGN)
        vcf_dataset_type = VCF_ANNOTATED_DATASET_TYPE
    else:
        vcf_dataset_type = VCF_DATASET_TYPE

    # Parse the resulting vcf.
    parse_alignment_group_vcf(alignment_group, vcf_dataset_type)


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

    # Get handles for each of the bam files.
    def _get_bam_location(sample_alignment):
        bam_dataset = get_dataset_with_type(sample_alignment, alignment_type)
        return bam_dataset.get_absolute_location()
    bam_files = [_get_bam_location(sample_alignment) for sample_alignment in
            sample_alignment_list if sample_alignment]

    # Keep only valid bam_files
    valid_bam_files = []
    for bam_file in bam_files:
        if bam_file is None:
            continue
        if not os.stat(bam_file).st_size > 0:
            continue
        valid_bam_files.append(bam_file)
    bam_files = valid_bam_files
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
        subprocess.check_call(full_command, stdout=fh)

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
