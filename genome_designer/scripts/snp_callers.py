"""
Functions for calling SNPs.
"""

import os
import re
import subprocess

import vcf

from main.celery_util import CELERY_ERROR_KEY
from main.celery_util import get_celery_worker_status
from main.models import clean_filesystem_location
from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.models import get_dataset_with_type
from scripts.alignment_pipeline import get_insert_size
from scripts.snpeff_util import run_snpeff
from scripts.util import fn_runner
from scripts.vcf_parser import parse_alignment_group_vcf
from settings import DEBUG_CONCURRENT
from settings import PWD
from settings import TOOLS_DIR

# TODO: These VCF types should be set somewhere else. snpeff_util and vcf_parser
# also use them, but where should they go? settings.py seems logical, but it
# cannot import from models.py... -dbg

# Dataset type to use for snp calling.
VCF_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES
# Dataset type to use for snp annotation.
VCF_ANNOTATED_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES_SNPEFF
# Dataset type for results of finding SVs.
VCF_SV_DATASET_TYPE = Dataset.TYPE.VCF_SV

def run_snp_calling_pipeline(alignment_group, concurrent=DEBUG_CONCURRENT):
    """Calls SNPs for all of the alignments in the alignment_group.
    """
    # Check whether Celery is running.
    if concurrent:
        celery_status = get_celery_worker_status()
        assert not CELERY_ERROR_KEY in celery_status, (
                celery_status[CELERY_ERROR_KEY])

    args = [alignment_group]
    fn_runner(run_snp_calling_pipeline_internal, alignment_group.reference_genome.project, args, concurrent)


def run_snp_calling_pipeline_internal(alignment_group):
    """Internal method to provide async interface.
    """
    run_analysis_pipeline(alignment_group, Dataset.TYPE.BWA_ALIGN)

    # For now, automatically run snpeff if a genbank annotation is available.
    # If no annotation, then skip it, and pass the unannotated vcf type.
    if alignment_group.reference_genome.is_annotated():
        run_snpeff(alignment_group, Dataset.TYPE.BWA_ALIGN)
        vcf_dataset_type = VCF_ANNOTATED_DATASET_TYPE
    else:
        vcf_dataset_type = VCF_DATASET_TYPE
    vcf_sv_dataset_type = VCF_SV_DATASET_TYPE

    # Parse the resulting vcfs.
    parse_alignment_group_vcf(alignment_group, vcf_dataset_type)
    parse_alignment_group_vcf(alignment_group, vcf_sv_dataset_type)


def run_analysis_pipeline(alignment_group, alignment_type):
    """ Get the input/output information for running freebayes
    and pindel/delly
    """
    # Grab the reference genome fasta for the alignment.
    fasta_ref = get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

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
    sv_vcf_dir = os.path.join(vcf_dir, 'sv')
    ensure_exists_0775_dir(sv_vcf_dir)
    vcf_sv_output_filename = os.path.join(
            sv_vcf_dir, alignment_type + '.vcf')

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

    run_freebayes(alignment_group, fasta_ref, bam_files, VCF_DATASET_TYPE, vcf_output_filename)
    run_sv(alignment_group, fasta_ref, bam_files, VCF_SV_DATASET_TYPE, sv_vcf_dir)

def _add_dataset(alignment_group, vcf_dataset_type, vcf_output_filename):
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


def run_freebayes(alignment_group, fasta_ref, bam_files, vcf_dataset_type, vcf_output_filename):
    """Run freebayes using the bam alignment files keyed by the alignment_type
    for all Genomes of the passed in ReferenceGenome.

    NOTE: If a Genome doesn't have a bam alignment file with this
    alignment_type, then it won't be used.
    """
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

    # add dataset for freebayes.
    _add_dataset(alignment_group, vcf_dataset_type, vcf_output_filename)


def _get_sample_uid(bam_file):
    """extract sample uid information from bam file"""
    # convert to readable sam file
    process = subprocess.Popen(['%s/samtools/samtools' % TOOLS_DIR, 'view', bam_file],
            stdout=subprocess.PIPE)
    # read only the first 10 lines, to avoid loading large bam files into memory
    samdata = subprocess.check_output(['head'], stdin=process.stdout)
    print 'samdata:', samdata
    # find the sample uid, which is found in the form RG:Z:[uid]
    match = re.search('RG:Z:(\S+)', samdata)
    print 'match:', match
    if not match:
        # should not happen if the bam file is structured properly
        assert 'Internal error: no sample uid found in bam file'
    return match.group(1)


def run_sv(alignment_group, fasta_ref, bam_files, vcf_dataset_type, sv_vcf_dir):
    """Run pindel and delly to find SVs."""
    # Create pindel config file
    pindel_config = os.path.join(sv_vcf_dir, 'pindel_config.txt')
    with open(pindel_config, 'w') as fh:
        for bam_file in bam_files:
            insert_size = get_insert_size(bam_file)
            sample_uid = _get_sample_uid(bam_file)
            fh.write('%s %s %s\n' % (bam_file, insert_size, sample_uid))

    # Build the full pindel command.
    print fasta_ref, pindel_config, sv_vcf_dir
    subprocess.check_call(['%s/pindel/pindel' % TOOLS_DIR,
        '-f', fasta_ref,
        '-i', pindel_config,
        '-c', 'ALL',
        '-o', os.path.join(sv_vcf_dir, 'pindel')
    ])

    # convert different types to vcf separately
    pindel_root = os.path.join(sv_vcf_dir, 'pindel')
    subprocess.check_call(['%s/pindel/pindel2vcf' % TOOLS_DIR,
        '-P', pindel_root,  # -P pindel output root
        '-r', fasta_ref,
        '-R', 'name',
        '-d', 'date'
    ])

    # add dataset for sv.
    _add_dataset(alignment_group, vcf_dataset_type, pindel_root + '.vcf')

