"""
Functions for calling SNPs.
"""

from celery import group
from celery import chain
from celery import task
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
from read_alignment import get_insert_size
from main.s3 import project_files_needed
from variant_effects import run_snpeff
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
VCF_PINDEL_TYPE = Dataset.TYPE.VCF_PINDEL
VCF_DELLY_TYPE = Dataset.TYPE.VCF_DELLY

# Delete this function eventually
# It is just a special case of find_variants (only runs freebayes), kept here so that test_snv_calling still works.
@task
@project_files_needed
def call_snvs(alignment_group):
    alignment_type = Dataset.TYPE.BWA_ALIGN
    common_params = {
            'alignment_group': alignment_group,
            'alignment_type': alignment_type,
            'fasta_ref': _get_fasta_ref(alignment_group),
            'output_dir': _create_output_dir(alignment_group),
            'bam_files': _find_valid_bam_files(alignment_group, alignment_type),
            }

    variant_tools = (
            ('freebayes', Dataset.TYPE.VCF_FREEBAYES, run_freebayes),
            )

    for variant_params in variant_tools:
        find_variants_with_tool(common_params, variant_params)

@task
@project_files_needed
def find_variants(alignment_group):
    """Calls SNVs for all of the alignments in the alignment_group.
    """
    alignment_type = Dataset.TYPE.BWA_ALIGN
    common_params = {
            'alignment_group': alignment_group,
            'alignment_type': alignment_type,
            'fasta_ref': _get_fasta_ref(alignment_group),
            'output_dir': _create_output_dir(alignment_group),
            'bam_files': _find_valid_bam_files(alignment_group, alignment_type),
            }

    variant_tools = (
            ('freebayes', Dataset.TYPE.VCF_FREEBAYES, run_freebayes),
            ('pindel', Dataset.TYPE.VCF_PINDEL, run_pindel),
            ('delly', Dataset.TYPE.VCF_DELLY, run_delly),
            )

    for variant_params in variant_tools:
        find_variants_with_tool(common_params, variant_params)

def _get_fasta_ref(alignment_group):
    # Grab the reference genome fasta for the alignment.
    return get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

def _create_output_dir(alignment_group):
    # Prepare a directory to put the output files.
    # We'll put them in /projects/<project_uid>/alignment_groups/vcf/<variant tool>/
    #     <alignment_type>.vcf
    # We'll save these for now, maybe it's not necessary later.
    vcf_dir = os.path.join(alignment_group.get_model_data_dir(), 'vcf')
    ensure_exists_0775_dir(vcf_dir)
    return vcf_dir

def _find_valid_bam_files(alignment_group, alignment_type):
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
    return bam_files

def find_variants_with_tool(common_params, variant_params):
    alignment_group = common_params['alignment_group']
    tool_name, vcf_dataset_type, tool_function = variant_params

    # Create subdirectory for this tool
    tool_dir = os.path.join(common_params['output_dir'], tool_name)
    ensure_exists_0775_dir(tool_dir)
    vcf_output_filename = os.path.join(tool_dir, common_params['alignment_type'] + '.vcf')

    # Run the tool
    tool_function(common_params['fasta_ref'],
            common_params['bam_files'], tool_dir, vcf_output_filename)

    # Add dataset
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

    # Do the following only for freebayes; right now just special if condition
    if tool_name == 'freebayes':
        # For now, automatically run snpeff if a genbank annotation is available.
        # If no annotation, then skip it, and pass the unannotated vcf type.
        if alignment_group.reference_genome.is_annotated():
            run_snpeff(alignment_group, Dataset.TYPE.BWA_ALIGN)
            vcf_dataset_type = VCF_ANNOTATED_DATASET_TYPE
        else:
            vcf_dataset_type = VCF_DATASET_TYPE

    # Parse the resulting vcf
    parse_alignment_group_vcf(alignment_group, vcf_dataset_type)


def run_freebayes(fasta_ref, bam_files, vcf_output_dir, vcf_output_filename):
    """Run freebayes using the bam alignment files keyed by the alignment_type
    for all Genomes of the passed in ReferenceGenome.

    NOTE: If a Genome doesn't have a bam alignment file with this
    alignment_type, then it won't be used.
    """
    vcf_dataset_type = VCF_DATASET_TYPE

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


def _get_sample_uid(bam_file):
    """extract sample uid information from bam file"""
    # convert to readable sam file
    process = subprocess.Popen(['%s/samtools/samtools' % TOOLS_DIR, 'view', bam_file],
            stdout=subprocess.PIPE)
    # read only the first 10 lines, to avoid loading large bam files into memory
    samdata = subprocess.check_output(['head'], stdin=process.stdout)
    # find the sample uid, which is found in the form RG:Z:[uid]
    match = re.search('RG:Z:(\S+)', samdata)
    if not match:
        # should not happen if the bam file is structured properly
        assert 'Internal error: no sample uid found in bam file'
    return match.group(1)


def _filter_small_variants(vcf_file, cutoff):
    """Go through each line of vcf, and remove small structural variants"""
    vcf_file_tmp = vcf_file + '.tmp'
    with open(vcf_file_tmp, 'w') as fh:
        for line in open(vcf_file):
            match = re.search('SVLEN=(-?[0-9]+);', line)
            # Check if SVLEN > cutoff
            if not match or abs(int(match.group(1))) > cutoff:
                fh.write(line)

    # move temporary file back to vcf_file path
    subprocess.check_call(['mv', vcf_file_tmp, vcf_file])


def run_pindel(fasta_ref, bam_files, vcf_output_dir, vcf_output_filename):
    """Run pindel to find SVs."""
    vcf_dataset_type = VCF_PINDEL_TYPE

    # Create pindel config file
    pindel_config = os.path.join(vcf_output_dir, 'pindel_config.txt')
    with open(pindel_config, 'w') as fh:
        for bam_file in bam_files:
            insert_size = get_insert_size(bam_file)
            sample_uid = _get_sample_uid(bam_file)
            fh.write('%s %s %s\n' % (bam_file, insert_size, sample_uid))

    # Build the full pindel command.
    pindel_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    subprocess.check_call(['%s/pindel/pindel' % TOOLS_DIR,
        '-f', fasta_ref,
        '-i', pindel_config,
        '-c', 'ALL',
        '-o', pindel_root
    ])

    # convert all different structural variant types to vcf
    subprocess.check_call(['%s/pindel/pindel2vcf' % TOOLS_DIR,
        '-P', pindel_root,
        '-r', fasta_ref,
        '-R', 'name',
        '-d', 'date'
    ])
    _filter_small_variants(vcf_output_filename, 10)


def run_delly(fasta_ref, bam_files, vcf_output_dir, vcf_output_filename):
    """Run delly to find SVs."""
    vcf_dataset_type = VCF_DELLY_TYPE

    delly_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    transformations = ['DEL', 'DUP', 'INV']
    vcf_outputs = map(lambda transformation:
            '%s_%s.vcf' % (delly_root, transformation), transformations)

    # run delly for each type of transformation
    for transformation, vcf_output in zip(transformations, vcf_outputs):
        # not checked_call, because delly errors if it doesn't find any SVs
        subprocess.call(['%s/delly/delly' % TOOLS_DIR,
            '-t', transformation,
            '-o', vcf_output,
            '-g', fasta_ref] + bam_files)

    # combine the separate vcfs for each transformation
    # TODO update sample uids
    vcf_outputs = filter(lambda file: os.path.exists(file), vcf_outputs)
    if vcf_outputs:
        with open(vcf_output_filename, 'w') as fh:
            subprocess.check_call(['vcf-concat'] + vcf_outputs, stdout=fh)
    else:
        # hack: create empty vcf
        subprocess.check_call(['touch', delly_root])
        subprocess.check_call(['%s/pindel/pindel2vcf' % TOOLS_DIR,
            '-p', delly_root,  # TODO does this work?
            '-r', fasta_ref,
            '-R', 'name',
            '-d', 'date'
        ])

