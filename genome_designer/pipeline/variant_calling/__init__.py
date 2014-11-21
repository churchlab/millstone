"""
Functions for calling Variants.
"""

import os
import shutil
import subprocess
from uuid import uuid4

from celery import task

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.model_utils import clean_filesystem_location
from main.model_utils import get_dataset_with_type
from main.s3 import project_files_needed
from pipeline.variant_effects import run_snpeff
from pipeline.variant_calling.common import get_common_tool_params
from pipeline.variant_calling.delly import run_delly
from pipeline.variant_calling.freebayes import run_freebayes
from pipeline.variant_calling.lumpy import run_lumpy
from pipeline.variant_calling.pindel import run_pindel
from utils.jbrowse_util import add_vcf_track
from utils import uppercase_underscore
from variants.vcf_parser import parse_alignment_group_vcf
from variants.variant_sets import add_variants_to_set_from_bed

# TODO: These VCF types should be set somewhere else. snpeff_util and
# vcf_parser also use them, but where should they go? settings.py seems
# logical, but it cannot import from models.py... -dbg

# Dataset type to use for snp calling.
VCF_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES
# Dataset type to use for snp annotation.
VCF_ANNOTATED_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES_SNPEFF

# Dataset type for results of finding SVs.
VCF_PINDEL_TYPE = Dataset.TYPE.VCF_PINDEL
VCF_DELLY_TYPE = Dataset.TYPE.VCF_DELLY

TOOL_FREEBAYES = 'freebayes'
TOOL_PINDEL = 'pindel'
TOOL_DELLY = 'delly'
TOOL_LUMPY = 'lumpy'


# Map from variant caller name to params.
VARIANT_TOOL_PARAMS_MAP = {
    TOOL_FREEBAYES: {
        'tool_name': TOOL_FREEBAYES,
        'dataset_type': Dataset.TYPE.VCF_FREEBAYES,
        'runner_fn': run_freebayes
    },
    TOOL_PINDEL: {
        'tool_name': TOOL_PINDEL,
        'dataset_type': Dataset.TYPE.VCF_PINDEL,
        'runner_fn': run_pindel
    },
    TOOL_DELLY: {
        'tool_name': TOOL_DELLY,
        'dataset_type': Dataset.TYPE.VCF_DELLY,
        'runner_fn': run_delly
    },
    TOOL_LUMPY: {
        'tool_name': TOOL_LUMPY,
        'dataset_type': Dataset.TYPE.VCF_LUMPY,
        'runner_fn': run_lumpy
    }
}


@task
@project_files_needed
def find_variants_with_tool(alignment_group, variant_params_dict):
    """Applies a variant caller to the alignment data contained within
    alignment_group.

    Args:
        alignment_group: AlignmentGroup with all alignments complete.
        variant_params_dict: Dictionary with keys:
            * tool_name
            * dataset_type
            * runner_fn

    Returns:
        Boolean indicating whether we made it through this entire function.
    """
    # TODO: More informative failure information.
    try:
        common_params = get_common_tool_params(alignment_group)
    except:
        alignment_group.status = AlignmentGroup.STATUS.FAILED
        alignment_group.save(update_fields=['status'])
        return False

    tool_name = variant_params_dict['tool_name']
    vcf_dataset_type = variant_params_dict['dataset_type']
    tool_function = variant_params_dict['runner_fn']

    # Finding variants means that all the aligning is complete, so now we
    # are VARIANT_CALLING.
    alignment_group.status = AlignmentGroup.STATUS.VARIANT_CALLING
    alignment_group.save()

    # Create subdirectory for this tool
    tool_dir = os.path.join(common_params['output_dir'], tool_name)
    ensure_exists_0775_dir(tool_dir)
    vcf_output_filename = os.path.join(tool_dir,
            uppercase_underscore(common_params['alignment_type']) + '.vcf')

    # Run the tool
    tool_succeeded = tool_function(
            vcf_output_dir=tool_dir,
            vcf_output_filename=vcf_output_filename,
            **common_params)
    if not tool_succeeded:
        return False

    # Add dataset
    # If a Dataset already exists, delete it, might have been a bad run.
    existing_set = Dataset.objects.filter(
            type=vcf_dataset_type,
            label=vcf_dataset_type,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    if len(existing_set) > 0:
        existing_set[0].delete()

    vcf_dataset = Dataset.objects.create(
            type=vcf_dataset_type,
            label=vcf_dataset_type,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    alignment_group.dataset_set.add(vcf_dataset)

    # Do the following only for freebayes; right now just special if condition
    if tool_name == TOOL_FREEBAYES:
        # For now, automatically run snpeff if a genbank annotation is
        # available.
        # If no annotation, then skip it, and pass the unannotated vcf type.
        if alignment_group.reference_genome.is_annotated():
            run_snpeff(alignment_group, Dataset.TYPE.BWA_ALIGN)
            vcf_dataset_type = VCF_ANNOTATED_DATASET_TYPE
        else:
            vcf_dataset_type = VCF_DATASET_TYPE

    sort_vcf(vcf_dataset.get_absolute_location())

    # Tabix index and add the VCF track to Jbrowse
    add_vcf_track(alignment_group.reference_genome, alignment_group,
        vcf_dataset_type)

    # Parse the resulting vcf, grab variant objects
    parse_alignment_group_vcf(alignment_group, vcf_dataset_type)

    flag_variants_from_bed(alignment_group, Dataset.TYPE.BED_CALLABLE_LOCI)

    return True


def flag_variants_from_bed(alignment_group, bed_dataset_type):
    sample_alignments = alignment_group.experimentsampletoalignment_set.all()
    for sample_alignment in sample_alignments:

        # If there is no callable_loci bed, skip the sample alignment.
        # TODO: Make this extensible to other BED files we might have
        callable_loci_bed = get_dataset_with_type(
                entity=sample_alignment,
                type=Dataset.TYPE.BED_CALLABLE_LOCI)

        if not callable_loci_bed:
            continue

        # need to add sample_alignment and bed_dataset here.
        add_variants_to_set_from_bed(
                sample_alignment=sample_alignment,
                bed_dataset=callable_loci_bed)


def sort_vcf(input_vcf_filepath):
    """Sorts a vcf file by chromosome and position.

    Overwrites the input.
    """
    temp_vcf = os.path.splitext(input_vcf_filepath)[0] + str(uuid4())[:8]
    assert not os.path.exists(temp_vcf)

    sort_cmd = (
            '(grep ^"#" {original_vcf}; grep -v ^"#" {original_vcf} | '
            'sort -k1,1 -k2,2n) > {sorted_vcf}'
    ).format(
            original_vcf=input_vcf_filepath,
            sorted_vcf=temp_vcf
    )
    subprocess.call(sort_cmd, shell=True)

    shutil.move(temp_vcf, input_vcf_filepath)
