"""
Functions for calling Variants.
"""

import os

from celery import task

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.s3 import project_files_needed
from pipeline.variant_effects import run_snpeff
from pipeline.variant_calling.common import add_vcf_dataset
from pipeline.variant_calling.common import get_common_tool_params
from pipeline.variant_calling.common import process_vcf_dataset
from pipeline.variant_calling.constants import TOOL_DELLY
from pipeline.variant_calling.constants import TOOL_FREEBAYES
from pipeline.variant_calling.constants import TOOL_LUMPY
from pipeline.variant_calling.constants import TOOL_PINDEL
from pipeline.variant_calling.delly import run_delly
from pipeline.variant_calling.freebayes import run_freebayes
from pipeline.variant_calling.lumpy import run_lumpy
from pipeline.variant_calling.pindel import run_pindel
from utils import uppercase_underscore


###############################################################################
# Constants and Maps
###############################################################################

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


###############################################################################
# Tasks
###############################################################################

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
            * tool_kwargs (optional, passed on top of common_params to runner_fn)

    Returns:
        Boolean indicating whether we made it through this entire function.
    """
    # TODO: More informative failure information.
    try:
        common_params = get_common_tool_params(alignment_group)
    except Exception as e:
        alignment_group.status = AlignmentGroup.STATUS.FAILED
        alignment_group.save(update_fields=['status'])
        raise e

    tool_name = variant_params_dict['tool_name']
    vcf_dataset_type = variant_params_dict['dataset_type']
    tool_function = variant_params_dict['runner_fn']
    tool_kwargs = variant_params_dict.get('tool_kwargs', {})

    is_parallel_tool = 'region_num' in tool_kwargs

    # Finding variants means that all the aligning is complete, so now we
    # are VARIANT_CALLING.
    alignment_group.status = AlignmentGroup.STATUS.VARIANT_CALLING
    alignment_group.save(update_fields=['status'])

    # Create subdirectory for this tool
    tool_dir = os.path.join(common_params['output_dir'], tool_name)
    ensure_exists_0775_dir(tool_dir)

    # Make vcf output filename
    if is_parallel_tool:
        vcf_output_filename = os.path.join(tool_dir,
                uppercase_underscore(common_params['alignment_type']) +
                '.partial.' + str(tool_kwargs['region_num']) + '.vcf')
    else:
        vcf_output_filename = os.path.join(tool_dir,
                uppercase_underscore(common_params['alignment_type']) +
                '.vcf')

    # Run the tool
    common_params.update(tool_kwargs)
    try:
        tool_succeeded = tool_function(
                vcf_output_dir=tool_dir,
                vcf_output_filename=vcf_output_filename,
                **common_params)
    except Exception as exc:
        alignment_group = AlignmentGroup.objects.get(id=alignment_group.id)
        alignment_group.status = AlignmentGroup.STATUS.FAILED
        alignment_group.save(update_fields=['status'])
        print 'Tool failed: ' + str(exc)
        return False  # failed

    if not tool_succeeded:
        return False

    # Certain tools require merging and so are responsible for their own
    # vcf processing.
    does_vcf_processing_occur_elsewhere = (
            is_parallel_tool and tool_name in [TOOL_FREEBAYES, TOOL_LUMPY])
    if does_vcf_processing_occur_elsewhere:
        return True

    # Otherwise proceed with processing vcf dataset.

    # add unannotated vcf dataset first
    add_vcf_dataset(alignment_group, vcf_dataset_type, vcf_output_filename)

    # If freebayes is not parallelized, then run snpeff now,
    # then update the vcf_output_filename and vcf_dataset_type.
    if (tool_name == TOOL_FREEBAYES and
            alignment_group.reference_genome.is_annotated()):
        vcf_output_filename = run_snpeff(alignment_group, Dataset.TYPE.BWA_ALIGN)
        vcf_dataset_type = Dataset.TYPE.VCF_FREEBAYES_SNPEFF
        add_vcf_dataset(alignment_group, vcf_dataset_type, vcf_output_filename)

    # Finally, generate variants on the final vcf.
    process_vcf_dataset(alignment_group, vcf_dataset_type)

    return True  # success
