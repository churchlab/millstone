"""
The alignment pipeline.

We start with the .fastq files and the reference for a particular genome
and carry out the gauntlet of steps to perform alignments, alignment
cleaning, snv calling, and effect prediction.
"""

from celery import chain
from celery import group
from django.conf import settings

from main.celery_util import CELERY_ERROR_KEY
from main.celery_util import get_celery_worker_status
from main.models import AlignmentGroup
from read_alignment import align_with_bwa_mem
from snv_calling import get_variant_tool_params
from snv_calling import find_variants_with_tool


# The default alignment function.
ALIGNMENT_FN = align_with_bwa_mem


def _assert_celery_running():
    """Make sure celery is running, unless settings.CELERY_ALWAYS_EAGER = True.
    """
    if settings.CELERY_ALWAYS_EAGER:
        return
    celery_status = get_celery_worker_status()
    assert not CELERY_ERROR_KEY in celery_status, celery_status[CELERY_ERROR_KEY]


def run_pipeline_multiple_ref_genomes(alignment_group_label, ref_genome_list,
        sample_list, test_models_only=False):
    """
    Runs the pipeline for each ref genome in the ref_genome_list.

    We create an AlignmentGroup for each ReferenceGenome and align all
    ExpermentSamples to each ReferenceGenome separately.

    Args:
        ref_genome_list: List of ReferenceGenome instances.
        sample_list: List of sample instances. Must belong to same project as
            ReferenceGenomes.
        test_models_only: If True, don't actually run alignments. Just create
            models.
    """

    assert len(alignment_group_label) > 0, "Name must be non-trivial string."
    assert len(ref_genome_list) > 0, (
            "Must provide at least one ReferenceGenome.")
    assert len(sample_list) > 0, (
            "Must provide at least one ExperimentSample.")
    _assert_celery_running()

    # Save the alignment group objects for returning if required.
    alignment_groups = {}

    for ref_genome in ref_genome_list:
        alignment_groups[ref_genome.uid] = run_pipeline(
                alignment_group_label + '_' + ref_genome.label,
                ref_genome, sample_list)

    # Return a dictionary of all alignment groups indexed by ref_genome uid.
    return alignment_groups


def run_pipeline(alignment_group_label, ref_genome, sample_list,
        test_models_only=False):
    """
    Creates an AlignmentGroup if not created and kicks off alignment for each one.

    Args:
        ref_genome: ReferenceGenome instance
        sample_list: List of sample instances. Must belong to same project as
            ReferenceGenomes.
        test_models_only: If True, don't actually run alignments. Just create
            models.
    """

    assert len(alignment_group_label) > 0, "Name must be non-trivial string."
    assert len(sample_list) > 0, (
            "Must provide at least one ExperimentSample.")
    _assert_celery_running()

    (alignment_group, created) = AlignmentGroup.objects.get_or_create(
            label=alignment_group_label,
            reference_genome=ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

    # Kick of the alignments concurrently.
    # TODO(gleb): Use this list to block on when integrating with
    # variant calling.

    # Since we don't want results to be passed as arguments in the
    # chain, use .si(...) and not .s(...)
    # http://stackoverflow.com/
    #       questions/15224234/celery-chaining-tasks-sequentially

    alignment_tasks = []
    for sample in sample_list:
        # create a task signature for this subtask
        align_task_signature = ALIGNMENT_FN.si(
                alignment_group, sample, None, test_models_only,
                project=ref_genome.project)

        alignment_tasks.append(align_task_signature)

    align_task_group = group(alignment_tasks)

    # create signatures for all variant tools
    variant_callers = []
    for variant_params in get_variant_tool_params():
        variant_caller_signature = find_variants_with_tool.si(
                alignment_group, variant_params, project=ref_genome.project)

        variant_callers.append(variant_caller_signature)

    variant_caller_group = group(variant_callers)

    whole_pipeline = chain(align_task_group, variant_caller_group)

    # now, run the whole pipeline
    whole_pipeline.apply_async()

    return alignment_group
