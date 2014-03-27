"""
The alignment pipeline.

We start with the .fastq files and the reference for a particular genome
and carry out the gauntlet of steps to perform alignments, alignment
cleaning, snv calling, and effect prediction.
"""

from datetime import datetime

from celery import chain
from celery import group
from celery import task

from main.celery_util import assert_celery_running
from main.models import AlignmentGroup
from main.models import Dataset
from read_alignment import align_with_bwa_mem
from snv_calling import get_variant_tool_params
from snv_calling import find_variants_with_tool


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

    # Make sure all samples are ready.
    relevant_datasets = Dataset.objects.filter(
            experimentsample__in=sample_list)
    for d in relevant_datasets:
        assert d.status == Dataset.STATUS.READY, (
                "Dataset %s for sample %s has status %s. Expected %s." % (
                        d.label, d.experimentsample_set.all()[0].label,
                        d.status, Dataset.STATUS.READY))

    assert_celery_running()

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
    assert_celery_running()

    (alignment_group, created) = AlignmentGroup.objects.get_or_create(
            label=alignment_group_label,
            reference_genome=ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

    # Kick of the alignments concurrently.

    # Since we don't want results to be passed as arguments in the
    # chain, use .si(...) and not .s(...)
    # http://stackoverflow.com/
    #       questions/15224234/celery-chaining-tasks-sequentially

    alignment_tasks = []
    for sample in sample_list:
        # create a task signature for this subtask
        align_task_signature = align_with_bwa_mem.si(
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

    pipeline_completion = pipeline_completion_tasks.s(
            alignment_group=alignment_group)

    whole_pipeline = chain(
            align_task_group, 
            variant_caller_group, 
            pipeline_completion)

    # now, run the whole pipeline
    whole_pipeline.apply_async()

    return alignment_group

@task
def pipeline_completion_tasks(variant_caller_group_result, alignment_group):
    """
    Things that happen when the pipeline completes that don't need
    to be run in parallel.
    """
    
    if hasattr(variant_caller_group_result, 'join'):
        variant_caller_group_result.join()

    print 'FINISHING PIPELINE.'

    # The alignment group is now officially complete.
    alignment_group.status = AlignmentGroup.STATUS.COMPLETED
    alignment_group.end_time = datetime.now()
    alignment_group.save()

    print alignment_group.status
