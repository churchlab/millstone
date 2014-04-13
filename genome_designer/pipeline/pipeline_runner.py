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
from main.models import ReferenceGenome
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from pipeline.read_alignment import align_with_bwa_mem
from pipeline.snv_calling import get_variant_tool_params
from pipeline.snv_calling import find_variants_with_tool
import settings


def run_pipeline_multiple_ref_genomes(alignment_group_label, ref_genome_list,
        sample_list):
    """
    Runs the pipeline for each ref genome in the ref_genome_list.

    We create an AlignmentGroup for each ReferenceGenome and align all
    ExpermentSamples to each ReferenceGenome separately.

    Args:
        ref_genome_list: List of ReferenceGenome instances.
        sample_list: List of sample instances. Must belong to same project as
            ReferenceGenomes.
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


def run_pipeline(alignment_group_label, ref_genome, sample_list):
    """
    Creates an AlignmentGroup if not created and kicks off alignment for each one.

    Args:
        ref_genome: ReferenceGenome instance
        sample_list: List of sample instances. Must belong to same project as
            ReferenceGenomes.
    """

    assert len(alignment_group_label) > 0, "Name must be non-trivial string."
    assert len(sample_list) > 0, (
            "Must provide at least one ExperimentSample.")
    assert_celery_running()

    alignment_group, _ = AlignmentGroup.objects.get_or_create(
            label=alignment_group_label,
            reference_genome=ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

    # The pipeline has two phases:
    # 1) Alignments - run in parallel.
    # 2) Variant calling - each variant caller runs in parallel, but waits
    #    for all alignments to be complete before starting.

    # NOTE: Since we don't want results to be passed as arguments in the
    # chain, use .si(...) and not .s(...)
    # See: http://stackoverflow.com/questions/15224234/celery-chaining-tasks-sequentially

    # First we create Models so that we can track status from the ui
    # immediately. The ui determines the status of each
    # ExperimentSampleToAlignment by looking at the status of its BWA dataset.
    sample_alignments_to_run = []
    for sample in sample_list:
        sample_alignment, _ = ExperimentSampleToAlignment.objects.get_or_create(
                alignment_group=alignment_group, experiment_sample=sample)

        # Get or create a Dataset to store the alignment result.
        sample_alignment_bwa_datasets = sample_alignment.dataset_set.filter(
                type=Dataset.TYPE.BWA_ALIGN)
        assert len(sample_alignment_bwa_datasets) <= 1
        if len(sample_alignment_bwa_datasets) == 1:
            bwa_dataset = sample_alignment_bwa_datasets[0]
        else:
            bwa_dataset = Dataset.objects.create(
                    label=Dataset.TYPE.BWA_ALIGN,
                    type=Dataset.TYPE.BWA_ALIGN,
                    status=Dataset.STATUS.NOT_STARTED)
            sample_alignment.dataset_set.add(bwa_dataset)
            sample_alignment.save()

        # Add it to the list of alignments to run, unless already done.
        if not bwa_dataset.status == Dataset.STATUS.READY:
            sample_alignments_to_run.append(sample_alignment)

    # Now we aggregate the alignments that need to be run, collecting their
    # signatures in a Celery group so that these alignments can be run in
    # parallel.
    # Before we do so, let's update the ref genome object.
    ref_genome = ReferenceGenome.objects.get(uid=ref_genome.uid)

    alignment_task_signatures = [align_with_bwa_mem.si(
                    alignment_group, sample_alignment,
                    project=ref_genome.project)
            for sample_alignment in sample_alignments_to_run]
    align_task_group = group(alignment_task_signatures)

    # Aggregate variant callers, which run in parallel once all alignments
    # are done.
    variant_caller_group = group([find_variants_with_tool.si(
                    alignment_group, variant_params,
                    project=ref_genome.project)
            for variant_params in get_variant_tool_params() if
                    variant_params[0] in settings.ENABLED_VARIANT_CALLERS])

    pipeline_completion = pipeline_completion_tasks.s(
            alignment_group=alignment_group)

    # Put together the whole pipeline.
    # We check for there being more than 0 alignments since celery
    # expects groups to have at least 1 task. Is there a better way to do this?
    if len(alignment_task_signatures) > 0:
        whole_pipeline = chain(
                align_task_group,
                variant_caller_group,
                pipeline_completion)
    else:
        whole_pipeline = chain(
                variant_caller_group,
                pipeline_completion)

    # Run the pipeline.
    ref_genome.save()
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
