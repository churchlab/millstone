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
from django.conf import settings

from main.celery_util import assert_celery_running
from main.models import AlignmentGroup
from main.models import ReferenceGenome
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from pipeline.read_alignment import align_with_bwa_mem
from pipeline.variant_calling import find_variants_with_tool
from pipeline.variant_calling import VARIANT_TOOL_PARAMS_MAP
from pipeline.variant_calling.freebayes import merge_freebayes_parallel
from pipeline.variant_calling.freebayes import freebayes_regions


def run_pipeline(alignment_group_label, ref_genome, sample_list,
        perform_variant_calling=True, alignment_options={}):
    """Runs the entire bioinformatics pipeline, including alignment and
    variant calling.

    Steps:
        * Create AlignmentGroup if not created
        * get_or_create ExperimentSampleToAlignments and respective Datasets
        * Kick off alignments
        * When all alignments are done, kick off variant calling

    Args:
        alignment_group_label: Name for this alignment.
        ref_genome: ReferenceGenome instance
        sample_list: List of sample instances. Must belong to same project as
            ReferenceGenomes.

    Returns:
        Tuple pair (alignment_group, async_result).
    """
    assert len(alignment_group_label) > 0, "Name must be non-trivial string."
    assert len(sample_list) > 0, (
            "Must provide at least one ExperimentSample.")
    assert_celery_running()

    # Make sure all samples are ready.
    relevant_datasets = Dataset.objects.filter(
            experimentsample__in=sample_list)
    for d in relevant_datasets:
        assert d.status == Dataset.STATUS.READY, (
                "Dataset %s for sample %s has status %s. Expected %s." % (
                        d.label, d.experimentsample_set.all()[0].label,
                        d.status, Dataset.STATUS.READY))

    alignment_group, _ = AlignmentGroup.objects.get_or_create(
            label=alignment_group_label,
            reference_genome=ref_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

    alignment_group.alignment_options.update(alignment_options)

    # The pipeline has two synchronous phases, each of whose components
    # maybe run in parallel:
    #     1) Alignments - run in parallel.
    #     2) Variant calling - each variant caller runs in parallel, but waits
    #         for all alignments to be complete before starting.

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

        # Add it to the list of alignments to run, unless already done.
        if not bwa_dataset.status == Dataset.STATUS.READY:
            sample_alignments_to_run.append(sample_alignment)

    # Before we continue, let's update the ref genome object. This is code
    # left over from when we were fighting a concurrency bug.
    # TODO: Revisit such calls and see if we can clean them up.
    ref_genome = ReferenceGenome.objects.get(uid=ref_genome.uid)

    # Now we aggregate the alignments that need to be run, collecting their
    # signatures in a Celery group so that these alignments can be run in
    # parallel.
    alignment_task_signatures = [align_with_bwa_mem.si(
                    alignment_group, sample_alignment,
                    project=ref_genome.project)
            for sample_alignment in sample_alignments_to_run]
    align_task_group = group(alignment_task_signatures)

    # Aggregate variant callers, which run in parallel once all alignments
    # are done.
    if perform_variant_calling:

        # If we are parallelizing freebayes, then we need to kick off
        # several tasks for each parallelized region:
        if settings.FREEBAYES_PARALLEL:

            # a list of variant callers that we will run in parallel, with
            # parameters.
            variant_param_list = []

            # determine the regions to send to each freebayes worker
            fb_regions = freebayes_regions(ref_genome)

            for tool in settings.ENABLED_VARIANT_CALLERS:

                params = VARIANT_TOOL_PARAMS_MAP[tool]

                #create a variant_param dictionary for each region
                if tool == 'freebayes':
                    for region_num, fb_region in enumerate(fb_regions):
                        this_region_num = region_num
                        region_params = dict(params)
                        region_params['tool_kwargs'] = {
                                    'region':fb_region,
                                    'region_num':this_region_num}
                        variant_param_list.append(region_params)
                else:
                    variant_param_list.append(params)

        # no freebayes parallel:
        else:
            variant_param_list = [VARIANT_TOOL_PARAMS_MAP[tool]
                    for tool in settings.ENABLED_VARIANT_CALLERS]

        variant_caller_group = group([
                find_variants_with_tool.si(
                        alignment_group, variant_params,
                        project=ref_genome.project)
                for variant_params in variant_param_list])
    else:
        variant_caller_group = None

    # We add a final task which runs only after all previous tasks are complete.
    merge_fb_parallel = perform_variant_calling and settings.FREEBAYES_PARALLEL
    pipeline_completion = pipeline_completion_tasks.s(
            alignment_group=alignment_group,
            merge_fb_parallel=merge_fb_parallel)

    # Put together the whole pipeline.
    # Since this method might be called after alignments have already been
    # complete and only variant calling needs to be done, it's possible
    # for there to be no alignment tasks to do, in which case we skip adding
    # an empty group to the chain else celery will throw an error.
    groups = []
    if len(alignment_task_signatures) > 0:
        groups.append(align_task_group)
    if variant_caller_group is not None:
        groups.append(variant_caller_group)
    groups.append(pipeline_completion)
    whole_pipeline = chain(*groups)

    # Run the pipeline.
    ref_genome.save()
    async_result = whole_pipeline.apply_async()

    return (alignment_group, async_result)

@task
def pipeline_completion_tasks(variant_caller_group_result, alignment_group,
        merge_fb_parallel):
    """Final set of synchronous steps after all alignments and variant callers
    are finished.

    If we ran freebayes, see if it was parallelized, and merge and add
    the resulting vcf.

    Otherwise, just check for completion, timestamp, and save
    the alignment group.
    """
    if hasattr(variant_caller_group_result, 'join'):
        variant_caller_group_result.join()

    # Get fresh copy of alignment_group.
    alignment_group = AlignmentGroup.objects.get(id=alignment_group.id)

    if alignment_group.status != AlignmentGroup.STATUS.VARIANT_CALLING:
        alignment_group.status = AlignmentGroup.STATUS.FAILED
    else:
        # if ran freebayes in parallel, merge the partial vcfs and process the
        # vcf dataset.
        if merge_fb_parallel:
            merge_freebayes_parallel(alignment_group)

        # The alignment group is now officially complete.
        alignment_group.status = AlignmentGroup.STATUS.COMPLETED

    alignment_group.end_time = datetime.now()
    alignment_group.save()
