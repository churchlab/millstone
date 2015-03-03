"""
The alignment pipeline.

We start with the .fastq files and the reference for a particular genome
and carry out the gauntlet of steps to perform alignments, alignment
cleaning, snv calling, and effect prediction.
"""

from datetime import datetime

from celery import chord
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


# List of variant callers to use. At time of writing, this was not hooked
# up to the ui and only used internally.
VARIANT_CALLING_OPTION__CALLER_OVERRIDE = 'enabled_variant_callers_override'


def run_pipeline(alignment_group_label, ref_genome, sample_list,
        skip_alignment=False, perform_variant_calling=True, alignment_options={},
        variant_calling_options={}):
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
        skip_alignment: If True, skip alignment.
        perform_variant_calling: Whether to run variant calling.
        alignment_options: Control aspects of alignment.
        variant_calling_options: Control aspects of calling variants.

    Returns:
        Tuple pair (alignment_group, async_result).
    """
    if not skip_alignment:
        _assert_pipeline_is_safe_to_run(alignment_group_label, sample_list)

    # Create AlignmentGroup, the entity which groups together the alignments
    # of individual samples, and results of variant calling which happens
    # for all samples together.
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
    sample_alignments_to_run = _get_or_create_sample_alignment_datasets(
            alignment_group, sample_list)

    # Before we continue, let's update the ref genome object. This is code
    # left over from when we were fighting a concurrency bug.
    # TODO: Revisit such calls and see if we can clean them up.
    ref_genome = ReferenceGenome.objects.get(uid=ref_genome.uid)

    # Now we aggregate the alignments that need to be run, collecting their
    # signatures in a Celery group so that these alignments can be run in
    # parallel.
    alignment_task_signatures = []
    for sample_alignment in sample_alignments_to_run:
        alignment_task_signatures.append(
                align_with_bwa_mem.si(
                        alignment_group, sample_alignment,
                        project=ref_genome.project))
    # Aggregate variant callers, which run in parallel once all alignments
    # are done.
    if perform_variant_calling:
        variant_caller_group = _construct_variant_caller_group(
                alignment_group, variant_calling_options)
    else:
        variant_caller_group = None

    # We add a final task which runs only after all previous tasks are complete.
    pipeline_completion = pipeline_completion_tasks.si(alignment_group)

    # Put together the whole pipeline.
    whole_pipeline = None
    if len(alignment_task_signatures) > 0:
        whole_pipeline = group(alignment_task_signatures)
    if variant_caller_group is not None:
        if whole_pipeline is None:
            whole_pipeline = variant_caller_group
        else:
            whole_pipeline = whole_pipeline | variant_caller_group
    whole_pipeline = whole_pipeline | pipeline_completion

    # TODO(gleb): We had this to deal with race conditions. Do we still need it?
    ref_genome.save()

    # HACK(gleb): Force ALIGNING so that UI starts refreshing. This should be
    # right, but I'm open to removing if it's not right for some case I
    # didn't think of.
    alignment_group.status = AlignmentGroup.STATUS.ALIGNING
    alignment_group.start_time = datetime.now()
    alignment_group.end_time = None
    alignment_group.save(update_fields=['status', 'start_time', 'end_time'])

    # Run the pipeline. This is a non-blocking call when celery is running so
    # the rest of code proceeds immediately.
    async_result = whole_pipeline.apply_async()

    return (alignment_group, async_result)


def _assert_pipeline_is_safe_to_run(alignment_group_label, sample_list):
    """Helper that checks that pipeline is ready to run.

    Raises:
        AssertionError if any problems.
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


def _get_or_create_sample_alignment_datasets(alignment_group, sample_list):
    """Creates Dataset models that allow tracking status of alignment from ui.

    Does not start alignments.

    Returns list of ExperimentSampleToAlignments.
    """
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
    return sample_alignments_to_run


def _construct_variant_caller_group(alignment_group, variant_calling_options):
    """Returns celery Group of variant calling tasks that can be run
    in parallel.
    """
    # Get fresh copy of ReferenceGenome to avoid potential issues with
    # race conditions.
    ref_genome = ReferenceGenome.objects.get(
            uid=alignment_group.reference_genome.uid)

    # List of tasks that can be run in parallel. These will be combined into a
    # single celery.group.
    parallel_tasks = []

    # Determine which variant callers to use.
    if VARIANT_CALLING_OPTION__CALLER_OVERRIDE in variant_calling_options:
        effective_variant_callers = variant_calling_options[
                VARIANT_CALLING_OPTION__CALLER_OVERRIDE]
    else:
        effective_variant_callers = settings.ENABLED_VARIANT_CALLERS

    # Iterate through tools and kick off tasks.
    for tool in effective_variant_callers:
        # Common params for this tool.
        tool_params = VARIANT_TOOL_PARAMS_MAP[tool]

        if settings.FREEBAYES_PARALLEL and tool == 'freebayes':
            # Special handling for freebayes if running parallel.
            # Create a chord consisting of the group of freebayes tasks
            # distributed across the regions of the genome and finally a call
            # to merge results into single vcf add corresponding Variant
            # models to db.
            fb_parallel_tasks = []

            # Break up ReferenceGenome into regions and create separate
            # job for each.
            fb_regions = freebayes_regions(ref_genome)

            for region_num, fb_region in enumerate(fb_regions):
                region_params = dict(tool_params)
                region_params['tool_kwargs'] = {
                    'region': fb_region,
                    'region_num': region_num
                }
                fb_parallel_tasks.append(find_variants_with_tool.si(
                        alignment_group, region_params,
                        project=ref_genome.project))

            # Create chord so that merge task occurs after all freebayes
            # tasks are complete.
            freebayes_chord = (group(fb_parallel_tasks) |
                    merge_freebayes_parallel.si(alignment_group))
            parallel_tasks.append(freebayes_chord)
        else:
            parallel_tasks.append(find_variants_with_tool.si(
                    alignment_group, tool_params,
                    project=ref_genome.project))

    return group(parallel_tasks)


@task
def pipeline_completion_tasks(alignment_group):
    """Final set of synchronous steps after all alignments and variant callers
    are finished.

    Sets end_time and status on alignment_group.
    """
    try:
        # Get fresh copy of alignment_group.
        alignment_group = AlignmentGroup.objects.get(id=alignment_group.id)
        assert AlignmentGroup.STATUS.VARIANT_CALLING == alignment_group.status

        # The alignment pipeline is now officially complete.
        alignment_group.status = AlignmentGroup.STATUS.COMPLETED
    except:
        # TODO(gleb): Failure logging.
        alignment_group.status = AlignmentGroup.STATUS.FAILED

    alignment_group.end_time = datetime.now()
    alignment_group.save(update_fields=['end_time', 'status'])
