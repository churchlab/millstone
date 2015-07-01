"""
The alignment pipeline.

We start with the .fastq files and the reference for a particular genome
and carry out the gauntlet of steps to perform alignments, alignment
cleaning, snv calling, and effect prediction.
"""

from datetime import datetime
import time

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
from pipeline.variant_calling import TOOL_FREEBAYES
from pipeline.variant_calling import TOOL_LUMPY
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
    #

    # NOTE: Nested chords in celery don't work so we need to break up the
    # pipeline into # two separate pipelines: 1) alignment and 2) variant
    # calling. This task is the first task in the variant calling pipeline which
    # polls the database until all alignments are complete before kicking off
    # parallel variant calling tasks.

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

    if len(alignment_task_signatures) > 0:
        alignment_task_group = group(alignment_task_signatures)
        alignment_task_group_async_result = alignment_task_group.apply_async()
    else:
        alignment_task_group_async_result = None

    # HACK(gleb): Force ALIGNING so that UI starts refreshing. This should be
    # right, but I'm open to removing if it's not right for some case I
    # didn't think of.
    alignment_group.status = AlignmentGroup.STATUS.ALIGNING
    alignment_group.start_time = datetime.now()
    alignment_group.end_time = None
    alignment_group.save(update_fields=['status', 'start_time', 'end_time'])

    # Aggregate variant callers, which run in parallel once all alignments
    # are done.
    if perform_variant_calling:
        variant_caller_group = _construct_variant_caller_group(
                alignment_group, variant_calling_options)
    else:
        variant_caller_group = None

    # Put together the whole pipeline.
    variant_calling_pipeline = start_variant_calling_pipeline_task.si(
            alignment_group)
    if variant_caller_group is not None:
        variant_calling_pipeline = (variant_calling_pipeline |
                variant_caller_group)

    # Add a final task which runs only after all previous tasks are complete.
    pipeline_completion = pipeline_completion_tasks.si(alignment_group)
    variant_calling_pipeline = variant_calling_pipeline | pipeline_completion

    # TODO(gleb): We had this to deal with race conditions. Do we still need it?
    ref_genome.save()

    # Run the pipeline. This is a non-blocking call when celery is running so
    # the rest of code proceeds immediately.
    variant_calling_async_result = variant_calling_pipeline.apply_async()

    return (
            alignment_group,
            alignment_task_group_async_result,
            variant_calling_async_result)


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

    # Determine which variant callers to use.
    if VARIANT_CALLING_OPTION__CALLER_OVERRIDE in variant_calling_options:
        effective_variant_callers = variant_calling_options[
                VARIANT_CALLING_OPTION__CALLER_OVERRIDE]
    else:
        effective_variant_callers = settings.ENABLED_VARIANT_CALLERS

    # List of tasks that can be run in parallel. These will be combined into a
    # single celery.group.
    parallel_tasks = []

    # Iterate through tools and kick off tasks.
    for tool in effective_variant_callers:
        # Common params for this tool.
        tool_params = VARIANT_TOOL_PARAMS_MAP[tool]

        parallel_tasks = []

        if settings.FREEBAYES_PARALLEL and tool == TOOL_FREEBAYES:
            # Special handling for freebayes if running parallel. Break up
            # ReferenceGenome into regions and create separate job for each.
            fb_regions = freebayes_regions(ref_genome)
            assert len(fb_regions) >= 0

            for region_num, fb_region in enumerate(fb_regions):
                region_params = dict(tool_params)
                region_params['tool_kwargs'] = {
                    'region': fb_region,
                    'region_num': region_num
                }
                parallel_tasks.append(find_variants_with_tool.si(
                        alignment_group, region_params,
                        project=ref_genome.project))
        elif tool == TOOL_LUMPY:
            sample_alignment_list = (
                    alignment_group.experimentsampletoalignment_set.all())
            # TODO: What if some alignments failed?
            for sa in sample_alignment_list:
                # Create separate lumpy task for each sample.
                per_sample_params = dict(tool_params)
                per_sample_params['tool_kwargs'] = {
                    'region_num': sa.id,
                    'sample_alignments': [sa]
                }
                parallel_tasks.append(find_variants_with_tool.si(
                        alignment_group, per_sample_params,
                        project=ref_genome.project))
        else:
            parallel_tasks.append(find_variants_with_tool.si(
                    alignment_group, tool_params, project=ref_genome.project))

    variant_calling_pipeline = (group(parallel_tasks) |
            merge_variant_data.si(alignment_group))
    return variant_calling_pipeline


@task
def start_variant_calling_pipeline_task(alignment_group):
    """First task in variant calling pipeline which waits for all alignments
    to be complete.

    Nested chords in celery don't work so we need to break up the pipeline into
    two separate pipelines: 1) alignment and 2) variant calling. This task is
    the first task in the variant calling pipeline which polls the database
    until all alignments are complete.
    """
    print 'START VARIANT CALLING PIPELINE. WAITING FOR ALIGNMENTS TO COMPLETE.'

    POLL_INTERVAL_SEC = 5

    sample_alignment_list = ExperimentSampleToAlignment.objects.filter(
            alignment_group=alignment_group)

    all_samples_ready = False
    failed = False
    while not all_samples_ready:
        all_samples_ready = True
        for sa in sample_alignment_list:
            sa_fresh = ExperimentSampleToAlignment.objects.get(id=sa.id)
            bwa_dataset = sa_fresh.dataset_set.get(label=Dataset.TYPE.BWA_ALIGN)
            if bwa_dataset.status == Dataset.STATUS.FAILED:
                failed = True
                # # DEBUG: Uncomment to see alignment error output
                # print '---------------ALIGNMENT ERROR OUTPUT---------------\n'
                # error_dataset = sa_fresh.dataset_set.get(
                #         label=Dataset.TYPE.BWA_ALIGN_ERROR)
                # with open(error_dataset.get_absolute_location(),'r') as fh:
                #     for line in fh:
                #         print line
                break

            if not bwa_dataset.status == Dataset.STATUS.READY:
                all_samples_ready = False
                break

        if failed:
            alignment_group = AlignmentGroup.objects.get(id=alignment_group.id)
            alignment_group.status = AlignmentGroup.STATUS.FAILED
            alignment_group.save(update_fields=['status'])
            raise Exception("Alignment failed.")

        if not all_samples_ready:
            time.sleep(POLL_INTERVAL_SEC)

    # All ready. Set VARIANT_CALLING.
    alignment_group.status = AlignmentGroup.STATUS.VARIANT_CALLING
    alignment_group.save(update_fields=['status'])


@task
def merge_variant_data(alignment_group):
    """Merges results of variant caller data after pipeline is complete.
    """
    merge_freebayes_parallel(alignment_group)


@task
def pipeline_completion_tasks(alignment_group):
    """Final set of synchronous steps after all alignments and variant callers
    are finished.

    Sets end_time and status on alignment_group.
    """
    print 'START PIPELINE COMPLETION'
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
