
"""
The alignment pipeline.

We start with the .fastq files and the reference for a particular genome
and carry out the gauntlet of steps to perform alignments, alignment 
cleaning, snv calling, and effect prediction.
"""

from celery import group
from celery import chain
from celery import task
import copy
import os
import subprocess
from subprocess import CalledProcessError
import sys

from main.celery_util import CELERY_ERROR_KEY
from main.celery_util import get_celery_worker_status
from main.models import clean_filesystem_location
from main.models import get_dataset_with_type
from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from scripts.import_util import add_dataset_to_entity
from scripts.jbrowse_util import add_bam_file_track
from scripts.jbrowse_util import prepare_jbrowse_ref_sequence
from snv_calling import get_variant_tool_params
from snv_calling import find_variants_with_tool
from read_alignment import align_with_bwa_mem
from settings import DEBUG_CONCURRENT
from settings import PWD
from settings import TOOLS_DIR
from settings import BASH_PATH

ALIGNMENT_FN = align_with_bwa_mem

def run_pipeline(
        alignment_group_label, ref_genome_list, sample_list,
        test_models_only=False):

    """Creates an AlignmentGroup and kicks off alignment for each one.

    We create an AlignmentGroup for each ReferenceGenome and align all
    ExpermentSamples to each ReferenceGenome separately.

    Args:
        ref_genome_list: List of ReferenceGenome instances.
        sample_list: List of sample instances. Must belong to same project as
            ReferenceGenomes.
        test_models_only: If True, don't actually run alignments. Just create
            models.
        concurrent: If True, run alignments in parallel.
    """

    assert len(alignment_group_label) > 0, "Name must be non-trivial string."
    assert len(ref_genome_list) > 0, (
            "Must provide at least one ReferenceGenome.")
    assert len(sample_list) > 0, (
            "Must provide at least one ExperimentSample.")

    # If running concurrently, check whether Celery is running.
    celery_status = get_celery_worker_status()
    assert not CELERY_ERROR_KEY in celery_status, celery_status[CELERY_ERROR_KEY]

    # Save the alignment group objects for returning if required.
    alignment_groups = {}

    for ref_genome in ref_genome_list:
        alignment_group = AlignmentGroup.objects.create(
                label=alignment_group_label + '_' + ref_genome.label,
                reference_genome=ref_genome,
                aligner=AlignmentGroup.ALIGNER.BWA)

        alignment_groups[ref_genome.uid] = alignment_group

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

    # Return a dictionary of all alignment groups indexed by ref_genome uid.
    return alignment_groups

def resume_alignment_group_alignment(alignment_group, sample_list,
        concurrent=DEBUG_CONCURRENT):
    """Resumes alignments that are not complete for an alignment group.

    Any samples that are not in the ready state are deleted and re-started.
    """

    # This should be a part of run_alignment_pipeline, which will 
    # intelligently decide to resume or restart all alignments and
    # re-call SNPs as necessary. For now, let's just not do this. -dbg
    raise NotImplementedError

    # If running concurrently, check whether Celery is running.
    if concurrent:
        celery_status = get_celery_worker_status()
        assert not CELERY_ERROR_KEY in celery_status, (
                celery_status[CELERY_ERROR_KEY])

    for sample in sample_list:
        # Skip samples that have already been aligned.
        maybe_sample_alignment_list = (
                alignment_group.experimentsampletoalignment_set.filter(
                    experiment_sample=sample))
        if len(maybe_sample_alignment_list):
            sample_alignment = maybe_sample_alignment_list[0]
            maybe_dataset_set = sample_alignment.dataset_set.filter(
                type=Dataset.TYPE.BWA_ALIGN)
            if len(maybe_dataset_set):
                bwa_dataset = maybe_dataset_set[0]
                if bwa_dataset.status == Dataset.STATUS.READY:
                    # Skip repeating alignment.
                    continue
                else:
                    # Delete the sample_alignment and its data and just start
                    # over.
                    sample_alignment.delete()

        # Carry out alignment.
        args = [alignment_group, sample, None, False]
        # fn_runner(ALIGNMENT_FN, alignment_group.reference_genome.project,
        #         args, concurrent=concurrent)


def wait_and_check_pipe(pipe, popenargs=None):
    """Similar to subprocess.check_call() except takes a running pipe
    as input.

    Args:
        pipe: The running pipe.
        popenargs: Optional list of args passed to Popen for error reporting.
    """
    retcode = pipe.wait()
    if retcode:
        if popenargs:
            cmd = popenargs[0]
        else:
            cmd = ''
        raise CalledProcessError(retcode, cmd)
    return 0
