"""
Functions for running custom SV methods based on assembly and detecting
low coverage.

Celery tasks are here. Implementations are in assembly.py.
"""

from celery import group
from celery import task

from genome_finish.assembly import clean_up_previous_runs_of_sv_calling_pipeline
from genome_finish.assembly import FAILURE_REPORT__CONTIG
from genome_finish.assembly import FAILURE_REPORT__DETECT_DELETION
from genome_finish.assembly import FAILURE_REPORT__PARSE_VARIANTS
from genome_finish.assembly import generate_contigs
from genome_finish.assembly import parse_variants_from_vcf
from genome_finish.celery_task_decorator import report_failure_stats
from genome_finish.celery_task_decorator import set_assembly_status
from genome_finish.detect_deletion import cov_detect_deletion_make_vcf
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from pipeline.read_alignment_util import ensure_bwa_index
from utils.jbrowse_util import compile_tracklist_json
from utils.jbrowse_util import prepare_jbrowse_ref_sequence


def run_de_novo_assembly_pipeline(sample_alignment_list,
        sv_read_classes={}, input_velvet_opts={},
        overwrite=True):
    """Kicks off Millstone's custom SV calling pipeline.

    NOTE: Despite the name, in addition to de novo assembly, our custom
    SV-calling pipeline also uses non-assembly based methods like low-coverage
    detection to call deletions.
    """
    # First, we delete any data from previous runs of this custom SV-calling
    # pipeline, and update the status of the sample alignments to indicate
    # that custom SV-calling is taking place.

    for sample_alignment in sample_alignment_list:
        set_assembly_status(
                sample_alignment,
                ExperimentSampleToAlignment.ASSEMBLY_STATUS.CLEARING,
                force=True)
        clean_up_previous_runs_of_sv_calling_pipeline(sample_alignment)
        set_assembly_status(
                sample_alignment,
                ExperimentSampleToAlignment.ASSEMBLY_STATUS.QUEUED, force=True)

    # Recompile the tracklist after deleting the indiv_tracks dirs for these
    # deleted contigs.
    ref_genome = sample_alignment_list[0].alignment_group.reference_genome
    compile_tracklist_json(ref_genome)

    # Next, we ensure reference genome fasta is indexed. We do this before
    # the async tasks.
    ref_genome_fasta = ref_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    ensure_bwa_index(ref_genome_fasta)
    prepare_jbrowse_ref_sequence(ref_genome)

    # Finally we assemble the async tasks that be parallelized.
    async_result = get_sv_caller_async_result(
            sample_alignment_list)

    return async_result


def single_sample_alignment_assembly(sample_alignment):
    """
    Run the assembly pipeline without celery on a single sample_alignment.
    Used for testing & debug.
    """

    set_assembly_status(
            sample_alignment,
            ExperimentSampleToAlignment.ASSEMBLY_STATUS.CLEARING,
            force=True)

    clean_up_previous_runs_of_sv_calling_pipeline(sample_alignment)

    set_assembly_status(
            sample_alignment,
            ExperimentSampleToAlignment.ASSEMBLY_STATUS.QUEUED, force=True)

    ref_genome = sample_alignment.alignment_group.reference_genome

    compile_tracklist_json(ref_genome)

    ref_genome_fasta = ref_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    ensure_bwa_index(ref_genome_fasta)
    prepare_jbrowse_ref_sequence(ref_genome)

    generate_contigs(sample_alignment)
    cov_detect_deletion_make_vcf(sample_alignment)
    parse_variants_from_vcf(sample_alignment)


@task
def _chordfinisher(*args, **kwargs):
    """
    Needs to run at the end of a chord to delay the variant parsing step.

    http://stackoverflow.com/questions/
    15123772/celery-chaining-groups-and-subtasks-out-of-order-execution
    """
    return "FINISHED VARIANT FINDING."


def get_sv_caller_async_result(sample_alignment_list):
    """Builds a celery chord that contains tasks for calling SVs for each
    ExperimentSampleToAlignment in sample_alignment_list in parallel. Each task
    generates vcfs, named according to the method use to call the contained
    variants. The callback to the chord is a chain of tasks (applied
    synchronously) that parse variants from vcfs.

    Returns an AsyncResult object.
    """
    generate_contigs_tasks = []
    cov_detect_deletion_tasks = []
    for sample_alignment in sorted(sample_alignment_list,
            key=lambda x: x.experiment_sample.label):

        # These tasks are based on de novo assembly.
        generate_contigs_tasks.append(
                generate_contigs_async.si(sample_alignment))

        # These tasks use coverage to call large deletions.
        cov_detect_deletion_tasks.append(
                cov_detect_deletion_make_vcf_async.si(sample_alignment))

    variant_finding = group(generate_contigs_tasks + cov_detect_deletion_tasks)

    sv_task_chain = (variant_finding |
            _chordfinisher.si() |
            parse_variants_for_sa_list_async.si(sample_alignment_list))

    return sv_task_chain()


@task(ignore_result=False)
@report_failure_stats(FAILURE_REPORT__CONTIG)
def generate_contigs_async(sample_alignment,
        sv_read_classes={}, input_velvet_opts={},
        overwrite=True):
    """
    Async wrapper for contigs function.
    """
    generate_contigs(sample_alignment,
            sv_read_classes={}, input_velvet_opts={},
            overwrite=True)


@task(ignore_result=False)
@report_failure_stats(FAILURE_REPORT__DETECT_DELETION)
def cov_detect_deletion_make_vcf_async(sample_alignment):
    """
    Async wrapper for deletion coverage function.
    """
    cov_detect_deletion_make_vcf(sample_alignment)


@task(ignore_result=False)
@report_failure_stats(FAILURE_REPORT__PARSE_VARIANTS)
def parse_variants_for_sa_list_async(sample_alignment_list):
    """
    Async wrapper for generation of vcf variants from SV calls.
    """
    for sample_alignment in sorted(sample_alignment_list,
            key=lambda x: x.experiment_sample.label):
        parse_variants_from_vcf(sample_alignment)
