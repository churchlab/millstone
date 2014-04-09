"""
Methods for aligning raw fastq reads to a reference genome.
"""

import copy
from datetime import datetime
import os
import re
import string
import subprocess

from celery import task

from main.models import get_dataset_with_type
from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from main.models import ReferenceGenome
from main.model_utils import clean_filesystem_location
from main.s3 import project_files_needed
from pipeline.read_alignment_util import ensure_bwa_index
from scripts.jbrowse_util import add_bam_file_track
from scripts.jbrowse_util import add_bed_file_track
from settings import TOOLS_DIR
from settings import BASH_PATH
from scripts.util import titlecase_spaces


SAMTOOLS_BINARY = '%s/samtools/samtools' % TOOLS_DIR


@task
@project_files_needed
def align_with_bwa_mem(alignment_group, sample_alignment):
    """
    REPLACES OLD BWA PIPELINE USING ALN AND SAMPE/SAMSE
    Aligns a sample to a reference genome using the bwa tool.

    Args:
        alignment_group: AlignmentGroup that this alignment is part of.
        sample_alignment: ExperimentSampleToAlignment. The respective dataset
            is assumed to have been created as well.
    """

    # Start by gettng fresh objects from database.
    sample_alignment = ExperimentSampleToAlignment.objects.get(
            id=sample_alignment.id)
    experiment_sample = sample_alignment.experiment_sample
    alignment_group = AlignmentGroup.objects.get(id=alignment_group.id)


    # Grab the reference genome fasta for the alignment.
    ref_genome_fasta = get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # Get the BWA Dataset and set it to computing.
    bwa_dataset = sample_alignment.dataset_set.get(
                type=Dataset.TYPE.BWA_ALIGN)
    bwa_dataset.status = Dataset.STATUS.COMPUTING
    bwa_dataset.save(update_fields=['status'])

    # Create a file that we'll write stderr to.
    error_path = os.path.join(sample_alignment.get_model_data_dir(),
            'bwa_align.error')
    error_output = open(error_path, 'w')

    # The alignment group is now officially ALIGNING.
    if alignment_group.status != AlignmentGroup.STATUS.ALIGNING:
        alignment_group.status = AlignmentGroup.STATUS.ALIGNING
        alignment_group.start_time = datetime.now()
        alignment_group.end_time = None
        alignment_group.save(update_fields=['status','start_time','end_time'])

    error_output.write(
            "==START OF ALIGNMENT PIPELINE FOR %s, (%s) ==" % (
            sample_alignment.experiment_sample.label, sample_alignment.uid))

    # We wrap the alignment logic in a try-except so that if an error occurs,
    # we record it and update the status of the Dataset to FAILED if anything
    # should fail.
    try:
        # Build index if the index doesn't exist.
        # NOTE: When aligning multiple samples to the same reference genome
        # concurrently, the build index method should be called once to completion
        # before starting the concurrent alignment jobs.
        ensure_bwa_index(ref_genome_fasta)

        # Grab the fastq sources, and determine whether we are doing paired ends.

        for dataset in experiment_sample.dataset_set.all():
            if dataset.type == Dataset.TYPE.FASTQ1:
                input_reads_1_fq = dataset.wrap_if_compressed()
                input_reads_1_fq_path = dataset.get_absolute_location()
            elif dataset.type == Dataset.TYPE.FASTQ2:
                input_reads_2_fq = dataset.wrap_if_compressed()
                input_reads_2_fq_path = dataset.get_absolute_location()

        assert input_reads_1_fq, "Must have at least one .fastq file."
        is_paired_end = not input_reads_2_fq is None

        # There are two steps to alignment with bwa.
        #     1. Call 'bwa aln' to generate SA coordinate indices. Note that
        #        we run a separate bwa aln for each end, and do so in parallel
        #        when paired end reads.
        #     2. Call 'bwa sampe' (paired-end) or 'bwa samse' (single-end) to
        #        generate SAM output.

        # 1. Generate SA coordinates.
        read_fq_1_path, read_fq_1_fn = os.path.split(input_reads_1_fq_path)

        align_input_args = ' '.join([
            '%s/bwa/bwa' % TOOLS_DIR,
            'mem',
            '-t', '1', # threads
            '-M', # picard compatibility
            ref_genome_fasta,
            input_reads_1_fq,
        ])

        if is_paired_end:
            read_fq_2_path, read_fq_2_fn = os.path.split(input_reads_2_fq_path)
            align_input_args += ' ' + input_reads_2_fq

        # To skip saving the SAM file to disk directly, pipe output directly to
        # make a BAM file.
        align_input_args += ' | ' + SAMTOOLS_BINARY + ' view -bS -'

        ### 2. Generate SAM output.
        output_bam = os.path.join(sample_alignment.get_model_data_dir(),
                'bwa_align.bam')

        error_output.write(align_input_args)

        with open(output_bam, 'w') as fh:
            subprocess.check_call(align_input_args,
                    stdout=fh, stderr=error_output,
                    shell=True, executable=BASH_PATH)

        # Do several layers of processing on top of the initial alignment.
        result_bam_file = process_sam_bam_file(sample_alignment,
                alignment_group.reference_genome, output_bam, error_output)

        # Add the resulting file to the dataset.
        bwa_dataset.filesystem_location = clean_filesystem_location(
                result_bam_file)
        bwa_dataset.save()

        # Add track to JBrowse.
        add_bam_file_track(alignment_group.reference_genome, sample_alignment,
                Dataset.TYPE.BWA_ALIGN)

        bwa_dataset.status = Dataset.STATUS.READY
        bwa_dataset.save()

        delete_redundant_files(sample_alignment.get_model_data_dir())

    # except subprocess.CalledProcessError as e:
    #     error_output.write(str(e))
    #     bwa_dataset.status = Dataset.STATUS.FAILED
    #     bwa_dataset.save()
    #     return
    except:
        import traceback
        error_output.write(traceback.format_exc())
        bwa_dataset.status = Dataset.STATUS.FAILED
        bwa_dataset.save()
        return
    finally:
        print error_path
        error_output.write('==END OF ALIGNMENT PIPELINE==\n')
        error_output.close()

        # Add the error Dataset to the object.
        error_dataset = Dataset.objects.create(
                label=Dataset.TYPE.BWA_ALIGN_ERROR,
                type=Dataset.TYPE.BWA_ALIGN_ERROR,
                filesystem_location=clean_filesystem_location(error_path))
        sample_alignment.dataset_set.add(error_dataset)
        sample_alignment.save()

        return sample_alignment


DEFAULT_PROCESSING_MASK = {
    'make_bam': True,
    'sort': True,
    'add_groups': True,
    'indel_realigner': True,
    'compute_insert_metrics': True,
    'index': True,
    'compute_callable_loci': True,
    'withmd': True,
}


def process_sam_bam_file(sample_alignment, reference_genome,
        sam_bam_file_location, error_output=None,
        opt_processing_mask=DEFAULT_PROCESSING_MASK):
    """Converts to bam, sorts, and creates index.

    Args:
        sample_alignment: The relationship between a sample and an alignment
        sam_bam_file_location: The full path to the .sam/.bam alignment output.
            If sam input, will be converted from sam to bam.
        error_output: File handle that can be passed as stderr to subprocess
            calls.
        opt_processing_mask: An optional mask that can be modified to
            specify which processes to run. This is useful when doing a
            re-run of a partially-complete run and previous steps completed.
            NOTE: Not for amateurs.

    Returns:
        The path of the final .bam file.
    """

    experiment_sample = sample_alignment.experiment_sample

    # For any keys missing from the processing mask, give them the values from
    # the default mask.
    effective_mask = copy.copy(DEFAULT_PROCESSING_MASK)
    effective_mask.update(opt_processing_mask)

    #1. Convert to bam file if necessary
    if os.path.splitext(sam_bam_file_location)[1] == '.sam':
        sam_file_location = sam_bam_file_location
        bam_file_location = os.path.splitext(sam_file_location)[0] + '.bam'

        if effective_mask['make_bam']:
            with open(bam_file_location, 'w') as fh:
                subprocess.check_call([
                    SAMTOOLS_BINARY,
                    'view',
                    '-bS',
                    sam_file_location
                ], stdout=fh, stderr=error_output)
    else:
        bam_file_location = sam_bam_file_location

    # TODO: Make this more robust, or better yet, pipe things together.
    assert os.path.splitext(bam_file_location)[1] == '.bam'

    # 2. Sort
    sorted_output_name = os.path.splitext(bam_file_location)[0] + '.sorted'
    sorted_bam_file_location = sorted_output_name + '.bam'
    if effective_mask['sort']:
        # 2a. Perform the actual sorting.
        subprocess.check_call([
            SAMTOOLS_BINARY,
            'sort',
            bam_file_location,
            sorted_output_name
        ], stderr=error_output)

        # 2b. Index the sorted result.
        subprocess.check_call([
            SAMTOOLS_BINARY,
            'index',
            sorted_bam_file_location,
        ], stderr=error_output)

    # 3. Compute insert size metrics
    # Subsequent steps screw up pairing info so this has to
    # be done here.
    if effective_mask['compute_insert_metrics']:
        compute_insert_metrics(sorted_bam_file_location, sample_alignment, error_output)

    # 4. Add groups
    grouped_output_name= (
            os.path.splitext(sorted_bam_file_location)[0] +
            '.grouped')
    grouped_bam_file_location = grouped_output_name + '.bam'
    if effective_mask['add_groups']:
        add_groups(experiment_sample, sorted_bam_file_location, grouped_bam_file_location,
                error_output)

    # 5. Perform realignment accounting for indels.
    realigned_bam_file_location = (
            os.path.splitext(grouped_bam_file_location)[0] +
            '.realigned.bam'
    )
    if effective_mask['indel_realigner']:
        # Make sure the previous result is indexed.
        subprocess.check_call([
                SAMTOOLS_BINARY,
                'index',
                grouped_bam_file_location,
        ], stderr=error_output)

        realign_given_indels(
                experiment_sample,
                reference_genome,
                grouped_bam_file_location,
                realigned_bam_file_location,
                error_output
        )

    # 6. Add back MD tags for visualization of mismatches by Jbrowse
    if effective_mask['withmd']:
        final_bam_location = (
                os.path.splitext(realigned_bam_file_location)[0] +
                '.withmd.bam'
        )

        ref_genome_fasta_location = get_dataset_with_type(
                reference_genome,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

        with open(final_bam_location, 'w') as fh:
            # Add MD tags for Jbrowse visualization
            subprocess.check_call([
                SAMTOOLS_BINARY,
                'fillmd', '-b',
                realigned_bam_file_location,
                ref_genome_fasta_location
            ], stderr=error_output, stdout=fh)

        # Re-index this new bam file.
        subprocess.check_call([
            SAMTOOLS_BINARY, 'index', final_bam_location],
            stderr=error_output)

    else:
        final_bam_location = realigned_bam_file_location

    # 7. Compute callable loci
    if effective_mask['compute_callable_loci']:
        compute_callable_loci(reference_genome, sample_alignment,
                final_bam_location, error_output)

    # 8. Create index.
    if effective_mask['index']:
        subprocess.check_call([
            SAMTOOLS_BINARY,
            'index',
            final_bam_location,
        ], stderr=error_output)

    return final_bam_location


def add_groups(experiment_sample, input_bam_file, bam_output_file, error_output):
    from main.models import ensure_exists_0775_dir
    # Create a temp directory in the experiment_sample data dir else we run out of space
    # with big alignment jobs.
    picard_tmp_dir = os.path.join(experiment_sample.get_model_data_dir(), 'picard_tmp')
    ensure_exists_0775_dir(picard_tmp_dir)

    prefix = experiment_sample.uid
    subprocess.check_call([
        'java', '-Xmx1024M',
        '-jar', '%s/picard/AddOrReplaceReadGroups.jar' % TOOLS_DIR,
        'I=' + input_bam_file,
        'O=' + bam_output_file,
        'RGPU=' + prefix,
        'RGLB=' + prefix,
        'RGID=' + prefix,
        'RGPL=illumina',
        'RGSM=' + prefix,
        'SORT_ORDER=coordinate',
        'TMP_DIR=' + picard_tmp_dir, # Write temp data locally to avoid exhausting space.
        'VALIDATION_STRINGENCY=LENIENT' # Prevent unmapped read problems
    ], stderr=error_output)


def realign_given_indels(
        experiment_sample,
        reference_genome,
        input_bam_file,
        output_bam_file,
        error_output):
    """Perform realignment accounting for indels. This should give a cleaner
    alignment.
    """
    temp_intervals_file = os.path.join(
            experiment_sample.get_model_data_dir(), 'realign.intervals')
    ref_genome_fasta_location = get_dataset_with_type(
            reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # Prepare realignment intervals file.
    subprocess.check_call([
        'java', '-Xmx1024M',
        '-jar', '%s/gatk/GenomeAnalysisTK.jar' % TOOLS_DIR,
        '-T', 'RealignerTargetCreator',
        '-I', input_bam_file,
        '-R', ref_genome_fasta_location,
        '-o', temp_intervals_file,
    ], stderr=error_output)

    # Perform realignment.
    subprocess.check_call([
        'java', '-Xmx1024M',
        '-jar', '%s/gatk/GenomeAnalysisTK.jar' % TOOLS_DIR,
        '-T', 'IndelRealigner',
        '-I', input_bam_file,
        '-R', ref_genome_fasta_location,
        '-targetIntervals',  temp_intervals_file,
        '-o', output_bam_file,
    ], stderr=error_output)


def compute_insert_metrics(bam_file_location, sample_alignment, stderr=None):
    output = _get_metrics_output_filename(bam_file_location)
    histogram_file = os.path.splitext(bam_file_location)[0] + '.histmet.txt'
    subprocess.check_call([
        'java', '-Xmx1024M',
        '-jar', '%s/picard/CollectInsertSizeMetrics.jar' % TOOLS_DIR,
        'I=' + bam_file_location,
        'O=' + output,
        'HISTOGRAM_FILE=' + histogram_file,
        'VALIDATION_STRINGENCY=LENIENT' # Prevent unmapped read problems
    ], stderr=stderr)

    # Add insert metrics file as a dataset
    insert_metrics = Dataset.objects.create(
            label=Dataset.TYPE.PICARD_INSERT_METRICS,
            type=Dataset.TYPE.PICARD_INSERT_METRICS,
            filesystem_location=clean_filesystem_location(output))

    sample_alignment.dataset_set.add(insert_metrics)
    sample_alignment.save()


def compute_callable_loci(reference_genome, sample_alignment,
            bam_file_location, stderr=None):

    try:
        ref_genome_fasta_location = get_dataset_with_type(
                reference_genome,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

        output = _get_callable_loci_output_filename(bam_file_location)
        summary_file = os.path.splitext(bam_file_location)[0] + '.loci_summary.txt'

        subprocess.check_call([
            'java', '-Xmx1024M',
            '-jar', '%s/gatk/GenomeAnalysisTK.jar' % TOOLS_DIR,
            '-T', 'CallableLoci',
            '-R', ref_genome_fasta_location,
            '-I', bam_file_location,
            '-summary', summary_file,
            '-o', output,
        ], stderr=stderr)

        # Add callable loci bed as dataset
        callable_loci_bed = Dataset.objects.create(
                label=Dataset.TYPE.BED_CALLABLE_LOCI,
                type=Dataset.TYPE.BED_CALLABLE_LOCI,
                filesystem_location=clean_filesystem_location(output))

        sample_alignment.dataset_set.add(callable_loci_bed)
        sample_alignment.save()

        callable_loci_bed_fn = callable_loci_bed.get_absolute_location()

        # Remove 'CALLABLE' rows - we assume no feature means callable
        # Also convert other features to title case

        output = subprocess.check_output(
                ['grep',  '-v', 'CALLABLE', callable_loci_bed_fn])

        with open(callable_loci_bed_fn, 'w') as callable_loci_bed_fh:
            for i, line in enumerate(output.split('\n')):
                try:
                    fields = line.split()
                    if len(fields) == 0: continue
                    chrom, start, end, feature = fields
                    if feature == 'CALLABLE':
                        continue
                    else:

                        feature = titlecase_spaces(feature)
                        # Bed feature can't have spaces =(
                        feature = feature.replace(' ','_')

                    print >> callable_loci_bed_fh, '\t'.join(
                            [chrom, start, end, feature])
                except Exception as e:
                   print >> stderr, (
                            'WARNING: Callable Loci line' +
                            '%d: (%s) couldn\'t be parsed: %s') % (
                                    i, line, str(e))
        # add it as a jbrowse track
        add_bed_file_track(reference_genome, sample_alignment, callable_loci_bed)

    except Exception as e:
        print >> stderr, 'WARNING: Callable Loci failed.'
        print >> stderr, str(e)


def _get_metrics_output_filename(bam_file_location):
    return os.path.splitext(bam_file_location)[0] + '.insertmet.txt'


def _get_callable_loci_output_filename(bam_file_location):
    return os.path.splitext(bam_file_location)[0] + '.callable_loci.bed'


def get_insert_size(sample_alignment):
    """Returns the average insert size for a bam_file.

    If the insert size can't be calculated, perhaps because of a bad alignment,
    return -1.
    """
    picard_metrics = get_dataset_with_type(sample_alignment,
            Dataset.TYPE.PICARD_INSERT_METRICS)

    if picard_metrics is None:
        bam_file = get_dataset_with_type(sample_alignment,
                Dataset.TYPE.BWA_ALIGN).get_absolute_location()
        compute_insert_metrics(bam_file, sample_alignment, stderr=None)
        picard_metrics = get_dataset_with_type(sample_alignment,
            Dataset.TYPE.PICARD_INSERT_METRICS)

    try:
        metric_filename = picard_metrics.get_absolute_location()
        assert os.path.exists(metric_filename)
    except:
        # if it doesn't work, return -1
        return -1

    # Read through the output file until you get to the data labels line,
    # then take the first value in the following line.
    with open(metric_filename) as fh:
        found_line = False
        insert_size = -1
        for line in fh:
            if found_line:
                insert_size = line.split()[0]
                break
            if re.match('MEDIAN_INSERT_SIZE', line):
                found_line = True
    return insert_size


##############################################################################
# Clean-ups
# TODO: This doesn't work cleanly with our mask strategy. Figure this out
#      when we implement support for users manipulating alignment options.
##############################################################################

FILES_TO_DELETE_AFTER_ALIGNMENT = set([
    'bwa_align.sam',
    'bwa_align.bam',
    'bwa_align.sorted.bam',
    'bwa_align.sorted.bam.bai',
    'bwa_align.sorted.grouped.bam',
    'bwa_align.sorted.grouped.bam.bai',
    'bwa_align.sorted.grouped.realigned.bam'
])


def delete_redundant_files(source_dir):
    """Delete redundant files upon alignment completion.

    This is sort of a hack for now to deal with the issue where we are
    generating a bunch of extra data during alignment.
    """
    for f in os.listdir(source_dir):
        if f in FILES_TO_DELETE_AFTER_ALIGNMENT:
            os.remove(os.path.join(source_dir, f))


def clean_alignment_group_data(alignment_group):
    """Utility method for cleaning redundant data manually.
    """
    for sample_alignment in alignment_group.experimentsampletoalignment_set.all():
        maybe_dataset_set = sample_alignment.dataset_set.filter(
                type=Dataset.TYPE.BWA_ALIGN)
        if len(maybe_dataset_set):
            bwa_dataset = maybe_dataset_set[0]
            if bwa_dataset.status == Dataset.STATUS.READY:
                delete_redundant_files(sample_alignment.get_model_data_dir())
