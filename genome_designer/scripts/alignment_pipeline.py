
"""
The alignment pipeline.

We start with the .fastq files and the reference for a particular genome
and carry out the gauntlet of steps to perform alignments and clean them up.
"""

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
from scripts.util import fn_runner
from settings import DEBUG_CONCURRENT
from settings import PWD
from settings import TOOLS_DIR
from settings import BASH_PATH

def create_alignment_groups_and_start_alignments(
        alignment_group_label, ref_genome_list, sample_list,
        test_models_only=False, concurrent=DEBUG_CONCURRENT):

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
    if concurrent:
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
        alignment_tasks = []
        for sample in sample_list:
            args = [alignment_group, sample, None, test_models_only]
            alignment_tasks.append(fn_runner(align_with_bwa_mem,
                alignment_group.reference_genome.project, args, concurrent))

            #alignment_tasks.append(fn_runner(align_with_bwa, args, concurrent))
            # TODO(DBG): Uncomment once bwa_mem works correctly.
            # See: https://github.com/churchlab/genome-designer-v2/issues/93
            # alignment_tasks.append(fn_runner(align_with_bwa_mem, args, concurrent))

    # Return a dictionary of all alignment groups indexed by ref_genome uid.
    return alignment_groups
    

def align_with_bwa_mem(alignment_group, experiment_sample=None,
        sample_alignment=None, test_models_only=False):
    """
    REPLACES OLD BWA PIPELINE USING ALN AND SAMPE/SAMSE
    Aligns a sample to a reference genome using the bwa tool.

    Args:
        alignment_group: AlignmentGroup that this alignment is part of.
        experiment_sample: ExperimentSample with handle to source fastq, etc.
            If not provided, then sample_alignment must be provided.
        sample_alignment: If provided, delete any previously contained
            alignment data and re-run alignment.
        test_models_only: If True, don't actually perform alignment, just
            create the basic models.
    """

    ### Validation
    if not experiment_sample:
        assert sample_alignment

    # Now continue with alignment.

    # First we check whether we are re-running a previously failed alignment.
    # In other words, if an alignment instance was passed as an argument.
    if sample_alignment:
        # TODO: Delete existing data?
        experiment_sample = sample_alignment.experiment_sample
    else:
        # Create the initial record.
        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=experiment_sample)

    # If we are only testing that the right models are created, then
    # we return at this point. This should never be True in production.
    if test_models_only:
        return

    # Grab the reference genome fasta for the alignment.
    ref_genome_fasta = get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # Write error to this file.
    error_path = os.path.join(sample_alignment.get_model_data_dir(),
            'bwa_align.error')
    error_output = open(error_path, 'w')

    # Create a Dataset object and set the state to COMPUTING.
    bwa_dataset = Dataset.objects.create(
            label=Dataset.TYPE.BWA_ALIGN,
            type=Dataset.TYPE.BWA_ALIGN,
            status=Dataset.STATUS.COMPUTING)
    sample_alignment.dataset_set.add(bwa_dataset)
    sample_alignment.save()

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
            ref_genome_fasta,
            input_reads_1_fq,
        ])

        if is_paired_end:
            read_fq_2_path, read_fq_2_fn = os.path.split(input_reads_2_fq_path)
            align_input_args += ' ' + input_reads_2_fq

        # To skip saving the SAM file to disk directly, pipe output directly to 
        # make a BAM file.
        align_input_args += ' | samtools view -bS -'

        # 2. Generate SAM output.
        output_bam = os.path.join(sample_alignment.get_model_data_dir(),
                'bwa_align.bam')

        with open(output_bam, 'w') as fh:
            subprocess.check_call(align_input_args, 
                    stdout=fh, stderr=error_output, 
                    shell=True, executable=BASH_PATH)

        # Do several layers of processing on top of the initial alignment.
        result_bam_file = process_sam_bam_file(experiment_sample,
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

    except subprocess.CalledProcessError as e:
        error_output.write(str(e))
        bwa_dataset.status = Dataset.STATUS.FAILED
        bwa_dataset.save()
        return
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


def align_with_bwa(alignment_group, experiment_sample=None,
        sample_alignment=None, test_models_only=False):
    """DEPRECATED: NOT USED BY CURRENT PIPELINE
    Aligns a sample to a reference genome using the bwa tool.

    Args:
        alignment_group: AlignmentGroup that this alignment is part of.
        experiment_sample: ExperimentSample with handle to source fastq, etc.
            If not provided, then sample_alignment must be provided.
        sample_alignment: If provided, delete any previously contained
            alignment data and re-run alignment.
        test_models_only: If True, don't actually perform alignment, just
            create the basic models.
    """
    ### Validation
    if not experiment_sample:
        assert sample_alignment

    # Now continue with alignment.

    # First we check whether we are re-running a previously failed alignment.
    # In other words, if an alignment instance was passed as an argument.
    if sample_alignment:
        # TODO: Delete existing data?
        experiment_sample = sample_alignment.experiment_sample
    else:
        # Create the initial record.
        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=experiment_sample)

    # If we are only testing that the right models are created, then
    # we return at this point. This should never be True in production.
    if test_models_only:
        return

    # Grab the reference genome fasta for the alignment.
    ref_genome_fasta = get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # Write error to this file.
    error_path = os.path.join(sample_alignment.get_model_data_dir(),
            'bwa_align.error')
    error_output = open(error_path, 'w')

    # Create a Dataset object and set the state to COMPUTING.
    bwa_dataset = Dataset.objects.create(
            label=Dataset.TYPE.BWA_ALIGN,
            type=Dataset.TYPE.BWA_ALIGN,
            status=Dataset.STATUS.COMPUTING)
    sample_alignment.dataset_set.add(bwa_dataset)
    sample_alignment.save()

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
        output_index_1 = os.path.join(
                sample_alignment.get_model_data_dir(),
                (os.path.splitext(read_fq_1_fn)[0] + '.sai'))
        print output_index_1

        with open(output_index_1, 'w') as input_reads_1_stdout_fh:
            align_input_1_args = ' '.join([
                '%s/bwa/bwa' % TOOLS_DIR,
                'aln',
                '-t', '1', # threads
                ref_genome_fasta,
                input_reads_1_fq,
            ])

            print align_input_1_args

            # NOTE: Non-blocking. See call to wait_and_check_pipe() below.
            align_input_1_proc = subprocess.Popen(align_input_1_args,
                    stdout=input_reads_1_stdout_fh, stderr=error_output,
                    shell=True, executable=BASH_PATH)

        if is_paired_end:
            read_fq_2_path, read_fq_2_fn = os.path.split(input_reads_2_fq_path)
            output_index_2 = os.path.join(
                    sample_alignment.get_model_data_dir(),
                    (os.path.splitext(read_fq_2_fn)[0] + '.sai'))

            with open(output_index_2, 'w') as input_reads_2_stdout_fh:
                align_input_2_args = ' '.join([
                    '%s/bwa/bwa' % TOOLS_DIR,
                    'aln',
                    '-t', '1', # threads
                    ref_genome_fasta,
                    input_reads_2_fq,
                ])
                align_input_2_proc = subprocess.Popen(align_input_2_args,
                        stdout=input_reads_2_stdout_fh, stderr=error_output,
                        shell=True, executable=BASH_PATH)
                wait_and_check_pipe(align_input_2_proc, align_input_2_args)

                print align_input_2_args

        # Block here until both calls to 'bwa aln' complete.
        wait_and_check_pipe(align_input_1_proc, align_input_1_args)

        # 2. Generate SAM output.
        output_sam = os.path.join(sample_alignment.get_model_data_dir(),
                'bwa_align.sam')

        if is_paired_end:
            with open(output_sam, 'w') as fh:
                subprocess.check_call(' '.join([
                    '%s/bwa/bwa' % TOOLS_DIR,
                    'sampe',
                    ref_genome_fasta,
                    output_index_1,
                    output_index_2,
                    input_reads_1_fq,
                    input_reads_2_fq,
                ]), stdout=fh, stderr=error_output,
                shell=True, executable=BASH_PATH)
        else:
            with open(output_sam, 'w') as fh:
                subprocess.check_call(' '.join([
                    '%s/bwa/bwa' % TOOLS_DIR,
                    'samse',
                    ref_genome_fasta,
                    output_index_1,
                    input_reads_1_fq,
                ]), stdout=fh, stderr=error_output,
                shell=True, executable=BASH_PATH)

        # Do several layers of processing on top of the initial alignment.
        result_bam_file = process_sam_bam_file(experiment_sample,
                alignment_group.reference_genome, output_sam, error_output)

        # Add the resulting file to the dataset.
        bwa_dataset.filesystem_location = clean_filesystem_location(
                result_bam_file)
        bwa_dataset.save()

        # Add track to JBrowse.
        add_bam_file_track(alignment_group.reference_genome, sample_alignment,
                Dataset.TYPE.BWA_ALIGN)

        bwa_dataset.status = Dataset.STATUS.READY
        bwa_dataset.save()

    except subprocess.CalledProcessError as e:
        error_output.write(str(e))
        bwa_dataset.status = Dataset.STATUS.FAILED
        bwa_dataset.save()
        return
    except:
        import traceback
        error_output.write(traceback.format_exc())
        bwa_dataset.status = Dataset.STATUS.FAILED
        bwa_dataset.save()
        return
    finally:
        print error_path
        error_output.write('test\n')
        error_output.close()

        # Add the error Dataset to the object.
        error_dataset = Dataset.objects.create(
                label=Dataset.TYPE.BWA_ALIGN_ERROR,
                type=Dataset.TYPE.BWA_ALIGN_ERROR,
                filesystem_location=clean_filesystem_location(error_path))
        sample_alignment.dataset_set.add(error_dataset)
        sample_alignment.save()

        return sample_alignment


def ensure_bwa_index(ref_genome_fasta, error_output=None):
    """Creates the bwa index if it doesn't exist already.

    We rely on the convention that the index file location is the fasta
    location with the extension '.bwt' appended to it.
    """
    if not os.path.exists(ref_genome_fasta + '.bwt'):
        build_bwa_index(ref_genome_fasta, error_output)

    # Also build the fasta index.
    if not os.path.exists(ref_genome_fasta + '.fai'):
        subprocess.check_call([
            '%s/samtools/samtools' % TOOLS_DIR,
            'faidx',
            ref_genome_fasta
        ], stderr=error_output)

    ref_genome_dict_path = os.path.splitext(ref_genome_fasta)[0] + '.dict'
    if not os.path.exists(ref_genome_dict_path):
        subprocess.check_call([
            'java', '-Xmx1024M',
            '-jar', '%s/picard/CreateSequenceDictionary.jar' % TOOLS_DIR,
            'R=', ref_genome_fasta,
            'O=', ref_genome_dict_path,
        ], stderr=error_output)


def build_bwa_index(ref_genome_fasta, error_output=None):
    """Calls the command that builds the bwa index required for alignment.

    This creates a file in the same directory as ref_genome_fasta, appending
    the extension '.bwt' to the name of the fasta.
    """
    subprocess.check_call([
        '%s/bwa/bwa' % TOOLS_DIR,
        'index',
        '-a',
        'is',
        ref_genome_fasta
    ], stderr=error_output)


DEFAULT_PROCESSING_MASK = {
    'make_bam': True,
    'sort': True,
    'add_groups': True,
    'indel_realigner': True,
    'compute_insert_metrics': True,
    'index': True,
}

def process_sam_bam_file(experiment_sample, reference_genome, sam_bam_file_location,
        error_output=None, opt_processing_mask=DEFAULT_PROCESSING_MASK):
    """Converts to bam, sorts, and creates index.

    Args:
        experiment_sample: The ExperimentSample for this alignment.
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
    # For any keys missing from the processing mask, give them the values from
    # the default mask.
    effective_mask = copy.copy(DEFAULT_PROCESSING_MASK)
    effective_mask.update(opt_processing_mask)

    #1. Convert to bam file if necessary
    if os.path.splitext(sam_bam_file_location)[1] == '.sam':
        sam_file_location = sam_bam_file_location
        bam_file_location = os.path.splitext(sam_file_location)[0] + '.bam'

        if effective_mask['make_bam']:
            fh = open(bam_file_location, 'w')
            subprocess.check_call([
                '%s/samtools/samtools' % TOOLS_DIR,
                'view',
                '-bS',
                sam_file_location
            ], stdout=fh, stderr=error_output)
            fh.close()
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
            '%s/samtools/samtools' % TOOLS_DIR,
            'sort',
            bam_file_location,
            sorted_output_name
        ], stderr=error_output)

        # 2b. Index the sorted result.
        subprocess.check_call([
            '%s/samtools/samtools' % TOOLS_DIR,
            'index',
            sorted_bam_file_location,
        ], stderr=error_output)

    # 3. Add groups
    grouped_output_name= (
            os.path.splitext(sorted_bam_file_location)[0] +
            '.grouped')
    grouped_bam_file_location = grouped_output_name + '.bam'
    if effective_mask['add_groups']:
        add_groups(experiment_sample, sorted_bam_file_location, grouped_bam_file_location,
                error_output)

    # 3. Perform realignment accounting for indels.
    realigned_bam_file_location = (
            os.path.splitext(grouped_bam_file_location)[0] +
            '.realigned.bam'
    )
    if effective_mask['indel_realigner']:
        # Make sure the previous result is indexed.
        subprocess.check_call([
                '%s/samtools/samtools' % TOOLS_DIR,
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

    # 4. Compute insert size metrics
    if effective_mask['compute_insert_metrics']:
        compute_insert_metrics(realigned_bam_file_location, error_output)

    # 5. Create index.
    if effective_mask['index']:
        subprocess.check_call([
            '%s/samtools/samtools' % TOOLS_DIR,
            'index',
            realigned_bam_file_location,
        ], stderr=error_output)

    return realigned_bam_file_location


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


def compute_insert_metrics(bam_file_location, stderr=None):
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


def _get_metrics_output_filename(bam_file_location):
    return os.path.splitext(bam_file_location)[0] + '.insertmet.txt'


def get_insert_size(bam_file_location):
    """Returns the average insert size for a bam_file."""
    output = _get_metrics_output_filename(bam_file_location)
    if not os.path.exists(output):
        compute_insert_metrics(bam_file_location)

    # Read through the output file until you get to the data labels line,
    # then take the first value in the following line.
    fh = open(output)
    found_line = False
    insert_size = -1
    for line in fh:
        if found_line:
            insert_size = line.split()[0]
            break
        if re.match('MEDIAN_INSERT_SIZE', line):
            found_line = True
    fh.close()
    return insert_size


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
