"""
Methods for aligning raw fastq reads to a reference genome.
"""

import copy
from datetime import datetime
import os
import subprocess
from subprocess import PIPE
from subprocess import Popen

from celery import task
from django.conf import settings

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSampleToAlignment
from main.model_utils import clean_filesystem_location
from main.model_utils import get_dataset_with_type
from main.s3 import project_files_needed
from pipeline.read_alignment_util import ensure_bwa_index
from pipeline.callable_loci import get_callable_loci
from pipeline.read_alignment_util import index_bam_file
from pipeline.read_alignment_util import extract_split_reads
from pipeline.read_alignment_util import extract_discordant_read_pairs
from utils.import_util import add_dataset_to_entity
from utils.jbrowse_util import add_bam_file_track
from utils.jbrowse_util import add_bed_file_track
from utils import titlecase_spaces


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
        alignment_group.save(update_fields=['status', 'start_time', 'end_time'])

    error_output.write(
            "==START OF ALIGNMENT PIPELINE FOR %s, (%s) ==\n" % (
                    sample_alignment.experiment_sample.label,
                    sample_alignment.uid))

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
        # First, grab fastq1, which must exist
        fq1_queryset = experiment_sample.dataset_set.filter(
                type=Dataset.TYPE.FASTQ1)
        assert fq1_queryset, "Must have at least one .fastq file"
        fq1_dataset = fq1_queryset[0]
        input_reads_1_fq = fq1_dataset.wrap_if_compressed()
        input_reads_1_fq_path = fq1_dataset.get_absolute_location()

        # Second, check if fastq2 exists and set is_paired_end
        fq2_queryset = experiment_sample.dataset_set.filter(
                type=Dataset.TYPE.FASTQ2)
        if fq2_queryset:
            is_paired_end = True
            fq2_dataset = fq2_queryset[0]
            input_reads_2_fq = fq2_dataset.wrap_if_compressed()
            input_reads_2_fq_path = fq2_dataset.get_absolute_location()
        else:
            is_paired_end = False

        # 1. Generate SA coordinates.
        read_fq_1_path, read_fq_1_fn = os.path.split(input_reads_1_fq_path)

        align_input_args = ' '.join([
            '%s/bwa/bwa' % settings.TOOLS_DIR,
            'mem',
            '-t', '1', # threads
            '-R', '"'+read_group_string(experiment_sample)+'"',
            # uncomment this to keep secondary alignments (for finding and marking paralogy regions)
            # But before we can uncomment we need to fix de novo assembly code
            '-a',
            ref_genome_fasta,
            input_reads_1_fq,
        ])

        if is_paired_end:
            read_fq_2_path, read_fq_2_fn = os.path.split(input_reads_2_fq_path)
            align_input_args += ' ' + input_reads_2_fq

        # To skip saving the SAM file to disk directly, pipe output directly to
        # make a BAM file.
        align_input_args += ' | ' + settings.SAMTOOLS_BINARY + ' view -bS -'

        ### 2. Generate SAM output.
        output_bam = os.path.join(sample_alignment.get_model_data_dir(),
                'bwa_align.bam')

        error_output.write(align_input_args)

        # Flush the output here so it gets written before the alignments.
        error_output.flush()

        with open(output_bam, 'w') as fh:
            subprocess.check_call(align_input_args,
                    stdout=fh, stderr=error_output,
                    shell=True, executable=settings.BASH_PATH)

        # Set processing mask to not compute insert metrics if reads are
        # not paired end, as the lumpy script only works on paired end reads
        opt_processing_mask = {}
        if not is_paired_end:
            opt_processing_mask['compute_insert_metrics'] = False

        # Do several layers of processing on top of the initial alignment.
        result_bam_file = process_sam_bam_file(sample_alignment,
                alignment_group.reference_genome, output_bam, error_output,
                opt_processing_mask=opt_processing_mask)

        # Add the resulting file to the dataset.
        bwa_dataset.filesystem_location = clean_filesystem_location(
                result_bam_file)
        bwa_dataset.save()

        # Isolate split and discordant reads for SV calling.
        get_discordant_read_pairs(sample_alignment)
        get_split_reads(sample_alignment)

        # Add track to JBrowse.
        add_bam_file_track(alignment_group.reference_genome, sample_alignment,
                Dataset.TYPE.BWA_ALIGN)

        bwa_dataset.status = Dataset.STATUS.READY
        bwa_dataset.save()

        delete_redundant_files(sample_alignment.get_model_data_dir())

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
    'rmdup': True,
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
                    settings.SAMTOOLS_BINARY,
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

    # We are only performing rmdup if we are sorting - if no sorting, then
    # skip this even if it is flagged. We probably should expose this to
    # the user, but for now let's just throw an error.
    if not effective_mask['sort'] and not effective_mask['rmdup']:
        raise Exception('Cannot remove duplicates without sorting!'
                ' Contents of processing mask: \n'+str(effective_mask))

    elif effective_mask['sort'] and effective_mask['rmdup']:
        # 2a. Perform the actual sorting and rmdup on stream.
        sort_rmdup_cmd = '|'.join([
                ' '.join([
                        settings.SAMTOOLS_BINARY,
                        'sort',
                        '-o',
                        bam_file_location,
                        sorted_output_name + '.tmp.bam']),
                ' '.join([
                        settings.SAMTOOLS_BINARY,
                        'rmdup',
                        '-',
                        sorted_bam_file_location])])

        subprocess.check_call(sort_rmdup_cmd, shell=True, stderr=error_output)

        # 2b. Index the sorted result.
        index_bam_file(sorted_bam_file_location, error_output)

    elif effective_mask['sort']:
        # 2a. Perform the actual sorting.
        subprocess.check_call([
            settings.SAMTOOLS_BINARY,
            'sort',
            bam_file_location,
            sorted_output_name
        ], stderr=error_output)

        # 2b. Index the sorted result.
        index_bam_file(sorted_bam_file_location, error_output)

    # 3. Compute insert size metrics
    # Subsequent steps screw up pairing info so this has to
    # be done here.
    if effective_mask['compute_insert_metrics']:
        compute_insert_metrics(sorted_bam_file_location,
                sample_alignment, error_output)

    # 4. Add back MD tags for visualization of mismatches by Jbrowse
    if effective_mask['withmd']:
        final_bam_location = (
                os.path.splitext(sorted_bam_file_location)[0] +
                '.withmd.bam'
        )

        ref_genome_fasta_location = get_dataset_with_type(
                reference_genome,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

        with open(final_bam_location, 'w') as fh:
            # Add MD tags for Jbrowse visualization
            subprocess.check_call([
                settings.SAMTOOLS_BINARY,
                'fillmd', '-b',
                sorted_bam_file_location,
                ref_genome_fasta_location
            ], stderr=error_output, stdout=fh)

        # Re-index this new bam file.
        index_bam_file(final_bam_location, error_output)

    else:
        final_bam_location = sorted_bam_file_location

    # 5. Compute callable loci
    if effective_mask['compute_callable_loci']:
        compute_callable_loci(reference_genome, sample_alignment,
                final_bam_location, error_output)

    # 6. Create index.
    if effective_mask['index']:
        index_bam_file(final_bam_location, error_output)

    return final_bam_location


def read_group_string(experiment_sample):
    """
    Generate a SAM header string for the read group.
    Passed to bwa with escaped tabs.
    """
    read_group_fields = [
            '@RG',
            'ID:'+experiment_sample.uid,
            'PL:illumina',
            'PU:'+experiment_sample.uid,
            'LB:'+experiment_sample.uid,
            'SM:'+experiment_sample.uid]

    return('\\t'.join(read_group_fields))


def compute_insert_metrics(bam_file, sample_alignment, stderr=None):
    """Computes read fragment insert size distribution.

    Creates a Dataset for each of:
        * histogram file
        * file with mean and stdev comma-separated

    Raises:
        ValueError if calculating paired-end distribution failed.
    """
    histo_file = os.path.splitext(bam_file)[0] + '.insert_size_histogram.txt'
    mean_stdev_file = (os.path.splitext(bam_file)[0] +
            '.insert_size_mean_stdev.txt')

    # First, we analyze the bam distribution.
    read_bam_cmd = [
            settings.SAMTOOLS_BINARY,
            'view',
            bam_file
    ]
    p1 = Popen(read_bam_cmd, stdout=PIPE, stderr=stderr)

    read_length = get_read_length(bam_file)

    pairend_distro_cmd = [
        settings.LUMPY_PAIREND_DISTRO_BIN,
        '-r', str(read_length),
        '-X', '4', # num stdevs from end to extend
        '-N', '10000', # number to sample
        '-o', histo_file
    ]
    p2 = Popen(pairend_distro_cmd, stdin=p1.stdout, stdout=PIPE, stderr=stderr)

    # Allow p1 to receive a SIGPIPE if p2 exits.
    p1.stdout.close()

    # Run the command and get mean, stdev
    mean_and_stdev_str = p2.communicate()[0]
    mean_and_stdev_parts = mean_and_stdev_str.split('\t')
    if len(mean_and_stdev_parts) != 2:
        raise ValueError(
            "Poor alignment. Perhaps you tried aligning to the wrong reference "
            "genome?")
    raw_mean, raw_stdev = mean_and_stdev_parts
    mean = int(float(raw_mean.split(':')[1].strip()))
    stdev = int(float(raw_stdev.split(':')[1].strip()))

    # Lumpy doesn't like stdev of 0.
    if stdev < 1:
        stdev = 1

    # Save the histogram file as a Dataset.
    add_dataset_to_entity(sample_alignment,
            Dataset.TYPE.LUMPY_INSERT_METRICS_HISTOGRAM,
            Dataset.TYPE.LUMPY_INSERT_METRICS_HISTOGRAM,
            filesystem_location=histo_file)

    # Write mean, stdev to another file and create another Dataset.
    with open(mean_stdev_file, 'w') as fh:
        fh.write("%d,%d" % (mean, stdev))
    add_dataset_to_entity(sample_alignment,
            Dataset.TYPE.LUMPY_INSERT_METRICS_MEAN_STDEV,
            Dataset.TYPE.LUMPY_INSERT_METRICS_MEAN_STDEV,
            filesystem_location=mean_stdev_file)


def compute_callable_loci(reference_genome, sample_alignment,
            bam_file_location, stderr=None):

    # Set output fn to None in case try fails.
    callable_loci_bed_fn = None

    try:
        ref_genome_fasta_location = get_dataset_with_type(
                reference_genome,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

        callable_loci_bed_fn = (
                _get_callable_loci_output_filename(bam_file_location))

        get_callable_loci(bam_file_location, callable_loci_bed_fn)

        # Add callable loci bed as dataset
        callable_loci_bed_dataset = Dataset.objects.create(
                label=Dataset.TYPE.BED_CALLABLE_LOCI,
                type=Dataset.TYPE.BED_CALLABLE_LOCI,
                filesystem_location=clean_filesystem_location(callable_loci_bed_fn))

        sample_alignment.dataset_set.add(callable_loci_bed_dataset)
        sample_alignment.save()

        clean_bed_fn = clean_bed_features(callable_loci_bed_dataset, stderr=stderr)

        # add it as a jbrowse track
        add_bed_file_track(
                reference_genome,
                sample_alignment,
                callable_loci_bed_dataset)


    except Exception as e:
        print >> stderr, 'WARNING: Callable Loci failed.'
        print >> stderr, str(e)
        clean_bed_fn = ''

    finally:
        return clean_bed_fn


def clean_bed_features(callable_loci_bed_dataset, stderr=None):
    """
    This helper function updates the callable loci bed to make the
    features titlecase and remove spaces. It reads the whole bed into
    memory, cleans it, and prints it back out to the original file.

    It does NOT update the jbrowse track, which needs to be done separately
    with:
        add_bed_file_track(
                reference_genome, sample_alignment, callable_loci_bed)
    """
    # update the fn with the clean filesystem location
    callable_loci_bed_fn = callable_loci_bed_dataset.get_absolute_location()

    try:
        # store the bed in memory so we can clean/update it
        output = subprocess.check_output(
                ['cat', callable_loci_bed_fn])

        # write the bed again from memory, cleaning fields
        with open(callable_loci_bed_fn, 'w') as callable_loci_bed_fh:
            for i, line in enumerate(output.split('\n')):
                try:
                    fields = line.split()
                    if len(fields) == 0:
                        continue
                    chrom, start, end, feature = fields

                    feature = titlecase_spaces(feature)
                    # Bed feature can't have spaces =(
                    feature = feature.replace(' ', '_')

                    print >> callable_loci_bed_fh, '\t'.join(
                            [chrom, start, end, feature])
                except Exception as e:
                    print >> stderr, (
                        'WARNING: Callable Loci line' +
                        '%d: (%s) couldn\'t be parsed: %s') % (
                                i, line, str(e))

    except Exception as e:
        print >> stderr, 'WARNING: Callable Loci clean failed.'
        print >> stderr, str(e)
    finally:
        return callable_loci_bed_fn


def _get_metrics_output_filename(bam_file_location):
    return os.path.splitext(bam_file_location)[0] + '.insertmet.txt'


def _get_callable_loci_output_filename(bam_file_location):
    return os.path.splitext(bam_file_location)[0] + '.callable_loci.bed'


def get_read_length(bam_file):
    """Returns the median read length for this sample alignment by looking
    at the, at most, first 1000 reads.
    """
    assert len(os.path.splitext(bam_file)[1]), "Invalid bam: %s" % bam_file

    p = subprocess.Popen([settings.SAMTOOLS_BINARY, 'view', bam_file],
        stdout=subprocess.PIPE)

    lines = 0
    base_sum = 0

    for line in p.stdout:
        try:
            base_sum += len(line.split('\t')[9])
        except:
            break
        lines += 1
        if lines > 1000:
            break

    p.stdout.close()

    return int(base_sum / lines)


def get_insert_size_mean_and_stdev(sample_alignment, stderr=None, _iteration=0):
    """Returns a tuple (mean, stdev) for insert sizes from the alignment.

    Calls the compute functoin if metrics don't exist.

    If the insert size can't be calculated, perhaps because of a bad alignment,
    returns (-1, -1).

    Args:
        sample_alignment: ExperimentSampleToAlignment we want metrics for.
        iteration: Used internally to avoid getting stuck in case where
            computation repeatedly fails.

    Returns:
        Tuple of ints (mean, stdev).
    """
    # Prevent getting stuck in case computation keeps failing.
    if _iteration >= 3:
        return (-1, -1)

    def _compute():
        """Calls compute function then recursively calls
        get_insert_size_mean_and_stdev().
        """
        bam_file = get_dataset_with_type(sample_alignment,
                Dataset.TYPE.BWA_ALIGN).get_absolute_location()
        compute_insert_metrics(bam_file, sample_alignment, stderr=stderr)
        return get_insert_size_mean_and_stdev(sample_alignment, stderr,
                _iteration=_iteration + 1)

    mean_stdev_dataset = get_dataset_with_type(sample_alignment,
            Dataset.TYPE.LUMPY_INSERT_METRICS_MEAN_STDEV)
    if not mean_stdev_dataset:
        return _compute()

    file_location = mean_stdev_dataset.get_absolute_location()
    if not os.path.exists(file_location):
        return _compute()

    with open(file_location) as fh:
        combined_str = fh.read().strip()
        parts = combined_str.split(',')
        if not len(parts) == 2:
            return _compute()
        else:
            return tuple([int(p) for p in parts])


def get_discordant_read_pairs(sample_alignment):
    """Isolate discordant pairs of reads from a sample alignment.
    """
    # First, check if completed dataset already exists.
    bwa_disc_dataset = get_dataset_with_type(
            sample_alignment, Dataset.TYPE.BWA_DISCORDANT)
    if bwa_disc_dataset is not None:
        if (bwa_disc_dataset.status == Dataset.STATUS.READY and
                os.path.exists(bwa_disc_dataset.get_absolute_location())):
            return bwa_disc_dataset
    else:
        bwa_disc_dataset = Dataset.objects.create(
                label=Dataset.TYPE.BWA_DISCORDANT,
                type=Dataset.TYPE.BWA_DISCORDANT)
        sample_alignment.dataset_set.add(bwa_disc_dataset)

    # If here, we are going to run or re-run the Dataset.
    bwa_disc_dataset.status = Dataset.STATUS.NOT_STARTED
    bwa_disc_dataset.save(update_fields=['status'])

    bam_dataset = get_dataset_with_type(sample_alignment, Dataset.TYPE.BWA_ALIGN)
    bam_filename = bam_dataset.get_absolute_location()

    assert os.path.exists(bam_filename), "BAM file '%s' is missing." % (
            bam_filename)

    # NOTE: This assumes the index just adds at .bai, w/ same path otherwise
    # - will this always be true?
    if not os.path.exists(bam_filename+'.bai'):
        index_bam_file(bam_filename)

    bam_discordant_filename = os.path.join(sample_alignment.get_model_data_dir(),
            'bwa_discordant_pairs.bam')

    try:
        bwa_disc_dataset.status = Dataset.STATUS.COMPUTING
        bwa_disc_dataset.save(update_fields=['status'])
        extract_discordant_read_pairs(bam_filename, bam_discordant_filename)

    except subprocess.CalledProcessError:
        bwa_disc_dataset.filesystem_location = ''
        bwa_disc_dataset.status = Dataset.STATUS.FAILED
    finally:
        bwa_disc_dataset.status = Dataset.STATUS.READY
        bwa_disc_dataset.filesystem_location = clean_filesystem_location(
                bam_discordant_filename)

    bwa_disc_dataset.save()

    return bwa_disc_dataset


def get_split_reads(sample_alignment):
    """Isolate split reads from a sample alignment.

    This uses a python script supplied with Lumpy that is run as a
    separate process.

    NOTE THAT THIS SCRIPT ONLY WORKS WITH BWA MEM.
    """
    bwa_split_dataset = get_dataset_with_type(
            sample_alignment, Dataset.TYPE.BWA_SPLIT)
    if bwa_split_dataset is not None:
        if (bwa_split_dataset.status == Dataset.STATUS.READY and
                os.path.exists(bwa_split_dataset.get_absolute_location())):
            return bwa_split_dataset
    else:
        bwa_split_dataset = Dataset.objects.create(
                label=Dataset.TYPE.BWA_SPLIT,
                type=Dataset.TYPE.BWA_SPLIT,
                status=Dataset.STATUS.NOT_STARTED)
        sample_alignment.dataset_set.add(bwa_split_dataset)

    # If here, we are going to run or re-run the Dataset.
    bwa_split_dataset.status = Dataset.STATUS.NOT_STARTED
    bwa_split_dataset.save(update_fields=['status'])

    bam_dataset = get_dataset_with_type(sample_alignment, Dataset.TYPE.BWA_ALIGN)
    bam_filename = bam_dataset.get_absolute_location()

    assert os.path.exists(bam_filename), "BAM file '%s' is missing." % (
            bam_filename)

    bam_split_filename = os.path.join(sample_alignment.get_model_data_dir(),
            'bwa_split_reads.bam')

    try:
        bwa_split_dataset.status = Dataset.STATUS.COMPUTING
        bwa_split_dataset.save(update_fields=['status'])
        extract_split_reads(bam_filename, bam_split_filename)

    except subprocess.CalledProcessError:
        # if there are no split reads, then fail.
        bwa_split_dataset.filesystem_location = ''
        bwa_split_dataset.status = Dataset.STATUS.FAILED
    finally:
        bwa_split_dataset.status = Dataset.STATUS.READY
        bwa_split_dataset.filesystem_location = clean_filesystem_location(
                bam_split_filename)


    bwa_split_dataset.save()

    return bwa_split_dataset

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
    'bwa_align.sorted.realigned.bam'
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
