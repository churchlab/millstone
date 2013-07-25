"""
The alignment pipeline.

We start with the .fastq files and the reference for a particular genome
and carry out the gauntlet of steps that we are interested in.
"""

import os
import re
import shutil
import subprocess
from subprocess import CalledProcessError

import jbrowse_util
import util as gd_script_util
from util import copy_and_add_dataset_source
from util import genbank_to_fasta


def import_reference_genome(
        project,
        label,
        source_fasta=None,
        source_gbk=None,
        reference_genome=None,
        move_source=False):
    """Creates a ReferenceGenome associated with the given Project.

    Args:
        project: The Project we're storing everyting relative to.
        label: The human-readable label for the ReferenceGenome.
        source_fasta: Location of the source_fasta file.
        source_gbk: Location of the source_gbk file.
        reference_genome: Use this object if already exists.
        move_source: If True, move the source files tothe destination.
            Otherwise default behavior is to make a copy.
    """
    from main.models import Dataset
    from main.models import ReferenceGenome

    assert source_fasta or source_gbk, "Must have either gbk or fasta."
    if source_fasta:
        assert os.path.exists(source_fasta)
    if source_gbk:
        assert os.path.exists(source_gbk)

    if not reference_genome:
        reference_genome = ReferenceGenome.objects.create(
                project=project,
                label=label)

    if not source_fasta:
        source_fasta = os.path.join(
                reference_genome.get_data_directory_path(),
                os.path.splitext(os.path.split(source_gbk)[1])[0] + '.fa')
        genbank_to_fasta(source_gbk, source_fasta)

    # Add the datasets
    fasta_dataset = copy_and_add_dataset_source(
            reference_genome,
            Dataset.TYPE.REFERENCE_GENOME,
            Dataset.TYPE.REFERENCE_GENOME,
            source_fasta,
            unique=True
    )

    if source_gbk:
        gbk_dataset = copy_and_add_dataset_source(
                reference_genome,
                Dataset.TYPE.REFERENCE_GENOME_GBK,
                Dataset.TYPE.REFERENCE_GENOME_GBK,
                source_gbk,
                unique=True
        )

    # Prepare indexes for alignment.
    # TODO: Figure out how to integrate indexers at the point of use
    # instead of here awkwardly.
    bwa_build_index(fasta_dataset.filesystem_location)
    bowtie2_build_index(fasta_dataset.filesystem_location)

    # TODO: snpEff stuff
    #     * Prepare .gbk > .gff for snpEff run
    #     * Make snpEffectPredictor.bin

    return reference_genome


def import_genome(
        project,
        reference_genome,
        label,
        source_fastq1,
        source_fastq2=None):
    """Add a Genome to the ReferenceGenome.
    """
    from main.models import Dataset
    from main.models import Genome

    genome = Genome.objects.create(
            project=project,
            reference_genome=reference_genome,
            label=label
    )

    # Add the datasets.
    copy_and_add_dataset_source(
            genome,
            Dataset.TYPE.FASTQ1,
            Dataset.TYPE.FASTQ1,
            source_fastq1,
            unique=True
    )

    if source_fastq2:
        copy_and_add_dataset_source(
                genome,
                Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2,
                source_fastq2,
                unique=True
        )

    # Add to the ReferenceGenome.
    reference_genome.genome_set.add(genome)
    reference_genome.save()
    return genome


def import_genome_deprecated(label, project, ref_genome_fasta_location, **kwargs):
    """Deprecated.

    Creates a new Genome object for the project.

    NOTE: In progress.
    """
    from main.models import Dataset
    from main.models import DatasetType
    from main.models import Genome
    from main.models import ReferenceGenome

    genome = Genome.objects.create(
        project=project,
        label=label
    )

    # Add barcode if included in kwargs.
    BARCODE_KEY = 'barcode'
    if BARCODE_KEY in kwargs:
        genome.barcode = kwargs[BARCODE_KEY]
        genome.save()

    # Create a db record for the reference Genome.
    ref_genome_label = label + '_ref'
    ref_genome_dataset = Dataset.objects.create(
            label=ref_genome_label,
            dataset_type=DatasetType.objects.get_or_create(name='reference')[0],
            filesystem_location=ref_genome_fasta_location)
    reference_genome = ReferenceGenome.objects.create(
            label=ref_genome_label,
            source_dataset=ref_genome_dataset)
    genome.reference_genome = reference_genome
    genome.save()

    return genome


def _get_bowtie2_index_prefix(ref_genome_location):
    return os.path.splitext(ref_genome_location)[0]


def bowtie2_build_index(ref_genome_location):
    """Build the index for bowtie2."""
    index_prefix = _get_bowtie2_index_prefix(ref_genome_location)

    # Create the bowtie2 index.
    subprocess.check_call([
        'bowtie2-build',
        ref_genome_location,
        index_prefix
    ])


def align_with_bowtie2(
        reference_genome,
        genome,
        input_reads_1_fq,
        input_reads_2_fq=None):
    """Perform alignment using bowtie2. Adds the resulting bam_file as a
    Dataset to the genome.

    Args:
        ref_genome_fasta: Full path to the ref_genome. The resulting index name
                will be the prefix (.fa) stripped.
        genome: The Genome object which we can store Datasets for.
        input_reads_1_fq: Full path to the fastq reads file.
        input_reads_2_fq: Full path to the fastq reads file. If omitted,
            the alignment will be done with input_reads_1_fq treated as
            unpaired reads.
    """
    from main.models import Dataset

    # First, delete any existing BOWTIE2 Datasets.
    delete_existing_datasets_of_type(genome, Dataset.TYPE.BOWTIE2_ALIGN)

    # Create the Dataset with status COMPUTING.
    bowtie2_dataset = Dataset.objects.create(
            label=Dataset.TYPE.BOWTIE2_ALIGN,
            type=Dataset.TYPE.BOWTIE2_ALIGN,
            status=Dataset.STATUS.COMPUTING)
    genome.dataset_set.add(bowtie2_dataset)
    genome.save()

    ref_genome_fasta = reference_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME).filesystem_location

    output_bam = os.path.join(
            genome.get_data_directory_path(), 'bowtie2_align.bam')

    error_output = open(os.path.join(genome.get_data_directory_path(),
            'bowtie2_align.error'), 'w')

    try:
        index_prefix = _get_bowtie2_index_prefix(ref_genome_fasta)

        # Determine whether we are doing paired-end or not.
        is_paired_end = not input_reads_2_fq is None

        # If the fastq files are zipped, we need to unzip them using process
        # substitution before bowtie2 will read them. This looks like:
        # BEFORE: filename.fastq.gz
        # AFTER:  <(zcat filename.fastq.gz)
        input_reads_1_fq = '<( zcat ' + input_reads_1_fq + ' )'
        if is_paired_end:
            input_reads_2_fq = '<( zcat ' + input_reads_2_fq + ' )'

        # Now run alignment.
        with open(output_bam, 'w') as fh:
            THREADS = '1'
            args = [
                'bowtie2',
                '--local',
                '-q',
                '-p', THREADS,
                '--phred33',
                '-k', '3', # Search for multiple distinct alignments
                index_prefix,
            ]
            if is_paired_end:
                args.extend([
                    '-1', input_reads_1_fq,
                    '-2', input_reads_2_fq
                ])
            else:
                args.extend([
                    '-U', input_reads_1_fq,
                ])

            # in order to enable process substitution for unzipping 
            # gzipped fastq files in memory, we need to turn the arg list
            # into a string, and enable shell=True, setting the shell executable
            # to /bin/bash, which supports process substutiton (sh does not). 
            args = ' '.join(args)

            #we will also add a bash pipe to samtools to convert sam -> bam on the fly
            args += ' | samtools view -bS -'

            subprocess.check_call(args, stdout=fh, stderr=error_output, 
                shell=True, executable='/bin/bash')

        # Do several layers of processing on top of the initial alignment.
        result_bam_file = process_sam_bam_file(genome, output_bam, error_output)

    except subprocess.CalledProcessError as e:
        error_output.write(str(e))
        bowtie2_dataset.status = Dataset.STATUS.FAILED
        bowtie2_dataset.save()
        return

    else:
        bowtie2_dataset.filesystem_location = result_bam_file
        bowtie2_dataset.status = Dataset.STATUS.READY
        bowtie2_dataset.save()

        # Create a JBrowse track.
        jbrowse_util.add_bam_file_track(
                reference_genome, genome, Dataset.TYPE.BOWTIE2_ALIGN)

    finally:
        # No more errors here.
        error_output.close()


def bwa_build_index(ref_genome_fasta):
    subprocess.check_call([
        'bwa',
        'index',
        '-a',
        'is',
        ref_genome_fasta
    ])


def align_with_bwa(
        reference_genome,
        genome,
        input_reads_1_fq,
        input_reads_2_fq=None):
    """Perform alignment using bwa. Adds the resulting bam_file as a
    Dataset to the genome.

    Args:
        ref_genome_fasta: Full path to the ref_genome. The resulting index name
                will be the prefix (.fa) stripped.
        genome: The Genome object which we can store Datasets for.
        input_reads_1_fq: Full path to the fastq reads file.
        input_reads_2_fq: Full path to the fastq reads file.
        output_filename: Full path to the output .sam file.
    """

    from main.models import Dataset

    # First, delete any existing BOWTIE2 Datasets.
    delete_existing_datasets_of_type(genome, Dataset.TYPE.BWA_ALIGN)

    # Create the Dataset and set status as computing.
    bwa_dataset = Dataset.objects.create(
            label=Dataset.TYPE.BWA_ALIGN,
            type=Dataset.TYPE.BWA_ALIGN,
            status=Dataset.STATUS.COMPUTING)
    genome.dataset_set.add(bwa_dataset)
    genome.save()

    ref_genome_fasta = reference_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME).filesystem_location

    output_sam = os.path.join(
            genome.get_data_directory_path(), 'bwa_align.sam')

    error_output = open(os.path.join(genome.get_data_directory_path(),
            'bwa_align.error'), 'w')

    try:
        # Determine whether we are doing paired-end or not.
        is_paired_end = not input_reads_2_fq is None

        # There are two steps to alignment with bwa.
        #     1. Call 'bwa aln' to generate SA coordinate indeces. Note that
        #        we run a separate bwa aln for each end, and do so in parallel
        #        when paired end reads.
        #     2. Call 'bwa sampe' (paired-end) or 'bwa samse' (single-end) to
        #        generate SAM output.

        # 1. Generate SA coordinates.
        output_index_1 = os.path.splitext(input_reads_1_fq)[0] + '.sai'
        with open(output_index_1, 'w') as input_reads_1_stdout_fh:
            align_input_1_args = [
                'bwa',
                'aln',
                '-t', '6', # threads
                ref_genome_fasta,
                input_reads_1_fq,
            ]

            # NOTE: Non-blocking. See call to wait_and_check_pipe() below.
            align_input_1_proc = subprocess.Popen(align_input_1_args,
                    stdout=input_reads_1_stdout_fh, stderr=error_output)

        if is_paired_end:
            output_index_2 = os.path.splitext(input_reads_2_fq)[0] + '.sai'
            with open(output_index_2, 'w') as input_reads_2_stdout_fh:
                align_input_2_args = [
                    'bwa',
                    'aln',
                    '-t', '6', # threads
                    ref_genome_fasta,
                    input_reads_2_fq,
                ]
                align_input_2_proc = subprocess.Popen(align_input_2_args,
                        stdout=input_reads_2_stdout_fh, stderr=error_output)
                wait_and_check_pipe(align_input_2_proc, align_input_2_args)

        # Block here until both calls to 'bwa aln' complete.
        wait_and_check_pipe(align_input_1_proc, align_input_1_args)

        # 2. Generate SAM output.
        if is_paired_end:
            with open(output_sam, 'w') as fh:
                subprocess.check_call([
                    'bwa',
                    'sampe',
                    ref_genome_fasta,
                    output_index_1,
                    output_index_2,
                    input_reads_1_fq,
                    input_reads_2_fq,
                ], stdout=fh, stderr=error_output)
        else:
            with open(output_sam, 'w') as fh:
                subprocess.check_call([
                    'bwa',
                    'samse',
                    ref_genome_fasta,
                    output_index_1,
                    input_reads_1_fq,
                ], stdout=fh, stderr=error_output)

        # Do several layers of processing on top of the initial alignment.
        result_bam_file = process_sam_bam_file(genome, output_sam, error_output)

    except subprocess.CalledProcessError as e:
        error_output.write(str(e))
        bwa_dataset.status = Dataset.STATUS.FAILED
        bwa_dataset.save()
        return

    else:
        bwa_dataset.filesystem_location = result_bam_file
        bwa_dataset.status = Dataset.STATUS.READY
        bwa_dataset.save()

        # Create a JBrowse track.
        jbrowse_util.add_bam_file_track(
                reference_genome, genome, Dataset.TYPE.BWA_ALIGN)

    finally:
        error_output.close()

DEFAULT_PROCESSING_MASK = {
    'make_bam': True,
    'sort': True,
    'add_groups': True,
    'indel_realigner': True,
    'compute_insert_metrics': True,
    'index': True,
}

def process_sam_bam_file(genome, sam_bam_file_location, error_output,
        opt_processing_mask=DEFAULT_PROCESSING_MASK):
    """Converts to bam, sorts, and creates index.

    Args:
        genome: The Genome object.
        sam_bam_file_location: The full path to the .sam/.bam alignment output.
            If sam input, will be converted from sam to bam.
        error_output: File handle that can be passed as stderr to subprocess
            calls.
        opt_processing_mask: An optional mask that can be modified to
            specify which processes to run. This is useful when doing a
            re-run of a partially-complete run and previous steps completed.
            NOTE: Not for amateurs.

    Returns the path of the final .bam file.
    """
    # Make sure the mask is valid.
    assert (set(opt_processing_mask.keys()) ==
            set(DEFAULT_PROCESSING_MASK.keys()))

    #1. Convert to bam file if necessary
    if os.path.splitext(sam_bam_file_location)[1] == '.sam':
        sam_file_location = sam_bam_file_location
        bam_file_location = os.path.splitext(sam_file_location)[0] + '.bam'

        if opt_processing_mask['make_bam']:
            fh = open(bam_file_location, 'w')
            subprocess.check_call([
                'samtools',
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
    if opt_processing_mask['sort']:
        # 2a. Perform the actual sorting.
        subprocess.check_call([
            'samtools',
            'sort',
            bam_file_location,
            sorted_output_name
        ], stderr=error_output)

        # 2b. Index the sorted result.
        subprocess.check_call([
            'samtools',
            'index',
            sorted_bam_file_location,
        ], stderr=error_output)

    # 3. Add groups
    grouped_output_name= (
            os.path.splitext(sorted_bam_file_location)[0] +
            '.grouped')
    grouped_bam_file_location = grouped_output_name + '.bam'
    if opt_processing_mask['add_groups']:
        add_groups(genome, sorted_bam_file_location, grouped_bam_file_location,
                error_output)

    # 3. Perform realignment accounting for indels.
    realigned_bam_file_location = (
            os.path.splitext(grouped_bam_file_location)[0] +
            '.realigned.bam'
    )
    if opt_processing_mask['indel_realigner']:
        # Make sure the previous result is indexed.
        subprocess.check_call([
                'samtools',
                'index',
                grouped_bam_file_location,
        ], stderr=error_output)

        realign_given_indels(
                genome,
                grouped_bam_file_location,
                realigned_bam_file_location,
                error_output
        )

    # 4. Compute insert size metrics
    if opt_processing_mask['compute_insert_metrics']:
        compute_insert_metrics(realigned_bam_file_location, error_output)

    # 5. Create index.
    if opt_processing_mask['index']:
        subprocess.check_call([
            'samtools',
            'index',
            realigned_bam_file_location,
        ], stderr=error_output)

    return realigned_bam_file_location


def add_groups(genome, input_bam_file, bam_output_file, error_output):
    from main.models import ensure_exists_0775_dir
    # Create a temp directory in the genome data dir else we run out of space
    # with big alignment jobs.
    picard_tmp_dir = os.path.join(genome.get_data_directory_path(), 'picard_tmp')
    ensure_exists_0775_dir(picard_tmp_dir)

    prefix = genome.uid
    subprocess.check_call([
        'java', '-jar', '/opt/picard/AddOrReplaceReadGroups.jar',
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
        genome,
        input_bam_file,
        output_bam_file,
        error_output):
    """Perform realignment accounting for indels. This should give a cleaner
    alignment.
    """
    temp_intervals_file = os.path.join(
            genome.get_data_directory_path(), 'realign.intervals')
    ref_genome_fasta_location = (
            genome.reference_genome.get_fasta_source().filesystem_location)

    # Prepare realignment intervals file.
    subprocess.check_call([
        'java', '-jar', '/opt/gatk/GenomeAnalysisTK.jar',
        '-T', 'RealignerTargetCreator',
        '-I', input_bam_file,
        '-R', ref_genome_fasta_location,
        '-o', temp_intervals_file,
    ], stderr=error_output)

    # Perform realignment.
    subprocess.check_call([
        'java', '-jar', '/opt/gatk/GenomeAnalysisTK.jar',
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
        'java', '-jar', '/opt/picard/CollectInsertSizeMetrics.jar',
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


def delete_existing_datasets_of_type(genome, dataset_type):
    """Helper function that deletes any existing Dataset with given type."""
    from main.models import Dataset

    type_error_msg = "alignment_type must be valid Dataset type"
    assert dataset_type in dir(Dataset.TYPE), type_error_msg

    deleted = True
    while deleted:
        deleted = False
        dataset_set = genome.dataset_set.filter(type=dataset_type)
        if len(dataset_set) > 0:
            dataset_set[0].delete()
            deleted = True


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


if __name__ == '__main__':
    gd_script_util.setup_django_env()

    from main.models import Dataset
    from main.models import Genome
    from main.models import ReferenceGenome

    reference_genome = ReferenceGenome.objects.get(uid='b4fae5bc')
    # genome = Genome.objects.get(uid='d676dbe2') # BW25
    genome = Genome.objects.get(uid='90980be7') # BW38

    alignment_type = Dataset.TYPE.BWA_ALIGN
    delete_existing_datasets_of_type(genome, alignment_type)

    # Create the Dataset and set status as computing.
    dataset = Dataset.objects.create(
            label=alignment_type,
            type=alignment_type,
            status=Dataset.STATUS.COMPUTING)
    genome.dataset_set.add(dataset)
    genome.save()

    try:
        sam_file_location = os.path.join(
            genome.get_data_directory_path(), 'bwa_align.sam')

        error_output = open(os.path.join(genome.get_data_directory_path(),
                'bwa_align.error'), 'a')

        mask = {
            'make_bam': False,
            'sort': True,
            'add_groups': True,
            'index': True,
        }

        result_bam_file = process_sam_bam_file(genome, sam_file_location, error_output, mask)

    except subprocess.CalledProcessError as e:
        error_output.write(str(e))
        dataset.status = Dataset.STATUS.FAILED
        dataset.save()

    else:
        dataset.filesystem_location = result_bam_file
        dataset.status = Dataset.STATUS.READY
        dataset.save()

        # Create a JBrowse track.
        jbrowse_util.add_bam_file_track(
                reference_genome, genome, alignment_type)

    finally:
        error_output.close()
