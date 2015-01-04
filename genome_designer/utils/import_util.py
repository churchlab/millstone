"""
Methods related to importing data.
"""

import copy
import csv
import os
import shutil
import re
import subprocess
from tempfile import mkdtemp
from tempfile import mkstemp
from tempfile import NamedTemporaryFile

from BCBio import GFF
from Bio import Entrez
from Bio import SeqIO
from celery import task
from django.conf import settings
from django.db import transaction

from main.celery_util import assert_celery_running
from main.exceptions import ValidationException
from main.models import Chromosome
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ReferenceGenome
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.model_utils import clean_filesystem_location
from main.model_utils import get_dataset_with_type
from main.s3 import project_files_needed
from pipeline.read_alignment_util import ensure_bwa_index
from pipeline.variant_effects import build_snpeff
from utils import generate_safe_filename_prefix_from_label
from utils import uppercase_underscore
from utils.jbrowse_util import prepare_jbrowse_ref_sequence
from utils.jbrowse_util import add_genbank_file_track
from variants.vcf_parser import get_or_create_variant


IMPORT_FORMAT_TO_DATASET_TYPE = {
    'fasta': Dataset.TYPE.REFERENCE_GENOME_FASTA,
    'genbank': Dataset.TYPE.REFERENCE_GENOME_GENBANK,
    'gff': Dataset.TYPE.REFERENCE_GENOME_GFF,
    'vcf_user': Dataset.TYPE.VCF_USERINPUT
}

SAMPLE_SERVER_COPY_KEY__SAMPLE_NAME = 'Sample_Name'
SAMPLE_SERVER_COPY_KEY__READ_1 = 'Read_1_Path'
SAMPLE_SERVER_COPY_KEY__READ_2 = 'Read_2_Path'

REQUIRED_SAMPLE_SERVER_COPY_HEADER = [
    SAMPLE_SERVER_COPY_KEY__SAMPLE_NAME,
    SAMPLE_SERVER_COPY_KEY__READ_1,
]

# Cols that we know about, to distinguish them from user-defined cols.
PRE_DEFINED_SAMPLE_SERVER_COPY_HEADER_PARTS = (
        REQUIRED_SAMPLE_SERVER_COPY_HEADER +
        [SAMPLE_SERVER_COPY_KEY__READ_2])

SAMPLE_BROWSER_UPLOAD_KEY__SAMPLE_NAME = 'Sample_Name'
SAMPLE_BROWSER_UPLOAD_KEY__READ_1 = 'Read_1_Filename'
SAMPLE_BROWSER_UPLOAD_KEY__READ_2 = 'Read_2_Filename'

REQUIRED_SAMPLE_UPLOAD_THROUGH_BROWSER_HEADER = [
    SAMPLE_BROWSER_UPLOAD_KEY__SAMPLE_NAME,
    SAMPLE_BROWSER_UPLOAD_KEY__READ_1,
]

# Cols that we know about, to distinguish them from user-defined cols.
PRE_DEFINED_SAMPLE_UPLOAD_THROUGH_BROWSER_PARTS = (
        REQUIRED_SAMPLE_UPLOAD_THROUGH_BROWSER_HEADER +
        [SAMPLE_BROWSER_UPLOAD_KEY__READ_2])

REQUIRED_VCF_HEADER_PART = ['CHROM', 'POS', 'ID', 'REF', 'ALT']


if settings.S3_ENABLED:
    from main.s3 import s3_temp_get, s3_get

    def import_reference_genome_from_s3(project, label, s3file, import_format):
        with s3_temp_get(s3file) as f:
            return import_reference_genome_from_local_file(
                    project, label, f, import_format)

    @project_files_needed
    def import_samples_from_s3(project, targets_file_rows, s3files):
        tmp_dir = mkdtemp()
        local_s3files_map = {}
        for s3file in s3files:
            filepath = os.path.join(tmp_dir, s3file.name)
            s3_get(s3file.key, filepath)
            local_s3files_map[s3file.name] = filepath

        for row in targets_file_rows:
            sample_label = row['Sample_Name']
            experiment_sample = ExperimentSample.objects.create(
                    project=project, label=sample_label)
            copy_and_add_dataset_source(experiment_sample, Dataset.TYPE.FASTQ1,
                    Dataset.TYPE.FASTQ1, local_s3files_map[row['Read_1_Path']])
            if 'Read_2_Path' in row and row['Read_2_Path']:
                copy_and_add_dataset_source(experiment_sample, Dataset.TYPE.FASTQ2,
                        Dataset.TYPE.FASTQ2, local_s3files_map[row['Read_2_Path']])
        shutil.rmtree(tmp_dir)


class DataImportError(Exception):
    """Exception thrown when there are errors in imported data.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return 'DataImportError: ' + str(self.msg)


@project_files_needed
def import_reference_genome_from_local_file(project, label, file_location,
        import_format, move=False):
    """Creates a ReferenceGenome associated with the given Project.

    Args:
        project: The Project we're storing everyting relative to.
        label: The human-readable label for the ReferenceGenome.
        file_location: Location of the genome on the server.
        import_format: Must be 'fasta' or 'genbank'.
        move: move instead of copy the original file_location - for instance,
        if we saved it to a temporary file. Moving is of course faster than
        copying.

    Returns:
        ReferenceGenome.
    """
    # Validate the input.
    assert import_format in ['fasta', 'genbank']

    # Validate the file.
    assert os.path.exists(file_location), "File %s doesn't exist." % (
            file_location)

    # Validate the input by parsing it with BioPython, while also
    # counting the number of chromosomes.
    num_bases = 0
    for genome_record in SeqIO.parse(file_location, import_format):
        num_bases += len(genome_record)

    # Make sure sequence exists.
    if not num_bases > 0:
        raise DataImportError("No sequence in file.")

    # Create the ReferenceGenome object.
    reference_genome = ReferenceGenome.objects.create(
            project=project,
            label=label)

    # Copy the source file to the ReferenceGenome data location.
    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE[import_format]
    copy_and_add_dataset_source(reference_genome, dataset_type,
            dataset_type, file_location)

    return reference_genome


def add_chromosomes(reference_genome, dataset):
    """ Makes a Chromosome for each unique SeqRecord.name in the dataset
    """
    
    chromosomes = [chrom.label for chrom \
            in Chromosome.objects.filter(reference_genome=reference_genome)]

    def _make_chromosome(seq_rec_iter):
        for seq_record in seq_rec_iter:
            if seq_record.name and seq_record.name not in chromosomes:
                Chromosome.objects.create(reference_genome=reference_genome,
                        label=seq_record.name, num_bases=len(seq_record))

    # Add chromosome ids
    dataset_path = dataset.get_absolute_location()

    if dataset.TYPE.REFERENCE_GENOME_FASTA:
        _make_chromosome(SeqIO.parse(dataset_path, "fasta"))
    elif dataset.TYPE.REFERENCE_GENOME_GENBANK:
        _make_chromosome(SeqIO.parse(dataset_path, "genbank"))     
    else:
        raise AssertionError("Unexpected Dataset type")


def generate_fasta_from_genbank(ref_genome):
    """If this reference genome has a genbank but not a FASTA, generate
    a FASTA from the genbank. """

    # If a FASTA already exists, then just return.
    if ref_genome.dataset_set.filter(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).exists():
        return

    # Check that a genbank exists.
    assert ref_genome.dataset_set.filter(
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).exists()

    # Get genbank path and filename components (for creating FASTA file name).
    genbank_path = get_dataset_with_type(
            ref_genome,
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()

    genbank_dir, genbank_filename = os.path.split(genbank_path)
    genbank_noext = os.path.splitext(genbank_filename)[0]

    # Put the fasta file in the same dir, just change the extension to .fa.
    fasta_filename = os.path.join(genbank_dir, (genbank_noext + '.fa'))

    # Get the individual records, each corresponding to a chromosome.
    genome_records = list(SeqIO.parse(genbank_path, 'genbank'))

    # SnpEFF takes the name attr, but the BioPython uses the id attr to make its
    # fasta file, so overwrite the id with the name when converting to fasta.
    for genome_record in genome_records:
        genome_record.id = genome_record.name

    SeqIO.write(genome_records, fasta_filename, 'fasta')

    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE['fasta']
    copy_and_add_dataset_source(ref_genome, dataset_type,
            dataset_type, fasta_filename)

def ensure_fasta_index(ref_genome_fasta):
    """
    Check if a fasta index is present w/ extension .fai. If not,
    use samtools to generate one.
    """
    if not os.path.exists(ref_genome_fasta + '.fai'):
        subprocess.check_call([
            settings.SAMTOOLS_BINARY,
            'faidx',
            ref_genome_fasta])


def generate_gff_from_genbank(ref_genome):
    """If this reference genome has a genbank but not a GFF, generate
    a GFF from the genbank. """

    # If a GFF already exists, then just return.
    if ref_genome.dataset_set.filter(
            type=Dataset.TYPE.REFERENCE_GENOME_GFF).exists():
        return

    # Check that a genbank exists.
    assert ref_genome.dataset_set.filter(
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).exists()

    # Get genbank path and filename components (for creating GFF file name).
    genbank_path = get_dataset_with_type(
            ref_genome,
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()

    genbank_dir, genbank_filename = os.path.split(genbank_path)
    genbank_noext = os.path.splitext(genbank_filename)[0]

    # Put the GFF file in the same dir, just change the extension to .gff.
    gff_filename = os.path.join(genbank_dir, (genbank_noext + '.gff'))

    # Get the individual records, each corresponding to a chromosome.
    genome_records = list(SeqIO.parse(genbank_path, 'genbank'))

    # SnpEFF takes the name attr, but the BioPython uses the id attr to make its
    # GFF file, so overwrite the id with the name when converting to GFF.

    for genome_record in genome_records:
        genome_record.id = genome_record.name

    GFF.write(genome_records, open(gff_filename, 'w'))

    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE['gff']
    copy_and_add_dataset_source(ref_genome, dataset_type,
            dataset_type, gff_filename)


def import_reference_genome_from_ncbi(project, label, record_id, import_format):
    """Imports a reference genome by accession from NCBI using efetch.
    """
    # Validate the input.
    assert import_format in ['fasta', 'genbank'], (
            'Import Format must be \'fasta\' or \'genbank\'')

    # Format keys for Efetch.
    # More info at:  http://www.ncbi.nlm.nih.gov/
    #       books/NBK25499/table/chapter4.chapter4_table1/?report=objectonly
    CONVERT_FORMAT = {
        'fasta': 'fa',
        'genbank': 'gbwithparts'
    }

    # What suffix to use for each input format
    # TODO: Should this be a property of the Dataset TYPE?
    FORMAT_SUFFIX = {
        'fasta': '.fa',
        'genbank': '.gb'
    }

    Entrez.email = settings.EMAIL
    handle = Entrez.efetch(
            db="nuccore",
            id=record_id,
            rettype=CONVERT_FORMAT[import_format],
            retmode="text")

    # Store results in temporary file.
    filename_prefix = generate_safe_filename_prefix_from_label(label) + '_'
    temp = NamedTemporaryFile(delete=False, prefix=filename_prefix,
            suffix=FORMAT_SUFFIX[import_format])
    temp.write(handle.read())
    handle.close()
    temp.close()

    # Create ref genome from this temporary file.
    reference_genome = import_reference_genome_from_local_file(
        project, label, temp.name, import_format, move=True)

    # Clean up.
    if os.path.isfile(temp.name):
        os.remove(temp.name)

    return reference_genome


def sanitize_record_id(record_id_string):
    """We want to grab only the first word-only part of each seqrecord in a
    FASTA/Genbank file, and use that as a consistent and readable id between
    genbank and FASTA.
    """
    return re.match( r'^\w{1,20}', record_id_string).group()


def _assert_sample_targets_file_size(targets_file):
    if hasattr(targets_file, "size"):
        assert targets_file.size < 1000000, (
                "Targets file is too large: %d" % targets_file.size)


@project_files_needed
def import_samples_from_targets_file(project, targets_file, options={}):
    """Uses the uploaded targets file to add a set of samples to the project.
    We need to check each line of the targets file for consistency before we
    do anything, however. Checking is moved to parse_targets_file() which parses
    targets_file and returns valid rows. parse_targets_file() will also be
    called from parse_targets_file_s3 in xhr_handlers in case of S3 uploading.

    It writes a copy of the uploaded targets file to a temporary file

    Args:
        project: The project we're storing everything relative to>
        targets_file: The UploadedFile django object that holds the targets
            in .tsv format.
        options: Dictionary of options. Currently a hack to allow different
            parts of the pipeline to not run during tests (e.g. FastQC).
    """
    assert_celery_running()

    parsed_rows = parse_experiment_sample_targets_file(
            project,
            targets_file,
            REQUIRED_SAMPLE_SERVER_COPY_HEADER,
            SAMPLE_SERVER_COPY_KEY__SAMPLE_NAME,
            SAMPLE_SERVER_COPY_KEY__READ_1,
            SAMPLE_SERVER_COPY_KEY__READ_2)

    # We perform the additional step of testing file locations, and for the
    # test data, switching out $GD_ROOT template variable.
    valid_rows = []
    for row in parsed_rows:
        updated_row = copy.copy(row)
        for field, value in row.iteritems():
            if field in (SAMPLE_SERVER_COPY_KEY__READ_1,
                    SAMPLE_SERVER_COPY_KEY__READ_2):
                updated_value = value.replace('$GD_ROOT', settings.PWD)
                with open(updated_value, 'rb') as test_file:
                    try:
                        test_file.read(8)
                    except:
                        raise AssertionError(
                                "Cannot read file at %s" % updated_value)
                updated_row[field] = updated_value
        valid_rows.append(updated_row)

    return create_samples_from_row_data(project, valid_rows, move=False,
            options=options)


def create_samples_from_row_data(
        project, data_source_list, move=False, options={}):
    """Creates ExperimentSample objects along with their respective Datasets.

    The data is copied to the entity location. We block until we've created the
    models, and then go async for actual copying.

    Args:
        project: Project these Samples should be added to.
        data_source_list: List of objects with keys:
            * Sample_Name
            * Read_1_Path
            * Read_2_Path (optional)
            * other metadata keys (optional)
            * ...
        move: Whether to move the source data. Else copy.
        options: Dictionary of options. Currently a hack to allow different
            parts of the pipeline to not run during tests (e.g. FastQC).

    Returns:
        List of ExperimentSamples.
    """
    experiment_samples = []
    for row in data_source_list:
        # Create ExperimentSample object and then store the data relative to
        # it.
        sample_label = row['Sample_Name']
        experiment_sample = ExperimentSample.objects.create(
                project=project, label=sample_label)

        # Create the Datasets before starting copying so we can show status in
        # the ui. This is a new pattern where we are moving copying to happen
        # asynchronously in Celery.
        _create_fastq_dataset(
                experiment_sample, row['Read_1_Path'], Dataset.TYPE.FASTQ1,
                Dataset.STATUS.QUEUED_TO_COPY)
        maybe_read2_path = row.get('Read_2_Path', '')
        if maybe_read2_path:
            _create_fastq_dataset(
                    experiment_sample, maybe_read2_path, Dataset.TYPE.FASTQ2,
                    Dataset.STATUS.QUEUED_TO_COPY)

        # Add extra metadata columns.
        _update_experiment_sample_data_for_row(experiment_sample, row,
                PRE_DEFINED_SAMPLE_SERVER_COPY_HEADER_PARTS)

        # Start the async job of copying.
        copy_experiment_sample_data.delay(
                project, experiment_sample, row, move=move, options=options)

        experiment_samples.append(experiment_sample)

    _update_experiment_sample_parentage(experiment_samples)

    return experiment_samples


def _create_fastq_dataset(experiment_sample, fastq_source, dataset_type,
        dataset_status):
    """Helper function for creating a Dataset that will point to a file.

    Since clients of this function are responsible for actuallying copying the
    data, this function sets the is_present bit to False on the Dataset.
    """
    fastq_dest = _get_copy_target_path(experiment_sample, fastq_source)
    reads_dataset = add_dataset_to_entity(experiment_sample,
            dataset_type, dataset_type, fastq_dest)
    reads_dataset.status = dataset_status
    reads_dataset.save()
    return reads_dataset


def _copy_dataset_data(experiment_sample, fastq_source, dataset_type,
            move=False, set_status=Dataset.STATUS.VERIFYING):
    """Helper to copy data and set status.
    """
    dataset = experiment_sample.dataset_set.get(type=dataset_type)
    dataset.status = Dataset.STATUS.COPYING
    dataset.save()
    copy_dataset_to_entity_data_dir(experiment_sample, fastq_source, move=move)
    dataset.status = set_status
    dataset.save()
    return dataset


@task
@project_files_needed
def copy_experiment_sample_data(
        project, experiment_sample, data, move=False,
        options={'skip_fastqc': False}):
    """Celery task that wraps the process of copying the data for an
    ExperimentSample.
    """

    # Copy read1.
    read1_dataset = _copy_dataset_data(experiment_sample, data['Read_1_Path'],
            Dataset.TYPE.FASTQ1, move=move)
    # Copy read2.
    maybe_read2_path = data.get('Read_2_Path', '')
    if maybe_read2_path:
        read2_dataset = _copy_dataset_data(experiment_sample, maybe_read2_path,
                Dataset.TYPE.FASTQ2, move=move)
    else:
        read2_dataset = None

    # Verification.
    if read2_dataset is not None:
        # Paired reads.
        if (read1_dataset.filesystem_location ==
                read2_dataset.filesystem_location):
            # Make sure the files are not the same.
            read1_dataset.status = Dataset.STATUS.FAILED
            read2_dataset.status = Dataset.STATUS.FAILED

            # TODO: Provide way for user to get an error message, similar to
            # how make an error link for alignments.
        else:
            read1_dataset.status = Dataset.STATUS.READY
            read2_dataset.status = Dataset.STATUS.READY
    else:
        # Unpaired.
        read1_dataset.status = Dataset.STATUS.READY
        read2_dataset.status = Dataset.STATUS.READY

    # Quality Control via FASTQC and save.
    read1_dataset.status = Dataset.STATUS.QC
    if not options.get('skip_fastqc', False):
        run_fastqc_on_sample_fastq(experiment_sample, read1_dataset)
    read1_dataset.status = Dataset.STATUS.READY
    read1_dataset.save()

    read2_dataset.status = Dataset.STATUS.QC
    if not options.get('skip_fastqc', False):
        run_fastqc_on_sample_fastq(experiment_sample, read2_dataset, rev=True)
    read2_dataset.status = Dataset.STATUS.READY
    read2_dataset.save()


def run_fastqc_on_sample_fastq(experiment_sample, fastq_dataset, rev=False):
    """Runs FASTQC on a fastq dataset object.

    Args:
        experiment_sample: The ExperimentSample for this fastq.
        fastq_dataset: Dataset that points to uploaded fastq file.

    Returns:
        New Dataset pointing to html file of FastQC results.
    """
    fastq_filename = fastq_dataset.get_absolute_location()

    # There's no option to pass the output filename to FastQC so we just
    # create the name that matches what FastQC outputs.
    if fastq_dataset.is_compressed():
        unzipped_fastq_filename = os.path.splitext(fastq_filename)[0]
    else:
        unzipped_fastq_filename = fastq_filename
    fastqc_filename = unzipped_fastq_filename + '_fastqc.html'

    if rev:
        dataset_type = Dataset.TYPE.FASTQC2_HTML
    else:
        dataset_type = Dataset.TYPE.FASTQC1_HTML

    # create the tmp dir if it doesn't exist
    if not os.path.exists(settings.TEMP_FILE_ROOT):
        os.mkdir(settings.TEMP_FILE_ROOT)

    command = [
            settings.FASTQC_BINARY,
            fastq_filename,
            '-o', experiment_sample.get_model_data_dir(),
            '-d', settings.TEMP_FILE_ROOT]

    fastqc_output = subprocess.check_output(
            command, stderr=subprocess.STDOUT)

    # Check that fastqc file has been made
    # TODO: We need proper error checking and logging probably, so that this
    # non-essential step doesn't destroy the whole import process.
    if not os.path.exists(fastqc_filename):
        print 'FastQC Failed for {}:\n{}'.format(
                    fastq_filename, fastqc_output)

    fastqc_dataset = add_dataset_to_entity(experiment_sample,
            dataset_type, dataset_type, fastqc_filename)
    fastqc_dataset.status = Dataset.STATUS.READY
    fastqc_dataset.save()

    return fastqc_dataset


def create_sample_models_for_eventual_upload(project, targets_file):
    """Parses the form to create sample placeholers that are awaiting
    data upload.

    Args:
        targets_file: The filled out form.

    Raises:
        ValidationException if validation fails.
    """
    try:
        valid_rows = parse_experiment_sample_targets_file(
                project,
                targets_file,
                REQUIRED_SAMPLE_UPLOAD_THROUGH_BROWSER_HEADER,
                SAMPLE_BROWSER_UPLOAD_KEY__SAMPLE_NAME,
                SAMPLE_BROWSER_UPLOAD_KEY__READ_1,
                SAMPLE_BROWSER_UPLOAD_KEY__READ_2)
    except AssertionError as e:
        raise ValidationException(e)

    for row in valid_rows:
        _create_sample_and_placeholder_dataset(project, row)


def _create_sample_and_placeholder_dataset(project, row):
    """Create Datasets but don't copy data.
    """
    # Parsing and validation.
    fastq1_filename = row['Read_1_Filename']
    maybe_fastq2_filename = row.get('Read_2_Filename', '')
    assert fastq1_filename != maybe_fastq2_filename

    # Now create the models.
    experiment_sample = ExperimentSample.objects.create(
            project=project, label=row['Sample_Name'])

    fastq1_filename = row['Read_1_Filename']
    _create_fastq_dataset(
            experiment_sample, fastq1_filename, Dataset.TYPE.FASTQ1,
            Dataset.STATUS.AWAITING_UPLOAD)

    if maybe_fastq2_filename:
        _create_fastq_dataset(
                experiment_sample, maybe_fastq2_filename,
                Dataset.TYPE.FASTQ2, Dataset.STATUS.AWAITING_UPLOAD)

    # Add extra metadata columns.
    _update_experiment_sample_data_for_row(experiment_sample, row,
            PRE_DEFINED_SAMPLE_UPLOAD_THROUGH_BROWSER_PARTS)


def _update_experiment_sample_data_for_row(experiment_sample, row, known_cols):
    """Updates the catch-all ExperimentSample.data field with user-defined
    fields.
    """
    for field, value in row.iteritems():
        if field not in known_cols:
            clean_field = uppercase_underscore(field)
            if not clean_field.startswith('SAMPLE_'):
                clean_field = 'SAMPLE_' + clean_field
            experiment_sample.data[clean_field] = str(value)

    experiment_sample.save(update_fields=['data'])


def _update_experiment_sample_parentage(experiment_samples):
    """
    Adds children/parent relations to ExperimentSamples according to a
    dictionary of parents parsed out of the targets file.
    """

    es_label_dict = dict([(es.label, es) for es in experiment_samples])

    # make a dict of child to parent labels.
    parent_dict = {}
    for es in experiment_samples:
        if 'SAMPLE_PARENTS' in es.data:
            parents = es.data['SAMPLE_PARENTS'].split('|')
            parent_dict[es.label] = parents

    # add the children of each parent to the model.
    for child, parents in parent_dict.items():
        assert child in es_label_dict.keys(), (
                'Child {} missing from ExperimentSamples'.format(child))
        for parent in parents:
            if not parent: continue #skip blank strings
            assert parent in es_label_dict.keys(), (
                   'Parent {} missing from ExperimentSamples'.format(parent))
            es_label_dict[parent].add_child(es_label_dict[child])

def parse_experiment_sample_targets_file(project,
        targets_filehandle_or_filename, required_header, sample_name_key,
        read_1_key, read_2_key):
    """Parses and validates the file.

    Returns:
        List of objects representing the rows.
    """
    _assert_sample_targets_file_size(targets_filehandle_or_filename)

    # The purpose of the next few lines of somewhat convoluted code is to
    # make sure we support weird template formats such as Excel on OsX might
    # output. In the end, we want to end up with the variable targets_file
    # being a File object that has been read in in universal mode,
    # open(..., 'rU'). This requirement is made slightly trickier by the fact
    # that the aptly named param targets_filehandle_or_filename is of ambiguous
    # type (because Python) and so the remaining code needs to work whether
    # it's a string filename, or a File object. One way we can solve all
    # these constraints is to write the contents of the file to a temporary
    # location, and then read it back in universal mode. I would welcome a more
    # elegant fix.
    if isinstance(targets_filehandle_or_filename, str):
        temp_file_location = targets_filehandle_or_filename
    else:
        # It's an open File object.
        if not os.path.exists(settings.TEMP_FILE_ROOT):
            os.mkdir(settings.TEMP_FILE_ROOT)
        _, temp_file_location = mkstemp(dir=settings.TEMP_FILE_ROOT)
        with open(temp_file_location, 'w') as temp_fh:
            temp_fh.write(targets_filehandle_or_filename.read())
    targets_file = open(temp_file_location, 'rU')
    os.remove(temp_file_location)

    # Now this works even if there are silly carriage return characters ^M.
    reader = csv.DictReader(targets_file, delimiter='\t')

    # Read the header / schema.
    targets_file_header = reader.fieldnames

    # Make sure all header cols are present.
    missing_header_cols = (set(required_header) - set(targets_file_header))
    assert 0 == len(missing_header_cols), (
            "Missing cols: %s" % ' '.join(missing_header_cols))

    # Query all relevant datasets to check for filename clashes.
    existing_sample_dataset_filename_set = set([
            os.path.split(ds.filesystem_location)[1]
            for ds in Dataset.objects.filter(
                    experimentsample__project=project)])

    # Set this to a boolean on the first iteration, and make sure all rows
    # are either paired or unpaired.
    is_paired_end = None

    # Initial aggregation and validation.
    valid_rows = []
    for raw_row_obj in reader:
        clean_row_obj = {}
        for key, value in raw_row_obj.iteritems():
            # Ignore rows of the form K/V pair {None: ''}
            if key is None:
                continue
            clean_row_obj[key] = value.strip()

        sample_name = clean_row_obj[sample_name_key]
        if not sample_name:
            # Null sample name, skip the row.
            continue

        assert len(targets_file_header) == len(clean_row_obj.keys()), (
                "Row %s has the wrong number of fields." % sample_name)

        # Determine whether paired-end data (first iteration only).
        if is_paired_end is None:
            is_paired_end = (read_2_key in clean_row_obj and
                    clean_row_obj[read_2_key])

        # Check filenames are present.
        assert clean_row_obj[read_1_key], (
                "No read 1 in row %s" % sample_name)
        if is_paired_end:
            assert clean_row_obj[read_2_key], (
                    "No read 2 in row %s" % sample_name)

        # Catch a common copy paste error where read1 matches read2.
        if is_paired_end:
            same = (clean_row_obj[read_1_key] == clean_row_obj[read_2_key])
            assert not same, "Read 1 filename is same as read 2 filename"

        # Make sure Dataset with that name doesn't exist.
        def _assert_not_filename_exists(filename_col):
            filename = os.path.basename(clean_row_obj[filename_col])
            assert not filename in existing_sample_dataset_filename_set, (
                "%s exists" % clean_row_obj[filename_col])
        _assert_not_filename_exists(read_1_key)
        if is_paired_end:
            _assert_not_filename_exists(read_2_key)

        valid_rows.append(clean_row_obj)

    # Make sure all the standard fields have unique values relative to each
    # other.
    def _assert_no_repeated_value(col):
        values = set([row[col] for row in valid_rows])
        assert len(values) == len(valid_rows), (
                "Non-unique %s detected." % col)
    _assert_no_repeated_value(sample_name_key)
    _assert_no_repeated_value(read_1_key)
    if is_paired_end:
        _assert_no_repeated_value(read_2_key)

    targets_file.close()

    return valid_rows


@transaction.commit_on_success
def import_variant_set_from_vcf(ref_genome, variant_set_name, variant_set_file):
    """Convert an uploaded VCF file into a new variant set object.

    Args:
        ref_genome: ReferenceGenome.
        variant_set_name: Name of the variant set (label).
        variant_set_file: Path to the variant set on disk.
    """
    # For now, variant set name must be unique even among diff ref genomes.
    variant_set_name_exists = bool(VariantSet.objects.filter(
            reference_genome=ref_genome,
            label=variant_set_name).count())
    assert not variant_set_name_exists, 'Variant set name must be unique.'

    # Create the VariantSet.
    variant_set = VariantSet.objects.create(
            reference_genome=ref_genome,
            label=variant_set_name)

    # First, save this vcf as a dataset, so we can point to it from the
    # new variant common_data_objs
    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE['vcf_user']
    dataset = copy_and_add_dataset_source(variant_set, dataset_type,
            dataset_type, variant_set_file)

    # Now read the variant set file.
    _read_variant_set_file_as_csv(variant_set_file, ref_genome, dataset,
            variant_set)

    # These actions invalidate the materialized view.
    ref_genome.invalidate_materialized_view()


def _read_variant_set_file_as_csv(variant_set_file, reference_genome,
        dataset, variant_set):
    """If reading the variant set file as a vcf fails (because we arent using
    all columns, as will usually be the case) then read it as a CSV and check
    manually for the required columns.

    Args:
        * variant_set_file: Path to vcf file.
        * reference_genome: ReferenceGenome object.
    """

    with open(variant_set_file) as fh:
         # Use this wrapper to skip the header lines
        # Double ##s are part of the header, but single #s are column
        # headings and must be stripped and kept.
        def remove_vcf_header(iterable):
            for line in iterable:
                if not line.startswith('##'):
                    if line.startswith('#'):
                        line = line.lstrip('#')
                    yield line
        vcf_noheader = remove_vcf_header(open(variant_set_file))

        reader = csv.DictReader(vcf_noheader, delimiter='\t')

        # Check that the required columns are present.
        assert (len(reader.fieldnames) >= len(REQUIRED_VCF_HEADER_PART)), (
            'Header for PseudoVCF %s is too short, should have [%s], has %s' % (
                    variant_set_file, ', '.join(REQUIRED_VCF_HEADER_PART),
                    ', '.join(reader.fieldnames)))

        for col, check in zip(reader.fieldnames[0:len(REQUIRED_VCF_HEADER_PART)],
                REQUIRED_VCF_HEADER_PART):
            assert col == check, (
                "Header column '%s' is missing or out of order; %s" % (check,
                    ', '.join(reader.fieldnames)))

        class PseudoVCF:
            """Pseudo wrapper class to satisfy interface of
            extract_raw_data_dict().
            """
            def __init__(self, **entries):
                self.__dict__.update(entries)
                self.__dict__['ALT'] = self.__dict__['ALT'].strip().split(',')
                self.__dict__['samples'] = []

        for record in reader:

            record = PseudoVCF(**record)

            # Get or create the Variant for this record.
            # NOTE: No samples so query_cache is not necessary.
            variant, alts = get_or_create_variant(
                    reference_genome, record, dataset, query_cache=None)

            # Create a link between the Variant and the VariantSet if
            # it doesn't exist.
            VariantToVariantSet.objects.get_or_create(
                    variant=variant,
                    variant_set=variant_set)


##############################################################################
# Helper Functions
##############################################################################

def copy_and_add_dataset_source(entity, dataset_label, dataset_type,
        original_source_location, move=False):
    """Copies the dataset to the entity location and then adds as
    Dataset. If the original_source_location is a file object, then
    it just read()s from the handle and writes to destination.

    If move is true, move instead of copying it. Good for files downloaded
    to a temp directory, since copying is slower.

    The model entity must satisfy the following interface:
        * property dataset_set
        * method get_model_data_dir()

    Returns:
        The Dataset object.
    """
    dest = copy_dataset_to_entity_data_dir(entity, original_source_location,
        move)
    dataset = add_dataset_to_entity(entity, dataset_label, dataset_type,
        dest)

    # First create the dataset and set copying status on it, dont set
    # the filesystem location.
    # dest = _get_copy_target_path(entity, original_source_location)
    # dataset = add_dataset_to_entity(entity, dataset_label, dataset_type, dest)
    # actual_dest = copy_dataset_to_entity_data_dir(
    #         entity, original_source_location, move=move)
    # assert actual_dest == dest, "If this fails, there's a bug."

    return dataset


def _get_copy_target_path(entity, original_source_location):
    """Returns the full path to the copy target.

    Args:
        entity: Model entity from which we determine the target dir.
        original_source_location: Original location from which we determine a
            filename.

    Returns:
        String describing full target path.
    """
    assert hasattr(entity, 'get_model_data_dir')
    source_name = os.path.split(original_source_location)[1]
    return os.path.join(entity.get_model_data_dir(), source_name)


def copy_dataset_to_entity_data_dir(entity, original_source_location,
        move=False):
    """If a file path, copy the data to the entity model data dir.
    If a handle, then just write it to the data dir.

    Returns:
        The destination to which the file was copied.
    """
    dest = _get_copy_target_path(entity, original_source_location)

    if not original_source_location == dest:
        try: #first try path
            if move:
                shutil.move(original_source_location, dest)
            else:
                shutil.copy(original_source_location, dest)
        except TypeError: #then try a handle
            open(dest,'w').write(
                original_source_location.read())
    return dest


def add_dataset_to_entity(entity, dataset_label, dataset_type,
        filesystem_location=None):
    """Helper function for adding a Dataset to a model.
    """
    dataset = Dataset.objects.create(
            label=dataset_label, type=dataset_type)

    if filesystem_location is not None:
        dataset.filesystem_location = clean_filesystem_location(
                filesystem_location)
        dataset.save()

    entity.dataset_set.add(dataset)
    entity.save()

    return dataset


def prepare_ref_genome_related_datasets(ref_genome, dataset):
    """Prepares data related to a ReferenceGenome.

    For example, if only Genbank exists, creates a Fasta Dataset.

    If related Datasets exists, this function is a no-op.

    Args:
        ref_genome: ReferenceGenome.
        dataset: A dataset pointing to a genome.

    Raises:
        AssertionError if dataset status is NOT_STARTED.
    """
    assert dataset.status != Dataset.STATUS.NOT_STARTED

    if dataset.type == Dataset.TYPE.REFERENCE_GENOME_FASTA:

        # make sure the fasta index is generated


        # Run jbrowse ref genome processing
        prepare_jbrowse_ref_sequence(ref_genome)

    elif dataset.type == Dataset.TYPE.REFERENCE_GENOME_GENBANK:
        # Run snpeff build after creating ReferenceGenome obj.
        build_snpeff(ref_genome)

        # These functions are NO-OPS if the respective Datasets exist.
        generate_fasta_from_genbank(ref_genome)
        generate_gff_from_genbank(ref_genome)

        # Run jbrowse genbank genome processing for genes
        add_genbank_file_track(ref_genome)

    # We create the bwa index once here, so that alignments running in
    # parallel don't step on each others' toes.
    ref_genome_fasta = get_dataset_with_type(ref_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()
    ensure_bwa_index(ref_genome_fasta)
