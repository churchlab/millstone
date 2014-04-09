"""
Methods related to importing data.
"""

import copy
import csv
import os
import shutil
import re
import vcf
import tempfile
import traceback
import shutil
from collections import namedtuple
from tempfile import NamedTemporaryFile

from BCBio import GFF
from Bio import Entrez
from Bio import SeqIO
from celery import task
from django.db import transaction
from django.conf import settings

from main.celery_util import assert_celery_running
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantSet
from main.models import VariantToVariantSet
from main.models import get_dataset_with_type
from main.model_utils import clean_filesystem_location
from main.s3 import project_files_needed
from scripts.vcf_parser import extract_raw_data_dict
from scripts.vcf_parser import get_or_create_variant
from scripts.util import uppercase_underscore
from settings import PWD
from settings import EMAIL


IMPORT_FORMAT_TO_DATASET_TYPE = {
    'fasta': Dataset.TYPE.REFERENCE_GENOME_FASTA,
    'genbank': Dataset.TYPE.REFERENCE_GENOME_GENBANK,
    'gff': Dataset.TYPE.REFERENCE_GENOME_GFF,
    'vcf_user': Dataset.TYPE.VCF_USERINPUT
}

REQUIRED_SAMPLE_HEADER_PART = ['Sample_Name', 'Read_1_Path', 'Read_2_Path']
REQUIRED_VCF_HEADER_PART = ['CHROM','POS','ID','REF','ALT']

if settings.S3_ENABLED:
    from main.s3 import s3_temp_get, s3_get

    def import_reference_genome_from_s3(project, label, s3file, import_format):
        with s3_temp_get(s3file) as f:
            return import_reference_genome_from_local_file(
                    project, label, f, import_format)

    @project_files_needed
    def import_samples_from_s3(project, targets_file_rows, s3files):
        tmp_dir = tempfile.mkdtemp()
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
    num_chromosomes = 0
    num_bases = 0
    for genome_record in SeqIO.parse(file_location, import_format):
        num_chromosomes += 1
        num_bases += len(genome_record)

    # Make sure sequence exists.
    if not num_bases > 0:
        raise DataImportError("No sequence in file.")

    # Create the ReferenceGenome object.
    reference_genome = ReferenceGenome.objects.create(
            project=project,
            label=label,
            num_chromosomes=num_chromosomes,
            num_bases=num_bases)

    # Copy the source file to the ReferenceGenome data location.
    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE[import_format]
    copy_and_add_dataset_source(reference_genome, dataset_type,
            dataset_type, file_location)

    return reference_genome


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

    return

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

    return

def import_reference_genome_from_ncbi(project, label, record_id, import_format):
    """
    Pull a reference genome by accession from NCBI using efetch.
    """
    # Validate the input.
    assert import_format in ['fasta', 'genbank'], (
            'Import Format must be \'fasta\' or \'genbank\'')

    # Format keys for Efetch.
    # More info at:  http://www.ncbi.nlm.nih.gov/
    #       books/NBK25499/table/chapter4.chapter4_table1/?report=objectonly
    CONVERT_FORMAT = {'fasta':'fa', 'genbank':'gbwithparts'}
    # What suffix to use for each input format
    # TODO: Should this be a property of the Dataset TYPE?
    FORMAT_SUFFIX = {'fasta':'.fa', 'genbank':'.gb'}

    Entrez.email = EMAIL
    handle = Entrez.efetch(
            db="nuccore", 
            id=record_id, 
            rettype=CONVERT_FORMAT[import_format], 
            retmode="text")

    temp = NamedTemporaryFile(delete=False, prefix=label+'_', 
            suffix=FORMAT_SUFFIX[import_format])
    temp.write(handle.read())
    handle.close()
    temp.close()
    reference_genome = import_reference_genome_from_local_file(
        project, label, temp.name, import_format, move=True)

    if os.path.isfile(temp.name):
        os.remove(temp.name)

    return reference_genome

def sanitize_record_id(record_id_string):
    """We want to grab only the first word-only part of each seqrecord in a
    FASTA/Genbank file, and use that as a consistent and readable id between
    genbank and FASTA.
    """
    return re.match( r'^\w{1,20}', record_id_string).group()


def parse_targets_file(project, targets_file, remove_directory_path=False):
    # The targets file shouldn't be over 1 Mb ( that would be ~3,000 genomes)
    if hasattr(targets_file, "size"):
        assert targets_file.size < 1000000, (
                "Targets file is too large: %d" % targets_file.size)

    # Detect the format.
    reader = csv.DictReader(targets_file, delimiter='\t')

    targets_file_header = reader.fieldnames

    REQUIRED_SAMPLE_HEADER_STRING = ', '.join(REQUIRED_SAMPLE_HEADER_PART)

    assert len(targets_file_header) >= len(REQUIRED_SAMPLE_HEADER_PART), (
        "Header is too short. Need columns: %s" % REQUIRED_SAMPLE_HEADER_STRING)

    header_idxes = {}

    for col in REQUIRED_SAMPLE_HEADER_PART:
        if not col in targets_file_header:
            raise ValueError("Required header column '%s' is missing." % col)

    # Validate all the rows.
    valid_rows = []
    for row_num, row in enumerate(reader):

        # Make a copy of the row so we can clean up the data for further
        # processing.
        clean_row = copy.copy(row)

        #TODO: Every row seems to have an empty K/V pair {None:''}, not sure
        #why. Here I remove it by hand:
        row = dict([(k, v) for k, v in row.iteritems() if k is not None])

        # Make sure the row has all the fields
        assert len(targets_file_header) == len(row.keys()), (
                "Row %d has the wrong number of fields." % row_num)

        for field_name, field_value in row.iteritems():
            if 'Path' not in field_name:
                #make sure each field is alphanumeric only
                assert re.match('^[\. \w-]*$', field_value) is not None, (
                        'Only alphanumeric characters and spaces are allowed, '
                        'except for the paths.\n(Row %d, "%s")' % (
                                row_num, field_name))
            else:
                if remove_directory_path:
                    # If it is a path, take the filename from path and return
                    # them as a list.
                    clean_field_value = os.path.basename(field_value)
                else:
                    # If it is a path, then try to open the file and read one byte.
                    # Replace the string '$GD_ROOT with the project path, so we
                    # can use the test data
                    clean_field_value = field_value.replace('$GD_ROOT', PWD)
                    with open(clean_field_value, 'rb') as test_file:
                        try:
                            test_file.read(8)
                        except:
                            raise Exception("Cannot read file at %s" %
                                    clean_field_value)
                clean_row[field_name] = clean_field_value

        # Save this row.
        valid_rows.append(clean_row)
    return valid_rows


@project_files_needed
def import_samples_from_targets_file(project, targets_file):
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
    """
    assert_celery_running()

    valid_rows = parse_targets_file(project, targets_file)

    def _create_fastq_dataset(experiment_sample, fastq_source, dataset_type):
        """Helper function for creating a Dataset to track the uploaded
        sample data.
        """
        fastq_dest = _get_copy_target_path(experiment_sample, fastq_source)
        reads_dataset = add_dataset_to_entity(experiment_sample,
                dataset_type, dataset_type, fastq_dest)
        reads_dataset.status = Dataset.STATUS.QUEUED_TO_COPY
        reads_dataset.save()
        return reads_dataset

    # Now create ExperimentSample objects along with their respective Datasets.
    # The data is copied to the entity location.
    # We do this asynchronously so catch all the results, so we can block on
    # them.
    experiment_samples = []
    for row in valid_rows:
        # Create ExperimentSample object and then store the data relative to
        # it.
        sample_label = row['Sample_Name']
        experiment_sample = ExperimentSample.objects.create(
                project=project, label=sample_label)

        # Create the Datasets before starting copying so we can show status in
        # the ui. This is a new pattern where we are moving copying to happen
        # asynchronously in Celery.
        _create_fastq_dataset(
                experiment_sample, row['Read_1_Path'], Dataset.TYPE.FASTQ1)
        maybe_read2_path = row.get('Read_2_Path', '')
        if maybe_read2_path:
            _create_fastq_dataset(
                experiment_sample, maybe_read2_path, Dataset.TYPE.FASTQ2)

        # Start the async job of copying.
        copy_experiment_sample_data.delay(project, experiment_sample, row)

        # Add extra metadata columns.
        for field, value in row.iteritems():
            if field not in REQUIRED_SAMPLE_HEADER_PART:
                clean_field = 'SAMPLE_'+uppercase_underscore(field)
                experiment_sample.data[clean_field] = str(value)

        experiment_sample.save()
        experiment_samples.append(experiment_sample)

    return experiment_samples


@task
@project_files_needed
def copy_experiment_sample_data(project, experiment_sample, data):
    """Celery task that wraps the process of copying the data for an
    ExperimentSample.
    """
    move = False

    def _copy_dataset_data(fastq_source, dataset_type):
        """Helper to copy data and set status to VERIFYING.
        """
        reads_dataset = experiment_sample.dataset_set.get(type=dataset_type)
        reads_dataset.status = Dataset.STATUS.COPYING
        reads_dataset.save()
        copy_dataset_to_entity_data_dir(experiment_sample, fastq_source,
                move=move)
        reads_dataset.status = Dataset.STATUS.VERIFYING
        reads_dataset.save()
        return reads_dataset

    # Copy read1.
    read1_dataset = _copy_dataset_data(data['Read_1_Path'],
            Dataset.TYPE.FASTQ1)
    # Copy read2.
    maybe_read2_path = data.get('Read_2_Path', '')
    if maybe_read2_path:
        read2_dataset = _copy_dataset_data(maybe_read2_path,
                Dataset.TYPE.FASTQ2)
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

    # Save the status.
    read1_dataset.save()
    read2_dataset.save()


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
