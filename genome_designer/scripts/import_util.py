"""
Methods related to importing data.
"""

import copy
import csv
import os
import shutil
import re
import vcf
from collections import namedtuple

from django.db import transaction
from main.models import clean_filesystem_location
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome
from main.models import Variant
from main.models import VariantCallerCommonData
from main.models import VariantSet
from main.models import VariantToVariantSet
from scripts.vcf_parser import extract_raw_data_dict
from scripts.vcf_parser import get_or_create_variant
from settings import PWD

from Bio import SeqIO


IMPORT_FORMAT_TO_DATASET_TYPE = {
    'fasta': Dataset.TYPE.REFERENCE_GENOME_FASTA,
    'genbank': Dataset.TYPE.REFERENCE_GENOME_GENBANK,
    'vcfu': Dataset.TYPE.VCF_USERINPUT
}

def import_reference_genome_from_local_file(project, label, file_location,
        import_format):
    """Creates a ReferenceGenome associated with the given Project.

    Args:
        project: The Project we're storing everyting relative to.
        label: The human-readable label for the ReferenceGenome.
        file_location: Location of the genome on the server.
        import_format: Must be 'fasta' or 'genbank'.
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
    genbank_path = ref_genome.dataset_set.get(
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()

    genbank_dir, genbank_filename = os.path.split(genbank_path)
    genbank_noext = os.path.splitext(genbank_filename)[0]

    # Put the fasta file in the same dir, just change the extension to .fa.
    fasta_filename = os.path.join(genbank_dir, (genbank_noext + '.fa'))

    # Get the individual records, each corresponding to a chromosome.
    genome_records = list(SeqIO.parse(genbank_path, 'genbank'))

    SeqIO.write(genome_records, fasta_filename, 'fasta')

    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE['fasta']
    copy_and_add_dataset_source(ref_genome, dataset_type,
            dataset_type, fasta_filename)

    return

def sanitize_record_id(record_id_string):
    """We want to grab only the first word-only part of each seqrecord in a
    FASTA/Genbank file, and use that as a consistent and readable id between
    genbank and FASTA.
    """
    return re.match( r'^\w{1,20}', record_id_string).group()


def import_samples_from_targets_file(project, targets_file):
    """Uses the uploaded targets file to add a set of samples to the project.
    We need to check each line of the targets file for consistency before we
    do anything, however.

    It writes a copy of the uploaded targets file to a temporary file

    Args:
        project: The project we're storing everything relative to>
        targets_file: The UploadedFile django object that holds the targets
            in .tsv format.
    """
    # The targets file shouldn't be over 1 Mb ( that would be ~3,000 genomes)
    assert targets_file.size < 1000000, "Targets file is too large."

    # Detect the format.
    reader = csv.DictReader(targets_file, delimiter='\t')

    # Validate the header.
    targets_file_header = reader.fieldnames

    assert len(targets_file_header) >= 6, "Bad header. Were columns removed?"

    REQUIRED_HEADER_PART = ['Sample_Name', 'Plate_or_Group', 'Well',
            'Read_1_Path', 'Read_2_Path','Parent_Samples']
    for col, check in zip(targets_file_header[0:len(REQUIRED_HEADER_PART)],
            REQUIRED_HEADER_PART):
        assert col == check, (
            "Header column '%s' is missing or out of order." % check)

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

    # Now create ExperimentSample objects along with their respective Datasets.
    # The data is copied to the entity location.
    for row in valid_rows:
        # Create ExperimentSample object and then store the data relative to
        # it.
        sample_label = row['Sample_Name']
        experiment_sample = ExperimentSample.objects.create(
                project=project, label=sample_label)
        copy_and_add_dataset_source(experiment_sample, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, row['Read_1_Path'])
        if 'Read_2_Path' in row and row['Read_2_Path']:
            copy_and_add_dataset_source(experiment_sample, Dataset.TYPE.FASTQ2,
                    Dataset.TYPE.FASTQ2, row['Read_2_Path'])


@transaction.commit_on_success
def import_variant_set_from_vcf(ref_genome, variant_set_name, variant_set_file):
    """Convert an uploaded VCF file into a new variant set object.

    Args:
        ref_genome: ReferenceGenome.
        variant_set_name: Name of the variant set (label).
        variant_set_file: Path to the variant set on disk.
    """
    # For now, variant set name must be unique even among diff ref genomes.
    variant_set_name_exists = len(VariantSet.objects.filter(
            label=variant_set_name)) > 0
    assert not variant_set_name_exists, 'Variant set name must be unique.'

    # Create the VariantSet.
    variant_set = VariantSet.objects.create(
            reference_genome=ref_genome,
            label=variant_set_name)

    # First, save this vcf as a dataset, so we can point to it from the
    # new variant common_data_objs
    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE['vcfu']
    dataset = copy_and_add_dataset_source(variant_set, dataset_type,
            dataset_type, variant_set_file)

    # Now read the variant set file.
    _read_variant_set_file_as_csv(variant_set_file, ref_genome, dataset,
            variant_set)


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
        REQUIRED_HEADER_PART = ['CHROM','POS','ID','REF','ALT']

        for col, check in zip(reader.fieldnames[0:len(REQUIRED_HEADER_PART)],
                REQUIRED_HEADER_PART):
            assert col == check, (
                "Header column '%s' is missing or out of order; %s" % (check,
                    ', '.join(reader.fieldnames)))

        class PseudoVCF:
            """Pseudo wrapper class to satisfy interface of
            extract_raw_data_dict().
            """
            def __init__(self, **entries):
                self.__dict__.update(entries)

                # Currently, we are storing alts as a string representing
                # the pickled/unpickled value parsed by pyvcf, so we
                # match that here.
                self.__dict__['ALT'] = '[' + self.__dict__['ALT'].strip() + ']'

        for record in reader:
            record = PseudoVCF(**record)

            # Build a dictionary of data for this record.
            raw_data_dict = extract_raw_data_dict(record)

            # Get or create the Variant for this record.
            variant = get_or_create_variant(reference_genome, raw_data_dict)

            # Create a common data object for this variant.
            common_data_obj = VariantCallerCommonData.objects.create(
                    variant=variant,
                    source_dataset=dataset,
                    data=raw_data_dict
            )

            # Create a link between the Variant and the VariantSet if
            # it doesn't exist.
            VariantToVariantSet.objects.get_or_create(
                    variant=variant,
                    variant_set=variant_set)


def copy_and_add_dataset_source(entity, dataset_label, dataset_type,
        original_source_location):
    """Copies the dataset to the entity location and then adds as
    Dataset.

    The model entity must satisfy the following interface:
        * property dataset_set
        * method get_model_data_dir()

    Returns the Dataset object.
    """
    dest = copy_dataset_to_entity_data_dir(entity, original_source_location)
    dataset = add_dataset_to_entity(entity, dataset_label, dataset_type,
        dest)

    return dataset


def copy_dataset_to_entity_data_dir(entity, original_source_location):
    """Copies the data to the entity model data dir.

    Returns:
        The destination to which the file was copied.
    """
    assert hasattr(entity, 'get_model_data_dir')
    source_name = os.path.split(original_source_location)[1]
    dest = os.path.join(entity.get_model_data_dir(), source_name)
    if not original_source_location == dest:
        shutil.copy(original_source_location, dest)
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
