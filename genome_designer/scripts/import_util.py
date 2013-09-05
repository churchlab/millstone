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
from scripts.vcf_parser import populate_common_data_from_vcf_record
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
    assert os.path.exists(file_location), "File doesn't exist."

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
def import_variant_set_from_vcf(project, ref_genome_uid, variant_set_name,
        variant_set_file):
    """Convert an uploaded VCF file into a new variant set object.

    Args:
        project: The Project this set is part of.
        ref_genome_uid: ReferenceGenome uid.
        variant_set_name: Name of the variant set (label).
        variant_set_file: Path to the variant set on disk.
    """

    ref_genome = ReferenceGenome.objects.get(
            project=project,
            uid=ref_genome_uid)

    # For now, variant set name must be unique even among diff ref genomes
    variant_set_name_exists = len(VariantSet.objects.filter(
            label=variant_set_name)) > 0

    assert not variant_set_name_exists, 'Variant set name must be unique.'

    variant_set = VariantSet.objects.create(
            reference_genome=ref_genome,
            label=variant_set_name)

    # First, save this vcf as a dataset, so we can point to it from the
    # new variant common_data_objs
    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE['vcfu']
    dataset = copy_and_add_dataset_source(variant_set, dataset_type,
            dataset_type, variant_set_file)

    # dbg 8/25/13 - this try/except fails, so I'm doing it only w/ csv below
    # First try with pyVCF, but if this fails then just use a csvReader
    # try:
    #     common_data_obj_iter = _read_variant_set_file_as_full_vcf(
    #             variant_set_file, ref_genome, dataset)
    # except IndexError:
    #     common_data_obj_iter = _read_variant_set_file_as_csv(
    #             variant_set_file, ref_genome, dataset)

    common_data_obj_iter = _read_variant_set_file_as_csv(
            variant_set_file, ref_genome, dataset)

    # Finally, create this variant if it doesn't already exist and add it
    # to the set.
    for common_data_obj in common_data_obj_iter:
        variant = get_or_create_variant(ref_genome, common_data_obj)
        VariantToVariantSet.objects.create(
                variant=variant,
                variant_set=variant_set)

def _read_variant_set_file_as_full_vcf(variant_set_file, reference_genome,
    dataset):
    """The importer will first try to read the VCF as an actual VCF, instead
    of a stripped down version with no samples, info, etc. We will use
    pyVCF for this purpose. If it fails, then we will resort to the CSV
    approach.

    Args:
        * variant_set_file - path to vcf file
        * reference_genome - django reference genome object
    Returns:
        an interator of common_data_objs used for creating or adding to
        variant objects
    """

    with open(variant_set_file) as fh:
       vcf_reader = vcf.Reader(fh)
       for record in vcf_reader:
           # Create a common data object for this variant.
           common_data_obj = VariantCallerCommonData.objects.create(
                   reference_genome=reference_genome,
                   source_dataset_id=dataset.id
           )
           # Extracting the data is somewhat of a process so we do
           # so in this function.
           populate_common_data_from_vcf_record(common_data_obj, record)

           yield common_data_obj

def _read_variant_set_file_as_csv(variant_set_file, reference_genome,
    dataset):
    """If reading the variant set file as a vcf fails (because we arent using
    all columns, as will usually be the case) then read it as a CSV and check
    manually for the required columns.

    Args:
        * variant_set_file - path to vcf file
        * reference_genome - django reference genome object
    Returns:
        an interator of common_data_objs used for creating or adding to
        variant objects
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

        #check that the required columns are present
        REQUIRED_HEADER_PART = ['CHROM','POS','ID','REF','ALT']

        for col, check in zip(reader.fieldnames[0:len(REQUIRED_HEADER_PART)],
                REQUIRED_HEADER_PART):
            assert col == check, (
                "Header column '%s' is missing or out of order; %s" % (check,
                    ', '.join(reader.fieldnames)))

        class PseudoVCF:
            def __init__(self, **entries):
                self.__dict__.update(entries)

        for record in reader:

            record = PseudoVCF(**record)

            # Create a common data object for this variant.
            common_data_obj = VariantCallerCommonData.objects.create(
                    reference_genome=reference_genome,
                    source_dataset_id=dataset.id
            )

            # Populate the common_data_obj with this row of the csv/vcf.
            # The function takes an object w/ attributes so we convert the
            # dict into one using a type() call.
            populate_common_data_from_vcf_record(common_data_obj, record)

            yield common_data_obj

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
        filesystem_location):
    """Helper function for adding a Dataset to a model.
    """
    dataset = Dataset.objects.create(
            label=dataset_label,
            type=dataset_type,
            filesystem_location=clean_filesystem_location(filesystem_location))

    entity.dataset_set.add(dataset)
    entity.save()

    return dataset
