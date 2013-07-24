"""
Methods related to importing data.
"""

import os
import shutil

from main.models import Dataset
from main.models import ReferenceGenome

from Bio import SeqIO


IMPORT_FORMAT_TO_DATASET_TYPE = {
    'fasta': Dataset.TYPE.REFERENCE_GENOME_FASTA,
    'genbank': Dataset.TYPE.REFERENCE_GENOME_GENBANK
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
    _copy_and_add_dataset_source(reference_genome, dataset_type,
            dataset_type, file_location)
    return reference_genome


def _copy_and_add_dataset_source(entity, dataset_label, dataset_type,
        original_source_location):
    """Copies the dataset to the entity location and then adds as
    Dataset.

    Returns the Dataset object.
    """
    source_name = os.path.split(original_source_location)[1]
    dest = os.path.join(entity.get_model_data_dir(), source_name)
    if not original_source_location == dest:
        shutil.copy(original_source_location, dest)
    return _add_dataset_to_genome(
            entity, dataset_label, dataset_type, dest)


def _add_dataset_to_genome(genome, dataset_label, dataset_type,
        filesystem_location):
    """Helper function for adding a Dataset to a ReferenceGenome.

    Returns the Dataset object.
    """
    dataset = Dataset.objects.create(
            label=dataset_label,
            type=dataset_type,
            filesystem_location=filesystem_location)
    genome.dataset_set.add(dataset)
    genome.save()
    return dataset
