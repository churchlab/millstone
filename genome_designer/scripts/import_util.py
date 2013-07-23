"""
Methods related to importing data.
"""

from main.models import Dataset
from main.models import ReferenceGenome


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

    # Create the ReferenceGenome object.
    reference_genome = ReferenceGenome.objects.create(project=project,
            label=label)

    # TODO: Finish Implementing.

    return reference_genome
