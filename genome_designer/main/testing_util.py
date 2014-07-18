"""
Utility methods for testing.
"""

from django.contrib.auth.models import User

from main.models import AlignmentGroup
from main.models import Chromosome
from main.models import ExperimentSample
from main.models import Project
from main.models import ReferenceGenome


TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'
TEST_PROJECT_NAME = 'recoli'
TEST_REF_GENOME_LABEL = 'mg1655'


def create_common_entities():
    """Creates the most common entities for testing.

    Returns at a User, Project, ReferenceGenome that are all
    related.
    """
    user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)

    project = Project.objects.create(
            title=TEST_PROJECT_NAME, owner=user.get_profile())

    reference_genome = ReferenceGenome.objects.create(
            project=project,
            label=TEST_REF_GENOME_LABEL)

    chromosome = Chromosome.objects.create(
            reference_genome=reference_genome,
            label='Chromosome',
            num_bases=9001)

    sample_1 = ExperimentSample.objects.create(
            project=project,
            label='es1')

    sample_2 = ExperimentSample.objects.create(
            project=project,
            label='es2')

    alignment_group_1 = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=reference_genome,
            aligner=AlignmentGroup.ALIGNER.BWA)

    return {
        'user': user,
        'project': project,
        'reference_genome': reference_genome,
        'chromosome': chromosome,
        'sample_1': sample_1,
        'sample_2': sample_2,
        'alignment_group_1': alignment_group_1
    }
