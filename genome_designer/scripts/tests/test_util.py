"""
Utility methods for testing.
"""

from django.contrib.auth.models import User

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
            label=TEST_REF_GENOME_LABEL,
            num_chromosomes=1,
            num_bases=777)

    return {
        'user': user,
        'project': project,
        'reference_genome': reference_genome
    }