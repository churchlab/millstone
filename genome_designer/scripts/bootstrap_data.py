#!/usr/bin/env python

"""
Script to setup some test data.

This is useful during development when we are continuously wiping the db
and want to get some new data in quickly.
"""

import os
import random
import shutil

from util import setup_django_env


# This is the directory where this bootstrap script is located.
PWD = os.path.dirname(os.path.realpath(__file__ ))


def bootstrap_fake_data():
    """Fill the database with fake data.
    """

    # Imports only work after the environment has been set up.
    from django.contrib.auth.models import User

    ### Get or create the user.
    TEST_USERNAME = 'gmcdev'
    TEST_PASSWORD = 'g3n3d3z'
    TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'
    try:
        user = User.objects.get(username=TEST_USERNAME)
    except User.DoesNotExist:
        user = User.objects.create_user(
                TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)

    ### Create some projects
    from main.models import Project
    TEST_PROJECT_NAME = 'recoli'
    (test_project, project_created) = Project.objects.get_or_create(
            title=TEST_PROJECT_NAME, owner=user.get_profile())
    (test_project_2, project_created) = Project.objects.get_or_create(
            title='project2', owner=user.get_profile())
    (test_project_3, project_created) = Project.objects.get_or_create(
            title='project3', owner=user.get_profile())

    ### Create some reference genomes
    from main.models import ReferenceGenome
    REF_GENOME_1_LABEL = 'mg1655'
    (ref_genome_1, ref_genome_created) = ReferenceGenome.objects.get_or_create(
            label=REF_GENOME_1_LABEL, project=test_project, num_chromosomes=1,
            num_bases=100)
    REF_GENOME_2_LABEL = 'c321D'
    (ref_genome_2, ref_genome_created) = ReferenceGenome.objects.get_or_create(
            label=REF_GENOME_2_LABEL, project=test_project, num_chromosomes=1,
            num_bases=200)

    ### Create some samples
    from main.models import ExperimentSample
    SAMPLE_1_LABEL = 'sample1'
    (sample_1, created) = ExperimentSample.objects.get_or_create(
            project=test_project,
            label=SAMPLE_1_LABEL)

    ### Add datasets to the samples.
    from main.models import Dataset

    dataset_1 = Dataset.objects.create(
            type=Dataset.TYPE.FASTQ1,
            label='sample1_fastq1',
            filesystem_location='blah')
    sample_1.dataset_set.add(dataset_1)

    dataset_2 = Dataset.objects.create(
            type=Dataset.TYPE.FASTQ2,
            label='sample1_fastq2',
            filesystem_location='blah2')
    sample_1.dataset_set.add(dataset_2)

    ### Create an alignment.
    from main.models import AlignmentGroup
    alignment_group_1 = AlignmentGroup.objects.create(
            label='Alignment 1',
            reference_genome=ref_genome_1,
            aligner=AlignmentGroup.ALIGNER.BWA)

    # Link it to a sample.
    from main.models import ExperimentSampleToAlignment
    ExperimentSampleToAlignment.objects.create(
            alignment_group=alignment_group_1,
            experiment_sample=sample_1)


    ### Create some fake variants
    from main.models import Variant
    for var_count in range(100):
        Variant.objects.create(
            type=Variant.TYPE.TRANSITION,
            reference_genome=ref_genome_1,
            chromosome='chrom',
            position=random.randint(1,ref_genome_1.num_bases),
            ref_value='A',
            alt_value='G')


def reset_database():
    """Deletes the old database and sets up a new one.

    For now, only works with the temp.db database to prevent
    accidentally deleting data down the line.
    """
    import settings

    ### Delete the old database if it exists.
    print 'Deleting old database ...'
    TEMP_DB_NAME = 'temp.db'
    temp_db_name = settings.DATABASES['default']['NAME']
    assert temp_db_name == TEMP_DB_NAME
    tempdb_dir = os.path.split(PWD)[0]
    tempdb_abs_path = os.path.join(tempdb_dir, TEMP_DB_NAME)
    if os.path.exists(tempdb_abs_path):
        os.remove(tempdb_abs_path)

    ### Run syncdb
    from django.core.management import call_command
    # NOTE: Remove interactive=False if you want to have the option of creating
    # a super user on sync.
    print 'Creating new database via syncdb ...'
    call_command('syncdb', interactive=False)

    ### Recreate the media root.
    if os.path.exists(settings.MEDIA_ROOT):
        shutil.rmtree(settings.MEDIA_ROOT)
    os.mkdir(settings.MEDIA_ROOT)


if __name__ == '__main__':
    setup_django_env()
    reset_database()
    bootstrap_fake_data()
