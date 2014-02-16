"""
Tests for pipeline.py
"""

import json
import os

from django.contrib.auth.models import User
from django.test import TestCase
from django.test.utils import override_settings

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import get_dataset_with_type
from main.models import ExperimentSample
from main.models import Project
from pipeline.pipeline import run_pipeline_multiple_ref_genomes
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import import_reference_genome_from_local_file
from scripts.jbrowse_util import prepare_jbrowse_ref_sequence
import settings


TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'

TEST_FASTA  = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')

class TestAlignmentPipeline(TestCase):

    def setUp(self):
        user = User.objects.create_user(TEST_USERNAME, password=TEST_PASSWORD,
                email=TEST_EMAIL)
        self.project = Project.objects.create(owner=user.get_profile(),
                title='Test Project')

        # Create a ref genome.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')

        # Create a sample.
        self.experiment_sample = ExperimentSample.objects.create(
                project=self.project, label='sample1')

        # Add fastq files to first sample.
        copy_and_add_dataset_source(self.experiment_sample, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, TEST_FASTQ1)
        copy_and_add_dataset_source(self.experiment_sample, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, TEST_FASTQ2)


    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS = True,
            CELERY_ALWAYS_EAGER = True, BROKER_BACKEND = 'memory')

    def test_run_pipeline(self):
        """Tests creating an alignment group.
        """
        ref_genome_list = [self.reference_genome]
        sample_list = [self.experiment_sample]

        run_pipeline_multiple_ref_genomes('name_placeholder',
                ref_genome_list, sample_list, test_models_only=False)

        alignment_group_obj_list = AlignmentGroup.objects.filter(
                reference_genome=self.reference_genome)
        self.assertEqual(1, len(alignment_group_obj_list))

        alignment_group_obj = alignment_group_obj_list[0]
        self.assertEqual(1,
                len(alignment_group_obj.experimentsampletoalignment_set.all()))

        # Make sure the initial JBrowse config has been created.
        jbrowse_dir = self.reference_genome.get_jbrowse_directory_path()
        self.assertTrue(os.path.exists(jbrowse_dir))
        self.assertTrue(os.path.exists(os.path.join(jbrowse_dir,
                'trackList.json')))
