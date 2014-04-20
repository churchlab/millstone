"""
Tests for pipeline_runner.py
"""

import os

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TransactionTestCase

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import Project
from pipeline.pipeline_runner import run_pipeline
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import import_reference_genome_from_local_file
from scripts.import_util import import_reference_genome_from_ncbi
from scripts.util import internet_on


TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'

TEST_FASTA = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')


class TestAlignmentPipeline(TransactionTestCase):
    """Tests for the pipeline.

    NOTE: We use TransactoinTestCase since we are transitioning to treating
    this as an integration test where we use celery, where it's necessary
    that the database changes are actually committed for the celery thread
    to actually see them.
    """

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


    def test_run_pipeline__samples_not_ready__fastq1(self):
        """Tests that the pipeline raises an AssertionError if samples aren't
        ready, fastq1.
        """
        fastq_dataset = self.experiment_sample.dataset_set.filter(
            type=Dataset.TYPE.FASTQ1)[0]
        fastq_dataset.status = Dataset.STATUS.QUEUED_TO_COPY
        fastq_dataset.save()

        sample_list = [self.experiment_sample]

        with self.assertRaises(AssertionError):
            run_pipeline('name_placeholder',
                    self.reference_genome, sample_list)


    def test_run_pipeline__samples_not_ready__fastq2(self):
        """Tests that the pipeline raises an AssertionError if samples aren't
        ready, fastq2.
        """
        fastq_dataset = self.experiment_sample.dataset_set.filter(
            type=Dataset.TYPE.FASTQ2)[0]
        fastq_dataset.status = Dataset.STATUS.QUEUED_TO_COPY
        fastq_dataset.save()

        sample_list = [self.experiment_sample]

        with self.assertRaises(AssertionError):
            run_pipeline('name_placeholder', self.reference_genome, sample_list)


    def test_run_pipeline__genbank_from_ncbi_with_spaces_in_label(self):
        """Tests the pipeline where the genome is imported from NCBI with
        spaces in the name.

        NOTE: This should really be an integration test.
        """
        if not internet_on():
            return
        MG1655_ACCESSION = 'NC_000913.3'
        MG1655_LABEL = 'mg1655 look a space'
        ref_genome = import_reference_genome_from_ncbi(self.project,
                MG1655_LABEL, MG1655_ACCESSION, 'genbank')
        sample_list = [self.experiment_sample]

        run_pipeline('name_placeholder', ref_genome, sample_list)

        alignment_group_obj_list = AlignmentGroup.objects.filter(
                reference_genome=ref_genome)
        self.assertEqual(1, len(alignment_group_obj_list))

        alignment_group_obj = alignment_group_obj_list[0]
        self.assertEqual(1,
                len(alignment_group_obj.experimentsampletoalignment_set.all()))
