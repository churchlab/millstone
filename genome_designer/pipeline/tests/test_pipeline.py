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
from main.models import Variant
from main.testing_util import FullVCFTestSet
from pipeline.pipeline_runner import run_pipeline
from utils.import_util import copy_and_add_dataset_source
from utils.import_util import import_reference_genome_from_local_file


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

        # Create a sample with a single fastq
        self.experiment_sample_single_fastq = ExperimentSample.objects.create(
                project=self.project, label='sample_single_fastq')

        # Add the fastq file to the sample
        copy_and_add_dataset_source(self.experiment_sample_single_fastq,
                Dataset.TYPE.FASTQ1, Dataset.TYPE.FASTQ1, TEST_FASTQ1)

    def test_run_pipeline(self):
        """End-to-end test of pipeline. Fails if any errors.
        """
        sample_list = [self.experiment_sample]
        result = run_pipeline(
                'name_placeholder', self.reference_genome, sample_list)
        alignment_group = result[0]
        alignment_async_result = result[1]
        variant_calling_async_result = result[2]
        alignment_async_result.get()
        variant_calling_async_result.get()
        alignment_group = AlignmentGroup.objects.get(uid=alignment_group.uid)
        self.assertEqual(AlignmentGroup.STATUS.COMPLETED,
                alignment_group.status)

        # Make sure some expected variants are found.
        variants = Variant.objects.filter(
                reference_genome=self.reference_genome)
        self.assertTrue(len(variants))
        v_1834 = Variant.objects.get(position=1834)
        v_1834_vccd = v_1834.variantcallercommondata_set.all()[0]
        v_1834_ve = v_1834_vccd.variantevidence_set.all()[0]
        self.assertFalse(v_1834_ve.data.get('IS_SV', False))

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

    def test_run_pipeline__single_fastq(self):
        """End-to-end test of pipeline with single fastq sample.
        Fails if any errors.
        """
        sample_list = [self.experiment_sample_single_fastq]
        result = run_pipeline(
                'name_placeholder', self.reference_genome, sample_list)
        alignment_group = result[0]
        alignment_async_result = result[1]
        variant_calling_async_result = result[2]
        alignment_async_result.get()
        variant_calling_async_result.get()
        alignment_group = AlignmentGroup.objects.get(uid=alignment_group.uid)
        self.assertEqual(AlignmentGroup.STATUS.COMPLETED,
                alignment_group.status)

    def test_run_pipeline__snps_with_effect__no_svs(self):
        """Tests pipeline with SNPs with effect, but no SVs called.
        """
        ref_genome = import_reference_genome_from_local_file(
                self.project, 'mg1655_tolC_through_zupT',
                FullVCFTestSet.TEST_GENBANK, 'genbank')

        sample_obj = ExperimentSample.objects.create(
                project=self.project,
                label='Sample %d' % 0)

        # Add raw reads to each sample.
        copy_and_add_dataset_source(sample_obj,
                Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1,
                FullVCFTestSet.FASTQ1[0])
        copy_and_add_dataset_source(sample_obj,
                Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2,
                FullVCFTestSet.FASTQ2[0])

        result = run_pipeline(
            'test_align', ref_genome, [sample_obj])

        alignment_group = result[0]
        alignment_async_result = result[1]
        variant_calling_async_result = result[2]
        alignment_async_result.get()
        variant_calling_async_result.get()
        alignment_group = AlignmentGroup.objects.get(uid=alignment_group.uid)
        self.assertEqual(AlignmentGroup.STATUS.COMPLETED,
                alignment_group.status)

        # Check that SnpEff worked.
        v_205 = Variant.objects.get(
                reference_genome=alignment_group.reference_genome, position=205)
        v_205_va = v_205.variantalternate_set.all()[0]
        self.assertEqual('tolC', v_205_va.data['INFO_EFF_GENE'])
