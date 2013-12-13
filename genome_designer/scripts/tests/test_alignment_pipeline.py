"""
Tests for alignment_pipeline.py
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
from scripts.alignment_pipeline import align_with_bwa
from scripts.alignment_pipeline import align_with_bwa_mem
from scripts.alignment_pipeline import create_alignment_groups_and_start_alignments
from scripts.alignment_pipeline import FILES_TO_DELETE_AFTER_ALIGNMENT
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

TEST_FASTQ1_GZ = os.path.join(settings.PWD, 'test_data', 'compressed_fastq',
        'sample0.simLibrary.1.fq.gz')
TEST_FASTQ2_GZ = os.path.join(settings.PWD, 'test_data', 'compressed_fastq',
        'sample0.simLibrary.2.fq.gz')


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

        # Create a sample for compressed fastq data.
        self.compressed_experiment_sample = ExperimentSample.objects.create(
                project=self.project, label='sample1')

        # Add fastq files to first sample.
        copy_and_add_dataset_source(self.experiment_sample, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, TEST_FASTQ1)
        copy_and_add_dataset_source(self.experiment_sample, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, TEST_FASTQ2)

        # Add compressed fastq files to second sample.
        copy_and_add_dataset_source(self.compressed_experiment_sample,
                Dataset.TYPE.FASTQ1, Dataset.TYPE.FASTQ1, TEST_FASTQ1_GZ)
        copy_and_add_dataset_source(self.compressed_experiment_sample,
                Dataset.TYPE.FASTQ2, Dataset.TYPE.FASTQ2, TEST_FASTQ2_GZ)


    def test_bwa_align(self):
        """Test a single BWA alignment.
        """
        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        # Run the alignment.
        experiment_sample_alignment = align_with_bwa(
                alignment_group, self.experiment_sample)

        # Check that the run was successful as indicatedby the Dataset status.

        # UNCOMMENT FOR DEBUG
        # error = get_dataset_with_type(
        #         experiment_sample_alignment,
        #         Dataset.TYPE.BWA_ALIGN_ERROR).get_absolute_location()
        # with open(error) as fh:
        #     print fh.read()
        # assert False

        bwa_align_dataset = get_dataset_with_type(
                experiment_sample_alignment, Dataset.TYPE.BWA_ALIGN)
        self.assertEqual(Dataset.STATUS.READY, bwa_align_dataset.status)

        # Check that the alignment data was saved to a valid destination.
        bwa_align_dataset_path = bwa_align_dataset.get_absolute_location()
        self.assertTrue(os.path.exists(bwa_align_dataset_path,), (
                "No file at location %s" % bwa_align_dataset_path))

        # Check that redundant data is gone.
        for f in os.listdir(experiment_sample_alignment.get_model_data_dir()):
            self.assertFalse(f in FILES_TO_DELETE_AFTER_ALIGNMENT)

        # Check that a JBrowse track for the alignment was created.
        jbrowse_dir = self.reference_genome.get_jbrowse_directory_path()
        with open(os.path.join(jbrowse_dir, 'trackList.json')) as fh:
            jbrowse_config_json = json.loads(fh.read())
            found_bam_track = False
            for track in jbrowse_config_json['tracks']:
                if ('storeClass' in track and track['storeClass'] ==
                        'JBrowse/Store/SeqFeature/BAM'):
                    found_bam_track = True

                    # Also check that the urlTemplate is correct.
                    EXPECTED_URL_TEMPLATE = os.path.join(
                            settings.JBROWSE_DATA_URL_ROOT,
                            bwa_align_dataset.filesystem_location)
                    self.assertEqual(EXPECTED_URL_TEMPLATE,
                            track['urlTemplate'])

                    break
            self.assertTrue(found_bam_track)


    def test_bwa_align_mem(self):
        """Test a single BWA alignment.
        """
        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        # Run the alignment.
        experiment_sample_alignment = align_with_bwa_mem(
                alignment_group, self.experiment_sample)

        # Check that the run was successful as indicatedby the Dataset status.

        # UNCOMMENT FOR DEBUG
        # error = get_dataset_with_type(
        #         experiment_sample_alignment,
        #         Dataset.TYPE.BWA_ALIGN_ERROR).get_absolute_location()
        # with open(error) as fh:
        #     print fh.read()
        # assert False

        bwa_align_dataset = get_dataset_with_type(
                experiment_sample_alignment, Dataset.TYPE.BWA_ALIGN)
        self.assertEqual(Dataset.STATUS.READY, bwa_align_dataset.status)

        # Check that the alignment data was saved to a valid destination.
        bwa_align_dataset_path = bwa_align_dataset.get_absolute_location()
        self.assertTrue(os.path.exists(bwa_align_dataset_path,), (
                "No file at location %s" % bwa_align_dataset_path))

        # Check that a JBrowse track for the alignment was created.
        jbrowse_dir = self.reference_genome.get_jbrowse_directory_path()
        with open(os.path.join(jbrowse_dir, 'trackList.json')) as fh:
            jbrowse_config_json = json.loads(fh.read())
            found_bam_track = False
            for track in jbrowse_config_json['tracks']:
                if ('storeClass' in track and track['storeClass'] ==
                        'JBrowse/Store/SeqFeature/BAM'):
                    found_bam_track = True

                    # Also check that the urlTemplate is correct.
                    # NOTE (Changping): This would fail if the project is 
                    # S3 backed and urlTemplate would point to S3 endpoint.

                    # EXPECTED_URL_TEMPLATE = os.path.join(
                    #         settings.JBROWSE_DATA_URL_ROOT,
                    #         bwa_align_dataset.filesystem_location)
                    # self.assertEqual(EXPECTED_URL_TEMPLATE,
                    #         track['urlTemplate'])

                    break
            self.assertTrue(found_bam_track)

    def test_compressed_bwa_align(self):
        """Test a single BWA alignment.
        """
        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', 
                reference_genome=self.reference_genome)

        # Run the alignment.
        experiment_sample_alignment = align_with_bwa_mem(
                alignment_group, self.compressed_experiment_sample)

        # Check that the run was successful as indicatedby the Dataset status.

        # UNCOMMENT FOR DEBUG
        # error = get_dataset_with_type(
        #         experiment_sample_alignment,
        #         Dataset.TYPE.BWA_ALIGN_ERROR).get_absolute_location()
        # with open(error) as fh:
        #     print fh.read()
        # assert False

        bwa_align_dataset = get_dataset_with_type(
                experiment_sample_alignment, Dataset.TYPE.BWA_ALIGN)
        self.assertEqual(Dataset.STATUS.READY, bwa_align_dataset.status)

        # Check that the alignment data was saved to a valid destination.
        bwa_align_dataset_path = bwa_align_dataset.get_absolute_location()
        self.assertTrue(os.path.exists(bwa_align_dataset_path,), (
                "No file at location %s" % bwa_align_dataset_path))

        # Check that a JBrowse track for the alignment was created.
        jbrowse_dir = self.reference_genome.get_jbrowse_directory_path()
        with open(os.path.join(jbrowse_dir, 'trackList.json')) as fh:
            jbrowse_config_json = json.loads(fh.read())
            found_bam_track = False
            for track in jbrowse_config_json['tracks']:
                if ('storeClass' in track and track['storeClass'] ==
                        'JBrowse/Store/SeqFeature/BAM'):
                    found_bam_track = True

                    # Also check that the urlTemplate is correct.
                    EXPECTED_URL_TEMPLATE = os.path.join(
                            settings.JBROWSE_DATA_URL_ROOT,
                            bwa_align_dataset.filesystem_location)
                    self.assertEqual(EXPECTED_URL_TEMPLATE,
                            track['urlTemplate'])

                    break
            self.assertTrue(found_bam_track)


    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS = True,
            CELERY_ALWAYS_EAGER = True, BROKER_BACKEND = 'memory')
    def test_create_alignment_groups_and_start_alignments(self):
        """Tests creating an alignment group.
        """
        ref_genome_list = [self.reference_genome]
        sample_list = [self.experiment_sample]
        create_alignment_groups_and_start_alignments('name_placeholder',
                ref_genome_list, sample_list, test_models_only=True,
                concurrent=True)

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
