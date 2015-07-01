"""
Tests for read_alignment.py
"""

import json
import os
import subprocess

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import get_dataset_with_type
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
from main.model_utils import clean_filesystem_location
from pipeline.read_alignment import align_with_bwa_mem
from pipeline.read_alignment import get_discordant_read_pairs
from pipeline.read_alignment import get_split_reads
from pipeline.read_alignment import get_read_length
from pipeline.read_alignment import get_insert_size_mean_and_stdev
from settings import TOOLS_DIR
from utils.import_util import copy_and_add_dataset_source
from utils.import_util import import_reference_genome_from_local_file
from utils.jbrowse_util import compile_tracklist_json


TEST_USERNAME = 'gmcdev'
TEST_PASSWORD = 'g3n3d3z'
TEST_EMAIL = 'gmcdev@genomedesigner.freelogy.org'

TEST_FASTA = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(settings.PWD, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')

TEST_FASTQ1_GZ = os.path.join(settings.PWD, 'test_data', 'compressed_fastq',
        'sample0.simLibrary.1.fq.gz')
TEST_FASTQ2_GZ = os.path.join(settings.PWD, 'test_data', 'compressed_fastq',
        'sample0.simLibrary.2.fq.gz')

TEST_DISC_SPLIT_BAM = os.path.join(settings.PWD, 'test_data',
        'discordant_split_reads',
        'bwa_align.bam')

SAMTOOLS_BINARY = '%s/samtools/samtools' % TOOLS_DIR

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


    def test_bwa_align_mem(self):
        """Test a single BWA alignment.
        """
        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        # Create the expected models.
        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=self.experiment_sample)
        bwa_dataset = Dataset.objects.create(
                    label=Dataset.TYPE.BWA_ALIGN,
                    type=Dataset.TYPE.BWA_ALIGN,
                    status=Dataset.STATUS.NOT_STARTED)
        sample_alignment.dataset_set.add(bwa_dataset)
        sample_alignment.save()

        # Run the alignment.
        experiment_sample_alignment = align_with_bwa_mem(
                alignment_group, sample_alignment, project=self.project)

        # Check that the run was successful as indicated by the Dataset status.

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

        # compile the tracklist
        compile_tracklist_json(self.reference_genome)

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

                    if not self.reference_genome.project.is_s3_backed():
                        EXPECTED_URL_TEMPLATE = os.path.join(
                                settings.JBROWSE_DATA_URL_ROOT,
                                bwa_align_dataset.filesystem_location)
                        self.assertEqual(EXPECTED_URL_TEMPLATE,
                                track['urlTemplate'])

                    break
            self.assertTrue(found_bam_track)

    def test_compressed_bwa_align(self):
        """Test a single BWA alignment.
        """
        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment',
                reference_genome=self.reference_genome)

        # Create the expected models.
        experiment_sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=self.experiment_sample)
        bwa_dataset = Dataset.objects.create(
                    label=Dataset.TYPE.BWA_ALIGN,
                    type=Dataset.TYPE.BWA_ALIGN,
                    status=Dataset.STATUS.NOT_STARTED)
        experiment_sample_alignment.dataset_set.add(bwa_dataset)
        experiment_sample_alignment.save()

        # Run the alignment.
        align_with_bwa_mem(alignment_group, experiment_sample_alignment,
                project=self.project)

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

        # compile the tracklist
        compile_tracklist_json(self.reference_genome)

        with open(os.path.join(jbrowse_dir, 'trackList.json')) as fh:
            jbrowse_config_json = json.loads(fh.read())
            found_bam_track = False
            for track in jbrowse_config_json['tracks']:
                if ('storeClass' in track and track['storeClass'] ==
                        'JBrowse/Store/SeqFeature/BAM'):
                    found_bam_track = True

                    # Also check that the urlTemplate is correct.
                    if not self.reference_genome.project.is_s3_backed():
                        EXPECTED_URL_TEMPLATE = os.path.join(
                                settings.JBROWSE_DATA_URL_ROOT,
                                bwa_align_dataset.filesystem_location)

                        self.assertEqual(EXPECTED_URL_TEMPLATE,
                                track['urlTemplate'])

                    break
            self.assertTrue(found_bam_track)


class TestAlignmentPieces(TestCase):

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


        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        # Create the expected models.
        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=self.experiment_sample)

        bwa_dataset = copy_and_add_dataset_source(
                sample_alignment,
                dataset_label=Dataset.TYPE.BWA_ALIGN,
                dataset_type=Dataset.TYPE.BWA_ALIGN,
                original_source_location=TEST_DISC_SPLIT_BAM)

        bwa_dataset.status = status=Dataset.STATUS.READY
        bwa_dataset.save()

        self.bwa_dataset = bwa_dataset
        self.sample_alignment = sample_alignment


    def test_get_discordant_read_pairs(self):

        bwa_disc_dataset = get_discordant_read_pairs(self.sample_alignment)
        bwa_disc_dataset_loc = bwa_disc_dataset.get_absolute_location()

        self.assertTrue(os.path.exists(bwa_disc_dataset_loc), (
                "No file at location %s" % bwa_disc_dataset_loc))

        self.assertEqual(bwa_disc_dataset_loc, get_dataset_with_type(
                self.sample_alignment,
                Dataset.TYPE.BWA_DISCORDANT).get_absolute_location())

        # check line counts
        p = subprocess.Popen([SAMTOOLS_BINARY, 'view', bwa_disc_dataset_loc],
                stdout=subprocess.PIPE)
        self.assertEqual(134, sum([1 for line in p.stdout]))


    def test_get_split_reads(self):

        bwa_sr_dataset = get_split_reads(self.sample_alignment)
        bwa_sr_dataset_loc = bwa_sr_dataset.get_absolute_location()

        self.assertTrue(os.path.exists(bwa_sr_dataset_loc), (
                "No file at location %s" % bwa_sr_dataset_loc))

        self.assertEqual(bwa_sr_dataset_loc, get_dataset_with_type(
                self.sample_alignment,
                Dataset.TYPE.BWA_SPLIT).get_absolute_location())

        # check line counts
        p = subprocess.Popen([SAMTOOLS_BINARY, 'view', bwa_sr_dataset_loc],
                stdout=subprocess.PIPE)
        self.assertEqual(3, sum([1 for line in p.stdout]))

    def test_get_read_length(self):
        bam_file = get_dataset_with_type(self.sample_alignment,
                Dataset.TYPE.BWA_ALIGN).get_absolute_location()
        self.assertEqual(70, get_read_length(bam_file))

    def test_get_insert_size(self):
        mean, stdev = get_insert_size_mean_and_stdev(self.sample_alignment)
        self.assertAlmostEqual(mean, 498, delta=2)
        self.assertAlmostEqual(stdev, 50, delta=1)
