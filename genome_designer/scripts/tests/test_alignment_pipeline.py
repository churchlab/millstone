"""
Tests for alignment_pipeline.py
"""

import os

from django.test import TestCase

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import get_dataset_with_type
from main.models import ExperimentSample
from main.models import Project
from scripts.alignment_pipeline import align_with_bwa
from scripts.bootstrap_data import bootstrap_fake_data
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import import_reference_genome_from_local_file
from settings import PWD as GD_ROOT


TEST_FASTA  = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')


class TestAlignmentPipeline(TestCase):
    def setUp(self):
        bootstrap_fake_data()

    def test_bwa_align(self):
        """Test a single BWA alignment.
        """
        # Grab any project + ref genome from the database.
        project = Project.objects.all()[0]
        reference_genome = import_reference_genome_from_local_file(
                project, 'ref_genome', TEST_FASTA, 'fasta')

        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=reference_genome)

        # Create a sample.
        experiment_sample = ExperimentSample.objects.create(
                project=project, label='sample1')
        copy_and_add_dataset_source(experiment_sample, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, TEST_FASTQ1)
        copy_and_add_dataset_source(experiment_sample, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, TEST_FASTQ2)

        # Run the alignment.
        experiment_sample_alignment = align_with_bwa(
                alignment_group, experiment_sample)

        # Check that the returned ExperimentSampleToAlignment object has the
        # data stored in its database.
        bwa_align_dataset = get_dataset_with_type(
                experiment_sample_alignment, Dataset.TYPE.BWA_ALIGN)
        bwa_align_dataset_path = bwa_align_dataset.get_absolute_location()
        self.assertTrue(os.path.exists(bwa_align_dataset_path,), (
                "No file at location %s" % bwa_align_dataset_path))
