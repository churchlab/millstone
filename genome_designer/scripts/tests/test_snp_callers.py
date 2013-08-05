"""
Tests for alignment_pipeline.py
"""

import os

from django.test import TestCase
import vcf

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import get_dataset_with_type
from main.models import Project
from scripts.snp_callers import run_snp_calling_pipeline
from scripts.bootstrap_data import bootstrap_fake_data
from scripts.import_util import add_dataset_to_entity
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import copy_dataset_to_entity_data_dir
from scripts.import_util import import_reference_genome_from_local_file
from settings import PWD as GD_ROOT


TEST_FASTA  = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_FASTQ1 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.1.fq')

TEST_FASTQ2 = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'test_genome_1.snps.simLibrary.2.fq')

TEST_BAM = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'bwa_align.sorted.grouped.realigned.bam')

TEST_BAM_INDEX = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        '38d786f2', 'bwa_align.sorted.grouped.realigned.bam.bai')


class TestSNPCallers(TestCase):

    def setUp(self):
        bootstrap_fake_data()

        # Grab a project.
        self.project = Project.objects.all()[0]

        # Create a ref genome.
        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')


    def test_run_snp_calling_pipeline(self):
        """Test running the pipeline that calls SNPS.
        """
        # Create a new alignment group.
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)

        # Create a sample.
        sample_1 = ExperimentSample.objects.create(
                project=self.project,
                label='sample1')
        ### Add the raw reads
        copy_and_add_dataset_source(sample_1, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, TEST_FASTQ1)
        copy_and_add_dataset_source(sample_1, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, TEST_FASTQ2)

        # Create alignment to the sample.
        sample_alignment = ExperimentSampleToAlignment.objects.create(
                alignment_group=alignment_group,
                experiment_sample=sample_1)
        ### Add alignment data. NOTE: Stored in sample model dir.
        copy_dest = copy_dataset_to_entity_data_dir(sample_1, TEST_BAM)
        copy_dataset_to_entity_data_dir(sample_1, TEST_BAM_INDEX)
        add_dataset_to_entity(sample_alignment, Dataset.TYPE.BWA_ALIGN,
                Dataset.TYPE.BWA_ALIGN, copy_dest)

        # Run the pipeline.
        run_snp_calling_pipeline(alignment_group)

        # Check that the alignment group has a freebayes vcf dataset associated
        # with it.
        vcf_dataset = get_dataset_with_type(alignment_group,
                Dataset.TYPE.VCF_FREEBAYES)
        self.assertIsNotNone(vcf_dataset)

        # Make sure the .vcf file actually exists.
        self.assertTrue(os.path.exists(vcf_dataset.get_absolute_location()))

        # Make sure the vcf is valid by reading it using pyvcf.
        with open(vcf_dataset.get_absolute_location()) as vcf_fh:
            try:
                reader = vcf.Reader(vcf_fh)
                reader.next()
            except:
                self.fail("Not valid vcf")
