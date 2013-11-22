"""
Test the full alignment to variant calling pipeline.

This is more of an integration tests to make sure that, say, what bwa_mem
produces can be sensibly use by freebayes.
"""

import os

from django.test import TestCase
from django.test.utils import override_settings

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import get_dataset_with_type
from main.models import Project
from main.models import User
from main.models import Variant
from scripts.alignment_pipeline import create_alignment_groups_and_start_alignments
from scripts.snp_callers import run_snp_calling_pipeline
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

class TestFullPipeline(TestCase):

    def setUp(self):
        user = User.objects.create_user('test_username', password='password',
                email='test@example.com')

        # Grab a project.
        self.project = Project.objects.create(title='test project',
                owner=user.get_profile())


    @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS = True,
            CELERY_ALWAYS_EAGER = True, BROKER_BACKEND = 'memory')
    def test_full_pipeline(self):
        """Test running the full pipeline.
        """
        EXPECTED_NUM_VARIANTS = 30

        # We created the test genome with these specs.
        EXPECTED_VARIANT_POSITIONS = [800] # one-indexed.
        while len(EXPECTED_VARIANT_POSITIONS) < EXPECTED_NUM_VARIANTS:
            EXPECTED_VARIANT_POSITIONS.append(
                    EXPECTED_VARIANT_POSITIONS[-1] + 20)

        KNOWN_SUBSTITUTIONS_ROOT = os.path.join(GD_ROOT, 'test_data',
                'test_genome_known_substitutions')

        TEST_GENOME_FASTA = os.path.join(KNOWN_SUBSTITUTIONS_ROOT,
                'test_genome_known_substitutions.fa')

        FAKE_READS_FASTQ1 = os.path.join(KNOWN_SUBSTITUTIONS_ROOT,
                'test_genome_known_substitutions_0.snps.simLibrary.1.fq')

        FAKE_READS_FASTQ2 = os.path.join(KNOWN_SUBSTITUTIONS_ROOT,
                'test_genome_known_substitutions_0.snps.simLibrary.2.fq')

        FAKE_READS_BAM = os.path.join(KNOWN_SUBSTITUTIONS_ROOT,
                'bwa_align.sorted.grouped.realigned.bam')

        FAKE_READS_BAM_INDEX = os.path.join(KNOWN_SUBSTITUTIONS_ROOT,
                'bwa_align.sorted.grouped.realigned.bam.bai')

        # Create a ref genome from the above.
        REFERENCE_GENOME = import_reference_genome_from_local_file(
                self.project, 'test_genome', TEST_GENOME_FASTA, 'fasta')

        # Make sure there are no variants before.
        self.assertEqual(0, len(Variant.objects.filter(
                reference_genome=REFERENCE_GENOME)))

        # Create a sample.
        sample_1 = ExperimentSample.objects.create(
                project=self.project,
                label='sample1')
        ### Add the raw reads
        copy_and_add_dataset_source(sample_1, Dataset.TYPE.FASTQ1,
                Dataset.TYPE.FASTQ1, FAKE_READS_FASTQ1)
        copy_and_add_dataset_source(sample_1, Dataset.TYPE.FASTQ2,
                Dataset.TYPE.FASTQ2, FAKE_READS_FASTQ2)

        # Perform the alignment.
        ref_genome_list = [REFERENCE_GENOME]
        sample_list = [sample_1]
        alignment_group_dict = create_alignment_groups_and_start_alignments(
                'name_placeholder', ref_genome_list, sample_list)
        alignment_group = alignment_group_dict[REFERENCE_GENOME.uid]

        # Run the pipeline.
        run_snp_calling_pipeline(alignment_group)

        # Grab the resulting variants.
        variants = Variant.objects.filter(reference_genome=REFERENCE_GENOME)

        # There should be a Variant object for each record.
        self.assertEqual(EXPECTED_NUM_VARIANTS, len(variants))

        # Check that each variant is accounted for.
        self.assertEqual(set(EXPECTED_VARIANT_POSITIONS),
                set([v.position for v in variants]))
