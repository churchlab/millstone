"""
Tests for alignment_pipeline.py
"""

import os

from django.test import TestCase

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import Project
from scripts.snp_callers import run_snp_calling_pipeline
from scripts.bootstrap_data import bootstrap_fake_data
from scripts.import_util import import_reference_genome_from_local_file
from settings import PWD as GD_ROOT


TEST_FASTA  = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')


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

        run_snp_calling_pipeline(alignment_group)
