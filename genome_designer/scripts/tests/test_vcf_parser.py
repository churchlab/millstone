"""
Tests for vcf_parser.py
"""

import os

from django.test import TestCase
import vcf

from main.models import AlignmentGroup
from main.models import Dataset
from main.models import Project
from main.models import Variant
from main.models import VariantCallerCommonData
from scripts.bootstrap_data import bootstrap_fake_data
from scripts.import_util import copy_and_add_dataset_source
from scripts.import_util import import_reference_genome_from_local_file
from scripts.vcf_parser import parse_alignment_group_vcf
from settings import PWD as GD_ROOT


TEST_FASTA  = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome.fa')

TEST_GENOME_SNPS = os.path.join(GD_ROOT, 'test_data', 'fake_genome_and_reads',
        'test_genome_snps.vcf')


class TestVCFParser(TestCase):

    def setUp(self):
        bootstrap_fake_data()

        self.project = Project.objects.all()[0]

        self.reference_genome = import_reference_genome_from_local_file(
                self.project, 'ref_genome', TEST_FASTA, 'fasta')


    def test_parser(self):
        """Basic tests for the parser.
        """
        VCF_DATATYPE = Dataset.TYPE.VCF_FREEBAYES
        alignment_group = AlignmentGroup.objects.create(
                label='test alignment', reference_genome=self.reference_genome)
        copy_and_add_dataset_source(alignment_group, VCF_DATATYPE,
                VCF_DATATYPE, TEST_GENOME_SNPS)

        # Count the number of records in the vcf file for testing.
        record_count = 0
        with open(TEST_GENOME_SNPS) as fh:
            for record in vcf.Reader(fh):
                record_count += 1

        # Parse the vcf
        parse_alignment_group_vcf(alignment_group, VCF_DATATYPE)

        # There should be one VariantCallerCommonData object for each record.
        self.assertEqual(record_count,
                len(VariantCallerCommonData.objects.filter(
                        reference_genome=self.reference_genome)))

        # There should also be one Variant object for each record.
        self.assertEqual(record_count, len(Variant.objects.filter(
                reference_genome=self.reference_genome)))

        # Spot-check a few variants.
        self.assertEqual(1, len(Variant.objects.filter(
                reference_genome=self.reference_genome,
                position=376)))
        self.assertEqual(1, len(Variant.objects.filter(
                reference_genome=self.reference_genome,
                position=453)))

        # Check false negatives.
        self.assertEqual(0, len(Variant.objects.filter(
                reference_genome=self.reference_genome,
                position=454)))
